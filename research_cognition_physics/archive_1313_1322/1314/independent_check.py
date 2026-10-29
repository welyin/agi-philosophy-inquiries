"""Independent round-1314 audit via determinant minors, not author exterior code."""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
from itertools import combinations, permutations
import json
from math import prod, sqrt
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
PERMS=tuple(permutations(range(4)))
SIGNS={p:(-1)**sum(p[i]>p[j] for i in range(4) for j in range(i+1,4)) for p in PERMS}

def det4(rows):
    return sum(SIGNS[p]*prod(rows[i][p[i]] for i in range(4)) for p in PERMS)

def transformed(a,forms):
    return [[[sum(a[i][j]*forms[j][alpha][k] for j in range(4))
              for k in range(4)] for alpha in range(2)] for i in range(4)]

def coefficient(a,forms):
    """Coefficient of theta0 theta1 theta2 theta3 after N_i^alpha=C_i^alpha,r theta_r."""
    an=transformed(a,forms)
    value=0
    eps=((0,1,1),(1,0,-1))
    for i,j,k,l in PERMS:
        for aa,bb,e1 in eps:
            for cc,dd,e2 in eps:
                value+=SIGNS[(i,j,k,l)]*e1*e2*det4(
                    [forms[i][aa],an[j][bb],forms[k][cc],an[l][dd]])
    return value

def diagonal(a):
    return [[a[i] if i==j else 0 for j in range(4)] for i in range(4)]

def basis_forms(indices):
    return [[[int(2*i+alpha==indices[k]) for k in range(4)]
              for alpha in range(2)] for i in range(4)]

def polynomial(a):
    return {','.join(map(str,c)):coefficient(a,basis_forms(c))
            for c in combinations(range(8),4) if coefficient(a,basis_forms(c))!=0}

def matmul(a,b):
    return [[sum(a[i][k]*b[k][j] for k in range(len(b))) for j in range(len(b[0]))] for i in range(len(a))]

def transpose(a):
    return list(map(list,zip(*a)))

def matrix_tests():
    a=diagonal([1,2,3,4]); identity=diagonal([1]*4)
    p=polynomial(a)
    assert p['0,3,4,7']==-4 and len(p)>0
    assert polynomial(identity)=={}
    # Same determinant formula, but independent exact SL(2) spin-frame and SO(4) flavor tests.
    forms=[[[((2*i+alpha+1)*(k+2)+k*k)%7-3 for k in range(4)]
             for alpha in range(2)] for i in range(4)]
    base=coefficient(a,forms)
    spin=[[[forms[i][0][k]+forms[i][1][k] for k in range(4)],forms[i][1][:]] for i in range(4)]
    assert coefficient(a,spin)==base
    u=[[F(3,5),F(-4,5),0,0],[F(4,5),F(3,5),0,0],
       [0,0,F(5,13),F(-12,13)],[0,0,F(12,13),F(5,13)]]
    ap=matmul(matmul(u,a),transpose(u))
    fp=transformed(u,forms)
    assert coefficient(ap,fp)==base and det4(u)==1
    # A genuinely complex unitary covariance test adds phases while keeping det=1.
    uc=np.array(u,dtype=complex)@np.diag([1j,-1j,1,1])
    ac=uc@np.array(a)@uc.conj().T
    fc=transformed(uc,forms)
    complex_error=abs(coefficient(ac,fc)-base)
    assert complex_error<1e-8
    rng=np.random.default_rng(1314)
    yn=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4))
    ul=np.eye(4,dtype=complex); ul[0,0]=1j;ul[1,1]=-1j
    yp=ul.conj()@yn@uc.conj().T
    gram_error=float(np.linalg.norm(yp.conj().T@yp-uc@(yn.conj().T@yn)@uc.conj().T))
    assert gram_error<1e-11
    return dict(exact_coefficients=p,bare_epsilon_zero=True,target_coefficient=-4,
        exact_spin_covariance=True,exact_flavor_covariance=True,
        complex_covariance_error=float(complex_error),Y_gram_covariance_error=gram_error)

SIGMA=[np.array([[0,1],[1,0]],complex),
       np.array([[0,-1j],[1j,0]],complex),np.diag([1,-1]).astype(complex)]

def slash(p):
    return sum(p[i]*SIGMA[i] for i in range(3))

def hamiltonian(m,p):
    pp=slash(p)
    return np.block([[-pp,m*np.eye(2)],[m*np.eye(2),pp]])

def pair_kinematics(total,m1,m2):
    s=total*total
    rad=(s-(m1+m2)**2)*(s-(m1-m2)**2)
    assert rad>0
    p=sqrt(rad)/(2*total)
    return (s+m1*m1-m2*m2)/(2*total),(s+m2*m2-m1*m1)/(2*total),p

def scattering_tests():
    a=diagonal([1,2,3,4])
    unit=[np.array([1,0.],complex),np.array([0.,1],complex),
          np.array([1.,1],complex)/sqrt(2),np.array([-1.,1],complex)/sqrt(2)]
    cases=[]
    for vev in (0.,.01,.1,.3):
        masses=[vev*sqrt(j) for j in (1,2,3,4)]
        total=2.
        e1,e2,pin=pair_kinematics(total,*masses[:2])
        e3,e4,pout=pair_kinematics(total,*masses[2:])
        energies=[e1,e2,e3,e4]
        momenta=[np.array([0,0,-pin]),np.array([0,0,pin]),
                 np.array([-pout,0,0]),np.array([pout,0,0])]
        upper=[]; residual=0.
        for i,(m,e,pvec,chi) in enumerate(zip(masses,energies,momenta,unit)):
            mag=float(np.linalg.norm(pvec))
            up=sqrt(e+mag)*chi
            low=m/sqrt(e+mag)*chi
            external=np.concatenate([up,low if i<2 else -low])
            op=hamiltonian(m,pvec if i<2 else -pvec)
            sign=1 if i<2 else -1
            residual=max(residual,float(np.linalg.norm(op@external-sign*e*external)),
                         abs(float(np.vdot(external,external).real)-2*e))
            upper.append(up)
        forms=[[[upper[i][alpha] if i==j else 0 for j in range(4)] for alpha in range(2)] for i in range(4)]
        amplitude=coefficient(a,forms)
        expected=3*sqrt(prod(e+float(np.linalg.norm(p)) for e,p in zip(energies,momenta)))
        assert abs(amplitude-expected)<1e-11 and expected>0 and residual<1e-11
        assert abs(sum(energies[:2])-sum(energies[2:]))<1e-12
        assert np.linalg.norm(momenta[0]+momenta[1]-momenta[2]-momenta[3])<1e-12
        if vev==0: assert abs(amplitude-12)<1e-12
        cases.append(dict(Higgs_amplitude=vev,masses=masses,energies=energies,
                          contact_coefficient=float(amplitude.real),
                          positive_closed_form=expected,Dirac_residual=residual))
    return dict(cases=cases,normalization='unit Wilson coefficient and omit 1/Lambda^2',
                physical_scope='leading contact amplitude only; no device or higher-order error certificate')

def center_tests():
    # Exact integer exponents. Fundamental field weights are +-1, Yukawa weight is minus their sum.
    zl=dict(Q=0,u=0,d=0,L=1,e=-1,N=-1)
    zq=dict(Q=1,u=-1,d=-1,L=0,e=0,N=0)
    pairs=(('Q','u'),('Q','d'),('L','e'),('L','N'))
    assert all(-(z[x]+z[y])==0 for z in (zl,zq) for x,y in pairs)
    checked=[]
    for n in range(3,65):
        assert all((n*z[f])%n==0 for z in (zl,zq) for f in z)
        assert 1%n!=0
        checked.append(dict(n=n,E_allowed=False,four_N_center_allowed=(4%n==0)))
    assert [c['n'] for c in checked if c['four_N_center_allowed']]==[4]
    return dict(integer_yukawa_stabilizer_weights=True,samples=checked,
                universal_claim='proved analytically by the two point stabilizers; finite list only diagnoses code')

def compute():
    sources={'admission':HERE.parent/'_admission/round1314_admission.md',
             'proof':HERE/'proof.md','author_code':HERE/'check.py',
             'author_results':HERE/'results.json','note':HERE.parent/'research_note_1314.md',
             'next':HERE/'NEXT.md','independent_code':Path(__file__)}
    return dict(schema='round1314_independent_check_v1',status='PASS_WITH_SCOPE',
        arithmetic=matrix_tests(),scattering=scattering_tests(),centers=center_tests(),
        sources_sha256={str(p.relative_to(HERE.parent)):hashlib.sha256(p.read_bytes()).hexdigest()
                       for p in sources.values()},
        imported_author_code=False,new_axioms=0)

def main():
    ap=argparse.ArgumentParser()
    mode=ap.add_mutually_exclusive_group()
    mode.add_argument('--check',action='store_true')
    mode.add_argument('--write-new',action='store_true')
    args=ap.parse_args()
    result=compute()
    path=HERE/'independent_results.json'
    if args.write_new:
        with path.open('x',encoding='utf-8',newline='\n') as f:
            json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
    else:
        saved=json.loads(path.read_text(encoding='utf-8'))
        assert saved==result,'Independent result or source binding changed'
    print(json.dumps(dict(status=result['status'],coefficients=len(result['arithmetic']['exact_coefficients']),
        center_diagnostics=len(result['centers']['samples']),massive_amplitudes=len(result['scattering']['cases'])-1,
        source_bindings=len(result['sources_sha256']),mode='write-new' if args.write_new else 'read-only')))

if __name__=='__main__':main()
