"""1190 independent review: polynomial directional derivatives and 2x2 covariance.
No author program is imported; no finite Fock diagonalization is used.
Default mode recomputes and compares the saved review receipt.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import cmath
import hashlib
import json
import math
import numpy as np

ROOT=Path(__file__).resolve().parent

def polynomial_value(poly,phi,scale):
    return sum(c*phi**i*scale**j for (i,j),c in poly.items())

def directional(poly,phi,scale,dphi,dscale):
    total=F(0)
    for (i,j),c in poly.items():
        if i: total+=c*i*phi**(i-1)*scale**j*dphi
        if j: total+=c*j*phi**i*scale**(j-1)*dscale
    return total

def calculate():
    # The reduced density, evaluated before defining any named source.
    # phi and tetrad scale are independent polynomial variables.
    zeta,v,m,q,X=F(5,2),F(3,2),F(7,5),F(4,7),F(-2,3)
    source_rows=[]
    for phi in (F(5,4),F(7,3),F(3)):
        for c in (F(0),F(2,5),F(1)):
            poly={(3,2):-zeta*X/2,
                  (3,4):-zeta*m*m*(1-c)*q*q/2,
                  (5,4):-zeta*m*m*c*q*q/(2*v*v)}
            # e=scale*e0: trace is the scale derivative divided by scale^3.
            trace=directional(poly,phi,F(1),F(0),F(1))
            je=directional(poly,phi,F(1),F(1),F(0))
            ja=directional(poly,phi,F(1),F(1),1/phi)
            fixed_h_density=directional(poly,phi,F(1),F(1),-1/(2*phi))
            assert ja==je+trace/phi
            assert fixed_h_density==ja-3*trace/(2*phi)
            # Fixed h has metric scale phi/v at this evaluation point.
            # Differentiate its Einstein-frame action with that metric held fixed.
            hs2=phi/v
            fh=1-c+c*phi*phi/(v*v)
            fprime=2*c*phi/(v*v)
            jh_direct=-F(1,2)*(2*zeta*v*phi*X/hs2
                +zeta*v*v*m*m*(fh+phi*fprime)*q*q)
            jh=(v/phi)**2*fixed_h_density
            assert jh==jh_direct
            source_rows.append({'phi':str(phi),'c':str(c),
                'trace':str(trace),'j_e':str(je),'j_A':str(ja),
                'j_h':str(jh),'residual':'0'})
    # Quadratic classical flow acts exactly on the full-Fock Gaussian covariance.
    phi=F(121,100)
    r=F(1331,1000) # M*omega; omega*t=pi/2 exactly.
    S=np.array([[0,1/float(r)],[-float(r),0]],dtype=float)
    J=np.array([[0,1],[-1,0]],dtype=float)
    cov=S@(.5*np.eye(2))@S.T
    n_cov=(float(np.trace(cov))-1)/2
    n_exact=(r-1/r)**2/4
    assert abs(n_cov-float(n_exact))<2e-15
    assert abs(float(np.linalg.det(cov))-.25)<2e-15
    assert np.max(np.abs(S.T@J@S-J))<2e-15
    # Integrate original Hermite-2 wavefunction against the evolved Gaussian.
    # psi(x)=(Re a/pi)^1/4 exp(-a*x*x/2), a obtained solely from covariance.
    a=complex(1,-2*float(cov[0,1]))/(2*float(cov[0,0]))
    amp2=a.real**.25/math.sqrt(8)*(4/(1+a)-2)*cmath.sqrt(2/(1+a))
    p2=abs(amp2)**2
    p2_n=float(n_exact)/(2*(1+float(n_exact))**1.5)
    assert abs(p2-p2_n)<2e-15
    # Exact positive squaring proves P2 > 3/100, no floating margin.
    strict_margin=n_exact*n_exact-F(3,50)**2*(1+n_exact)**3
    assert strict_margin>0
    energy_form=np.diag([float(phi),1/float(phi)**2])
    e0=float(np.trace(energy_form@(.5*np.eye(2)))/2)
    e1=float(np.trace(energy_form@cov)/2)
    assert abs(e1-e0)<2e-15
    # Same eigenfrequency with original basis held fixed is a pure rotation.
    S_diag=np.array([[0,1],[-1,0]],dtype=float)
    cov_diag=S_diag@(.5*np.eye(2))@S_diag.T
    assert np.array_equal(cov_diag,.5*np.eye(2))
    # Source on instantaneous number states: q^2 and p^2 matrix elements.
    # q^2_20=+sqrt(2)/2, p^2_20=-sqrt(2)/2 at phi=1.
    q20=math.sqrt(2)/2
    p20=-q20
    J20=-p20+q20/2
    assert abs(J20-3*math.sqrt(2)/4)<1e-15
    family=[]
    for c in (F(0),F(1,2),F(1)):
        alpha=c-F(1,2)
        # f=1-c+c*phi^2 is positive on phi>0 for every c in [0,1].
        assert 0<=c<=1
        source02=(F(3,2)+c)*q20
        family.append({'c':str(c),'H_at_phi_1_same':True,
            'gap_at_phi_1':'1','d_log_gap_at_phi_1':str(alpha),
            'source_02':float(source02),'positivity_argument':'convex sum of 1 and phi^2'})
    author=json.loads((ROOT/'results.json').read_text(encoding='utf-8-sig'))
    assert abs(p2-author['full_Fock_probability']['P2'])<2e-14
    assert str(n_exact)==author['full_Fock_probability']['occupation_exact']
    assert abs(J20-author['complete_quantum_source']['zero_to_two'])<2e-14
    assert [r['d_log_gap_at_phi_1'] for r in family]==[
        r['d_log_gap_at_reference'] for r in author['allowed_material_family']]
    hashes={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
        for name in ('admission.md','proof.md','sources.md','check.py','results.json','reviewer_check.py')}
    return {'status':'PASS_LIMITED_MECHANISM','independent_algorithm':
        'Exact polynomial directional derivatives; 2x2 symplectic covariance; Gaussian-Hermite overlap; no author code imports',
        'double_source_rows':source_rows,
        'quantum_diagnostic':{'symplectic_residual':float(np.max(np.abs(S.T@J@S-J))),
            'covariance':cov.tolist(),'n_exact':str(n_exact),'n_covariance':n_cov,
            'P2_Hermite_integral':p2,'P2_Bogoliubov_comparison':p2_n,
            'rational_strict_margin':str(strict_margin),'P2_strict_lower':'3/100',
            'frozen_original_basis_P2':'0','fixed_H_energy_residual':abs(e1-e0),
            'source_02':J20},
        'material_family':family,
        'scope':{'fixed_background_quantum_diagnostic':True,
            'parameter_switching_implemented':False,'actual_clock_instrument_certified':False,
            'closed_geometry_backreaction_solved':False,'all_axioms_countermodel':False},
        'reviewed_sha256':hashes}

def compare(a,b,path='root'):
    if isinstance(a,dict):
        assert isinstance(b,dict) and set(a)==set(b),path
        for k in a:compare(a[k],b[k],path+'.'+k)
    elif isinstance(a,list):
        assert isinstance(b,list) and len(a)==len(b),path
        for i,(x,y) in enumerate(zip(a,b)):compare(x,y,path+f'[{i}]')
    elif isinstance(a,float):
        assert isinstance(b,(int,float)) and math.isclose(a,b,rel_tol=1e-11,abs_tol=2e-14),path
    else:assert a==b,path

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    data=calculate()
    out=ROOT/'review_checks.json'
    if args.write:
        with out.open('x',encoding='utf-8',newline='\n') as f:
            json.dump(data,f,ensure_ascii=False,indent=2);f.write('\n')
        print('PASS / EXCLUSIVE_REVIEW_RECEIPT_CREATED')
    else:
        compare(data,json.loads(out.read_text(encoding='utf-8-sig')))
        print('PASS / READ_ONLY_REVIEW_MATCH')
