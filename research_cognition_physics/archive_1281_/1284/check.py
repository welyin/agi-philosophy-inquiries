#!/usr/bin/env python3
"""1284: direct Legendre, constraint drift and local electric completion.
No author-independent checker is imported. --write is exclusive first output;
the default recomputes and compares the stored result without writing.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
I = np.eye(2, dtype=complex)
X = np.array([[0,1],[1,0]], dtype=complex)
Y = np.array([[0,-1j],[1j,0]], dtype=complex)
Z = np.diag([1,-1]).astype(complex)
R = np.array([X,Y,Z])/2
ZERO = np.zeros_like(R)
K_GOOD = np.array([[(a@b+b@a)/2 for b in R] for a in R])
FD = 1.e-5


def real(x):
    assert abs(np.imag(x)) < 1.e-8, x
    return float(np.real(x))


def combine(a, mats):
    return np.einsum('a,aij->ij', a, mats)


def k00(a, k):
    return np.einsum('a,b,abij->ij', a, a, k)


def charge(phi, pi, mats):
    return np.array([real(1j*(pi@m@phi-phi.conj()@m@pi.conj())) for m in mats])


def temporal_l(phi, vel, a, c, s, k):
    c0, s0 = combine(a,c), combine(a,s)
    return real(np.vdot(vel,vel)+vel.conj()@(s0-1j*c0)@phi
                +phi.conj()@(s0+1j*c0)@vel+phi.conj()@k00(a,k)@phi)


def original_temporal_l(phi, vel, a, adot, c, s, k):
    return temporal_l(phi,vel,a,c,ZERO,k)-sum(
        adot[j]*real(np.vdot(phi,s[j]@phi)) for j in range(3))


def momentum_from_l(phi, vel, a, c, s, k):
    pi = np.zeros(2, dtype=complex)
    for j in range(2):
        d = np.zeros(2, dtype=complex); d[j] = FD
        dr = (temporal_l(phi,vel+d,a,c,s,k)-temporal_l(phi,vel-d,a,c,s,k))/(2*FD)
        di = (temporal_l(phi,vel+1j*d,a,c,s,k)-temporal_l(phi,vel-1j*d,a,c,s,k))/(2*FD)
        pi[j] = (dr-1j*di)/2
    return pi


def legendre(phi, pi, a, c, s, k):
    vel = pi.conj()-(combine(a,s)-1j*combine(a,c))@phi
    return 2*float(np.real(pi@vel))-temporal_l(phi,vel,a,c,s,k)


def h_expanded(phi, pi, a, c, s, k):
    c0,s0=combine(a,c),combine(a,s)
    qs=np.array([2*float(np.real(pi@m@phi)) for m in s])
    mm=(s0+1j*c0)@(s0-1j*c0)-k00(a,k)
    return real(pi@pi.conj())+float(a@(charge(phi,pi,c)-qs))+real(np.vdot(phi,mm@phi))


def primary_from_l(phi, pi, a, c, s, k):
    # Put the initial data on the prescribed Gauss surface: DE = Q_R.
    de=charge(phi,pi,R)
    ans=[]
    for j in range(3):
        da=np.zeros(3); da[j]=FD
        hp=legendre(phi,pi,a+da,c,s,k)-float((a+da)@de)
        hm=legendre(phi,pi,a-da,c,s,k)-float((a-da)@de)
        ans.append((hp-hm)/(2*FD))
    return np.array(ans)


def primary_formula(phi, pi, a, c, s, k):
    c0,s0=combine(a,c),combine(a,s)
    ans=charge(phi,pi,c)-charge(phi,pi,R)
    for j in range(3):
        dm=(s[j]+1j*c[j])@(s0-1j*c0)+(s0+1j*c0)@(s[j]-1j*c[j])
        dm-=2*np.einsum('b,bij->ij',a,k[j])
        ans[j]+=-2*float(np.real(pi@s[j]@phi))+real(np.vdot(phi,dm@phi))
    return ans


def spatial_energy(phi, pi, grad, aa, c, s, k, mass):
    val=real(pi@pi.conj())+real(np.vdot(phi,mass@phi))+0.1*real(np.vdot(phi,phi))**2
    for d in range(3):
        v=grad[d]; cd=combine(aa[d],c); sd=combine(aa[d],s)
        val+=real(np.vdot(v,v)+1j*(phi.conj()@cd@v-v.conj()@cd@phi)
                  +phi.conj()@k00(aa[d],k)@phi
                  +v.conj()@sd@phi+phi.conj()@sd@v)
    return val


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def calculate():
    rng=np.random.default_rng(1284)
    residuals={k:0.0 for k in ['momentum_from_original_l','legendre_transform',
         'boundary_lagrangian','boundary_one_form','primary_derivative',
         'correct_primary','correct_spatial_direction','electric_gauss']}
    def keep(k,x): residuals[k]=max(residuals[k],float(np.max(np.abs(x))))
    def herm():
        z=(rng.normal(size=(2,2))+1j*rng.normal(size=(2,2)))/5
        return (z+z.conj().T)/2
    for case in range(12):
        phi=(rng.normal(size=2)+1j*rng.normal(size=2))/3
        vel=(rng.normal(size=2)+1j*rng.normal(size=2))/3
        a=rng.normal(size=3)/3; adot=rng.normal(size=3)/3
        c=np.array([herm() for _ in range(3)])
        s=np.array([herm() for _ in range(3)])
        k=np.empty((3,3,2,2),dtype=complex)
        for j in range(3):
            for l in range(j,3): k[j,l]=k[l,j]=herm()
        c0,s0=combine(a,c),combine(a,s)
        pi=vel.conj()+phi.conj()@(s0+1j*c0)
        keep('momentum_from_original_l',momentum_from_l(phi,vel,a,c,s,k)-pi)
        keep('legendre_transform',legendre(phi,pi,a,c,s,k)-h_expanded(phi,pi,a,c,s,k))
        ob=np.array([real(np.vdot(phi,m@phi)) for m in s])
        db=float(adot@ob)+2*float(np.real(phi.conj()@s0@vel))
        keep('boundary_lagrangian',temporal_l(phi,vel,a,c,s,k)-original_temporal_l(phi,vel,a,adot,c,s,k)-db)
        # Direct canonical one-form comparison, including the old A0 momentum.
        dp=(rng.normal(size=2)+1j*rng.normal(size=2))/3
        da=rng.normal(size=3)/3
        pi_old=vel.conj()+1j*phi.conj()@c0
        theta_change=2*float(np.real((pi-pi_old)@dp))+float(ob@da)
        def boundary(aa,pp): return real(np.vdot(pp,combine(aa,s)@pp))
        db_fd=(boundary(a+FD*da,phi+FD*dp)-boundary(a-FD*da,phi-FD*dp))/(2*FD)
        keep('boundary_one_form',theta_change-db_fd)
        keep('primary_derivative',primary_from_l(phi,pi,a,c,s,k)-primary_formula(phi,pi,a,c,s,k))
        keep('correct_primary',primary_from_l(phi,pi,a,R,ZERO,K_GOOD))
        # Local gauge jet direction; this includes grad(epsilon), not a constant rotation alone.
        grad=(rng.normal(size=(3,2))+1j*rng.normal(size=(3,2)))/3
        aa=rng.normal(size=(3,3))/3; eps=rng.normal(size=3)/3
        deps=rng.normal(size=(3,3))/3; er=combine(eps,R)
        dphi=1j*er@phi; dpi=-1j*pi@er
        dgrad=np.array([1j*combine(deps[d],R)@phi+1j*er@grad[d] for d in range(3)])
        daa=np.array([deps[d]+np.cross(aa[d],eps) for d in range(3)])
        def es(t):
            return spatial_energy(phi+t*dphi,pi+t*dpi,grad+t*dgrad,aa+t*daa,R,ZERO,K_GOOD,I)
        keep('correct_spatial_direction',(es(FD)-es(-FD))/(2*FD))
    # Three distinct primary failures; each datum is locally Gauss-completable.
    phi=np.array([1,0],dtype=complex); az=np.zeros(3)
    c_bad=R.copy(); c_bad[0]+=.25*I
    s_bad=ZERO.copy(); s_bad[0]=.25*I
    k_bad=K_GOOD.copy(); k_bad[0,0]+=.25*I
    witnesses={
        'wrong_C':primary_from_l(phi,np.array([1j,0]),az,c_bad,ZERO,K_GOOD)[0],
        'wrong_S':primary_from_l(phi,np.array([1,0],dtype=complex),az,R,s_bad,K_GOOD)[0],
        'wrong_K':primary_from_l(phi,np.zeros(2,dtype=complex),np.array([1.,0,0]),R,ZERO,k_bad)[0]
    }
    for value in witnesses.values(): assert abs(value+.5)<1.e-8, value
    # A positive, noninvariant mass potential causes secondary Gauss drift.
    phiv=np.array([1,1j],dtype=complex); mv=np.diag([1.,2.])
    pdot=-phiv.conj()@mv
    gauss_drift=(charge(phiv,FD*pdot,R)-charge(phiv,-FD*pdot,R))/(2*FD)
    dphiv=1j*R[0]@phiv
    dv=2*float(np.real(phiv.conj()@mv@dphiv))
    assert abs(dv-1)<1.e-12 and abs(gauss_drift[0]+1)<1.e-12
    witnesses['wrong_positive_V_Gauss_drift']=float(gauss_drift[0])
    witnesses['wrong_positive_V_gauge_variation']=dv
    # Explicit non-Abelian local completion with a nonzero connection and field energy.
    aconn=1/3; q=np.array([1.,2.,3.]); avec=np.array([0.,0.,aconn])
    pic=np.array([-3j,-2-1j]); assert np.max(abs(charge(phi,pic,R)-q))<1.e-12
    def e_and_derivative(x):
        c=math.cos(aconn*x); s=math.sin(aconn*x)
        e=np.array([(q[0]*s+q[1]*(1-c))/aconn,(q[1]*s-q[0]*(1-c))/aconn,q[2]*x])
        ep=np.array([q[0]*c+q[1]*s,q[1]*c-q[0]*s,q[2]])
        return e,ep
    for x in np.linspace(0,1,17):
        ee,ep=e_and_derivative(float(x))
        keep('electric_gauss',ep+np.cross(avec,ee)-q)
        assert np.linalg.norm(ee)<=np.linalg.norm(q)*x+1.e-12
    nodes,weights=np.polynomial.legendre.leggauss(24)
    energy=sum(w*np.dot(e_and_derivative((x+1)/2)[0],e_and_derivative((x+1)/2)[0])/4 for x,w in zip(nodes,weights))
    assert 0<float(energy)<7/3
    # Deleted quantifiers: all defects disappear on the zero matter trajectory,
    # while this does not remove any of the nonzero local witnesses above.
    keep_zero=primary_formula(np.zeros(2,dtype=complex),np.zeros(2,dtype=complex),np.ones(3),c_bad,s_bad,k_bad)
    assert np.array_equal(keep_zero,np.zeros(3))
    for key,value in residuals.items(): assert value<2.e-8,(key,value)
    srcs=[
      'archive_1009_1043/research_note_1016.md',
      'archive_1009_1043/research_note_1027.md',
      'archive_1009_1043/research_note_1029.md',
      'archive_894_922/research_note_920.md']
    return {
      'schema':'research_1284_author_v1','round':1284,'verdict':'PASS_WITH_SCOPE','author_groups':1,
      'method':'direct original-L and boundary Legendre checks; independent electric ODE completion',
      'diagnostics':{'noncommuting_vertex_cases':12,'finite_difference_step':FD,
        'maximum_residuals':residuals,'finite_witnesses':witnesses,
        'witness_exact_targets':{'primary_each':'-1/2','secondary_V':'-1'},
        'electric_completion':{'connection':[0,0,'1/3'],'charge':[1,2,3],
             'unit_box_energy_upper':'7/3','energy_quadrature_diagnostic':float(energy),
             'gauss_points':17,'energy_upper_is_analytic':True},
        'zero_matter_trajectory_blindness':True},
      'scope':{'R_representation_is_input':True,'full_primary_secondary_chain':True,
        'arbitrary_local_classical_Gauss_data':True,'fixed_boundary_canonical_dictionary':True,
        'not_general_Dirac_conjecture':True,'no_new_axiom_adoption':True,
        'no_arbitrary_global_boundary_completion':True,'no_quantum_jet_preparation_certified':True,
        'no_new_device':True,'floating_diagnostics_do_not_replace_analytic_quantifiers':True},
      'sources_sha256':{name:digest(ROOT/name) for name in srcs}}


def equal(a,b,path='root'):
    if isinstance(a,dict):
        assert isinstance(b,dict) and a.keys()==b.keys(),path
        for k in a: equal(a[k],b[k],path+'.'+k)
    elif isinstance(a,list):
        assert isinstance(b,list) and len(a)==len(b),path
        for j,(x,y) in enumerate(zip(a,b)): equal(x,y,path+f'[{j}]')
    elif isinstance(a,float):
        assert isinstance(b,(float,int)) and math.isfinite(a) and abs(a-b)<=1.e-10*max(1,abs(a),abs(b)),(path,a,b)
    else: assert a==b,(path,a,b)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args(); result=calculate(); output=HERE/'results.json'
    if args.write:
        with output.open('x',encoding='utf-8',newline='\n') as f:
            json.dump(result,f,ensure_ascii=False,indent=2,sort_keys=True); f.write('\n')
        mode='FIRST_WRITE'
    else:
        equal(result,json.loads(output.read_text(encoding='utf-8'))); mode='READ_ONLY_MATCH'
    print(json.dumps({'verdict':result['verdict'],'mode':mode,'author_groups':1,
          'max_residual':max(result['diagnostics']['maximum_residuals'].values()),
          'primary_witnesses':'C,S,K each -1/2','secondary_witness':'-1',
          'electric_energy_upper':'7/3'},ensure_ascii=False,sort_keys=True))

if __name__=='__main__': main()
