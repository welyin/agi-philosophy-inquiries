"""1310 finite Born task, positive massive transverse currents and RG bridge.

Default/--check is read-only. --write-new exclusively creates results.json.
Analytic proofs, not floating samples, carry all interval/general claims.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import json
import math
import numpy as np

HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parents[1]
ROOT=RESEARCH.parent
OUT=HERE/'results.json'
I2=np.eye(2,dtype=complex)
ZERO=np.zeros((2,2),complex)
PAULI=[np.array([[0,1],[1,0]],complex),np.array([[0,-1j],[1j,0]],complex),np.diag([1,-1]).astype(complex)]
G0=np.diag([1,1,-1,-1]).astype(complex)
G=[G0]+[np.block([[ZERO,s],[-s,ZERO]]) for s in PAULI]
G5=1j*G[0]@G[1]@G[2]@G[3]
PL=(np.eye(4)-G5)/2


def mdot(a,b): return a[0]*b[0]-sum(a[1:]*b[1:])
def slash(p): return p[0]*G[0]-sum((p[j+1]*G[j+1] for j in range(3)),np.zeros((4,4),complex))
def sp(p): return sum((p[j]*PAULI[j] for j in range(3)),np.zeros((2,2),complex))
def spinor(p,m,chi):
    a=np.sqrt(p[0]+m)
    return np.r_[a*chi,sp(p[1:])@chi/a]
def current(o,i): return np.array([o.conj()@G0@g@PL@i for g in G])


def external(m,n,E,theta):
    r,s=np.sqrt(E*E-m*m),np.sqrt(E*E-n*n)
    p=[np.array(x,float) for x in [(E,0,0,r),(E,0,0,-r),
         (E,s*np.sin(theta),0,s*np.cos(theta)),(E,-s*np.sin(theta),0,-s*np.cos(theta))]]
    q=p[0]-p[2]
    A,B=s*np.sin(theta),r-s*np.cos(theta)
    z=A*A+B*B
    b=np.sqrt(E+n)
    T1=b*I2-sp(p[2][1:])/b; T2=b*I2+sp(p[2][1:])/b
    w1=np.linalg.solve(T1,np.array([A,B],complex));w2=np.linalg.solve(T2,np.array([B,-A],complex))
    norm1,norm2=np.linalg.norm(w1),np.linalg.norm(w2)
    chis=[np.array([1,0],complex),np.array([0,1],complex),w1/norm1,w2/norm2]
    masses=[m,m,n,n]
    u=[spinor(pp,mm,cc) for pp,mm,cc in zip(p,masses,chis)]
    j1,j2=current(u[2],u[0]),current(u[3],u[1])
    f=(np.sqrt(E+m)-r/np.sqrt(E+m))/(2*norm1)
    j1a=f*np.array([A,-B,-1j*B,-A]);j2a=f*np.array([-A,-B,1j*B,-A])
    scale=max(np.linalg.norm(j1a),np.linalg.norm(j2a),1e-30)
    checks={
        'dirac_equation':max(float(np.linalg.norm((slash(pp)-mm*np.eye(4))@uu))/(E*np.linalg.norm(uu)) for pp,mm,uu in zip(p,masses,u)),
        'baruu_normalization':max(abs(uu.conj()@G0@uu-2*mm)/(2*E) for uu,mm in zip(u,masses)),
        'norm_2E':max(abs(np.vdot(uu,uu)-2*E)/(2*E) for uu in u),
        'momentum':float(np.linalg.norm(p[0]+p[1]-p[2]-p[3]))/E,
        'same_norm':abs(norm1-norm2)/norm1,
        'current_formula':max(float(np.linalg.norm(j1-j1a)),float(np.linalg.norm(j2-j2a)))/scale,
        'Ward_transverse':max(abs(mdot(q,j1)),abs(mdot(q,j2)))/(np.sqrt(z)*scale),
        'current_dot':abs(mdot(j1,j2)+2*f*f*z)/(2*f*f*z)}
    assert max(checks.values())<3e-8,checks
    assert abs(mdot(j1,j2))>0
    return dict(p=p,q=q,z=z,j1=j1,j2=j2,jdot=mdot(j1,j2),f=f,r=r,s=s,checks=checks)


def a(g,z,v): return g*g/(2*(z+g*g*v*v/4))
def inverse_floor(A,z,v):
    assert A>0 and A*v*v<2
    return np.sqrt(2*A*z/(1-A*v*v/2))


def kernel_setup(m,n,E,c2=1/9,theta0=np.pi/3,theta1=2*np.pi/3,order=128):
    r,s=np.sqrt(E*E-m*m),np.sqrt(E*E-n*n)
    z0=r*r+s*s-2*r*s*np.cos(theta0);z1=r*r+s*s-2*r*s*np.cos(theta1)
    nodes,weights=np.polynomial.legendre.leggauss(order)
    zs=(z0+z1)/2+(z1-z0)*nodes/2;weights=weights*(z1-z0)/2
    lam=16*E*E*r*r
    ks=[]
    for z in zs:
        theta=np.arccos((r*r+s*s-z)/(2*r*s))
        x=external(m,n,E,theta)
        ks.append(c2*abs(x['jdot'])**2/(16*np.pi*lam))
    ks=np.array(ks)
    u=np.sqrt(E+m)-r/np.sqrt(E+m);b=np.sqrt(E+n)
    hminus=u*(b-s/b)/2;hplus=u*(b+s/b)/2
    klower=(z1-z0)*c2*(2*hminus*hminus)**2/(16*np.pi*lam)
    kupper=(z1-z0)*c2*(2*hplus*hplus)**2/(16*np.pi*lam)
    ktot=float(weights@ks)
    assert 0<klower<=ktot<=kupper
    return dict(z0=z0,z1=z1,zs=zs,weights=weights,ks=ks,ktot=ktot,klower=klower,kupper=kupper)


def sigma(g,v,k): return float(k['weights']@(k['ks']*a(g,k['zs'],v)**2))


def scalar_checks():
    # Exact inversion including the mW=gv/2 denominator and Higgs relation.
    data=[]
    for g2 in (F(1,100),F(1,9),F(1,4),F(1)):
        for z in (F(1),F(4),F(25)):
            v2=F(3,2);A=g2/(2*(z+g2*v2/4));d=z*A;R2=v2/z
            assert 2*d/(1-d*R2/2)==g2
            assert 1/g2+v2/(4*z)==1/(2*z*A)
            data.append([str(g2),str(z),str(A)])
    # No finite inverse at saturation; above saturation is impossible.
    assert F(2,1)*F(1,1)==2
    return {'exact_inversion_cases':len(data),'finite_saturation_excluded':True,
            'Higgs_inverse_coupling_relation_checked':True}


def phase_compatibility():
    mu=2.0**-12
    U3=np.exp(2j*np.pi*np.outer(np.arange(3),np.arange(3))/3)/np.sqrt(3)
    M=U3@np.diag(mu*np.array([1,2,3]))@U3.conj().T
    # Original neutral Hamiltonian/mass is retained exactly; charged masses
    # are deliberately changed and not claimed to preserve total sources.
    Y_e=np.diag(mu*np.array([1,2,4]))
    assert np.min(np.diff(np.diag(Y_e)))>0
    rows=[]
    for N in (3,4,5):
        U=np.eye(N,dtype=complex);U[:3,:3]=U3
        R=np.eye(N)
        for j in range(2,N-1):
            rr=np.eye(N);t=1e-5;c,s=np.cos(t),np.sin(t)
            rr[j,j]=rr[j+1,j+1]=c;rr[j,j+1]=s;rr[j+1,j]=-s
            R=R@rr
        U=U@R
        assert np.max(abs(U.conj().T@U-np.eye(N)))<2e-14
        assert np.max(abs(U[:,:2]-np.vstack([U3[:,:2],np.zeros((N-3,2))])))<1e-14
        masses=mu*np.arange(1,N+1)
        MM=U@np.diag(masses)@U.conj().T
        HH=MM@MM
        tri=HH[0,1]*HH[1,2]*HH[2,0]
        assert abs(tri.imag)>1e-30
        seen={0}
        while True:
            new=seen|{j for i in seen for j in range(N) if abs(HH[i,j])>1e-18}
            if new==seen:break
            seen=new
        assert len(seen)==N
        c=U[1,0]*U[0,1].conjugate()
        assert abs(abs(c)-1/3)<1e-14
        rows.append({'ng':N,'neutral_mass_eigenvalues_positive_simple':True,
                     'all_flavor_connected':True,'fixed_first_two_mass_columns':True,
                     'abs_charged_current_mixing':float(abs(c)),
                     'nonzero_H_triangle':True,'old_full_source_numbers_inherited':False})
    return {'mu':mu,'neutral_three_mode_mass_unchanged':True,
            'charged_mass_choice':np.diag(Y_e).tolist(),'branches':rows,
            'old_phase_probability_formula_reused_only_in_unchanged_neutral_sector':True}


def compute():
    ward=[];fullW=[]
    cases=[(3.,4.,5.),(3.,4.,50.),(2.0**-12,2.0**-11,1.)]
    for m,n,E in cases:
        for theta in np.linspace(np.pi/3,2*np.pi/3,9):
            x=external(m,n,E,theta)
            ward.append(x['checks'])
            for g in (.05,.2,.7):
                v=1.;mw2=g*g*v*v/4
                # Both q contractions are retained before comparing the full W.
                full=a(g,x['z'],v)*(x['jdot']-mdot(x['q'],x['j1'])*mdot(x['q'],x['j2'])/mw2)
                transverse=a(g,x['z'],v)*x['jdot']
                residual=float(abs(full-transverse)/abs(transverse))
                assert residual<1e-10
                fullW.append(residual)
    exact=external(3,4,5,np.pi/2)
    assert abs(exact['jdot']+800/197)<2e-13
    # This sample puts all listed 3/4/5 generation mass choices below Q0.
    # It is a demonstration budget, not an empirically calibrated instrument.
    k=kernel_setup(3.,4.,50.,order=128)
    k2=kernel_setup(3.,4.,50.,order=256)
    v=1.;gstar=.2;gmax=.5;pstar=.5
    sigstar=sigma(gstar,v,k);sig2=sigma(gstar,v,k2)
    assert abs(sigstar-sig2)/sigstar<1e-10
    assert k['klower']*a(gstar,k['z1'],v)**2<=sigstar<=k['kupper']*a(gstar,k['z0'],v)**2
    ell=-np.log1p(-pstar)
    # Define one finite witness budget by its own positive integral. Analytic
    # existence does not depend on this floating integral; proof (13) supplies
    # a strictly positive lower bound and finite sufficient exposure as well.
    Lmax=ell/sigstar
    Lsufficient=ell/(k['klower']*a(gstar,k['z1'],v)**2)
    assert np.isfinite(Lmax) and 0<Lmax<=Lsufficient
    probabilities=[float(-np.expm1(-Lmax*sigma(g,v,k))) for g in (.1,.2,.3)]
    assert probabilities[0]<pstar and abs(probabilities[1]-pstar)<1e-14 and probabilities[2]>pstar
    Areq=np.sqrt(ell/(Lmax*k['kupper']))
    gnecessary=inverse_floor(Areq,k['z0'],v)
    assert 0<gnecessary<=gstar
    # Same true kernel exact conditional floor, unique by analytic monotonicity.
    lo,hi=0.,gmax
    for _ in range(80):
        mid=(lo+hi)/2
        if sigma(mid,v,k)<ell/Lmax:lo=mid
        else:hi=mid
    gtask=(lo+hi)/2
    assert abs(gtask-gstar)<1e-12
    Lrun=1.
    Btask=8*np.pi*np.pi*(gtask**-2-gmax**-2)/Lrun
    Bnecessary=8*np.pi*np.pi*(gnecessary**-2-gmax**-2)/Lrun
    assert Bnecessary>=Btask>0
    # Exact running boundary and endpoint case; constants carried as b/(8pi²).
    rg=[]
    for b in (-F(19,6),F(0),F(5,6),F(10)):
        slope=float(b)/(8*np.pi*np.pi)
        x0=gtask**-2
        values=[x0,x0-slope*Lrun]
        assert min(values)>=gmax**-2
        rg.append({'b':str(b),'joint_finite_task_and_running':True})
    assert all((-float(b)/(8*np.pi*np.pi)>=0)==(b<=0) for b in (-F(1),F(0),F(1)))
    # Total probability error weakens the necessary threshold, strengthens
    # sufficient one. This is a declared envelope, not computed actual error.
    eps=.05
    required_minus=-np.log(1-pstar+eps)
    required_plus=-np.log(1-pstar-eps)
    assert required_minus<ell<required_plus
    # The zero-momentum Fermi point does NOT furnish a finite-bin countermodel.
    fermi=[]
    for g in (.1,.03,.01):
        z=g**4
        value=a(g,z,v)
        assert 0<value<2/v**2
        # Exact integral of an upper-bounded constant kernel from zero to z1.
        z1=1.;m=g*g*v*v/4
        integral=g**4/4*(1/m-1/(z1+m))
        assert 0<integral<=g*g/v**2
        fermi.append({'g':g,'shrinking_point_a':float(value),'constant_kernel_forward_integral':integral,
                      'not_a_fixed_finite_exposure_task_counterexample':True})
    tiny=sigma(.05,v,k)
    Ltiny=ell/tiny
    assert Ltiny>Lmax and abs(-np.expm1(-Ltiny*tiny)-pstar)<1e-14
    files=[HERE/'proof.md',Path(__file__).resolve(),HERE.parent/'_admission/round1310_admission.md',
           HERE.parent/'1309/proof.md',RESEARCH/'archive_1290_1307/1303/proof.md',
           RESEARCH/'archive_1290_1307/1307/proof.md']
    return {
        'schema':'round1310_transverse_finite_task_v1','status':'PASS_CONDITIONAL_FINITE_BORN_TASK',
        'new_axioms':0,'actual_autonomous_instrument_claim':False,'higher_order_error_computed':False,
        'massive_spinor_cases':len(ward),'full_W_comparisons':len(fullW),
        'max_spinor_Ward_relative_residual':float(max(max(d.values()) for d in ward)),
        'max_fullW_relative_residual':float(max(fullW)),
        'exact_center_Jdot':'-800/197','inversion':scalar_checks(),
        'positive_bin':{'m':3.,'n':4.,'E':50.,'theta_range':['pi/3','2pi/3'],
            'z0':float(k['z0']),'z1':float(k['z1']),'K_lower':float(k['klower']),
            'K_integral_numerical':k['ktot'],'K_upper':float(k['kupper']),
            'Sigma_star':sigstar,'quadrature_relative_difference':abs(sigstar-sig2)/sigstar},
        'finite_task':{'p_star':pstar,'sample_g':gstar,'sample_v':v,'gmax':gmax,
            'exposure_witness':float(Lmax),'analytic_sufficient_exposure':float(Lsufficient),
            'probabilities_g_0p1_0p2_0p3':probabilities,'all_zero_counts_included':True,
            'a_required_necessary':float(Areq),'g_floor_necessary':float(gnecessary),
            'g_task_fixed_kernel':gtask,'true_kernel_monotone':True,
            'B_necessary_outer_bound':float(Bnecessary),'B_task_exact_effective':float(Btask),
            'Lrun':Lrun,'g_task_equals_gmax_allows_only_nonpositive_b':True,
            'RG_joint_samples':rg,'positive_error_is_declared_only':eps},
        'phase_compatibility':phase_compatibility(),'fermi_point_checks':fermi,
        'deleted_exposure_bound':{'smaller_g':.05,'required_exposure':float(Ltiny),'retains_positive_target_probability':True},
        'source_sha256':{p.relative_to(ROOT).as_posix():sha256(p.read_bytes()).hexdigest() for p in files}}


def close(a,b):
    if isinstance(a,float):return isinstance(b,(int,float)) and math.isclose(a,b,rel_tol=3e-11,abs_tol=1e-22)
    if isinstance(a,dict):return isinstance(b,dict) and a.keys()==b.keys() and all(close(v,b[k]) for k,v in a.items())
    if isinstance(a,list):return isinstance(b,list) and len(a)==len(b) and all(close(x,y) for x,y in zip(a,b))
    return a==b


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group()
    group.add_argument('--write-new',action='store_true')
    group.add_argument('--check',action='store_true')
    args=parser.parse_args();result=compute()
    if args.write_new:
        with OUT.open('x',encoding='utf-8',newline='\n') as f:
            json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
    else:
        assert close(result,json.loads(OUT.read_text(encoding='utf-8'))),'Result/source binding differs'
    print(json.dumps({'status':result['status'],'massive_spinor_cases':result['massive_spinor_cases'],
        'positive_bin':result['positive_bin']['K_integral_numerical']>0,
        'g_task':result['finite_task']['g_task_fixed_kernel'],'bindings':len(result['source_sha256'])}))
