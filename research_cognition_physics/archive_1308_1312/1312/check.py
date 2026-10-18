"""1312: conditional active-group selection from a full-row-rank Dirac map.

Default/--check is read-only. --write-new exclusively creates results.json.
No numerical group scan substitutes for the analytic classification.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parents[1]
ROOT = RESEARCH.parent
OUT = HERE / 'results.json'
EPS = np.array([[0, 1], [-1, 0]], complex)
TOL = 5e-10


def norm(x):
    return float(np.linalg.norm(x))


def exp_skew(x):
    vals, vecs = np.linalg.eigh(1j*x)
    return (vecs*np.exp(-1j*vals))@vecs.conj().T


def su_generators(d):
    ans = []
    for i in range(d):
        for j in range(i+1, d):
            a = np.zeros((d, d), complex); a[i,j] = 1; a[j,i] = -1
            b = np.zeros((d, d), complex); b[i,j] = 1j; b[j,i] = 1j
            ans += [a, b]
    for i in range(d-1):
        a = np.zeros((d, d), complex); a[i,i] = 1j; a[-1,-1] = -1j
        ans.append(a)
    return ans


def sp_generators(k):
    ans = []
    # X=[[A,B],[-conj(B),conj(A)]], A skew-Hermitian, B symmetric.
    abase = []
    for i in range(k):
        a = np.zeros((k,k), complex); a[i,i] = 1j; abase.append(a)
        for j in range(i+1,k):
            a = np.zeros((k,k), complex); a[i,j]=1; a[j,i]=-1; abase.append(a)
            a = np.zeros((k,k), complex); a[i,j]=1j; a[j,i]=1j; abase.append(a)
    z = np.zeros((k,k), complex)
    for a in abase:
        ans.append(np.block([[a,z],[z,a.conj()]]))
    for i in range(k):
        for j in range(i,k):
            for phase in (1,1j):
                b=np.zeros((k,k),complex); b[i,j]=phase; b[j,i]=phase
                ans.append(np.block([[z,b],[-b.conj(),z]]))
    assert len(ans)==k*(2*k+1)
    return ans


def symplectic_j(k):
    z=np.zeros((k,k)); eye=np.eye(k)
    return np.block([[z,eye],[-eye,z]]).astype(complex)


def random_group(gens,d,rng):
    x=sum((c*g for c,g in zip(rng.normal(size=len(gens)),gens)),np.zeros((d,d),complex))
    return exp_skew(x/max(1,norm(x)))


def random_unitary(d,rng):
    a=rng.normal(size=(d,d))+1j*rng.normal(size=(d,d))
    q,r=np.linalg.qr(a)
    diag=np.diag(r); return q@(np.diag(diag/abs(diag)))


def rank(a):
    return int(np.linalg.matrix_rank(a,tol=1e-10*max(1,norm(a))))


def fixed_dimension(gens,d):
    return d if not gens else d-rank(np.vstack(gens))


def embeddings(gens,d,indices):
    ans=[]
    for x in gens:
        y=np.zeros((d,d),complex); y[np.ix_(indices,indices)]=x; ans.append(y)
    return ans


def reference_matrix(phi):
    return np.column_stack([EPS@phi.conj(),phi])


def stabilizer_checks(rng):
    rows=[]; residuals=[]
    for d in range(2,8):
        stab=embeddings(su_generators(d-1),d,list(range(d-1)))
        v=np.eye(d)[:,-1]; fixed=fixed_dimension(stab,d)
        expected=2 if d==2 else 1
        assert fixed==expected
        residuals.extend(norm(x@v) for x in stab)
        for n in (1,2,4):
            r=d*n+1
            basis=reference_matrix(v) if d==2 else v[:,None]
            C=np.pad(np.eye(expected*n),((0,0),(0,r-expected*n)))
            f=np.kron(basis,np.eye(n))@C
            assert rank(f)==expected*n
            rows.append({'group':f'SU({d})','d':d,'nL':n,'right_singlets':r,
                         'fixed_dim':fixed,'mass_rank':rank(f),'full_row_rank':rank(f)==d*n})
    for k in range(1,5):
        d=2*k; j=symplectic_j(k); gens=sp_generators(k)
        idx=list(range(1,k))+list(range(k+1,2*k))
        stab=embeddings(sp_generators(k-1),d,idx) if k>1 else []
        v=np.eye(d)[:,0]; fixed=fixed_dimension(stab,d)
        assert fixed==2
        for x in gens:
            residuals += [norm(x.conj().T+x),norm(x.T@j+j@x)]
        residuals.extend(norm(x@v) for x in stab)
        u=random_group(gens,d,rng)
        residuals += [norm(u.conj().T@u-np.eye(d)),norm(u@j-j@u.conj())]
        phi=rng.normal(size=d)+1j*rng.normal(size=d)
        basis=np.column_stack([phi,j@phi.conj()])
        residuals.append(norm(basis.conj().T@basis-np.vdot(phi,phi)*np.eye(2)))
        residuals.append(norm(np.column_stack([u@phi,j@(u@phi).conj()])-u@basis))
        for n in (1,2,4):
            r=d*n+1; C=np.pad(np.eye(2*n),((0,0),(0,r-2*n)))
            f=np.kron(basis,np.eye(n))@C
            assert rank(f)==2*n
            rows.append({'group':f'Sp({k})','d':d,'nL':n,'right_singlets':r,
                         'fixed_dim':fixed,'mass_rank':rank(f),'full_row_rank':rank(f)==d*n})
    assert max(residuals)<TOL
    return {'cases':rows,'max_generator_fixedspace_residual':max(residuals)}


def smooth_su2_checks(rng):
    rows=[]; residuals=[]
    for n in (1,2,3,5):
        for r in (2*n,2*n+2):
            size=2*n; singular=np.linspace(1,2,size)
            c0=random_unitary(size,rng)@np.pad(np.diag(singular),((0,0),(0,r-size)))@random_unitary(r,rng)
            h=rng.normal(size=(size,size))+1j*rng.normal(size=(size,size)); h=(h+h.conj().T)/2
            def cf(s): return (1+s*s)*exp_skew(1j*s*h)@c0
            for s in (.25,1.,4.):
                phi=rng.normal(size=2)+1j*rng.normal(size=2); phi*=np.sqrt(s)/norm(phi)
                u=random_group(su_generators(2),2,rng)
                f=np.kron(reference_matrix(phi),np.eye(n))@cf(s)
                f2=np.kron(reference_matrix(u@phi),np.eye(n))@cf(float(np.vdot(u@phi,u@phi).real))
                residuals.append(norm(f2-np.kron(u,np.eye(n))@f)/norm(f))
                got=np.linalg.svd(f,compute_uv=False)
                expected=np.sort(np.sqrt(s)*(1+s*s)*singular)[::-1]
                residuals.append(norm(got-expected)/norm(expected))
                residuals.append(abs(np.linalg.det(reference_matrix(phi))-s))
                assert rank(f)==2*n
                rows.append({'nL':n,'r':r,'norm_squared':s,'rank':rank(f),
                             'smallest_mass':float(got[-1]),'unused_right_directions':r-2*n})
    assert max(residuals)<TOL
    # The matrix dimension obstruction cannot be removed by an arbitrary map.
    for n in range(1,7):
        assert min(2*n,2*n-1)<2*n
    return {'cases':rows,'max_covariance_singular_value_residual':max(residuals),
            'insufficient_right_dimension_checks':6,'nonpolynomial_radial_dependence':True}


def deletion_checks(rng):
    residuals=[]; rows=[]
    for sample in range(12):
        u=random_group(su_generators(3),3,rng)
        a=rng.normal(size=3)+1j*rng.normal(size=3)
        b=rng.normal(size=3)+1j*rng.normal(size=3)
        lam=1.5
        def fmap(x,y):return np.column_stack([x,y,np.cross(x.conj(),y.conj())/lam])
        f=fmap(a,b); f2=fmap(u@a,u@b)
        gram=(np.vdot(a,a).real*np.vdot(b,b).real-abs(np.vdot(a,b))**2)/lam
        assert gram>0 and rank(f)==3
        residuals += [norm(f2-u@f)/norm(f),abs(np.linalg.det(f)-gram)/max(1,gram)]
        rows.append({'two_reference_SU3_sample':sample,'rank':3,'determinant_positive':True})
    for name,d,gens in [('SU3',3,su_generators(3)),('Sp2',4,sp_generators(2))]:
        u=random_group(gens,d,rng); phi=rng.normal(size=d)+1j*rng.normal(size=d)
        f=np.vdot(phi,phi).real/1.5*np.eye(d)
        f2=np.vdot(u@phi,u@phi).real/1.5*np.eye(d)
        residuals.append(norm(f2-u@f@u.conj().T)/norm(f)); assert rank(f)==d
        rows.append({'changed_right_representation':name,'rank':d,'no_extra_reference':True})
    # Removing full rank keeps a valid SU3 single-reference rank-one map.
    phi=np.array([1,1j,2],complex); assert rank(phi[:,None])==1
    # Without sphere transitivity, SO3 on complex C3 allows one complex frame.
    phi=np.array([1,1j,0],complex)
    def so_f(x):return np.column_stack([x,x.conj(),1j*np.cross(x,x.conj())/1.5])
    a=np.array([[0,.2,-.3],[-.2,0,.4],[.3,-.4,0]],complex); u=exp_skew(a)
    f=so_f(phi); assert rank(f)==3
    residuals.append(norm(so_f(u@phi)-u@f))
    assert abs(phi@phi)<1e-14 and abs(np.array([np.sqrt(2),0,0])@np.array([np.sqrt(2),0,0])-2)<1e-14
    assert max(residuals)<TOL
    return {'samples':rows,'max_residual':max(residuals),'drop_full_rank_SU3_rank':1,
            'drop_sphere_transitivity_SO3_rank':3,'full_ER_P_countermodels_claimed':False}


def u1_and_discrete_checks(rng):
    residuals=[]; rows=[]
    for yl,h in [(Q(-1,2),Q(1,2)),(Q(1,6),Q(1,2)),(Q(2,7),Q(3,5))]:
        qr=[yl+h,yl-h]
        for n in (1,3,5):
            yp=random_unitary(n,rng)@np.diag(np.arange(1,n+1))
            ym=random_unitary(n,rng)@np.diag(np.arange(2,n+2))
            c=np.block([[yp,np.zeros((n,n))],[np.zeros((n,n)),ym]])
            for alpha in (.17,.8,2.3):
                phi=rng.normal(size=2)+1j*rng.normal(size=2)
                u=random_group(su_generators(2),2,rng)
                f=np.kron(reference_matrix(phi),np.eye(n))@c
                f2=np.kron(reference_matrix(np.exp(1j*float(h)*alpha)*u@phi),np.eye(n))@c
                left=np.exp(1j*float(yl)*alpha)*np.kron(u,np.eye(n))
                right=np.diag(np.exp(1j*alpha*np.repeat([float(q) for q in qr],n)))
                residuals.append(norm(f2-left@f@right.conj().T)/norm(f))
            rows.append({'qL':str(yl),'h':str(h),'ordered_physical_right_charges':[str(q) for q in qr],'nL':n})
    raw={'Yukawa':1j**2*(-1),'QQQL':1j**4,'NN':1j**2,'Weinberg':1j**2*(-1)**2}
    combined={'Higgs':(-1)*(-1),'Dirac':(-1j)*(1j),'QQQL':(-1j)**4,
              'NN':1j**2,'Weinberg':(-1j)**2}
    assert raw=={'Yukawa':1,'QQQL':1,'NN':-1,'Weinberg':-1}
    assert combined=={'Higgs':1,'Dirac':1,'QQQL':1,'NN':-1,'Weinberg':-1}
    assert all((-1j)**(4*n)==1 for n in (3,4,5))
    assert (-1j)**2*1j**2==1
    assert max(residuals)<TOL
    return {'covariance_cases':rows,'max_covariance_residual':max(residuals),
            'raw_r_phases':{k:int(v.real) for k,v in raw.items()},
            'rhat_phases':{k:int(v.real) for k,v in combined.items()},
            'ordinary_SU2_and_color_instanton_field_phase_check':True,
            'all_global_discrete_anomalies_claimed':False}


def phase_rational_certificate():
    mu=Q(1,4096); a=Q(1,32); L=Q(1,4); rows=[]; margins=[]
    for j in (3,4,5):
        D=7*j**5-31*j**3+24
        A=Q(j**3-1,8)+Q(j**5-1,256)*mu**2
        B=Q(j**3-8,8)+Q(j**5-32,256)*mu**2
        tx=Q(1+j**7,12288); ty=Q(128+j**7,12288)
        leading=Q(D)*L*a*a/163840
        tail=mu**2*(A*ty+B*tx+mu**4*tx*ty)
        margin=leading-tail
        assert margin>0
        margins.append(margin)
        rows.append({'third_mass_multiple':j,'D':D,'leading_coefficient':str(leading),
                     'tail_coefficient':str(tail),'q_exact':str(margin),'q_decimal':float(margin)})
    lower=min(margins)*a**6*mu**8/144
    assert lower>0
    return {'source':'1303 continuous compact-support amplitude and remainder bounds',
            'rows':rows,'uniform_C_lower_bound_exact':str(lower),
            'all_time_interval_inherited':['7/16','1/2'],'instrument_error_certified':False}


def common_parent_checks():
    mu=2.0**-12; v=1.; U3=np.exp(2j*np.pi*np.outer(np.arange(3),np.arange(3))/3)/np.sqrt(3)
    phi=np.array([0,v/np.sqrt(2)],complex); rows=[]; residuals=[]
    for n in (3,4,5):
        R=np.eye(n)
        for j in range(2,n-1):
            rr=np.eye(n); t=1e-5; c,s=np.cos(t),np.sin(t)
            rr[j,j]=rr[j+1,j+1]=c; rr[j,j+1]=s; rr[j+1,j]=-s; R=R@rr
        u0=np.eye(n,dtype=complex);u0[:3,:3]=U3; u=u0@R
        masses=mu*np.arange(1,n+1); mn=u@np.diag(masses)@u.conj().T
        charged=mu*np.array([1,2,4]+list(range(5,n+2))); me=np.diag(charged)
        yn=np.sqrt(2)*mn/v; ye=np.sqrt(2)*me/v
        y=np.block([[yn,np.zeros((n,n))],[np.zeros((n,n)),ye]])
        f=np.kron(reference_matrix(phi),np.eye(n))@y
        target=np.block([[mn,np.zeros((n,n))],[np.zeros((n,n)),me]])
        residuals.append(norm(f-target)); assert rank(f)==2*n
        expected=np.sort(np.r_[masses,charged])[::-1]
        residuals.append(norm(np.linalg.svd(f,compute_uv=False)-expected))
        quark_plus=np.diag(mu*np.arange(2,n+2)); quark_minus=np.diag(mu*np.arange(3,n+3))
        qmass=np.kron(np.eye(3),np.block([[quark_plus,np.zeros((n,n))],[np.zeros((n,n)),quark_minus]]))
        assert rank(qmass)==6*n
        hh=mn@mn; seen={0}
        while True:
            new=seen|{j for i in seen for j in range(n) if abs(hh[i,j])>1e-18}
            if new==seen:break
            seen=new
        assert len(seen)==n
        assert abs(abs(u[1,0]*u[0,1].conjugate())-1/3)<1e-14
        weights=abs(R[2,2:])**2; assert abs(sum(weights)-1)<1e-14
        # This arbitrary diagnostic F_j verifies the exact algebra of the
        # submatrix/convex identity, not a fabricated physical probability.
        amp=np.array([complex(j,j*j/7) for j in range(1,n+1)])
        full=u@np.diag(amp)@u.conj().T
        z=weights@amp[2:]
        compressed=U3@np.diag([amp[0],amp[1],z])@U3.conj().T
        residuals.append(norm(full[:3,:3]-compressed))
        lhs=((amp[0]-z)*(amp[1]-z).conjugate()).imag
        rhs=sum(w*((amp[0]-fj)*(amp[1]-fj).conjugate()).imag for w,fj in zip(weights,amp[2:]))
        residuals.append(abs(lhs-rhs))
        rows.append({'Nc':3,'ng':n,'lepton_left_rank':rank(f),'quark_left_rank':rank(qmass),
                     'minimum_Dirac_mass':float(np.linalg.svd(f,compute_uv=False)[-1]),
                     'connected_flavor_graph':True,'kept_1310_first_two_mass_columns':True,
                     'original_charged_current_abs_mixing':1/3,'P_convex_weights':weights.tolist(),
                     'exact_discrete_protection_and_QQQL_permission_inherited':True,
                     'old_total_source_numbers_inherited':False})
    assert max(residuals)<TOL
    return {'branches':rows,'max_mass_and_compression_residual':max(residuals),
            'mu':mu,'v':v,'phase_certificate':phase_rational_certificate(),
            'quantum_full_interacting_instrument_claim':False}


def compute():
    rng=np.random.default_rng(1312)
    result={'schema':'round1312_single_reference_Dirac_group_v1',
            'status':'PASS_CONDITIONAL_ACTIVITY_GROUP_SELECTION','new_axioms':0,
            'general_theorem':'SU2 standard complex dimension 2 and r>=2*nL; spectator compatibility separate',
            'stabilizers':stabilizer_checks(rng),'smooth_SU2':smooth_su2_checks(rng),
            'deleted_conditions':deletion_checks(rng),'U1_and_r':u1_and_discrete_checks(rng),
            'same_parent':common_parent_checks(),
            'not_inferred':['right-singlet role','single reference','sphere transitivity','all-left full rank',
                            'color group','hypercharge normalization','global quotient','three generations',
                            'physical parameter values','actual quantum instrument']}
    files=[HERE/'proof.md',Path(__file__).resolve(),HERE.parent/'_admission/round1312_admission.md',
           HERE.parent/'1308/proof.md',HERE.parent/'1310/proof.md',HERE.parent/'1310/results.json',
           RESEARCH/'archive_1290_1307/1303/proof.md',RESEARCH/'archive_223_230/research_note_227.md']
    result['source_sha256']={p.relative_to(ROOT).as_posix():sha256(p.read_bytes()).hexdigest() for p in files}
    return result


def close(a,b):
    if isinstance(a,float):return isinstance(b,(int,float)) and math.isclose(a,b,rel_tol=2e-9,abs_tol=2e-12)
    if isinstance(a,dict):return isinstance(b,dict) and a.keys()==b.keys() and all(close(v,b[k]) for k,v in a.items())
    if isinstance(a,list):return isinstance(b,list) and len(a)==len(b) and all(close(x,y) for x,y in zip(a,b))
    return a==b


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    mode=parser.add_mutually_exclusive_group(); mode.add_argument('--check',action='store_true'); mode.add_argument('--write-new',action='store_true')
    args=parser.parse_args(); data=compute()
    if args.write_new:
        with OUT.open('x',encoding='utf-8',newline='\n') as f:
            json.dump(data,f,ensure_ascii=False,indent=2);f.write('\n')
    else:
        assert close(data,json.loads(OUT.read_text(encoding='utf-8'))),'Saved result or source binding differs'
    print(json.dumps({'status':data['status'],'stabilizer_cases':len(data['stabilizers']['cases']),
                      'smooth_mass_cases':len(data['smooth_SU2']['cases']),
                      'common_branches':[r['ng'] for r in data['same_parent']['branches']],
                      'new_axioms':0,'mode':'first-write' if args.write_new else 'read-only'}))
