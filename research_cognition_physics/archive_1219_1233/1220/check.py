"""1220: finite quantum source, mass-changing detector and two field modes.

Default read-only replay. --write-results writes only this round's results.json.
This is a registered effective Hamiltonian, not a local relativistic QFT.
Actual record permissions and their work account are conditional R1220 inputs.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
EPS=np.finfo(float).eps
PARAM=dict(K=14,N=6,m=25.,mass=20.,kappa=2.,g=.15)
T=5*math.pi/8
GAP=2/3
NUMERIC_RESERVE=1e-5
CONTRACT=.001
TARGET=.005
STEPS=4096


def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def exp_neg(x):
    """exp(-i*x), Taylor on |x/2**s|<=1/8 then repeated squaring.

    Degree 18 tail < exp(1/8)*(1/8)**19/19!. Complex arithmetic
    gets a deliberately loose 512*eps per reduced evaluation, then
    (2**s)*(512*eps+16*s*eps) after squaring. See proof section 6.
    """
    x=np.asarray(x,float)
    maximum=float(np.max(abs(x))) if x.size else 0.
    s=max(0,int(math.ceil(math.log2(maximum/.125)))) if maximum>.125 else 0
    z=-1j*x/(2**s); term=np.ones_like(z); out=term.copy()
    for j in range(1,19): term=term*z/j;out+=term
    for _ in range(s): out=out*out
    return out


def states(K,N,hin):
    return [(h,p,a,b) for h in (0,1) for p in range(-K,K+1)
            for a in range(N+1) for b in range(N+1-a) if (h+a+b)%2==hin]


def hops(x,kappa,g):
    h,p,a,b=x
    for dp in (-1,1): yield (h,p+dp,a,b),-kappa/2,'source'
    yield (1-h,p-1,a+1,b),g*math.sqrt(a+1),'field'
    yield (1-h,p+1,a,b+1),g*math.sqrt(b+1),'field'
    if a: yield (1-h,p+1,a-1,b),g*math.sqrt(a),'field'
    if b: yield (1-h,p-1,a,b-1),g*math.sqrt(b),'field'


def model(hin,gap=GAP,kappa=2.,g=.15,rel=True,K=14,N=6):
    bs=states(K,N,hin); ix={x:i for i,x in enumerate(bs)}; d=len(bs)
    H=np.zeros((d,d)); Vs=np.zeros_like(H); Vf=np.zeros_like(H)
    diag=np.zeros(d); nr=np.zeros(d); q=np.zeros(d); labels=[]; outside={}
    for i,x in enumerate(bs):
        h,p,a,b=x; qq=-p-a+b; q[i]=qq; mm=25.+gap*h
        exact=math.sqrt(mm*mm+p*p)-25.+math.sqrt(400.+qq*qq)-20.
        nonrel=gap*h+p*p/(2*mm)+qq*qq/40
        diag[i]=(exact if rel else nonrel)+a+b+kappa
        nr[i]=nonrel-exact
        H[i,i]=diag[i]
        labels.append((h,int(np.sign(b-a)),int(qq>0)))
        for target,amp,kind in hops(x,kappa,g):
            if amp==0: continue
            j=ix.get(target)
            if j is None:
                outside.setdefault(target,np.zeros(d))[i]+=amp
            else:
                H[j,i]+=amp
                (Vs if kind=='source' else Vf)[j,i]+=amp
            # Conservation holds on every actual creation/annihilation/force hop.
            hh,pp,aa,bb=target; qs=-pp-aa+bb
            if kind=='source': assert aa==a and bb==b and qs==qq-(pp-p)
            else: assert qs==qq and hh==1-h
    R=np.array(list(outside.values()))
    assert np.linalg.norm(H-H.T)<1e-12
    assert np.linalg.norm(R)>0
    return dict(bs=bs,ix=ix,H=H,Vs=Vs,Vf=Vf,R=R,q=q,labels=labels,nr=nr)


def diagonalize(H):
    w,V=np.linalg.eigh(H)
    n=len(w); vf=float(np.linalg.norm(V,'fro')); hf=float(np.linalg.norm(H,'fro'))
    orth=float(np.linalg.norm(V.T@V-np.eye(n),'fro'))
    eig=float(np.linalg.norm(H@V-V*w,'fro'))
    # Generous real-operation rounding bounds for verification products,
    # matrix construction, and spectral evaluation; not an empirical error.
    gamma=16*n*EPS/(1-16*n*EPS)
    # Factor two also covers rounding in Frobenius norms, gamma, and their
    # positive sums. It is not inferred from a repeated eigensolver run.
    orth_upper=2*(orth+gamma*vf*vf+1e-14)
    eig_upper=2*(eig+gamma*hf*vf+256*EPS*n*max(1.,float(np.max(abs(H))))+1e-14)
    s=max(0,int(math.ceil(math.log2(max(1.,T*float(np.max(abs(w))))/.125))))
    phase_upper=2*2**s*(512*EPS+16*s*EPS+1.14*.125**19/math.factorial(19))
    numerical_upper=50*(1+T)*(orth_upper+eig_upper+phase_upper)
    assert np.max(np.sum(abs(H),axis=1))<64 and orth_upper<1e-6
    assert numerical_upper<NUMERIC_RESERVE
    return w,V,dict(orthogonality_fro=orth,eigen_residual_fro=eig,
                    verified_product_rounding_upper=numerical_upper,
                    polynomial_phase_scalar_error_upper=phase_upper,
                    reserved_state_norm=NUMERIC_RESERVE)


def initial(m,hin):
    Z=np.zeros((len(m['bs']),3),complex)
    for j,p in enumerate((-1,0,1)): Z[m['ix'][(hin,p,0,0)],j]=1
    alpha=np.array([.5j,1/math.sqrt(2),-.5j])
    return Z,alpha


def probability_data(m,w,V,Z,alpha,t=T):
    Y=V@(exp_neg(t*w)[:,None]*(V.T@Z))
    psi=Y@alpha; coherent=abs(psi)**2
    mixed=(abs(Y)**2)@(abs(alpha)**2)
    bs=m['bs']; hin=bs[np.argmax(abs(Z[:,0]))][0]
    event=np.array([h!=hin and b>a and q>0 for (h,p,a,b),q in zip(bs,m['q'])])
    flip=np.array([h!=hin for h,p,a,b in bs]); plus=np.array([a>b for h,p,a,b in bs]); minus=np.array([b>a for h,p,a,b in bs])
    photon=np.array([a+b for h,p,a,b in bs]); pvals=np.array([p for h,p,a,b in bs])
    psi0=Z@alpha
    E0=float(np.vdot(psi0,m['H']@psi0).real)
    E0mix=float(np.dot(abs(alpha)**2,np.sum(Z.conj()*(m['H']@Z),axis=0).real))
    Ef=float(np.vdot(psi,m['H']@psi).real)
    assert abs(E0-E0mix)<1e-11 and abs(Ef-E0)<1e-10
    labels=np.array(m['labels'])
    keep=np.all(labels[:,None,:]==labels[None,:,:],axis=2)
    # The reference is 12-block pinching, not full basis dephasing.
    work_s=np.where(keep,0.,-m['Vs']); work_f=np.where(keep,0.,-m['Vf'])
    ws=float(np.vdot(psi,work_s@psi).real); wf=float(np.vdot(psi,work_f@psi).real)
    assert np.linalg.norm(work_f+m['Vf'])<1e-12
    # Every fixed h,n+,- fibre loses just the q=0 <-> 1 source edge.
    assert np.max(np.sum(abs(work_s),axis=0))<=1+1e-12
    q_negative=m['q']<0
    return dict(joint_coherent=float(coherent@event),joint_dephased=float(mixed@event),
                joint_difference=float((coherent-mixed)@event),
                flip_probability=float(coherent@flip),dephased_flip_probability=float(mixed@flip),
                flip_positive_field=float(coherent@(flip&plus)),flip_negative_field=float(coherent@(flip&minus)),
                initial_reduced_energy=E0,initial_dephased_reduced_energy=E0mix,
                energy_conservation_residual=abs(Ef-E0),
                projected_detector_p_mean=float(coherent@pvals),projected_source_q_mean=float(coherent@m['q']),
                projected_field_momentum=float(coherent@np.array([a-b for h,p,a,b in bs])),
                projected_field_number=float(coherent@photon),
                projected_reference_pinching_work_source=ws,projected_reference_pinching_work_field=wf,
                projected_reference_pinching_work_total=ws+wf,
                initial_force_detector=float((1j*(np.vdot(psi0,m['H']@(pvals*psi0))-np.vdot(psi0,pvals*(m['H']@psi0)))).real),
                initial_force_source=float((1j*(np.vdot(psi0,m['H']@(m['q']*psi0))-np.vdot(psi0,m['q']*(m['H']@psi0)))).real),
                pointer_event_definition='h!=hin AND n_minus>n_plus AND source_q>0'),Y


def residual_certificate(m,w,V,Z):
    H,R=m['H'],m['R']
    rnorm=math.sqrt(float(np.max(np.sum(abs(R),axis=0)))*float(np.max(np.sum(abs(R),axis=1))))
    C2=rnorm*float(np.linalg.norm(H@(H@Z),'fro'))
    B=R@V; C=V.T@Z; dt=T/STEPS
    midpoint_sum=derivative_sum=0.
    for start in range(0,STEPS,128):
        times=(np.arange(start,min(start+128,STEPS))+.5)*dt
        phase=exp_neg(w[:,None]*times)
        norms=np.zeros(len(times)); dnorms=np.zeros(len(times))
        for j in range(3):
            amps=phase*C[:,j,None]
            v=B@amps; dv=B@(-1j*w[:,None]*amps)
            norms+=np.sum(abs(v)**2,axis=0); dnorms+=np.sum(abs(dv)**2,axis=0)
        midpoint_sum+=float(np.sum(np.sqrt(norms)))
        derivative_sum+=float(np.sum(np.sqrt(dnorms)))
    value=dt*midpoint_sum+dt*dt*derivative_sum/4+T*dt*dt*C2/24
    state=value+NUMERIC_RESERVE
    # Pure/entangled preparations in range(Z); deliberately loose event bound.
    event=2*state+state*state
    return dict(outside_rows=R.shape[0],uniform_preparation_dimension=3,
                intervals=STEPS,R_norm_upper=rnorm,second_derivative_upper=C2,
                midpoint_integral=dt*midpoint_sum,linear_envelope=dt*dt*derivative_sum/4,
                quadratic_envelope=T*dt*dt*C2/24,
                projection_state_norm_upper=value,numeric_state_reserve=NUMERIC_RESERVE,
                total_state_norm_upper=state,event_probability_upper=event)


def finite_nr_comparison(m,w,V,Z,alpha,data):
    # Difference of the two registered finite effective kinetic Hamiltonians.
    # The norm bound is state/preparation dependent, not full Lorentz accuracy.
    D=m['nr']; C=V.T@Z; steps=512; dt=T/steps
    C2=float(np.max(abs(D)))*float(np.linalg.norm(m['H']@(m['H']@Z),'fro'))
    s0=s1=0.
    for start in range(0,steps,64):
        times=(np.arange(start,min(start+64,steps))+.5)*dt
        phase=exp_neg(w[:,None]*times); v0=np.zeros(len(times));v1=np.zeros(len(times))
        for j in range(3):
            a=phase*C[:,j,None]
            z=D[:,None]*(V@a);dz=D[:,None]*(V@(-1j*w[:,None]*a))
            v0+=np.sum(abs(z)**2,axis=0);v1+=np.sum(abs(dz)**2,axis=0)
        s0+=float(np.sum(np.sqrt(v0)));s1+=float(np.sum(np.sqrt(v1)))
    bound=dt*s0+dt*dt*s1/4+T*dt*dt*C2/24+NUMERIC_RESERVE
    hn=m['H']+np.diag(D);wn,Vn,_=diagonalize(hn)
    mn=dict(m,H=hn)
    nd,_=probability_data(mn,wn,Vn,Z,alpha)
    # Add already-certified projection residual for the full NR comparator;
    # outside interaction terms are unchanged.
    return dict(finite_projection_state_norm_upper=bound,
                finite_joint_probability=nd['joint_coherent'],
                finite_joint_difference_from_sqrt=nd['joint_coherent']-data['joint_coherent'],
                declaration='Comparison of these two effective kinetic choices only; no bound to full SR or QFT.')


def dependencies():
    names=['archive_1106_1119/1111/check.py','archive_1106_1119/1118/proof.md',
           'archive_1106_1119/1119/proof.md','archive_1159_1160/1160/proof.md',
           'archive_1161_1200/1169/proof.md','archive_1219_/1219/candidate_proof.md',
           'archive_1219_/1220/admission.md']
    return {name:sha(ROOT/name) for name in names}


def old_clock():
    p=ROOT/'archive_1106_1119/1111/check.py'
    spec=importlib.util.spec_from_file_location('old1111',p); mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    old=mod.model(); tau=math.pi/4+T; W=old['w'](tau)
    state=W.reshape(4,8,3,3,2)
    probs=np.sum(abs(state[:,[1,3,5,7],...])**2,axis=(0,1,2,3))
    gram=float(np.linalg.norm(W.conj().T@W-np.eye(2)))
    assert gram<1e-10 and abs(probs[0]-probs[1])<1e-10
    return dict(start_age=math.pi/4,elapsed=T,finish_age=tau,Mr=probs.tolist(),
                unknown_Q_isometry_residual=gram,
                limitation='Rest lab clock only, independent of accelerator; no accelerated proper-clock certification.')


def calculate():
    out=dict(round=1220,status='author candidate; independent acceptance required',
             model='translation-invariant circle, two finite-band boson modes, finite quantum source',
             parameters=dict(PARAM,gap=GAP,T=T,holdout=True,contract_each=CONTRACT,target=TARGET),
             source_sha256=dependencies(),code_sha256=sha(__file__),groups={})
    cases={}
    for hin in (0,1):
        m=model(hin);w,V,num=diagonalize(m['H']);Z,alpha=initial(m,hin)
        data,Y=probability_data(m,w,V,Z,alpha)
        cert=residual_certificate(m,w,V,Z)
        nr=finite_nr_comparison(m,w,V,Z,alpha,data)
        lower=data['joint_difference']-2*(cert['event_probability_upper']+CONTRACT)
        assert lower>TARGET
        # Matched zero-force comparison: only registered kappa changes.
        mz=model(hin,kappa=0);wz,Vz,_=diagonalize(mz['H']);zz,az=initial(mz,hin)
        zero,_=probability_data(mz,wz,Vz,zz,az)
        assert abs(zero['joint_difference'])<1e-10
        # With g=0 the internal value cannot flip, exactly, on the full model.
        assert all(target[0]==x[0] for x in m['bs'] for target,amp,kind in hops(x,2,0) if amp!=0)
        Nbound=2*(data['initial_reduced_energy']+1e-9+4*.15**2)
        work_bound=1+2*.15*math.sqrt(2*Nbound)
        assert abs(data['projected_reference_pinching_work_total'])<work_bound
        cases[str(hin)]=dict(dimension=len(m['bs']),response=data,certificate=cert,
               numerical_checks=num,nonrelativistic_comparator=nr,
               zero_force_joint=zero['joint_coherent'],zero_force_difference=zero['joint_difference'],
               source_off_probability_difference_bound='exact zero for coherent-vs-dephased on this event',
               zero_field_flip_probability=0.,certified_joint_difference_lower=lower,
               total_model_energy_above_m_plus_M=data['initial_reduced_energy']+2*.15**2,
               full_model_field_number_bound=Nbound,
               reference_12_block_pinching_mean_work_absolute_bound=work_bound)
    out['groups']['same_kernel_source_field_mass_and_joint_holdout']=cases
    out['groups']['source_and_field_controls']={h:dict(zero_force_difference=c['zero_force_difference'],
        zero_field_flip_probability=c['zero_field_flip_probability']) for h,c in cases.items()}
    out['groups']['projection_and_numerical_certificate']={h:dict(
        event_probability_upper=c['certificate']['event_probability_upper'],
        verified_rounding_upper=c['numerical_checks']['verified_product_rounding_upper']) for h,c in cases.items()}
    out['groups']['reference_record_work_account']={h:dict(
        reference_mean_work=c['response']['projected_reference_pinching_work_total'],
        full_model_absolute_mean_work_bound=c['reference_12_block_pinching_mean_work_absolute_bound']) for h,c in cases.items()}
    out['groups']['original_lab_clock_and_unknown_Q']=old_clock()
    out['check_groups']=len(out['groups'])
    assert out['check_groups']==5
    out['scope']=dict(
        field_reference='Two working-mode vacuum |0,0>, independent of source, not Minkowski vacuum.',
        initial_force='coherent detector/source = -sqrt(2), +sqrt(2); dephased = 0,0',
        conserved_quantity='P=p+q+n_plus-n_minus, exact zero on the chosen sector',
        positive_energy='H_full=H_reduced+m+M+2g^2; lower bound m+M via square completion',
        reference_record='12-block (h,sign(n_minus-n_plus),1[q>0]) pinching; finite positive-mass actual apparatus is R1220 input',
        readout_work='Neither QND nor energy free; source-cross-bin potential and field coupling work included. Mean-work bound is not total apparatus cost or deterministic work capacity.',
        preparation='Finite-support joint source/material state; dephased version retains a three-valued preparation label. Reference/coherence/record preparation is input.',
        clock='Existing rest D is retained, no moving thermometer clock supplied.',
        no_claim=['Unruh temperature','KMS reconstruction','full Lorentz-covariant interaction','causal local QFT','Einstein dynamics','manufactured instrument','full stage completion'])
    return out


def compare(a,b,path=''):
    if isinstance(a,dict):
        assert a.keys()==b.keys(),path
        for k in a: compare(a[k],b[k],path+'/'+k)
    elif isinstance(a,list):
        assert len(a)==len(b),path
        for i,(x,y) in enumerate(zip(a,b)): compare(x,y,path+str(i))
    elif isinstance(a,float): assert math.isclose(a,b,rel_tol=2e-7,abs_tol=1e-10),(path,a,b)
    else: assert a==b,(path,a,b)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true');args=parser.parse_args()
    result=calculate();path=HERE/'results.json'
    if args.write_results:
        if path.exists(): raise FileExistsError('Refusing to overwrite existing results.json; default replay is read-only.')
        path.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:compare(result,json.loads(path.read_text(encoding='utf-8')))
    print(json.dumps(dict(status='PASS',round=1220,read_only=not args.write_results,
        check_groups=result['check_groups'],checks=list(result['groups']),
        certified_lower={h:c['certified_joint_difference_lower'] for h,c in result['groups']['same_kernel_source_field_mass_and_joint_holdout'].items()}),ensure_ascii=False))
