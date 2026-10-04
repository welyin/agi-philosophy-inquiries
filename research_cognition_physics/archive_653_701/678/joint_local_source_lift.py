"""678: full local constrained Grassmann lift of677's original source functional.
Auxiliary variables and regulator determinant are representation data.
No reflection positivity or identification with original physical HF is claimed.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_rational_physical_limit as prior
base=prior.base
HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_local_source_lift_results.json'


def err(a):
    return float(np.max(np.abs(a),initial=0))


def logpf(matrix):
    """615 pivot convention, logarithmic magnitude to retain huge bulk factors."""
    a=np.array(matrix,dtype=complex,copy=True);phase=1+0j;log=0.
    assert len(a)%2==0 and err(a+a.T)<2e-11
    for k in range(0,len(a),2):
        j=k+1+int(np.argmax(abs(a[k,k+1:])))
        assert abs(a[k,j])>1e-27,('unresolved pivot',k)
        if j!=k+1:
            a[[k+1,j],:]=a[[j,k+1],:]
            a[:,[k+1,j]]=a[:,[j,k+1]]
            phase=-phase
        pivot=a[k,k+1];mag=abs(pivot)
        phase*=pivot/mag;log+=np.log(mag)
        x=a[k,k+2:].copy();y=a[k+1,k+2:].copy()
        a[k+2:,k+2:]+=(np.outer(y,x)-np.outer(x,y))/pivot
    return phase,float(log)


def blocks(x,g5,a,layers):
    n=len(x);eye=np.eye(n);b=2*eye+a*x
    ap=b+a*g5@x;am=b-a*g5@x
    m=(layers+1)*n
    k=np.zeros((m,m),complex);k[:n,:n]=b;k[:n,-n:]=b
    for s in range(1,layers+1):
        k[s*n:(s+1)*n,(s-1)*n:s*n]=-am
        k[s*n:(s+1)*n,s*n:(s+1)*n]=ap
    inject=np.zeros((m,n),complex);inject[:n]=eye
    out=np.zeros((n,m),complex);out[:,:n]=b;out[:,-n:]=-b
    return k,inject,out


def lift(x,mat,a=.23,layers=1,lam=.37):
    jm,jp,mass,bar,pair=mat;n=len(x);r=n//2
    qp=jp@jp.conj().T;qm=jm@jm.conj().T;g5=qp-qm
    k,inject,out=blocks(x,g5,a,layers);m=len(k);total=2*n+2*m
    psi=slice(0,n);chi=slice(n,2*n);z=slice(2*n,2*n+m);eta=slice(2*n+m,total)
    y=np.zeros((n,total),complex)
    y[:r,psi]=.5*jm.conj().T;y[:r,z]=.5*jm.conj().T@out
    y[r:,psi]=-.75*jp.T@mass;y[r:,chi]=jp.T;y[r:,z]=.25*jp.T@mass@out
    nc=np.zeros((total,total),complex)
    nc[psi,psi]=mass@qp;nc[chi,chi]=bar@qm
    nc[chi,psi]=.5*np.eye(n);nc[psi,chi]=-.5*np.eye(n)
    nc[chi,z]=.5*g5@out;nc[z,chi]=-nc[chi,z].T
    nc[eta,z]=k;nc[z,eta]=-k.T
    nc[eta,psi]=-inject;nc[psi,eta]=inject.T
    full=nc+lam*y.T@pair@y
    f=np.linalg.solve(k,inject)
    pull=np.zeros((total,2*n),complex);pull[:2*n]=np.eye(2*n);pull[z,:n]=f
    eps=out@f
    reference=prior.soft(eps,mat,lam)
    errors=dict(skew=err(full+full.T),
        full_action_pullback=err(pull.T@full@pull-reference['N']),
        physical_source_pullback=err(y@pull-reference['Phi']))
    assert max(errors.values())<2e-12
    return dict(N=full,Y=y,K=k,out=out,F=f,eps=eps,reference=reference,errors=errors,
        n=n,m=m,a=a,L=layers)


def augmented(n,y,z):
    c=y.T@z;k=z.shape[1]
    return np.block([[n,c],[-c.T,np.zeros((k,k))]])


def full_source_check():
    links,e,phis=prior.fixture();mat=base.fixed_matrices(e,phis)
    _,_,_,h,_=base.kernel(links)
    g5=mat[1]@mat[1].conj().T-mat[0]@mat[0].conj().T
    d=lift(g5@h,mat)
    eps,_,_=prior.regulate(h,g5,d['a'],d['L'])
    assert err(d['eps']-eps)<1e-13
    n=d['n'];m=d['m'];q=d['reference']
    sk,lk=np.linalg.slogdet(d['K'])
    nu=(-1)**(m*(m+1)//2);assert nu==1
    rng=np.random.default_rng(67810)
    z=rng.normal(size=(n,4))+1j*rng.normal(size=(n,4))
    z/=np.linalg.norm(z,axis=0)
    rows=[]
    for count in (0,2,4):
        zz=z[:,:count]
        large=augmented(d['N'],d['Y'],zz)
        small=augmented(q['N'],q['Phi'],zz)
        p,l=logpf(large)
        target=base.pf(small)
        ps,ls=logpf(small)
        assert abs(ps*np.exp(ls)/target-1)<2e-12
        ratio=(p/(nu*sk*target/abs(target)))*np.exp(l-lk-np.log(abs(target)))
        assert abs(ratio-1)<4e-10
        rows.append(dict(source_count=count,full_dimension=len(large),
            full_log_abs_Pf=l,bulk_log_abs_det=float(lk),
            original_coefficient=base.old.cpair(target),
            corrected_bulk_coefficient=base.old.cpair(p/(nu*sk)*np.exp(l-lk)),
            phase_preserving_relative_error=float(abs(ratio-1))))
    # Independent auxiliary-block inverse check, including nonzero mass zz block.
    zblock=d['N'][2*n:2*n+m,2*n:2*n+m]
    ki=np.linalg.inv(d['K'])
    qi=np.block([[np.zeros((m,m)),ki],[-ki.T,ki.T@zblock@ki]])
    qq=d['N'][2*n:,2*n:];mix=d['N'][:2*n,2*n:]
    effective=d['N'][:2*n,:2*n]+mix@qi@mix.T
    effy=d['Y'][:,:2*n]+d['Y'][:,2*n:]@qi@mix.T
    errors=dict(inverse=err(qq@qi-np.eye(2*m)),
        Schur_action=err(effective-q['N']),Schur_source=err(effy-q['Phi']),
        no_source_contact=err(d['Y'][:,2*n:]@qi@d['Y'][:,2*n:].T))
    assert max(errors.values())<2e-11
    return dict(original_full_nonflat_16_channels=True,a=d['a'],L=d['L'],
        n=n,bulk_dimension=m,original_mass_and_four_distinct_phi=True,
        direct_identity_errors=d['errors'],auxiliary_Schur_errors=errors,
        bare_berezin_sign=nu,rows=rows,no_division_by_original_physical_weight_in_theorem=True)


def normalization_and_source_check():
    links,e,phis=prior.fixture();mat=base.fixed_matrices(e,phis)
    g5=mat[1]@mat[1].conj().T-mat[0]@mat[0].conj().T
    _,_,_,h,_=base.kernel(links)
    a=.23;layers=3
    k,_,_=blocks(g5@h,g5,a,layers)
    sk,lk=np.linalg.slogdet(k)
    gram=k.conj().T@k
    sg,lg=np.linalg.slogdet(gram)
    mineig=float(np.linalg.eigvalsh(gram)[0])
    assert mineig>0 and abs(sg-1)<2e-12
    cancellation=sk*sk.conjugate()/sg*np.exp(2*lk-lg)
    assert abs(cancellation-1)<3e-10
    step=1e-4;rows=[];ks=[]
    for angle in (-step,step):
        moved=links.copy()
        moved[1,0]=base.prior.rep(np.eye(3),np.eye(2),np.exp(1j*angle))@moved[1,0]
        _,_,_,ht,_=base.kernel(moved)
        kt,_,_=blocks(g5@ht,g5,a,layers);ks.append(kt)
        st,lt=np.linalg.slogdet(kt)
        sgt,lgt=np.linalg.slogdet(kt.conj().T@kt)
        rows.append(dict(angle=angle,numerator_aux_logdet=float(lt),
            compensating_fermion_logdet=float(lt),boson_log_integral=float(-lgt),
            corrected_log_factor=float(2*lt-lgt)))
    dk=(ks[1]-ks[0])/(2*step)
    trace=np.trace(np.linalg.solve(k,dk))
    fd=(rows[1]['numerator_aux_logdet']-rows[0]['numerator_aux_logdet'])/(2*step)
    corrected=(rows[1]['corrected_log_factor']-rows[0]['corrected_log_factor'])/(2*step)
    assert abs(trace.real-fd)<2e-5 and abs(trace.imag)<1e-10
    assert abs(corrected)<1e-6 and abs(fd)>1
    # Actual676 obstruction flux: full local action/source identity still works.
    flux=np.tile(np.eye(16,dtype=complex),(2,4,1,1))
    for i,(t,x) in enumerate(base.prior.SITES):
        for mu,z in enumerate(((-1.)**t if x==1 else 1.,1j if x==1 else 1.)):
            flux[mu,i]=base.prior.rep(np.eye(3),np.eye(2),z)
    _,_,_,hf,_=base.kernel(flux);df=lift(g5@hf,mat,layers=1)
    return dict(original_676_flux_lift_errors=df['errors'],
        L=layers,bulk_dimension=len(k),boson_gram_minimum_eigenvalue=mineig,
        gaussian_logdet_identity_error=float(abs(lg-2*lk)),
        complex_phase_cancellation_error=float(abs(cancellation-1)),
        original_link_source_rows=rows,bulk_source_trace=base.old.cpair(trace),
        bulk_source_finite_difference=float(fd),corrected_total_source=float(corrected),
        local_boson_uses_K_adjoint_K=True,compensating_fermion_uses_K_adjoint=True,
        positive_boson_gaussian_is_not_reflection_positivity=True,
        finite_source_derivative_only=True)


def run():
    deps=('joint_rational_physical_limit.py','research_note_659.md','research_note_673.md',
          'research_note_676.md','research_note_677.md','round678_drafts/local_block_probe_results.json')
    return dict(date='2026-10-02',round=678,tests_run=2,failures=0,errors=0,
        complete_local_physical_source_identity=full_source_check(),
        normalization_and_original_gauge_source=normalization_and_source_check(),
        dependency_hashes={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in deps},
        scope='Exact finite original-candidate Grassmann lift for all physical sources, original mass and signed Pfaffian; original local Wilson recurrence, explicit determinant and local convergent boson/fermion compensation. After auxiliary integration677 full-average limit applies. Auxiliary fields are not new particles or cognitive units. No finite regulator reflection positivity, normalized original physical state, HF/time identity, physical continuum or quantum GR claim.',
        all_checks_passed=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true')
    args=p.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(round=678,tests_run=2,all_checks_passed=True,
        physical=result['complete_local_physical_source_identity'],
        normalization=result['normalization_and_original_gauge_source'])))

