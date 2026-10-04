"""659: original free auxiliary Pfaffian from the specified local domain wall.

Finite-volume boundary mapping, signed Grassmann weights and retained bulk
normalization. The positive transfer here is in the auxiliary fifth direction,
not a physical-time transfer or a reconstruction of the original 32CAR H.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_spatial_auxiliary_geometry as old

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_domain_wall_boundary_mapping_results.json'


def op(a):return float(np.linalg.norm(a,2))


def data(kappa=1.,nx=3,nt=2,a=.1):
    f=old.frame(kappa,nx,nt);n=len(f['H']);eye=np.eye(n)
    g=np.kron(np.eye(len(f['sites'])),old.old.spin.G5)
    x=g@f['H'];pp=(eye+g)/2;pm=eye-pp
    up=np.kron(np.eye(len(f['sites'])),old.old.VP)
    um=np.kron(np.eye(len(f['sites'])),old.old.VM)
    change=np.column_stack((up,um));r=n//2
    xc=change.conj().T@x@change
    w=np.eye(r)+a*xc[:r,:r];z=a*xc[:r,r:]
    wb=np.eye(r)+a*xc[r:,r:]
    wi=np.linalg.inv(w)
    tc=np.block([[wi,-wi@z],[-z.conj().T@wi,wb+z.conj().T@wi@z]])
    mc=np.block([[wi,-wi@z],[np.zeros_like(wi),np.eye(r)]])
    transfer=change@tc@change.conj().T
    boundary=change@mc@change.conj().T
    ht=g@(a*x)@np.linalg.inv(2*eye+a*x)
    ev,v=np.linalg.eigh(ht);minus=v[:,ev<0];plus=v[:,ev>0]
    assert len(minus.T)==r and min(np.linalg.eigvalsh(w))>0
    assert op(transfer-transfer.conj().T)<1e-12
    return dict(f=f,n=n,eye=eye,g=g,X=x,A=eye+a*x,pp=pp,pm=pm,
        um=um,T=transfer,M=boundary,HT=ht,ev=ev,minus=minus,plus=plus,a=a)


def schur_frame(d,length):
    s=d['A'].copy()
    for _ in range(1,length):s=d['A']-d['pp']@np.linalg.solve(s,d['pm'])
    _,sing,vh=np.linalg.svd(d['pp']@s)
    r=vh.conj().T[:,d['n']//2:]
    assert max(sing[d['n']//2:])<2e-12
    return r


def modified_bulk(d,length,ap=False):
    n=d['n'];s=np.diag(np.ones(length-1),1)
    if ap:s[-1,0]=-1
    q=np.kron(np.eye(length),d['A'])-np.kron(s,d['pm'])-np.kron(s.T,d['pp'])
    if not ap:q[-n:,-n:]=d['pp']@d['A']
    return q


def full_pfaffian(d,length,e):
    """Original 16 internal modes, original B/T; bar-E fixed at e0.

    Fixing bar-E is harmless for E ratios: it saturates only the decoupled
    left-null rows. Its independent sphere average is an E-independent factor.
    """
    n=d['n'];dim=16*n*length
    m=np.zeros((dim,dim),complex);bar=np.zeros_like(m)
    localpm=(np.eye(4)-old.old.spin.G5)/2
    for i,ei in enumerate(e):
        sl=slice(16*(n*(length-1)+4*i),16*(n*(length-1)+4*i+4))
        m[sl,sl]=-np.kron(old.old.B,sum(ea*t for ea,t in zip(ei,old.old.T)))
        bar[sl,sl]=-np.kron(localpm@old.old.B,old.old.T[0].conj().T)
    bulk=np.kron(modified_bulk(d,length),np.eye(16))
    return old.old.pfaffian(np.block([[m,-bulk.T/2],[bulk/2,bar]]))


def ratio(d,r,e):return old.ratio(dict(d['f'],V=r),e)


def fields(count,seed=65901):
    rng=np.random.default_rng(seed)
    e=rng.normal(size=(count,10))*.08;e[:,0]=1
    return e/np.linalg.norm(e,axis=1)[:,None]


def boundary_identity_check():
    rows=[]
    for nx,nt,a,length in ((2,2,.2,2),(3,2,.1,4),(2,4,.15,3)):
        d=data(nx=nx,nt=nt,a=a);i=d['eye'];x=d['X'];t=d['T'];m=d['M']
        errors=[op((2*i+a*x)@m-i-t),op(d['g']@(a*x)@m-i+t),
            op(d['HT']-(i-t)@np.linalg.inv(i+t))]
        r=schur_frame(d,length);pr=r@r.conj().T
        advanced=m@np.linalg.matrix_power(t,length-1)@d['um']
        qa,_=np.linalg.qr(advanced);errors.append(op(pr-qa@qa.conj().T))
        bulk=modified_bulk(d,length);_,s,vh=np.linalg.svd(bulk)
        null=vh.conj().T[:,s<1e-11];qn,_=np.linalg.qr(null[-d['n']:])
        errors.extend([op(pr-qn@qn.conj().T),op(bulk@null),
            op(d['f']['B'].conj().T@pr.conj()@d['f']['B']-pr)])
        assert max(errors)<2e-12
        # Check the infinite-depth subspace and the indispensable boundary map.
        qinf,_=np.linalg.qr(np.linalg.solve(2*i+a*x,d['minus']))
        long=schur_frame(d,512);plong=long@long.conj().T
        mapped_error=op(plong-qinf@qinf.conj().T)
        naive_error=op(plong-d['minus']@d['minus'].conj().T)
        assert mapped_error<2e-12 and naive_error>1e-4
        rows.append(dict(nx=nx,nt=nt,a5=a,L5=length,
            exact_matrix_and_nullspace_error=max(errors),mapped_limit_error=mapped_error,
            unmapped_projection_error=naive_error,smallest_fifth_transfer_eigenvalue=float(min(np.linalg.eigvalsh(t)))))
    return dict(rows=rows,all_original_internal_matrices_retained=True)


def grassmann_check():
    d=data(nx=2,nt=2,a=.2);length=2;e=fields(4);e0=np.zeros_like(e);e0[:,0]=1
    r=schur_frame(d,length);full=full_pfaffian(d,length,e);base=full_pfaffian(d,length,e0)
    actual=full/base;reduced=ratio(d,r,e);err=float(abs(actual-reduced))
    assert err<2e-12 and abs(base)>1e-100
    return dict(nx=2,nt=2,L5=2,a5=.2,full_Nambu_dimension=1024,
        boundary_fields=e.tolist(),full_signed_ratio=[float(actual.real),float(actual.imag)],
        reduced_signed_ratio=[float(reduced.real),float(reduced.imag)],identity_error=err,
        full_reference_Pfaffian=[float(base.real),float(base.imag)],
        normalization_not_replaced_by_modulus=True)


def limit_check():
    rows=[];bounds=[]
    target=old.frame(1,3,2);e=fields(6,65902);target_ratio=old.ratio(target,e)
    for a in (.2,.1,.05,.025):
        d=data(a=a);i=d['eye'];m=d['M'];ht=d['HT'];ev=d['ev'];un=d['minus'];up=d['plus']
        qinf,_=np.linalg.qr(np.linalg.solve(2*i+a*d['X'],un));pinf=qinf@qinf.conj().T
        hh=op(target['H']);gap=float(min(abs(target['ev'])));tt=a*hh/2
        delta=a*hh*hh/(2*(1-tt));assert delta<gap
        projector_bound=delta/(2*(gap-delta))+tt/(1-tt)
        assert op(pinf-target['P'])<=projector_bound+1e-12
        coeff=op((up.conj().T@d['um'])@np.linalg.inv(un.conj().T@d['um']))
        positive=ev[ev>0];negative=ev[ev<0]
        q=float(max((1-positive)/(1+positive))/min((1-negative)/(1+negative)))
        assert 0<q<1
        for length in (8,32,128,512):
            r=schur_frame(d,length);p=r@r.conj().T
            filter_bound=float(np.linalg.cond(m)*coeff*q**(length-1))
            actual_error=op(p-pinf);assert actual_error<filter_bound+3e-12
            rows.append(dict(a5=a,L5=length,boundary_to_original_projector_error=op(p-target['P']),
                finite_depth_to_fixed_a_limit_error=actual_error,filter_bound=filter_bound))
        rr=ratio(d,r,e)
        bounds.append(dict(a5=a,L5=512,boundary_limit_to_original_error=op(pinf-target['P']),
            rigorous_bias_bound=projector_bound,original_signed_ratio=[float(target_ratio.real),float(target_ratio.imag)],
            mapped_signed_ratio=[float(rr.real),float(rr.imag)],ratio_error=float(abs(rr-target_ratio))))
    assert bounds[-1]['ratio_error']<bounds[0]['ratio_error']
    return dict(depth_rows=rows,joint_limit_rows=bounds,
        joint_limit_requires_a5_to_zero_and_a5_times_L5_to_infinity=True,
        no_uniform_time_volume_or_continuum_limit_claim=True)


def normalization_source_check():
    # The original spatial Wilson coefficient is a source, not auxiliary E.
    length=2;a=.2;h=2e-4;e0=np.zeros((4,10));e0[:,0]=1
    logs=[];bulk_logs=[]
    for k in (1-h,1+h):
        d=data(kappa=k,nx=2,nt=2,a=a)
        logs.append(float(np.log(abs(full_pfaffian(d,length,e0)))))
        bulk_logs.append(float(16*np.linalg.slogdet(modified_bulk(d,length,ap=True))[1]))
    numerator=(logs[1]-logs[0])/(2*h);bulk_fd=(bulk_logs[1]-bulk_logs[0])/(2*h)
    d=data(nx=2,nt=2,a=a);dx=d['g']@d['f']['dH']
    dd=np.kron(np.eye(length),a*dx)
    bulk_exact=float(16*np.trace(np.linalg.solve(modified_bulk(d,length,ap=True),dd)).real)
    assert abs(bulk_exact-bulk_fd)<2e-6 and abs(bulk_exact)>1
    remainder=numerator-bulk_exact;assert abs(remainder)>1
    # Response of the original E weight on two physical times also converges.
    response=[];angles=np.array([.12,-.17,.08,.21,-.09,.03])
    def logratio(k,a5=None):
        f=old.frame(k,3,2);e,_=old.plane(f,angles)
        if a5 is None:return float(np.log(abs(old.ratio(f,e))))
        q=data(kappa=k,a=a5);return float(np.log(abs(ratio(q,schur_frame(q,512),e))))
    target=(logratio(1+h)-logratio(1-h))/(2*h)
    for a5 in (.1,.05,.025):
        actual=(logratio(1+h,a5)-logratio(1-h,a5))/(2*h)
        response.append(dict(a5=a5,finite_difference_source=actual,original_source=target,error=abs(actual-target)))
    assert response[-1]['error']<response[0]['error']
    return dict(nx=2,nt=2,a5=a,L5=length,source='original spatial Wilson coefficient kappa',
        full_reference_log_Pfaffian_derivative=float(numerator),
        AP_bulk_log_determinant_derivative=bulk_exact,AP_bulk_finite_difference=float(bulk_fd),
        AP_derivative_error=float(abs(bulk_exact-bulk_fd)),
        retained_reference_log_weight_derivative=float(remainder),
        response_rows=response,source_dependent_reference_factor_cannot_be_discarded=True,
        general_geometry_derivative_limit_not_claimed=True)


def run():
    deps=('joint_spatial_auxiliary_geometry.py','joint_subgroup_measure_source.py',
          'research_note_612.md','research_note_656.md','research_note_658.md',
          'round659_drafts/temporal_memory_probe_results.json')
    return dict(date='2026-10-02',round=659,tests_run=4,failures=0,errors=0,
        exact_boundary=boundary_identity_check(),full_Grassmann=grassmann_check(),
        controlled_boundary_limit=limit_check(),common_normalization_source=normalization_source_check(),
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope='Specified free finite-volume domain-wall boundary maps to original signed S9 Pfaffian, with controlled auxiliary-depth/spacing limit and explicit source-dependent subtraction. No positive common physical-time semigroup, physical 32CAR identification, interacting gauge construction or quantum GR.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args();r=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:json.dump(r,f,ensure_ascii=False,indent=2)
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps({k:r[k] for k in ('round','tests_run','failures','errors')}))
