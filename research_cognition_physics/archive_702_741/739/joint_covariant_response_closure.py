"""739: common neutral/metric response and finite-band causal kernel identity.

Two checks: a covariant twelve-component symbol including original contacts,
and the exact spatial-frequency shift of full mass/shear spectral kernels.
No nonlinear solution, physical polarization count or UV-uniform bound.
"""
import argparse
import hashlib
import json
from functools import lru_cache
from pathlib import Path
import numpy as np
import joint_mixed_neutral_response as mixed
import joint_local_source_normalization as massmap

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_covariant_response_closure_results.json'


def maxabs(a):return float(np.max(np.abs(a)))


@lru_cache(maxsize=1)
def basis():
    B=[]
    for i in range(4):
        e=np.zeros((4,4));e[i,i]=1.;B.append(e)
    for i in range(4):
        for j in range(i):
            e=np.zeros((4,4));e[i,j]=e[j,i]=1/np.sqrt(2);B.append(e)
    return np.array(B)


@lru_cache(maxsize=1)
def local_data():
    old=mixed.old;x=mixed.original_data()[0]
    phi=np.array([0.,np.sqrt(old.U0[0]),0.,0.,np.sqrt(old.U0[1])])
    F=old.matter.original.F(phi)
    invJ=np.linalg.inv(massmap.jacobian(phi))[np.ix_([1,4],[1,4])]
    Z=invJ.T@old.matter.original.metric(phi)[np.ix_([1,4],[1,4])]@invJ
    _,_,Wh=old.W_jets(old.U0)
    D=np.diag(2*np.sqrt(old.U0))
    Hphi=D@(old.L/2+Wh)@D/F**2
    H=invJ.T@Hphi@invJ
    return x,Z,H


def projectors(p):
    w=p@p;theta=np.eye(4)-np.outer(p,p)/w;B=basis()
    tv=np.einsum('aij,ij->a',B,theta)
    P0=np.outer(tv,tv)/3
    P2=np.array([[np.sum(a*(theta@b@theta)) for b in B] for a in B])-P0
    return w,theta,tv,P0,P2


def direct_local_metric(p,f,alpha,beta):
    """Fierz-Pauli plus independently contracted linear curvature tensors."""
    B=basis();w=p@p;tr=np.trace(B,axis1=1,axis2=2)
    q=np.einsum('i,aij,j->a',p,B,p);D=np.einsum('aij,j->ia',B,p)
    fp=f/4*(w*np.eye(10)-w*np.outer(tr,tr)+np.outer(tr,q)+np.outer(q,tr)-2*D.T@D)
    scalar=w*tr-q;ricci=[];riemann=[]
    for e in B:
        d=e@p
        ricci.append((w*e+np.outer(p,p)*np.trace(e)-np.outer(p,d)-np.outer(d,p))/2)
        R=np.empty((4,4,4,4),complex)
        for a in range(4):
            for b in range(4):
                for c in range(4):
                    for d0 in range(4):
                        R[a,b,c,d0]=(-p[c]*p[b]*e[a,d0]-p[d0]*p[a]*e[b,c]
                                      +p[d0]*p[b]*e[a,c]+p[c]*p[a]*e[b,d0])/2
        riemann.append(R)
    ricci=np.array(ricci).reshape(10,-1);riemann=np.array(riemann).reshape(10,-1)
    # Analytic bilinear symbols: transpose, NOT Hermitian conjugation.
    Cgram=riemann@riemann.T-2*ricci@ricci.T+np.outer(scalar,scalar)/3
    return fp+2*alpha*np.outer(scalar,scalar)+2*beta*Cgram,scalar


def full_symbol(p,f,ell,Z,H,alpha,beta):
    x=local_data()[0];w,theta,tv,P0,P2=projectors(p);Q=np.sqrt(w+0j)
    Rm=w*mixed.kernel(Q)
    _,_,_,m,d=mixed.original_data()
    R2=mixed.tensor.tensor_dispersion(1j*Q,m,d)
    local,scalar=direct_local_metric(p,f,alpha,beta)
    full=np.zeros((12,12),complex)
    full[:2,:2]=w*Z+H-Rm
    cross=-np.outer(ell,scalar)/2-np.outer(Rm@x,tv)/6
    full[:2,2:]=cross;full[2:,:2]=cross.T
    full[2:,2:]=local-(x@Rm@x)*np.outer(tv,tv)/36-R2*P2/4
    projected_local=f*w/4*(P2-2*P0)+6*alpha*w*w*P0+beta*w*w*P2
    return full,Rm,R2,maxabs(local-projected_local)


def covariant_symbol_check():
    x,Z0,H=local_data();errors=[];ward=[];active=[];nulls=[]
    fixtures=[(1.,np.zeros(2),Z0,0.,0.),
              (1.03,np.array([.02,-.01]),Z0+np.diag([.04,.01]),.02,-.007)]
    for f,ell,Z,alpha,beta in fixtures:
        ae=float(x@Z@x+6*ell@x-6*f)
        assert abs(ae)>1
        nulls.append(dict(f=f,f_gradient=ell.tolist(),Z=Z.tolist(),alpha=alpha,beta=beta,
                          complementary_second_order_coefficient=ae,
                          complementary_fourth_order_coefficient=72*alpha))
        for p in (np.array([.8,.2,-.3,.4]),np.array([1.7,-.5,.6,.1]),
                  np.array([1.1+.3j,.2,-.3,.4])):
            full,Rm,R2,err=full_symbol(p,f,ell,Z,H,alpha,beta);errors.append(err)
            for xi in np.eye(4):
                gauge=np.outer(p,xi)+np.outer(xi,p)
                vec=np.r_[np.zeros(2),np.einsum('aij,ij->a',basis(),gauge)]
                ward.append(maxabs(full@vec))
            if np.iscomplexobj(p):continue
            w=p@p;n=p/np.sqrt(w);v=np.eye(4)[0]-n
            O=np.eye(4)-2*np.outer(v,v)/(v@v)
            T=O[:,1:];theta=T@T.T
            e=np.array([T@b@T.T for b in mixed.tensor.stress_basis()[1:]])
            lift=np.zeros((12,8));lift[:2,:2]=np.eye(2);lift[:2,2]=-x
            lift[2:,2]=2*np.einsum('aij,ij->a',basis(),theta)
            lift[2:,3:]=2*np.einsum('aij,bij->ab',basis(),e)
            predicted=np.zeros((8,8),complex)
            predicted[:2,:2]=w*Z+H-Rm
            mixed_local=-(Z@x+3*ell)*w-H@x
            predicted[:2,2]=mixed_local;predicted[2,:2]=mixed_local
            predicted[2,2]=ae*w+x@H@x+72*alpha*w*w
            predicted[3:,3:]=np.eye(5)*(f*w+4*beta*w*w-R2)
            active.append(maxabs(lift.T@full@lift-predicted))
    assert max(errors)<2e-13 and max(ward)<2e-13 and max(active)<3e-13
    return dict(independent_curvature_projector_error=max(errors),
                full_twelve_component_diffeomorphism_Ward_error=max(ward),
                original_active_eight_component_pullback_error=max(active),
                original_matched_mass_coordinate_potential_Hessian=H.tolist(),
                declared_finite_coefficient_fixtures=nulls,
                fixture_coefficients_not_predictions=True,
                no_full_internal_gauge_or_nonlinear_constraint_claim=True)


def channel_terms(channel):
    if channel=='neutral':return mixed.terms()
    _,_,_,m,d=mixed.original_data()
    return [(f'shear{i}',2*mi,0.,'T',np.array([[di/(80*np.pi**2)]]))
            for i,(mi,di) in enumerate(zip(m,d))]


def profile(x,r,kind):
    if kind=='T':return (np.maximum(0.,1-x*x))**1.5*(1+2*x*x/3)
    return mixed.profile(x,r,kind)


def invariant_kernel(s,channel):
    if channel=='neutral':return mixed.kernel(s)
    if abs(s)==0:return np.zeros((1,1),complex)
    _,_,_,m,d=mixed.original_data()
    return np.array([[-mixed.tensor.tensor_dispersion(1j*s,m,d)/s**4]])


def shifted_decomposition(s,k,channel):
    x,q=mixed.quad(192);terms=channel_terms(channel)
    size=terms[0][-1].shape[0]
    C=np.zeros((size,size));static=invariant_kernel(k,channel)
    I=static.copy();B=static.copy();R=np.zeros_like(I)
    for _,a,r,kind,Cj in terms:
        threshold=np.hypot(a,k)
        arg=a*x/np.sqrt(threshold**2-k*k*x*x)
        g=profile(arg,r,kind);gap=1-g
        C+=Cj
        I+=Cj*np.dot(q,s*s*x*g/(threshold**2+s*s*x*x))
        B+=Cj*(-np.log(threshold)-np.dot(q,gap/x))
        R+=Cj*(.5*np.log(1+threshold**2/s**2)+
               np.dot(q,gap*threshold**2/(x*(threshold**2+s*s*x*x))))
    return I,C,B,R


def spatial_memory_check():
    rows=[]
    for channel in ('neutral','shear'):
        shift_error=0.;constant_error=0.;split_error=0.
        _,C,B0,_=shifted_decomposition(1.2+.3j,0.,channel)
        for k,s in ((.2,1.1+.3j),(.7,.8+.9j),(1.,2.3-.5j)):
            I,_,B,R=shifted_decomposition(s,k,channel)
            actual=invariant_kernel(np.sqrt(s*s+k*k),channel)
            shift_error=max(shift_error,maxabs(actual-I))
            constant_error=max(constant_error,maxabs(B-B0))
            split_error=max(split_error,maxabs(I-C*np.log(s)-B0-R))
        assert max(shift_error,constant_error,split_error)<3e-12
        eig,V=np.linalg.eigh(C);W=(V*(1/np.sqrt(eig)))@V.T
        K=1.;T=.01;memory=0.
        for _,a,r,kind,Cj in channel_terms(channel):
            AK=np.hypot(a,K);assert AK*T<1
            memory+=np.linalg.norm(W@Cj@W,2)*AK*AK*T*T*(7/4+np.log(1/(AK*T)))
        bnorm=float(np.linalg.norm(W@B0@W,2))
        rate=mixed.causal.inverse_norm_upper(T)*(bnorm+memory)
        assert rate<1
        rows.append(dict(channel=channel,exact_spatial_shift_error=shift_error,
                         unchanged_constant_matrix_error=constant_error,
                         exact_shifted_memory_decomposition_error=split_error,
                         declared_spatial_band=K,time_window=T,
                         normalized_memory_uniform_L1_upper=float(memory),
                         normalized_constant_norm=bnorm,pure_kernel_contraction_upper=float(rate)))
    return dict(rows=rows,all_masses_from_the_same632_background=True,
                finite_band_uniform_bound_not_K_to_infinity=True,
                local_gravity_terms_not_included_in_reported_numeric_contraction=True)


def run():
    names=('covariant_symbol_check','spatial_memory_check')
    result={name:globals()[name]() for name in names}
    deps=('research_note_583.md','research_note_601.md','research_note_630.md',
          'research_note_631.md','research_note_632.md','research_note_732.md',
          'research_note_735.md','research_note_737.md','research_note_738.md',
          'joint_mixed_neutral_response.py','joint_background_contact_matching.py',
          'joint_tensor_stress_spectrum.py','joint_local_source_normalization.py')
    return dict(round=739,tests_run=2,failures=0,errors=0,checks=list(names),results=result,
                dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
                scope='Common stationary neutral/metric one-fermion-loop response including covariant local jets, linear diffeomorphism Ward identities and finite-spatial-band causal memory. Conditional retarded linear closure modulo coordinate gauge; no all-scale bound, all internal gauge dynamics, nonlinear integrability or physical stability theorem.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true')
    args=p.parse_args();r=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
