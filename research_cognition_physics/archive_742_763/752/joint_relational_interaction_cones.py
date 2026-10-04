"""752: conditional compatibility of the original derivative reference and a new local interaction.

Original 573/647 fields are reused. The added action is an explicit candidate,
not the old diagnostic observable and not a realized quantum instrument.
"""
import argparse, hashlib, json
from pathlib import Path
import numpy as np
import joint_dynamic_continuum_reference as inherited

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_relational_interaction_cones_results.json'
ETA=np.diag([-1.,1.,1.,1.])
EPS=.02
WIDTH=np.array([.1,.1])

def plateau(x):
    q=float(x*x)
    if q<=1: return 1.
    if q>=4: return 0.
    t=(4-q)/3
    a=np.exp(-1/t);b=np.exp(-1/(1-t))
    return float(a/(a+b))

def data():
    p=inherited.point();d=inherited.initial();i=p['index']
    h,s=p['phi'][[1,4]];F=2-(h*h+s*s)/6
    K=np.eye(2)/F+np.outer([h,s],[h,s])/(6*F*F)
    # Covectors in the original point's orthonormal frame; no new flat solution.
    q=np.zeros((2,4));q[:,0]=p['v'][[1,4]]
    q[0,1:]=d['f']['dh'][i]/d['psi'][i]**2
    q[1,1:]=d['f']['ds'][i]/d['psi'][i]**2
    AB=np.einsum('ai,ij,aj->a',q,ETA,q)
    return dict(h=float(h),s=float(s),F=float(F),K=K,q=q,AB=AB,
                W=float(.5+np.sin(s)/4))

def window(AB,center,rho=0.):
    z=(AB-center)/WIDTH
    C=np.array([[1.,rho],[rho,1.]])
    f=np.exp(-.5*z@C@z)*plateau(z[0])*plateau(z[1])
    # Derivative checks use the open plateau, not a nonsmooth cutoff.
    return float(f)

def window_jet(AB,center,rho=0.):
    z=(AB-center)/WIDTH
    assert max(abs(z))<1.
    C=np.array([[1.,rho],[rho,1.]])
    f=window(AB,center,rho);d=-C@z/WIDTH
    return f,f*d,f*(np.outer(d,d)-C/np.outer(WIDTH,WIDTH))

def lagrangian(q,gi,K,W,center,rho):
    AB=np.einsum('ai,ij,aj->a',q,gi,q)
    return float(-.5*np.einsum('ab,ai,ij,bj',K,q,gi,q)-EPS*W*window(AB,center,rho))

def principal(q,gi,K,W,center,rho):
    AB=np.einsum('ai,ij,aj->a',q,gi,q)
    f,df,ddf=window_jet(AB,center,rho);up=q@gi
    P=np.einsum('ab,ij->aibj',K+2*EPS*W*np.diag(df),gi)
    for a in range(2):
        for b in range(2):
            P[a,:,b,:]+=4*EPS*W*ddf[a,b]*np.outer(up[a],up[b])
    return P

def fd_hessian(fun,x,step):
    x=np.asarray(x,float);n=x.size;out=np.zeros((n,n));base=fun(x)
    unit=np.eye(n)*step
    for i in range(n):
        out[i,i]=(fun(x+unit[i])-2*base+fun(x-unit[i]))/step**2
        for j in range(i):
            out[i,j]=out[j,i]=(fun(x+unit[i]+unit[j])-fun(x+unit[i]-unit[j])
                              -fun(x-unit[i]+unit[j])+fun(x-unit[i]-unit[j]))/(4*step**2)
    return out

def source_principal_check(d):
    q=d['q'];K=d['K'];W=d['W'];center=d['AB']+WIDTH*np.array([.21,-.17]);rho=.37
    P=principal(q,ETA,K,W,center,rho).reshape(8,8)
    fun=lambda x:lagrangian(x.reshape(2,4),ETA,K,W,center,rho)
    steps=(2e-4,1e-4,5e-5)
    errs=[float(np.max(abs(fd_hessian(fun,q.ravel(),t)+P))) for t in steps]
    assert errs[-1]<2e-7 and errs[-1]<errs[0]/6
    f,df,ddf=window_jet(d['AB'],center,rho)
    T=-EPS*W*f*ETA+2*EPS*W*sum(df[a]*np.outer(q[a],q[a]) for a in range(2))
    # Independently vary sqrt(-det g)*Delta L with respect to inverse metric.
    def density(gi):
        return -EPS*W*window(np.einsum('ai,ij,aj->a',q,gi,q),center,rho)/np.sqrt(-np.linalg.det(gi))
    metric_errors=[]
    for t in (1e-4,5e-5,2.5e-5):
        err=0.
        for i in range(4):
            for j in range(i,4):
                E=np.zeros((4,4));E[i,j]=1;E[j,i]=1
                deriv=(density(ETA+t*E)-density(ETA-t*E))/(2*t)
                predicted=-.5*np.sum(T*E)
                err=max(err,abs(deriv-predicted))
        metric_errors.append(float(err))
    assert metric_errors[-1]<2e-9
    # Coordinate tensor covariance, including nonzero mixed window derivatives.
    L=np.array([[1.2,.1,0,.03],[.2,1.1,.1,0],[0,.05,.9,.02],[.01,0,.04,1.03]])
    gi=L@ETA@L.T;qnew=q@np.linalg.inv(L)
    Pnew=principal(qnew,gi,K,W,center,rho)
    transported=np.einsum('mi,aibj,nj->ambn',L,P.reshape(2,4,2,4),L)
    cov_error=float(np.max(abs(Pnew-transported)))
    assert cov_error<2e-14
    assert abs(ddf[0,1])>1
    return dict(hessian_steps=list(steps),gradient_hessian_max_errors=errs,
                metric_source_max_errors=metric_errors,tensor_covariance_max_error=cov_error,
                mixed_window_second_derivative=float(ddf[0,1]),
                added_stress_covariant=T.tolist(),
                scope='Local source and principal tensor of the explicitly added action; not a solved new Einstein history.')

def cone_check(d):
    P=principal(d['q'],ETA,d['K'],d['W'],d['AB'],0.)
    kinetic=-P[:,0,:,0]
    spatial=P[:,1:,:,1:].reshape(6,6)
    kt=np.linalg.eigvalsh(kinetic);sp=np.linalg.eigvalsh(spatial)
    assert min(kt)>0 and min(sp)>0
    # The original point has zero h_x,s_x (up to cos(pi/2) roundoff).
    xblock=P[:,1,:,1];mixed=P[:,0,:,1]+P[:,1,:,0]
    assert np.max(abs(mixed))<1e-16
    ev,U=np.linalg.eigh(kinetic);isqrt=(U/np.sqrt(ev))@U.T
    speed2=np.linalg.eigvalsh(isqrt@xblock@isqrt)
    assert 0<speed2[0]<speed2[1]<1
    null=np.array([-1.,1.,0.,0.])
    null_P=np.einsum('aibj,i,j->ab',P,null,null)
    assert np.linalg.eigvalsh(null_P).max()<-1e-6
    roots_error=0.;roots_max_imag=0.;minimum_root_abs=2.
    for j in range(21):
        theta=.27+j*.19
        n=np.array([np.cos(theta)*np.sin(.4+theta),np.sin(theta)*np.sin(.4+theta),np.cos(.4+theta)])
        A=P[:,0,:,0]
        B=-np.einsum('abi,i->ab',P[:,0,:,1:]+P[:,1:,:,0].transpose(0,2,1),n)
        C=np.einsum('aibj,i,j->ab',P[:,1:,:,1:],n,n)
        companion=np.block([[np.zeros((2,2)),np.eye(2)],[-np.linalg.solve(A,C),-np.linalg.solve(A,B)]])
        roots=np.linalg.eigvals(companion)
        roots_max_imag=max(roots_max_imag,float(np.max(abs(roots.imag))))
        minimum_root_abs=min(minimum_root_abs,float(np.min(abs(roots))))
        for z in roots:
            roots_error=max(roots_error,abs(np.linalg.det(A*z*z+B*z+C)))
    assert roots_max_imag<1e-12 and roots_error<1e-12 and minimum_root_abs>0
    # Positive class: affine A/B interactions, no Hessian of window.
    v=np.array([.3,-.2]);Keff=d['K']+2*EPS*d['W']*np.diag(v)
    assert np.linalg.eigvalsh(Keff).min()>0
    Paff=np.einsum('ab,ij->aibj',Keff,ETA)
    aff_null=float(np.max(abs(np.einsum('aibj,i,j->ab',Paff,null,null))))
    assert aff_null==0
    return dict(epsilon=EPS,widths=WIDTH.tolist(),kinetic_eigenvalues=kt.tolist(),
                spatial_energy_eigenvalues=sp.tolist(),speed_squared_in_original_x_frame=speed2.tolist(),
                speeds_in_original_x_frame=np.sqrt(speed2).tolist(),
                original_metric_null_principal_eigenvalues=np.linalg.eigvalsh(null_P).tolist(),
                sampled_directions=21,maximum_characteristic_root_imaginary_part=roots_max_imag,
                maximum_characteristic_determinant_residual=float(roots_error),
                affine_candidate_target_eigenvalues=np.linalg.eigvalsh(Keff).tolist(),
                affine_common_null_symbol_error=aff_null,
                scope='Original on-shell initial jet used to test a new action locally. No full perturbed constrained solution or instrument claimed.')

def run():
    d=data();a=source_principal_check(d);b=cone_check(d)
    deps=('research_note_352.md','research_note_573.md','research_note_647.md','research_note_731.md',
          'joint_gravity_material_coordinates.py','joint_dynamic_continuum_reference.py')
    return dict(round=752,tests_run=2,failures=0,errors=0,
                original_jet=dict(h=d['h'],s=d['s'],F=d['F'],W=d['W'],covectors=d['q'].tolist(),
                                  derivative_labels=d['AB'].tolist(),radial_target_metric=d['K'].tolist()),
                common_source_and_principal=a,conditional_common_cone=b,
                dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
                scope='Conditional incompatibility of promoting the original derivative relational window into a local action while demanding both radial polarizations retain the exact metric cone throughout its support. Affine derivative windows are the allowed positive-target class; a nonzero compact four-label window is excluded. A known k-essence principal-symbol method is applied, not a new general theorem. This is neither failure of relational observables nor quantum-gravity construction, and it does not replace the frozen original action.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args()
    r=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8') as f:json.dump(r,f,ensure_ascii=False,indent=2)
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps({k:r[k] for k in ('round','tests_run','failures','errors')}))

