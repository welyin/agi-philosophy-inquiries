"""618: original five-scalar nonminimal boundary momenta and common gluing.

Classical two-derivative Lorentz branch; timelike smooth boundary, spacelike
unit normal. This does not quantize gravity or derive the Einstein action.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_curved_quantum_source as old

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_geometric_boundary_matching_results.json'
M=old.M
ETA=np.diag([-1.,1.,1.])

def data(phi):
    return float(old.F(phi)),-np.asarray(phi)/3

def jordan_momenta(phi,h,Kcov,v):
    F,f=data(phi);hi=np.linalg.inv(h);tr=float(np.einsum('ab,ab',hi,Kcov))
    Ku=hi@Kcov@hi
    P=F*(tr*hi-Ku)+(f@v)*hi  # 2 Pi^{ab}/sqrt|h|
    j=f*tr-v                # pi_a/sqrt|h|
    return P,j

def invert_momenta(phi,h,P,j):
    F,f=data(phi);p=float(np.einsum('ab,ab',h,P))
    tr=(p+3*f@j)/(2*M)
    v=f*tr-j
    Ptf=P-p*np.linalg.inv(h)/3
    Kcov=tr*h/3-h@Ptf@h/F
    return Kcov,v

def einstein_data(phi,h,Kcov,v):
    F,f=data(phi)
    hE=F*h
    KE=np.sqrt(F)*(Kcov+h*(f@v)/(2*F))
    vE=v/np.sqrt(F)
    hi=np.linalg.inv(hE);tr=float(np.einsum('ab,ab',hi,KE))
    PE=tr*hi-hi@KE@hi
    jE=-old.metric(phi)@vE
    return hE,KE,vE,PE,jE

def momentum_check():
    rng=np.random.default_rng(618);rows=[];worst=0.
    for fraction in (.05,.35,.7,.95):
        direction=rng.normal(size=5);direction/=np.linalg.norm(direction)
        phi=direction*np.sqrt(6*M)*fraction
        T=np.eye(3)+.15*rng.normal(size=(3,3));h=T.T@ETA@T
        K=rng.normal(size=(3,3));K=(K+K.T)/2
        v=rng.normal(size=5)*.2
        F,f=data(phi);P,j=jordan_momenta(phi,h,K,v)
        Ki,vi=invert_momenta(phi,h,P,j)
        hE,KE,vE,PE,jE=einstein_data(phi,h,K,v)
        root=np.sqrt(abs(np.linalg.det(h)));rootE=np.sqrt(abs(np.linalg.det(hE)))
        Pi=.5*root*P;pi=root*j;PiE=.5*rootE*PE;piE=rootE*jE
        # Cotangent chain rule includes scalar work from changing h_E=F h_J.
        piframe=piE+np.einsum('ab,ab',PiE,h)*f
        dh=rng.normal(size=(3,3));dh=(dh+dh.T)/2
        dphi=rng.normal(size=5)
        dhE=F*dh+(f@dphi)*h
        workJ=np.einsum('ab,ab',Pi,dh)+pi@dphi
        workE=np.einsum('ab,ab',PiE,dhE)+piE@dphi
        errors=[np.max(abs(K-Ki)),np.max(abs(v-vi)),np.max(abs(Pi-F*PiE)),
                np.max(abs(pi-piframe)),abs(workJ-workE),abs(F+1.5*(f@f)-M)]
        # Same Schur quantity underlies the old target metric/inverse.
        kinetic=np.eye(5)/F+1.5*np.outer(f,f)/F**2
        errors.extend([np.max(abs(kinetic-old.metric(phi))),
                       np.max(abs(kinetic@old.inverse(phi)-np.eye(5)))])
        err=float(max(errors));worst=max(worst,err);assert err<2e-11
        Hess=np.block([[np.array([[6*F]]),3*f[None,:]],[3*f[:,None],-np.eye(5)]])
        assert abs(np.linalg.det(Hess)+6*M)<2e-12
        rows.append(dict(radial_fraction=fraction,F=F,Schur=F+1.5*float(f@f),
            normal_trace_scalar_Hessian_determinant=float(np.linalg.det(Hess)),max_error=err,
            omitted_scalar_cotangent_term_norm=float(np.linalg.norm(piframe-piE))))
    return dict(rows=rows,max_error=worst,full_five_scalar_inverse=True,
        Hessian_is_not_claimed_positive=True,common_boundary_work_preserved=True)

def vacuum():
    _,u,_=old.lattice.scalar.parameters()
    return np.array([0.,np.sqrt(u[0]),0.,0.,np.sqrt(u[1])])

# Coefficients of a polynomial q(x), q0+q1*x+q2*x^2; independent curvature eval.
A=np.array([.03,.12,-.07])
PHI=np.stack([vacuum()+np.array([.06,-.04,.03,.02,-.05]),
              np.array([.08,.11,-.03,.04,.07]),np.array([-.03,.02,.04,-.01,.03])])

def profiles(x,variation=0.):
    # Endpoint variations are independent of normal derivatives.
    da=np.array([.13,-.21,0.])
    dp=np.array([[.03,.02,-.04,.01,.05],[-.06,.01,.02,-.03,.02],[0.,0.,0.,0.,0.]])
    aa=A+variation*da;pp=PHI+variation*dp
    a=aa[0]+aa[1]*x+aa[2]*x*x;ap=aa[1]+2*aa[2]*x
    app=np.full_like(x,2*aa[2])
    phi=pp[0]+x[...,None]*pp[1]+x[...,None]**2*pp[2]
    v=pp[1]+2*x[...,None]*pp[2];vp=np.broadcast_to(2*pp[2],phi.shape)
    return a,ap,app,phi,v,vp

def potential_gradient(phi):
    L,u,_=old.lattice.scalar.parameters()
    delta=np.stack([np.sum(phi[...,:4]**2,axis=-1)-u[0],phi[...,4]**2-u[1]],axis=-1)
    w=delta@L
    return np.concatenate([phi[...,:4]*w[...,:1],phi[...,4:]*w[...,1:]],axis=-1)

def action(variation=0.,order=64):
    nodes,weights=np.polynomial.legendre.leggauss(order);x=(nodes+1)/2;weights=weights/2
    a,w,wp,phi,v,vp=profiles(x,variation)
    F=old.F(phi);f=-phi/3
    Fp=np.einsum('...i,...i->...',f,v)
    Fpp=-np.sum(v*v,axis=-1)/3+np.einsum('...i,...i->...',f,vp)
    V=old.node_potential(phi)*F**2
    RJ=-6*wp-12*w*w
    bulkJ=float(weights@(np.exp(3*a)*(F*RJ/2-np.sum(v*v,axis=-1)/2-V)))
    reduced=float(weights@(np.exp(3*a)*(3*F*w*w+3*Fp*w-np.sum(v*v,axis=-1)/2-V)))
    N=np.sqrt(F);b=a+np.log(F)/2
    bp=w+Fp/(2*F);bpp=wp+Fpp/(2*F)-Fp**2/(2*F**2)
    RE=-6/F*(bpp-bp*Fp/(2*F)+2*bp*bp)
    kin=np.einsum('...i,...ij,...j->...',v,old.metric(phi),v)
    bulkE=float(weights@(N*np.exp(3*b)*(RE/2-kin/(2*F)-V/F**2)))
    ae,we,_,pe,ve,_=profiles(np.array([0.,1.]),variation)
    Fe=old.F(pe);Fpe=-np.sum(pe*ve,axis=-1)/3
    endpointJ=3*np.exp(3*ae)*Fe*we
    endpointE=3*np.exp(3*ae)*Fe*(we+Fpe/(2*Fe))
    boundaryJ=float(endpointJ[1]-endpointJ[0]);boundaryE=float(endpointE[1]-endpointE[0])
    return dict(J=bulkJ+boundaryJ,E=bulkE+boundaryE,reduced=reduced,
        bulkJ=bulkJ,bulkE=bulkE,boundaryJ=boundaryJ,boundaryE=boundaryE,
        min_F=float(F.min()))

def variation_rhs():
    nodes,weights=np.polynomial.legendre.leggauss(64);x=(nodes+1)/2;weights/=2
    a,w,wp,phi,v,vp=profiles(x);F=old.F(phi);f=-phi/3
    Fp=np.einsum('...i,...i->...',f,v)
    Fpp=-np.sum(v*v,axis=-1)/3+np.einsum('...i,...i->...',f,vp)
    V=old.node_potential(phi)*F**2
    ell=3*F*w*w+3*Fp*w-np.sum(v*v,axis=-1)/2-V
    ELa=np.exp(3*a)*(3*ell-3*w*(6*F*w+3*Fp)-6*Fp*w-6*F*wp-3*Fpp)
    ELphi=np.exp(3*a)[...,None]*(
        vp+3*w[...,None]*v-3*f*wp[...,None]-6*f*w[...,None]**2-potential_gradient(phi))
    da=.13-.21*x
    dp=np.array([.03,.02,-.04,.01,.05])+x[:,None]*np.array([-.06,.01,.02,-.03,.02])
    bulk=float(weights@(ELa*da+np.einsum('...i,...i->...',ELphi,dp)))
    x=np.array([0.,1.]);a,w,wp,phi,v,vp=profiles(x)
    F=old.F(phi);f=-phi/3;Fp=np.sum(f*v,axis=-1)
    pa=np.exp(3*a)*(6*F*w+3*Fp)
    pi=np.exp(3*a)[:,None]*(3*f*w[:,None]-v)
    da=.13-.21*x
    dp=np.array([.03,.02,-.04,.01,.05])+x[:,None]*np.array([-.06,.01,.02,-.03,.02])
    bw=pa*da+np.einsum('...i,...i->...',pi,dp)
    boundary=float(bw[1]-bw[0])
    return bulk,boundary

def action_check():
    base=action();higher=action(order=96)
    err=max(abs(base['J']-base['E']),abs(base['J']-base['reduced']),
            abs(base['J']-higher['J']))
    bulk,boundary=variation_rhs();rows=[]
    for step in (2e-3,1e-3,5e-4):
        derivative=(action(step)['J']-action(-step)['J'])/(2*step)
        error=abs(derivative-bulk-boundary)
        rows.append(dict(step=step,action_derivative=derivative,
                         variation_identity_error=float(error)))
    assert err<5e-13 and rows[-1]['variation_identity_error']<2e-8
    assert rows[-1]['variation_identity_error']<.28*rows[-2]['variation_identity_error']
    assert abs(base['bulkJ']-base['bulkE'])>1e-3 and abs(boundary)>1e-3
    return dict(original_nonperiodic_action=base,action_equality_error=float(err),
        variation_bulk=bulk,variation_boundary=boundary,derivative_rows=rows,
        omitted_GHY_frame_defect=abs(base['bulkJ']-base['bulkE']),
        actual_old_potential_retained=True,off_shell_action_identity_not_global_solution=True)

def interface_check():
    phi=vacuum();F,f=data(phi);h=ETA;hi=ETA;tau=.07
    # A test surface action -tau int sqrt|h_J|. Extra diagnostic input only.
    P=tau*hi;j=np.zeros(5)
    K,v=invert_momenta(phi,h,P,j)
    hE,KE,vE,PE,jE=einstein_data(phi,h,K,v)
    surfacePE=-tau*F**(-1.5)*np.linalg.inv(hE)
    surfacejE=1.5*tau*F**(-2.5)*f
    errors=[np.max(abs(PE+surfacePE)),np.max(abs(jE+surfacejE))]
    Pback,jback=jordan_momenta(phi,h,K,v)
    errors.extend([np.max(abs(P-Pback)),np.max(abs(j-jback))])
    # Matched h and phi alone do not match normal momenta.
    Kbad=.04*h;vbad=np.zeros(5)
    Pbad,jbad=jordan_momenta(phi,h,Kbad,vbad)
    # No-source inverse is exactly zero for the original nondegenerate branch.
    Kzero,vzero=invert_momenta(phi,h,np.zeros((3,3)),np.zeros(5))
    assert np.max(abs(Kzero))==0 and np.max(abs(vzero))==0
    assert max(errors)<2e-13 and np.linalg.norm(surfacejE)>1e-3
    assert np.linalg.norm(Pbad)>1e-2 and np.linalg.norm(jbad)>1e-3
    return dict(original_vacuum_phi=phi.tolist(),original_F=F,diagnostic_surface_tension=tau,
        outward_extrinsic_trace_sum=float(np.einsum('ab,ab',hi,K)),
        outward_scalar_derivative_sum=v.tolist(),
        transformed_scalar_surface_source=surfacejE.tolist(),
        omitted_transformed_scalar_source_norm=float(np.linalg.norm(surfacejE)),
        source_matching_error=float(max(errors)),
        continuous_fields_only_metric_momentum_defect=float(np.linalg.norm(Pbad)),
        continuous_fields_only_scalar_momentum_defect=float(np.linalg.norm(jbad)),
        no_source_matching_is_C1_in_adapted_coordinates=True,
        shell_is_new_diagnostic_not_a_postulated_constituent=True,
        global_shell_solution_not_constructed=True)

def run():
    deps=('research_note_358.md','research_note_574.md','research_note_583.md',
          'research_note_600.md','research_note_601.md','research_note_617.md',
          'joint_curved_quantum_source.py','joint_frame_hessian_matching.py')
    return dict(round=618,tests_run=3,failures=0,errors=0,
        momenta=momentum_check(),action=action_check(),interface=interface_check(),
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope=dict(original_F_target_and_boundary_Schur_jointly_fixed=True,
            nonminimal_metric_and_scalar_momenta_cannot_be_chosen_independently=True,
            classical_timelike_two_derivative_scalar_gravity_branch=True,
            gauge_fermion_surface_terms_null_corners_and_higher_orders_not_completed=True,
            Einstein_action_is_input_and_quantum_GR_still_open=True))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true')
    args=p.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps(result,ensure_ascii=False))

