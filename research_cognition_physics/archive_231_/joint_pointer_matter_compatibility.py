"""576: the autonomous pointer changes the common matter spectrum and amplitudes.
Same two-derivative Einstein-frame effective action, flat constant vacuum, tree level.
This is not a derivation of a quantum continuum limit or an experimental fit.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_autonomous_pointer_geometry as old

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_pointer_matter_compatibility_results.json'
M=old.geo.M02
MAT,U0,_=old.lat.scalar.parameters()
VAC=np.sqrt(U0)


def base_metric(x):
    f=M-np.dot(x,x)/6
    return np.eye(2)/f+np.outer(x,x)/(6*f*f)


def potential(x):
    f=M-np.dot(x,x)/6
    d=x*x-U0
    return d@MAT@d/(4*f*f)


def connection(x,g=1.,shift='none'):
    h,s=x;h0,s0=VAC
    a=np.array([g*s*h,0.],dtype=np.result_type(x,g))
    if shift=='normal':
        # f = g*s0*(T-T0) + (g/2)*(s-s0)*(T-T0), globally invariant.
        a-=np.array([g*s0*h+g*(s-s0)*h/2,g*(h*h-h0*h0)/4])
    elif shift=='polynomial':
        # An arbitrary smooth fibre coordinate change, not a new interaction.
        a-=np.array([.4*h-.31*s+.17,-.31*h+.14*s-.09])
    elif shift!='none':raise ValueError(shift)
    return a


def metric(x,g=1.,shift='none'):
    k=base_metric(x);a=connection(x,g,shift)
    return np.block([[k+np.outer(a,a),a[:,None]],[a[None,:],np.ones((1,1))]])


def metric_derivatives(x,g=1.,shift='none'):
    out=np.zeros((3,3,3))  # derivative coordinate, matrix row, matrix column
    for k in range(2):
        z=np.array(x,dtype=complex);z[k]+=1e-30j
        out[k]=metric(z,g,shift).imag/1e-30
    return out


def christoffel(x,g=1.,shift='none'):
    inverse=np.linalg.inv(metric(x,g,shift));dg=metric_derivatives(x,g,shift)
    gamma=np.zeros((3,3,3))
    for a in range(3):
        for i in range(3):
            for j in range(3):
                gamma[a,i,j]=sum(inverse[a,l]*(dg[i,j,l]+dg[j,i,l]-dg[l,i,j])/2 for l in range(3))
    return gamma


def generalized_modes(hess,kin):
    invT=np.linalg.inv(np.linalg.cholesky(kin).T)
    vals,rot=np.linalg.eigh(invT.T@hess@invT)
    return vals,invT@rot


def vacuum(g=1.,shift='none'):
    h,s=VAC;f=M-np.dot(VAC,VAC)/6
    kin=base_metric(VAC)
    hess=2*np.diag(VAC)@MAT@np.diag(VAC)/f**2
    masses,B=generalized_modes(hess,kin)
    if np.linalg.det(B)<0:B[:,0]*=-1
    a=connection(VAC,g,shift)
    E=np.zeros((3,3));E[:2,:2]=B;E[2,:2]=-a@B;E[2,2]=1
    full_hessian=np.zeros((3,3));full_hessian[:2,:2]=hess
    c=g*h/np.sqrt(np.linalg.det(kin))
    return dict(F=f,K=kin,H=hess,m2=masses,B=B,E=E,G=metric(VAC,g,shift),
                full_hessian=full_hessian,c=c,a=a)


def vacuum_spectrum_check():
    q=vacuum();h,s=VAC;phi=np.array([0.,h,0.,0.,s])
    G6,_=old.target_blocks(phi)
    H6=np.zeros((6,6));H6[np.ix_([1,4],[1,4])]=q['H']
    vals6,_=generalized_modes(H6,G6)
    assert np.min(np.linalg.eigvalsh(MAT))>0 and q['F']>0 and abs(potential(VAC))<1e-28
    assert max(abs(q['E'].T@q['G']@q['E']-np.eye(3)).flat)<1e-14
    assert np.max(abs(q['E'].T@q['full_hessian']@q['E']-np.diag([*q['m2'],0.])))<1e-14
    assert np.max(abs(vals6[:4]))<1e-14 and np.max(abs(vals6[4:]-q['m2']))<1e-14
    hessian_errors=[]
    for step in (2e-4,1e-4,5e-5):
        hessian_errors.append(float(np.max(abs(old.lat.scalar.hessian(potential,VAC,step)-q['H']))))
    assert hessian_errors[-1]<hessian_errors[0]/12 and hessian_errors[-1]<1e-7
    assert abs(np.linalg.det(q['K'])-M/q['F']**3)<1e-15
    # Higgs angular tangent directions remain orthogonal to the neutral fibre connection.
    tangent=np.eye(5)[[0,2,3]]
    assert np.max(abs(tangent@old.connection(phi)))==0
    return dict(vacuum=VAC.tolist(),F=q['F'],radial_metric=q['K'].tolist(),old_mass_squared=q['m2'].tolist(),
        full_six_scalar_unreduced_eigenvalues=vals6.tolist(),radial_plus_pointer_mass_squared=[*q['m2'].tolist(),0.],
        numerical_Hessian_errors=hessian_errors,extra_neutral_scalar_modes=1,
        three_Higgs_angular_zeros_not_extra_physical_scalars=True,
        flat_vacuum_not_inhomogeneous_CMC_source=True)


def momenta(m2,direction):
    ml2,mh2=m2;mh=np.sqrt(mh2);p=(mh2-ml2)/(2*mh);el=(mh2+ml2)/(2*mh)
    n=np.asarray(direction);n=n/np.linalg.norm(n)
    # All incoming, order light, heavy, pointer; Minkowski signature +---.
    return np.array([[-el,*(-p*n)],[mh,0.,0.,0.],[-p,*(p*n)]])


def direct_amplitude(q,g=1.,shift='none',direction=(0.,0.,1.)):
    p=momenta(q['m2'],direction);dot=(p*np.array([1.,-1.,-1.,-1.]))@p.T
    dg=metric_derivatives(VAC,g,shift)
    vertex=np.einsum('kij,ia,jb,kc->abc',dg,q['E'],q['E'],q['E'])
    # The partial potential third derivative with a vertical external pointer is zero.
    amplitude=-(vertex[0,1,2]*dot[0,1]+vertex[0,2,1]*dot[0,2]+vertex[1,2,0]*dot[1,2])
    return float(amplitude)


def covariant_amplitude(q,g=1.,shift='none'):
    gamma=christoffel(VAC,g,shift);H=q['full_hessian']
    # At a critical point, U_;chi,j,k has no partial derivative contribution.
    tensor=np.zeros((3,3))
    for j in range(3):
        for k in range(3):
            tensor[j,k]=-sum(gamma[l,2,j]*H[l,k]+gamma[l,2,k]*H[j,l] for l in range(3))
    return -float(q['E'][:,0]@tensor@q['E'][:,1])


def cubic_amplitude_check():
    q=vacuum();expected=-q['c']*np.diff(q['m2'])[0]/2
    direct=direct_amplitude(q);covariant=covariant_amplitude(q)
    errors=[abs(direct-expected),abs(covariant-expected)]
    assert max(errors)<1e-14
    qn=vacuum(shift='normal')
    assert np.max(abs(qn['a']))<1e-15
    dg=metric_derivatives(VAC,shift='normal')
    v=np.einsum('kij,ia,jb,kc->abc',dg,qn['E'],qn['E'],qn['E'])
    assert abs(v[0,2,1]-q['c']/2)<1e-14
    assert abs(v[1,2,0]+q['c']/2)<1e-14
    # Normal-coordinate Christoffels carry the same mass-difference coupling.
    E=qn['E'];frame_gamma=np.einsum('ai,ijk,jb,kc->abc',np.linalg.inv(E),
        christoffel(VAC,shift='normal'),E,E)
    assert abs(frame_gamma[1,2,0]+q['c']/2)<1e-14
    assert abs(frame_gamma[0,2,1]-q['c']/2)<1e-14
    return dict(curvature_coupling=q['c'],mass_gap_squared=float(np.diff(q['m2'])[0]),
        predicted_signed_amplitude=expected,raw_metric_vertex=direct,covariant_potential_vertex=covariant,
        independent_amplitude_errors=errors,normal_connection_cubic_coefficients=[v[0,2,1],v[1,2,0]])


def coordinate_invariance_check():
    rows=[];maxerr=0.
    for g in (0.,1.,-.7):
        for shift in ('none','normal','polynomial'):
            q=vacuum(g,shift);ref=-q['c']*np.diff(q['m2'])[0]/2
            amp=direct_amplitude(q,g,shift,(.3,-.4,.8));cov=covariant_amplitude(q,g,shift)
            maxerr=max(maxerr,abs(amp-ref),abs(cov-ref))
            rows.append(dict(g=g,fibre_coordinates=shift,amplitude=amp,covariant_amplitude=cov))
    assert maxerr<1e-14
    # No mass splitting closes this decay, but does not remove the extra field.
    q=vacuum();q['m2']=np.array([.1,.1])
    assert abs(direct_amplitude(q))<1e-14
    return dict(cases=rows,maximum_error=maxerr,degeneracy_closes_only_this_channel=True)


def decay_width_check():
    q=vacuum();ml,mh=np.sqrt(q['m2']);delta=q['m2'][1]-q['m2'][0]
    width=q['c']**2*delta**3/(64*np.pi*mh**3)
    # Independent two-body phase-space angular quadrature using the raw metric vertex.
    xs,weights=np.polynomial.legendre.leggauss(7);angles=np.arange(11)*2*np.pi/11
    integral=0.;spread=[];p=delta/(2*mh)
    for x,w in zip(xs,weights):
        for phi in angles:
            n=(np.sqrt(1-x*x)*np.cos(phi),np.sqrt(1-x*x)*np.sin(phi),x)
            amp=direct_amplitude(q,direction=n);spread.append(amp)
            integral+=w*(2*np.pi/len(angles))*amp**2
    integrated=p*integral/(32*np.pi**2*mh**2)
    assert abs(integrated-width)<1e-16 and max(spread)-min(spread)<1e-14
    assert width>0 and np.sqrt(q['m2'][1])>np.sqrt(q['m2'][0])
    return dict(light_mass=ml,heavy_mass=mh,pointer_mass=0.,tree_partial_width=width,
                direct_phase_space_width=integrated,width_to_heavy_mass=width/mh,
                raw_vertex_angular_spread=max(spread)-min(spread),
                model_units_only_no_pole_mass_or_experimental_matching=True)


def zero_flow_check():
    data=old.source(16);C=data['C'];psi=data['psi'];geo=old.geo
    pdot=sum(geo.derivative(psi**2*C[...,i],i) for i in range(3))
    x,y,_=np.moveaxis(data['q']['grid'],-1,0)
    h=data['q']['f']['h'];s=data['q']['f']['s'];h0,s0=VAC
    dh=.05*h0*np.cos(x+y);ds=-.06*s0*np.sin(y);ddh=-.05*h0*np.sin(x+y)
    analytic=ds*h*dh+2*s*dh**2+2*s*h*ddh
    flat=sum(geo.derivative(C[...,i],i) for i in range(3))
    assert np.max(abs(flat-analytic))<1e-13
    assert np.max(abs(pdot))>1e-3 and abs(float(np.mean(pdot)))<1e-15
    # dC = ds wedge dT; a zero one-form C on an open set would require dC=0 there.
    curl_xy=-ds*h*dh
    assert np.max(abs(curl_xy))>1e-5
    return dict(max_local_pointer_momentum_derivative=float(np.max(abs(pdot))),
        total_pointer_charge_derivative=float(geo.VOL*np.mean(pdot)),
        rms_pointer_momentum_derivative=float(np.sqrt(np.mean(pdot*pdot))),
        flat_divergence_independent_error=float(np.max(abs(flat-analytic))),
        max_pulled_back_connection_curvature=float(np.max(abs(curl_xy))),
        initial_Pchi_zero_not_invariant_local_constraint=True)


def massive_completion_check():
    q=vacuum();a=q['a'];G=q['G'];H=q['full_hessian']
    low_limit,_=generalized_modes(q['H'],q['K']+np.outer(a,a))
    rows=[];errors=[]
    for mu in (1.,4.,16.):
        bare=H.copy();bare[2,2]+=mu*mu
        bare_eigs,_=generalized_modes(bare,G)
        # Correlated positive periodic pinning: mu^2[1-cos(chi+f(phi))], df*=a*.
        v=np.r_[a,1.]
        matched=H+mu*mu*np.outer(v,v)
        matched_eigs,_=generalized_modes(matched,G)
        matched_error=float(np.max(abs(matched_eigs-np.r_[q['m2'],mu*mu])))
        assert matched_error<2e-12
        errors.append(float(np.max(abs(bare_eigs[:2]-low_limit))))
        # Schur complement of the same quadratic inverse propagator.
        z=-.03
        block=z*G-bare
        eliminated=block[:2,:2]-np.outer(block[:2,2],block[2,:2])/block[2,2]
        formula=z*q['K']-q['H']+z*mu*mu/(mu*mu-z)*np.outer(a,a)
        assert np.max(abs(eliminated-formula))<1e-13
        rows.append(dict(mu=mu,bare_chi_pin_mass_squared=bare_eigs.tolist(),
                         correlated_pin_mass_squared=matched_eigs.tolist(),
                         correlated_pin_error=matched_error))
    assert errors[-1]<errors[0]/100
    assert np.max(abs(low_limit-q['m2']))>1e-3
    return dict(cases=rows,bare_pin_heavy_limit_radial_masses_squared=low_limit.tolist(),
        heavy_limit_errors=errors,original_radial_masses_squared=q['m2'].tolist(),
        bare_pin_changes_original_radial_physics=True,
        correlated_pin_is_an_explicit_alternative_with_additional_potential=True,
        record_and_inhomogeneous_geometry_not_revalidated_for_pinned_models=True)


def run():
    checks=('vacuum_spectrum_check','cubic_amplitude_check','coordinate_invariance_check',
            'decay_width_check','zero_flow_check','massive_completion_check')
    evidence={name:globals()[name]() for name in checks}
    deps=('joint_autonomous_pointer_geometry.py','joint_curved_quantum_source.py',
          'joint_scalar_propagation_matching.py','research_round_575_checks.json')
    return dict(round=576,tests_run=len(checks),failures=0,errors=0,checks=list(checks),evidence=evidence,
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope='same explicitly extended two-derivative effective matter action, constant flat classical vacuum and tree amplitudes; extra neutral massless mode, nonzero heavy-to-light pointer channel, current-source truncation failure and qualified massive alternatives; not a quantum continuum proof, realistic decay fit, full apparatus or quantum gravity')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps({k:result[k] for k in ('round','tests_run','failures','errors')}))
