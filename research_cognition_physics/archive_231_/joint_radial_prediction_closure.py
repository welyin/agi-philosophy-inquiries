"""564: same radial density and total mean energy, different radial futures.

Full Gauss reduction is a polar-coordinate identity, not emergent spacetime.
All dynamics and pointer inputs retain the scope of 557--563.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import numpy as np
import joint_radial_link_alignment as old
import joint_gauge_link_reference as link
import joint_finite_time_gauge_probe as inherited

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_radial_prediction_closure_results.json'
A=old.A; B=old.B; K=old.K


def radial_grid(n):
    x,w=old.legendre(n);r=1+.4*x
    f=np.exp(-1/(1-x*x));weights=.4*w*r**3*f*f;weights/=weights.sum()
    lp=-2*x/(.4*(1-x*x)**2)
    lpp=-2/(.4**2*(1-x*x)**2)-8*x*x/(.4**2*(1-x*x)**3)
    return r,weights,lp,lpp


def class_grid(n):
    x,w=old.legendre(n);chi=np.pi*(x+1)/2
    return np.cos(chi),w*np.sin(chi)**2  # normalized Haar class measure after dchi


def basis_checks(n=320):
    z,w=class_grid(n);basis=[np.ones_like(z),2*z]
    for j in range(2,9):basis.append(2*z*basis[-1]-basis[-2])
    basis=np.array(basis);gram=(basis*w)@basis.T
    multiplication=(basis*(w*z))@basis.T
    target=np.diag(np.full(8,.5),1)+np.diag(np.full(8,.5),-1)
    residual=max(float(np.max(abs(gram-np.eye(9)))),float(np.max(abs(multiplication-target))))
    assert residual<2e-13
    # Differentiate independent polynomial coefficient recurrences.
    polys=[np.array([1.]),np.array([0.,2.])]
    for j in range(2,9):
        q=np.polynomial.polynomial.polymul([0.,2.],polys[-1])
        polys.append(np.polynomial.polynomial.polysub(q,polys[-2]))
    eigen_error=0.
    for j,p in enumerate(polys):
        d=np.polynomial.polynomial.polyder(p);dd=np.polynomial.polynomial.polyder(d)
        actual=-(1-z*z)*np.polynomial.polynomial.polyval(z,dd)+3*z*np.polynomial.polynomial.polyval(z,d)
        eigen_error=max(eigen_error,float(np.max(abs(actual-j*(j+2)*basis[j]))))
    assert eigen_error<2e-10
    return dict(orthogonality_and_z_matrix_residual=residual,
        class_Casimir_eigenvalue_residual=eigen_error,first_eigenvalues=[j*(j+2) for j in range(5)],
        constant_channel_to_first_channel_z_matrix_element=float(multiplication[0,1]))


def full_reduction_check():
    rng=np.random.default_rng(4564);L,u,_=inherited.constants()
    def bump(r):
        x=(r-1)/.4
        if abs(x)>=1:return 0.
        return math.exp(-1/(1-x*x))
    def radial_lap(r):
        x=(r-1)/.4
        lp=-2*x/(.4*(1-x*x)**2)
        lpp=-2/(.4**2*(1-x*x)**2)-8*x*x/(.4**2*(1-x*x)**3)
        return lpp+lp*lp+3*lp/r
    def relative(q,U):
        n1=q[:4]/np.linalg.norm(q[:4]);n2=q[5:9]/np.linalg.norm(q[5:9])
        conj=n1.copy();conj[1:]*=-1
        return link.left(conj)@U@n2
    def alpha(w):
        return 1+.3*w[0]+.17j*w[1]+.12*w[2]*w[3]+.08*w[0]**2
    def sphere_lap(w):
        return -3*(.3*w[0]+.17j*w[1])-8*.12*w[2]*w[3]+.08*(2-8*w[0]**2)
    def wave(q,U):
        r1=np.linalg.norm(q[:4]);r2=np.linalg.norm(q[5:9])
        return bump(r1)*bump(r2)*np.exp(-(q[4]**2+q[9]**2)/2)*alpha(relative(q,U))
    errors=[[],[],[]];gauge=0.
    for _ in range(12):
        q=.5*rng.normal(size=10)
        for sl in (slice(0,4),slice(5,9)):q[sl]*=rng.uniform(.9,1.1)/np.linalg.norm(q[sl])
        quat=rng.normal(size=(3,4));quat/=np.linalg.norm(quat,axis=1)[:,None]
        U,G,H=map(link.left,quat)
        moved=q.copy();moved[:4]=G@q[:4];moved[5:9]=H@q[5:9]
        gauge=max(gauge,float(np.max(abs(relative(q,U)-relative(moved,G@U@H.T)))))
        r1=np.linalg.norm(q[:4]);r2=np.linalg.norm(q[5:9]);w=relative(q,U)
        W=0.
        for o,r in ((0,r1),(5,r2)):
            d=np.array([r*r,q[o+4]**2])-u;W+=d@L@d/4
        F=K*np.sum((q[:4]-U@q[5:9])**2)/2
        base=-A*(radial_lap(r1)+radial_lap(r2))+A*(2-q[4]**2-q[9]**2)+W+F
        angular_coefficient=A*(1/r1**2+1/r2**2)+B/4
        pref=bump(r1)*bump(r2)*np.exp(-(q[4]**2+q[9]**2)/2)
        target=pref*(base*alpha(w)-angular_coefficient*sphere_lap(w))
        psi=wave(q,U)
        for index,h in enumerate((.008,.004,.002)):
            lap=0j;glap=0j
            for i in range(10):
                d=np.zeros(10);d[i]=h
                lap+=(-wave(q+2*d,U)+16*wave(q+d,U)-30*psi+16*wave(q-d,U)-wave(q-2*d,U))/(12*h*h)
            for j in link.J:
                vals=[]
                for n in (-2,-1,1,2):
                    rot=math.cos(n*h/2)*np.eye(4)+2*math.sin(n*h/2)*j
                    vals.append(wave(q,rot@U))
                glap+=(-vals[0]+16*vals[1]-30*psi+16*vals[2]-vals[3])/(12*h*h)
            actual=-A*lap-B*glap+(W+F)*psi
            errors[index].append(float(abs(actual-target)/max(1.,abs(target))))
    maxima=[max(row) for row in errors]
    assert gauge<2e-14
    assert maxima[-1]<1e-7
    assert maxima[1]<maxima[0]/8 and maxima[2]<maxima[1]/8
    return dict(full_Gauss_relative_quaternion_residual=gauge,
        noncentral_full_ten_plus_three_coordinate_H_steps=[.008,.004,.002],
        reduced_H_step_residuals=maxima)


def equal_energy_sources(n=160,nz=320):
    r,w,lp,lpp=radial_grid(n);r1=r[:,None];r2=r[None,:];wt=w[:,None]*w[None,:]
    delta=r1-r2;R=delta*delta/2;Ac=A*(1/r1**2+1/r2**2)+B/4
    Abar=float(np.sum(wt*Ac));c=K*float(np.sum(wt*r1*r2))
    theta=c/(3*Abar)  # alpha=(e0+theta e1)/sqrt(1+theta^2), an internal amplitude
    mean_z=theta/(1+theta*theta);mean_casimir=3*theta*theta/(1+theta*theta)
    z,wz=class_grid(nz);alpha=(1+2*theta*z)/math.sqrt(1+theta*theta)
    Kalpha=6*theta*z/math.sqrt(1+theta*theta)
    numerical_norm=float(wz@(alpha*alpha));numerical_z=float(wz@(z*alpha*alpha))
    numerical_casimir=float(wz@(alpha*Kalpha))
    assert max(abs(numerical_norm-1),abs(numerical_z-mean_z),abs(numerical_casimir-mean_casimir))<1e-13
    assert abs(Abar*mean_casimir-c*mean_z)<1e-14
    L,u,_=inherited.constants()
    Wavg=np.zeros_like(R)
    # Singlet Gaussian <s^2>=1/2, <s^4>=3/4, in the original complete W.
    for rr in (r1,r2):
        d=rr*rr-u[0]
        Wavg=Wavg+.25*(L[0,0]*d*d+2*L[0,1]*d*(.5-u[1])+L[1,1]*(.75-u[1]+u[1]**2))
    radial_lap=lpp+lp*lp+3*lp/r
    Hbase=-A*(radial_lap[:,None]+radial_lap[None,:])+A+Wavg+K*(r1*r1+r2*r2)/2
    EH0=float(np.sum(wt*Hbase))
    EHtheta=float(np.sum(wt*(Hbase*numerical_norm+Ac*numerical_casimir-K*r1*r2*numerical_z)))
    direct_form=2*A*float(w@(lp*lp))+A+float(np.sum(wt*Wavg))+K*float(w@(r*r))
    assert abs(EH0-EHtheta)<2e-12 and abs(EH0-direct_form)<2e-11
    sigma=.5;y0=.4
    standardized=(R-y0)/sigma
    G=np.array([.5*(1+math.erf(t/math.sqrt(2))) for t in standardized.ravel()]).reshape(R.shape)
    Gp=np.exp(-standardized**2/2)/(sigma*math.sqrt(2*math.pi))
    Gpp=-(R-y0)*Gp/sigma**2
    lapR=8-3*(r2/r1+r1/r2)
    HGminusGH=-A*(4*R*Gpp+Gp*lapR+2*Gp*delta*(lp[:,None]-lp[None,:]))
    # Independent full H form identity: r'' = -2 Re <H psi, [H,G] psi>.
    second0=-2*float(np.sum(wt*Hbase*HGminusGH))
    second_theta=-2*float(np.sum(wt*(Hbase+Ac*numerical_casimir-K*r1*r2*numerical_z)*HGminusGH))
    grad_G_grad_A=-2*A*Gp*delta*(r1**-3-r2**-3)
    assert np.min(grad_G_grad_A)>-1e-15
    I=float(np.sum(wt*R*Gp));J=float(np.sum(wt*grad_G_grad_A))
    force_difference=-2*A*J*mean_casimir-4*A*K*I*mean_z
    assert force_difference<0
    assert abs(second_theta-second0-force_difference)<1e-11
    probability=float(np.sum(wt*G))
    # The n=0 source leaves that angular channel at first order in state vector.
    leakage_coefficient=K*K/4*float(np.sum(wt*r1*r1*r2*r2))
    return dict(Abar=Abar,mean_link_product_times_k=c,theta=theta,
        angle_mean=mean_z,angle_Casimir=mean_casimir,
        angular_norm=numerical_norm,
        energy0_full_H=EH0,energy_theta_full_H=EHtheta,energy0_positive_form=direct_form,
        same_complete_radial_singlet_density_by_factorization=True,
        same_mean_full_H_not_just_common_upper_bound=True,
        same_initial_radial_current_zero=True,initial_actual_event_probability=probability,
        event_second_derivative0=second0,event_second_derivative_theta=second_theta,
        event_second_difference_full_H=second_theta-second0,
        event_second_difference_force=force_difference,
        radial_RGprime_mean=I,radial_gradient_coefficient=J,
        short_time_probability_difference_coefficient=force_difference/2,
        n0_leakage_probability_t2_coefficient=leakage_coefficient,
        higher_energy_moments_not_claimed_equal=True)


def run():
    basis=basis_checks();reduction=full_reduction_check()
    pair=equal_energy_sources();other=equal_energy_sources(240,480)
    numerical_fields=[key for key,value in pair.items() if type(value) is float]
    residual=max(abs(pair[key]-other[key]) for key in numerical_fields)
    assert residual<2e-10
    checks=['complete_Gauss_H_reduction_including_noncentral_sources',
        'class_Haar_Casimir_and_multiplication_matrix',
        'identical_radial_density_and_identical_complete_mean_energy_sources',
        'actual_bounded_radial_event_second_derivative_full_H_and_force_identity',
        'independent_quadrature_and_noninvariant_constant_angular_channel']
    names=('joint_radial_link_alignment.py','joint_radial_link_alignment_results.json',
           'joint_gauge_link_reference.py','joint_finite_time_gauge_probe.py',
           'joint_direct_radial_record.py','joint_direct_radial_record_results.json',
           'round564_drafts/radial_closure_source_probe.py','round564_drafts/radial_closure_source_probe_results.json')
    return dict(round=564,tests_run=len(checks),failures=0,errors=0,checks=checks,
        dependency_hashes={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in names},
        full_Gauss_reduction=reduction,class_basis=basis,equal_energy_pair=pair,
        independent_quadrature_max_difference=residual,
        scope=dict(only_radial_density_plus_total_mean_energy_prediction_contract_is_excluded=True,
            same_source_H_and_actual_Gaussian_threshold=True,
            complete_radial_angular_joint_state_retained_gives_exact_reduction=True,
            class_z_reduction_only_for_declared_class_function_sector=True,
            mean_z_or_separate_angular_marginal_not_proved_sufficient=True,
            finite_pointer_duration_by_strong_limit_existence_only=True,
            no_uniform_duration_or_optimal_cost_claim=True,
            graph_group_probe_and_boundary_domains_are_inputs=True,
            no_spacetime_dimension_gravity_or_unified_completion=True,
            floating_diagnostics_not_machine_interval_certificates=True))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(round=564,tests=result['tests_run'],reduction=result['full_Gauss_reduction'],
        pair=result['equal_energy_pair'],quadrature=result['independent_quadrature_max_difference']),ensure_ascii=False))
