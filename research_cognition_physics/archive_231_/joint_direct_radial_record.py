"""563: a direct radial record without angular alignment in the unchanged source.

The all-source impulse energy identity and the regular-source finite-time
existence certificate are separate. New RP coupling and a compensated terminal
quadrature are explicit instruments, not autonomous matter-derived devices.
"""
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import numpy as np
import joint_matter_energy_moment_control as poly
import joint_gauge_link_reference as link

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_direct_radial_record_results.json'
spec=importlib.util.spec_from_file_location('radial563',HERE/'round563_drafts/radial_pointer_source_probe.py')
probe=importlib.util.module_from_spec(spec);spec.loader.exec_module(probe)


def evaluate(P,q):
    return sum(c*math.prod(x**n for x,n in zip(q,powers)) for powers,c in P.items())


def radial_source():
    # Independent representation in (r1,r2,z,s1,s2), including the Haar law of z.
    m=probe.source.build();L,u,_=probe.source.old.constants()
    a,b,k,eta,omega=m['a'],m['b'],m['k'],m['eta'],m['omega']
    r1,r2,z,s1,s2=[poly.var(i,5) for i in range(5)]
    one=poly.const(1,5);r12=poly.mul(r1,r1);r22=poly.mul(r2,r2)
    S=poly.add(r12,r22);c=poly.mul(poly.mul(r1,r2),z)
    F=poly.scale(poly.add(S,poly.scale(c,-2)),k/2)
    R=poly.scale(poly.mul(poly.add(r1,poly.scale(r2,-1)),poly.add(r1,poly.scale(r2,-1))),.5)
    W={}
    for rr,ss in ((r12,s1),(r22,s2)):
        d=[poly.add(rr,poly.const(-u[0],5)),poly.add(poly.mul(ss,ss),poly.const(-u[1],5))]
        W=poly.add(W,poly.scale(poly.quadratic(d,L),.25))
    hlocal=poly.add(poly.const(10*a*omega,5),
        poly.scale(poly.add(S,poly.mul(s1,s1),poly.mul(s2,s2)),-a*omega**2),W,F)
    P=poly.add(one,poly.scale(poly.add(F,poly.const(-4*k*m['variance'],5)),1j*eta))
    HP=poly.add(poly.mul(hlocal,P),poly.scale(poly.add(poly.const(-8*a*k,5),
        poly.scale(F,4*a*omega),poly.scale(c,-3*b*k/4)),1j*eta))
    inv1={(-1,0,0,0,0):1.};inv2={(0,-1,0,0,0):1.}
    deltaR=poly.add(poly.const(8,5),poly.scale(poly.mul(r2,inv1),-3),poly.scale(poly.mul(r1,inv2),-3))
    def J(Q):
        radial_derivative=poly.mul(poly.add(r1,poly.scale(r2,-1)),
            poly.add(poly.deriv(Q,0),poly.scale(poly.deriv(Q,1),-1)))
        return poly.scale(poly.add(poly.mul(poly.add(deltaR,poly.scale(R,-4*omega)),Q),
            poly.scale(radial_derivative,2)),1j*a)
    def expect(Q):
        total=0j;t=m['variance']
        for powers,cq in Q.items():
            nr,ns,nz,n1,n2=powers
            if nz%2:continue
            j=nz//2
            angular=math.factorial(2*j)/(4**j*math.factorial(j)*math.factorial(j+1))
            assert nr>-4 and ns>-4
            radial=(2*t)**((nr+ns)/2)*math.gamma(2+nr/2)*math.gamma(2+ns/2)
            total+=cq*radial*angular*poly.normal_moment(n1,t)*poly.normal_moment(n2,t)
        return total
    assert abs(expect(poly.mul(poly.conj(P),P))-m['norm'])<1e-12
    exact={}
    for name,Q in [('Hpsi',HP),('J_R_psi',J(P)),('J_R_Hpsi',J(HP)),
                   ('Gamma_R_psi',poly.scale(poly.mul(R,P),4*a)),
                   ('Gamma_R_Hpsi',poly.scale(poly.mul(R,HP),4*a))]:
        squared=poly.real(expect(poly.mul(poly.conj(Q),Q)))/m['norm']
        assert squared>0
        exact[name]=math.sqrt(squared)
    return dict(model=m,P=P,HP=HP,J=J,R=R,expect=expect,exact_norms=exact)


def source_crosscheck(radial,certificate):
    rng=np.random.default_rng(563);model=radial['model'];HP=model['H'](model['P']);worst=0.
    for _ in range(24):
        q=rng.normal(size=10);r1=np.linalg.norm(q[:4]);r2=np.linalg.norm(q[5:9]);z=q[:4]@q[5:9]/(r1*r2)
        reduced=np.array([r1,r2,z,q[4],q[9]])
        a=evaluate(HP,q);b=evaluate(radial['HP'],reduced)
        worst=max(worst,abs(a-b)/max(1.,abs(a)))
    assert worst<1e-12
    rows={}
    for name,actual in radial['exact_norms'].items():
        upper=certificate['full_source_norm_upper'][name]
        assert actual<upper
        rows[name]=dict(independent_Laurent_Haar_Gaussian_norm=actual,Cartesian_triangle_upper=upper)
    return dict(polynomial_H_change_of_variables_residual=worst,norms=rows)


def differential_check(radial):
    rng=np.random.default_rng(2563);model=radial['model'];a=model['a'];b=model['b'];k=model['k'];eta=model['eta']
    L,u,_=probe.source.old.constants();theta=.37;errors=[];cross=0.
    step_ratios=(.04,.02,.01);step_errors=[[] for _ in step_ratios]
    def wave(q,U):
        r1=np.linalg.norm(q[:4]);r2=np.linalg.norm(q[5:9]);R=(r1-r2)**2/2
        F=k*np.sum((q[:4]-U@q[5:9])**2)/2
        return np.exp(-q@q/2-1j*theta*R)*(1+1j*eta*(F-1.6))
    for radius in (.025,.1,.7,1.3):
        for _ in range(3):
            q=.5*rng.normal(size=10);q[:4]*=radius/np.linalg.norm(q[:4])
            quat=rng.normal(size=4);quat/=np.linalg.norm(quat);U=link.left(quat)
            r1=np.linalg.norm(q[:4]);r2=np.linalg.norm(q[5:9]);z=q[:4]@U@q[5:9]/(r1*r2)
            R=(r1-r2)**2/2;reduced=np.array([r1,r2,z,q[4],q[9]])
            P=evaluate(radial['P'],reduced)
            target=np.exp(-q@q/2-1j*theta*R)*(evaluate(radial['HP'],reduced)
                +theta*evaluate(radial['J'](radial['P']),reduced)+theta**2*4*a*R*P)
            n1=q[:4]/r1;n2=q[5:9]/r2
            gradR=np.r_[(r1-r2)*n1,0.,-(r1-r2)*n2,0.]
            gradF=np.r_[k*(q[:4]-U@q[5:9]),0.,k*(q[5:9]-U.T@q[:4]),0.]
            cross=max(cross,abs(gradR@gradF-2*k*R*(1+z)))
            assert abs(gradR@gradR-4*R)<1e-12
            psi=wave(q,U)
            W=0.
            for o,r in ((0,r1),(5,r2)):
                d=np.array([r*r,q[o+4]**2])-u;W+=d@L@d/4
            F=k*np.sum((q[:4]-U@q[5:9])**2)/2
            # Fourth-order stencil over three step sizes: avoid subtractive
            # cancellation from the discarded 1e-4*r second-order stencil.
            for index,ratio in enumerate(step_ratios):
                h=ratio*min(r1,r2,1.);lap=glap=0j
                for i in range(10):
                    shift=np.zeros(10);shift[i]=h
                    lap+=(-wave(q+2*shift,U)+16*wave(q+shift,U)-30*psi
                          +16*wave(q-shift,U)-wave(q-2*shift,U))/(12*h*h)
                hg=ratio
                for j in link.J:
                    samples=[]
                    for multiple in (-2,-1,1,2):
                        angle=multiple*hg
                        rotation=math.cos(angle/2)*np.eye(4)+2*math.sin(angle/2)*j
                        samples.append(wave(q,rotation@U))
                    glap+=(-samples[0]+16*samples[1]-30*psi+16*samples[2]-samples[3])/(12*hg*hg)
                actual=-a*lap-b*glap+(W+F)*psi
                step_errors[index].append(float(abs(actual-target)/max(1.,abs(target))))
            errors.append(step_errors[-1][-1])
    assert max(errors)<5e-6 and cross<2e-14
    maxima=[max(row) for row in step_errors]
    assert maxima[1]<maxima[0]/5 and maxima[2]<maxima[1]/5
    return dict(full_ten_coordinate_and_three_group_phase_stencil_residual=max(errors),
        radii_tested=[.025,.1,.7,1.3],cross_gradient_residual=cross,
        fourth_order_step_ratios=list(step_ratios),step_maximum_residuals=maxima)


def threshold(n):
    x,w=np.polynomial.legendre.leggauss(n)
    r=4*(x+1);weight=4*w*2*r**3*np.exp(-r*r)
    r1=r[:,None];r2=r[None,:];wt=weight[:,None]*weight[None,:]
    S=r1*r1+r2*r2;R=(r1-r2)**2/2;k=.8;eta=.3;a=.5;N=1.1152
    f2=k*k/4*((S-4)**2+r1*r1*r2*r2)
    z=(R-.4)/.5
    G=np.array([.5*(1+math.erf(t/math.sqrt(2))) for t in z.ravel()]).reshape(R.shape)
    Gprime=np.exp(-z*z/2)/(.5*math.sqrt(2*math.pi))
    prob=float(np.sum(wt*(1+eta*eta*f2)*G))/N
    current=4*a*k*eta/N*float(np.sum(wt*R*Gprime))
    # Full-H density derivative after the original Gaussian angular average.
    schrodinger=4*a*k*eta/N*float(np.sum(wt*(S-4)*G))
    mean=float(np.sum(wt*R*(1+eta*eta*f2)))/N
    return dict(instantaneous_event_probability=prob,slope_current=current,
        slope_full_H_density=schrodinger,source_R_mean=mean,radial_integration_cutoff=8.)


def impulse_check(certificate):
    eta=.3;k=.8;N=1+2*eta*eta*k*k
    closed=(2-9*math.pi/16+eta*eta*k*k/4*(24-441*math.pi/64))/N
    assert abs(closed-certificate['source_radial_mean'])<2e-14
    q,w=np.polynomial.hermite.hermgauss(28);w/=math.sqrt(math.pi)
    worst=0.
    for R,nu in ((0.,.5),(.2391764717569,.5),(2.,.2)):
        noise=math.sqrt(2)*nu*q;y=R+noise
        # Integral of the squared weak gradient multiplier of the Kraus function.
        density_score=noise/(2*nu*nu)
        values=np.array([w@y,w@(y*y),w@density_score,w@(density_score*density_score)])
        target=np.array([R,R*R+nu*nu,0.,1/(4*nu*nu)])
        worst=max(worst,float(np.max(abs(values-target))))
    assert worst<1e-12
    return dict(closed_radial_mean=closed,kernel_moment_residual=worst,
        weak_gradient_energy_identity='Delta E_matter = a <R> / nu^2',
        nonselective_impulse_electric_energy_change=0.,
        conditional_or_finite_time_zero_electric_change_not_claimed=True)


def budget_check(certificate):
    model=probe.source.build();c=model['c'];E=certificate['same_source_E1']
    k=model['k'];a=model['a'];ell=c['ell'];v=c['v'];Ar=c['uh'];Br2=2/(v*ell)
    D0=2*Ar*Ar;D1=2*Br2;tau=certificate['sufficient_duration'];M=certificate['pointer_mass']
    rng=np.random.default_rng(3563);margin=[]
    L,u,_=probe.source.old.constants()
    for _ in range(80):
        q=rng.normal(size=10);r1=np.linalg.norm(q[:4]);r2=np.linalg.norm(q[5:9]);R=(r1-r2)**2/2
        W=0.
        for o,r in ((0,r1),(5,r2)):
            d=np.array([r*r,q[o+4]**2])-u;W+=d@L@d/4
        assert R<=Ar+math.sqrt(Br2*W)+1e-12
        assert R*R<=D0+D1*W+1e-12
        margin.append(D0+D1*W-R*R)
    assert abs(1-M*D1/(2*tau*tau)-.5)<1e-14
    mu1=math.sqrt(2/math.pi);mu3=2*mu1
    # R<=H/k, unlike the old F<=H; the k factors must remain here.
    energy_change=3*E*mu1/(k*tau)+2*(E/k+Ar)/(k*tau*tau)+Br2*mu3/(k*tau**3)
    mean_bound=min(E/k,Ar+math.sqrt(Br2*E))
    return dict(R_bound_constants=dict(Ar=Ar,Br_squared=Br2,D0=D0,D1=D1),
        sampled_square_bound_minimum_margin=float(min(margin)),
        source_R_mean_upper=mean_bound,
        impulse_source_energy_increase_upper=a*mean_bound/.5**2,
        impulse_second_record_moment_upper=D0+D1*E+.5**2,
        finite_time_energy_change_sufficient_upper=energy_change,
        whole_apparatus_remaining_W_coercivity=.5,
        total_H_lower=-D0/(2*D1),pointer_initial_energy=1/(2*M),
        fixed_derivative_readout_noise_cost='nu_R=epsilon^2*sigma_D; Delta E=a<R>/(epsilon^4*sigma_D^2)',
        switch_on_mean_work=0.,switch_off_work_equals_source_energy_change=True)


def run():
    certificate=probe.run();assert certificate==json.loads(probe.TARGET.read_text('utf8'))
    radial=radial_source();norms=source_crosscheck(radial,certificate)
    differential=differential_check(radial);impulse=impulse_check(certificate)
    row=threshold(160);other=threshold(220)
    residual=max(abs(row[k]-other[k]) for k in row)
    assert residual<1e-11
    assert abs(row['slope_current']-row['slope_full_H_density'])<1e-11
    assert abs(row['source_R_mean']-certificate['source_radial_mean'])<1e-11
    assert row['slope_current']>certificate['actual_threshold_derivative_positive_lower']
    budget=budget_check(certificate)
    checks=['same_full_source_independent_radial_Haar_representation',
        'inverse_radius_L2_norms_and_independent_exact_moments',
        'full_H_phase_conjugation_including_near_axis_points',
        'actual_impulse_record_and_nonselective_energy_identity',
        'bounded_radial_threshold_current_and_full_H_density_response',
        'finite_time_positive_record_and_common_source_pointer_budget']
    deps=('round563_drafts/radial_pointer_source_probe.py','round563_drafts/radial_pointer_source_probe_results.json',
        'round561_drafts/short_time_source_probe.py','joint_gauge_link_reference.py',
        'joint_matter_energy_moment_control.py','joint_singlet_common_mass_rg_results.json')
    return dict(round=563,tests_run=len(checks),failures=0,errors=0,checks=checks,
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        finite_time_certificate=certificate,independent_source_norms=norms,
        full_H_radial_conjugation=differential,impulse_record=impulse,
        actual_radial_threshold=row,threshold_quadrature_difference=residual,common_budget=budget,
        scope=dict(no_angular_alignment_or_new_source_H_required=True,
            same_560_regular_nonstationary_Gauss_source_used=True,
            new_RP_instrument_and_compensated_terminal_reader_are_inputs=True,
            only_impulse_nonselective_electric_energy_is_unchanged=True,
            weak_Cartesian_derivatives_not_flat_half_line_momentum=True,
            finite_duration_is_conservative_existence_not_optimal_cost=True,
            radial_edge_response_not_four_coordinates_or_continuum=True,
            no_three_dimensions_GR_or_unified_completion=True,
            floating_evaluations_not_machine_interval_certificates=True))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(round=563,tests=result['tests_run'],norms=result['independent_source_norms'],
        differential=result['full_H_radial_conjugation'],threshold=result['actual_radial_threshold'],
        budget=result['common_budget']),ensure_ascii=False))
