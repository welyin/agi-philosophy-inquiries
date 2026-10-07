"""Local spatial-reference readout in a supplied finite harmonic lattice.

Geometry, Dirichlet sources, Gaussian preparation, and finite-noise final
instruments are inputs. This is neither a continuum limit nor an Einstein solution.
"""
import argparse
import itertools
import json
import math
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
TARGET=HERE/'local_reference_probe_results.json'


def lattice(n,d):
    sites=np.asarray(list(itertools.product(range(1,n+1),repeat=d)),int)
    lookup={tuple(x):k for k,x in enumerate(sites)}
    k=2*d*np.eye(len(sites))
    boundary=np.zeros((len(sites),d))
    for j,x in enumerate(sites):
        for axis in range(d):
            for step in (-1,1):
                y=x.copy(); y[axis]+=step
                if tuple(y) in lookup:
                    k[j,lookup[tuple(y)]]=-1
                else:
                    boundary[j]+=y
    assert np.max(abs(k@sites-boundary))<1e-14
    return sites,k


def spectral(mat,fun):
    ev,u=np.linalg.eigh(mat)
    return (u*fun(ev))@u.T


def setup(n,d,coupling=.5,mass2=1.):
    sites,k=lattice(n,d)
    N=len(sites); bmat=k+mass2*np.eye(N)
    a=np.block([[k,coupling*np.eye(N)],[coupling*np.eye(N),bmat]])
    w0=np.zeros_like(a); p0=np.zeros_like(a)
    w0[:N,:N]=spectral(k,lambda s:.5/np.sqrt(s))
    w0[N:,N:]=spectral(bmat,lambda s:.5/np.sqrt(s))
    p0[:N,:N]=spectral(k,lambda s:.5*np.sqrt(s))
    p0[N:,N:]=spectral(bmat,lambda s:.5*np.sqrt(s))
    gamma=np.block([[w0,np.zeros_like(a)],[np.zeros_like(a),p0]])
    omega=np.block([[np.zeros_like(a),np.eye(2*N)],[-np.eye(2*N),np.zeros_like(a)]])
    return sites,k,a,gamma,omega


def flow(a,t):
    vals,u=np.linalg.eigh(a)
    assert vals[0]>0
    w=np.sqrt(vals)
    c=(u*np.cos(t*w))@u.T
    s=(u*(np.sin(t*w)/w))@u.T
    # Stable 1-cos avoids cancellation in the forced mean.
    f=(u*(2*np.sin(t*w/2)**2/vals))@u.T
    return np.block([[c,s],[-a@s,c]]),f,s


def independent_series(a,force,t,terms=42):
    m=len(a); gen=np.block([[np.zeros_like(a),np.eye(m)],[-a,np.zeros_like(a)]])
    s=np.eye(2*m); term=s.copy()
    mean=np.zeros((2*m,force.shape[1]))
    v=np.vstack([np.zeros_like(force),-force])*t
    mean+=v
    for j in range(1,terms+1):
        term=term@(t*gen)/j; s+=term
        v=(t*gen)@v/(j+1); mean+=v
    return s,mean


def meter_check(gamma,mean,a,force,nu,N):
    """Explicit Gaussian meter: q_meter += q_probe, p_probe -= p_meter."""
    modes=3*N
    full=np.zeros((2*modes,2*modes))
    sysidx=np.r_[np.arange(2*N),modes+np.arange(2*N)]
    full[np.ix_(sysidx,sysidx)]=gamma
    full[2*N:3*N,2*N:3*N]=nu*nu*np.eye(N)
    full[5*N:6*N,5*N:6*N]=np.eye(N)/(4*nu*nu)
    r=np.eye(modes); r[2*N:3*N,N:2*N]=np.eye(N)
    shear=np.block([[r,np.zeros_like(r)],[np.zeros_like(r),np.linalg.inv(r).T]])
    after=shear@full@shear.T
    reduced=after[np.ix_(sysidx,sysidx)]
    expected=gamma.copy(); expected[3*N:4*N,3*N:4*N]+=np.eye(N)/(4*nu*nu)
    assert np.max(abs(reduced-expected))<1e-12
    output_cov=after[2*N:3*N,2*N:3*N]
    assert np.max(abs(output_cov-gamma[N:2*N,N:2*N]-nu*nu*np.eye(N)))<1e-12
    added_energy=.5*np.trace(reduced[2*N:,2*N:]-gamma[2*N:,2*N:])
    assert abs(added_energy-N/(8*nu*nu))<1e-11
    # Conditional reference covariance from this same joint meter model.
    refidx=np.r_[np.arange(N),modes+np.arange(N)]
    meteridx=np.arange(2*N,3*N)
    cross=after[np.ix_(refidx,meteridx)]
    conditional=after[np.ix_(refidx,refidx)]-cross@np.linalg.solve(output_cov,cross.T)
    om=np.block([[np.zeros((N,N)),np.eye(N)],[-np.eye(N),np.zeros((N,N))]])
    positivity=float(np.linalg.eigvalsh(conditional+.5j*om)[0])
    assert positivity>-1e-10
    return dict(one_axis_meter_energy_increment=float(added_energy),
                conditional_reference_uncertainty_min_eigenvalue=positivity,
                reference_record_cross_covariance_norm=float(np.linalg.norm(cross)))


def example(n,d,b=8192.,coupling=.5,mass2=1.,nu=.25,alpha=.01):
    x,k,a,gamma0,omega=setup(n,d,coupling,mass2)
    N=len(x); m=2*N
    ellmin=4*d*math.sin(math.pi/(2*(n+1)))**2
    assert coupling**2<ellmin*(ellmin+mass2)
    spec=np.linalg.eigvalsh(a)
    formula_min=ellmin+mass2/2-math.sqrt(mass2**2/4+coupling**2)
    assert abs(spec[0]-formula_min)<1e-12
    row_bound=4*d+mass2+abs(coupling)
    t=1/math.sqrt(n*row_bound)
    s,f,sine=flow(a,t)
    force=np.vstack([np.zeros((N,d)),coupling*b*x])
    mean=np.vstack([-f@force,-sine@force])
    gamma=s@gamma0@s.T
    series,series_mean=independent_series(a,force,t)
    series_error=float(np.max(abs(series-s)))
    mean_error=float(np.max(abs(series_mean-mean)))
    symplectic=float(np.max(abs(s@omega@s.T-omega)))
    assert series_error<1e-11 and mean_error<1e-8 and symplectic<1e-11
    assert np.linalg.eigvalsh(gamma+.5j*omega)[0]>-1e-11
    qmean,pmean=mean[:m],mean[m:]
    quantum_energy=lambda g:.5*np.trace(g[m:,m:]+a@g[:m,:m])
    energy_fluct_error=abs(quantum_energy(gamma)-quantum_energy(gamma0))
    energy_mean=.5*np.sum(pmean*pmean,axis=0)+.5*np.sum(qmean*(a@qmean),axis=0)+np.sum(force*qmean,axis=0)
    assert energy_fluct_error<1e-10 and np.max(abs(energy_mean))<1e-7
    gain=coupling*b*t*t/2
    calibrated=-qmean[N:]/gain
    bias=float(np.max(abs(calibrated-x)))
    z=t*math.sqrt(row_bound)
    beta=2*(math.cosh(z)-1)/z**2-1
    uniform_bias=2*math.cosh(1)-3
    assert bias<=n*beta+1e-12 and n*beta<=uniform_bias+1e-12 and uniform_bias<.125
    gap=min(float(np.max(abs(calibrated[i]-calibrated[j]))) for i in range(N) for j in range(i)) if N>1 else None
    if gap is not None: assert gap>=1-2*uniform_bias
    vmax_bound=1/(2*math.sqrt(ellmin))+t*t*math.sqrt(4*d+mass2)/2+nu*nu
    record_cov=gamma[N:2*N,N:2*N]+nu*nu*np.eye(N)
    vmax=float(np.max(np.diag(record_cov)))
    assert vmax<=vmax_bound
    needed_b=math.sqrt(512*vmax_bound*math.log(4*d*N/alpha))/(abs(coupling)*t*t)
    assert b>=needed_b
    noise_failure_bound=2*d*N*math.exp(-gain*gain/(128*vmax_bound))
    assert noise_failure_bound<=alpha/2
    delta=gain/4
    cutoff=gain*(n+uniform_bias)+math.sqrt(2*vmax_bound*math.log(4*d*N/alpha))
    cap=math.ceil(cutoff/delta)
    symbols=2*cap+2  # finite bins -cap,...,+cap and one overflow record
    bits=math.ceil(math.log2(symbols))
    # Invert each finite recorded bin, then integrate its exact Gaussian law.
    # Gaussian marginal integration does not assume independence across sites.
    coordinate_failures=[]
    for i in range(d):
        for j in range(N):
            good=[r for r in range(-cap,cap+1) if int(np.rint(-r*delta/gain))==x[j,i]]
            assert good and good==list(range(min(good),max(good)+1))
            lo=(min(good)-.5)*delta; hi=(max(good)+.5)*delta
            mu=qmean[N+j,i]; sigma=math.sqrt(record_cov[j,j])
            err=.5*math.erfc((mu-lo)/(math.sqrt(2)*sigma))+.5*math.erfc((hi-mu)/(math.sqrt(2)*sigma))
            coordinate_failures.append(err)
    meter=meter_check(gamma,mean,a,force,nu,N)
    # The reference gradient energy includes edges to fixed boundary values.
    gradient_per_axis=.5*b*b*(n+1)*n**(d-1)
    edge_energy=0.
    for axis in range(d):
        # For one reference axis, exactly (n+1)n^(d-1) oriented edges carry b.
        edge_energy+=gradient_per_axis
    fluct_energy=float(d*quantum_energy(gamma0))
    # b changes displacement but not the exact quadratic covariance.
    half_mean=np.vstack([-f@(force/2),-sine@(force/2)])
    assert np.max(abs(mean-2*half_mean))==0
    return dict(d=d,n=n,sites=N,reference_fields=d,probe_modes=d*N,
        inputs=dict(b=b,coupling=coupling,probe_mass_squared=mass2,read_noise_sd=nu,alpha=alpha,time=t),
        positive_joint_min_eigenvalue=float(spec[0]),
        checks=dict(symplectic_residual=symplectic,independent_series_matrix_error=series_error,
            independent_series_mean_error=mean_error,fluctuation_energy_residual=float(energy_fluct_error),
            mean_energy_residual=float(np.max(abs(energy_mean)))),
        coordinates=dict(max_calibration_bias=bias,finite_n_bias_bound=n*beta,
            uniform_bias_bound=uniform_bias,minimum_distinct_mean_linf_gap=gap,
            actual_record_variance_max=vmax,analytic_record_variance_bound=vmax_bound,
            sufficient_single_shot_b=needed_b,chosen_b=b,
            all_coordinate_noise_failure_upper=noise_failure_bound,
            all_record_overflow_failure_upper=alpha/2,
            all_coordinates_single_shot_failure_upper=noise_failure_bound+alpha/2,
            finite_bin_step=delta,finite_bin_cap=cap,symbols_per_local_axis=symbols,
            bits_per_local_axis=bits,total_record_bits=d*N*bits,
            exact_gaussian_marginal_failure_union_diagnostic=float(sum(coordinate_failures))),
        backreaction=dict(reference_mean_displacement_max=float(np.max(abs(qmean[:N]))),
            reference_mean_momentum_max=float(np.max(abs(pmean[:N]))),
            reference_covariance_change_norm=float(np.linalg.norm(gamma[:N,:N]-gamma0[:N,:N])),
            **meter),
        resources=dict(reference_gradient_energy=float(edge_energy),
            initial_fluctuation_energy=fluct_energy,
            initial_total_energy=float(edge_energy+fluct_energy),
            final_read_energy_increment=float(d*N/(8*nu*nu)),
            independent_source_preparations=1,full_read_rounds=1),
        first_local_mean_coordinate=[float(v) for v in calibrated[0]])


def negative_controls():
    n,d=3,2
    x,k,a,gamma0,omega=setup(n,d)
    N=len(x)
    # Single x-axis reference: transverse reflection preserves its entire
    # Hamiltonian, source, initial Gaussian law, and declared local read.
    perm=[int(np.flatnonzero(np.all(x==[z[0],n+1-z[1]],axis=1))[0]) for z in x]
    p=np.eye(N)[perm]
    assert np.max(abs(p@k-k@p))==0 and np.max(abs(p@x[:,0]-x[:,0]))==0
    max_mean=max_variance=0.
    for t in (.1,.7,2.):
        s,f,si=flow(a,t)
        j=np.r_[np.zeros(N),.5*4096*x[:,0]]
        means=-(f@j)[N:]
        gamma=s@gamma0@s.T
        variance=np.diag(gamma[N:2*N,N:2*N])+.25**2
        max_mean=max(max_mean,float(np.max(abs(means-p@means))))
        max_variance=max(max_variance,float(np.max(abs(variance-p@variance))))
    assert max_mean<1e-10 and max_variance<1e-12
    # One analytically identified unstable larger lattice, not a parameter scan.
    x6,k6,a6,g6,o6=setup(6,1)
    unstable=float(np.linalg.eigvalsh(a6)[0])
    assert unstable<0
    critical=(math.sqrt(2)-1)/2
    assert 10/49<critical  # pi^2<10 gives ell_min(n=6)<10/49
    # Decoupling removes the coordinate mean entirely.
    x0,k0,a0,g0,o0=setup(3,2,coupling=0)
    s0,f0,si0=flow(a0,.7)
    source0=np.zeros((len(a0),2))
    assert np.max(abs(f0@source0))==0
    return dict(single_axis_transverse_reflection=dict(n=n,d=d,
         example_distinct_sites=[[1,1],[1,3]],
         sampled_mean_difference=max_mean,sampled_variance_difference=max_variance,
         exact_argument='transverse permutation symmetry; one local final read only'),
        fixed_coupling_large_domain_failure=dict(n=6,d=1,coupling=.5,mass_squared=1.,
         joint_quadratic_min_eigenvalue=unstable,critical_laplacian_eigenvalue=critical),
        zero_coupling_coordinate_signal=0.)


def run():
    rows=[example(3,d) for d in (1,2,3)]
    controls=negative_controls()
    return dict(date='2026-09-30',diagnostic_groups=6,failures=0,errors=0,
        numbered_round_created=False,numbered_test_increment=0,
        diagnostics=dict(finite_lattice_joint_models=rows,negative_controls=controls),
        scope=dict(supplied_geometry=True,dirichlet_reference_boundary_input=True,
            specified_gaussian_initial_family=True,one_shot_finite_record_certificate=True,
            all_unknown_quantum_inputs=False,quantum_clock_control_included=False,
            continuum_or_einstein_backreaction_proved=False,three_dimensions_selected=False))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true')
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    answer=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf-8',newline='\n') as out:
            out.write(json.dumps(answer,ensure_ascii=False,indent=2)+'\n')
    if args.check:
        assert json.loads(TARGET.read_text('utf-8'))==json.loads(json.dumps(answer))
    print(json.dumps(answer,ensure_ascii=False,indent=2))
