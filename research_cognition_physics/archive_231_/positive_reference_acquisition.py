"""Blank-probe coordinate acquisition under the fixed positive alignment rule.

Finite supplied lattice, fixed boundary sources, known Gaussian preparation,
and a final noisy instrument are inputs. Stable acquisition is not equilibration.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import numpy as np
import local_reference_probe_model as old

HERE=Path(__file__).resolve().parent
TARGET=HERE/"positive_reference_acquisition_results.json"


def build(n,d,b=8192.,m2=1.,r=.5):
    x,k,a0,gamma0,omega=old.setup(n,d,coupling=m2*r,mass2=m2)
    N=len(x)
    a=a0.copy(); a[:N,:N]+=m2*r*r*np.eye(N)
    lam=m2*r
    force=np.vstack([lam*r*b*x,lam*b*x])
    mu2=m2*(1+r*r)
    return x,k,a,gamma0,omega,force,mu2


def certificate(x,k,a,gamma0,omega,force,mu2,b,m2,r,nu=.25,alpha=.01):
    N,d=x.shape; n=int(np.max(x)); lam=m2*r; modes=2*N
    t=1/math.sqrt(n*(4*d+mu2))
    s,f,sine=old.flow(a,t)
    mean=np.vstack([-f@force,-sine@force])
    gamma=s@gamma0@s.T
    # An independent massive-mode expression uses the corrected force on BOTH fields.
    am=k+mu2*np.eye(N)
    sm,fm,sinem=old.flow(am,t)
    probe=-lam*b*fm@x
    probe_p=-lam*b*sinem@x
    direct=np.vstack([r*probe,probe,r*probe_p,probe_p])
    assert np.max(abs(mean-direct))<1e-8
    series,series_mean=old.independent_series(a,force,t)
    assert np.max(abs(s-series))<1e-11
    assert np.max(abs(mean-series_mean))<1e-8
    assert np.max(abs(s@omega@s.T-omega))<1e-11
    assert np.linalg.eigvalsh(gamma+.5j*omega)[0]>-1e-11
    gain=lam*b*t*t/2
    calibrated=-probe/gain
    beta1=2*math.cosh(1)-3
    z=t*math.sqrt(4*d+mu2)
    finite_bias=n*(2*(math.cosh(z)-1)/z**2-1)
    bias=float(np.max(abs(calibrated-x)))
    assert bias<=finite_bias+1e-12 and finite_bias<=beta1+1e-12
    assert beta1<.125
    ellmin=4*d*math.sin(math.pi/(2*(n+1)))**2
    vstar=1/(2*math.sqrt(ellmin))+t*t*math.sqrt(4*d+m2)/2+nu*nu
    record_cov=gamma[N:2*N,N:2*N]+nu*nu*np.eye(N)
    assert np.max(np.diag(record_cov))<=vstar
    bmin=math.sqrt(512*vstar*math.log(4*d*N/alpha))/(abs(lam)*t*t)
    assert b>=bmin
    noise_upper=2*d*N*math.exp(-gain*gain/(128*vstar))
    assert noise_upper<=alpha/2
    delta=abs(gain)/4
    cap=math.ceil((abs(gain)*(n+beta1)+math.sqrt(2*vstar*math.log(4*d*N/alpha)))/delta)
    symbols=2*cap+2
    # Direct integration of finite record bins supplements, not replaces, the bound.
    marginal_union=0.
    for i in range(d):
        for j in range(N):
            good=[z for z in range(-cap,cap+1) if int(np.rint(-z*delta/gain))==x[j,i]]
            assert good==list(range(min(good),max(good)+1))
            lo,hi=(min(good)-.5)*delta,(max(good)+.5)*delta
            sd=math.sqrt(record_cov[j,j]); center=probe[j,i]
            marginal_union+=.5*math.erfc((center-lo)/(math.sqrt(2)*sd))+.5*math.erfc((hi-center)/(math.sqrt(2)*sd))
    qmean,pmean=mean[:modes],mean[modes:]
    mean_energy=.5*np.sum(pmean*pmean,axis=0)+.5*np.sum(qmean*(a@qmean),axis=0)+np.sum(force*qmean,axis=0)
    scale=.5*np.sum(pmean*pmean,axis=0)+.5*np.sum(abs(qmean*(a@qmean)),axis=0)+np.sum(abs(force*qmean),axis=0)
    assert np.max(abs(mean_energy)/np.maximum(scale,1.))<1e-12
    fluct=lambda g:.5*float(np.trace(g[modes:,modes:]+a@g[:modes,:modes]))
    fluct_res=abs(fluct(gamma)-fluct(gamma0))
    assert fluct_res<1e-10
    egrad=d*b*b*(n+1)*n**(d-1)/2
    counterterm_mean=m2*r*r*b*b*float(np.sum(x*x))/2
    counterterm_variance=d*m2*r*r*float(np.trace(gamma0[:N,:N]))/2
    meter=old.meter_check(gamma,mean,a,force,nu,N)
    # Correct common reference mean is constant; the relative mean is NOT aligned.
    common_displacement=(qmean[:N]-r*qmean[N:])/math.sqrt(1+r*r)
    assert np.max(abs(common_displacement))<1e-8
    common_initial=b*x/math.sqrt(1+r*r)
    relative_initial=r*b*x/math.sqrt(1+r*r)
    relative_stationary=np.linalg.solve(am,k@relative_initial)
    relative=relative_initial+(r*qmean[:N]+qmean[N:])/math.sqrt(1+r*r)
    relative_p=(r*pmean[:N]+pmean[N:])/math.sqrt(1+r*r)
    displacement0=relative_initial-relative_stationary
    excess0=.5*np.sum(displacement0*(am@displacement0))
    displacement=relative-relative_stationary
    excess=.5*np.sum(relative_p*relative_p)+.5*np.sum(displacement*(am@displacement))
    assert excess0>0 and abs(excess-excess0)/excess0<1e-12
    # Keeping the old force after changing A is a different, incorrectly sourced model.
    wrong_force=force.copy(); wrong_force[:N]=0.
    wrong=-f@wrong_force
    wrong_probe_gap=float(np.max(abs(wrong[N:]-probe)))
    assert wrong_probe_gap>1e-4
    return dict(n=n,d=d,sites=N,inputs=dict(b=b,mass2=m2,r=r,coupling=lam,nu=nu,alpha=alpha,time=t),
        checks=dict(normal_mode_mean_error=float(np.max(abs(mean-direct))),
            independent_series_error=float(np.max(abs(s-series))),
            independent_mean_error=float(np.max(abs(mean-series_mean))),
            symplectic_residual=float(np.max(abs(s@omega@s.T-omega))),
            mean_energy_absolute_residual=float(np.max(abs(mean_energy))),
            mean_energy_scaled_residual=float(np.max(abs(mean_energy)/np.maximum(scale,1.))),
            fluctuation_energy_residual=fluct_res),
        coordinates=dict(max_bias=bias,bias_upper=finite_bias,sufficient_b=bmin,
            record_variance_max=float(np.max(np.diag(record_cov))),variance_upper=vstar,
            all_coordinate_failure_upper=noise_upper+alpha/2,
            local_axis_symbols=symbols,total_bits=d*N*math.ceil(math.log2(symbols)),
            exact_marginal_failure_union_diagnostic=float(marginal_union)),
        source_and_backreaction=dict(common_mean_displacement=float(np.max(abs(common_displacement))),
            reference_mean_displacement_max=float(np.max(abs(qmean[:N]))),
            probe_mean_magnitude_max=float(np.max(abs(probe))),
            old_force_probe_error=wrong_probe_gap,
            relative_initial_norm=float(np.linalg.norm(relative_initial)),
            relative_stationary_norm=float(np.linalg.norm(relative_stationary)),
            relative_at_read_norm=float(np.linalg.norm(relative)),
            relative_oscillation_excess_energy=float(excess0),
            excess_energy_relative_residual=float(abs(excess-excess0)/excess0),
            **meter),
        resources=dict(gradient_energy=egrad,counterterm_mean_energy=counterterm_mean,
            counterterm_fluctuation_energy=counterterm_variance,
            total_initial_energy=egrad+counterterm_mean+d*fluct(gamma0),
            final_read_energy_increment=d*N/(8*nu*nu),
            source_preparations=1,read_rounds=1))


def run():
    b,m2,r=8192.,1.,.5
    rows=[]
    for n,d in ((3,1),(3,2),(3,3),(6,1)):
        x,k,a,g0,om,j,mu2=build(n,d,b,m2,r)
        assert np.linalg.eigvalsh(a)[0]>0
        ell=4*d*math.sin(math.pi/(2*(n+1)))**2
        assert abs(np.linalg.eigvalsh(a)[0]-ell)<1e-12
        rows.append(certificate(x,k,a,g0,om,j,mu2,b,m2,r))
    # The old n=6 failure now has the SAME fixed coupling and a strictly positive lowest mode.
    stability=[]
    for n in (6,127):
        ell=4*math.sin(math.pi/(2*(n+1)))**2
        a=ell*np.eye(2)+m2*np.array([[r*r,r],[r,1.]])
        ev=np.linalg.eigvalsh(a)
        assert abs(ev[0]-ell)<1e-14
        stability.append(dict(n=n,d=1,min_eigenvalue=float(ev[0]),coupling=m2*r))
    return dict(date="2026-09-30",diagnostic_groups=6,failures=0,errors=0,
        numbered_round_created=False,numbered_test_increment=0,
        dependencies={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                      for name in ("local_reference_probe_model.py","reference_alignment_completion.py")},
        diagnostics=dict(blank_probe_examples=rows,fixed_rule_large_box_spectrum=stability),
        scope=dict(fixed_positive_rule_all_finite_box_stability=True,
            same_rule_blank_probe_coordinate_readout=True,known_gaussian_source_family=True,
            all_unknown_input_diamond_guarantee=False,permanent_relative_alignment_formed=False,
            size_independent_reading_budget=False,quantum_clock_control=False,
            dynamic_einstein_source_solved=False,dimension_selected=False,complete_unification=False))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write-results",action="store_true")
    p.add_argument("--check",action="store_true")
    args=p.parse_args(); result=run()
    if args.write_results:
        with TARGET.open("x",encoding="utf8",newline="\n") as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    if args.check:
        assert json.loads(TARGET.read_text("utf8"))==json.loads(json.dumps(result))
    print(json.dumps(result,ensure_ascii=False,indent=2))
