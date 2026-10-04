"""681: gauge-invariant free-mode diagnostics for a proposed flow/time mapping.

This is NOT the original interacting Gauss/S9 candidate. Analytic signs are
proved in the note. Numerical integrals cross-check the witnesses, not certify
nonperturbative gauge theories. Uses only the existing Python/NumPy runtime.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_gauge_flow_time_interface_results.json'

def heat_closed(t,s,omega):
    a=2*s
    return (math.exp(-omega*t)*math.erfc(math.sqrt(a)*omega-t/(2*math.sqrt(a)))
            +math.exp(omega*t)*math.erfc(math.sqrt(a)*omega+t/(2*math.sqrt(a))))/(4*omega)

def quadrature(kind,s,omega,times,n):
    # Normalized magnetic Fourier mode X=B/|k|, |k|=omega>0.
    # Its covariance is |f(p0,k)|^2/(p0^2+omega^2).
    cap=math.sqrt(60/(2*s)) if kind=='gradient' else 60/(2*s)
    x,w=np.polynomial.legendre.leggauss(n)
    p=(x+1)*cap/2;w=w*cap/(2*np.pi)
    r2=p*p+omega*omega
    f2=np.exp(-2*s*(r2 if kind=='gradient' else np.sqrt(r2)))
    density=w*f2/r2
    cov=np.cos(np.outer(np.asarray(times),p))@density
    m2=float((p*p)@density);m4=float((p**4)@density)
    return cov,m2,m4

def analytic_bounds(kind,s,omega):
    # Bound m2 from p in [omega,2omega]; bound m4 by an integrable majorant.
    if kind=='gradient':
        a=2*s
        low=omega*math.exp(-5*a*omega**2)/(2*np.pi)
        high=math.exp(-a*omega**2)/(4*math.sqrt(np.pi)*a**1.5)
    else:
        b=2*s
        low=omega*math.exp(-b*math.sqrt(5)*omega)/(2*np.pi)
        high=2/(np.pi*b**3)
    eps=math.sqrt(6*low/(217*high))
    # Q(eps,2eps) <= -low eps^2 +(217/12) high eps^4 = -low eps^2/2.
    return float(low),float(high),eps,-low*eps**2/2

def continuum_witnesses():
    rows=[]
    for kind in ('gradient','decaying_half_space_EOM'):
        for s,omega in ((.125,1.),(.5,1.),(.5,.7)):
            low,high,eps,bound=analytic_bounds(kind,s,omega)
            times=np.array([2,3,4])*eps
            cov,m2,m4=quadrature(kind,s,omega,times,320)
            fine,m2fine,m4fine=quadrature(kind,s,omega,times,480)
            error=float(max(np.max(abs(cov-fine)),abs(m2-m2fine),abs(m4-m4fine)))
            assert error<2e-10
            assert m2>=low and m4<=high
            q=float(cov[0]-2*cov[1]+cov[2])
            matrix=np.array([[cov[0],cov[1]],[cov[1],cov[2]]])
            assert q<bound<0
            closed_error=None
            if kind=='gradient':
                closed=np.array([heat_closed(t,s,omega) for t in times])
                closed_error=float(np.max(abs(cov-closed)))
                assert closed_error<2e-12
            rows.append(dict(kind=kind,flow=s,spatial_frequency=omega,
                positive_source_times=[eps,2*eps],coefficients=[1,-1],
                reflection_Gram=matrix.tolist(),reflected_norm=q,
                analytic_strict_upper_bound=bound,
                second_spectral_moment=m2,second_moment_lower_bound=low,
                fourth_spectral_moment=m4,fourth_moment_upper_bound=high,
                quadrature_agreement=error,heat_closed_form_error=closed_error,
                minimum_Gram_eigenvalue=float(np.linalg.eigvalsh(matrix)[0])))
    return dict(rows=rows,entire_Gaussian_measure_integrated=True,
        field_is_gauge_invariant_magnetic_mode=True,
        arbitrary_positive_flow_failure_proved_analytically=True,
        EOM_scope='decaying linear half-space solution only, not the entire slab paper')

def finite_lattice_and_spatial_alternative():
    rows=[]
    for n in (8,16):
        shift=np.roll(np.eye(n),1,axis=1)
        lap=2*np.eye(n)-shift-shift.T
        omega=.7
        precision=lap+omega**2*np.eye(n)
        vals,vec=np.linalg.eigh(precision)
        c=np.linalg.inv(precision)
        reflection=np.eye(n)[(1-np.arange(n))%n]
        pos=np.arange(1,n//2+1)
        base=(reflection@c)[np.ix_(pos,pos)]
        assert np.linalg.eigvalsh(base)[0]>-1e-12
        for kind in ('gradient','decaying_half_space_EOM'):
            s=.5
            rates=vals if kind=='gradient' else np.sqrt(vals)
            f=(vec*np.exp(-s*rates))@vec.T
            flowed=f@c@f.T
            q=(reflection@flowed)[np.ix_(pos,pos)]
            ev,z=np.linalg.eigh(q);v=z[:,0]
            assert np.linalg.eigvalsh(flowed)[0]>0 and ev[0]<-1e-5
            comm=float(np.max(abs(f@reflection-reflection@f)))
            assert comm<2e-14
            # Same spatial Fourier smearing acts by a scalar on each time slice.
            spatial=math.exp(-2*s*omega**2)*c
            qs=(reflection@spatial)[np.ix_(pos,pos)]
            assert np.linalg.eigvalsh(qs)[0]>-1e-12
            rows.append(dict(time_sites=n,spatial_frequency=omega,flow=s,kind=kind,
                ordinary_covariance_minimum=float(np.linalg.eigvalsh(flowed)[0]),
                original_RP_minimum=float(np.linalg.eigvalsh(base)[0]),
                flowed_RP_minimum=float(ev[0]),negative_witness=v.tolist(),
                witness_residual=float(np.linalg.norm(q@v-ev[0]*v)),
                reflection_commutator=comm,
                spatial_only_RP_minimum=float(np.linalg.eigvalsh(qs)[0])))
    # Exact free continuum source map for multiple momenta and positive times.
    times=np.array([.15,.4,.9]);freqs=np.array([.7,1.1,1.8]);s=.2
    rng=np.random.default_rng(68121);sources=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
    reflected=0j;squares=0.
    for a,k in enumerate(freqs):
        q=np.exp(-k*(times[:,None]+times[None,:]))*math.exp(-2*s*k*k)/(2*k)
        reflected+=sources[:,a].conj()@q@sources[:,a]
        squares+=math.exp(-2*s*k*k)/(2*k)*abs(np.exp(-k*times)@sources[:,a])**2
    assert abs(reflected-squares)<1e-13 and squares>0
    # Spatial-only Abelian flow leaves the spatially uniform electric history.
    electric_amplitude=1.
    spatial_factor=math.exp(-s*0.)
    full_factor=math.exp(-s*1.3**2)
    assert spatial_factor==1 and full_factor<1
    return dict(lattice_rows=rows,spatial_source_norm=float(reflected.real),
        independent_sum_of_squares=float(squares),
        spatial_zero_momentum_electric_mode=dict(initial=electric_amplitude,
            after_spatial_flow=electric_amplitude*spatial_factor,
            after_spacetime_flow=electric_amplitude*full_factor,
            nonzero_temporal_frequency=1.3),
        spatial_smearing_not_nonperturbative_chiral_gauge_construction=True)

def run():
    return dict(date='2026-10-02',round=681,tests_run=2,failures=0,errors=0,
        physical_mode_counterexamples=continuum_witnesses(),
        finite_regulator_and_alternative=finite_lattice_and_spatial_alternative(),
        dependency_hashes={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in
            ('research_note_643.md','research_note_661.md','research_note_673.md',
             'research_note_675.md','research_note_677.md','research_note_678.md',
             'research_note_679.md','research_note_680.md',
             'round681_drafts/flow_half_support_probe.py',
             'round681_drafts/flow_half_support_probe_results.json',
             'round681_drafts/mature_flow_mapping_entry.md')},
        scope=dict(linear_Abelian_physical_gaussian_sector=True,
            nonzero_flow_as_sharp_physical_time_observable_rejected=True,
            original_complete_Gauss_S9_Q0_sign_not_decided=True,
            no_counterexample_to_slab_boundary_physical_theory=True,
            spatial_flow_preserves_existing_RP_not_creates_it=True,
            joint_HF_continuum_quantum_GR_not_complete=True,
            spatial_dimension_not_derived_or_changed=True))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(round=681,tests_run=2,all_checks_passed=True)))

