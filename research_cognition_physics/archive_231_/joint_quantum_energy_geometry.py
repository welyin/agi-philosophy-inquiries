"""578: original curved-target kinetic gap, finite volume and gravitational matching.
Spectral geometry is inherited mathematics; the application constrains one unrenormalized
finite-graph-to-geometry prescription. No continuum quantum gravity or vacuum energy fit.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_curved_quantum_source as original

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_quantum_energy_geometry_results.json'
R=np.sqrt(6.)
GAP=2/3
SCALAR_CURVATURE=-10/3


def quadrature(left,right,n=240):
    q,w=np.polynomial.legendre.leggauss(n)
    return (left+right)/2+(right-left)*q/2,(right-left)*w/2


def radial_form_check():
    rows=[]
    for ell,left,right,phase in ((0,.2,2.6,0.),(1,.7,4.2,.37),(2,2.,8.,-.23)):
        r,w=quadrature(left,right)
        z=(2*r-left-right)/(right-left);dz=2/(right-left)
        amp=np.exp(-1/(1-z*z));derivative=amp*(-2*z)*dz/(1-z*z)**2
        u=amp*np.exp(1j*phase*r)
        du=(derivative+1j*phase*amp)*np.exp(1j*phase*r)
        S=R*np.sinh(r/R)
        norm=float(np.dot(w,abs(u)**2))
        # Direct original weighted radial derivative after f=u/S^2.
        direct=float(np.dot(w,abs(du-2/(R*np.tanh(r/R))*u)**2
                               +ell*(ell+3)*abs(u)**2/S**2)/norm)
        transformed=float(np.dot(w,abs(du)**2
            +(GAP+(2+ell*(ell+3))/S**2)*abs(u)**2)/norm)
        error=abs(direct-transformed)
        assert error<2e-12 and direct>GAP
        rows.append(dict(ell=ell,interval=[left,right],phase=phase,
                         direct_rayleigh=direct,transformed_rayleigh=transformed,identity_error=error))
    # Verify that the inherited five-field target indeed has this polar metric.
    r=np.linspace(.04,7.,17);coordinate_radius=np.sqrt(6*original.M)*np.tanh(r/R)
    metric_errors=[]
    for rho,x in zip(r,coordinate_radius):
        phi=np.array([x,0.,0.,0.,0.]);K=original.metric(phi)
        drho_dx=np.sqrt(original.M)/original.F(phi)
        metric_errors.append(max(abs(K[0,0]-drho_dx**2),
                                 abs(K[1,1]*x*x-6*np.sinh(rho/R)**2)))
    assert max(metric_errors)<2e-10
    return dict(rows=rows,original_target_polar_metric_error=max(metric_errors),
                scalar_curvature=SCALAR_CURVATURE,kinetic_gap=GAP)


def sharp_kinetic_gap_check():
    rows=[]
    for length in (4.,8.,16.,32.):
        r,w=quadrature(length,2*length)
        u=np.sin(np.pi*(r-length)/length)
        du=np.pi/length*np.cos(np.pi*(r-length)/length)
        norm=np.dot(w,u*u)
        value=float(np.dot(w,du*du+(GAP+1/(3*np.sinh(r/R)**2))*u*u)/norm)
        upper=GAP+np.pi**2/length**2+1/(3*np.sinh(length/R)**2)
        assert GAP<value<=upper+1e-12
        rows.append(dict(length=length,rayleigh=value,excess_over_gap=value-GAP,
                         analytic_upper=float(upper)))
    assert rows[-1]['excess_over_gap']<rows[0]['excess_over_gap']/50
    return dict(rows=rows,
                trial_is_compact_H1_radial_function_approximable_by_smooth_ones=True,
                sharpness_only_for_free_target_laplacian_not_full_matter_H=True)


def graph_volume_bound_check():
    hbar=.07;volume=2.3;rows=[]
    for side in (2,4,8,16):
        count=side**3
        angle=2*np.pi*np.arange(count)/count
        weights=1+.21*np.sin(angle)+.09*np.cos(3*angle)
        weights=weights/np.sum(weights)*volume
        original_bound=float(hbar*hbar/3*np.sum(1/weights))
        volume_bound=float(hbar*hbar*count**2/(3*volume))
        assert original_bound>=volume_bound*(1-1e-14)
        epsilon=1/side
        rows.append(dict(side=side,nodes=count,volume=float(np.sum(weights)),
                         exact_weight_bound=original_bound,volume_only_bound=volume_bound,
                         scaled_volume_bound=volume_bound*epsilon**6))
    scaled=[q['scaled_volume_bound'] for q in rows]
    assert max(scaled)-min(scaled)<1e-14
    # Fixed geometry cap |R|+tau^2 <= L, Einstein convention R-|K|²+tau²=2 kappa rho.
    L=3.;kappa=1.;max_nodes=np.sqrt(3*L/(2*kappa))*volume/hbar
    compatible=[bool(q['nodes']<=max_nodes) for q in rows]
    assert compatible[0] and not compatible[-1]
    return dict(hbar=hbar,rows=rows,curvature_trace_cap=L,gravitational_coefficient=kappa,
                necessary_node_bound=float(max_nodes),passes_necessary_bound=compatible,
                source_energy_identification_is_an_additional_matching_condition=True,
                bound_applies_to_all_states_not_just_product_sources=True)


def metric_subtraction_check():
    count=27;hbar=.07;volume0=2.3
    A=hbar*hbar*count**2/3
    def gap_energy(volume):return A/volume
    rows=[]
    for volume in (1.2,volume0,4.7):
        step=volume*1e-4
        # Independent four-point variation of the specified subtraction functional.
        derivative=(gap_energy(volume-2*step)-8*gap_energy(volume-step)
                   +8*gap_energy(volume+step)-gap_energy(volume+2*step))/(12*step)
        pressure=-derivative;density=gap_energy(volume)/volume
        assert abs(pressure-density)<2e-10
        rows.append(dict(volume=volume,gap_energy=gap_energy(volume),
                         subtraction_functional_pressure=float(pressure),density=density,
                         pressure_over_density=float(pressure/density)))
    lam=A/(volume0*volume0)
    matched_energy=gap_energy(volume0)-lam*volume0
    residual_pressure=A/volume0**2+lam
    assert abs(matched_energy)<1e-14 and residual_pressure>0
    xi=.2
    shifted_gap=GAP+xi*SCALAR_CURVATURE
    assert abs(shifted_gap)<1e-14
    return dict(rows=rows,constant_density_matched_at_volume=volume0,
                matched_energy_residual=matched_energy,
                residual_pressure_of_gap_minus_constant_density=residual_pressure,
                target_curvature_ordering_xi=xi,shifted_free_kinetic_gap=shifted_gap,
                varying_fixed_node_count=True,not_total_state_pressure_or_actual_vacuum_energy=True,
                cancelling_this_lower_bound_does_not_prove_continuum_limit=True)


def run():
    checks=('radial_form_check','sharp_kinetic_gap_check','graph_volume_bound_check','metric_subtraction_check')
    evidence={name:globals()[name]() for name in checks}
    deps=('joint_curved_quantum_source.py','research_note_574.md','research_note_572.md',
          'research_note_577.md','research_round_577_checks.json')
    return dict(round=578,tests_run=len(checks),failures=0,errors=0,checks=list(checks),
        evidence=evidence,dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope='original positive Laplace-Beltrami fixed-graph Hamiltonian: all-state kinetic lower bound; conditional obstruction to unrenormalized fixed-hbar bounded-volume regular Einstein matching; fixed-geometry energy shifts are not automatically the same metric source; no full ground-state energy, continuum limit, spacetime dimension or general unification no-go')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        payload=json.dumps(result,ensure_ascii=False,indent=2)+'\n'
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(payload)
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps({k:result[k] for k in ('round','tests_run','failures','errors')}))
