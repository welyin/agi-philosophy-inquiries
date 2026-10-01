"""630: original whole-generation mass, Lorentzian cuts and common Weyl source.

A continuum free-fermion vacuum / one-fermion-loop background calculation.
No finite-graph continuum limit, full interacting SM, or Einstein derivation.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_fermion_scalar_loop_matching as old
import joint_fermion_gauss_completion as finite

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_continuum_source_spectrum_results.json'
QSTAR=.27


def quadrature(n=192):
    x,w=np.polynomial.legendre.leggauss(n)
    return (x+1)/2,w/2


def original_data():
    M0,direction=old.physical_mass_direction()
    f=np.sqrt(6)*np.sinh(QSTAR/np.sqrt(6))
    mass=M0*f
    # All-left symmetric matrix: half weight per singular value.
    # Dirac pairs occur twice; no extra spin/Nambu multiplicity.
    masses=np.linalg.svd(mass,compute_uv=False)
    weights=np.full(16,.5)
    phi=np.sqrt(6*old.original.M)*np.tanh(QSTAR/np.sqrt(6))*direction
    return mass,masses,weights,phi,direction


def density(omega,masses,weights):
    omega=np.asarray(omega,float)
    out=np.zeros_like(omega)
    for m,d in zip(masses,weights):
        active=omega>2*m
        w=omega[active]
        out[active]+=d*m*m*w*w/(4*np.pi)*(1-(2*m/w)**2)**1.5
    return out


def finite_subtracted(z,masses,weights,n=192):
    # z^4/pi int rho(w)/(w^3 (w^2-z^2)) dw,
    # x=2m/w removes the infinite endpoint. Valid below all thresholds here.
    x,w=quadrature(n)
    total=0j
    for m,d in zip(masses,weights):
        total+=d/(16*np.pi**2)*np.dot(w,x*(1-x*x)**1.5/(1-(z*x/(2*m))**2))
    return z**4*total


def mass_geometry_check():
    mass,masses,weights,phi,direction=original_data()
    h,delta=finite.mass_matrices(phi)
    BdG=np.block([[h,delta],[-delta.conj(),-h.T]])
    exact=np.linalg.eigvalsh(BdG)[32:]
    mass_error=float(np.max(abs(exact-np.sort(np.repeat(masses,2)))))
    assert mass_error<3e-15
    a=mass.conj().T@mass
    S2=float(np.dot(weights,masses**2));S4=float(np.dot(weights,masses**4))
    assert abs(S2-np.trace(a).real/2)<1e-15
    assert abs(S4-np.trace(a@a).real/2)<1e-15
    radial_error=0.;compensated_error=0.
    for q in (.11,.2,.27,.4,.65):
        p=np.sqrt(6*old.original.M)*np.tanh(q/np.sqrt(6))*direction
        hp,dp=finite.mass_matrices(p)
        lam=np.sinh(q/np.sqrt(6))/np.sinh(QSTAR/np.sqrt(6))
        radial_error=max(radial_error,float(np.max(abs(hp-lam*h))),float(np.max(abs(dp-lam*delta))))
    for sigma in (-.3,-.1,.12,.31):
        q=np.sqrt(6)*np.arcsinh(np.exp(-sigma)*np.sinh(QSTAR/np.sqrt(6)))
        p=np.sqrt(6*old.original.M)*np.tanh(q/np.sqrt(6))*direction
        hp,dp=finite.mass_matrices(p)
        compensated_error=max(compensated_error,float(np.max(abs(np.exp(sigma)*hp-h))),
                              float(np.max(abs(np.exp(sigma)*dp-delta))))
    b=float(1/np.sqrt(6)/np.tanh(QSTAR/np.sqrt(6)))
    probe=float(4*np.max(masses))
    rho=float(density(np.array([probe]),masses,weights)[0])
    v=np.array([b,1.]);kernel=rho*np.outer(v,v)
    null=np.array([1.,-b])
    joint=float(null@kernel@null)
    deleted=float(null@np.diag(np.diag(kernel))@null)
    assert max(radial_error,compensated_error)<2e-15
    assert abs(joint)<1e-14 and deleted>.01
    return dict(q_star=QSTAR,phi=phi.tolist(),F=float(old.original.F(phi)),
        Weyl_components=16,Dirac_equivalent_weight=float(weights.sum()),
        masses_sorted=np.sort(masses).tolist(),minimum_pair_threshold=float(2*min(masses)),
        S2=S2,S4=S4,old_BdG_mass_error=mass_error,radial_mass_scaling_error=radial_error,
        compensated_metric_mass_error=compensated_error,
        scalar_coordinate_source_factor=b,probe_frequency=probe,
        nonlocal_source_matrix=kernel.tolist(),compensated_quadratic_kernel=joint,
        erase_cross_kernel_value=deleted,
        not_a_gauge_symmetry_of_full_gravity_scalar_action=True)


def spectrum_matching_check():
    mass,masses,weights,phi,direction=original_data()
    S2=np.dot(weights,masses**2);S4=np.dot(weights,masses**4)
    # Independently use canonical free Dirac projectors to evaluate the cut.
    I=np.eye(2);Z=np.zeros((2,2));sigma=np.diag([1.,-1.])
    alpha=np.block([[Z,sigma],[sigma,Z]]);beta=np.diag([1.,1.,-1.,-1.])
    trace_errors=[]
    for m in masses:
        for ratio in (.21,.8,2.4):
            k=m*ratio;E=np.hypot(m,k);H=alpha*k+beta*m
            pplus=(np.eye(4)+H/E)/2;pminus=np.eye(4)-pplus
            transition=np.trace(pminus@(m*beta)@pplus@(m*beta)).real
            trace_errors.append(abs(transition-2*m*m*k*k/(E*E)))
    assert max(trace_errors)<1e-15
    # Boundary derivative of the spectral cutoff integral, not a log fit.
    cutoff=500*max(masses);rho=float(density(np.array([cutoff]),masses,weights)[0])
    leading=S2*cutoff**2/(4*np.pi)
    observed_mass_log=(rho-leading)/np.pi
    expected_mass_log=-3*S4/(2*np.pi**2)
    observed_kinetic_log=rho/(np.pi*cutoff**2)
    expected_kinetic_log=S2/(4*np.pi**2)
    assert abs(observed_mass_log/expected_mass_log-1)<3e-5
    assert abs(observed_kinetic_log/expected_kinetic_log-1)<3e-5
    # 599: b4_mass = 2 S4 lambda^4 + 2 S2 (d lambda)^2.
    # Gamma_div=b4/(32 pi^2 epsilon); 1/epsilon -> 2 log cutoff.
    # Gamma_E''=-Pi_E, so mass and Lorentz-frequency kinetic signs differ.
    heat_mass_log=-(24*S4)/(16*np.pi**2)
    heat_kinetic_log=4*S2/(16*np.pi**2)
    assert abs(heat_mass_log-expected_mass_log)<1e-15
    assert abs(heat_kinetic_log-expected_kinetic_log)<1e-15
    # Exact zero-frequency dispersive integral with k-cutoff is the
    # second derivative of the same Dirac-sea energy, with minus sign.
    x,w=quadrature(256);K=3*max(masses);ks=K*x
    sea_curvature=0.;dispersion=0.
    for m,d in zip(masses,weights):
        sea_curvature-=d*m*m/np.pi**2*K*np.dot(w,ks**4/(ks*ks+m*m)**1.5)
        upper=2*np.hypot(K,m)
        om=2*m+(upper-2*m)*x
        one=density(om,np.array([m]),np.array([d]))
        dispersion+=(upper-2*m)*np.dot(w,one/(np.pi*om))
    assert abs(sea_curvature+dispersion)<2e-11
    return dict(max_Dirac_projector_cut_error=float(max(trace_errors)),
        ultraviolet_cutoff=float(cutoff),mass_log_spectral=float(observed_mass_log),
        mass_log_from_original_heat_kernel=float(heat_mass_log),
        kinetic_log_spectral=float(observed_kinetic_log),
        kinetic_log_from_original_heat_kernel=float(heat_kinetic_log),
        sea_energy_second_derivative=float(sea_curvature),
        zero_frequency_spectral_integral=float(dispersion),
        static_identity_error=float(abs(sea_curvature+dispersion)),
        continuum_input='Canonical Lorentzian fermions, constant original scalar background, flat Einstein metric and vacuum; gauge/scalar/metric loops omitted.',
        quadratic_power_divergence_scheme_dependent=True,
        finite_contact_coefficients_not_fixed_by_spectral_cut=True)


def shared_window_check():
    mass,masses,weights,phi,direction=original_data()
    gap=2*min(masses)
    L4=float(weights.sum()/(80*np.pi**2))
    L6=float(np.sum(weights/masses**2)/(1120*np.pi**2))
    rows=[]
    for ratio in (.1,.35,.65,.9):
        z=ratio*gap
        value=finite_subtracted(z,masses,weights)
        fine=finite_subtracted(z,masses,weights,384)
        upper=z**4*L4/(1-ratio**2)
        remainder_bound=z**6*L6/(1-ratio**2)
        assert abs(value-fine)<2e-13
        assert 0<=value.real<=upper*(1+1e-12) and abs(value.imag)<1e-16
        assert 0<=value.real-z**4*L4<=remainder_bound*(1+1e-12)
        rows.append(dict(frequency=float(z),fraction_of_first_threshold=ratio,
                         twice_subtracted_kernel=float(value.real),
                         leading_fourth_order=float(z**4*L4),
                         bound=float(upper),sixth_order_remainder_bound=float(remainder_bound),
                         quadrature_change=float(abs(value-fine))))
    # Vacuum Gaussian source pulses, normalized per spatial volume.
    # Source is an external diagnostic; no finite apparatus has been derived.
    pulses=[]
    for tau in (1.,3.,8.):
        x,w=quadrature(320);work=0.;prob=0.
        for m,d in zip(masses,weights):
            # w=2m+y/tau; y in [0,12], exponential tail < exp(-144).
            om=2*m+12*x/tau
            rr=density(om,np.array([m]),np.array([d]))
            base=12/tau*w*rr*np.exp(-(tau*om)**2)*tau*tau
            prob+=float(np.sum(base));work+=float(np.dot(base,om))
        assert work>0 and prob>0
        pulses.append(dict(tau=tau,excitation_density_per_epsilon_squared=prob,
                           absorbed_energy_density_per_epsilon_squared=work))
    # A local real polynomial has no positive-frequency cut or vacuum noise.
    frequencies=gap*np.array([.5,1.01,1.5,4.,20.])
    spectral=density(frequencies,masses,weights)
    assert spectral[0]==0 and np.all(spectral[1:]>0)
    return dict(L4=L4,L6=L6,below_threshold_rows=rows,
        probe_frequencies=frequencies.tolist(),commutator_spectral_density=spectral.tolist(),
        symmetric_vacuum_noise=(spectral/2).tolist(),Gaussian_pulse_rows=pulses,
        local_real_action_alone_cannot_match_open_pair_channels=True,
        infinite_space_vacuum_not_original_finite_graph_Gibbs_state=True,
        actual_record_instrument_and_graph_continuum_map_not_constructed=True)


def run():
    deps=('research_note_599.md','research_note_600.md','research_note_602.md',
          'research_note_620.md','research_note_625.md','research_note_626.md',
          'research_note_627.md','research_note_628.md','research_note_629.md',
          'joint_fermion_scalar_loop_matching.py','joint_fermion_gauss_completion.py',
          'joint_curved_quantum_source.py')
    return dict(round=630,tests_run=3,failures=0,errors=0,
        original_mass_and_metric=mass_geometry_check(),
        spectral_and_UV_matching=spectrum_matching_check(),
        common_effective_window=shared_window_check(),
        dependency_hashes={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in deps},
        scope='One generation of the original nonzero mass data, same radial ray as 599, canonical 3+1 continuum free fermions and vacuum. Exact connected quadratic source cut, twice-subtracted dispersion and nonlocal Weyl/scalar-source relation in the one-fermion-loop background sector. Local counterterms and trace anomaly retained as separate data. No full interacting theory, finite-graph matching, physical parameter fit, actual apparatus or GR derivation.')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:
        assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps(result,ensure_ascii=False,indent=2))

