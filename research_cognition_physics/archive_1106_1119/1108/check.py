"""1108: a shared physical-translation contract ties clock and inertial energy.

The general representation statements are analytical, not finite-matrix tests.
No material deployment, exact rest packet, or complete SR model is certified.
"""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OLD = ROOT / 'research_cognition_physics/archive_923_934/925/record_resource_exchange.py'


def norm(a):
    return float(np.linalg.norm(a, 2))


def fmat(h, f):
    e, v = np.linalg.eigh(h)
    return (v * f(e)) @ v.conj().T


def evolve(h, t):
    return fmat(h, lambda e: np.exp(-1j * t * e))


def run():
    spec = importlib.util.spec_from_file_location('source925_for1108', OLD)
    old = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(old)
    z = old.model()
    eye = np.eye(z['d'], dtype=complex)
    indices = [z['idx'](loc, q, 3 + (q if loc else 0))
               for q in range(3) for loc in range(2)]
    w = eye[:, indices]
    h = w.conj().T @ z['h'] @ w
    h0 = w.conj().T @ z['h0'] @ w
    interaction = h - h0
    identity = np.eye(6)
    source = identity[:, [0, 2, 4]]
    receiver = identity[:, [1, 3, 5]]
    projector = receiver @ receiver.conj().T
    m = 100.
    mass0 = m * identity + h0
    mass = m * identity + h
    T = np.pi / (2 * z['g'])
    phase_target = evolve(h, T) @ source
    source_checks = dict(
        closed_subspace_error=norm(z['h'] @ w - w @ h),
        original_dimension=z['d'], retained_dimension=6,
        bare_interaction_commutator=norm(h0 @ interaction - interaction @ h0),
        internal_eigenvalues=np.linalg.eigvalsh(h).tolist(),
        rest_transfer_isometry_error=norm(projector @ phase_target - phase_target))
    assert source_checks['closed_subspace_error'] < 1e-13
    assert source_checks['bare_interaction_commutator'] < 1e-13
    assert source_checks['rest_transfer_isometry_error'] < 1e-12

    rows = []
    for speed in (0., .3, .6, .8):
        gamma = 1 / np.sqrt(1 - speed**2)
        t, x = gamma * T, speed * gamma * T
        good_e, good_p = gamma * mass, gamma * speed * mass
        bad_p = gamma * speed * mass0
        # Independent evaluation of the proposed bad dispersion on its
        # common-velocity spectral preparation; all matrices commute here.
        bad_e = fmat(bad_p @ bad_p + mass0 @ mass0, np.sqrt) + interaction
        good_phase = t * good_e - x * good_p
        bad_phase = t * bad_e - x * bad_p
        ug, ub = evolve(good_phase, 1.), evolve(bad_phase, 1.)
        prob = lambda u: (source.conj().T @ u.conj().T @ projector @ u @ source)
        pg, pb = prob(ug), prob(ub)
        expected_bad = float(np.sin(gamma * np.pi / 2)**2)
        boost_e_error = norm(bad_e - gamma * mass)
        boost_p_error = norm(bad_p - gamma * speed * mass)
        expected_shell_defect = 2 * (gamma - 1) * mass0 @ interaction
        shell_defect = bad_e @ bad_e - bad_p @ bad_p - mass @ mass
        row = dict(speed=speed, gamma=gamma, coordinate_time=t,
                   good_all_input_probability_error=norm(pg - np.eye(3)),
                   bad_probability=float(np.trace(pb).real / 3),
                   bad_analytic_probability=expected_bad,
                   bad_all_input_probability_identity_error=norm(pb - expected_bad * np.eye(3)),
                   rest_energy_generator_difference=norm((mass0 + interaction) - mass),
                   bad_energy_boost_defect=boost_e_error,
                   bad_momentum_boost_defect=boost_p_error,
                   bad_mass_shell_defect=norm(shell_defect),
                   shell_defect_formula_error=norm(shell_defect - expected_shell_defect),
                   good_path_phase_identity_error=norm(good_phase - T * mass),
                   bad_path_phase_identity_error=norm(bad_phase - T * (mass0 + gamma * interaction)),
                   fixed_lab_generator_conservation_error=norm(ub.conj().T @ bad_e @ ub - bad_e))
        assert row['good_all_input_probability_error'] < 3e-12
        assert row['bad_all_input_probability_identity_error'] < 3e-12
        assert row['shell_defect_formula_error'] < 2e-10
        assert row['good_path_phase_identity_error'] < 3e-12
        assert row['bad_path_phase_identity_error'] < 3e-12
        assert row['fixed_lab_generator_conservation_error'] < 3e-11
        assert abs(boost_e_error - (gamma - 1) * z['g']) < 2e-12
        assert abs(boost_p_error - gamma * speed * z['g']) < 2e-12
        rows.append(row)
    assert abs(rows[2]['bad_probability'] - (2 + np.sqrt(2)) / 4) < 2e-12

    # Common clock h does not select the absolute rest-energy offset.
    # At rest both offset choices have precisely the same internal channel.
    p = 20.
    offsets = []
    for offset in (100., 101.):
        M = offset * identity + h
        E = fmat(p*p*identity + M @ M, np.sqrt)
        velocity = p * np.linalg.inv(E)
        offsets.append(dict(offset=offset,
            rest_channel_phase_removed_error=norm(np.exp(1j*offset*T)*evolve(M,T)-evolve(h,T)),
            common_momentum_velocity_eigenvalues=np.linalg.eigvalsh(velocity).tolist()))
        assert offsets[-1]['rest_channel_phase_removed_error'] < 3e-12
    assert abs(offsets[0]['common_momentum_velocity_eigenvalues'][0] -
               offsets[1]['common_momentum_velocity_eigenvalues'][0]) > 1e-4

    # A finite-width velocity packet: discrete quadrature checks the proven
    # support-uniform bound. It is NOT a finite-dimensional Poincare rep.
    eigenvalues = np.linalg.eigvalsh(h)
    center = (eigenvalues[0] + eigenvalues[-1]) / 2
    h_centered = h - center * identity
    packet_rows = []
    for radius in (.02, .05, .1):
        rapidities = np.linspace(-radius, radius, 9)
        weights = np.array([1.,2.,3.,4.,5.,4.,3.,2.,1.])
        weights /= weights.sum()
        exact, ideal = [], []
        for r, weight in zip(rapidities, weights):
            k = np.cosh(r)
            exact.append(np.sqrt(weight) * evolve(mass, T*k))
            ideal.append(np.sqrt(weight) * np.exp(-1j*(m+center)*T*k) * evolve(h_centered,T))
        exact, ideal = np.vstack(exact), np.vstack(ideal)
        bound = float(T * norm(h_centered) * (np.cosh(radius)-1))
        actual = norm(exact - ideal)
        assert norm(exact.conj().T @ exact - identity) < 1e-12
        assert actual <= bound + 1e-12
        packet_rows.append(dict(relative_rapidity_support=radius,
            quadrature_isometry_error=actual, analytic_complete_reference_bound=bound,
            finite_quadrature_certifies_continuum_representation=False,
            physical_localized_readout_certified=False))

    # The complete algebra clause has content. Population-only comparisons
    # cannot distinguish these two rest generators; one coherent test can.
    hz = np.diag([0., 1.]); rz = np.diag([100., 102.])
    basis = np.eye(2); plus = np.ones(2)/np.sqrt(2)
    pop_error = max(norm(evolve(hz,T) @ np.outer(b,b) @ evolve(hz,T).conj().T -
                         evolve(rz,T) @ np.outer(b,b) @ evolve(rz,T).conj().T) for b in basis)
    p_clock = float(abs(np.vdot(plus, evolve(hz,np.pi) @ plus))**2)
    p_mass = float(abs(np.vdot(plus, evolve(rz,np.pi) @ plus))**2)
    assert pop_error < 1e-12 and p_clock < 1e-12 and abs(p_mass-1) < 1e-12

    sources = [OLD,
        ROOT/'research_cognition_physics/archive_935_955/research_note_941.md',
        ROOT/'research_cognition_physics/archive_1086_1105/_admission/dynamic_contract_audit.md',
        ROOT/'research_cognition_physics/认知时空操作公理与洛伦兹运动学_阶段论文.md']
    return dict(round=1108,status='PASS_CONDITIONAL_TRANSLATION_BRIDGE',
        groups=dict(original_process=source_checks,
            omitted_interaction_inertia=rows,absolute_mass_freedom=offsets,
            finite_velocity_width=packet_rows,
            restricted_readout=dict(population_error=pop_error,
                coherent_clock_probability=p_clock,coherent_mass_probability=p_mass)),
        source_sha256={str(q.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(q.read_bytes()).hexdigest() for q in sources},
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        new_adopted_axioms=0,scientific_calibration_increment=0,
        empirical_data=False,actual_clock_port_implementation=False,
        full_SR_joint_model=False,goal_complete=False)


def compare(old, new):
    if isinstance(new,dict):
        assert old.keys()==new.keys()
        for k in new:compare(old[k],new[k])
    elif isinstance(new,list):
        assert len(old)==len(new)
        for x,y in zip(old,new):compare(x,y)
    elif isinstance(new,float):
        assert np.isclose(old,new,atol=1e-11,rtol=1e-10),(old,new)
    else:assert old==new,(old,new)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args(); result=run()
    if args.write_results:
        with (HERE/'results.json').open('x',encoding='utf8') as f:
            json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
    else:compare(json.loads((HERE/'results.json').read_text(encoding='utf8')),result)
    print(json.dumps(dict(round=1108,status=result['status'],groups=len(result['groups']),
                         saved_results_match=not args.write_results)))
