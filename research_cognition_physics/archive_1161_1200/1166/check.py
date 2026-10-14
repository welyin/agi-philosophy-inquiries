"""1166: actual old clock/port identification and flat-response non-entailment.

Default replays two frozen actual entry points and compares this saved result.
--write-once creates only this round's first results.json. No null probe,
clock reset, new geometry instrument or new physical mother permission.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import subprocess
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
ROOT = BASE.parent
TOL = 3e-9
ETA = np.diag([1., -1., -1., -1.])
CAPTURE = math.pi / 4
COORDINATE_DURATION = 1.2
SOURCES = {
    '02_可组合认知结构与有限维量子理论_阶段论文.md': '381dfab12a958a68bd5b36bf0353053c9cc5e0f24fcfd4f69f0cddfe6396e953',
    '03_认知操作五公理与三维关系空间_阶段论文.md': 'd14a56d9d9dd138cf39dbb652cd18b94427efa4ed22e9b8760d5af144f8a7d1d',
    '04_认知时空操作公理与洛伦兹运动学_阶段论文.md': 'ac72400f9043b72ac4344ea350ee7a10226a205db81c319d02e76a306a6c3b86',
    'archive_1106_1119/_shared/1116_1119/adoption_and_dependencies.md': '1c1395c16cf2392be0c6e7275d660f454113290bd0e23d6892e08b3aa41b4a0e',
    'archive_1120_1158/_shared/1157_1158/adoption_and_dependencies.md': 'f7af93c680768e4e786edf9aad5e1f95c8c6759e69f2ac94137219960bfe95af',
    'archive_1106_1119/1112/check.py': '6e49d71ec5535fd6999432481938663a5543e2704db993ec7a00b8f475b6cd9f',
    'archive_1106_1119/1114/proof.md': '36529e5f56971a8baf66065928b5b7d1899d6444f257b1cc7008a01fd42be9a2',
    'archive_1106_1119/1115/proof.md': 'b7b2239946e44cb1aade62d1282525d998297deb6dcadc36b3005dc7c454a8da',
    'archive_1106_1119/1115/results.json': 'e17b36f5835e44564be721722b9242490e8727314f89d2c0e0b663cba5eebc68',
    'archive_1106_1119/1119/proof.md': 'b4a97154a8302a49c4f4f6788f615909b4849820bfa42e3bf11574c16808008d',
    'archive_1106_1119/1119/check.py': 'efe8f4aabafb94d146e13f773337ddc453afb26b7705d3eacd92f56244f7a863',
    'archive_1106_1119/1119/results.json': 'fc48bf4935d25f975eae2b2337fe87d474f976e325b4da3847e3ec47c13c5d27',
    'archive_1120_1158/1158/proof.md': '8d05617c800014afd8965eaf51874fa3235e16cf2e49a2918ad00ebe5a8e7608',
    'archive_1120_1158/1158/check.py': 'd0f56b05cad8f576ab608896a0165d807ea9b4c128e8cbe642ee65af2cd0d3fc',
    'archive_1120_1158/1158/results.json': '24187911743d4cb8660c0c717fd74bbdd835874a482e7a1870a9011d17056929',
    'archive_1120_1158/1158/acceptance.json': 'c8d9b32be4c04552b02beed4c02437f9e86732f258e7acadfac0418b6fecaf0c',
    'archive_1106_1119/1111/results.json': '7e5b27990183dd54a75ce9eda35340ce79ae3d1c24b31b25dfb2064bc9d257b3',
}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def small(value, label):
    value = float(np.linalg.norm(value))
    assert value < TOL, (label, value)
    return value


def read_json(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def compare(actual, saved):
    if isinstance(saved, dict):
        assert actual.keys() == saved.keys()
        for key in saved:
            compare(actual[key], saved[key])
    elif isinstance(saved, list):
        assert len(actual) == len(saved)
        for a, b in zip(actual, saved):
            compare(a, b)
    elif isinstance(saved, float):
        assert math.isclose(actual, saved, rel_tol=2e-9, abs_tol=3e-10), (actual, saved)
    else:
        assert actual == saved, (actual, saved)


def replay():
    rows = []
    runner = ROOT / 'scripts/run_research_gr_active.py'
    for name, marker in [
        ('archive_1106_1119/1119/check.py', 'saved results matched'),
        ('archive_1120_1158/1158/check.py', 'READ_ONLY_MATCH'),
    ]:
        run = subprocess.run([sys.executable, '-B', '-X', 'utf8', str(runner),
                              '--script', str(BASE/name)], cwd=ROOT,
                             capture_output=True, text=True, encoding='utf-8', timeout=120)
        assert run.returncode == 0, (name, run.stdout, run.stderr)
        assert marker in run.stdout, (name, run.stdout)
        rows.append(dict(script=name, exit_code=run.returncode, result_comparison=marker))
    return rows


PAIRS = [(i, j) for i in range(4) for j in range(i, 4)]


def design(vectors):
    return np.array([[w[i]*w[j]*(1 if i == j else 2) for i, j in PAIRS] for w in vectors])


def matrix(coefficients):
    out = np.zeros((4, 4))
    for (i, j), value in zip(PAIRS, coefficients):
        out[i, j] = out[j, i] = value
    return out


def group_clock(material_module):
    material = material_module.Material(1)
    hd = material.baseline*np.eye(288)+material.hfk
    hin = material.baseline*np.eye(288)+np.kron(material.hc, np.eye(72))
    eigenvalues, eigenvectors = np.linalg.eigh(hd)
    # The actual 1106 incoming ready code, including arbitrary qutrit columns.
    joint = np.eye(3, dtype=complex)[:, None, :]
    incoming = material.encode_joint(joint, 0.).reshape(288, 3)
    capture = material.gamma @ material_module.exp_h(hin, CAPTURE)
    captured = capture @ incoming
    captured_eigenbasis = eigenvectors.conj().T @ captured
    cap_error = small(hd@capture-capture@hin, 'same complete positive capture')
    a = .3
    velocities = [np.zeros(3)]
    velocities += [sign*a*np.eye(3)[i] for i in range(3) for sign in (-1., 1.)]
    velocities += [a*(np.eye(3)[i]+np.eye(3)[j]) for i in range(3) for j in range(i+1, 3)]
    rows, probabilities, free_times = [], [], []
    max_clock_error = max_gram = 0.
    for velocity in velocities:
        duration = COORDINATE_DURATION*math.sqrt(1-float(velocity@velocity))
        state = eigenvectors @ (np.exp(-1j*duration*eigenvalues)[:, None]*captured_eigenbasis)
        reshaped = state.reshape(4, 8, 3, 3, 3)
        p = np.sum(abs(reshaped[:, [1, 3, 5, 7], ...])**2, axis=(0, 1, 2, 3))
        expected = math.sin((CAPTURE+duration)/2)**6
        max_clock_error = max(max_clock_error, small(p-expected, 'original Mr including capture age'))
        max_gram = max(max_gram, small(state.conj().T@state-np.eye(3), 'full output isometry'))
        inferred = 2*math.asin(float(p[0])**(1/6))-CAPTURE
        small(inferred-duration, 'free proper interval, not total clock age')
        probabilities.append(float(p[0]))
        free_times.append(duration)
        rows.append(dict(velocity=velocity.tolist(), coordinate_duration=COORDINATE_DURATION,
                         free_proper_time=duration, total_clock_age=CAPTURE+duration,
                         actual_Mr_all_qutrit_probabilities=p.tolist(), inferred_free_time=inferred))
    vectors = np.column_stack([np.ones(10), np.array(velocities)])
    A = design(vectors)
    assert np.linalg.matrix_rank(A) == 10
    inferred_times = 2*np.arcsin(np.array(probabilities)**(1/6))-CAPTURE
    y = (inferred_times/COORDINATE_DURATION)**2
    recovered = matrix(np.linalg.solve(A, y))
    reconstruction_error = small(recovered-ETA, 'metric from actual full-clock probabilities')
    return dict(full_clock_dimension=288, open_qutrit_columns=3, capture_duration=CAPTURE,
                capture_energy_intertwiner_error=cap_error, maximum_probability_error=max_clock_error,
                complete_output_gram_error=max_gram, calibration_cases=rows, design_rank=10,
                recovered_bilinear_form=recovered.tolist(), reconstruction_error=reconstruction_error,
                no_reset_or_intermediate_clock_measurement=True), (A, vectors, y, np.array(probabilities), np.array(free_times))


def group_nontrivial(old):
    network = old['same_entangled_process_and_active_future']
    assert abs(network['CHSH']-2*math.sqrt(2)) < TOL
    assert network['noncommuting_active_future_residual'] < TOL
    gap = network['coherent_Bell_recovery_probability']-network['EA_dephased_different_process_Bell_recovery_probability']
    small(gap-.5, 'actual retained versus dephased record future')
    # This is exactly 1158's already-certified complete recoil menu.
    pairs = [(np.array([5., 4, 0, 0]), np.array([5., -4, 0, 0])),
             (np.array([5., 0, 4, 0]), np.array([5., 0, -4, 0]))]
    sigma = 1.2
    end_rows = []
    for label, momenta in zip(('xx', 'yy'), pairs):
        total = sum(momenta)
        small(total-np.array([10., 0, 0, 0]), 'same full pair charge')
        for body, p in zip(('A', 'B'), momenta):
            delta_z = sigma*p/3
            error = small(delta_z@ETA@delta_z-sigma*sigma, 'same interval despite changed endpoint')
            end_rows.append(dict(branch=label, body=body, four_momentum=p.tolist(),
                                 actual_free_displacement=delta_z.tolist(), proper_square_error=error))
    endpoint_difference = float(np.linalg.norm(sigma*(pairs[1][0]-pairs[0][0])/3))
    assert endpoint_difference > 2
    delta_a = pairs[1][0]-pairs[0][0]
    delta_b = pairs[1][1]-pairs[0][1]
    small(delta_a+delta_b, 'full exchange keeps recoil partner')
    saved_geo = old['actual_causal_lines_four_P_and_original_recoil']
    small(delta_a-np.array(saved_geo['individual_recoil_A']), 'same actual recoil identity')
    assert saved_geo['boosted_individual_exchange_A'][0] == 3.
    T0 = 23/15+math.pi
    speeds = [math.tanh(r+math.atanh(1/T0)) for r in (0., 1., 3., 6.)]
    assert all(0 < v < 1 for v in speeds)
    assert speeds == sorted(speeds) and speeds[-1] > .9999
    return dict(complete_record_output_dimension=network['complete_output_dimension'],
                open_input_dimension=network['open_input_dimension'], CHSH=network['CHSH'],
                record_dephasing_changes_actual_future_probability_by=gap,
                active_future_identity_error=network['noncommuting_active_future_residual'],
                retained_arbitrary_reference_identity_error=network['reference_residual'],
                original_recoil_free_legs=end_rows, branch_endpoint_difference=endpoint_difference,
                boosted_actual_individual_energy_exchange=[3., -3.],
                finite_boost_task_speeds=speeds, common_metric=ETA.tolist(),
                exact_null_probe_claimed=False)


def group_finite_and_dictionary(data):
    A, vectors, y, probabilities, free_times = data
    lo, hi, epsilon = .2, .5, 1e-6
    assert probabilities.min() > lo and probabilities.max() < hi
    L = 1/(3*lo**(5/6)*math.sqrt(1-hi**(1/3)))
    signs = np.array([1., -1., 1., -1., 1., -1., 1., -1., 1., -1.])
    perturbed = probabilities+epsilon*signs
    inferred = 2*np.arcsin(perturbed**(1/6))-CAPTURE
    reconstructed = matrix(np.linalg.solve(A, (inferred/COORDINATE_DURATION)**2))
    actual_error = float(np.max(abs(reconstructed-ETA)))
    inverse_norm = float(np.linalg.norm(np.linalg.inv(A), np.inf))
    bound = inverse_norm*(2*max(free_times)*L*epsilon+L*L*epsilon**2)/COORDINATE_DURATION**2
    assert actual_error <= bound
    q = .02
    wrong = np.diag([1., -math.exp(2*q), -math.exp(-2*q), -1.])
    wrong_defect = float(np.max(abs(np.einsum('ni,ij,nj->n', vectors, wrong, vectors)-y)))
    assert wrong_defect > .003
    S = np.diag([1., math.exp(q), math.exp(-q), 1.])
    inverse = np.linalg.inv(S)
    transported = inverse.T@ETA@inverse
    moved = vectors@S.T
    common_error = small(np.einsum('ni,ij,nj->n', moved, transported, moved)-y,
                         'common coordinate and readout dictionary')
    recovered_moved = matrix(np.linalg.solve(design(moved), y))
    moved_error = small(recovered_moved-transported, 'same flat form in transported coordinates')
    return dict(probability_calibration_window=[lo, hi], probability_error=epsilon,
                inverse_clock_lipschitz_bound=L, inverse_design_infinity_norm=inverse_norm,
                observed_coefficient_error=actual_error, certified_coefficient_error_bound=bound,
                wrong_changed_metric_fixed_endpoints_defect=wrong_defect,
                common_coordinate_transport_residual=common_error,
                recovered_transported_form_error=moved_error,
                endpoint_time_and_capture_calibration_errors_not_included=True,
                these_math_counterchecks_do_not_grant_new_geometry_instruments=True)


def calculate():
    for name, expected in SOURCES.items():
        assert digest(BASE/name) == expected, name
    old_replay = replay()
    material = load_module('r1112_for_1166', BASE/'archive_1106_1119/1112/check.py')
    clock, data = group_clock(material)
    old = read_json(BASE/'archive_1120_1158/1158/results.json')
    result = dict(schema='round1166_flat_effective_nonentailment_v1', status='PASS', check_groups=4,
        code_sha256=digest(Path(__file__)), source_sha256=SOURCES,
        group1_actual_frozen_replay=old_replay,
        group2_same_complete_clock_and_operational_metric=clock,
        group3_nontrivial_complete_process_flat_law=group_nontrivial(old),
        group4_finite_accuracy_and_dictionary=group_finite_and_dictionary(data),
        scope=dict(theory='S_eff: exact FUCP parent plus declared finite material domain and common realization',
                   conclusion='nonzero certified common geometric response is not forced in this semantics',
                   single_layer_all_physical_mother_C_not_certified=True,
                   no_distributional_stress_or_Einstein_equation_counterexample=True,
                   calibrated_quadratic_interval_class_only=True,
                   all_finite_context_quantifiers_proved_by_generator_induction=True,
                   new_adopted_axioms=0, new_empirical_calibration=0, GR_hard_deliverables_completed=0))
    for name, expected in SOURCES.items():
        assert digest(BASE/name) == expected, ('source changed during run', name)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-once', action='store_true')
    args = parser.parse_args()
    result = calculate()
    path = HERE/'results.json'
    if args.write_once:
        with path.open('x', encoding='utf-8') as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write('\n')
    else:
        compare(result, read_json(path))
    print(json.dumps(dict(status='PASS', mode='WRITE_ONCE' if args.write_once else 'READ_ONLY_MATCH',
                          groups=4, source_files=len(SOURCES), old_entry_points_replayed=2,
                          metric_rank=10,
                          metric_reconstruction_error=result['group2_same_complete_clock_and_operational_metric']['reconstruction_error'],
                          GR_hard_deliverables_completed=0), ensure_ascii=False))


if __name__ == '__main__':
    main()
