"""1234: same complete material, ideal thermal role, unchanged actual domain.

No extra physical preparation or energy-shell gate is granted by this check.
The comparison group is mathematical TR evidence, not an apparatus model.
Default is read-only comparison. --write-once creates the result exclusively.
--replay-history additionally replays 1166 through the current archive view.
"""
from __future__ import annotations
import argparse
from fractions import Fraction
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
TOL = 2e-10
ETA = np.diag([1., -1., -1., -1.])
SOURCES = [
    'archive_1106_1119/1112/check.py',
    'archive_1106_1119/1115/proof.md',
    'archive_1106_1119/1119/proof.md',
    'archive_1120_1158/1158/proof.md',
    'archive_1161_1200/1166/proof.md',
    'archive_1161_1200/1166/results.json',
    'archive_1161_1200/1172/proof.md',
    'archive_1201_1218/_shared/thermal_roles_adoption_after1218.md',
    'archive_1201_1218/1209/proof.md',
    'archive_1201_1218/1210/proof.md',
    'archive_1201_1218/1215/proof.md',
    'archive_1201_1218/1218/proof.md',
    'archive_1219_1233/1219/candidate_proof.md',
    'archive_1234_/1234/proof.md',
    'archive_1234_/1234/NEXT.md',
    'archive_1234_/1234/check.py',
]


def small(a, name, tol=TOL):
    value = float(np.linalg.norm(a))
    assert value < tol, (name, value)
    return value


def entropy_eigs(values):
    values = np.asarray(values, float)
    assert values.min(initial=0) > -TOL
    positive = values[values > 1e-15]
    return float(-np.dot(positive, np.log(positive)))


def entropy(a):
    return entropy_eigs(np.linalg.eigvalsh((a+a.conj().T)/2))


def trace_distance(a):
    return float(np.abs(np.linalg.eigvalsh((a+a.conj().T)/2)).sum()/2)


def fixture():
    spec = importlib.util.spec_from_file_location('r1112_for1234',
        BASE/'archive_1106_1119/1112/check.py')
    old = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(old)
    m = old.Material(1)
    I = np.eye(288, dtype=complex)
    K = m.hfk
    M = m.baseline*I+K
    eigs, V = np.linalg.eigh(K)
    spectrum_error = small(eigs-np.repeat(np.arange(4), 72), 'complete spectrum')
    projectors = []
    for j in range(4):
        P = I.copy()
        for k in range(4):
            if j != k:
                P = P@(K-k*I)/(j-k)
        projectors.append(P)
    r = sum(w*p for w, p in zip([8., 4., 2., 1.], projectors))/1080
    r_spec = (V*(2.**(-eigs)/np.sum(2.**(-eigs))))@V.conj().T
    thermal_error = small(r-r_spec, 'same-H polynomial/exponential')
    assert np.linalg.eigvalsh(r).min() > 0
    return m, M, r, dict(
        complete_dimension=288,
        baseline=m.baseline,
        eigenvalue_multiplicity=72,
        beta=math.log(2),
        spectrum_error=spectrum_error,
        thermal_identity_error=thermal_error,
        trace_error=small(np.trace(r)-1, 'normalization'),
        stationarity_error=small(M@r-r@M, 'same physical generator'),
        shell_probabilities=[float(np.trace(p@r).real) for p in projectors],
        exact_shell_probabilities=['8/15', '4/15', '2/15', '1/15'],
        mean_energy=float(np.trace(M@r).real),
        mean_energy_error=small(np.trace(M@r).real-m.baseline-11/15, 'same energy'),
        largest_eigenvalue=1/135,
        smallest_eigenvalue=1/1080,
    )


def ready_gap(m, M, r):
    delta = math.pi/4
    columns = np.eye(3, dtype=complex)[:, None, :]
    J = m.encode_joint(columns, delta).reshape(288, 3)
    isometry_error = small(J.conj().T@J-np.eye(3), 'old captured code')
    compression_error = small(J.conj().T@r@J-np.eye(3)/320,
                              'original code/Gibbs overlap')
    p = float(np.trace(J.conj().T@r@J).real)
    assert abs(p-3/320) < TOL
    assert Fraction(1)-3*Fraction(1, 135) == Fraction(44, 45)
    assert Fraction(1)-Fraction(3, 320) == Fraction(317, 320)
    ready = J@np.diag([.5, .3, .2])@J.conj().T
    distance = trace_distance(ready-r)
    assert distance >= 317/320-TOL
    return dict(
        actual_code_dimension=3,
        isometry_error=isometry_error,
        thermal_compression_error=compression_error,
        gibbs_probability_of_parent_code_effect=p,
        exact_parent_code_probability='3/320',
        all_original_code_states_trace_distance_lower_bound='317/320',
        spectrum_only_lower_bound='44/45',
        one_mixed_code_trace_distance=distance,
        no_new_code_effect_instrument=True,
        no_actual_gibbs_preparation_certificate=True,
    )


def comparison_group():
    # Full D=A_level(4) x degeneracy(72); active block a=0 only.
    # All other 71 blocks are retained explicitly by their identical spectrum.
    q = np.array([8., 4., 2., 1.])/15
    block = np.diag(q)/72
    a = np.array([1., 1j], complex)/math.sqrt(2)
    b = np.outer(a, a.conj())
    theta = .37
    W = np.eye(8, dtype=complex)
    i, j = 1, 2  # |0,1> and |1,0>, same total energy.
    W[np.ix_([i, j], [i, j])] = [[math.cos(theta), -math.sin(theta)],
                                [math.sin(theta), math.cos(theta)]]
    K = np.diag(np.repeat(np.arange(4), 2)+np.tile(np.arange(2), 4))
    unitary_error = small(W.conj().T@W-np.eye(8), 'full block unitary')
    conservation_error = small(W@K-K@W, 'exact mathematical energy shell')
    incoming = np.kron(block, b)
    outgoing = W@incoming@W.conj().T
    rho = outgoing.reshape(4, 2, 4, 2)
    d0 = np.einsum('aibi->ab', rho)
    a0 = np.einsum('aiaj->ij', rho)
    aout = a0+(71/72)*b
    sin = entropy_eigs(np.tile(q/72, 72))
    sout = entropy(d0)+71*entropy(block)
    saout = entropy(aout)
    sjointout = entropy(outgoing)+71*entropy(incoming)
    mutual = sout+saout-sjointout
    de_d = float(np.dot(np.arange(4), np.diag(d0-block).real))
    de_a = float((aout-b)[1, 1].real)
    assert abs(de_d+de_a) < TOL
    drel = float(-sout-np.dot(np.diag(d0).real+71*q/72, np.log(q/72)))
    identity_rhs = -math.log(2)*de_a+saout-mutual
    identity_error = small(drel-identity_rhs, 'full relative entropy identity')
    dtr = trace_distance(d0-block)
    delta = abs(de_a)
    kappa = abs(saout)
    epsilon = trace_distance(aout-b)
    bound = math.sqrt((math.log(2)*delta+kappa)/2)
    assert drel > 0 and dtr > 0 and dtr <= bound+TOL
    assert mutual > -TOL
    return dict(
        full_material_dimension=288,
        full_auxiliary_dimension=2,
        retained_inactive_degeneracy_blocks=71,
        active_energy_shell=[1, 2],
        unitary_error=unitary_error,
        energy_shell_error=conservation_error,
        total_entropy_error=small(sjointout-sin, 'whole entropy'),
        delta_material_energy=de_d,
        delta_auxiliary_energy=de_a,
        auxiliary_entropy_change=saout,
        mutual_information=mutual,
        relative_entropy=drel,
        full_relative_entropy_identity_error=identity_error,
        three_budgets=dict(delta=delta, kappa=kappa, epsilon=epsilon),
        material_trace_distance=dtr,
        pinsker_uniform_bound=bound,
        mathematical_comparison_only=True,
        initial_auxiliary_purity_and_coherence_retained=True,
    )


def actual_geometry(m, M, r):
    columns = np.eye(3, dtype=complex)[:, None, :]
    delta = math.pi/4
    J = m.encode_joint(columns, delta).reshape(288, 3)
    vals, V = np.linalg.eigh(M)
    mask = np.tile(np.repeat(np.arange(8) % 2, 9), 4)
    rows = []
    for s, u in [(.31, np.array([1., 0., 0., 0.])),
                 (.77, np.array([1.25, .75, 0., 0.]))]:
        P = [u[i]*M for i in range(4)]
        dz = s*u
        interval = float(dz@ETA@dz)
        expected = math.sin((delta+s)/2)**6
        U = (V*np.exp(-1j*s*vals))@V.conj().T
        output = U@J
        effect = output.conj().T@(mask[:, None]*output)
        rows.append(dict(
            proper_duration=s,
            complete_clock_readout=expected,
            full_input_readout_error=small(effect-expected*np.eye(3), 'old Mr effect'),
            interval_error=small(interval-s*s, 'same eta interval'),
            proper_energy_error=small(sum(ETA[i, i]*u[i]*P[i] for i in range(4))-M,
                                      'same observer and energy'),
            full_phase_error=small(sum(ETA[i, i]*dz[i]*P[i] for i in range(4))-s*M,
                                   'complete material phase'),
            equilibrium_clock_stationarity_error=small(U@r@U.conj().T-r,
                                                      'ideal reference is stationary'),
        ))
    return dict(
        rows=rows,
        actual_preparations_operations_and_futures_unchanged=True,
        all_actual_domain_statement_proved_by_1166_structural_induction=True,
        finite_rows_do_not_prove_the_universal_statement=True,
        universal_curvature_or_einstein_claim=False,
    )


def compare(actual, saved):
    if isinstance(saved, dict):
        assert actual.keys() == saved.keys()
        for k in saved:
            compare(actual[k], saved[k])
    elif isinstance(saved, list):
        assert len(actual) == len(saved)
        for a, b in zip(actual, saved):
            compare(a, b)
    elif isinstance(saved, float):
        assert math.isclose(actual, saved, rel_tol=2e-8, abs_tol=3e-10), (actual, saved)
    else:
        assert actual == saved, (actual, saved)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write-once', action='store_true')
    parser.add_argument('--replay-history', action='store_true')
    args = parser.parse_args()
    if args.replay_history:
        runner = ROOT/'scripts/run_research_gr_response_active.py'
        command = [sys.executable, '-B', '-X', 'utf8', str(runner), '--script',
                   str(BASE/'archive_1161_1200/1166/check.py')]
        replay = subprocess.run(command, cwd=ROOT, text=True, capture_output=True,
                                encoding='utf-8', timeout=180)
        assert replay.returncode == 0, (replay.stdout, replay.stderr)
        assert 'READ_ONLY_MATCH' in replay.stdout, replay.stdout
        print(json.dumps(dict(history_replay='PASS', output=replay.stdout.strip()),
                         ensure_ascii=False))
    m, M, r, first = fixture()
    results = dict(
        status='PASS_IN_IDEAL_ROLE_EXTENSION_SCOPE',
        round=1234,
        scientific_check_groups=3,
        inherited_interface_checks=1,
        new_axioms_adopted=0,
        GR_hard_deliverables_completed=0,
        source_sha256={p: hashlib.sha256((BASE/p).read_bytes()).hexdigest() for p in SOURCES},
        same_complete_material_reference=first,
        actual_ready_state_distinction=ready_gap(m, M, r),
        mathematical_TR_budget_comparison=comparison_group(),
        unchanged_actual_geometric_interfaces=actual_geometry(m, M, r),
        scope=dict(
            ideal_reference_department_nonempty=True,
            new_actual_thermal_preparation=False,
            stronger_physical_expansion_proved=False,
            TW_optional_not_universal=True,
            TC_not_universal_q_positive=True,
            entire_GR_goal_completed=False,
        ),
    )
    path = HERE/'results.json'
    if args.write_once:
        with path.open('x', encoding='utf-8', newline='\n') as f:
            json.dump(results, f, ensure_ascii=False, indent=2)
            f.write('\n')
        mode = 'WRITE_ONCE'
    else:
        compare(results, json.loads(path.read_text(encoding='utf-8')))
        mode = 'READ_ONLY_MATCH'
    print(json.dumps(dict(status=results['status'], mode=mode,
                         groups=3, inherited_checks=1, sources=len(SOURCES),
                         complete_dimension=288,
                         actual_code_thermal_distance_lower_bound='317/320',
                         GR_hard_deliverables_completed=0), ensure_ascii=False))


if __name__ == '__main__':
    main()
