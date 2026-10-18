"""1289: finite checks of local composite-vector overlap classification.

The all-N and all-smooth-local-jet statements are analytic in proof.md.
No quantum detector, preparation, source budget or nonperturbative error is
inferred from these finite linear-algebra checks. Default execution is read-only;
--write exclusively creates a previously absent results.json.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parents[1]
ROOT = RESEARCH.parent
TARGET = HERE / "results.json"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(error, tolerance=3e-10):
    error = float(error)
    assert error < tolerance, (error, tolerance)
    return error


def rank(matrix, tolerance=2e-10):
    singular = np.linalg.svd(np.asarray(matrix), compute_uv=False)
    if not len(singular):
        return 0
    return int(np.count_nonzero(singular > tolerance * max(1., singular[0])))


def generators(n):
    """Tr(Ta Tb)=delta_ab/2; N=2 is the ordinary Pauli/2 basis."""
    out = []
    for i in range(n):
        for j in range(i + 1, n):
            symmetric = np.zeros((n, n), complex)
            symmetric[i, j] = symmetric[j, i] = .5
            antisymmetric = np.zeros((n, n), complex)
            antisymmetric[i, j], antisymmetric[j, i] = -.5j, .5j
            out.extend([symmetric, antisymmetric])
    for k in range(1, n):
        diagonal = np.zeros(n)
        diagonal[:k], diagonal[k] = 1, -k
        out.append(np.diag(diagonal) / math.sqrt(2 * k * (k + 1)))
    return np.asarray(out)


def matrix_units(n):
    return np.eye(n * n, dtype=complex).reshape(n * n, n, n)


def invariant_tensor_dimensions(ts):
    n = ts.shape[1]
    units = matrix_units(n)
    mixed = np.column_stack([
        np.concatenate([(t @ b - b @ t).ravel() for t in ts]) for b in units])
    holomorphic = np.column_stack([
        np.concatenate([(t.T @ b + b @ t).ravel() for t in ts]) for b in units])
    return n * n - rank(mixed), n * n - rank(holomorphic)


def stabilizer_dimension(ts):
    n = ts.shape[1]
    if n == 2:
        return len(ts), 0.
    hs = []
    for small in generators(n - 1):
        h = np.zeros((n, n), complex)
        h[:n-1, :n-1] = small
        hs.append(h)
    columns = []
    for t in ts:
        column = np.concatenate([(h @ t - t @ h).ravel() for h in hs])
        columns.append(np.concatenate([column.real, column.imag]))
    linear_map = np.column_stack(columns)
    residual = require(np.linalg.norm(linear_map[:, -1]))
    return len(ts) - rank(linear_map), residual


EPSILON = np.array([[0., 1.], [-1., 0.]], complex)


def data(ts, nvec, g, v):
    vectors = np.array([t @ nvec for t in ts])
    mass = g*g*v*v * (vectors.conj() @ vectors.T).real
    neutral = g*v*v * np.array([np.vdot(nvec, t @ nvec).real for t in ts])
    if len(nvec) == 2:
        hol = -.5j*g*v*v * np.array([nvec.T @ EPSILON @ t @ nvec for t in ts])
        overlap = np.array([2*hol.imag, 2*hol.real, neutral])
    else:
        overlap = neutral[None, :]
    return mass, overlap, neutral


def change_representative(ts, nvec, mass, overlap, g, v):
    """A gauge change of description, not an applied physical channel."""
    generator = sum((.11 / (i + 1)) * t for i, t in enumerate(ts))
    ev, eigenvectors = np.linalg.eigh(generator)
    u = (eigenvectors * np.exp(1j * ev)) @ eigenvectors.conj().T
    rotated = np.array([u @ t @ u.conj().T for t in ts])
    adjoint = np.array([[2*np.trace(a @ b).real for b in rotated] for a in ts])
    mp, cp, _ = data(ts, u @ nvec, g, v)
    checks = dict(
        unitary=require(np.linalg.norm(u.conj().T @ u - np.eye(len(nvec)))),
        determinant_one=require(abs(np.linalg.det(u)-1)),
        adjoint_orthogonal=require(np.linalg.norm(adjoint @ adjoint.T - np.eye(len(ts)))),
        mass_transport=require(np.linalg.norm(mp-adjoint @ mass @ adjoint.T)),
        overlap_transport=require(np.linalg.norm(cp-overlap @ adjoint.T)),
    )
    if len(nvec) == 2:
        checks["epsilon_invariant"] = require(np.linalg.norm(u.T @ EPSILON @ u - EPSILON))
    return checks


def check_case(n, g, v):
    ts = generators(n)
    nvec = np.zeros(n, complex)
    nvec[-1] = 1
    mass, c, neutral = data(ts, nvec, g, v)
    dof = n*n-1
    gram = np.array([[np.trace(a @ b) for b in ts] for a in ts])
    normalization = require(np.linalg.norm(gram-.5*np.eye(dof)))
    m_off, m_diag = g*g*v*v/4, g*g*v*v*(n-1)/(2*n)
    expected = np.sort([0.] * (n*n-2*n) + [m_off]*(2*n-2) + [m_diag])
    values, vectors = np.linalg.eigh(mass)
    mass_error = require(np.linalg.norm(values-expected))
    assert rank(mass) == 2*n-1
    expected_rank = 3 if n == 2 else 1
    assert rank(c) == expected_rank
    assert rank(neutral[None, :]) == 1
    residue_rows = []
    pole_groups = [(m_off, 3, 3)] if n == 2 else [(0., n*n-2*n, 0), (m_off, 2*n-2, 0), (m_diag, 1, 1)]
    for pole, multiplicity, expected_residue_rank in pole_groups:
        selected = np.abs(values-pole) < 2e-10
        assert int(np.count_nonzero(selected)) == multiplicity
        projector = vectors[:, selected] @ vectors[:, selected].conj().T
        residue = c @ projector @ c.T
        residue_rank = rank(residue)
        assert residue_rank == expected_residue_rank
        residue_rows.append(dict(mass_squared=float(pole), bare_multiplicity=multiplicity,
                                 residue_rank=residue_rank, residue_trace=float(np.trace(residue).real)))
    propagator_residuals = []
    for p2 in (.3, 1.1, 4.):
        actual = c @ np.linalg.solve(p2*np.eye(dof)+mass, c.T)
        predicted = (g*g*v**4/4)/(p2+m_off)*np.eye(3) if n == 2 else np.array([
            [(g*g*v**4*(n-1)/(2*n))/(p2+m_diag)]])
        propagator_residuals.append(require(np.linalg.norm(actual-predicted)))
    mixed, hol = invariant_tensor_dimensions(ts)
    assert mixed == 1 and hol == int(n == 2)
    center_factor = np.exp(4j*np.pi/n)
    if n >= 3:
        assert abs(center_factor-1) > .5
    stabilizer, stab_error = stabilizer_dimension(ts)
    assert stabilizer == expected_rank
    assert rank(c) == stabilizer
    c_gram = c @ c.T
    expected_gram = (g*g*v**4/4)*np.eye(3) if n == 2 else np.array([[g*g*v**4*(n-1)/(2*n)]])
    overlap_residual = require(np.linalg.norm(c_gram-expected_gram))
    gauge_checks = change_representative(ts, nvec, mass, c, g, v)
    alpha = .37
    mass_phase, c_phase, neutral_phase = data(ts, np.exp(1j*alpha)*nvec, g, v)
    phase_checks = dict(mass_invariant=require(np.linalg.norm(mass_phase-mass)),
                        neutral_invariant=require(np.linalg.norm(neutral_phase-neutral)))
    if n == 2:
        rotation = np.array([[math.cos(2*alpha), math.sin(2*alpha), 0.],
                             [-math.sin(2*alpha), math.cos(2*alpha), 0.], [0., 0., 1.]])
        phase_checks["charged_pair_covariant"] = require(np.linalg.norm(c_phase-rotation @ c))
        phase_checks["single_current_neutral_subspace_dimension"] = 3-rank(rotation-np.eye(3))
        assert phase_checks["single_current_neutral_subspace_dimension"] == 1
        require(np.linalg.norm(c-(g*v*v/2)*np.diag([1., 1., -1.])))
    else:
        phase_checks["single_current_neutral_subspace_dimension"] = 1
    return dict(N=n, g=g, v=v, number_of_generators=dof,
                generator_normalization_residual=normalization,
                mass_eigenvalues=[float(x) for x in values], mass_spectrum_residual=mass_error,
                bare_mass_rank=rank(mass), complete_menu_overlap_rank=rank(c),
                global_neutral_single_current_rank=rank(neutral[None, :]),
                mixed_invariant_tensor_complex_dimension=mixed,
                holomorphic_invariant_tensor_complex_dimension=hol,
                center_holomorphic_factor=[float(center_factor.real), float(center_factor.imag)],
                residual_stabilizer_centralizer_dimension=stabilizer,
                diagonal_centralizer_residual=stab_error,
                overlap_gram_residual=overlap_residual,
                transverse_propagator_residuals=propagator_residuals,
                distinct_pole_groups=residue_rows,
                gauge_representative_change=gauge_checks, global_phase=phase_checks)


def run():
    rows = [check_case(n, g, v) for n in (2, 3, 4, 5) for g, v in ((1., 1.), (.7, 1.4), (1.7, .6))]
    zero_checks = []
    for n in (2, 3, 4, 5):
        ts = generators(n)
        nvec = np.eye(n, dtype=complex)[:, -1]
        for g, v in ((0., 1.), (1., 0.)):
            mass, c, _ = data(ts, nvec, g, v)
            assert np.linalg.norm(mass) == np.linalg.norm(c) == 0
            zero_checks.append(dict(N=n, g=g, v=v, bare_mass_rank=rank(mass), overlap_rank=rank(c)))
    source_files = [HERE / "proof.md", Path(__file__), HERE.parent / "_shared/goal.md",
        RESEARCH / "archive_554_584/research_note_568.md",
        RESEARCH / "archive_554_584/research_note_569.md",
        RESEARCH / "archive_923_934/research_note_934.md",
        RESEARCH / "archive_1009_1043/research_note_1016.md",
        RESEARCH / "archive_1009_1043/research_note_1030.md",
        RESEARCH / "archive_1009_1043/research_note_1039.md",
        RESEARCH / "archive_1281_1287/research_note_1284.md"]
    return dict(round=1289, date="2026-10-18", all_scientific_checks_passed=True,
                finite_cases=rows, zero_coupling_or_amplitude_controls=zero_checks,
                scope=dict(
                    all_N_classification_proved_analytically=True,
                    all_smooth_local_finite_jets_bound_proved_at_trivial_background=True,
                    complete_degree2_one_covariant_derivative_menu=True,
                    additional_noninvariant_color_spurion=False,
                    rank_counts_internal_vector_directions_not_spatial_dimension=True,
                    distinct_mass_poles_not_identified_with_channel_rank=True,
                    global_charge_not_gauge_charge=True,
                    gauge_representative_change_is_not_physical_control=True,
                    high_order_or_nonperturbative_particles_excluded=False,
                    any_SU_N_group_excluded=False,
                    SU2_uniquely_selected_by_cognitive_principles=False,
                    actual_unknown_physical_code_certified=False,
                    actual_composite_readout_permission_granted=False,
                    source_complete_probe_built=False,
                    reference_preparation_or_transport_built=False,
                    new_GR_source_or_geometry_validation=False,
                    nonperturbative_pole_error_bound=False,
                    new_axioms=0, new_common_devices=0, stage_goal_completed=False),
                source_hashes={str(p.relative_to(ROOT)).replace('\\', '/'):sha(p) for p in source_files})


def compare(actual, saved, location="root"):
    if isinstance(actual, dict):
        assert set(actual) == set(saved), location
        for key, value in actual.items():
            compare(value, saved[key], location + "." + key)
    elif isinstance(actual, list):
        assert len(actual) == len(saved), location
        for i, (a, b) in enumerate(zip(actual, saved)):
            compare(a, b, location + f"[{i}]")
    elif isinstance(actual, float):
        assert math.isclose(actual, saved, rel_tol=2e-11, abs_tol=2e-10), (location, actual, saved)
    else:
        assert actual == saved, (location, actual, saved)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    if args.write and TARGET.exists():
        raise FileExistsError("results.json exists; refusing to overwrite scientific output")
    result = run()
    if args.write:
        with TARGET.open("x", encoding="utf-8") as handle:
            json.dump(result, handle, ensure_ascii=False, indent=2)
            handle.write("\n")
        status = "CREATED"
    else:
        compare(result, json.loads(TARGET.read_text(encoding="utf-8-sig")))
        status = "MATCH"
    print(json.dumps(dict(status=status, round=1289, all_scientific_checks_passed=True,
        cases=len(result["finite_cases"]),
        bare_vs_menu=[dict(N=r["N"], mass_rank=r["bare_mass_rank"], menu_rank=r["complete_menu_overlap_rank"],
            neutral_rank=r["global_neutral_single_current_rank"])
            for r in result["finite_cases"] if r["g"] == r["v"] == 1.]), ensure_ascii=False))
