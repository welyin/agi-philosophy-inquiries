"""1288 independent CAR / finite cyclic-parent audit.

No imports of author or historical calculation implementations. Exact integer
matrix checks carry the finite diagnostics; all-N and all-time statements are
proved in the separately reviewed author text. --write is exclusive creation.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in HERE.parents if (p / 'research_cognition_physics').is_dir() and (p / 'scripts').is_dir())
TARGET = HERE / 'independent_results.json'
OLD = ROOT / 'research_cognition_physics/archive_956_989'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def local_car():
    basis = tuple(n for n in range(16) if n.bit_count() == 2)
    index = {n: i for i, n in enumerate(basis)}
    q = np.array([(n & 3).bit_count() - 1 for n in basis], dtype=np.int64)
    raising = np.zeros((6, 6), dtype=np.int64)
    # Apply c_i^dagger c_(i+2), computing both fermionic signs explicitly.
    for j, bits in enumerate(basis):
        for i in (0, 1):
            destroy = i + 2
            if not ((bits >> destroy) & 1):
                continue
            sign = (-1) ** (bits & ((1 << destroy) - 1)).bit_count()
            rest = bits ^ (1 << destroy)
            if (rest >> i) & 1:
                continue
            sign *= (-1) ** (rest & ((1 << i) - 1)).bit_count()
            raising[index[rest | (1 << i)], j] += sign
    assert np.array_equal((q[:, None] - q[None, :]) * raising, raising)
    assert not np.any(raising[:, q == 1])
    assert not np.any(raising.T[:, q == -1])
    return basis, q, raising


def material():
    basis, q, raising = local_car()
    ident = np.eye(6, dtype=np.int64)
    qa, qb = np.repeat(q, 6), np.tile(q, 6)
    ta, tb = np.kron(raising, ident), np.kron(ident, raising)
    sv = -(ta + ta.T + tb + tb.T)
    su = np.diag(qa * qa + qb * qb)
    sk = np.diag(qa * qb)
    h10 = sv + 10 * su + 2 * sk
    return basis, qa, qb, ta, tb, h10, (su, sk, sv)


def parent(n, qa, qb, ta, tb):
    # Residue basis 0,...,n-1. Principal representatives agree with +/-1 for n>=3.
    e = np.array([j if j <= n // 2 else j - n for j in range(n)], dtype=np.int64)
    dim = 36 * n * n
    ea = np.tile(np.repeat(e, n), 36)
    eb = np.tile(np.tile(e, n), 36)
    qaa, qbb = np.repeat(qa, n*n), np.repeat(qb, n*n)
    sv = np.zeros((dim, dim), dtype=np.int64)
    for raising, leg in ((ta, 0), (tb, 1)):
        rows, cols = np.nonzero(raising)
        for s, t in zip(rows, cols):
            for a in range(n):
                for b in range(n):
                    aa, bb = ((a + 1) % n, b) if leg == 0 else (a, (b + 1) % n)
                    src, dst = t*n*n + a*n + b, s*n*n + aa*n + bb
                    sv[dst, src] -= raising[s, t]
                    sv[src, dst] -= raising[s, t]
    su, sk = np.diag(ea*ea + eb*eb), np.diag(ea*eb)
    h10 = sv + 10*su + 2*sk
    jidx = np.arange(36)*n*n + (qa % n)*n + qb % n
    ga, gb = (ea-qaa) % n, (eb-qbb) % n
    assert np.array_equal(np.flatnonzero((ga == 0) & (gb == 0)), np.sort(jidx))
    rr, cc = np.nonzero(h10)
    assert np.all(ga[rr] == ga[cc]) and np.all(gb[rr] == gb[cc])
    # It is genuinely a Z_n parent: continuous E-Q is not a full-kinematic symmetry.
    continuous_defect = max(int(np.max(np.abs((ea-qaa)[rr]-(ea-qaa)[cc]))),
                            int(np.max(np.abs((eb-qbb)[rr]-(eb-qbb)[cc]))))
    assert continuous_defect == n
    return h10, (su, sk, sv), jidx, continuous_defect


def transport_columns(operator, jidx, oldoperator):
    actual = operator[:, jidx]
    expected = np.zeros_like(actual)
    expected[jidx, :] = oldoperator
    assert np.array_equal(actual, expected)
    return int(np.max(np.abs(actual - expected)))


def dressed_apply(m, amplitudes, n, qa, qb):
    """Controlled cyclic translations; no compression to the physical sector."""
    inp = amplitudes.reshape(36, n, n, -1)
    out = np.zeros_like(inp, dtype=np.complex128)
    for s, t in zip(*np.nonzero(m)):
        shifted = np.roll(inp[t], (int(qa[s]-qa[t]), int(qb[s]-qb[t])), axis=(0, 1))
        out[s] += m[s, t] * shifted
    return out.reshape(amplitudes.shape)


def operator_tests(n, qa, qb, jidx):
    # Two complete Kraus branches, each a different local permutation after
    # a diagonal filter. These are tests of the old representation theorem,
    # not a claim that new arbitrary operations have been physically supplied.
    ids = np.arange(36)
    perm0 = ((ids // 6 + 1) % 6)*6 + ids % 6
    perm1 = (ids // 6)*6 + 5-ids % 6
    w0 = np.where(ids % 2 == 0, 3, 4)
    w1 = np.where(ids % 2 == 0, 4, 3)
    k0, k1 = np.zeros((36, 36)), np.zeros((36, 36))
    k0[perm0, ids], k1[perm1, ids] = w0/5, w1/5
    assert np.array_equal(w0*w0+w1*w1, np.full(36, 25))
    # Generic deterministic, entangled material-reference vector (R dimension 4).
    i, j = np.indices((36, 4))
    psi = ((3*i+2*j) % 7-3) + 1j*((i-j) % 5-2)
    psi = psi / np.linalg.norm(psi)
    assert np.linalg.matrix_rank(psi) == 4
    lifted = np.zeros((36*n*n, 4), dtype=complex)
    lifted[jidx] = psi
    branch_error, total = 0.0, 0.0
    for k in (k0, k1):
        actual = dressed_apply(k, lifted, n, qa, qb)
        expected = np.zeros_like(actual)
        expected[jidx] = k @ psi
        branch_error = max(branch_error, float(np.max(np.abs(actual-expected))))
        total += float(np.linalg.norm(actual)**2)
    assert branch_error < 1e-14 and abs(total-1) < 1e-14
    # A vector outside Gauss=0 also verifies product and adjoint identities.
    a, b = np.indices((36*n*n, 2))
    trial = ((a+3*b) % 11-5) + 1j*((2*a+b) % 7-3)
    trial = trial / np.linalg.norm(trial)
    left = dressed_apply(k0 @ k1, trial, n, qa, qb)
    right = dressed_apply(k0, dressed_apply(k1, trial, n, qa, qb), n, qa, qb)
    product_error = float(np.max(np.abs(left-right)))
    completeness = dressed_apply(k0.T, dressed_apply(k0, trial, n, qa, qb), n, qa, qb)
    completeness += dressed_apply(k1.T, dressed_apply(k1, trial, n, qa, qb), n, qa, qb)
    complete_error = float(np.max(np.abs(completeness-trial)))
    assert product_error < 1e-14 and complete_error < 1e-14
    return dict(reference_dimension=4, actual_branch_amplitude_max_error=branch_error,
                branch_probability_sum=total, full_kinematic_product_max_error=product_error,
                full_kinematic_kraus_completeness_max_error=complete_error)


def unitary(h, time):
    vals, vecs = np.linalg.eigh(h.astype(float))
    return (vecs*np.exp(-1j*time*vals)) @ vecs.T


def calculate():
    basis, qa, qb, ta, tb, h10, sources = material()
    rows = []
    for n in (3, 4, 5):
        hn, sn, jidx, defect = parent(n, qa, qb, ta, tb)
        herr = transport_columns(hn, jidx, h10)
        serr = [transport_columns(a, jidx, b) for a, b in zip(sn, sources)]
        tests = operator_tests(n, qa, qb, jidx)
        row = dict(N=n, kinematic_dimension=36*n*n, gauss_dimension=len(jidx),
                   integer_10H_intertwining_error=herr,
                   exact_source_errors=dict(zip(('U', 'kappa', 'v'), serr)),
                   nonzero_continuous_gauss_charge_jump=defect, **tests)
        if n == 3:
            # One small, unrelated diagnostic time. The original long record
            # window is inherited analytically, never refitted or recomputed.
            up, um = unitary(hn/10, 3/7), unitary(h10/10, 3/7)
            expected = np.zeros((len(hn), 36), dtype=complex)
            expected[jidx] = um
            err = float(np.max(np.abs(up[:, jidx]-expected)))
            assert err < 1e-12
            row['one_time_full_parent_propagator_max_error'] = err
        rows.append(row)
    # N=2 collision: (+1,+1) and (-1,+1) have the identical two link residues.
    # No field-only function, not just the chosen quadratic polynomial, can
    # match both 2U+kappa and 2U-kappa closer than |kappa|.
    kappa = F(1, 5)
    desired = (F(2)+kappa, F(2)-kappa)
    minimax = (max(desired)-min(desired))/2
    assert minimax == kappa
    source_minimax = (F(1)-F(-1))/2
    assert source_minimax == 1
    old = json.loads((OLD/'958/capacitive_material_write_results.json').read_text('utf-8-sig'))
    assert old['parameters']['U'] == 1 and old['parameters']['hopping'] == .1 and old['parameters']['kappa'] == .2
    paths = [Path(__file__), OLD/'research_note_959.md', OLD/'research_note_958.md',
             OLD/'958/capacitive_material_write_results.json', OLD/'959/gauss_material_parent_results.json']
    return dict(round=1288, status='PASS_WITH_SCOPE', independent_groups=1,
        implementation='independent occupation-bit CAR, integer matrices, cyclic shifts; no author imports',
        local_two_particle_basis=list(basis), material_dimension=36,
        parameters=dict(U='1', v='1/10', kappa='1/5'), cyclic_parent_cases=rows,
        N2_obstruction=dict(link_residues_for_both_states=[1, 1],
            required_electric_values=[str(x) for x in desired],
            all_field_only_energy_uniform_error_lower=str(minimax),
            all_field_only_kappa_source_uniform_error_lower=str(source_minimax),
            scope='fixed nonzero-kappa pure-link electrostatic dictionary; not all Z2 theories'),
        old_958_budget_inherited_without_new_record_simulation={k: old['exact_interface'][k] for k in
            ('all_time_unknown_input_error_upper', 'controlled_phase_time', 'read_window_half_width', 'analytic_receiver_contrast_lower')},
        scope=dict(all_N_ge_3_by_analytic_argument=True, all_times_by_exact_intertwiner=True,
            arbitrary_passive_partner=True, entire_fixed_number_Gauss_sector=True,
            original_three_parameter_sources_only=True, actual_instruments_only_if_previously_licensed=True,
            arbitrary_new_instrument_permission=False, full_kinematic_U1_equivalence=False,
            full_3plus1_stress_tensor=False, propagating_photon_or_Lorentz_recovery=False,
            entire_cognitive_framework_countermodel=False, group_uniqueness_selected=False),
        source_sha256={p.relative_to(ROOT).as_posix(): digest(p) for p in paths})


def same(saved, current, location='root'):
    if isinstance(current, dict):
        assert saved.keys() == current.keys(), location
        for k in current: same(saved[k], current[k], location+'.'+k)
    elif isinstance(current, list):
        assert len(saved) == len(current), location
        for i, (a,b) in enumerate(zip(saved,current)): same(a,b,location+f'[{i}]')
    elif isinstance(current, float):
        assert math.isclose(saved, current, rel_tol=2e-11, abs_tol=1e-12), location
    else:
        assert saved == current, location


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = calculate()
    if args.write:
        with TARGET.open('x', encoding='utf-8', newline='\n') as f:
            json.dump(result, f, ensure_ascii=False, indent=2)
            f.write('\n')
        mode = 'EXCLUSIVE_WRITE'
    else:
        same(json.loads(TARGET.read_text('utf-8')), result)
        mode = 'READ_ONLY_MATCH'
    print(json.dumps(dict(status=result['status'], mode=mode,
        tested_N=[r['N'] for r in result['cyclic_parent_cases']],
        N2_energy_lower=result['N2_obstruction']['all_field_only_energy_uniform_error_lower']), ensure_ascii=False))
