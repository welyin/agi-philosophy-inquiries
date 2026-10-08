"""Round 1011: physical Yukawa matrices constrain generation charges.

Finite calibration of the analytic classification in research_note_1011.md.
No fitted CKM data, apparatus model, anomaly approximation, or new dynamics.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
OUT = HERE / "flavor_charge_selection_results.json"


def comm(a, b):
    return a @ b - b @ a


def op(a):
    return float(np.linalg.norm(a, ord=2))


def psqrt(a):
    vals, vecs = np.linalg.eigh(a)
    assert min(vals) > 0
    return (vecs * np.sqrt(vals)) @ vecs.conj().T


def hermitian_basis(n):
    basis = []
    for i in range(n):
        a = np.zeros((n, n), dtype=complex)
        a[i, i] = 1
        basis.append(a)
    for i in range(n):
        for j in range(i + 1, n):
            a = np.zeros((n, n), dtype=complex)
            a[i, j] = a[j, i] = 1 / math.sqrt(2)
            basis.append(a)
            a = np.zeros((n, n), dtype=complex)
            a[i, j], a[j, i] = 1j / math.sqrt(2), -1j / math.sqrt(2)
            basis.append(a)
    return basis


def commutant_dimension(*matrices):
    cols = []
    for b in hermitian_basis(len(matrices[0])):
        vals = np.concatenate([comm(b, a).ravel() for a in matrices])
        cols.append(np.concatenate([vals.real, vals.imag]))
    svals = np.linalg.svd(np.array(cols).T, compute_uv=False)
    tol = 1e-10 * max(1., float(svals[0]))
    return int(len(svals) - np.count_nonzero(svals > tol))


def components(a, tol=1e-12):
    # Used only for explicit finite matrices; exact graph criterion is proved
    # analytically in the note, not inferred from this numerical threshold.
    unseen, out = set(range(len(a))), []
    while unseen:
        todo = [min(unseen)]
        found = set()
        while todo:
            i = todo.pop()
            if i in found:
                continue
            found.add(i)
            unseen.discard(i)
            todo.extend(j for j in unseen if abs(a[i, j]) > tol)
        out.append(sorted(found))
    return out


def matmul(a, b):
    return [[sum(x*y for x,y in zip(row, col)) for col in zip(*b)] for row in a]


def det3(a):
    return (a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])
            -a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])
            +a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]))


def fraction_rank(rows):
    a = [[F(x) for x in row] for row in rows]
    if not a:
        return 0
    r = 0
    for c in range(len(a[0])):
        pivot = next((j for j in range(r, len(a)) if a[j][c]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        v = a[r][c]
        a[r] = [x/v for x in a[r]]
        for j in range(len(a)):
            if j != r:
                v = a[j][c]
                a[j] = [x-v*y for x,y in zip(a[j], a[r])]
        r += 1
        if r == len(a):
            break
    return r


def charge_lift(ql, yu, yd, h):
    eye = np.eye(len(ql))
    qu = np.linalg.solve(yu, ql @ yu) + h*eye
    qd = np.linalg.solve(yd, ql @ yd) - h*eye
    return qu, qd


def ckm(t12, t13, t23, delta):
    c12,s12,c13,s13,c23,s23 = (
        math.cos(t12),math.sin(t12),math.cos(t13),math.sin(t13),
        math.cos(t23),math.sin(t23))
    e = np.exp(1j*delta)
    return np.array([
        [c12*c13, s12*c13, s13/e],
        [-s12*c23-c12*s23*s13*e,
         c12*c23-s12*s23*s13*e, s23*c13],
        [s12*s23-c12*c23*s13*e,
         -c12*s23-s12*c23*s13*e, c23*c13],
    ])


def compare(a, b):
    if isinstance(a, dict):
        assert a.keys() == b.keys()
        for key in a:
            compare(a[key], b[key])
    elif isinstance(a, list):
        assert len(a) == len(b)
        for x, y in zip(a, b):
            compare(x, y)
    elif isinstance(a, float):
        assert math.isclose(a, b, rel_tol=2e-9, abs_tol=2e-11), (a, b)
    else:
        assert a == b, (a, b)


def run():
    hu = np.diag([1., 4., 9.])
    yu = np.diag([1., 2., 3.])
    h = .5
    hd_real = np.array([[4., 1., 1.], [1., 5., 1.], [1., 1., 6.]])
    hd_complex = np.array([[4., 1., .5j], [1., 5., 1.], [-.5j, 1., 6.]])
    hd_split = np.array([[4., 1., 0.], [1., 5., 0.], [0., 0., 6.]])
    hd_chain = np.array([[4., 1., 0.], [1., 5., 1.], [0., 1., 6.]])
    cases = []
    max_intertwiner = max_hermitian = max_basis_covariance = 0.
    rng = np.random.default_rng(1011)
    for label, hd, qs in [
        ("connected_real_cp_conserving", hd_real, [1/6]*3),
        ("connected_complex", hd_complex, [1/6]*3),
        ("two_components", hd_split, [1., 1., 3.]),
        ("three_components", np.diag([4., 5., 6.]), [1., 2., 3.]),
    ]:
        yd, ql = psqrt(hd), np.diag(qs)
        qu, qd = charge_lift(ql, yu, yd, h)
        errors = [op(ql@yu-yu@(qu-h*np.eye(3))),
                  op(ql@yd-yd@(qd+h*np.eye(3)))]
        herm = max(op(qu-qu.conj().T), op(qd-qd.conj().T))
        ncomp = len(components(hd))
        dim = commutant_dimension(hu, hd)
        assert dim == ncomp
        assert max(errors) < 1e-12 and herm < 1e-12
        max_intertwiner = max(max_intertwiner, *errors)
        max_hermitian = max(max_hermitian, herm)
        for _ in range(4):
            ws = []
            for j in range(3):
                w, _ = np.linalg.qr(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))
                ws.append(w)
            w, ru, rd = ws
            yu1, yd1 = w@yu@ru.conj().T, w@yd@rd.conj().T
            ql1 = w@ql@w.conj().T
            qu1, qd1 = charge_lift(ql1, yu1, yd1, h)
            err = max(op(qu1-ru@qu@ru.conj().T),
                      op(qd1-rd@qd@rd.conj().T))
            assert err < 2e-12
            assert commutant_dimension(yu1@yu1.conj().T, yd1@yd1.conj().T) == dim
            max_basis_covariance = max(max_basis_covariance, err)
        c = comm(hu, hd)
        prod = hd[0,1]*hd[1,2]*hd[2,0]
        vand = (1-4)*(4-9)*(9-1)
        determinant_formula = 2j*vand*np.imag(prod)
        assert abs(np.linalg.det(c)-determinant_formula) < 1e-11
        cases.append(dict(case=label,components=components(hd),
            hermitian_charge_dimension=dim,
            down_eigenvalues=np.linalg.eigvalsh(hd).tolist(),
            commutator_det=[float(np.linalg.det(c).real),float(np.linalg.det(c).imag)]))

    # Exact graph ranks for all 8 three-vertex edge subsets. The proof in the
    # note works for all n; exhaustive finite graph calibration is not that proof.
    graph_checks = []
    for mask in range(8):
        edges = [(i,j) for k,(i,j) in enumerate([(0,1),(1,2),(0,2)]) if mask & (1<<k)]
        rows = []
        for i,j in edges:
            row = [0]*3
            row[i], row[j] = 1, -1
            rows.append(row)
        a = np.diag([4., 5., 6.])
        for i,j in edges:
            a[i,j] = a[j,i] = 1
        rank = fraction_rank(rows)
        comp = len(components(a))
        assert rank == 3-comp
        assert commutant_dimension(hu, a) == comp
        graph_checks.append(dict(mask=mask,components=comp,exact_incidence_rank=rank))

    # Rational Sylvester certificate and real 3x3 skew determinant: no floating
    # threshold is used to establish the CP-conserving counterexample.
    b = [[F(x) for x in row] for row in [[4,1,1],[1,5,1],[1,1,6]]]
    minors = [b[0][0], b[0][0]*b[1][1]-b[0][1]*b[1][0], det3(b)]
    real_comm = [[F([1,4,9][i]-[1,4,9][j])*b[i][j] for j in range(3)] for i in range(3)]
    exact_real_det = det3(real_comm)
    assert all(x > 0 for x in minors) and exact_real_det == 0

    # Non-Hermitian right charge is the obstruction, even though solving the
    # intertwiner alone is possible for every full-rank Y.
    q_bad = np.diag([0., 1., 2.])
    qu_bad, qd_bad = charge_lift(q_bad, yu, psqrt(hd_real), h)
    bad_hermitian = op(qd_bad-qd_bad.conj().T)
    assert op(comm(q_bad, hd_real)) > 1 and bad_hermitian > .1
    assert op(q_bad@psqrt(hd_real)-psqrt(hd_real)@(qd_bad+h*np.eye(3))) < 1e-12

    # A degenerate up spectrum does not fix an eigenbasis. A connected graph in
    # an arbitrary such basis need not yield a scalar commutant.
    degenerate_dim = commutant_dimension(np.eye(3), hd_real)
    assert degenerate_dim == 3
    assert op(comm(hd_real, np.eye(3))) == 0

    # Full-rank omission leaves a right-handed null-space generator free, even
    # with a scalar left charge forced by the two left mass matrices.
    yu_zero = np.diag([0., 2., 3.])
    ql_zero = np.eye(3)/6
    qu_zero = (1/6+h)*np.eye(3)+np.diag([.75,0.,0.])
    qd_zero = (1/6-h)*np.eye(3)
    yd_zero = psqrt(hd_real)
    rank_zero_error = max(op(ql_zero@yu_zero-yu_zero@(qu_zero-h*np.eye(3))),
                          op(ql_zero@yd_zero-yd_zero@(qd_zero+h*np.eye(3))))
    assert rank_zero_error < 1e-12
    assert commutant_dimension(yu_zero@yu_zero.T, hd_real) == 1

    # Fixed-basis finite calibration bound. This is not a permission to violate
    # an exact gauge symmetry or anomaly cancellation in a physical theory.
    approximate_rows = []
    for step in [.002, .02, .2]:
        ql = np.diag([0., step, 2*step])
        eps, bmin = op(comm(ql, hd_chain)), 1.
        spread, bound = 2*step, 2*eps/bmin
        assert spread <= bound + 1e-14
        approximate_rows.append(dict(step=step,commutator_norm=eps,
                                     graph_edge_lower_bound=bmin,charge_spread=spread,bound=bound))
    weak_rows = []
    for delta in [1e-2, 1e-4, 1e-6]:
        a = np.array([[4.,1.,0.],[1.,5.,delta],[0.,delta,6.]])
        ql = np.diag([0.,0.,1.])
        eps = op(comm(ql, a))
        assert math.isclose(eps, delta, rel_tol=1e-12)
        weak_rows.append(dict(edge=delta,commutator_norm=eps,charge_spread=1.))

    # Physical-spectrum/mixing bridge, using synthetic fixed-scale data rather
    # than current observational numbers. Nonzero invariant is sufficient only.
    mixing_rows = []
    for delta in [0., .7]:
        v = ckm(.23,.07,.31,delta)
        a = v@np.diag([2.,5.,11.])@v.conj().T
        cdet = np.linalg.det(comm(hu, a))
        j = np.imag(v[0,0]*v[1,1]*v[0,1].conjugate()*v[1,0].conjugate())
        vd = (2-5)*(5-11)*(11-2)
        assert abs(cdet-2j*j*120*vd) < 2e-11
        assert commutant_dimension(hu, a) == 1
        mixing_rows.append(dict(cp_phase=delta,J=float(j),
            commutator_det=[float(cdet.real),float(cdet.imag)],
            hermitian_charge_dimension=1))

    sources = [BASE/"archive_923_934"/f"research_note_{n}.md" for n in [924,926,927]]
    sources += [BASE/"archive_531_553"/f"research_note_{n}.md" for n in [531,532,536,537]]
    sources += [BASE/"archive_956_989/981/drafts/common_parent_contract_v1.md"]
    return dict(round=1011,calibration_groups=1,cumulative_test_groups=3789,
        scope="finite flavor-matrix calibration; analytic proof covers the stated class",
        graph_cases=graph_checks,physical_charge_cases=cases,
        largest_intertwiner_residual=max_intertwiner,
        largest_hermiticity_residual=max_hermitian,
        largest_basis_covariance_residual=max_basis_covariance,
        real_cp_conserving_certificate=dict(leading_principal_minors=[str(x) for x in minors],
            exact_commutator_determinant=str(exact_real_det),commutant_dimension=1),
        invalid_noncommuting_left_charge_right_hermiticity_defect=bad_hermitian,
        degenerate_up_spectrum_counterexample_dimension=degenerate_dim,
        rank_deficient_right_charge=dict(null_space_shift=.75,intertwiner_residual=rank_zero_error,
            full_anomaly_test_claimed=False),
        finite_precision=approximate_rows,weak_edge_counterexample=weak_rows,
        synthetic_mixing=mixing_rows,
        all_scientific_calibrations_passed=True,
        new_cognitive_axioms=0,actual_experiment=False,all_physics_generated=False,
        historical_source_sha256={str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    if args.write:
        with OUT.open("x", encoding="utf8") as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
    else:
        compare(result, json.loads(OUT.read_text(encoding="utf8")))
    print(json.dumps({k:v for k,v in result.items() if k != "historical_source_sha256"},
                     ensure_ascii=False, indent=2))
