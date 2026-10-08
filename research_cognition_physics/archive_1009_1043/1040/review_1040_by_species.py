"""Independent reviewer: rational Maxwell polynomial and covariant metric variation.

Does not import the author's science or verification implementation.
"""
from fractions import Fraction as Q
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "review_1040_by_species_results.json"


def mul(a, b):
    return [[sum((x*y for x, y in zip(row, col)), Q(0))
             for col in zip(*b)] for row in a]


def rank(a):
    a = [list(row) for row in a]
    piv = 0
    for col in range(len(a[0])):
        p = next((i for i in range(piv, len(a)) if a[i][col]), None)
        if p is None:
            continue
        a[piv], a[p] = a[p], a[piv]
        q = a[piv][col]
        a[piv] = [x/q for x in a[piv]]
        for i in range(len(a)):
            if i != piv:
                q = a[i][col]
                a[i] = [x-q*y for x, y in zip(a[i], a[piv])]
        piv += 1
    return piv


def rational_maxwell():
    ans = []
    for R in [Q(5, 4), Q(3, 2), Q(2), Q(5, 2)]:
        inv_eps = [R, R, 1/R]
        inv_mu = [1/R, 1/R, 1/R**3]
        for k in [(Q(1), Q(0), Q(0)), (Q(0), Q(0), Q(1)),
                  (Q(1, 3), Q(2, 3), Q(2, 3)),
                  (Q(-2, 7), Q(3, 7), Q(6, 7))]:
            x, y, z = k
            K = [[Q(0), -z, y], [z, Q(0), -x], [-y, x, Q(0)]]
            N = [[Q(0) for _ in range(6)] for _ in range(6)]
            for i in range(3):
                for j in range(3):
                    N[i][j+3] = -K[i][j]*inv_mu[j]
                    N[i+3][j] = K[i][j]*inv_eps[j]
            w2 = z*z+(x*x+y*y)/R**2
            n3 = mul(mul(N, N), N)
            assert n3 == [[w2*v for v in row] for row in N]
            assert rank(N) == 4
            gauss = [list(k)+[Q(0)]*3, [Q(0)]*3+list(k)]
            assert rank(gauss) == 2
            assert mul(gauss, N) == [[Q(0)]*6]*2
            assert sum(N[i][i] for i in range(6)) == 0
            # Three distinct roots, rank four and zero trace => each nonzero
            # root has multiplicity two. Image equals the Gauss kernel.
            ans.append(dict(R=str(R), k=[str(v) for v in k],
                            omega_squared=str(w2), matrix_rank=4,
                            minimal_polynomial_identity=True,
                            physical_polarization_multiplicity=2))
    return ans


def metric_density(g, tensor, T):
    det = np.linalg.det(g)
    assert det < 0
    volume = np.sqrt(-det)
    inverse = np.linalg.inv(g)
    raised = inverse @ tensor @ inverse.T
    S = -np.sum(tensor*raised)/4
    pfaffian = (tensor[0, 1]*tensor[2, 3]
                -tensor[0, 2]*tensor[1, 3]
                +tensor[0, 3]*tensor[1, 2])
    P = -pfaffian/volume
    return volume*T*(1-np.sqrt(1-2*S/T-P*P/T**2))


def covariant_metric_variation():
    rng = np.random.default_rng(104039)
    eta = np.diag([-1., 1., 1., 1.])
    maxerr = 0.
    for _ in range(18):
        E, B = rng.uniform(-.25, .25, size=(2, 3))
        T = rng.uniform(.7, 2.5)
        ex, ey, ez = E
        bx, by, bz = B
        tensor = np.array([[0, -ex, -ey, -ez], [ex, 0, bz, -by],
                           [ey, -bz, 0, bx], [ez, by, -bx, 0]])
        S, P = (E@E-B@B)/2, E@B
        root = np.sqrt(1-2*S/T-P*P/T**2)
        L, ls, lp = T*(1-root), 1/root, P/(T*root)
        lower = ls*(tensor @ eta @ tensor.T)+eta*(L-P*lp)
        upper = eta @ lower @ eta
        h = .001
        for i in range(4):
            for j in range(i, 4):
                direction = np.zeros((4, 4))
                direction[i, j] = direction[j, i] = 1.
                f = lambda t: metric_density(eta+t*direction, tensor, T)
                derivative = (f(-2*h)-8*f(-h)+8*f(h)-f(2*h))/(12*h)
                expected = .5*np.sum(upper*direction)
                maxerr = max(maxerr, abs(float(derivative-expected)))
    assert maxerr < 1e-8
    return dict(backgrounds=18, covariant_metric_directions=180,
                method="five-point real derivative in g_covariant",
                maximum_stress_variation_residual=maxerr)


def threshold_controls():
    cases = []
    for B2, delta in [(Q(1, 3), Q(1, 12)), (Q(7, 2), Q(2, 7))]:
        boundary = B2*(1-delta)/delta
        assert B2/(boundary+B2) == delta
        assert B2/(Q(99, 100)*boundary+B2) > delta
        assert B2/(Q(101, 100)*boundary+B2) < delta
        cases.append(dict(B_squared=str(B2), delta=str(delta),
                          T_boundary=str(boundary), both_sides_checked=True))
    a, b, B2 = Q(1, 100), Q(7, 400), Q(1, 4)
    ls = 1-a*B2
    eps = [ls, ls, ls+2*b*B2]
    mu = [ls, ls, ls-2*a*B2]
    assert all(x > 0 for x in eps+mu)
    vy2, vz2 = mu[2]/eps[1], mu[1]/eps[2]
    assert Q(0) < vz2 < vy2 < 1
    return dict(threshold_cases=cases, positive_birefringence=True,
                two_quartic_speeds_squared=[str(vy2), str(vz2)])


def run():
    receipt = json.loads((HERE/"research_round_1040_checks.json").read_text("utf8"))
    for name, digest in receipt["owned_sha256"].items():
        assert hashlib.sha256((HERE/name).read_bytes()).hexdigest() == digest
    return dict(round=1040, independent_reviewer="round1036_species_selection",
                passed=True, scientific_groups_added=0,
                rational_six_component_cases=rational_maxwell(),
                independent_metric_variation=covariant_metric_variation(),
                finite_controls=threshold_controls(),
                author_owned_sha256=receipt["owned_sha256"],
                reviewed_receipt_sha256=hashlib.sha256(
                    (HERE/"research_round_1040_checks.json").read_bytes()).hexdigest())


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    data = run()
    if args.write:
        with OUT.open("x", encoding="utf8") as f:
            json.dump(data, f, ensure_ascii=False, indent=2)
            f.write("\n")
    else:
        assert data == json.loads(OUT.read_text("utf8"))
    print(json.dumps(dict(passed=True, rational_six_component_cases=16,
                          **data["independent_metric_variation"]),
                     ensure_ascii=False))
