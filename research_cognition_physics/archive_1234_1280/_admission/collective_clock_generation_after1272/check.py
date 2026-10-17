"""Exact finite collective admission; default replay never rewrites saved output.

Polynomial coefficient identities prove the displayed commuting solution. Matrix
substitution at rational times is a separate diagnostic, not a continuum proof.
No finite-N quantum instrument or geometric metric is certified here.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import json
import numpy as np

HERE = Path(__file__).resolve().parent
ALPHA, U0 = F(8, 195), F(65, 8)
R, B = F(400, 10003), F(9997, 10003)
END = F(1, 100)


def trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def add(*ps):
    out = [F(0)] * max(map(len, ps))
    for p in ps:
        for i, x in enumerate(p):
            out[i] += x
    return trim(out)


def scale(p, x):
    return trim([x*y for y in p])


def mul(p, q):
    out = [F(0)] * (len(p)+len(q)-1)
    for i, x in enumerate(p):
        for j, y in enumerate(q):
            out[i+j] += x*y
    return trim(out)


def deriv(p):
    return trim([i*p[i] for i in range(1, len(p))] or [F(0)])


def val(p, t):
    ans = F(0)
    for x in reversed(p):
        ans = ans*t + x
    return ans


def zero(p):
    assert all(x == 0 for x in p), p


def encode(x):
    return {"fraction": str(x), "decimal": float(x)}


def branch(sign):
    s0 = -2*B-sign*R
    d = s0*s0+ALPHA*U0*U0
    den, ns = [F(1), -4*s0, 4*d], [s0, -2*d]
    p = [[F(1), v, (v*v+1)/4] for v in (8*B+sign*12*R, 16*B)]
    return den, ns, p


def eigen_to_site(plus, minus):
    return np.array([[(plus+minus)/2, (plus-minus)/2],
                     [(plus-minus)/2, (plus+minus)/2]], dtype=object)


def run():
    assert B*B+F(3, 4)*R*R == 1
    assert ALPHA*U0*3 == 1
    poly_checks = 0
    for sign in (1, -1):
        den, ns, ps = branch(sign)
        q = add([F(1)], *ps)
        qp, qpp, dp = deriv(q), deriv(deriv(q)), deriv(den)
        identities = [
            add(q, scale(den, -3)),
            add(scale(dp, -1), scale(ns, -4)),
            add(mul(deriv(ns), den), scale(mul(ns, dp), -1),
                [2*ALPHA*U0*U0], scale(mul(ns, ns), -2)),
            add(scale(q, ALPHA*U0), scale(den, -1)),
            add(mul(ns, q), scale(mul(den, qp), F(1, 4))),
            add(mul(mul(ns, ns), q), scale(mul(mul(ns, den), qp), F(1, 2)),
                scale(mul(mul(den, den), qpp), F(1, 8)), scale(den, -U0)),
        ]
        for p in ps:
            identities.append(add(scale(mul(p, deriv(deriv(p))), 2),
                                  scale(mul(deriv(p), deriv(p)), -1), [F(-1)]))
        for identity in identities:
            zero(identity)
            poly_checks += 1

    I = np.array([[F(1), F(0)], [F(0), F(1)]], dtype=object)
    X = np.array([[F(0), F(1)], [F(1), F(0)]], dtype=object)
    s = -2*B*I-R*X
    t1, t2 = R*X, B*I-R*X/2
    u = s@s+4*t1@t1+4*t2@t2+I/8
    assert np.all(u == U0*I)
    assert np.all(-u+ALPHA*u@(3*I)@u == 0)
    assert np.all(s+2*(t1+t2) == 0)
    matrix_checks = []
    for diagonal in ((1, 1), (1, 2), (2, 3)):
        theta = np.diag([F(x) for x in diagonal])
        for kappa in (F(0), F(1), F(5, 4), F(2)):
            y = kappa*(-2*theta@s)+(1-kappa)*(-4*s@theta+2*theta@s)
            rr = 4*t1@theta+y
            assert np.all(rr-rr.T == 2*(5-4*kappa)*R*(X@theta-theta@X))
            expected = 8*B*theta+6*R*(X@theta+theta@X)
            assert np.all(rr+rr.T == expected)
            matrix_checks.append({"lapse": list(diagonal), "kappa": str(kappa),
                                  "rotation_01": str((rr-rr.T)[0, 1]),
                                  "p1_tangent_01": str(expected[0, 1])})

    # Independent direct substitution into all scalar eigenbranch original ODEs.
    # The polynomial checks above, not this grid, certify the entire interval.
    substitutions = 0
    endpoint_branches = []
    for sign in (1, -1):
        den, ns, ps = branch(sign)
        for j in range(11):
            tau = END*F(j, 10)
            dv, nv = val(den, tau), val(ns, tau)
            assert dv > 0
            uv, sv = U0/dv, nv/dv
            up = -U0*val(deriv(den), tau)/dv**2
            sp = (val(deriv(ns), tau)*dv-nv*val(deriv(den), tau))/dv**2
            pp = [val(p, tau) for p in ps]
            pv = [val(deriv(p), tau) for p in ps]
            pa = [val(deriv(deriv(p)), tau) for p in ps]
            ts = [sv/2+v/(8*p) for p, v in zip(pp, pv)]
            tv = [sp/2+(acc*p-v*v)/(8*p*p) for p, v, acc in zip(pp, pv, pa)]
            assert all(p >= 1 for p in pp)
            reconstructed_u = sv*sv+sum(4*p*t*t+1/(16*p) for p, t in zip(pp, ts))
            qv = 1+sum(pp)
            assert reconstructed_u == uv
            assert -uv+ALPHA*uv*qv*uv == 0
            assert sv+2*sum(p*t for p, t in zip(pp, ts)) == 0
            assert up == 4*sv*uv
            assert sp == -2*ALPHA*uv*uv+2*sv*sv
            for p, v, t, tp in zip(pp, pv, ts, tv):
                assert v == 8*p*t-4*p*sv
                assert tp == -4*t*t-ALPHA*uv*uv+1/(16*p*p)+4*sv*t
            substitutions += 1
        endpoint_branches.append({"sign": sign, "p1": encode(pp[0]),
                                  "p2": encode(pp[1]), "U": encode(uv), "s": encode(sv)})

    c = [F(1), 8*B, F(65, 4)+24*R*R]
    off = [F(0), 12*R, 48*B*R]
    plus, _, pplus = branch(1)
    minus, _, pminus = branch(-1)
    zero(add(c, scale(add(pplus[0], pminus[0]), F(-1, 2))))
    zero(add(off, scale(add(pplus[0], scale(pminus[0], -1)), F(-1, 2))))
    assert 8*B-12*R > 7 and 8*B > 7 and 16*B > 15
    assert B < 1 and R < F(1, 25)
    # Positive coefficients and lower clock slope certify these on all [0, END].
    off_upper = F(12, 25)*END+F(48, 25)*END*END
    assert off_upper == F(78, 15625) < F(1, 200)
    assert val(off, END) > 0
    assert (8*B)**2 > 49
    return {
        "status": "PASS_EXACT_COLLECTIVE_CLOCK_SCOPE",
        "parameters": {"alpha": str(ALPHA), "b": str(B), "r": str(R), "U0": str(U0)},
        "exact_polynomial_identities": poly_checks+2,
        "rational_branch_ode_substitutions": substitutions,
        "frame_tests": matrix_checks,
        "finite_interval": ["0", str(END)],
        "p1_site_clock_derivative_lower_bound": encode(8*B),
        "initial_diagonal_lapse_clock_jacobian": "8*b times I2",
        "initial_clock_jacobian_determinant": encode((8*B)**2),
        "p1_min_eigen_initial_slope": encode(8*B-12*R),
        "p2_uniform_clock_initial_slope": encode(16*B),
        "p1_site_clock_endpoint": encode(val(c, END)),
        "p1_offdiagonal_endpoint": encode(val(off, END)),
        "normalized_p1_offdiagonal_endpoint": encode(val(off, END)/val(c, END)),
        "uniform_offdiagonal_upper_bound": encode(off_upper),
        "endpoint_branches": endpoint_branches,
        "scope": {
            "classical_collective_equations_only": True,
            "finite_N_physical_state_or_instrument": False,
            "three_dimensional_metric_or_GR_verified": False,
            "kappa_is_newton_coupling": False,
            "full_CT_or_lapse_independence_verified": False,
            "new_axiom": False,
            "new_numbered_round": False,
        },
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="Create result once; never overwrite.")
    args = parser.parse_args()
    out = run()
    dest = HERE/"results.json"
    if args.write:
        with dest.open("x", encoding="utf-8") as f:
            f.write(json.dumps(out, ensure_ascii=False, indent=2)+"\n")
        print("SAVED_EXACT_COLLECTIVE_CLOCK_SCOPE")
    else:
        assert out == json.loads(dest.read_text(encoding="utf-8"))
        print("READ_ONLY_MATCH: exact constraints, polynomial ODE certificates and clock window")
