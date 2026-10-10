"""Count0 exact checks: old G6 characteristic lattice and a compact axion period.

Default is read-only. --save-exclusive creates results.json once.
This is algebra calibration, not a construction of bundles or a quantum detector.
"""
from fractions import Fraction as F
import json
from pathlib import Path
import sys


def instantons(u, v, w):
    return (F(2, 3)*u-v, F(1, 2)*u-w, F(1, 36)*u)


def integer_character(k3, k2, k1, kt=0):
    return ((24*k3+18*k2+k1)/F(36), -F(k3), -F(k2), F(kt))


def periodic(k3, k2, k1, kt=0):
    return all(v.denominator == 1 for v in integer_character(k3, k2, k1, kt))


def pair_integral(x, y):
    # H^2(S2 x S2): basis x,y with x^2=y^2=0 and integral xy=1.
    return F(x[0])*y[1]+F(x[1])*y[0]


def calculate():
    a = (1, 1)
    roots3 = (a, (0, 0), (0, 0))
    roots2 = ((-1, -1), (0, 0))
    a_square = pair_integral(a, a)
    trace3 = sum(pair_integral(x, x) for x in roots3)/2
    trace2 = sum(pair_integral(x, x) for x in roots2)/2
    # E3=L+1+1 and E2=L^-1+1; detE3 detE2 is trivial.
    nu_from_roots = (trace3-a_square/6, trace2-a_square/4, a_square/72)
    assert nu_from_roots == (F(2, 3), F(1, 2), F(1, 36))
    assert nu_from_roots == instantons(F(1), F(0), F(0))

    cases = (
        ("color_unit_only_rejected", (F(1), F(0), F(0)), False),
        ("color_unit_hypercharge_completion", (F(1), F(0), F(12)), True),
        ("color_three_only_allowed", (F(3), F(0), F(0)), True),
        ("fractional_color_rejected_by_v_basis", (F(1, 2), F(0), F(24)), False),
        ("fractional_weak_rejected_by_w_basis", (F(0), F(1, 2), F(27)), False),
    )
    rows = []
    for label, k, expected in cases:
        c = integer_character(*k)
        assert periodic(*k) == expected
        kgamma = k[1]/2+k[2]/36
        assert c[0] == 2*k[0]/3+kgamma
        rows.append({"case": label, "K3_K2_K1": list(map(str, k)),
                     "integer_basis_coefficients_u_v_w_t": list(map(str, c)),
                     "passes_period": expected, "bare_k_gamma": str(kgamma),
                     "u_basis_phase_turns_mod_1": str(c[0] % 1)})
    assert integer_character(1, 0, 0)[0] % 1 == F(2, 3)
    # t=0 on the witness: no gravitational coefficient can alter this failure.
    assert [sum(c*b for c, b in zip(integer_character(1, 0, 0, kt), (1,0,0,0)))
            for kt in (F(0), F(2), F(1, 7))] == [F(2, 3)]*3

    # Independent lattice coordinates: any integer (z,K3,K2) determines
    # K1=36z-24K3-18K2. Basis vectors suffice for this linear identity.
    lattice = []
    for z, k3, k2 in ((1,0,0), (0,1,0), (0,0,1), (2,1,-1)):
        k1 = 36*z-24*k3-18*k2
        assert integer_character(k3,k2,k1) == (F(z),F(-k3),F(-k2),F(0))
        assert k1 % 6 == 0
        lattice.append([k3,k2,k1])

    # Old629 chiral-index rows, including gravity. The overall sign depends
    # on the declared passive field-redefinition convention.
    A = (
        (3,-2,-3,-12), (2,-1,0,-6), (1,-1,0,-6),
        (1,0,-1,-4), (1,0,0,-2), (0,0,0,-2)
    )
    rephasings = []
    for au,av,aw,at in A:
        delta = (-av,-aw,36*au+24*av+18*aw,at)
        assert integer_character(*delta) == tuple(map(F,(au,av,aw,at)))
        assert periodic(*delta)
        rephasings.append(list(delta))

    # Same local fa, different fundamental period, no statement about a UV model.
    fa = F(7)
    local_examples = []
    for k3,k2,k1 in ((1,0,12),(3,0,0)):
        f = k3*fa
        assert F(k3)/f == 1/fa
        local_examples.append({"K3_K2_K1":[k3,k2,k1],"fundamental_f":str(f),
                               "local_fa":str(f/F(k3)),
                               "bare_k_gamma":str(F(k2,2)+F(k1,36))})

    return {
        "passed": True,
        "new_scientific_groups": 0,
        "new_empirical_groups": 0,
        "new_cognitive_axioms": 0,
        "scope": "closed ordinary-spin G6 constant-period consistency; exact Fraction algebra only",
        "root_bundle_witness_nu3_nu2_nu1": list(map(str,nu_from_roots)),
        "cases": rows,
        "integer_lattice_generators_and_control": lattice,
        "old629_index_rephasing_delta_K3_K2_K1_Kt": rephasings,
        "same_local_fa_distinct_period_examples": local_examples,
        "not_certified": [
            "spacetime axion winding, defects, boundaries, or full QFT construction",
            "physical topology preparation or axion detection",
            "global gauge group or fundamental axion period selected by cognition",
            "all non-spin tangential structures or extra topological sectors"
        ]
    }


def main():
    result = calculate()
    target = Path(__file__).with_name("results.json")
    if sys.argv[1:] == ["--save-exclusive"]:
        with target.open("x", encoding="utf-8", newline="\n") as out:
            json.dump(result, out, ensure_ascii=False, indent=2)
            out.write("\n")
    elif sys.argv[1:]:
        raise SystemExit("only --save-exclusive is supported; default is read-only")
    else:
        assert json.loads(target.read_text(encoding="utf-8")) == result
    print("axion global-period exact checks passed (count0)")


if __name__ == "__main__":
    main()
