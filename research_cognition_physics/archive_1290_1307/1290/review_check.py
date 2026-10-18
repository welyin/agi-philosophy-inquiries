"""Independent source-moment audit; no imports of the author's calculation.

Default: recompute and compare with the frozen local JSON. --write creates it
exclusively. All floating quadratures below are diagnostics, not error bounds.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "review_results.json"
ETA = (F(-1), F(1), F(1), F(1))


def shell(m, u):
    u = tuple(map(F, u))
    s = sum(x*x for x in u)
    assert s < 1
    return (-m*(1+s)/(1-s), *(2*m*x/(1-s) for x in u))


def dot(p, q):
    return sum(ETA[i]*p[i]*q[i] for i in range(4))


def vertex(p, q, m, keep_mass=True):
    # Independent differentiation of the two creation plane waves.
    return [[-(p[i]*q[j]+p[j]*q[i])/2
             + (ETA[i]*(dot(p,q)-(m*m if keep_mass else 0))/2 if i == j else 0)
             for j in range(4)] for i in range(4)]


def ward(t, p, q):
    return [sum(ETA[i]*(p[i]+q[i])*t[i][j] for i in range(4))
            for j in range(4)]


def sf(x):
    x = F(x)
    return f"{x.numerator}/{x.denominator}"


def angular_integral(z):
    """Integral_-1^1 (1/20+c^2/60) exp(-z*c) dc, analytic in c."""
    ans = np.empty_like(z)
    small = z < 0.2
    zs = z[small]
    val = np.zeros_like(zs)
    for k in range(9):
        val += 2*zs**(2*k)/math.factorial(2*k) * (
            1/(20*(2*k+1)) + 1/(60*(2*k+3)))
    ans[small] = val
    zz = z[~small]
    i0 = 2*np.sinh(zz)/zz
    i2 = 2*((zz*zz+2)*np.sinh(zz)-2*zz*np.cosh(zz))/(zz**3)
    ans[~small] = i0/20+i2/60
    return ans


def gaussian_diagnostic(m, order):
    # Integrate c analytically, unlike a three-dimensional sampled integral.
    nodes, weights = np.polynomial.legendre.leggauss(order)
    p = 4*(nodes+1)  # finite diagnostic interval [0,8]
    w = 4*weights
    pp, qq = p[:,None], p[None,:]
    ep, eq = np.sqrt(pp*pp+m*m), np.sqrt(qq*qq+m*m)
    z = 2*pp*qq
    integrand = pp**4*qq**4/(ep*eq) * np.exp(-(ep+eq)**2-pp*pp-qq*qq)
    integrand *= angular_integral(z)
    return float(np.sum(w[:,None]*w[None,:]*integrand)/(16*math.pi**4))


def run():
    pars = [(0,0,0), (F(1,4),0,0), (F(1,5),F(1,6),F(1,7)),
            (F(-1,3),F(1,8),F(-1,9))]
    cases = 0
    max_ratio = F(0)
    bad_mass_defect = F(0)
    for m in (F(1), F(3,2), F(7,3)):
        for u in pars:
            p = shell(m,u)
            assert dot(p,p) == -m*m
            for v in pars:
                q = shell(m,v)
                t = vertex(p,q,m)
                assert ward(t,p,q) == [F(0)]*4
                for i in range(4):
                    for j in range(4):
                        assert t[i][j] == t[j][i]
                        ratio = abs(t[i][j])/(p[0]*q[0])
                        max_ratio = max(max_ratio, ratio)
                        assert ratio <= 3
                broken = ward(vertex(p,q,m,False),p,q)
                bad_mass_defect = max(bad_mass_defect,max(map(abs,broken)))
                cases += 1
            q_opposite = (p[0], *(-x for x in p[1:]))
            assert vertex(p,q_opposite,m)[0][0] == 0
    assert bad_mass_defect > 0

    # Wick normalization from the two possible Kronecker contractions.
    # Symmetric real pair coefficient; sum over BOTH ordered mode labels.
    c = [[F(2,3),F(-3,7)],[F(-3,7),F(5,11)]]
    norm_by_contractions = sum(c[i][j]*c[k][l] * (
        int(i==k and j==l)+int(i==l and j==k))
        for i in range(2) for j in range(2) for k in range(2) for l in range(2))
    norm_by_fock_states = 2*c[0][0]**2 + 4*c[0][1]**2 + 2*c[1][1]**2
    assert norm_by_contractions == norm_by_fock_states
    assert norm_by_contractions == 2*sum(v*v for row in c for v in row)

    # Reconstruct the constant independently before inserting pi.
    radial_prefactor_without_pi4 = F(2,4) * F(16,64) * 9
    assert radial_prefactor_without_pi4 == F(9,8)
    beta44 = sum(F((-1)**k * math.comb(3,k),4+k) for k in range(4))
    assert beta44 == F(1,140)
    constants = {}
    for r in (5,6,7,8):
        c_r_no_pi = radial_prefactor_without_pi4*beta44/F(2*r-8)
        assert c_r_no_pi == F(9,1120*(2*r-8))
        constants[str(r)] = {"C_r_times_pi4":sf(c_r_no_pi),
                             "doubling_mass_bound_ratio":sf(F(1,2**(2*r-8)))}

    # Compact six-box convolution followed by a fixed normalized smooth bump.
    box_h = F(1,8)
    derivative_l1_upper = (2/box_h)**6
    support_radius = 6*box_h/2+F(1,8)
    assert support_radius == F(1,2)
    mass = 2**20
    compact_upper = F(9,1120*(12-8)*3**4) * derivative_l1_upper**2/F((2*mass)**4)
    assert compact_upper == F(1,40320*2**36)
    assert compact_upper < F(1,10**15)

    # Sharp-time shear: q=-p, p=(R,R,0), complete t12=R^2.
    # Its pair density grows like R^2, and radial measure supplies R^2 dR.
    sharp = []
    for radius in (1.,10.,100.):
        ep2 = 2*radius**2+1
        pair_over_radius2 = radius**4/ep2/radius**2
        sharp.append({"radius":radius,"pair_weight_divided_by_R2":pair_over_radius2})
    assert sharp[-1]["pair_weight_divided_by_R2"] > .4999

    diagnostics = []
    for m in (.5,1.,2.):
        lo, hi = gaussian_diagnostic(m,48), gaussian_diagnostic(m,72)
        assert lo > 0 and hi > 0
        assert abs(lo-hi) <= 2e-9*hi+1e-17
        # Hermite Cauchy--Schwarz gives B6^2 <= 6!, independently of M.
        rough = 9*math.factorial(6)/(1120*math.pi**4*4)*(2*m)**(-4)
        assert hi < rough
        diagnostics.append({"M":m,"order48":lo,"order72":hi,
                            "absolute_grid_difference":abs(lo-hi),
                            "rigorous_rough_bound":rough})

    return {
        "status":"PASS_WITH_SCOPE",
        "exact_vertex_checks":{"on_shell_pairs":cases,
            "ward_identity":True,"opposite_spatial_t00_zero":True,
            "tested_max_vertex_ratio":sf(max_ratio),
            "mass_term_deletion_ward_defect":sf(bad_mass_defect)},
        "wick":{"two_contraction_norm":sf(norm_by_contractions),
                "factor_two_confirmed":True},
        "constant_derivation":{"beta44":sf(beta44),"orders":constants},
        "compact_witness":{"support_radius":sf(support_radius),
            "derivative_L1_upper":sf(derivative_l1_upper),"M":mass,
            "variance_strict_upper_using_pi_gt_3":sf(compact_upper),
            "below_1e_minus15":True},
        "sharp_time_diagnostic":sharp,
        "gaussian_analytic_angle_diagnostics":diagnostics,
        "scope":{"fixed_finite_spacetime_test_menu":True,
            "full_minimal_mass_vertex":True,"mass_varying_vacuum_subtraction_declared":True,
            "quadrature_is_rigorous_certificate":False,
            "actual_stress_instrument_certified":False,
            "curved_background_or_HG_extension_certified":False,
            "nonzero_portal_full_history_certified":False,
            "all_axiom_species_freedom_proved":False}
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--write",action="store_true")
    args=ap.parse_args()
    payload=run()
    text=json.dumps(payload,ensure_ascii=False,indent=2,sort_keys=True)+"\n"
    if args.write:
        with OUT.open("x",encoding="utf-8",newline="\n") as f:
            f.write(text)
        print("CREATED review_results.json; independent checks PASS_WITH_SCOPE")
    else:
        expected=json.loads(OUT.read_text(encoding="utf-8"))
        if expected != payload:
            raise AssertionError("independent saved payload differs")
        print("MATCH; independent checks PASS_WITH_SCOPE")


if __name__ == "__main__":
    main()
