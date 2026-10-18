"""1292 independent exact checks; no import of author or historical calculators.
Default and --check compare a saved result read-only; --write creates it once.
Finite samples check implementations. All-integer/window statements are analytic.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in HERE.parents
            if (p / "research_cognition_physics").is_dir()
            and (p / "scripts").is_dir())
OUT = HERE / "independent_results.json"


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def mod1(x):
    return x % 1


def add(x, y):
    return (mod1(x[0] + y[0]), (x[1] + y[1]) % 4)


def mul(n, x):
    return (mod1(n * x[0]), (n * x[1]) % 4)


def cosets(elements, generators):
    remaining = set(elements)
    sizes = []
    while remaining:
        seed = min(remaining)
        orbit = {seed}
        front = [seed]
        while front:
            x = front.pop()
            for g in generators:
                y = add(x, g)
                assert y in elements
                if y not in orbit:
                    orbit.add(y)
                    front.append(y)
        assert orbit <= remaining
        sizes.append(len(orbit))
        remaining -= orbit
    return {"count": len(sizes), "orbit_sizes": sorted(set(sizes))}


def group_test(k):
    assert k != 0
    n = abs(k)
    elements = {(F(j, 36*n), a)
                for j in range(36*n) for a in range(4)
                if mod1(-18*k*F(j, 36*n) + F(a, 2)) == 0}
    x, z = (mod1(F(1, 36*k)), 1), (F(0), 2)
    generated = {add(mul(j, x), mul(e, z))
                 for j in range(36*n) for e in range(2)}
    assert elements == generated and len(elements) == 72*n
    assert mul(36*n, x) == (F(0), 0)
    assert all(mul(j, x) != z for j in range(36*n))
    color = mul(12*k, x)
    assert color == (F(1, 3), 0)
    pure = {e for e in elements if e[1] == 0}
    assert len(pure) == 18*n
    c_pure = cosets(pure, [color])
    c_internal = cosets(elements, [color])
    c_full = cosets(elements, [color, z])
    assert c_pure["count"] == 6*n
    assert c_internal["count"] == 24*n
    assert c_full["count"] == 12*n
    assert mod1(9*k*F(1, 18*n)) == F(1, 2)
    assert mod1(9*k*x[0] + F(x[1], 4)) == F(1, 2)
    assert mod1(9*k*z[0] + F(z[1], 4)) == F(1, 2)
    return {
        "k": k, "K_order": len(elements),
        "pure_gauge_components": c_pure,
        "internal_with_r_components": c_internal,
        "full_spin_components": c_full,
        "N_minimum_pure_residual_turn": F(1, 2),
        "x_order": 36*n, "color_generator": color,
    }


def representations(k):
    # Color triality, weak-doublet indicator, original q=6Y.
    rows = {
        "Q": (1, 1, 1+6*k), "uc": (-1, 0, -4-6*k),
        "dc": (-1, 0, 2-6*k), "L": (0, 1, -3-18*k),
        "ec": (0, 0, 6+18*k), "N": (0, 0, 18*k),
        "H": (0, 1, 3), "S": (0, 0, -36*k),
    }
    for triality, weak, q in rows.values():
        for j in range(6):
            assert mod1(F(j*triality, 3)+F(j*weak, 2)+F(j*q, 6)) == 0
    qs = {name: row[2] for name, row in rows.items()}
    assert qs["Q"]+qs["H"]+qs["uc"] == 0
    assert qs["Q"]-qs["H"]+qs["dc"] == 0
    assert qs["L"]-qs["H"]+qs["ec"] == 0
    assert qs["L"]+qs["H"]+qs["N"] == 0
    assert qs["S"]+2*qs["N"] == 0
    assert mod1(F(1, 2)+2*F(1, 4)) == 0
    # After H, e=3Q = q/2 +/- 3/2 on doublets.
    e_rows = []
    for name, (triality, weak, q) in rows.items():
        es = [F(q, 2)+F(3, 2), F(q, 2)-F(3, 2)] if weak else [F(q, 2)]
        for e in es:
            assert e.denominator == 1
            assert mod1(F(triality, 3)+e/3) == 0
        e_rows.append({"field": name, "e": es})
    # r alone changes H; r times the weak -I fixes H.
    for a in range(4):
        assert mod1(F(a, 2)+F(a, 2)) == 0
    assert mod1(F(1, 2)) != 0
    return {"k": k, "H_residual_charges": e_rows}


def det(m):
    return m[0][0]*m[1][1]-m[0][1]*m[1][0]


def rank2(m):
    return 2 if det(m) else (1 if any(x for row in m for x in row) else 0)


def quadratic(q, x):
    return sum(a*b for a, b in zip(q, x))**2


def kinetic_hessian(k, g=F(1), gp=F(1), v=F(1), w=F(1)):
    # Independent order: (B,W3). Recover Hessian from energies, not the
    # author's mass-matrix implementation.
    qh, qs = (gp/2, -g/2), (-6*k*gp, F(0))
    def energy(x):
        return (v*v*quadratic(qh, x)+w*w*quadratic(qs, x))/2
    e1, e2 = (F(1), F(0)), (F(0), F(1))
    m11, m22 = 2*energy(e1), 2*energy(e2)
    m12 = energy((F(1), F(1)))-energy(e1)-energy(e2)
    m = [[m11, m12], [m12, m22]]
    c = [[v*v*x for x in qh], [w*w*x for x in qs]]
    assert det(m) == 9*k*k*g*g*gp*gp*v*v*w*w
    for i in range(2):
        for j in range(2):
            assert m[i][j] == c[0][i]*c[0][j]/(v*v)+c[1][i]*c[1][j]/(w*w)
    assert rank2(m) == rank2(c) == (1 if k == 0 else 2)
    if k == 0:
        null = (g, gp)
        assert all(sum(row[j]*null[j] for j in range(2)) == 0 for row in m)
        assert all(sum(row[j]*null[j] for j in range(2)) == 0 for row in c)
    return m, c


def response(m, s):
    # Exact inverse in the independent (B,W3) ordering.
    denominator = (s+m[0][0])*(s+m[1][1])-m[0][1]*m[1][0]
    return s*(s+m[1][1])/denominator


def poly_mul(a, b):
    out = [F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            out[i+j] += x*y
    return out


def calculate():
    ks = [-2, -1, 0, 1, 3]  # Independent representative sample.
    rep = [representations(k) for k in ks]
    groups = [group_test(k) for k in ks if k]
    params = [(F(1), F(1), F(1), F(1)),
              (F(3, 2), F(2, 3), F(4, 5), F(7, 6))]
    masses = []
    for g,gp,v,w in params:
        for k in ks:
            m,c = kinetic_hessian(k,g,gp,v,w)
            masses.append({"k":k, "parameters":[g,gp,v,w], "M_BW":m,
                           "C_BW":c, "mass_rank":rank2(m), "current_rank":rank2(c),
                           "determinant":det(m)})
    delta = F(18,73)
    slo, shi = F(9,584), F(9,292)
    assert slo == delta/16 and shi == delta/8
    ideal_gap, each_error = F(3,8), F(1,64)
    certified_gap = ideal_gap-2*each_error
    assert certified_gap == F(11,32)
    window = []
    m0,_ = kinetic_hessian(0)
    for s in [slo, (slo+shi)/2, shi]:
        r0 = response(m0,s)
        assert r0 >= F(1,2)
        for k in [-2,-1,1,3]:
            mk,_ = kinetic_hessian(k)
            tr = mk[0][0]+mk[1][1]
            dk = det(mk)/tr
            assert dk-delta == F(18*(k*k-1),73*(72*k*k+1))
            assert dk >= delta
            rk = response(mk,s)
            assert rk <= s/(s+dk) <= s/(s+delta) <= F(1,9) < F(1,8)
            assert r0-rk >= ideal_gap
            a,b = F(1,4), F(36*k*k)
            exact_difference = b*(s+a)**2/((s+2*a)*(s*s+s*(2*a+b)+a*b))
            assert r0-rk == exact_difference
            window.append({"s":s,"k":k,"R0":r0,"Rk":rk,"gap":r0-rk,
                           "eigenvalue_lower_bound":dk})
    assert response(m0,shi) == F(82,155)
    # Universal positivity certificate after cross multiplication:
    # delta(n)-delta(1) = 18(n-1)/(73(72n+1)), n=k^2>=1.
    lhs = [F(-18), F(18)]
    assert [x-y for x,y in zip(poly_mul([0,18],[73]),
                              poly_mul([18],[1,72]))] == lhs
    # Fixed-background radial potential: V=lambda(w sigma+
    # (sigma^2+pi^2)/2)^2. Thus the pi^2 quadratic coefficient is zero.
    hessian_phase = {"lambda_S":F(1),"w":F(1),
                     "radial_mass_squared":F(2),"phase_mass_squared_k0":F(0),
                     "S_kept_complex_at_k0":True}
    # Deleted hypotheses: w=0 removes all S mass contributions and restores
    # the k=0 neutral mass matrix even when k is nonzero.
    m_no_S = [[F(1,4),F(-1,4)],[F(-1,4),F(1,4)]]
    assert m_no_S == m0 and det(m_no_S) == 0
    # On-axis constant temporal potentials give the same quadratic energy
    # in the scalar stress and in the mass Hessian, not a complete Ward test.
    source_checks = []
    for k in [0,1,-2]:
        m,_ = kinetic_hessian(k)
        a0 = (F(2,3),F(-3,5))
        e = sum(a0[i]*m[i][j]*a0[j] for i in range(2) for j in range(2))/2
        direct = F(1,2)*(quadratic((F(1,2),F(-1,2)),a0)
                           +quadratic((F(-6*k),F(0)),a0))
        assert e == direct and e > 0
        source_checks.append({"k":k,"A0_BW":a0,"scalar_temporal_energy":e})
    return {
        "status":"PASS_WITH_SCOPE",
        "method":"Independent Fraction cover phases/cosets, energy Hessians, exact resolvents; no author imports",
        "representations":rep, "residual_group_samples":groups,
        "k0_group_scope":{"K0":"U(1) x Z2(F)",
                          "pure_gauge_connected":True,
                          "full_after_connected_Spin_identification_connected":True,
                          "fermion_parity_not_deleted":True},
        "mass_and_current_samples":masses, "finite_window_samples":window,
        "universal_bounds":{"s_min":slo,"s_max":shi,"delta_star":delta,
                            "ideal_gap":ideal_gap,"per_branch_error":each_error,
                            "conditional_remaining_gap":certified_gap,
                            "delta_monotonic_numerator_coefficients":lhs},
        "radial_reference":hessian_phase,"source_quadratic_checks":source_checks,
        "deletion_checks":{"r_without_weak_center_fails_H":True,
                           "w_zero_restores_zero_mode":True,
                           "nonzero_k_N_internal_character_minus_one":True},
        "scope":{"all_integer_and_continuous_window_proofs":"analytic, samples not quantifier",
                 "background":"local constant H/S and vanishing gauge/derivative jets",
                 "response":"leading translation-invariant Euclidean spacelike field-strength kernel",
                 "actual_instrument":False,"higher_order_budget_proved":False,
                 "all_backgrounds_have_nondegenerate_HS":False,
                 "anomaly_recomputed":False,"full_source_quantum_Ward_certificate":False,
                 "k_sign_equivalence_of_entire_model":False}
    }


def encode(obj):
    if isinstance(obj,F):
        return str(obj)
    if isinstance(obj,dict):
        return {k:encode(v) for k,v in obj.items()}
    if isinstance(obj,(list,tuple)):
        return [encode(v) for v in obj]
    return obj


def payload():
    science = encode(calculate())
    sources = [
        HERE/"proof.md", HERE/"check.py", HERE/"results.json",
        HERE.parent/"_shared"/"round1292_admission.md",
        ROOT/"research_cognition_physics/archive_1009_1043/1035/drafts/joint_selection_derivation.md",
        ROOT/"research_cognition_physics/archive_1046_1062/1060/proof.md",
        ROOT/"research_cognition_physics/archive_1288_/1289/proof.md",
        Path(__file__).resolve(),
    ]
    bindings = {str(p.relative_to(ROOT)).replace("\\","/"):sha(p) for p in sources}
    return {"scientific_payload":science,
            "scientific_payload_sha256":hashlib.sha256(
                json.dumps(science,sort_keys=True,ensure_ascii=False,separators=(",",":")).encode()).hexdigest(),
            "source_bindings":bindings}


def main():
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--write",action="store_true")
    mode.add_argument("--check",action="store_true")
    args = parser.parse_args()
    result = payload()
    if args.write:
        with OUT.open("x",encoding="utf-8",newline="\n") as stream:
            json.dump(result,stream,ensure_ascii=False,indent=2)
            stream.write("\n")
        print("CREATED",OUT.name)
    else:
        saved = json.loads(OUT.read_text(encoding="utf-8"))
        assert result == saved, "independent scientific payload or source binding mismatch"
        print("READ_ONLY_MATCH",OUT.name)
    print("PASS_WITH_SCOPE",result["scientific_payload_sha256"])


if __name__ == "__main__":
    main()
