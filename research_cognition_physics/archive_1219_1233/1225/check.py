"""1225: six original full-Gauss record configurations and area dictionary.

Numerics evaluate target-space preparation coefficients, finite-star incidence,
and inverse error transport. They do NOT simulate the full-H propagator, compute
its third moments, or certify an actual clock/area device. The analytic rank
bound holds at the fixed nonzero packet widths independently of quadrature.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parents[1]
PREVIOUS = HERE.parent / "1224" / "check.py"
spec = importlib.util.spec_from_file_location("round1224_coefficients", PREVIOUS)
old = importlib.util.module_from_spec(spec)
spec.loader.exec_module(old)
RESULTS = HERE / "results.json"
SOURCES = [
    "archive_585_628/research_note_589.md",
    "archive_585_628/research_note_598.md",
    "archive_629_652/research_note_644.md",
    "archive_1219_/1221/proof.md",
    "archive_1219_/1224/proof.md",
    "archive_1219_/1224/check.py",
    "archive_1219_/1224/results.json",
    "archive_1219_/1224/NEXT.md",
    "archive_1219_/_admission/after1224_joint_geometry_thermal_screen.md",
]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def analytic_bounds():
    c = old.constants()
    b0, bk = c["gradient_bound"], c["common_center_radius"]
    h, q1, a0 = c["hessian_bound"], c["divergence_L1_bound"], c["A_origin_component_lower"]
    lf = old.M / 4 + math.sqrt(old.M) / (2 * old.R)
    a_lower = b0 * (a0 - bk * h * q1 - lf * bk)
    b_upper = b0 * (h * h * q1 * bk + b0 * h)
    s_lower = c["S_lower"]
    inverse_inf = max(1 / s_lower, (1 + 2 * b_upper / s_lower) / a_lower)
    inverse_frob = math.sqrt(3 / a_lower**2 + 6 * b_upper**2 / (a_lower**2 * s_lower**2) + 3 / s_lower**2)
    assert a_lower > 0 and s_lower > 0 and c["eta"] > 0
    return dict(fixed_probe_width=old.WIDTH, fixed_narrow_width=bk,
                source_amplitude=c["eta"], Hess_F_global_upper=lf,
                a_lower=a_lower, a_upper=b0*b0, abs_b_upper=b_upper,
                S_lower=s_lower, S_upper=b0*b0,
                inverse_infinity_norm_upper=inverse_inf,
                singular_value_lower=1/inverse_frob,
                bounds_are_analytic_not_quadrature=True,
                width_limit_used=False)


def menu_incidence():
    side = 5
    e = np.eye(3, dtype=int)
    origin = (0, 0, 0)
    rows = []
    for kind, directions in [("adjacent", (i,)) for i in range(3)] + [
            ("mixed", (i, j)) for i in range(3) for j in range(i+1, 3)]:
        j = e[directions[0]] if kind == "adjacent" else -e[directions[0]] + e[directions[1]]
        source = tuple((j % side).tolist())
        common = []
        for k in np.ndindex(side, side, side):
            ends = [tuple(((np.array(k)+d) % side).tolist()) for d in e]
            star = [k] + ends
            if origin in star and source in star:
                common.append(dict(center=list(k), probe_role=star.index(origin), source_role=star.index(source)))
        assert len(common) == 1
        if kind == "adjacent":
            assert common[0]["center"] == [0, 0, 0] and common[0]["probe_role"] == 0
        else:
            assert common[0]["center"] == ((-e[directions[0]]) % side).tolist()
        rows.append(dict(kind=kind, directions=[i+1 for i in directions],
                         source=list(source), common=common))
    return dict(configurations=rows, fixed_configurations=6, source_signs_per_configuration=2,
                terminal_effects_per_preparation=2, graph_side=side,
                probe_center_narrow_packet_swapped_by_fixed_recipe=True,
                all_full_H_terms_retained=True)


def preparation_coefficients(order):
    c = old.constants()
    R, M = old.R, old.M
    rr, wr = old.radial_nodes(order, old.WIDTH)
    u, wu = old.legendre(order, -1, 1)
    wu *= (3/4)*(1-u*u)
    j = np.arange(1, order+1)
    ang = np.cos(j*np.pi/(order+1))
    wa = 2*np.sin(j*np.pi/(order+1))**2/(order+1)
    r, polar, rel = np.meshgrid(rr, u, ang, indexing="ij")
    weight = (wr[:, None, None]*wu[None, :, None]*wa[None, None, :]).reshape(-1)
    y = np.zeros((r.size, 5))
    y[:, 0] = (r*np.sqrt(1-polar**2)*rel).reshape(-1)
    y[:, 1] = (r*np.sqrt(1-polar**2)*np.sqrt(1-rel**2)).reshape(-1)
    y[:, -1] = (r*polar).reshape(-1)
    Y = old.embed(y)
    U = old.effect_gradient_ambient(Y)
    rb, wb = old.radial_nodes(order, c["common_center_radius"])
    a = b = S = tangent = 0.
    for radius, rw in zip(rb, wb):
        for polar_b, uw in zip(u, wu):
            B = old.embed(np.array([radius*math.sqrt(1-polar_b**2), 0, 0, 0, radius*polar_b]))
            logs, ders = old.log_and_differential(B, Y, U)
            A = weight @ ders
            L = weight @ logs
            # Remaining S2 transverse Higgs orientations were analytically integrated.
            A[2:5] = 0
            L[2:5] = 0
            grad = old.effect_gradient_ambient(B)
            sqrtD = B[0]/R
            s = math.sqrt(M)*B[-1]/sqrtD
            fold = M/sqrtD**2
            grad_s_sq = fold*(1-s*s/(6*M))
            lap_s = -fold*s/(3*M)
            lap_f = (math.cos(s)*lap_s-math.sin(s)*grad_s_sq)/4
            dr_f = math.cos(s)*math.sqrt(M)/4/math.cosh(radius/R)**2*polar_b
            rf_over_kappa = 4*old.N*radius/(c["common_center_radius"]**2-radius**2)*dr_f-lap_f
            w = float(rw*uw)
            a += w*float(old.mink(grad, A))
            b -= w*rf_over_kappa*float(old.mink(A, L))
            S += w*float(old.mink(A, A))
            tangent = max(tangent, abs(float(old.mink(B, A))), abs(float(old.mink(B, L))))
    s = math.sqrt(M)*Y[:, -1]/(Y[:, 0]/R)
    fold = M/(Y[:, 0]/R)**2
    I_F = float(weight @ (np.cos(s)*fold*(1-s*s/(6*M))))
    v0 = old.HBAR*old.PHASE*I_F/4
    bnd = analytic_bounds()
    assert bnd["a_lower"] <= a <= bnd["a_upper"]
    assert abs(b) <= bnd["abs_b_upper"] and bnd["S_lower"] <= S <= bnd["S_upper"]
    assert tangent < 1e-12 and v0 > 0
    return dict(order=order, a=a, b=b, S=S, I_F=I_F, v0=v0,
                maximum_tangent_residual=tangent)


def matrix_from(row):
    a, b, S = (row[k] for k in ("a", "b", "S"))
    incidence = np.array([[1,1,0],[1,0,1],[0,1,1]], dtype=float)
    out = np.zeros((6,6))
    out[:3,:3] = a*np.eye(3)
    out[:3,3:] = b*incidence
    out[3:,3:] = -S*np.eye(3)
    return out


def stable_inverse(row):
    L = matrix_from(row)
    inv = np.linalg.inv(L)
    a,b,S = (row[k] for k in ("a", "b", "S"))
    determinant = -a**3*S**3
    assert abs(np.linalg.det(L)-determinant) < 1e-17
    bounds = analytic_bounds()
    smin = float(np.linalg.svd(L, compute_uv=False)[-1])
    assert smin >= bounds["singular_value_lower"]
    assert np.linalg.norm(inv, np.inf) <= bounds["inverse_infinity_norm_upper"]
    # Two previously used shapes, plus one strictly diagonal equal-determinant
    # shape: demonstration of a fixed inverse, not a parameter scan or full H.
    shapes = [np.eye(3), np.array([[4/3,2/3,1/2],[2/3,4/3,0],[1/2,0,1]]),
              np.diag([2.,.5,1.])]
    volumes = [.6, .8, 1.2]
    examples=[]
    for G,w in zip(shapes,volumes):
        C=np.linalg.inv(G)
        P=w**(4/3)*C
        p=np.array([P[0,0],P[1,1],P[2,2],P[0,1],P[0,2],P[1,2]])
        z=L@p
        recovered=inv@z
        assert np.max(abs(recovered-p)) < 1e-14
        # Real relational coefficients are eta/v0^2 times these z.
        raw=old.constants()["eta"]/row["v0"]**2*z
        A=np.sqrt(np.diag(P))
        h=w**(2/3)*G
        independent=np.array([math.sqrt(np.linalg.det(h[np.ix_([j for j in range(3) if j!=i], [j for j in range(3) if j!=i])])) for i in range(3)])
        assert np.max(abs(A-independent)) < 1e-14
        examples.append(dict(volume=w, shape=G.tolist(), normalized_records=z.tolist(),
                             relational_quadratic_records=raw.tolist(), face_areas=A.tolist(),
                             inverse_error=float(np.max(abs(recovered-p)))))
    return dict(matrix=L.tolist(), determinant=determinant,
                singular_value_minimum=smin, inverse_infinity_norm=float(np.linalg.norm(inv,np.inf)),
                examples=examples, areas_are_independent_of_entropy=True,
                examples_are_algebraic_not_full_H_simulations=True)


def error_contract(row):
    # A finite numeric inverse-transport check ONLY. The probability certificate
    # still needs full-H moment bounds; these are not assigned arbitrary values.
    lower, upper, gminus, gplus = .5, 1.5, .25, 4.
    pmin=lower**(4/3)/gplus
    pmax=upper**(4/3)/gminus
    K=analytic_bounds()["inverse_infinity_norm_upper"]
    delta_L=1e-6
    delta_z=1e-7
    delta_P=K/(1-K*delta_L)*(delta_z+delta_L*pmax)
    assert K*delta_L<1 and delta_P<pmin
    area_error=delta_P/(math.sqrt(pmin-delta_P)+math.sqrt(pmin))
    L=matrix_from(row)
    Lhat=L+delta_L/6*np.ones((6,6))
    P=np.array([.8,.9,1.1,.1,-.1,.2])
    zhat=L@P+delta_z*np.array([1,-1,1,-1,1,-1])
    phat=np.linalg.solve(Lhat,zhat)
    actual=float(np.max(abs(phat-P)))
    assert actual<=delta_P
    return dict(declared_compact_range=dict(w_min=lower,w_max=upper,g_min=gminus,g_max=gplus),
                illustrative_coefficient_error=delta_L,illustrative_normalized_record_error=delta_z,
                cofactor_entry_error_bound=delta_P,area_error_bound=area_error,
                inverse_perturbation_check=actual,
                coefficient_errors_are_assumed_not_interval_certified=True,
                actual_probability_error_contract="(v0^2/eta)[2 R_rel d + (2 epsilon_p + 4 M1_B epsilon_c/v_min)/d^2]",
                dictionary_error_contract="K/(1-K delta_L) [delta_z+delta_L P_max]",
                full_H_moments_uniformly_finite=True,full_H_moments_numerically_computed=False,
                specific_wait_time_or_sample_count_certified=False,
                fixed_finite_waiting_grid_uniform_Hoeffding=True,
                absolute_meter_calibration_or_autonomous_clock=False)


def run():
    rows=[preparation_coefficients(n) for n in (12,18,24)]
    differences={k:abs(rows[-1][k]-rows[-2][k]) for k in ("a","b","S","I_F")}
    assert max(differences.values())<2e-10
    return dict(round=1225,scientific_groups=4,source_sha256={p:sha(RESEARCH/p) for p in SOURCES},
                menu=menu_incidence(),fixed_width_rank_bounds=analytic_bounds(),
                target_integrals=dict(rows=rows,last_refinement_difference=differences,
                                     not_interval_quadrature=True,not_full_H_evolution=True),
                stable_dictionary=stable_inverse(rows[-1]),
                finite_error_contract=error_contract(rows[-1]),
                scope="Six full-H finite-width geometric record configurations; stable cofactor/area inverse conditional on declared task resources")


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--write",action="store_true")
    parser.add_argument("--verify",action="store_true")
    args=parser.parse_args()
    data=run()
    data["code_sha256"]=sha(Path(__file__))
    if args.write:
        RESULTS.write_text(json.dumps(data,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
        print("1225: wrote 4 groups; fixed-width six-row inverse, no device budget claimed")
    else:
        saved=json.loads(RESULTS.read_text("utf-8"))
        assert saved==data,"1225 saved result, source or code differs"
        print("1225: PASS 4 groups; full-H moment values remain uncomputed")


if __name__=="__main__":
    main()
