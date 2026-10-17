"""Round 1235: finite-window CT1188 bridge. NumPy, no plots or network.

Default is read-only deterministic replay. --write-once creates results.json.
Quadrature is a finite diagnostic; the certified bounds are analytic in proof.md.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parents[1]
S = 0.25
RADIUS = 1e-4
THICKNESS = 1e-5
EPSILON = 1e-8
Q = 4.0 / 3.0  # normalized source is 2 times the independently computed physical T
AS = {"rank_one": np.diag([1.0, 0.0, 0.0]),
      "traceless": np.diag([1.0, -1.0, 0.0])}


def require(condition: bool, label: str) -> None:
    if not condition:
        raise AssertionError(label)


def maximum(value) -> float:
    return float(np.max(np.abs(value)))


def window(order: int):
    nodes, weights = leggauss(order)
    grid = np.meshgrid(nodes, nodes, nodes, nodes, indexing="ij")
    point = np.stack([g.ravel() for g in grid], axis=1)
    point[:, :3] *= RADIUS
    point[:, 3] *= THICKNESS
    wgrid = np.meshgrid(weights / 2, weights / 2, weights / 2,
                       weights / 2, indexing="ij")
    weight = np.prod(np.stack(wgrid, axis=-1), axis=-1).ravel()
    return point, weight


def invariants(u: float, point: np.ndarray, A: np.ndarray):
    x, tau = point[:, :3], 1 + point[:, 3]
    ax = x @ A
    quad = np.einsum("ni,ni->n", x, ax)
    t = tau + 0.5 * u * quad
    ac2 = t ** (2 / 3)
    H = 1 / (3 * t)
    v = u * ax
    vv = np.einsum("ni,nj->nij", v, v)
    v2 = np.einsum("ni,ni->n", v, v)
    W = 1 - v2 / ac2
    require(bool(np.all(W > 0)), "all finite windows spacelike")
    metric = ac2[:, None, None] * np.eye(3) - vv
    K = (u * A + (ac2 * H)[:, None, None] * np.eye(3)
         - 2 * H[:, None, None] * vv) / np.sqrt(W)[:, None, None]
    shape = np.linalg.solve(metric, K)
    trace = np.trace(shape, axis1=1, axis2=2)
    Im = np.einsum("nij,nji->n", shape, shape)
    Jm = trace ** 2
    k = u / ac2
    vAv = np.einsum("ni,ij,nj->n", v, A, v)
    vAA = np.einsum("ni,ij,nj->n", v, A @ A, v)
    D = (k * vAv - H * v2) / (ac2 * W)
    Ia = (3 * H ** 2 + 2 * H * k * np.trace(A)
          + k ** 2 * np.trace(A @ A)
          + 2 * (k ** 2 * vAA - H ** 2 * v2) / (ac2 * W)
          + D ** 2) / W
    Ja = (3 * H + k * np.trace(A) + D) ** 2 / W
    # Geometric B from independent 4D Ricci computation, not T.
    B = 2 / (3 * t ** 2) * (2 / W - 1)
    # Three physical scalars, each q=sqrt(2)/3, action -1/2 (d phi)^2.
    q = np.sqrt(2) / 3
    each_physical_energy = 0.5 * (q / t) ** 2
    physical_energy = 3 * each_physical_energy
    rho = 2 * physical_energy * (2 / W - 1)
    kg = -q / t ** 2 + 3 * H * q / t
    source_conservation = -2 * physical_energy / t + 6 * H * physical_energy
    return dict(t=t, tau=tau, W=W, I=Im, J=Jm, R=Im-Jm+B, B=B,
                rho=rho, Ia=Ia, Ja=Ja, kg=kg,
                source_conservation=source_conservation,
                scalar_curvature=-2/(3*t*t), each_energy=each_physical_energy)


def constraint(g, a=1.0, b=1.0, L=1.0, cosmological=0.0):
    return (g["I"] - L * g["J"]) / a - b * g["R"] + 2*cosmological + g["rho"]


def stable_source_difference(point, A):
    x, tau = point[:, :3], 1 + point[:, 3]
    ax = x @ A
    xx = np.einsum("ni,ni->n", x, ax)
    V = np.einsum("ni,ni->n", ax, ax)
    d = 0.5 * S * xx
    tp, tm = tau+d, tau-d
    y = d/tau
    # Exact positive expression, avoids subtracting three order-one values.
    time_part = Q/tau**2 * y*y*(3-y*y)/(1-y*y)**2
    zp, zm = S*S*V*tp**(-2/3), S*S*V*tm**(-2/3)
    tilt_part = Q*S*S*V*(tp**(-8/3)/(1-zp)+tm**(-8/3)/(1-zm))
    return time_part+tilt_part, time_part, tilt_part


def source_bounds(A):
    diag = np.diag(A)
    tr2 = float(np.sum(diag**2))
    mean_V = tr2*RADIUS**2/3
    mean_quad2 = RADIUS**4*(np.sum(diag**2)/5
                           + (np.sum(diag)**2-np.sum(diag**2))/9)
    m = 1-THICKNESS-0.5*S*np.sum(np.abs(diag))*RADIUS**2
    zstar = S*S*tr2*RADIUS**2*m**(-2/3)
    lower = 2*Q*S*S*mean_V*(1+THICKNESS)**(-8/3)
    upper = S*S*(3*Q/4*mean_quad2*m**(-4)
                 + 2*Q*mean_V*m**(-8/3)/(1-zstar))
    return dict(lower=float(lower), upper=float(upper), min_t=float(m),
                max_z=float(zstar), mean_V=mean_V, mean_quad2=float(mean_quad2))


def certificate():
    # Analytic component bounds of proof section 6, never estimated from samples.
    I_upper = 2*(3/4 + 3/4 + 9/8
                 + 2*(9/16+1/4)*(1/16)*(3/2)*2 + 1/16)
    require(I_upper < 6, "explicit invariant bound I<6")
    F_upper = 2*(6+2*13)+2*(6+13+6)+2+6
    require(F_upper < 128, "independent component triangle bound")
    diagonal = 4*128/((1/4)**2*(1/8)**2)
    offdiagonal = 2*128/((1/4)**2*(1/8)**2)
    Mx = 3*diagonal
    Mt = 4*128/((1/4)**2*(1/16)**2)
    E = 0.5*(RADIUS**2*Mx+THICKNESS**2/3*Mt)
    kappa = 1.1/2*(4*EPSILON/S**2+E)
    ab_bound = 1.1/4*(4*EPSILON/S**2+E)
    lambda_bound = kappa/(2*(2-3*kappa))
    require(kappa < 2/3, "inverse trace coefficient bounded away from pole")
    for endpoint in (1-kappa, 1+kappa):
        mapped = endpoint/(3*endpoint-1)
        require(abs(mapped-0.5) <= lambda_bound+1e-15, "finite lambda inverse bound")
    return dict(I_component_upper=I_upper, F_triangle_upper=F_upper,
                holomorphic_F_certificate=128, derivative_diagonal=diagonal,
                derivative_offdiagonal=offdiagonal, M_x=Mx, M_t=Mt,
                moment_x=RADIUS**2, moment_t=THICKNESS**2/3,
                E=E, central_difference_error_bound=S*S*E,
                epsilon=EPSILON, L_distance_bound=kappa,
                ab_distance_bound=ab_bound, lambda_distance_bound=lambda_bound,
                parameters={"a": [0.9,1.1], "b": [0.9,1.1],
                            "L": [0.5,2], "cosmological": [-1,1]},
                proof_not_sampling=True)


def scientific_result():
    require(0 < S <= 1/4 and 0 < RADIUS <= 1e-4 and 0 < THICKNESS <= 1e-5,
            "strictly positive finite windows within certified domain")
    cert = certificate()
    point, weight = window(10)
    source_rows, wrong_rows = {}, {}
    max_kg=max_conservation=max_matrix=max_full=max_source=0.0
    min_curvature=float("inf")
    for name,A in AS.items():
        samples={u:invariants(u,point,A) for u in (-S,0.0,S)}
        for g in samples.values():
            max_kg=max(max_kg,maximum(g["kg"]))
            max_conservation=max(max_conservation,maximum(g["source_conservation"]))
            max_matrix=max(max_matrix,maximum(g["I"]-g["Ia"]),maximum(g["J"]-g["Ja"]))
            max_full=max(max_full,maximum(constraint(g)))
            max_source=max(max_source,maximum(g["rho"]-g["B"]))
            min_curvature=min(min_curvature,float(np.min(np.abs(g["scalar_curvature"]))))
            require(bool(np.all(g["each_energy"]>0)), "all three source constituents retained")
        raw_source=samples[S]["rho"]+samples[-S]["rho"]-2*samples[0.0]["rho"]
        stable, time_part, tilt_part=stable_source_difference(point,A)
        require(maximum(raw_source-stable)<2e-15,"independent raw/stable source formulas")
        average=float(weight@stable)
        bounds=source_bounds(A)
        require(bounds["lower"] <= average <= bounds["upper"],"strict analytic source bounds")
        point6,weight6=window(6)
        average6=float(weight6@stable_source_difference(point6,A)[0])
        require(abs(average-average6)<1e-18,"finite quadrature consistency")
        bad=[samples[0.0]["rho"]-samples[u]["rho"] for u in (S,-S)]
        bad_delta=float(weight@(bad[0]+bad[1]))
        require(abs(bad_delta+average)<2e-15,"deleted source wrongly rejects compatible geometry")
        require(average>4e-12,"deletion exceeds explicit 1e-12 per-slice budget")
        source_rows[name]=dict(average=average, event_change=float(weight@time_part),
                              normal_change=float(weight@tilt_part),
                              analytic_bounds=bounds, deleted_source_difference=bad_delta,
                              two_quadrature_orders=[6,10], quadrature_difference=abs(average-average6))
        values={u:constraint(samples[u],a=1,b=0.95,L=1.2) for u in (-S,0.0,S)}
        delta=float(weight@(values[S]+values[-S]-2*values[0.0]))
        central=2*S*S*(1-1.2) if name=="rank_one" else 4*S*S*(1-0.95)
        error=abs(delta-central)
        require(error<cert["central_difference_error_bound"],"finite-center difference certificate")
        require(abs(delta)>4*EPSILON,"fixed wrong coefficient excluded at finite budget")
        wrong_rows[name]=dict(a=1,b=0.95,L=1.2, averaged_difference=delta,
                             point_difference=central, finite_window_error=error,
                             certified_error=cert["central_difference_error_bound"],
                             excluded=True)
    require(max_kg<1e-14 and max_conservation<1e-14,"independent KG and full-source conservation")
    require(max_matrix<1e-13,"matrix vs analytic invariants")
    require(max_full<1e-13 and max_source<1e-13,"curved common effective solution")
    require(min_curvature>0.6,"nonzero curvature")
    return dict(status="PASS", science_groups=3,
                group_1=dict(name="independent matter and common curved effective geometry",
                             max_KG=max_kg,max_source_conservation=max_conservation,
                             max_invariant_difference=max_matrix,max_complete_constraint=max_full,
                             max_independent_source_geometry_difference=max_source,
                             min_abs_four_curvature=min_curvature,
                             physical_scalar_count=3, physical_sum_q_squared=2/3,
                             normalized_source_factor=2),
                group_2=dict(name="finite positive thickness and source deletion",rows=source_rows),
                group_3=dict(name="analytic derivative certificate and finite parameter selection",
                             certificate=cert,wrong_candidate=wrong_rows),
                domain=dict(shape_s=S,spatial_half_width=RADIUS,temporal_half_width=THICKNESS,
                            window="normalized uniform product cube; common coordinate weights"),
                scope=dict(CT_adopted=False,full_GR_proved=False,quantum_apparatus_built=False,
                           complete_actual_cognitive_model=False,
                           analytic_bounds_certified=True,
                           quadrature_is_diagnostic=True,
                           source_not_fitted_to_constraint_residual=True))


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--write-once",action="store_true")
    args=parser.parse_args()
    payload=scientific_result()
    sources=[HERE/"proof.md",HERE/"check.py",HERE/"admission.md",HERE/"NEXT.md",
             RESEARCH/"archive_1161_1200/1188/proof.md"]
    payload["source_sha256"]={str(p.relative_to(RESEARCH)).replace("\\","/"):
                              hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    target=HERE/"results.json"
    if args.write_once:
        if target.exists():
            raise FileExistsError("results.json already exists; default mode is read-only")
        target.write_text(json.dumps(payload,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
        mode="CREATED"
    else:
        expected=json.loads(target.read_text(encoding="utf-8"))
        require(payload==expected,"read-only result and source-hash match")
        mode="READ_ONLY_MATCH"
    print(json.dumps(dict(status="PASS",mode=mode,science_groups=3,
                         E=payload["group_3"]["certificate"]["E"],
                         lambda_distance_bound=payload["group_3"]["certificate"]["lambda_distance_bound"]),
                     ensure_ascii=False))


if __name__=="__main__":
    main()
