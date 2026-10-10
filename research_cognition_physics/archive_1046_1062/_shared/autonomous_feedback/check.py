"""Count-0 diagnostic of the adopted Koski four-charge-state model.

This is not a fit to the device's nonuniform electron temperatures, nor a
simulation of cotunneling, thermometry, reservoir depletion or quantum memory.
Default: recompute and compare without writes. --save: first exclusive creation.
"""
from pathlib import Path
import argparse
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
E = 1.602176634e-19
KB = 1.380649e-23
STATES = [(0, 0), (0, 1), (1, 0), (1, 1)]


def rate(delta, resistance, temperature):
    x = delta / (KB * temperature)
    return KB * temperature / (E * E * resistance) * (
        1.0 if x == 0.0 else x / math.expm1(x))


def evaluate(voltage, ts=0.070, td=0.070):
    j, rs, rd = KB * 0.350, 580000.0, 43000.0
    energy = np.array([j * (2*n-1)*(2*N-1)/2 for n, N in STATES])
    channels = []
    generator = np.zeros((4, 4))
    for src, (n, N) in enumerate(STATES):
        for label, mu, res, temp in [
                ("L", E*voltage/2, rs, ts),
                ("R", -E*voltage/2, rs, ts), ("D", 0.0, rd, td)]:
            dst = STATES.index((n, 1-N) if label == "D" else (1-n, N))
            dn = 0 if label == "D" else 1-2*n
            dh = float(energy[dst] - energy[src])
            de = dh - mu*dn
            gamma = rate(de, res, temp)
            channels.append(dict(src=src, dst=dst, label=label, dn=dn,
                                 dh=dh, de=de, gamma=gamma, mu=mu, temp=temp))
            generator[dst, src] += gamma
            generator[src, src] -= gamma
    scale = float(np.max(-np.diag(generator)))
    matrix = generator / scale
    matrix[-1] = 1.0
    rhs = np.array([0.0, 0.0, 0.0, 1.0])
    p = np.linalg.solve(matrix, rhs)
    pn = [sum(p[i] for i, s in enumerate(STATES) if s[0] == n) for n in (0, 1)]
    pd = [sum(p[i] for i, s in enumerate(STATES) if s[1] == n) for n in (0, 1)]
    point_information = np.array([math.log(p[i]/(pn[n]*pd[N]))
                                  for i, (n, N) in enumerate(STATES)])
    lookup = {(c["label"], c["src"], c["dst"]): c for c in channels}
    heat = dict(S=0.0, D=0.0)
    h_flow = dict(S=0.0, D=0.0)
    info = dict(S=0.0, D=0.0)
    partial_ep = dict(S=0.0, D=0.0)
    current = dict(L=0.0, R=0.0)
    work, ldb_error = 0.0, 0.0
    for c in channels:
        a, b = c["src"], c["dst"]
        sector = "D" if c["label"] == "D" else "S"
        flux = p[a]*c["gamma"]
        rev = lookup[c["label"], b, a]
        revflux = p[b]*rev["gamma"]
        heat[sector] -= flux*c["de"]
        h_flow[sector] += flux*c["dh"]
        info[sector] += flux*(point_information[b]-point_information[a])
        partial_ep[sector] += 0.5*(flux-revflux)*math.log(flux/revflux)
        ldb_error = max(ldb_error, abs(math.log(c["gamma"]/rev["gamma"])
                                      + c["de"]/(KB*c["temp"])))
        work += c["mu"]*flux*c["dn"]
        if c["label"] in current:
            current[c["label"]] += E*flux*c["dn"]

    # Independent two-sector solution and the paper's equations (3)--(5).
    jp, jm = j+E*voltage/2, j-E*voltage/2
    gs = lambda de: rate(de, rs, ts)
    gd = lambda de: rate(de, rd, td)
    gc, cg = gs(jm)+gs(jp)+gd(j), gs(-jp)+gs(-jm)+gd(-j)
    pg, pc = cg/(gc+cg), gc/(gc+cg)
    p_closed = np.array([pc/2, pg/2, pg/2, pc/2])
    qs_closed = -(jm*gs(jm)+jp*gs(jp))*pg + (jm*gs(-jm)+jp*gs(-jp))*pc
    qd_closed = -j*gd(j)*pg+j*gd(-j)*pc
    i_closed = E/2*((gs(jm)-gs(jp))*pg+(gs(-jp)-gs(-jm))*pc)
    id_closed = math.log(pc/pg)*gd(j)*pg+math.log(pg/pc)*gd(-j)*pc

    # Retaining n,N but suppressing the physical L/R junction label.
    coarse = np.zeros((4, 4))
    for c in channels:
        if c["label"] != "D":
            coarse[c["dst"], c["src"]] += c["gamma"]
    coarse_entropy, coarse_ep = 0.0, 0.0
    for a in range(4):
        for b in range(4):
            if coarse[b, a] > 0:
                f, r = p[a]*coarse[b, a], p[b]*coarse[a, b]
                coarse_entropy += f*math.log(coarse[b, a]/coarse[a, b])
                coarse_ep += 0.5*(f-r)*math.log(f/r)
    entropy_s, entropy_d = heat["S"]/(KB*ts), heat["D"]/(KB*td)
    marginal = np.array([[pn[n]*pd[N] for n, N in STATES]])[0]
    checks = dict(
        stationarity_scaled=float(np.max(abs(generator @ p))/scale),
        column_sum_scaled=float(np.max(abs(generator.sum(axis=0)))/scale),
        closed_probability=float(np.max(abs(p-p_closed))),
        local_detailed_balance=ldb_error,
        energy_ledger_scaled=abs(heat["S"]+heat["D"]-work)/(j*scale),
        stationary_H_scaled=abs(sum(h_flow.values()))/(j*scale),
        current_conservation_scaled=abs(sum(current.values()))/(E*scale),
        power_IV_scaled=abs(work-current["L"]*voltage)/(j*scale),
        information_balance_scaled=abs(sum(info.values()))/scale,
        entropy_S_identity_scaled=abs(partial_ep["S"]-(entropy_s-info["S"]))/scale,
        entropy_D_identity_scaled=abs(partial_ep["D"]-(entropy_d-info["D"]))/scale,
        coarse_entropy_identity_scaled=abs(coarse_ep-(coarse_entropy-info["S"]))/scale,
        closed_heat_S_scaled=abs(heat["S"]-qs_closed)/(j*scale),
        closed_heat_D_scaled=abs(heat["D"]-qd_closed)/(j*scale),
        closed_current_scaled=abs(current["L"]-i_closed)/(E*scale),
        closed_information_scaled=abs(info["D"]-id_closed)/scale)
    assert min(p) > 0 and abs(sum(p)-1) < 1e-14
    assert max(checks.values()) < 1e-12
    assert min(partial_ep.values()) > -1e-8
    assert coarse_ep > -1e-8 and entropy_s-coarse_entropy > -1e-8
    if voltage == 0.0:
        assert max(abs(heat["S"]), abs(heat["D"]), abs(work)) < 1e-28
        assert abs(current["L"]) < 1e-25
        gibbs = np.exp(-energy/(KB*ts))
        gibbs /= sum(gibbs)
        assert np.max(abs(p-gibbs)) < 1e-14
    else:
        assert heat["S"] < 0 < heat["D"]
        assert info["D"] > 0 and entropy_s-coarse_entropy > 0
    return dict(
        voltage_microvolt=voltage*1e6, temperature_S_K=ts, temperature_D_K=td,
        state_order=[list(s) for s in STATES], stationary_probability=p.tolist(),
        stationary_G_probability=pg, current_fA=current["L"]*1e15,
        heat_S_aW=heat["S"]*1e18, heat_D_aW=heat["D"]*1e18,
        supplied_IV_aW=work*1e18,
        interaction_energy_flow_S_aW=h_flow["S"]*1e18,
        interaction_energy_flow_D_aW=h_flow["D"]*1e18,
        mutual_information_nats=float(p @ point_information),
        information_D_per_s=info["D"], information_S_per_s=info["S"],
        heat_entropy_S_per_s=entropy_s, heat_entropy_D_per_s=entropy_d,
        entropy_production_S_per_s=partial_ep["S"],
        entropy_production_D_per_s=partial_ep["D"],
        coarse_visible_entropy_S_per_s=coarse_entropy,
        coarse_entropy_production_S_per_s=coarse_ep,
        hidden_junction_entropy_per_s=entropy_s-coarse_entropy,
        product_marginal_stationarity_residual_scaled=float(
            np.max(abs(generator @ marginal))/scale),
        checks=checks)


def calculate():
    return dict(
        new_science_groups=0, new_empirical_groups=0, new_cognitive_axioms=0,
        scope="Adopted sequential-tunneling model diagnostics; no experiment fit.",
        experimental_fit=False, cotunneling_included=False,
        complete_unknown_quantum_reference_instrument=False,
        finite_reservoir_depletion_model=False,
        constants_SI=dict(e=E, k_B=KB),
        adopted_device_scale=dict(J_over_kB_K=0.350, R_s_ohm=580000,
                                  R_d_ohm=43000, E_s_over_kB_K=1.7,
                                  E_d_over_kB_K=0.810, n_g=0.5, N_g=0.5),
        diagnostic_temperature_choice="Both sectors 0.070 K, not measured-fit temperatures.",
        primary_sources=["https://arxiv.org/html/1507.00530",
                         "https://arxiv.org/pdf/1509.08288"],
        nonequilibrium=evaluate(20e-6), equilibrium_control=evaluate(0.0),
        all_checks_passed=True)


def compare(a, b, path=""):
    if isinstance(a, dict):
        assert a.keys() == b.keys(), path
        for key in a:
            compare(a[key], b[key], path+"/"+key)
    elif isinstance(a, list):
        assert len(a) == len(b), path
        for i, (x, y) in enumerate(zip(a, b)):
            compare(x, y, path+f"/{i}")
    elif isinstance(a, float):
        # Residuals and equilibrium zeros are roundoff-level diagnostics.
        atol = 1e-8 if path.endswith("_per_s") else 1e-10
        assert math.isclose(a, b, rel_tol=2e-11, abs_tol=atol), (path, a, b)
    else:
        assert a == b, (path, a, b)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--save", action="store_true")
    args = parser.parse_args()
    result = calculate()
    target = HERE/"results.json"
    if args.save:
        with target.open("x", encoding="utf-8") as handle:
            json.dump(result, handle, ensure_ascii=False, indent=2, allow_nan=False)
            handle.write("\n")
        print("First diagnostic results saved; adopted science count = 0.")
    else:
        compare(result, json.loads(target.read_text(encoding="utf-8")))
        print("Read-only four-state, energy, information and coarse-channel checks passed.")
