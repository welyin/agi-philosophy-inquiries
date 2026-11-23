"""Boundary-value charged-flavour transport for two split complete fives."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import runpy

import numpy as np

HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE / path).read_bytes())


CFG = read("inputs.json")
BASE = runpy.run_path(str(HERE / CFG["gauge_and_low_inputs"]))
PION = runpy.run_path(str(HERE / CFG["pion_response"]))
LOOP = 16 * math.pi**2
SHAPES = ((3, 3), (3, 5), (5, 3), (2, 5), (5, 2), (2, 2))
BOUNDARIES = np.cumsum([np.prod(s) for s in SHAPES])[:-1]


def pack(*matrices):
    return np.concatenate([m.ravel() for m in matrices])


def unpack(state):
    return [a.reshape(s) for a, s in zip(np.split(state, BOUNDARIES), SHAPES)]


def dagger(a):
    return a.conj().T


def norm2(a):
    return float(np.vdot(a, a).real)


def ckm_params(v):
    s13 = abs(v[0, 2])
    c13 = math.sqrt(1 - s13**2)
    s12, s23 = abs(v[0, 1]) / c13, abs(v[1, 2]) / c13
    c12, c23 = np.sqrt(1 - np.array([s12, s23])**2)
    j = (v[0, 0]*v[1, 1]*v[0, 1].conjugate()*v[1, 0].conjugate()).imag
    sine = j/(s12*c12*s23*c23*s13*c13**2)
    cosine = (abs(v[1, 0])**2-s12**2*c23**2-c12**2*s23**2*s13**2)/(2*s12*c12*s23*c23*s13)
    return np.array([s12, s13, s23, math.atan2(sine, cosine)])


def null_and_row(row):
    right = dagger(np.linalg.svd(row, full_matrices=True)[2])
    return right[:, row.shape[0]:], right[:, :row.shape[0]]


def initial(params):
    ys = np.exp(params[:9])
    mix = BASE["FLOW"]["mixing"](np.arcsin(np.exp(params[9:12])), params[12])
    down = mix * ys[3:6]
    yb, ye, ym, yt = ys[5:9]
    sm = math.sqrt((yb**-2-yt**-2)/(ym**-2-yt**-2))
    cm = math.sqrt(1 - sm**2)
    el = mix @ np.array([[1, 0, 0], [0, cm, -sm], [0, sm, cm]]) @ np.diag([ye, ym, yt])
    v = np.linalg.solve(el, down[:, 2])
    assert abs(norm2(v)-1) < 1e-11
    k, _ = null_and_row(v.conj()[None, :])
    color = np.column_stack((down, el@k))
    weak = np.column_stack((el, down[:, :2])).T
    eye = np.eye(5)
    weak_light = eye[:, 2:3] @ dagger(v[:, None]) + eye[:, 3:] @ dagger(k)
    weak_basis = np.column_stack((weak_light, eye[:, :2]))
    assert np.linalg.norm(color @ weak_basis - weak.T) < 1e-12
    r = eye[3:, :]
    return pack(np.diag(ys[:3]), color, weak, r, r.T, np.eye(2)), mix, weak_basis


def descend(label, params, step):
    data = BASE["inputs"](label)
    point = read(CFG["two_pair_spectra"])[label]["two_split_complete_pairs"]
    fields = BASE["thresholds"](point, data["MZ"])
    high, t_color, t_weak = np.log(np.array([point["matching_scale_GeV"],
        point["M_T_GeV"][0], point["M_D_GeV"][0]])/data["MZ"])
    y, mix, weak_basis = initial(params)
    damp = BASE["FLOW"]["GAUGE_DAMPING"]

    def rhs(t, state):
        u, d, e, rt, rd, higgs_transfer = unpack(state)
        uu, dd, ee = u@dagger(u), d@dagger(d), e@dagger(e)
        trace = 3*(norm2(u)+norm2(d))+norm2(e)
        g = BASE["couplings"](data, fields, t)
        ku, kd, ke = trace-damp@g
        weak, hyper = g[1], .6*g[0]
        fourth = 3*(norm2(uu)+norm2(dd))+norm2(ee)
        # The fractional linear transfer integrates the full Higgs Riccati equation.
        linear = 4*trace-9*weak-3*hyper
        constant = -2*fourth+3*(2*weak**2+(weak+hyper)**2)/8
        transfer = np.array([[linear, constant], [-24, 0]])@higgs_transfer
        return pack(1.5*(uu-dd)@u+ku*u, 1.5*(dd-uu)@d+kd*d, 1.5*ee@e+ke*e,
                    rt@dagger(d)@d, .5*ee@rd, transfer)/LOOP

    portals = []
    for a, b, species in [(high, t_color, "color"), (t_color, t_weak, "weak"), (t_weak, 0, None)]:
        n = math.ceil((a-b)/step)
        h = (b-a)/n
        for i in range(n):
            t = a+i*h
            k1 = rhs(t, y)
            k2 = rhs(t+h/2, y+h*k1/2)
            k3 = rhs(t+h/2, y+h*k2/2)
            k4 = rhs(t+h, y+h*k3)
            y = y+h*(k1+2*k2+2*k3+k4)/6
        u, d, e, rt, rd, transfer = unpack(y)
        if species == "color":
            t_light, heavy = null_and_row(rt)
            portals.append(norm2(d@heavy))
            physical = d@t_light
            d.fill(0)
            d[:, :3] = physical
            rt.fill(0)
        if species == "weak":
            l_light, heavy = null_and_row(dagger(rd))
            portals.append(norm2(dagger(heavy)@e))
            physical = dagger(l_light)@e
            e.fill(0)
            e[:3] = physical
            rd.fill(0)
    ul, uv, uvr = np.linalg.svd(u)
    dl, dv, dvr = np.linalg.svd(d[:, :3])
    _, ev, evr = np.linalg.svd(e[:3])
    v = dagger(ul[:, ::-1])@dl[:, ::-1]
    obs = np.r_[np.log(np.r_[uv[::-1], dv[::-1], ev[::-1]]), ckm_params(v)]
    obs[9:12] = np.log(obs[9:12])
    # Fixed gauge-basis source only, external-line projections with running masses.
    external = dagger(weak_basis@l_light.conj())@t_light@dagger(dvr)[:, ::-1]
    # The 10 contains left-handed Q as well as u^c and e^c.
    up, down = ul[:, -1], dl[:, -1]
    uc = uvr.T[:, -1]
    ec = evr.T[:, ::-1]
    pion = dagger(ec)@(np.vdot(uc, up)*down+np.vdot(uc, down)*up)
    (a, b), (c, d) = transfer.real
    high_quartic = (b-d*data["lam"])/(c*data["lam"]-a)
    return obs, portals, np.sum(abs(external)**2, axis=0), high_quartic, abs(pion)**2


def shoot(label, step):
    data = BASE["inputs"](label)
    low = ckm_params(data["v"])
    target = np.r_[np.log(data["y"]), np.log(low[:3]), low[3]]
    old = read(CFG["initial_guess_only"])[label]["two_split_complete_pairs"]["point"]
    hp = np.r_[.5*np.log(old["squared_yukawas_at_Q"]), target[9:]]
    inv = np.eye(13)
    previous = None
    for _ in range(16):
        observed, portals, external, quartic, pion = descend(label, hp, step)
        err = observed-target
        if max(abs(err)) < 2e-10:
            return hp, portals, external, max(abs(err)), quartic, pion
        if previous:
            x0, r0 = previous
            change = hp-x0
            delta = err-r0
            inv += np.outer(change-inv@delta, delta)/(delta@delta)
        previous = hp.copy(), err.copy()
        hp -= inv@err
    raise ValueError(f"Charged fit residual {err}")


def point(label):
    h = CFG["step_log_mu"]
    coarse = shoot(label, h)[0]
    p, portals, external, error, quartic, pion = shoot(label, h/2)
    assert max(abs(p-coarse)) < 2e-8
    old = read(CFG["initial_guess_only"])[label]["two_split_complete_pairs"]["point"]
    ys = np.exp(p[:9])
    yb, ye, ym, yt = ys[5:9]
    pom = [(yb**-2-yt**-2)/(a**-2-yt**-2) for a in (ye, ym)]
    p_tau = 1-pom[1]
    _, mix, _ = initial(p)
    q = math.sqrt(1+3*abs(mix[0, 0])**2)
    overlap = abs(mix[0, 0]*mix[0, 2])
    tau_bound = math.sqrt(p_tau)*overlap + math.sqrt(1-p_tau)*math.sqrt(q*q-overlap**2)
    strength = q*q-tau_bound**2
    mass_bounds, example_bounds = {}, {}
    for tag, envelope in [("central", False), ("one_width_shape_box", True)]:
        k = np.array([PION["response"](i, envelope)[0] for i in range(2)])
        widths = CFG["hbar_GeV_s"]/(CFG["seconds_per_year"]*np.array(CFG["lifetime_lower_years_90"]))
        mass_bounds[tag] = (strength/sum(widths/k))**.25
        example_bounds[tag] = max((pion[:2]*k/widths)**.25)
    gram = [sum(ys[6:]**2)-yb**2, sum(ys[3:5]**2)]
    light_trace = 3*sum(ys[:6]**2)+sum(ys[6:]**2)
    gauge_ratio = math.sqrt(old["g_squared_at_Q"][0])*mass_bounds["one_width_shape_box"]/old["Q_GeV"]
    return dict(Q_GeV=old["Q_GeV"], D_T_poles_GeV=[t[0] for t in old["physical_D_T_GeV"]],
        high_yukawas=ys.tolist(), high_ckm_sines_delta=np.r_[np.exp(p[9:12]),p[12]].tolist(),
        maximum_Y_change_from_no_portal=float(max(abs(ys/np.sqrt(old["squared_yukawas_at_Q"])-1))),
        down_aligned_inverse_norm=yb/yt, up_aligned_inverse_norm=np.linalg.norm(yb*mix[:, 2]/ys[6:]),
        ten_e_mu_max_probabilities=pom, ten_tau_min=p_tau, fivebar_tau_max=(yb/yt)**2,
        ten_pion_e_mu_strength_min=strength, Lambda_geV_lower=mass_bounds,
        MX_Q_lower_if_A2_1_same_g=gauge_ratio,
        minimum_relative_same_chirality_vector_change_if_MX_Q_and_A2_1=1-gauge_ratio**-2,
        high_heavy_T_D_Higgs_gram_traces=gram, at_threshold_heavy_T_D_Higgs_gram_traces=portals,
        fixed_high_current_low_external_10_pion_e_mu_tau=pion.tolist(),
        low_external_witness_Lambda_lower=example_bounds,
        high_quartic_active=quartic, change_from_no_portal_quartic=quartic-old["lambda_observed"],
        residual_if_equal_kinetic_trace=old["S"]-light_trace-3*gram[0]-gram[1],
        threshold_summed_partial_widths_over_mass=[portals[0]/(16*math.pi),portals[1]/(32*math.pi)],
        fixed_high_current_low_external_d_s_b=external.tolist(),
        fit_residual_upper=2e-10, verified_log_parameter_step_error_upper=2e-8)


def verify(a, b):
    if isinstance(b, dict):
        assert a.keys() == b.keys()
        for k in b:
            verify(a[k], b[k])
    elif isinstance(b, list):
        assert len(a) == len(b)
        for x, y in zip(a, b):
            verify(x, y)
    else:
        assert math.isclose(a, b, rel_tol=1e-7, abs_tol=1e-13), (a, b)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    outfile = HERE/"results.json"
    values = {tag:point(tag) for tag in CFG["gauge_fits"]}
    if parser.parse_args().check:
        verify(values, json.loads(outfile.read_bytes()))
    else:
        outfile.write_bytes((json.dumps(values, indent=2, allow_nan=False)+"\n").encode())
    print("1811 main", hashlib.sha256(outfile.read_bytes()).hexdigest()[:12])
