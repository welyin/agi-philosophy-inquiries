"""Round 1557. Signed flavor transport with running thresholds, Python + NumPy."""
import argparse
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
HIGH, LOW = 1e16, 91.1876
LOOP = 16 * np.pi**2
GAUGE_C = np.array([[.85, 2.25, 8], [.25, 2.25, 8],
                    [2.25, 2.25, 0], [.45, 2.25, 0]])


def load(path):
    return json.loads(path.read_text(encoding="utf-8"))


def modes(matrix):
    left, values, _ = np.linalg.svd(matrix)
    return values[::-1], left[:, ::-1]


def ratio(state):
    d = modes(state[1])[0]
    e = modes(state[2])[0]
    return d[0] / d[1], e[1] / e[0], d[0] / d[1] * e[1] / e[0]


def model(t, singlets):
    h = np.diag(t["H_diag"])
    f = np.zeros_like(h, dtype=complex)
    ix = np.triu_indices(3)
    f[ix] = np.array(t["F_real_upper"]) + 1j * np.array(t["F_imag_upper"])
    f = f + np.triu(f, 1).T
    s = complex(t["s_real"], t["s_imag"])
    r = t["r"]
    v = 174.
    state = np.zeros((6, 3, 3), dtype=complex)
    state[:3] = [r*(h+s*f), h+f, h-3*f]
    if singlets:
        state[3] = r*(h-3*s*f)
        state[4] = -v**2*r**2*1e9*f/(t["t_eV"]*HIGH)
    return state


def rhs(log_down, state, g_mz):
    # Rows are left-handed weak flavors, columns are right-handed mass ports.
    u, d, e, n, mass, _ = state
    gram = state[:4] @ state[:4].conj().transpose(0, 2, 1)
    trace = np.trace(3*(gram[0]+gram[1])+gram[2]+gram[3]).real
    log_up = np.log(HIGH/LOW)-log_down
    g2 = 1/(1/g_mz**2 - np.array([41/10, -19/6, -7])*2*log_up/LOOP)
    deriv = np.zeros_like(state)
    deriv[:4] = ((GAUGE_C@g2-trace)[:, None, None]*state[:4]
                 -1.5*(gram-gram[[1, 0, 3, 2]])@state[:4])/LOOP
    right = n.conj().T @ n
    deriv[4] = -(right.T@mass+mass@right)/LOOP
    yd, ld = modes(d)
    ye, le = modes(e)
    up_d = np.diag(ld.conj().T@gram[0]@ld).real
    nu_e = np.diag(le.conj().T@gram[3]@le).real
    deriv[5].flat[:4] = 1.5/LOOP*np.array([
        yd[1]**2-yd[0]**2, up_d[0]-up_d[1],
        ye[0]**2-ye[1]**2, nu_e[1]-nu_e[0]])
    return deriv


def rk5(fun, x, y, h):
    # Cash--Karp fifth-order method. Its six stages differ from classical RK4.
    k1 = fun(x, y)
    k2 = fun(x+h/5, y+h*k1/5)
    k3 = fun(x+3*h/10, y+h*(3*k1+9*k2)/40)
    k4 = fun(x+3*h/5, y+h*(3*k1/10-9*k2/10+6*k3/5))
    k5 = fun(x+h, y+h*(-11*k1/54+5*k2/2-70*k3/27+35*k4/27))
    k6 = fun(x+7*h/8, y+h*(1631*k1/55296+175*k2/512+575*k3/13824
                            +44275*k4/110592+253*k5/4096))
    return y+h*(37*k1/378+250*k3/621+125*k4/594+512*k6/1771)


def cross(state, active):
    mass = state[4, :active, :active]
    u, m, _ = np.linalg.svd(mass)
    takagi = u.conj()
    diag = np.diag(takagi.T@mass@takagi)
    takagi = takagi*np.exp(-.5j*np.angle(diag))
    np.testing.assert_allclose(takagi.T@mass@takagi, np.diag(m), atol=1e-16)
    n = state[3, :, :active]
    # Check the coefficient of the light Weinberg operator at the instant
    # of matching; its further running does not enter the dimension-four flow.
    pre = n @ np.linalg.solve(mass, n.T)
    turned = n@takagi
    heavy = np.outer(turned[:, 0], turned[:, 0])/m[0]
    kept = (turned[:, 1:] / m[1:]) @ turned[:, 1:].T
    mismatch = np.linalg.norm(pre-heavy-kept)/np.linalg.norm(pre)
    assert mismatch < 5e-9
    state[3:5] = 0.
    state[3, :, :active-1] = turned[:, 1:]
    state[4, :active-1, :active-1] = np.diag(m[1:])
    return float(m[0]*HIGH)


def transport(tensors, g_mz, singlets, cells):
    y = model(tensors, singlets)
    initial = ratio(y)
    active = 3 if singlets else 0
    thresholds = []
    x, end = 0., np.log(HIGH/LOW)
    fun = lambda position, state: rhs(position, state, g_mz)
    while x < end-1e-12:
        h = min(end/cells, end-x)
        trial = rk5(fun, x, y, h)
        def event(position, state):
            mass = state[4, :active, :active]
            return np.log(np.linalg.svd(mass, compute_uv=False)[0])+position
        if active and event(x+h, trial) >= 0:
            a, b = 0., h
            for _ in range(36):
                mid = (a+b)/2
                if event(x+mid, rk5(fun, x, y, mid)) >= 0:
                    b = mid
                else:
                    a = mid
            h = (a+b)/2
            trial = rk5(fun, x, y, h)
            thresholds.append(cross(trial, active))
            active -= 1
        x, y = x+h, trial
    assert active == 0
    final = ratio(y)
    terms = y[5].real.flat[:4]
    assert abs(sum(terms)-np.log(final[-1]/initial[-1])) < 4e-9
    return {"high_ratios": list(initial), "low_ratios": list(final),
            "log_terms_down_self_up_on_down_lepton_self_nu_on_lepton": terms.tolist(),
            "thresholds_GeV_descending": thresholds}


def quark_bound(tensors, gauge_upper):
    """Orientation-free upper bound on low-energy d/s, for arbitrary singlet Yukawas."""
    y = model(tensors, False)
    kup, kdown = (np.linalg.norm(y[i], 2)**2 for i in (0, 1))
    k = np.array([kup, kdown])
    b = np.array([41/10, -19/6, -7])
    length = np.log(HIGH/LOW)
    cells = 256
    dl = length/cells
    gauge_path = lambda x: 1/(1/gauge_upper**2-2*b*(length-x)/LOOP)
    integral = np.zeros(2)
    for i in range(cells):
        g2 = np.maximum(gauge_path(i*dl), gauge_path((i+1)*dl))
        force = GAUGE_C[:2] @ g2
        exponent = 2*force*dl/LOOP
        term = 4.5*k/force*np.expm1(exponent)
        integral += LOOP/9*np.log1p(term)
        k = k*np.exp(exponent)/(1+term)
    transfer = float(np.exp(1.5*integral.sum()/LOOP))
    return {"integrals_up_down_upper": integral.tolist(),
            "d_over_s_max_factor": transfer,
            "d_over_s_at_MZ_upper": float(ratio(y)[0]*transfer)}


def get_results():
    archived = load(HERE.parent/"1556/results.json")["candidates"]
    weak = load(HERE.parent/"1386/inputs.json")["yukawa"]
    g = load(HERE.parent/"1369/inputs.json")
    gauge = np.array([g[k][0] for k in ("g1_GUT", "g2", "g3")])
    gauge_hi = np.array([g[name][0]+2*g[name][-1] for name in ("g1_GUT", "g2", "g3")])
    d, e = np.array(weak["down"]), np.array(weak["charged_leptons"])
    ds_lo = (d[0,0]-2*d[0,1])/(d[1,0]+2*d[1,-1])
    double_lo = ds_lo*(e[1,0]-2*e[1,1])/(e[0,0]+2*e[0,-1])
    data = {"input_MZ_lower_d_over_s_and_double_ratio": [ds_lo, double_lo],
            "cases": {}}
    for branch in ("A", "B"):
        t = archived["release_down_"+branch]["tensors"]
        full = transport(t, gauge, True, 160)
        fine = transport(t, gauge, True, 320)
        sm = transport(t, gauge, False, 160)
        for key in full:
            np.testing.assert_allclose(full[key], fine[key], rtol=8e-8, atol=3e-10)
        ratio_lo = ds_lo/fine["low_ratios"][0]
        double_factor = double_lo/fine["low_ratios"][-1]
        bound = quark_bound(t, gauge_hi)
        assert bound["d_over_s_at_MZ_upper"] < ds_lo
        assert ratio_lo > 1 and double_factor > 1
        data["cases"][branch] = {
            "full": fine, "same_UV_without_singlets": sm,
            "quark_only_comparison_arbitrary_singlet_Yukawa": bound,
            "minimum_down_singular_ratio_gain": ratio_lo,
            "minimum_joint_matching_condition_product": double_factor,
            "minimum_common_left_or_right_norm_relative_error_if_one_per_sector":
                (np.sqrt(double_factor)-1)/(np.sqrt(double_factor)+1)}
    return data


def compare(saved, actual):
    if isinstance(saved, dict):
        assert saved.keys() == actual.keys()
        for key in saved:
            compare(saved[key], actual[key])
    else:
        np.testing.assert_allclose(saved, actual, rtol=2e-7, atol=5e-10)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--check", action="store_true")
    group.add_argument("--write-new", action="store_true")
    args = parser.parse_args()
    answer = get_results()
    dest = HERE/"independent_results.json"
    if args.write_new:
        with dest.open("x", encoding="utf-8", newline="\n") as out:
            json.dump(answer, out, indent=2, allow_nan=False)
            out.write("\n")
    else:
        compare(load(dest), answer)
    print("1557 independent signed transport and matching verified")
