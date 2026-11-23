"""Resolve the two physical charged flavour lines with controllable Yukawa couplings."""
import argparse
import json
import math
from pathlib import Path
import runpy

import numpy as np

HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE / name).read_bytes())


cfg = read("inputs.json")


def upper_matrices():
    model = runpy.run_path(str(HERE / cfg["endpoint_program"]))
    ref = read(cfg["endpoint_data"])["same_1708_point"][cfg["nominal_case"]]
    fl, sim, state = model["fl"], model["sim"], model["state"]
    f = model["charged"]
    fl["cfg"]["direct_semileptonic_moduli"]["Vub"][0] = f["Vub"][0]-np.linalg.norm(f["Vub"][1:])
    state["D0"] = model["new_dirac"](ref["Dirac_ratio_squared"])
    rep = sim["graph"](1, -1, "unit_third", .20)[1]
    y = np.zeros((4,3,4), complex)
    y[0,:,:3] = np.diag(fl["central"]("up"))
    down = fl["central"]("down")
    down[0] -= fl["low"]["yukawa"]["down"][0][1]
    y[1,:,:3] = fl["low_ckm"](rep["gamma_MZ_deg"], fl["angles"]("direct")) @ np.diag(down)
    y[2,:,:3] = sim["zer"]["unpack"](rep["charged_left_MZ"]) @ np.diag(fl["central"]("charged_leptons"))
    poles = model["spectra"]["neutral_poles_GeV"]
    for start, susy, dt, gauge in fl["segments"](.2):
        if start >= poles[2]:
            return y[:2, :, :3].copy(), ref
        assert not susy
        for j, pole in enumerate(poles):
            if start == pole:
                y[3,:,j] = fl["physics"]["Y0"][:,j]
        for i in range(0, len(gauge)-1, 2):
            beta = fl["beta"]
            k1 = beta(y,gauge[i],False)
            k2 = beta(y+dt*k1/2,gauge[i+1],False)
            k3 = beta(y+dt*k2/2,gauge[i+1],False)
            k4 = beta(y+dt*k3,gauge[i+2],False)
            y += dt*(k1+2*k2+2*k3+k4)/6
    raise ValueError("upper cluster missing")


def majority(p, n):
    return sum(math.comb(n,k)*p**k*(1-p)**(n-k) for k in range(n//2+1,n+1))


def filter_rows(y, kind, tol):
    bright, values, _ = np.linalg.svd(y)
    v = values[::-1]
    bright = bright[:,::-1]
    ea, eb = np.pi/(2*v[1]), np.pi/(2*v[2])
    # The first pulse checks identical chirality (up) or changed chirality (down).
    # The second checks identical chirality, regardless of which side it started on.
    if kind == "up":
        per = np.cos(ea*v)**2*np.cos(eb*v)**2
        wrong = max(np.sin(np.pi*tol/2)**2, (ea*ea+eb*eb)*((1+tol)*v[0])**2)
        target = 0
    else:
        per = np.sin(ea*v)**2*np.cos(eb*v)**2
        wrong = max(np.sin(np.pi*tol/2)**2+(eb*(1+tol)*v[1])**2, (ea*(1+tol)*v[0])**2)
        target = 1
    return dict(y=v, basis=bright, target=target, gap=v[1]-v[0], gap2=v[1]**2-v[0]**2,
                per=per, wrong=wrong, pulse_sum=ea+eb)


def chirality_direct(y, out):
    block = np.block([[0*y,y],[y.conj().T,0*y]])
    d,w = np.linalg.eigh(block)
    t,u = np.pi/(2*out["y"][1]),np.pi/(2*out["y"][2])
    left, right = np.diag([1,1,1,0,0,0]), np.diag([0,0,0,1,1,1])
    evol = lambda t: (w*np.exp(-1j*t*d))@w.conj().T
    middle = right if out["target"] == 1 else left
    success = middle @ evol(u) @ middle @ evol(t) @ left
    actual = (success.conj().T@success)[:3,:3]
    expected = (out["basis"]*out["per"])@out["basis"].conj().T
    return np.linalg.norm(actual-expected,2)


def main():
    matrices, ref = upper_matrices()
    tol = cfg["relative_spectral_calibration_corridor"]
    data = [filter_rows(y, k, tol) for y,k in zip(matrices, ("up","down"))]
    charged = read("../1702/inputs.json")["direct_semileptonic_moduli"]
    width = np.linalg.norm(charged["Vus"][1:])
    limit = cfg["amplitude_error_fraction_of_Kl3_width"]*width
    levels = {n:[majority(d["wrong"],n) for d in data] for n in range(1,12,2)}
    n = next(n for n,err in levels.items() if sum(err)<limit)
    assert n==7
    filters, projectors, rows = [], [], {}
    for kind, y, d, eps in zip(("up","down"), matrices, data, levels[n]):
        v, U, idx = d["y"], d["basis"], d["target"]
        ideal = np.outer(U[:,idx],U[:,idx].conj())
        op = (U*majority(d["per"],n)) @ U.conj().T
        filters.append(op)
        projectors.append(ideal)
        residual = chirality_direct(y,d)
        assert residual<1e-10
        rows[kind] = {"Yukawa_at_7171GeV":v.tolist(), "nearest_linear_gap":d["gap"],
                      "nearest_squared_gap":d["gap2"], "unamplified_error_bound":d["wrong"],
                      "uniform_amplitude_error":eps,"nominal_operator_error":np.linalg.norm(op-ideal,2),
                      "direct_chiral_unitary_test_error":residual,
                      "necessary_YYdagger_area":(1-2*eps)/d["gap2"],
                      "necessary_chiral_Y_area":(1-2*eps)/d["gap"],
                      "implemented_chiral_Y_area":2*n*d["pulse_sum"]}
    A,B = filters
    pu,ps = projectors
    us = np.linalg.norm(pu@ps,2)
    np.testing.assert_allclose(us,ref["Vus_at_upper_neutral_cluster"],rtol=1e-9)
    p_obs=(np.trace(B@A@A@B)/np.trace(B@B)).real
    eu,es = levels[n]
    bound = 2*eu+2*es**2/(1-es)**2
    amp_error = eu+es
    interval_r = [2*(us+sign*amp_error)**2 for sign in (-1,1)]
    def tree_ratio(r):
        trace_over_product = (4-16*r+121*r*r)/(6*r)
        return (trace_over_product+math.sqrt(trace_over_product**2-4))/2
    assert interval_r[1] < 2/11
    # Normalized strange-state preparation needs two failed dimensions in its tail.
    prepared = np.trace(B@B).real/3
    return {"underlying_threshold_Vus":us,"K_width":width,"target_amplitude_error":limit,
            "majority_error_by_repetitions":levels,"selected_odd_repetitions":n,
            "ports": rows,"uniform_product_amplitude_error":sum(levels[n]),
            "uniform_sequential_probability_error":bound,"nominal_mixed_port_preparation_probability":prepared,
            "allowed_1709_ratio_if_only_modulus_transfers":interval_r,
            "allowed_tree_mass_ratio_at_same_inputs":[tree_ratio(r) for r in reversed(interval_r)],
            "operator_norm_actual_open_port":np.linalg.norm(A@B,2),
            "nominal_conditional_Cabibbo_probability":p_obs,
            "phase_blind_instrument_squared_amplitude":us*us,
            "unconditional_strange_preparation_floor":(1-es)**2/3}


def tidy(x):
    if isinstance(x,dict): return {str(k):tidy(v) for k,v in x.items()}
    if isinstance(x,list): return [tidy(v) for v in x]
    if isinstance(x,(float,np.floating)): return float(f"{x:.10g}")
    return x


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    out = tidy(main())
    if p.parse_args().check:
        assert read("results.json")==out
        print("1711 verified")
    else:
        (HERE/"results.json").write_text(json.dumps(out,indent=2)+"\n",encoding="utf8")
