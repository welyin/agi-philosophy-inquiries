"""Match a common pair of colour poles to each commuting light-family assignment."""
import argparse
import itertools
import json
import math
import runpy
from functools import lru_cache
from pathlib import Path
from types import SimpleNamespace

import numpy as np


HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE / name).read_bytes())


PORT = runpy.run_path(str(HERE / "../1755/joint.py"))
AMP = runpy.run_path(str(HERE / "../1738/independent_check.py"))["amplitudes"]
CLOCK = read("../1736/independent_results.json")["data"]
STRONG = read("../1755/inputs.json")
TABLE = read("../1738/inputs.json")["table_IV"]
C_S = read("../1765/inputs.json")["mixing_coefficient_ps_inv"][0]
TARGET_TBS = CLOCK["delta_ms"][0] / C_S
GOLDEN_S = CLOCK["golden_S"][0]
ORIGINAL = "c1h"


@lru_cache(maxsize=24)
def load_mother(us, ub, cb, step):
    """Use the measured heavy/charged spectrum and new weak trial, with no d/s calibration."""
    module = runpy.run_path(str(HERE / "../1744/check.py"))

    def load(path):
        model = runpy.run_path(str(path))
        direct = model["DIRECT"]
        spec = [("us", us), ("ub", ub), ("cb", cb)]
        weak = {key: (val-direct["V"+key][0])/math.hypot(*direct["V"+key][1:]) for key,val in spec}
        model["CFG"]["nearby_measured_inputs_in_marginal_widths"][module["CFG"]["parent_point"]] = {
            "weak": weak, "g3": 0., "top": 0., "up": 0., "charm": 0.}
        return model

    module["mother"].__globals__["runpy"] = SimpleNamespace(run_path=load)
    bg = module["mother"](step)
    F, z = PORT["matter_filter"](bg)
    return bg, F, z


def split_evolution(bg, F, assignment, gain):
    model, pt = bg["model"], bg["pt"]
    flow = pt["flow"]
    y = bg["yW"].copy()
    y[2] = gain*F.conj().T @ y[2]
    y[3] = F.conj().T @ y[3]
    poles = flow["fields"]["neutral_poles_GeV"]
    transport = np.eye(3, dtype=complex)
    for start, susy, dt, grid in reversed(bg["segments"]):
        assert not susy
        y, transport = model["advance"](y, transport, -dt, grid[::-1], flow)
        if start == bg["Mcolor"]:
            u, ys, rr = np.linalg.svd(y[1, :, :3])
            lu = np.linalg.svd(y[0, :, :3])[0]
            p = abs(lu[:, -1].conj() @ u[:, ::-1])**2
            c = np.sqrt(p[0])
            amps = {"c": c, "1": 1., "h": np.sqrt((1+c)/2)}
            gains = np.array([amps[k] for k in assignment])
            before = ys[::-1]
            after = before * gains
            assert np.all(np.diff(after)>0)
            y[1, :, :3] = (u * after[::-1]) @ rr
            # Two normalized rows on three ordinary and two extra lines furnish these eigenvalues.
            X = np.diag(np.sqrt(1/gains**2-1))[gains!=1]
            W = (1 / np.sqrt(1+np.diag(X@X.T)))[:,None] * np.column_stack((X, np.eye(2)))
            np.testing.assert_allclose(W @ W.T, np.eye(2), atol=1e-14)
            null = np.vstack((np.eye(3), -X))*gains
            np.testing.assert_allclose(W @ null, 0., atol=1e-14)
            np.testing.assert_allclose(null.T @ null, np.eye(3),atol=1e-14)
            color = dict(c=c, threshold_ds=after[0]/after[1], threshold_bs=after[2]/after[1],
                         gain_ds=gains[0]/gains[1], gain_bs=gains[2]/gains[1],
                         parent_zero_residual=(-p[0]*before[0]+p[1]*before[1]-p[2]*before[2])/before[0])
        for j, mass in enumerate(poles):
            if mass == start:
                y[3, :, j] = 0.
    q, lep, gamma, _, _ = flow["extraction"](y)
    mod = dict(zip(("Vus","Vub","Vcb"), np.sin(q[3:])*[np.cos(q[4]),1,np.cos(q[4])]))
    weak = {k[1:]:(v-model["DIRECT"][k][0])/math.hypot(*model["DIRECT"][k][1:]) for k,v in mod.items()}
    matching = model["MATCH"]
    factors, _ = matching["thresholds"](dict(moduli_MZ=mod, weak_offsets=weak, gamma=gamma),
                                        matching["source_thresholds"]())
    _, _, beta, norms = AMP(mod, gamma)
    return dict(moduli=mod, gamma=gamma, S=np.sin(2*beta), tbs_squared=norms[1]**2,
                td_ts=norms[0]/norms[1], ds=q[0]/q[1]/factors[0], bs=q[2]/q[1]/factors[2],
                lepton_spectrum=lep, **color)


def eval_point(log_moduli, assignment, step):
    bg, F, z = load_mother(*np.exp(log_moduli), step)
    gain = 1.
    target = bg["model"]["JOINT"]["LEPTONS"][:,0]
    for _ in range(3):
        row = split_evolution(bg, F, assignment, gain)
        gain *= target[2]/row["lepton_spectrum"][2]
    assert abs(row["lepton_spectrum"][2]/target[2]-1)<2e-11
    row["relative_lepton_masses"] = (row.pop("lepton_spectrum")/target).tolist()
    return dict(**row, tau_gain=gain, z=z)


def values(row):
    return np.log([row["gamma"],row["S"],row["tbs_squared"]])


def solve(target, assignment, start, jac, step):
    x = start.copy()
    for _ in range(10):
        row = eval_point(x, assignment, step)
        residual = values(row)-target
        if max(abs(residual))<3e-10:
            row["max_log_residual"] = max(abs(residual))
            return row, x
        x -= np.linalg.solve(jac, residual)
    raise ValueError((assignment, residual))


def lattice_distances(ds, bs):
    centers = np.array(TABLE["means"][:3])
    widths = np.array(TABLE["sigmas"][:3])
    cov = np.outer(widths,widths)*np.array(TABLE["correlation"])[:3,:3]
    contrast = np.array([0.,-ds,1.])
    ratios = np.array([[0.,-ds,1.],[1.,-bs,0.]])
    off = ratios @ centers
    return dict(MILC_ds=(ds*centers[1]-centers[2])/np.sqrt(contrast@cov@contrast),
                MILC_joint_ds_bs=np.sqrt(off@np.linalg.solve(ratios@cov@ratios.T,off)),
                QCDSF_sd=(1/ds-STRONG["QCDSF_s_d"][0])/STRONG["QCDSF_s_d"][1])


def clean(x):
    if isinstance(x, dict):
        return {k:clean(v) for k,v in x.items()}
    if isinstance(x,(list,tuple,np.ndarray)):
        return [clean(v) for v in x]
    return float(f"{x:.10g}")


def calculate():
    original = read("../1765/results.json")["rows"]["archived_time_integrated_angle"]
    start = np.log([original["moduli"][k] for k in ("Vus","Vub","Vcb")])
    step = .25
    base = eval_point(start, ORIGINAL, step)
    assert abs(base["gamma"]-original["gamma"]) < 2.e-4
    jac = np.column_stack([(values(eval_point(start+d, ORIGINAL, step))-values(base))/2e-4
                           for d in np.eye(3)*2e-4])
    target = np.log([read("../1707/inputs.json")["time_integrated_gamma_deg"]["central"],GOLDEN_S,TARGET_TBS])
    anchor, x = solve(target, ORIGINAL, start, jac, step)
    same_trial, refits = {}, {}
    for perm in itertools.permutations(ORIGINAL):
        name = "".join(perm)
        row = eval_point(x, name, step)
        same_trial[name] = dict(ds=row["ds"], bs=row["bs"], gamma_difference_deg=row["gamma"]-anchor["gamma"],
                                S_difference=row["S"]-anchor["S"], relative_tbs_squared=row["tbs_squared"]/anchor["tbs_squared"]-1)
        got, _ = solve(target,name,x,jac,step)
        got["strong_comparison"] = lattice_distances(got["ds"],got["bs"])
        refits[name] = got
    # Keep both B measurements at their centres. This angle is a retrospective
    # feasibility check inside the updated time-integrated interval.
    alt_target = target.copy()
    alt_target[0] = np.log(64.)
    _, alt_x = solve(alt_target, ORIGINAL, x, jac, step)
    alternate = {}
    for perm in itertools.permutations(ORIGINAL):
        name = "".join(perm)
        row, _ = solve(alt_target, name, alt_x, jac, step)
        alternate[name] = {k:row[k] for k in ("gamma","S","tbs_squared","td_ts","ds","bs","moduli")}
        alternate[name]["strong_comparison"] = lattice_distances(row["ds"], row["bs"])
    assert refits["1ch"]["ds"] > 1.05 * refits[ORIGINAL]["ds"]
    assert all(row["strong_comparison"]["MILC_ds"]>4.48 for row in refits.values())
    assert all(abs(row["strong_comparison"]["QCDSF_sd"])<2 for row in refits.values())
    assert alternate[ORIGINAL]["strong_comparison"]["MILC_joint_ds_bs"]<1
    assert min(row["strong_comparison"]["MILC_joint_ds_bs"] for name,row in alternate.items() if name!=ORIGINAL)>6
    refined = eval_point(x, ORIGINAL, step/2)
    fine = {k:refined[k]-anchor[k] for k in ("ds","bs","gamma","S","tbs_squared")}
    assert max(abs(v) for v in fine.values()) < 4e-7
    return clean(dict(same_trial=same_trial, refits=refits, exploratory_64_deg=alternate, original_finer_step_change=fine,
                      primary_local_jacobian_at_start=jac))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--check",action="store_true")
    args = ap.parse_args()
    results = calculate()
    if args.check:
        assert results == read("independent_results.json")
        print("1766 six colour assignments, charged recalibration, covariance and two-threshold flow verified")
    else:
        (HERE/"independent_results.json").write_text(json.dumps(results,indent=2)+"\n",encoding="utf-8")
