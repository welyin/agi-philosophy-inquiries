"""Propagate experimental calibration precision through the charged texture, including nonlinear examples."""
import argparse
import itertools
import json
from pathlib import Path
import runpy

import numpy as np


HERE = Path(__file__).resolve().parent
model = runpy.run_path(str(HERE/"check.py"))
read = model["read"]


def calculate():
    ref = read("results.json")
    cfg = read("inputs.json")
    lattice = model["STRONG"]["QCDSF_s_d"]
    angles = read(cfg["secondary_tree_angle"])["time_integrated_gamma_deg"]
    # For a symmetric local diagnostic use the larger of the asymmetric widths.
    angle = [angles["central"], max(angles["central"]-angles["68"][0], angles["68"][1]-angles["central"])]
    coeff = cfg["mixing_coefficient_ps_inv"]
    clock = model["CLOCKS"]["delta_ms"]
    # Independent local variations in log(s/d), log(gamma), log(C_s), log(DeltaMs).
    D = np.array([[-1.,0,0,0], [0,1.,0,0], [0,0,-1.,1.]])
    sigma = np.array([record[1]/record[0] for record in [lattice,angle,coeff,clock]])
    central_moduli = ref["rows"]["archived_time_integrated_angle"]["moduli"]
    x = np.log([central_moduli[k] for k in ("Vus","Vub","Vcb")])
    initial_jac = np.array(ref["fractional_observable_vs_trial_moduli_jacobian"])
    mean = np.array([item[0] for item in [lattice,angle,coeff,clock]])
    def target(values):
        sd,g,c,dm = values
        return np.log([1/sd,g,dm/c])
    step=.25
    central, x = model["solve"](target(mean),x,initial_jac,step)
    raw_jac = model["jacobian"](x,step,model["logs"](central))
    def low_moduli(u):
        row = model["evaluate"](u,step)
        return np.log([row["moduli"][k] for k in ("Vus","Vub","Vcb")])
    delta = 1e-4
    low_jac = np.column_stack([(low_moduli(x+d)-low_moduli(x-d))/(2*delta)
                               for d in delta*np.eye(3)])
    response = low_jac @ np.linalg.solve(raw_jac,D)
    modes = response * sigma
    covariance = modes@modes.T
    corners = []
    for signs in itertools.product((-1,1),repeat=3):
        # The tiny clock uncertainty is included in the linearization; hold it fixed here.
        shifts = np.array([*signs,0])*sigma
        if signs[1] > 0:
            shifts[1] = (angles["68"][1]-angle[0])/angle[0]
        failed=[]
        for factor in (1.,.5,.25):
            changed = mean*(1+factor*shifts)
            try:
                row, _ = model["solve"](target(changed), x, raw_jac, step)
            except RuntimeError as exc:
                if "shooting did not close" not in str(exc):
                    raise
                failed.append(factor)
                continue
            assert max(abs(model["logs"](row)-target(changed))) < 1e-8
            corners.append(dict(signs=list(signs), factor=factor, failed_trials=failed,
                                moduli=row["moduli"], golden_S=row["sin2beta"]))
            break
        else:
            corners.append(dict(signs=list(signs), failed_trials=failed, factor=None))
    return dict(calibrating_angle=angle[0], fractional_width_per_input=modes.tolist(),
                local_fractional_standard_errors=np.sqrt(np.diag(covariance)).tolist(),
                covariance_log_moduli=covariance.tolist(),sampled_corners=corners)


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    data=json.loads(json.dumps(calculate()),parse_float=lambda x:float(f"{float(x):.8g}"))
    target_file=HERE/"precision_results.json"
    if args.check:
        assert data==json.loads(target_file.read_bytes())
        print("1765: input uncertainty projection and eight independently re-solved corners verified")
    else:
        with target_file.open("x",encoding="utf-8") as output:
            json.dump(data,output,indent=2,allow_nan=False)
            output.write("\n")
