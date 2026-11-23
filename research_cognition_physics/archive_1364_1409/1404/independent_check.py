"""Direct CKM-ratio inversion with the rare-decay measurement held out from calibration."""
import argparse
import importlib.util
import json
from math import hypot, radians, sqrt
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OLD = HERE.parents[1] / "archive_1323_1352" / "1351"
spec = importlib.util.spec_from_file_location("charged_tree_1351", OLD / "check.py")
prior = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prior)


def calculate():
    # Physical neutral f_Bd, not f_B+; held fixed in the central comparison.
    fd, fs = 190.5, 230.7
    md, ms, mu = 5279.58, 5366.90, 105.6584
    beta = lambda mass: sqrt(1-(2*mu/mass)**2)
    factor = (fd/fs)**2 * (md/ms) * (beta(md)/beta(ms)) * 1.517/1.616
    obs = .035
    low, high = obs - hypot(.012,.001), obs + hypot(.014,.002)
    k, w, rho = .2157/.9704, .0351/.93, .083
    preds = {}
    for gamma in (65,245):
        v2, matrix = prior.ratio(k,w,rho,radians(gamma))
        assert np.isclose(abs(matrix[2,0]/matrix[2,1])**2,v2)
        preds[str(gamma)] = factor * v2
        assert low < preds[str(gamma)] < high
    # Present external central scales. The uncertainty calculation conditions on these values.
    v = sqrt(obs/factor)
    return dict(fixed_hadron_clock_factor=factor,measured_R_working_box=[low,high],
                rare_inferred_abs_Vtd_over_Vts=v,ratio_box_fixed_other_inputs=[sqrt(low/factor),sqrt(high/factor)],
                gamma_to_predicted_R=preds)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    want=p.parse_args().check
    result=calculate()
    if want:
        saved=json.loads((HERE/"independent_results.json").read_text(encoding="utf8"))
        assert result.keys()==saved.keys()
        for k,v in result.items():
            if isinstance(v,dict):
                assert v.keys()==saved[k].keys() and all(np.isclose(a,saved[k][b]) for b,a in v.items())
            else:
                assert np.allclose(v,saved[k])
        print("1404 independent rare-decay CKM inversion and both phase-branch witnesses verified")
    else:
        print(json.dumps(result,indent=2))
