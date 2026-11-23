"""Inverse neutral-current normalization using physical proton and cesium targets."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np


BASE = Path(__file__).resolve().parent
Z, N = 55, 78
ZPOLE_REF = .23129
THOMSON_REF = .23868


def loaded(path):
    return json.loads((BASE/path).read_text(encoding="utf-8-sig"))


def one_loop(k, threshold=None):
    # Inputs are the reference central non-weak anchors adopted by the author.
    inv_em, inv_s, mz = 127.934, 1/.11823, 91.1876
    sm = np.array([41/6, -19/6, -7.])
    high = sm if threshold is None else np.array([11.,1.,-3.])
    ell = 0 if threshold is None else math.log(threshold/mz)/(2*math.pi)
    subtract = (sm-high)*ell
    unified, t = np.linalg.solve([[k+1, high[0]+high[1]], [1, high[2]]],
                                 [inv_em-sum(subtract[:2]),inv_s-subtract[2]])
    inv_at_z = unified*np.array([k,1,1])+t*high+subtract
    assert t > ell and np.all(inv_at_z > 0)
    assert abs(sum(inv_at_z[:2])-inv_em) < 2e-12 and abs(inv_at_z[2]-inv_s)<2e-12
    return {"s_Z":inv_at_z[1]/inv_em, "MX_GeV":mz*math.exp(2*math.pi*t)}


def limits(center, width, n):
    return [center-n*width, center+n*width]


def rho_and_angle(p,cs):
    rho = (Z*p-cs)/N
    angle = (1-p/rho)/4
    assert rho > 0
    np.testing.assert_allclose([rho*(1-4*angle),rho*(Z*(1-4*angle)-N)],[p,cs])
    return rho, angle


def minimum_transfer(z, p, cs):
    # Minkowski enlargement of the experimental rectangle along each nuisance axis.
    assert 0<z<N/Z
    vals = [(N-Z*z)*x+z*y for x,y in itertools.product(p,cs)]
    gap = max(min(vals),-max(vals),0.)
    eps = gap/(N-Z*z+(Z+N)*z)
    # Check the resulting ray hits both expanded observation intervals.
    p_rho = [(p[0]-eps)/z,(p[1]+eps)/z]
    c_rho = sorted([(cs[0]-(Z+N)*eps)/(Z*z-N),(cs[1]+(Z+N)*eps)/(Z*z-N)])
    assert max(p_rho[0],c_rho[0],0) <= min(p_rho[1],c_rho[1])+1e-12
    return eps


def compute():
    historical = loaded("../1385/inputs.json")["proton"]
    atom = loaded("../1462/inputs.json")
    qp, errp = historical["pv_electron_scattering_only"]
    qpref, referrp = historical["SM"]
    qcref, referrc = atom["baseline_Cs"]
    low_shift = THOMSON_REF-ZPOLE_REF
    p_ref = 1-4*THOMSON_REF
    c_ref = Z*p_ref-N
    pred = {"pure_SM_k53":one_loop(5/3),
            "one_TeV_MSSM_k53":one_loop(5/3,1000.),
            "pure_SM_k43":one_loop(4/3),
            "pure_SM_k32_over25":one_loop(32/25)}
    results = {}
    for name,(qc, e1, e2) in atom["atomic_branches"].items():
        summaries = {}
        for n in (0,1,2):
            p = limits(qp-qpref+p_ref, errp+referrp, n)
            cs = limits(qc-qcref+c_ref, e1+e2+referrc, n)
            calculated = np.array([rho_and_angle(x,y) for x,y in itertools.product(p,cs)])
            rho_range, angle_range = np.array([calculated.min(axis=0),calculated.max(axis=0)]).T
            raw_p = limits(qp,errp,n)
            raw_c = limits(qc,e1+e2,n)
            scenarios = {}
            for key,item in pred.items():
                sz = item["s_Z"]
                z = 1-4*(sz+low_shift)
                scenarios[key] = {
                    "required_low_minus_pole":(angle_range-sz).tolist(),
                    "min_abs_charge_per_nucleon_raw":minimum_transfer(z,raw_p,raw_c),
                    "min_extra_per_nucleon_after_named_matching":minimum_transfer(z,p,cs)
                }
            summaries[str(n)] = {"range_rho":rho_range.tolist(),
                                 "angle_range_low":angle_range.tolist(),"scenarios":scenarios}
        results[name] = summaries
    return {"GUT_central_independent_solution":pred,"published_reference_shift":low_shift,
            "independent_proton_Cs":results}


def equal(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for key in a:
            equal(a[key],b[key])
    else:
        np.testing.assert_allclose(a,b,rtol=1e-10,atol=1e-13)


if __name__ == "__main__":
    cli=argparse.ArgumentParser()
    mode=cli.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write-new",action="store_true")
    mode.add_argument("--check",action="store_true")
    args=cli.parse_args()
    data=compute()
    file=BASE/"independent_results.json"
    if args.write_new:
        payload=json.dumps(data,indent=2,allow_nan=False)+"\n"
        with file.open("x",encoding="utf-8") as writer:
            writer.write(payload)
    else:
        equal(data,loaded(file))
        print("1568 independent gauge, weak-charge and normalization checks passed.")
