"""Finite-Higgs form factor and short-distance Yukawa response of the same prompt pair."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def scalar_triangle(mass, mh):
    tau = (mh/(2*mass))**2
    if tau <= 1:
        f = math.asin(math.sqrt(tau))**2
    else:
        v = math.sqrt(1-1/tau)
        f = -(math.log((1+v)/(1-v))-1j*math.pi)**2/4
    return 1.5*(tau+(tau-1)*f)/tau**2


def examine(masses, s, config, v, lam):
    # In the down sector take b and two heavy axes; d,s are spectators.
    c = math.sqrt(1-s*s)
    n = np.array([1.,1.])/math.sqrt(2)
    rotation = np.eye(3)
    rotation[0,0] = c
    rotation[0,1:] = s*n
    rotation[1:,0] = -s*n
    rotation[1:,1:] += (c-1)*np.outer(n,n)
    assert np.allclose(rotation@rotation.T, np.eye(3))
    projector = rotation.T @ np.diag([1.,0.,0.]) @ rotation
    yb = config["low_yb_matching"]
    d = np.diag([yb,*(math.sqrt(2)*np.array(masses)/v)])
    higgs = projector @ d
    # Trace difference without subtracting two almost coincident numbers.
    avg = sum(m*m for m in masses)/v**2
    shift = s*s*(avg-yb*yb)
    shift4 = shift*(2*yb**2+shift)
    y = higgs@higgs.T
    assert np.isclose(np.trace(y)-yb**2, shift, rtol=1e-7, atol=1e-16)
    assert np.isclose(np.trace(y@y)-yb**4, shift4, rtol=1e-7, atol=1e-19)
    h = config["kinematic_GeV"]["H"]
    b = config["kinematic_GeV"]["b"]
    factor = sum(scalar_triangle(m,h) for m in masses)/2-scalar_triangle(b,h)
    assert np.isclose(np.trace(projector),1.)
    assert np.allclose(projector.diagonal()[1:],s*s/2)
    for relative_v in [.4,1.001,1.7]:
        matrix = np.diag([relative_v,1,1]) @ rotation @ d
        assert np.isclose(np.linalg.det(matrix),relative_v*np.linalg.det(d))
    return {"max_trY2_shift":shift,"max_trY4_shift":shift4,
            "local_beta_lambda_shift_at_reference":(12*lam*shift-6*shift4)/(16*math.pi**2),
            "hgg_amplitude_shift_per_s_squared_real_imag":[factor.real,factor.imag],
            "hgg_heavy_unit_absolute_change":abs(factor)*s*s}


def evaluate():
    config = json.loads((HERE/"inputs.json").read_text(encoding="utf-8"))
    lam = json.loads((HERE/config["rg"]).read_text(encoding="utf-8"))["reference"]["lambda"]
    outputs = {}
    for masses in config["heavy_examples_GeV"]:
        outputs[str(masses[0])] = examine(masses,max(config["mixing_s_window"]),config,
                                          config["v_GeV"],lam)
    assert all(d["max_trY2_shift"]<1.2e-6 and d["hgg_heavy_unit_absolute_change"]<1.2e-8
               for d in outputs.values())
    return {"window_s_upper":max(config["mixing_s_window"]),"spectra":outputs}


if __name__=="__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    args = p.parse_args()
    output = HERE/"independent_results.json"
    result = evaluate()
    if args.check:
        assert result==json.loads(output.read_text(encoding="utf-8"))
    else:
        with output.open("x",encoding="utf-8") as stream:
            stream.write(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print("1610 independent: prompt sharing, traces and finite-mass Higgs loop verified")
