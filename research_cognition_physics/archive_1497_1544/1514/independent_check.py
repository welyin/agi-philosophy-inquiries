"""Invert the physical spectra of two hyperplane restrictions of one mass map."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def amplitudes(poles, measured):
    eig = np.square(measured)
    weights = [np.prod(pole - eig) / np.prod(pole - np.delete(poles, i))
               for i, pole in enumerate(poles)]
    np.testing.assert_allclose(sum(weights), 1., atol=2e-12)
    assert min(weights) > -1e-13
    return np.sqrt(np.maximum(weights, 0.))


def readout(parent, line, target):
    bare = np.column_stack((np.zeros(3), np.diag(parent)))
    gram = bare @ (np.eye(4) - np.outer(line, line.conj())) @ bare.T
    np.testing.assert_allclose(np.linalg.eigvalsh(gram), np.square(target),
                               rtol=2e-7, atol=1e-17)


def calculation():
    rows = json.loads((HERE/"../1402/inputs.json").read_text(encoding="utf8"))["scales"]["1e16"]
    vals = np.array(rows)
    d, e = np.split(vals[:,0], 2)
    examples = {}
    for scale in (1, 2, 10):
        t = np.maximum(d,e)
        t[2] *= scale
        poles = np.r_[0., t*t]
        a,b = (amplitudes(poles, y) for y in (d,e))
        components = a*b
        hi = sum(components)
        lo = max(0., 2*max(components)-hi)
        readout(t,a,d)
        readout(t,b,e)
        # Both endpoint phase choices are attained in these examples (<=2 shared components).
        signed = b.copy()
        signed[components.argmax()] *= -1
        np.testing.assert_allclose(abs(a @ signed), lo, atol=1e-12)
        readout(t,signed,e)
        examples[str(scale)] = {"Y_frobenius":float(np.linalg.norm(t)),
                                "X_singular_value_range":[float(lo),float(hi)],
                                "least_angle_deg":math.degrees(math.acos(min(1.,hi)))}
    def minimum_parent_overlap(z):
        x,y = np.split(z,2)
        return np.prod(np.minimum(x,y)/np.maximum(x,y))
    corners = [minimum_parent_overlap(vals[:,0]+2*vals[:,1]*sign)
               for sign in itertools.product((-1,1), repeat=6)]
    np.testing.assert_allclose(minimum_parent_overlap(vals[:,0]), examples["1"]["X_singular_value_range"])
    return {"conditions":"clean transfer, common 3x4 Yukawa, ranked hyperplanes at GUT scale",
            "pairings_fixed_throughout_box":True,
            "minimal_frobenius_parent_cosine_two_width":[float(min(corners)),float(max(corners))],
            "fixed_parent_cases":examples}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write-new",action="store_true")
    mode.add_argument("--check",action="store_true")
    args = parser.parse_args()
    out = HERE/"independent_results.json"
    result = calculation()
    if args.write_new:
        with out.open("x",encoding="utf8") as f:
            json.dump(result,f,indent=2)
            f.write("\n")
    else:
        stored = json.loads(out.read_text(encoding="utf8"))
        for k in result["fixed_parent_cases"]:
            for field,v in result["fixed_parent_cases"][k].items():
                np.testing.assert_allclose(v,stored["fixed_parent_cases"][k][field],rtol=1e-11)
        np.testing.assert_allclose(result["minimal_frobenius_parent_cosine_two_width"],
                                   stored["minimal_frobenius_parent_cosine_two_width"],rtol=1e-11)
    print("Verified compressed mass spectra and attainable shared gauge overlaps.")
