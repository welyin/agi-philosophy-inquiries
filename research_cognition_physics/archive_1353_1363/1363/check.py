"""Verify transfer spectra and flavor readouts in the same finite evolution."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


def evolution(h, t):
    e, u = np.linalg.eigh(h)
    return (u * np.exp(-1j * e * t)) @ u.conj().T


def chain(delta, end=0.0, second=1.0):
    return np.array([[0.,1.,0.],[1.,delta,second],[0.,second,end]])


def calculate():
    delta = math.sqrt(8/3)
    time = 4*math.pi/math.sqrt(32/3)
    result = {"time": time, "site": {}, "flavor": {}}
    for name,d in [("resonant",0.),("detuned",delta)]:
        h = chain(d)
        rows = []
        for t in [time/2,time]:
            prop = evolution(h,t)
            prob = abs(prop[:,0])**2
            assert abs(sum(prob)-1) < 1e-12
            rows.append(prob.tolist())
        result["site"][name] = rows
    assert abs(result["site"]["detuned"][1][2]-1) < 1e-12
    nonresonant = chain(.8,end=1.,second=.7)
    p = abs(evolution(nonresonant,.1)[2,0])**2
    assert p > 0
    result["unequal_ends_heralded_rate_at_point_one"] = float(p)
    identity = np.eye(2)
    for angle in (0.,math.pi/6):
        rot = np.array([[math.cos(angle),-math.sin(angle)],[math.sin(angle),math.cos(angle)]])
        mass = rot @ np.diag([2.,3.]) @ rot.T
        full = np.kron(chain(delta),identity)+np.kron(np.eye(3),mass)
        rows = []
        for t in (math.pi/2, math.pi):
            prop = evolution(full,t).reshape(3,2,3,2)
            block = prop[2,:,0,:]
            free = evolution(mass,t)
            scalar = evolution(chain(delta),t)[2,0]
            assert np.allclose(block, scalar*free)
            charge_read = abs(scalar)**2
            conditional = abs(block[1,0])**2/charge_read
            assert abs(conditional - math.sin(2*angle)**2*math.sin(t/2)**2) < 1e-12
            rows.append({"all_flavors_received":float(charge_read),"second_flavor_given_received":float(conditional)})
        result["flavor"][str(round(angle,6))] = rows
    return result


def same(actual, saved):
    if isinstance(actual,dict):
        assert actual.keys()==saved.keys()
        for key in actual:
            same(actual[key],saved[key])
    elif isinstance(actual,list):
        assert len(actual)==len(saved)
        for x,y in zip(actual,saved):
            same(x,y)
    else:
        np.testing.assert_allclose(actual,saved,atol=1e-13,rtol=1e-10)


if __name__ == "__main__":
    cli=argparse.ArgumentParser()
    cli.add_argument("--write-new",action="store_true")
    opts=cli.parse_args()
    data=calculate()
    path=Path(__file__).with_name("results.json")
    if opts.write_new:
        with path.open("x",encoding="utf8") as out:
            json.dump(data,out,indent=2,allow_nan=False)
            out.write("\n")
    else:
        same(data,json.loads(path.read_bytes()))
    print(json.dumps(data))
