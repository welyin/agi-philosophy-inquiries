"""Finite work store: check the relative clock, both switching works and full states."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

from check import compare, unitary


def run(pulsed, time, stop_store):
    labels = [(n,m,b) for n in range(2) for m in range(3) for b in range(3)]
    pos = {row: j for j,row in enumerate(labels)}
    ef = np.array([n*n+int(m != 0) for n,m,b in labels])
    eb = np.array([b for n,m,b in labels])
    idx = lambda n,m: pos[n,m,2-n*n-int(m != 0)]
    select = [idx(n,m) for n in range(2) for m in range(3)]
    eye = np.eye(len(labels), dtype=complex)
    magnetic = np.zeros_like(eye)
    swap = eye.copy()
    for m in range(3):
        i, j = idx(0, m), idx(0, (m+1)%3)
        magnetic[i,j] = magnetic[j,i] = 1.
        i,j = idx(0,m),idx(1,(m-1)%3)
        swap[i,i] = swap[j,j] = 0
        swap[i,j] = swap[j,i] = 1
    bare = np.diag(ef+eb)
    assert np.allclose(bare@magnetic, magnetic@bare)
    assert np.allclose(bare@swap, swap@bare)
    psi = eye[:,idx(0,0)]
    prepared = unitary(magnetic, math.pi/3) @ psi if pulsed else psi
    ref, path = prepared/math.sqrt(2), swap@prepared/math.sqrt(2)
    frozen_energy = float(np.sum(eb*(abs(ref)**2+abs(path)**2)))
    free_energy = ef if stop_store else ef+eb
    phase = np.exp(-1j*time*free_energy)
    a,b = phase*ref, swap.conj().T@(phase*path)
    bins = [(a+b)/math.sqrt(2), (a-b)/math.sqrt(2)]
    weights = [float(np.vdot(v,v).real) for v in bins]
    field_in_both = float(sum(np.sum(ef*abs(v)**2) for v in bins))
    store_in_both = float(sum(np.sum(eb*abs(v)**2) for v in bins))
    assert np.isclose(store_in_both+field_in_both, 2.)
    assert all(np.isclose(np.sum(abs(v[select])**2),w) for v,w in zip(bins,weights))
    if stop_store:
        expected = ((1/3,4/9) if pulsed else (1.,0.))[int(time>2)]
    else:
        expected = 0.
    assert np.isclose(weights[1], expected)
    joint = []
    # All nonzero output basis states are |0,m,b=2-I_(m!=0)>.
    for v in bins:
        x=v[[idx(0,m) for m in range(3)]]
        rho=np.outer(x,x.conj())
        joint.append({"real":rho.real.tolist(),"imag":rho.imag.tolist()})
    return {"plus_minus":weights,"joint_field_store_densities":joint,
            "freeze_restore_external_work": [-frozen_energy, frozen_energy] if stop_store else [0.,0.],
            "final_field_energy":field_in_both,"final_store_energy":store_in_both}


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--write",action="store_true")
    args=parser.parse_args()
    result={}
    for prep in (False,True):
        for stopped in (False,True):
            result[f"pulse_{prep}_store_clock_frozen_{stopped}"]=[
                run(prep,t,stopped) for t in [math.pi/2,math.pi]]
    path=Path(__file__).with_name("battery_results.json")
    if args.write:
        with path.open("x",encoding="utf8") as out:
            json.dump(result,out,indent=2,allow_nan=False)
            out.write("\n")
    else:
        compare(result,json.loads(path.read_text(encoding="utf8")))
    print({k:[r["plus_minus"][1] for r in rows] for k,rows in result.items()})
