"""Independent real-quadrature test for the sine branch of entangled neutral B decays."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


HERE=Path(__file__).resolve().parent


def band(center,errors):
    delta=2*math.sqrt(sum(e*e for e in errors))
    return max(-1.,center-delta),min(1.,center+delta)


def calculate():
    plus=band(1.37/2,[.14/2,.06/2])
    minus=band(1.17/2,[.18/2,.11/2])
    s=(max(plus[0],minus[0]), min(plus[1],minus[1]))
    ds=band(.80,[.14,.06,.03])
    dc=band(.91,[.22,.09,.07])
    phi=(max(math.asin(ds[0]),math.acos(dc[1])),min(math.asin(ds[1]),math.acos(dc[0])))
    near=(math.asin(s[0]),math.asin(s[1]))
    far=(math.pi-near[1],math.pi-near[0])
    gap=far[0]-phi[1]
    assert max(near[0],phi[0]) < min(near[1],phi[1]) and gap>1.
    # Two fictive but physical Dalitz bins with fixed magnitude and strong-phase differences.
    phase=.74
    z=np.array([.72+.31j,.33-.71j])
    def readout(angle,strong):
        return (np.exp(1j*angle)*strong).imag
    np.testing.assert_allclose(readout(phase,z),readout(math.pi-phase,z.conj()))
    gram=np.real(np.outer(z,z.conj()))
    assert np.allclose(gram,np.real(np.outer(z.conj(),z)))
    assert np.linalg.norm(readout(phase,z)-readout(math.pi-phase,z))>.5
    return {"round":1427,"BaBar_common_sine_rectangle":s,"independent_Dalitz_sine_box":ds,
            "independent_Dalitz_cosine_box":dc,"allowed_small_angle_rad":near,
            "second_B_sine_branch_rad":far,"D_filter_oriented_arc_rad":phi,
            "minimum_phase_mismatch_for_second_branch_degrees":math.degrees(gap),
            "unoriented_charm_bins_leave_real_quadrature_ambiguous":True}


def same(a,b):
    return a.keys()==b.keys() and all(np.allclose(a[k],b[k],rtol=1e-10,atol=1e-12) for k in a)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    args=p.parse_args()
    result=calculate()
    if args.check:
        assert same(result,json.loads((HERE/"independent_results.json").read_text(encoding="utf-8")))
        print("1427 independent Dalitz quadrature and strong-phase mirror verified")
    else:
        print(json.dumps(result,indent=2))
