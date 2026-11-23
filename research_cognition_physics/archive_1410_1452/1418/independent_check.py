"""Independent inverse amplitudes for unresolved CsI plus Ar and a coherent readout."""
import argparse
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def calculate():
    # Ideal tagged-flavor forward-bin response. Mass numbers approximate nuclear masses.
    argon = np.array([18., 22.])
    crystals = np.array([[55., 78.], [53., 74.]])
    weights = np.array([133., 127.])/260
    matrix = crystals.T @ (weights[:, None]*crystals)
    base = np.array([.5-2*.2386, -.5])
    null = np.array([argon[1], -argon[0]])
    partner = base - 2*(null@matrix@base)/(null@matrix@null)*null
    assert np.allclose(argon@base, argon@partner)
    assert np.isclose(base@matrix@base, partner@matrix@partner)
    silicon = np.array([14., 14.])
    ratio = ((silicon@partner)/(silicon@base))**2
    assert ratio > 15 and not np.allclose(base, partner)
    # Inclusive recoil traces over a surviving coherent neutrino, which a later CC may read.
    coherent = np.array([1., 1j])/np.sqrt(2)
    readout = coherent.copy()
    no_flip = np.eye(2)
    flip_e = np.diag([-1., 1.])
    recoils = [float(np.linalg.norm(q@coherent)**2) for q in (no_flip, flip_e)]
    joint = [float(abs(np.vdot(readout, q@coherent))**2) for q in (no_flip, flip_e)]
    assert np.allclose(recoils, [1.,1.]) and np.allclose(joint, [1.,0.])
    return {"round":1418, "illustrative_tree_sin2thetaW":.2386, "original_proton_neutron":base.tolist(),
            "second_proton_neutron":partner.tolist(), "argon_signed_charge":float(argon@base),
            "unresolved_CsI_response":float(base@matrix@base),
            "isoscalar_Si28_rate_ratio":float(ratio), "medium_forward_ratio":float(partner.sum()/base.sum()),
            "ideal_recoil_only":recoils, "ideal_recoil_then_coherent_projection":joint}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    opt=parser.parse_args()
    result=calculate()
    if opt.check:
        saved=json.loads((HERE/"independent_results.json").read_text(encoding="utf-8"))
        assert result.keys()==saved.keys() and all(np.allclose(result[k],saved[k]) for k in result)
        print("1418 independent target inversion and conditional final-flavor readout verified")
    else:
        print(json.dumps(result,indent=2))
