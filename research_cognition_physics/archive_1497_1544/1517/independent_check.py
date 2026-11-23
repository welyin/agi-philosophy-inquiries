"""Independent partial traces, exact four-state spectra and signed spectral gap."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

BASE = Path(__file__).resolve().parent


def tensor_test(em, eh):
    a = np.zeros((3, 4, 3))
    for k in range(3):
        a[k, k, k] = math.sqrt(1-em)
        a[(k+1) % 3, 3, k] = math.sqrt(em+eh)
    # Rotate the source coordinates so the null direction is not a standard axis.
    frame = np.array([[1,1,1,1], [1,-1,1,-1], [1,1,-1,-1], [1,-1,-1,1]]) / 2
    a = np.einsum("ijk,jl->ilk", a, frame)
    ell = frame[3]
    l = math.sqrt(1-3*(em+eh))*ell
    mx = np.vstack([a[:, :, k] for k in range(3)])
    hx = a.reshape((12,3))
    av = mx.T @ mx
    dm = av + np.outer(l,l) - np.eye(4)
    dh = hx.T @ hx - np.eye(3)
    spent = np.linalg.norm(mx @ ell)**2
    np.testing.assert_allclose([np.linalg.norm(dm,2),np.linalg.norm(dh,2),spent],
                               [em,eh,3*(em+eh)], atol=3e-15)
    return a, ell, {"epsilon_M":em,"epsilon_H":eh,"projected_cost":round(spent,12)}


def full_grams(y, m, ell, phi):
    matrices, dets = [], []
    for c in (2,-3):
        line = m + c*np.exp(1j*phi)*ell
        mass = np.vstack([y,line])
        matrices.append(mass @ mass.conj().T)
        dets.append(np.linalg.det(mass))
    return matrices, dets


def exact_and_perturbed():
    a, ell, exact = tensor_test(0.,0.)
    m = np.array([1.25,-.4,.8,.61])
    weights = np.array([2.,3.,6.])/7
    phi = 0.61
    y = a @ weights
    (d,e), dets = full_grams(y,m,ell,phi)
    assert np.max(abs(d.imag)) < 2e-15 and np.max(abs(e.imag)) < 2e-15
    assert abs(dets[0].imag)>1e-3
    ds,es=np.linalg.eigvalsh(d),np.linalg.eigvalsh(e)
    assert np.all(ds<=es+1e-12) or np.all(es<=ds+1e-12)
    assert np.linalg.matrix_rank(d-e,tol=1e-12)==1
    random = np.random.default_rng(130519)
    worst=0.
    for _ in range(180):
        y = random.normal(size=(3,4))+1j*random.normal(size=(3,4))
        m = random.normal(size=4)+1j*random.normal(size=4)
        ell_t = random.normal(size=4)+1j*random.normal(size=4)
        (d,e),_ = full_grams(y,m,ell_t,phi)
        off = 5*np.linalg.norm(y@ell_t.conj())
        eig = np.linalg.eigvalsh(d-e)
        np.testing.assert_allclose(eig[[1,2]],0,atol=5e-13)
        np.testing.assert_allclose(-eig[0]*eig[-1],off**2,rtol=2e-13,atol=1e-12)
        gaps=(np.linalg.eigvalsh(d)-np.linalg.eigvalsh(e))[:3]
        ratio=max(0.,max(gaps))*max(0.,-min(gaps))/off**2
        assert ratio<=1.+2e-12
        worst=max(worst,ratio)
    return {"exact_grams":exact,"exact_full_d_squared_spectrum":ds.tolist(),
            "exact_full_e_squared_spectrum":es.tolist(),
            "saturated_down_det_angle_mod_pi":float(np.angle(dets[0]) % math.pi),
            "finite_complex_checks":180,"max_observed_ratio":float(worst)}


def audit():
    tab=np.array(json.loads((BASE/"../1402/inputs.json").read_text(encoding="utf8"))["scales"]["1e16"])
    boxes={}
    for width in (0,1,2):
        low,high=tab[:,0]-width*tab[:,1],tab[:,0]+width*tab[:,1]
        plus=max(0.,max(low[:3]**2-high[3:]**2))
        minus=max(0.,max(low[3:]**2-high[:3]**2))
        assert plus>0 and minus>0
        boxes[str(width)]={"P_plus":float(plus),"P_minus":float(minus),
                           "min_norm_Y_ell_dimensionless":math.sqrt(plus*minus)/5}
    ads=json.loads((BASE/"../1380/inputs.json").read_text(encoding="utf8"))["lhcb"]
    assert abs(ads["asymmetry"])>2*(ads["stat"]+ads["syst"])
    return {"same_scale_yukawa_boxes":boxes,"complete_finite_mass":exact_and_perturbed(),
            "unequal_residual_sharp_check":tensor_test(.013,.027)[2],
            "sole_left_weak_CP_source_ADS":{"prediction":0.,"observed":ads["asymmetry"],
            "stat":ads["stat"],"syst":ads["syst"]}}

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    m = p.add_mutually_exclusive_group()
    m.add_argument("--write-new", action="store_true")
    m.add_argument("--check", action="store_true")
    args = p.parse_args()
    path = BASE/"independent_results.json"
    out = json.loads(json.dumps(audit(), default=float))
    if args.write_new:
        with path.open("x", encoding="utf-8") as f:
            f.write(json.dumps(out, indent=2)+"\n")
    elif args.check:
        assert out == json.loads(path.read_text(encoding="utf-8"))
        print("Independent tests passed.")
    else:
        print(json.dumps(out, indent=2))
