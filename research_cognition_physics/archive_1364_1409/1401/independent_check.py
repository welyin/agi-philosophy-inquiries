"""Absolute core gain certificate for two invertible three-family maps."""
import argparse
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(rel):
    return json.loads((HERE / rel).read_text(encoding="utf-8-sig"))


def mixing(s12, s13, s23, delta):
    c12, c13, c23 = np.sqrt(1 - np.square([s12, s13, s23]))
    z = np.exp(1j * delta)
    r12 = np.array([[c12,s12,0],[-s12,c12,0],[0,0,1]])
    r13 = np.array([[c13,0,s13/z],[0,1,0],[-s13*z,0,c13]])
    r23 = np.array([[1,0,0],[0,c23,s23],[0,-s23,c23]])
    return r23 @ r13 @ r12


def evaluate(up, down, angles, g, radians):
    up, down = np.array(up), np.array(down)
    a, b, c, phase = np.array(angles)[:,0]
    if radians:
        a, b, c = np.sin([a,b,c])
    yu = np.diag(up[:,0])
    yd = mixing(a,b,c,phase) @ np.diag(down[:,0])
    left, s, rhs = np.linalg.svd(np.linalg.solve(yu, yd))
    w, t, vh = np.linalg.svd((yu @ left) * np.sqrt(s))
    profile = ((w * t) @ w.conj().T) / g
    unitary = w @ vh
    cu = g * (unitary / np.sqrt(s)) @ left.conj().T
    cd = g * (unitary * np.sqrt(s)) @ rhs
    primitive = np.hstack([np.linalg.svd(x/g, compute_uv=False) for x in (cu,cd)])
    p = np.sqrt(max(s[0],1/s[-1]))
    residual = max(np.linalg.norm(yu-profile@cu), np.linalg.norm(yd-profile@cd))
    assert np.allclose([min(primitive),max(primitive)],[1/p,p])
    assert 0 < min(t)/g <= max(t)/g < 1 and residual < 1e-12
    ratio = (np.prod(up[:,0])/np.prod(down[:,0]))**(1/3)
    lower = (np.prod(up[:,0]-2*up[:,1])/np.prod(down[:,0]+2*down[:,1]))**(1/3)
    upper = (np.prod(up[:,0]+2*up[:,1])/np.prod(down[:,0]-2*down[:,1]))**(1/3)
    assert 1 < lower < ratio < upper
    return dict(g=g,relative_singular=s.tolist(),sharp_absolute_band_factor=p,
                sharp_log_radius=np.log(p),unit_volume_up_down_gain_ratio=ratio,
                unit_volume_ratio_double_width_box=[lower,upper],
                profile_eigenvalue_range=[min(t)/g,max(t)/g],
                max_reconstruction_error=float(residual))


def evaluate_all():
    low = read("../1386/inputs.json")
    high = read("../1397/inputs.json")["high_scale"]
    return {"MZ":evaluate(low["yukawa"]["up"],low["yukawa"]["down"],
                          low["quark_sines_and_phase"],read("../1369/inputs.json")["g2"][0],False),
            "1e16_GeV_SM_running":evaluate(high["up"],high["down"],
                          high["angles_theta12_theta13_theta23_delta"],high["g2"][0],True)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = evaluate_all()
    if args.check:
        saved = read("independent_results.json")
        assert saved.keys() == result.keys()
        for label, section in result.items():
            assert section.keys() == saved[label].keys()
            for key, value in section.items():
                assert np.allclose(value,saved[label][key],rtol=1e-10,atol=1e-14), (label,key)
        print("1401 fixed absolute gain, contraction reconstruction and volume checks passed")
    else:
        print(json.dumps(result, indent=2, allow_nan=False))
