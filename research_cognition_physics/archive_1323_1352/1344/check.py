"""Shared-record Ward identity, finite cones, and the soft/detector transfer size."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


def color_gram(lam):
    return np.array([[1., lam/2-1, -lam/2], [lam/2-1, 1., -lam/2],
                     [-lam/2, -lam/2, lam]])


def cone_kernel(a, beta):
    return np.log1p(a*np.tan(beta)**2/2)/np.log1p(np.tan(beta)**2)


def integrate(a, beta):
    relative_angle = np.arccos(1-a)
    nodes, weights = np.polynomial.legendre.leggauss(24)
    x = (1+np.cos(beta) + (1-np.cos(beta))*nodes[:, None])/2
    sine = np.sqrt(1-x*x)
    phi = np.arange(96)*2*np.pi/96
    w = a / ((1-sine*np.cos(phi))*(1-sine*np.cos(phi-relative_angle)))
    return weights @ np.mean(w, axis=1) * (1-np.cos(beta))/2


def ratio(a, lam, beta):
    # Three angles are opposite jets in descending energy order.
    # For gluon i the opposite pair supplies the q-qbar opening angle.
    x = 2*np.sqrt(np.maximum(0, a*(2-a)))
    x /= np.sum(x)/2
    raw = np.array([(1-x[i])*(x[(i+1)%3]**2+x[(i+2)%3]**2) for i in range(3)])
    probs = raw / np.sum(raw)
    f = cone_kernel(a, beta)
    per_tag = lam/2*np.sum(f) + (1-lam)*f
    return float(probs @ per_tag)


def run():
    beta = math.pi/6
    cap = 1.5*cone_kernel(1.5, beta)
    rng = np.random.default_rng(1344)
    for lam in (1., 2.25, .31, 3.7):
        g = color_gram(lam)
        assert np.linalg.norm(g @ np.ones(3)) < 1e-14
        eig = min(2-lam/2,1.5*lam)
        for _ in range(10):
            tags = rng.normal(size=(3,4))+1j*rng.normal(size=(3,4))
            defect = np.trace(tags.conj().T @ g @ tags).real
            assert defect+1e-12 >= eig*np.linalg.norm(tags-tags.mean(axis=0))**2
    for a in (.01, .6, 1., 1.7, 2.):
        assert abs(cone_kernel(a,beta)-integrate(a,beta)/integrate(2.,beta)) < 1e-12
    samples = []
    for t2,t3,n3,stat in [(113,135,.882,.027),(121,135,.856,.023),(131,148,.759,.024)]:
        angles = np.radians([360-t2-t3,t2,t3])
        a = 1-np.cos(angles)
        records = dict(angles_degrees=[360-t2-t3,t2,t3], corrected_charged_n3=n3,
                       stat_only_n3=stat, measured_ratio=n3/.570,
                       su3_central_soft=ratio(a,2.25,beta),so3_central_soft=ratio(a,1.,beta),
                       transfer_needed_vs_global_so3=n3/.570/cap-1)
        assert records["so3_central_soft"] <= cap
        samples.append(records)
    testrange = np.linspace(.02, .48, 75)*np.pi
    for i in testrange:
        for j in testrange:
            # Exterior angles are positive and add up to 2*pi for interior directions
            angles = np.pi-np.array([i,j,np.pi-i-j])
            a = 1-np.cos(angles)
            assert np.sum(cone_kernel(a,beta))/2 <= cap+1e-13
    return dict(round=1344, cap_point=9/8, cap_cone30=float(cap),angular_examples=samples,
                tagged_ward_nonzero_eigenvalues_in_CF={"su3":[.875,3.375],"so3":[1.5,1.5]})


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", action="store_true")
    result = run()
    path = Path(__file__).with_name("results.json")
    if ap.parse_args().write:
        with path.open("x",encoding="utf8") as f:
            json.dump(result,f,indent=2)
            f.write("\n")
    else:
        saved = json.loads(path.read_text(encoding="utf8"))
        assert np.isclose(result["cap_cone30"], saved["cap_cone30"])
        for row, ref in zip(result["angular_examples"],saved["angular_examples"], strict=True):
            assert row.keys() == ref.keys()
            for k in row:
                np.testing.assert_allclose(row[k],ref[k],rtol=1e-12,atol=1e-12)
        assert all(saved[k]==result[k] for k in ("round","cap_point","tagged_ward_nonzero_eigenvalues_in_CF"))
    print(json.dumps(result))
