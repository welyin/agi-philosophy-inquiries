"""Two endpoint success rates, mixed-state visibility, and a transfer-phase witness."""
from argparse import ArgumentParser
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def measurements(vector):
    x, y, z = vector
    return np.array([[1+z, x-1j*y], [x+1j*y, 1-z]]) / 2


def quark_matrix(lam, c, a, b):
    edges = np.zeros((3,3), dtype=complex)
    edges[0,1] = b*lam**3
    edges[1,2] = a*lam**2
    return edges + edges.conj().T + np.diag([c*abs(lam)**4,b*abs(lam)**2,1])


def observe(pair):
    vals, frames, kernels = [], [], []
    for m in pair:
        spec, vec = np.linalg.eigh(m)
        order = np.argsort(abs(spec))
        vals.extend(abs(spec[order]))
        frames.append(vec[:,order])
        kernels.append(m @ m.conj().T)
    v = frames[0].conj().T @ frames[1]
    invariant = (v[0,0]*v[1,1]*np.conj(v[0,1]*v[1,0])).imag
    loop = kernels[0][0,2]*kernels[1][2,0]
    return vals, abs(v)**2, float(invariant), float(np.angle(loop))


def calculate():
    cos_half = 1/np.sqrt(2)
    optimal = np.array([0., -cos_half, cos_half])
    ideal = measurements(optimal)
    epsilon = 0.01
    k = cos_half - 2*epsilon
    visibility = k/np.sqrt(1-k*k)
    mixed = measurements([0., -k, k])
    pure = measurements([np.sqrt(1-2*k*k),-k,k])
    assert min(np.linalg.eigvalsh(mixed)) > 0
    assert abs(np.linalg.det(pure)) < 1e-15
    mixed_visibility = abs(mixed[0,1])/np.sqrt((mixed[0,0]*mixed[1,1]).real)
    fidelity_min = 1-epsilon/cos_half
    trace_dist_max = np.sqrt(1-fidelity_min)
    assert np.isclose(visibility,mixed_visibility)
    assert np.isclose(np.trace(mixed@ideal).real,fidelity_min)
    assert np.isclose(sum(abs(np.linalg.eigvalsh(ideal-pure)))/2,trace_dist_max)
    p = json.loads((HERE/"inputs.json").read_text(encoding="utf8"))["HS_rounded_fit_table1"]
    down = p["lambda0"]*np.cos(np.pi/8)
    up = -1j*p["lambda0"]*np.sin(np.pi/8)
    hd = quark_matrix(down,p["c_down"],p["A0"],p["b"])
    direct, twisted = (observe([quark_matrix(s*up,p["c_up"],p["A0"],p["b"]),hd])
                       for s in [1,-1])
    assert np.allclose(direct[0],twisted[0]) and np.allclose(direct[1],twisted[1])
    assert direct[2] > 0 and np.isclose(direct[2],-twisted[2])
    assert np.isclose(direct[3],-np.pi/2) and np.isclose(twisted[3],np.pi/2)
    return {"endpoint_probability_loss":epsilon,"fidelity_lower":float(fidelity_min),
            "trace_distance_upper":float(trace_dist_max), "normalized_visibility_lower":float(visibility),
            "original_and_twisted_J":[direct[2],twisted[2]],"reference_and_population_data_same":True}


if __name__=="__main__":
    parser=ArgumentParser()
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument("--write-new",action="store_true")
    mode.add_argument("--check",action="store_true")
    args=parser.parse_args()
    out=calculate()
    dest=HERE/"independent_results.json"
    if args.write_new:
        with dest.open("x",encoding="utf8") as f:
            json.dump(out,f,indent=2)
            f.write("\n")
    else:
        assert out==json.loads(dest.read_text(encoding="utf8"))
    print(out)
