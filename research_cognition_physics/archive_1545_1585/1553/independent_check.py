"""Doubly-charged tau rate, weighted right-current row, and exact six pole traces."""
import argparse
import itertools
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
DATA = HERE.parent/"1402"/"inputs.json"


def cone_bound(d, e):
    # The right-hand side is the squared sine-complement of the enclosing ball.
    down, strange, bottom = d
    electron, muon, tau = e
    radius2 = (bottom+muon)**2+(strange+electron)**2+down**2
    return max(0.,1-radius2/tau**2)


def spectrum(f, d):
    zero = np.zeros((3,3),complex)
    mass = np.block([[zero,d],[d.T,f]])
    _, vals, vh = np.linalg.svd(mass)
    right = vh.conj().T[3:,:]
    weights = abs(right)**2
    strength = np.vdot(f,f).real
    tr_right = np.sum(weights*vals**2)
    np.testing.assert_allclose(2*tr_right-np.sum(vals**2),strength,rtol=1e-10)
    mean_tau = weights[2]@vals
    assert mean_tau+1e-12 >= abs(f[2,2])
    return vals, weights, float(strength)


def random_unitary(rng):
    z = rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
    return np.linalg.qr(z)[0]


def evaluate():
    raw = json.loads(DATA.read_text(encoding="utf-8"))
    table = np.array(raw["scales"]["1e16"])
    d, e = table[:,0].reshape(2,3)
    corners = np.array([table[:,0]+table[:,1]*sgn for sgn in
                        itertools.product((-2,2), repeat=6)])
    beta = min(cone_bound(row[:3],row[3:]) for row in corners)
    # Check the theorem also on genuinely complex unified flavor axes.
    rng = np.random.default_rng(91553)
    for _ in range(60):
        v = random_unitary(rng)
        f = (v*d)@v.T-np.diag(e)
        vertex_frac = abs(f[2,2])**2/np.linalg.norm(f)**2
        assert vertex_frac >= cone_bound(d,e)-1e-13
        # Independent Dirac couplings are allowed in the spectral identity.
        neutral_d = 0.001*random_unitary(rng)@np.diag([.1,.4,.7])
        spectrum(f,neutral_d)
    # Independent hierarchy and completeness budgets.
    eps = 0.001
    eta = 0.01
    ratio = 0.2
    norm2 = 1-2*eta-3*eps**2
    pmin = (np.sqrt(beta*norm2)-ratio)/(1-ratio)
    f = -np.diag(d+e)
    m, weights, _ = spectrum(f,np.diag([1e-7,2e-7,3e-7]))
    assert m[1]/m[0]<ratio and m[3]/m[0]<eps
    assert min(weights[:,:3].sum(axis=0))>1-eta
    assert weights[2,0]>pmin
    # Keep charged spectra but independently choose a hierarchical R with the
    # largest mass on e. All masses below are merely a well-conditioned example.
    independent = np.diag([0.02,0.002,0.001])
    ms, ws, _ = spectrum(independent,np.diag([1e-7,2e-7,3e-7]))
    assert ms[1]/ms[0]<ratio and ws[2,0] < 1e-12
    return {"tau_tau_cone_floor_two_width":beta, "neutral_mass_squared_trace_floor":norm2,
            "right_tau_probability_heaviest_floor":float(pmin),
            "shared_hierarchical_example_heaviest_weight":float(weights[2,0]),
            "independent_R_heaviest_weight":float(ws[2,0])}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--check",action="store_true")
    group.add_argument("--write-new",action="store_true")
    args = parser.parse_args()
    result = evaluate()
    target = HERE/"independent_results.json"
    if args.write_new:
        with target.open("x",encoding="utf-8") as stream:
            json.dump(result,stream,indent=2,allow_nan=False)
            stream.write("\n")
    else:
        saved = json.loads(target.read_text(encoding="utf-8"))
        assert result.keys()==saved.keys()
        for key in result:
            np.testing.assert_allclose(saved[key],result[key],atol=1e-12,rtol=1e-10)
        print("1553: charged cone bound, exact neutral trace and heaviest tau-tag tested")
