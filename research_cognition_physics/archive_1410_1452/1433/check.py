"""Spectral and orientation tests for a single rank-one quark-lepton splitting."""
import argparse
import json
from itertools import product
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def spectral_interval(down, lepton):
    # Both sets are in increasing mass order.
    return [float(max(lepton[:-1] / down[1:])), float(min(lepton[1:] / down[:-1]))]


def rank_one_witness(down, lepton):
    # Common first two rows. Complete the Gram matrix by one row for either spectrum.
    d, e = down[::-1]**2, lepton[::-1]**2
    low = np.maximum(d[1:],e[1:])
    high = np.minimum(d[:-1],e[:-1])
    assert np.all(low < high)
    poles = np.r_[np.sqrt(low*high), 0.]
    def complete(spectrum):
        weights = [-np.prod(pole-spectrum)/np.prod(pole-np.delete(poles,i))
                   for i,pole in enumerate(poles)]
        assert min(weights) > 0
        rows = np.diag(np.sqrt(poles))[:2]
        return np.vstack([rows,np.sqrt(weights)])
    return complete(d),complete(e)


def median_normalization(d, e, w):
    vals = np.linalg.svd(w*e[None,:]/d[:,None],compute_uv=False)
    return float(vals[1])


def invariant_from_rates(d, e, rates, k):
    # Normalize the generalized Gram eigenvalues to avoid determinants of tiny Yukawas.
    prod = float(np.prod(e/d)**2)
    tr = float(np.sum(rates*(e[None,:]/d[:,None])**2))
    tr_inv = float(np.sum(rates*(d[:,None]/e[None,:])**2))
    x = k*k
    return x**3-tr*x**2+prod*tr_inv*x-prod


def run():
    cfg = json.loads((HERE/"inputs.json").read_text(encoding="utf-8"))
    frozen = json.loads((HERE/cfg["prior"]).read_text(encoding="utf-8"))
    values = np.array(frozen["scales"][cfg["scale"]])
    centers = values[:,0]
    d,e = centers[:3],centers[3:]
    corners = [centers + 2*values[:,1]*s for s in product([-1,1],repeat=6)]
    allowed = np.array([spectral_interval(z[:3],z[3:]) for z in corners])
    aligned = [median_normalization(z[:3],z[3:],np.eye(3)) for z in corners]
    yd,ye = rank_one_witness(d,e)
    b = (yd-ye)/4
    a = (3*yd+ye)/4
    assert np.linalg.matrix_rank(b) == 1
    assert np.allclose(a+b,yd) and np.allclose(a-3*b,ye)
    ud,sd,_ = np.linalg.svd(yd)
    ue,se,_ = np.linalg.svd(ye)
    assert np.allclose(sd,d[::-1]) and np.allclose(se,e[::-1])
    w = (ud.T@ue)[::-1,::-1]  # Both left mass bases in ascending order.
    assert abs(median_normalization(d,e,w)-1)<1e-9
    rng = np.random.default_rng(1433)
    for _ in range(50):
        test,_ = np.linalg.qr(rng.normal(size=(3,3)) + 1j*rng.normal(size=(3,3)))
        # Well-conditioned test of all three polynomial roots from rates.
        q=np.array([1.,2.,3.]); ell=np.array([0.8,2.2,4.])
        roots=np.linalg.svd(test*ell[None,:]/q[:,None],compute_uv=False)
        for k in roots:
            assert abs(invariant_from_rates(q,ell,abs(test)**2,k)) < 1e-10
        assert np.isclose(median_normalization(q,ell,test.conj()),roots[1])
    return {"round":1433,"scale_GeV":1e16,
            "kappa_spectra_only_central":spectral_interval(d,e),
            "kappa_spectra_only_common_to_entire_box":[max(allowed[:,0]),min(allowed[:,1])],
            "light_to_heavy_aligned_e_over_d":list(e/d),
            "kappa_required_aligned_center":median_normalization(d,e,np.eye(3)),
            "kappa_required_aligned_two_widths":[min(aligned),max(aligned)],
            "kappa_one_example": {"left_rate_matrix":np.round(abs(w)**2,10).tolist(),
                                  "rank_B":int(np.linalg.matrix_rank(b))}}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    flags=parser.add_mutually_exclusive_group()
    flags.add_argument("--check",action="store_true")
    flags.add_argument("--write-new",action="store_true")
    args=parser.parse_args()
    data=run()
    path=HERE/"results.json"
    if not args.write_new:
        assert data == json.loads(path.read_text(encoding="utf-8"))
        print("1433 two-sided spectral feasibility and modulus-only orientation test verified")
    else:
        with path.open("x",encoding="utf-8") as stream:
            json.dump(data,stream,indent=2)
            stream.write("\n")
