"""Independent 1447 check using eliminated seed and analytic cubic in kappa."""
import argparse
from itertools import permutations
import json
import math
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent


def all_slopes(bottom, ratio):
    # Eliminate the two seed masses from kappa, directly solve polynomial in kappa_s.
    h=(bottom-1)/(3-bottom)**3*ratio**2
    roots=np.roots([-h,9*h,-27*h-1,27*h+1])
    return sorted(z.real for z in roots if abs(z.imag)<1e-7)


def test():
    inp=json.loads((HERE/"inputs.json").read_text(encoding="utf-8"))
    def older(key):
        return json.loads((HERE/inp[key]).read_text(encoding="utf-8"))
    rates=older("direct_ratio")["direct_ATLAS_ratio_boxes"]["bb_ZZ"]["2"]["normalized_rate_interval"]
    bk=[math.sqrt(rates[0]/1.1), math.sqrt(rates[1]/.9)]
    y=older("mass_calibration")["yukawa"]["down"]
    ratio=(y[1][0]+2*y[1][1])/(y[2][0]-2*y[2][1])*(1+inp["max_relative_transport_mZ_to_h"])
    ordinary=[min(all_slopes(k,ratio),key=lambda s:abs(s-1)) for k in bk]
    branches=all_slopes(7/9,.019)
    assert len(branches)==3 and branches[0]<-350 and branches[2]>350
    for b in branches:
        t=(b-1)/(3-b)
        assert np.isclose(t*(1+t)**2,-.1*.9**2*.019**2,rtol=1e-6)
    # Direct SVD/left-eigenvector invariance for a complex nonsymmetric tensor.
    rng=np.random.default_rng(1447)
    frame = lambda: np.linalg.qr(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))[0]
    u,v=frame(),frame()
    seed=np.array([.007,.075,.63])
    matrix=(u*seed)@v.conj().T
    deform=matrix+0.83*matrix@matrix.conj().T@matrix
    proj=u.conj().T@deform@deform.conj().T@u
    assert np.max(abs(proj-np.diag((seed+.83*seed**3)**2)))<1e-12
    # A radial counterterm neutral at v, with flavour polynomial vanishing on bottom.
    # At v=1, (h^2-1)h x(x^2-b^2) has h derivative 2x(x^2-b^2).
    shift=(1.4-1)/(2*(seed[1]**2-seed[2]**2))
    kap=1+2*shift*(seed**2-seed[2]**2)
    assert np.allclose(kap[1:],[1.4,1])
    # A phase-carrying Hermitian example, exact polynomial identity by 27 paths.
    a=np.array([[.001,.003+.002j,0],[.003-.002j,-.04,.087j],[0,-.087j,1]])
    multi=sum(a[0,i]*a[i,j]*a[j,2] for i in range(3) for j in range(3))
    assert np.isclose(multi,np.trace(a)*a[0,1]*a[1,2])
    return {"kappa_b":bk,"strange_to_bottom_max":ratio,"positive_connected_kappa_s":ordinary,
         "all_three_signed_kappa_s_same_bottom_and_mass_ratio":branches,
         "radial_counterterm_strange_and_bottom_at_same_masses":kap[1:].tolist(),
         "largest_svd_residual":float(np.max(abs(proj-np.diag((seed+.83*seed**3)**2))))}


if __name__=="__main__":
    p=argparse.ArgumentParser()
    ctl=p.add_mutually_exclusive_group()
    ctl.add_argument("--write",action="store_true"); ctl.add_argument("--check",action="store_true")
    args=p.parse_args()
    res=test()
    res["largest_svd_residual"]=round(res["largest_svd_residual"],12)
    out=json.dumps(res,indent=2,allow_nan=False)+"\n"
    target=HERE/"independent_results.json"
    if args.write:
        with target.open("x",encoding="utf-8") as f: f.write(out)
    elif out!=target.read_text(encoding="utf-8"):
        raise SystemExit("1447 independent mismatch")
    else: print("1447 independent OK")
