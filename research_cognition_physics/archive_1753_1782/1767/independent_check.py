"""Invert neutral-meson observables exactly and test the inherited common weak-current factor."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def reconstruct(dm, dg, asymmetry):
    """Convention: dg = Gamma_L - Gamma_H and a = (1-|q/p|^4)/(1+|q/p|^4)."""
    r = ((1-asymmetry)/(1+asymmetry))**.25
    split = complex(dm, dg/2)
    h12, h21 = split/(2*r), r*split/2
    dispersive = (h12+h21.conjugate())/2
    absorptive = 1j*(h12-h21.conjugate())
    hamiltonian = np.array([[-.33j, dispersive-1j*absorptive/2],
                            [dispersive.conjugate()-1j*absorptive.conjugate()/2, -.33j]])
    eigen = np.linalg.eigvals(hamiltonian)
    np.testing.assert_allclose(eigen[np.argmax(eigen.real)]-eigen[np.argmin(eigen.real)], split, atol=1e-13)
    np.testing.assert_allclose((abs(h12)**2-abs(h21)**2)/(abs(h12)**2+abs(h21)**2), asymmetry,atol=1e-14)
    rate = (absorptive/dispersive)
    return dict(twice_M_over_dm=2*abs(dispersive)/dm, abs_Gamma12_ps_inv=abs(absorptive),
                Gamma_over_M=[rate.real,rate.imag])


def solve():
    clocks = json.loads((HERE.parent/"1736/independent_results.json").read_bytes())["data"]
    dm, edm = clocks["delta_ms"]
    dg, edg = .087, math.hypot(.012,.009)
    a, ea = .0039, math.hypot(.0026,.0020)
    cases = {}
    for n in (1,2):
        upper_a = abs(a)+n*ea
        s = (dg+n*edg)/(2*(dm-n*edm))
        upper = math.sqrt(1+.5*(1+s*s)*math.expm1(-.5*math.log1p(-upper_a**2)))
        checked = reconstruct(dm-n*edm, dg+n*edg, upper_a)["twice_M_over_dm"]
        assert abs(checked-upper)<3e-15
        cases[str(n)]=dict(max_abs_asymmetry=upper_a, deltaM_approx_relative_cap=upper-1)
    s = (dg+2*edg)/(2*(dm-2*edm))
    target = 1+370/10180
    required = math.sqrt(1-(1+2*(target*target-1)/(1+s*s))**-2)
    theory_dm, theory_dg, theory_a = 18.77, .091, .0000206
    matched=reconstruct(dm,dm*theory_dg/theory_dm,theory_a)
    baseline=reconstruct(theory_dm,theory_dg,theory_a)
    np.testing.assert_allclose(matched["Gamma_over_M"],baseline["Gamma_over_M"])
    return dict(central=reconstruct(dm,dg,a),corners=cases, required_abs_asymmetry_for_370_over_10180=required,
                common_rescaling_example=dict(dg_at_measured_dm=dm*theory_dg/theory_dm, a_sl=theory_a,
                                             Gamma_over_M=matched["Gamma_over_M"]))


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args=parser.parse_args()
    out=json.loads(json.dumps(solve()),parse_float=lambda x:float(f"{float(x):.10g}"))
    path=HERE/"independent_results.json"
    if args.check:
        assert out==json.loads(path.read_bytes())
        print("1767 independent: exact mixing inverse, two-width envelope, common-factor invariance")
    else:
        with path.open("x",encoding="utf8") as handle:
            json.dump(out,handle,indent=2,allow_nan=False)
            handle.write("\n")
