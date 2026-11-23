"""Sharp high-order tail bounds for electroweak-background dependent fermion masses."""
import argparse
import json
from math import pi, sqrt
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE/name).read_text(encoding="utf-8"))


def extremal_coefficients(y, k, eps, cap, direct):
    # Every n >= 2 projected coefficient is -cap. All series below converge for eps < 1.
    tail = -cap*eps**2/(1-eps)
    first = (y-direct-tail)/eps
    mass = direct+first*eps+tail
    radial = mass+2*(first*eps-cap*(eps/(1-eps)**2-eps))
    assert np.isclose(mass,y) and np.isclose(radial,k*y)
    return first


def run():
    cfg = read("inputs.json")
    ydata = read(cfg["mass_data"])["yukawa"]
    w = cfg["upper_branch_widths"]
    v = cfg["v_GeV"]
    out = {}
    for key, sector in (("bb","down"),("tautau","charged_leptons")):
        r = cfg["amplitude_data"][key]
        k_hi = sqrt((r["raw"]+w*r["plus"])/((r["ref"]-w*r["ref_err"])*cfg["rate_response_min"]))
        center,err=ydata[sector][2]
        y = center - cfg["lower_mass_widths"]*err
        strength = (3-k_hi)*y/2
        eps = v*v/(2*cfg["scale_example_GeV"]**2)
        per_cap={}
        for cap in cfg["cap_examples"]:
            threshold = sqrt(strength/cap)
            smallest_eps = threshold/(1+threshold)
            largest_scale = v/sqrt(2*smallest_eps)
            first_null = extremal_coefficients(y,k_hi,smallest_eps,cap,0)
            assert 0 < first_null <= cap
            direct_bound = strength-cap*(eps/(1-eps))**2
            first_one_tev = extremal_coefficients(y,k_hi,eps,cap,direct_bound)
            # Verify the envelope at this specified scale whenever it claims a positive floor.
            if direct_bound>0:
                assert max(abs(direct_bound),abs(first_one_tev)) <= cap
            per_cap[str(round(cap,6))] = {
                "no_direct_largest_scale_GeV":largest_scale,
                "no_direct_minimum_epsilon":smallest_eps, "first_coefficient_at_boundary":first_null,
                "one_TeV_direct_net_fraction_floor":max(0.,direct_bound/y)}
        out[key]={"kappa_upper":k_hi,"y_lower":y,
                  "one_TeV_minimum_uniform_cap_if_no_direct":strength*((1-eps)/eps)**2,
                  "cap_examples":per_cap}
    # Embed a saturating scalar block in falsely misaligned three-flavour bases.
    rng=np.random.default_rng(1438)
    def unitary():
        return np.linalg.qr(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))[0]
    for _ in range(12):
        l,r=unitary(),unitary()
        coeff=[l@np.diag([.01,-.02,-1 if n>=2 else .3])@r.conj().T for n in range(7)]
        x=.08
        ymat=sum(m*x**n for n,m in enumerate(coeff))
        hmat=sum((2*n+1)*m*x**n for n,m in enumerate(coeff))
        left,mass,right=np.linalg.svd(ymat)
        z=left[:,0].conj()@hmat@right.conj().T[:,0]
        c0=left[:,0].conj()@coeff[0]@right.conj().T[:,0]
        assert c0.real >= (3*mass[0]-z.real)/2-(x/(1-x))**2-1e-12
    return out


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    flags=cli.add_mutually_exclusive_group()
    flags.add_argument("--check",action="store_true")
    flags.add_argument("--write-new",action="store_true")
    args=cli.parse_args()
    data=run()
    if args.write_new:
        with (HERE/"results.json").open("x",encoding="utf-8") as stream:
            json.dump(data,stream,indent=2)
            stream.write("\n")
    else:
        assert data == read("results.json")
        print("1438 Wilson tail envelope, saturation and actual bottom/tau windows verified")
