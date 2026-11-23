"""Test a finite-momentum CP-odd gluon/photon relation and a soft axion counterexample."""
import argparse
import cmath
import json
import math
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent


def form(s,m):
    z=s/(4*m*m)
    if z<=1:
        return math.asin(math.sqrt(z))**2/z
    beta=math.sqrt(1-1/z)
    return -(2*math.atanh(beta)-1j*math.pi)**2/(4*z)


def calc():
    inputs=json.loads((HERE/"inputs.json").read_text(encoding="utf-8-sig"))
    tar=json.loads((HERE.parent/"1532/inputs.json").read_text(encoding="utf-8-sig"))["benchmark"]
    mf=inputs["higgs_mass_GeV"]
    s=mf**2
    t=form(s,tar["top_MS_GeV"])
    ell=form(s,inputs["tau_mass_GeV"])
    endpoints=[]
    for mass in inputs["b_mass_work_range_GeV"]:
        b=form(s,mass)
        g=t-b
        photon=4*t/3-b/3+ell
        for local_g,local_gamma in [(0.,0.),(.7,-1.3),(-g.real,-photon.real)]:
            np.testing.assert_allclose((photon+local_gamma-(g+local_g)/3).imag,ell.imag,atol=1e-14)
        q=4*mass**2/s
        beta=math.sqrt(1-q)
        # d[q atanh(sqrt(1-q))]/dq = atanh(beta)-1/(2 beta).
        assert math.atanh(beta)>1/(2*beta)
        endpoints.append(dict(mb=mass,real_gluon=g.real,imag_gluon=g.imag,
                              absolute_gluon=abs(g),real_photon=photon.real,imag_photon=photon.imag))
    assert t.imag==0 and ell.imag>0
    # QCD trace cancellation at zero momentum, retaining the first nontrivial term.
    soft_s=1e-5
    b0=inputs["b_mass_work_range_GeV"][0]
    slope=(1/tar["top_MS_GeV"]**2-1/b0**2)/12
    np.testing.assert_allclose((form(soft_s,tar["top_MS_GeV"])-form(soft_s,b0))/soft_s,slope,rtol=1e-6)
    # A further colour-singlet, charge-one pseudoscalar Yukawa can cancel the tau cut.
    extra_mass=3.
    extra=-ell.imag/form(s,extra_mass).imag
    assert abs((ell+extra*form(s,extra_mass)).imag)<1e-14
    # Euclidean momentum: constant QCD susceptibility only in its low-energy range.
    axion={str(x):x/(1+x) for x in (0.,.01,1.,100.)}   # x=Q^2/m_a^2
    return dict(unit_heavy_fermion_normalization=True, mass_points=endpoints, tau_im=ell.imag,
                clean_absorptive_ratio_to_eta=ell.imag,
                ratio_real_contact_independent_amplitude_floor=3*ell.imag/4,
                extra_uncoloured_CPodd_Yukawa_per_eta=extra,example_axion_response=axion)


def equal(a,b):
    if isinstance(a,dict):
        return a.keys()==b.keys() and all(equal(x,b[k]) for k,x in a.items())
    if isinstance(a,list):
        return len(a)==len(b) and all(equal(x,y) for x,y in zip(a,b))
    return math.isclose(a,b,rel_tol=2e-11,abs_tol=1e-14)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    flags=parser.add_mutually_exclusive_group(required=True)
    flags.add_argument("--write",action="store_true")
    flags.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=calc()
    path=HERE/"independent_results.json"
    if args.write:
        with path.open("x",encoding="utf-8") as f:
            f.write(json.dumps(result,indent=2,allow_nan=False)+"\n")
    else:
        assert equal(result,json.loads(path.read_text(encoding="utf-8")))
    print("1534 independent on-shell and soft CP checks passed")
