"""Independent triplet audit: basis-invariant exterior power, EW threshold, spacelike current."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(relative):
    return json.loads((HERE / relative).read_text(encoding="utf-8-sig"))


def gauge_shift(masses, g0, mu0):
    """Sum the FPS (m²=2mu_H²) weak pole-matching entries at their own scales."""
    low, high = sorted(masses)
    out = []
    for m in (low, high):
        inv = 1/g0**2 + (19/6*math.log(m/mu0) - 4/3*math.log(m/low))/(8*math.pi**2)
        out.append(-6*m*m / ((4*math.pi)**4*inv*inv))
    return sum(out)


def inverse_alpha_delta(qlo, qhi, mass, charge_trace=1.):
    # Massive Dirac loop. log1p preserves the difference even at PeV scales.
    z, weights = np.polynomial.legendre.leggauss(24)
    x = (1+z)/2
    w = x*(1-x)
    integrand = w*np.log1p(w*(qhi-qlo)/(mass**2+w*qlo))
    return float(charge_trace/math.pi * np.dot(weights,integrand))


def check_two_columns(masses_eV, P, v):
    masses = np.array(masses_eV)*1e-9
    ratio = (masses[1]/masses[0])**(1/3)
    # Give the smaller light mass to the higher heavy mass.
    heavy = P*np.array([math.sqrt(ratio), 1/math.sqrt(ratio)])
    lower = 4*P**3*math.sqrt(masses[0]*masses[1])/v**2
    coeffs = np.array([[1.,3.,2.],[3.,-1.,0.],[.2,-.7,4.]])
    unitary = np.linalg.qr(coeffs*np.exp(1j*np.arange(9).reshape(3,3)**2))[0][:,:2]
    ans = []
    for angle in (0j, .54+.8j):
        c, s = np.cos(angle), np.sin(angle)
        R = np.array([[c,s],[-s,c]])
        B = (unitary*np.sqrt(masses)) @ R
        Y = (B * np.sqrt(2*heavy)) / v
        gram = Y.conj().T @ Y
        light = B @ B.T
        sv = np.linalg.svd(light, compute_uv=False)
        trace = float(np.dot(heavy**2, np.diag(gram).real))
        det_ratio = float(np.linalg.det(gram).real / (4*P**2*masses.prod()/v**4))
        assert abs(det_ratio-1)<2e-13 and np.allclose(sv[:2], masses[::-1], rtol=5e-13,atol=1e-24)
        assert trace >= lower*(1-1e-12)
        ans.append({"trace_GeV2":float(trace), "determinant_to_exterior_formula":det_ratio})
    assert math.isclose(ans[0]["trace_GeV2"], lower, rel_tol=1e-13)
    return {"ratio_heavy_for_least_trace":ratio,"minimum_trace_GeV2":lower,
            "equal_heavy_minimum_trace_over_optimum":sum(masses)/(2*math.sqrt(masses.prod())),
            "positive_and_complex_examples":ans}


def run():
    r=read("../1611/results.json")
    inp=read("../1608/inputs.json")
    oscillation=read("../1370/inputs.json")
    lo,hi=r["whole_triplet_geometric_mean_GeV"]
    mid=r["central_triplet_geometric_mean_GeV"]
    g0=inp["reference"]["g2"]
    mu0=inp["matching_scale"]
    mh_upper=inp["reference"]["mh"]+2*inp["marginal_width"]["mh"]
    shifts={s:abs(gauge_shift((p,p),g0,mu0)) for s,p in (("low",lo),("central",mid),("high",hi))}
    gauge = {"at_degenerate_pev_masses_GeV2":shifts,"lower_over_mh2":shifts["low"]/mh_upper**2,
             "relative_precision_of_independent_counterterm":mh_upper**2/shifts["low"],
             "scale_over_M_for_zero_explicit_gauge_term":math.exp(-1/12)}
    for p in (lo,mid,hi):
        for split in np.linspace(0,math.log(p/mu0),501):
            val=abs(gauge_shift((p*math.exp(-split),p*math.exp(split)),g0,mu0))
            assert val>=abs(gauge_shift((p,p),g0,mu0))*(1-1e-12)
    src=oscillation["NO"]
    no=(math.sqrt(src["solar"][0]),math.sqrt(src["atm_abs"][0]))
    src=oscillation["IO"]
    io=(math.sqrt(src["atm_abs"][0]-src["solar"][0]),math.sqrt(src["atm_abs"][0]))
    current={}
    for a,b in [(12.25,3434.),(1800.,21600.)]:
        key=f"{a:g}_{b:g}_GeV2"
        ratio_at_end=inverse_alpha_delta(a,b,mu0)+inverse_alpha_delta(a,b,lo**2/mu0)
        quark_bound=sum(inverse_alpha_delta(a,b,m,1/3) for m in (1100,2100))
        current[key]={"minimum_triplet_inverse_alpha_change":2*inverse_alpha_delta(a,b,hi),
                      "maximum_triplet_change":ratio_at_end,"maximum_two_quark_change":quark_bound,
                      "extension_upper_at_this_order":ratio_at_end+quark_bound}
    return {"gauge":gauge,"two_columns_at_central_P":{
        "NO":check_two_columns(no,mid,246.22),"IO":check_two_columns(io,mid,246.22)},
        "inverse_alpha_differences":current,
        "L3_2000_combined_marginal_error":math.hypot(.61,1.14)}


def compare(a,b):
    if isinstance(a, dict):
        assert a.keys()==b.keys()
        for k in a:
            compare(a[k],b[k])
    elif isinstance(a, list):
        assert len(a)==len(b)
        for x,y in zip(a,b):
            compare(x,y)
    else:
        assert math.isclose(a,b,rel_tol=2e-10,abs_tol=1e-18),(a,b)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    out=run()
    path=HERE/"independent_results.json"
    if args.check:
        compare(out,read("independent_results.json"))
        print("1612 independent OK: exterior determinant, full split bound and spacelike endpoints")
    else:
        with path.open("x",encoding="utf-8") as stream:
            json.dump(out,stream,indent=2)
            stream.write("\n")
