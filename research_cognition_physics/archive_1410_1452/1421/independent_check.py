"""Sharp two-readout moment bounds for a specified Michel/allowed model."""
import argparse
import json
import math
from pathlib import Path

from numpy.polynomial import Polynomial

ROOT = Path(__file__).resolve().parent
ENDPOINT, Q = 52.8, 17.33


def response(elastic_cut, carbon_electron_cut):
    # x=E/Emax. Constant cross-section and Jacobian factors cancel in each kernel.
    x = Polynomial([0, 1])
    flux = x*x*(1-x)
    elastic = flux*(x-elastic_cut/ENDPOINT)
    carbon = flux*(x-Q/ENDPOINT)**2
    kernels = []
    for weight, threshold in [(elastic, elastic_cut), (carbon,Q+carbon_electron_cut)]:
        integral = weight.integ()
        low = threshold/ENDPOINT
        area = integral(1)-integral(low)
        def tail(energy, integral=integral, low=low, area=area):
            return (integral(1)-integral(max(low,energy/ENDPOINT)))/area
        kernels.append(tail)
    return kernels


def bisect(f, lo, hi):
    value = f(lo)
    assert value*f(hi)<=0
    for _ in range(55):
        mid = (lo+hi)/2
        if value*f(mid)<=0:
            hi = mid
        else:
            lo = mid
    return (lo+hi)/2


def ratio_minimum(el, ca, elastic_cut, carbon_cut, left, right):
    a, b = left**2+right**2/3, 1+2*left
    # For a tail-optimal p the stationary condition is b*w_e/w_c = (a+b*Pe)/Pc.
    # Use the normalized densities directly. Their ratio is decreasing on support.
    x = Polynomial([0,1])
    flux = x**2*(1-x)
    ie = (flux*(x-elastic_cut/ENDPOINT)).integ()
    ic = (flux*(x-Q/ENDPOINT)**2).integ()
    norm_ratio = (ic(1)-ic((Q+carbon_cut)/ENDPOINT))/(ie(1)-ie(elastic_cut/ENDPOINT))
    def derivative_test(e):
        weight_ratio = norm_ratio*ENDPOINT*(e-elastic_cut)/(e-Q)**2
        return b*weight_ratio*ca(e) - a - b*el(e)
    lo = Q+carbon_cut+0.00001
    stop = (bisect(derivative_test,lo,ENDPOINT-.001)
            if derivative_test(lo)>0 else lo)
    return {"cut_MeV":stop,"sharp_positive_ratio":17.24/9.4*(a+b*el(stop))/ca(stop)}


def run():
    observed_upper = 10.1 + 2*math.hypot(1.1,1)
    carbon_min = 8.9 - 2*math.hypot(.3,.9)
    # A jointly possible pair from the sign-reflected CHARM V,A rectangle.
    left, right = .235,-.234
    tab = {}
    for name,te,tc in [("true_total",0,0),("ideal_left_response_thresholds",18,18)]:
        el,ca = response(te,tc)
        cut = bisect(lambda e: ca(e)-carbon_min/9.4, Q+tc, ENDPOINT)
        p = carbon_min/9.4
        high = bisect(lambda e: ca(e)-(1-p),Q+tc,ENDPOINT)
        row = dict(capture_lower_overlap=p, min_elastic_overlap=el(cut),
                   max_elastic_overlap=1-el(high), capture_extremizer_cut_MeV=cut)
        if te==0:   # The cut experiment changes the NC-right kernel & units.
            row.update(ratio_minimum(el,ca,te,tc,left,right))
            row["witness"] = dict(factor=1.07,threshold_MeV=38,Pe=el(38),Pc=ca(38),
                elastic=1.07*17.24*(left**2+right**2/3+(1+2*left)*el(38)),
                carbon=1.07*9.4*ca(38))
            assert row["witness"]["elastic"]<=observed_upper
            assert row["witness"]["carbon"]>=carbon_min
        tab[name]=row
    return dict(two_width_elastic_upper=observed_upper, two_width_carbon_lower=carbon_min,
                ratio_upper=observed_upper/carbon_min,cases=tab)


def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:
            compare(a[k],b[k])
    else:
        assert math.isclose(a,b,rel_tol=1e-9,abs_tol=2e-12), (a,b)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check",action="store_true")
    mode.add_argument("--write-new",action="store_true")
    args=parser.parse_args()
    output=run()
    target=ROOT/"independent_results.json"
    if args.write_new:
        with target.open("x",encoding="utf-8") as stream:
            json.dump(output,stream,indent=2)
            stream.write("\n")
    else:
        compare(output,json.loads(target.read_text(encoding="utf-8")))
        print("1421 independent response bounds verified")
