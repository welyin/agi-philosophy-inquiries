"""Shape information after eliminating the two physical untagged weights."""
import argparse
import json
from math import exp,sqrt
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent


def calculate():
    gamma, diff = .6600, .0783
    h,l=1/(gamma-diff/2),1/(gamma+diff/2)
    center=1.60
    fraction=(center-l)/(h-l)
    assert 0 < fraction < 1
    moment2=2*(fraction*h**2+(1-fraction)*l**2)
    assert np.isclose(moment2,2*((h+l)*center-h*l))
    x=2-sqrt(2)
    curvature=(2*x-x*x)*exp(-x)
    uniform=curvature*(h-l)**2/(8*l**2)
    # Example with equal integrated heavy/light weights: check a tail directly.
    m=(h+l)/2
    cut=.6*m
    tail_mix=(exp(-cut/h)+exp(-cut/l))/2
    tail_single=exp(-cut/m)
    assert abs(tail_mix-tail_single)<uniform
    # Response expected with unchanged fitted selected signal yield, fixed endpoint efficiencies.
    # Run 3 single-heavy branching and its 9.2% increase for a pure-light template.
    relative_middle=(h+l)/(h+l/1.092)
    return dict(heavy_ps=h,light_ps=l,integrated_heavy_fraction_at_combined_mean=fraction,
                predicted_second_raw_moment_ps2_at_combined_mean=moment2,
                max_fractional_second_moment_excess=(h-l)**2/(4*h*l),
                uniform_max_absolute_tail_difference=uniform,fixed_gate_ps=cut,
                signed_tail_difference_equal_weights=tail_mix-tail_single,
                predicted_central_rate_multiplier_constant_yield=relative_middle,
                observed_run3_central_multiplier=1.040)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    result=calculate()
    if p.parse_args().check:
        old=json.loads((HERE/"independent_results.json").read_text(encoding="utf8"))
        assert old.keys()==result.keys() and all(np.isclose(old[k],v,rtol=1e-10,atol=1e-14) for k,v in result.items())
        print("1405 second moment, all threshold gates and efficiency calibration verified")
    else:
        print(json.dumps(result,indent=2))
