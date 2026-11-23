"""Color threshold secants and their Higgs derivative, independently by interval integration."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np


def inverse_difference(low, high, masses, kappas, background=0.0):
    # Integrate the changed beta function on its actual constant pieces.
    poles = [m*math.exp(k*background) for m,k in zip(masses, kappas)]
    points = sorted([low, high] + [m for m in poles if low < m < high])
    area = 0
    for left, right in zip(points, points[1:]):
        active = sum(m < math.sqrt(left*right) for m in poles)
        area -= (2/3)*active*math.log(right/left)/(2*math.pi)
    return area


def curvature(energies, inverse_differences):
    low, mid, high = energies
    first, second = inverse_differences
    return first/math.log(mid/low)-second/math.log(high/mid)


def evaluate():
    energies = [474.,664.,896.]
    central = [.0936,.0894,.0889]
    width = [.0041,.0031,.0034]
    def measured(data):
        inv = 1/np.array(data)
        return float(curvature(energies,np.diff(inv)))
    corners = [measured([c+sign*w for c,sign,w in zip(central,signs,width)])
               for signs in itertools.product((-1,1),repeat=3)]
    def steps(poles,kappas):
        return [inverse_difference(lo,hi,poles,kappas) for lo,hi in zip(energies,energies[1:])]
    peak = curvature(energies,steps([664.,664.],[1,1]))
    np.testing.assert_allclose(peak, 2/(3*math.pi))
    for m in [200.,474.,896.,2000.]:
        np.testing.assert_allclose(curvature(energies,steps([m],[1])),0,atol=2e-15)
    # Bracket the two fundamental-color poles; the lower and higher spectator drop out.
    all_masses = [300.,520.,800.,1200.]
    all_kappa = [.3,.75,1.,-1.]
    def diff(eta):
        return inverse_difference(474.,896.,all_masses,all_kappa,eta)
    eps = 1e-4
    slope = (diff(eps)-diff(-eps))/(2*eps)
    expected = (2/3)*sum(all_kappa[1:3])/(2*math.pi)
    np.testing.assert_allclose(slope, expected, rtol=1e-9)
    # These same two Dirac thresholds have the same geometric mean.
    equal, split = [560.,560.], [500.,560.0**2/500.]
    d0, d1 = steps(equal,[0,0]), steps(split,[0,0])
    np.testing.assert_allclose(d0,d1,atol=2e-15)
    midread = [inverse_difference(474.,540.,p,[0,0]) for p in (equal,split)]
    assert abs(midread[0]-midread[1]) > .005
    return {"CMS_inverse_secant_difference":measured(central),
            "CMS_marginal_rectangle_range":[min(corners),max(corners)],
            "two_fundamental_maximum_one_loop_curvature":peak,
            "windowed_background_derivative":slope,
            "heavy_hgg_weight_over_alpha_s_over_8pi":(2/3)*sum(all_kappa[1:3]),
            "equal_and_split_poles_GeV":[equal,split], "shared_inverse_steps": d0,
            "extra_inverse_readout_474_to_540":midread}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    path = Path(__file__).with_name("independent_results.json")
    result = evaluate()
    if args.write:
        with path.open("x", encoding="utf-8") as out:
            json.dump(result, out, indent=2, allow_nan=False)
            out.write("\n")
    else:
        saved = json.loads(path.read_text(encoding="utf-8"))
        assert saved.keys()==result.keys()
        for k,v in result.items():
            np.testing.assert_allclose(v, saved[k], rtol=1e-9, atol=2e-14)
        print("1379: color slope triangle, threshold derivative and two spectra verified.")
