"""Mass determinant orders and one-loop bounds for a single charged spurion."""
import argparse
from fractions import Fraction
import itertools
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent


def data(name):
    return json.loads((HERE/name).read_text(encoding="utf-8-sig"))


def minor_orders(rows, cols):
    # Each injective assignment appears once; independent coefficients are allowed.
    return [sum(abs(rows[i]+c) for i, c in zip(ix, cols))
            for ix in itertools.permutations(range(len(rows)),len(cols))]


def threshold(gauge, sign2, sign3, widths):
    g2 = gauge["g2"][0]+widths*sign2*gauge["g2"][1]
    g3 = gauge["g3"][0]+widths*sign3*gauge["g3"][1]
    a2, a3 = 4*math.pi/g2**2, 4*math.pi/g3**2
    b2,b3 = -19/6,-7
    length = 2*math.pi*(a2-a3)/(b2-b3)
    return gauge["mu_ref_GeV"]*math.exp(length), a2-b2*length/(2*math.pi)


def run():
    cfg = data("inputs.json")
    g = data(cfg["threshold"])
    charges = cfg["charges"]
    p, s = charges["X_perp"], charges["X_on"]
    matrix_charges = {"15": ([p,s,-s],[-s]*3),
                      "24": ([p,p,s,s,-s,-s],[p,p,s,s,-s,-s]),
                      "10": ([p]*3+[s]*3,[p,s,-s])}
    orders = {rep: min(minor_orders(*ends)) for rep,ends in matrix_charges.items()}
    assert orders == {"15":35,"24":16,"10":25}
    q_max = math.gcd(abs(p-s),abs(2*s))
    anomaly = 13*p+13*s-17*s
    slopes = {"15":Fraction(14,3),"24":Fraction(10,3),"10":Fraction(2)}
    weighted = sum(slopes[rep]*k for rep,k in orders.items())
    assert q_max==1 and anomaly==0 and weighted==Fraction(800,3)
    mz_to_match = threshold(g,0,0,cfg["widths"])
    corners = [threshold(g,i,j,cfg["widths"]) for i,j in itertools.product((-1,1),repeat=2)]
    windows=[]
    for C in cfg["entry_caps"]:
        factors = {"15":3**1.5*C**3,"24":math.factorial(6)*C**6,"10":6*math.sqrt(20)*C**3}
        log_B = sum(float(slopes[rep])*math.log(val) for rep,val in factors.items())
        for inv in cfg["alpha_parent_inverse_floor"]:
            def eps_floor(a23):
                return math.exp(-(2*math.pi*(a23-3*inv)+log_B)/float(weighted))
            center = eps_floor(mz_to_match[1])
            floor_list = [eps_floor(a23) for _,a23 in corners]
            # A sparse simultaneously realizable spectrum has all matching coefficients 1:
            # M15~(1,eps^9,eps^26), M24~(1,1,1,1,eps^8,eps^8),
            # M10~(1,eps^8,eps^17), and M5~(1,1,1).
            sparse_cross = math.exp(-2*math.pi*(mz_to_match[1]-3*inv)/float(weighted))
            windows.append({"coefficient_cap":C,"min_inverse_parent":inv,
                            "necessary_epsilon_floor_central":center,
                            "necessary_epsilon_floor_gauge_box":[min(floor_list),max(floor_list)],
                            "sparse_unit_coefficient_crossing_central":sparse_cross,
                            "sparse_mass_floor_at_crossing_GeV":mz_to_match[0]*sparse_cross**26})
    edm = data(cfg["edm"])["edm"]
    theta_floor = edm["abs_slope_e_fm"][0]-2*sum(edm["abs_slope_e_fm"][1:])
    return {"full_gapping_charge_divisor":q_max,"mass_product_minimal_powers":orders,
            "weighted_log_threshold_coefficient":str(weighted),"SU15_anomaly":anomaly,
            "M23_GeV_and_SM_inverse_alpha":list(mz_to_match),"threshold_windows":windows,
            "theta_upper_if_other_neutron_dipoles_vanish":
                edm["limit_90_ecm"]/(edm["fm_in_cm"]*theta_floor)}


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--write-new",action="store_true")
    group.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=run()
    out=HERE/"results.json"
    if args.write_new:
        with out.open("x",encoding="utf-8") as handle:
            json.dump(result,handle,indent=2)
            handle.write("\n")
    else:
        assert json.loads(out.read_text())==result
    print(json.dumps(result))
