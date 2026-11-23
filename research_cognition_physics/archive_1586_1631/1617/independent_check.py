"""Independent audit. Columns, inclusive target and crossed LEP2 template."""
import argparse
import itertools
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
D = json.loads((HERE / "inputs.json").read_text(encoding="utf8"))


def make_result():
    n = D["residual"]
    widths = n["widths"]
    p_mu = (1 + math.cos(2 * math.pi * n["irrep_j"] / n["d"])) / 2
    theta13 = [math.asin(math.sqrt(a)) / 2
               for a in (n["fast"][0]-widths*n["fast"][1],
                         n["fast"][0]+widths*n["fast"][1])]
    gf = D["fermi_GeV_minus_2"]
    sigma_fb = (4 * gf**2 * 209**2 / (3 * math.pi) * 3.8937937e11
                * D["mu3e_ratio_90"])
    measured = {"emu": 22., "etau": 78.}
    ans = {}
    for order in ("solar_no", "solar_io_prior"):
        sun, err = n[order]
        es = [math.cos(h)**2*math.cos(s)**2
              for h,s in itertools.product(theta13,
                  [math.asin(math.sqrt(sun+sgn*widths*err)) for sgn in (-1,1)])]
        ratios = [(1-p_mu-e)/p_mu for e in es]
        tau = D["ordinary_tau_e"]
        correction = tau["combined"]*widths
        upper = [sigma_fb, sigma_fb*max(ratios)]
        ans[order] = {"electron_probability": [min(es),max(es)],
             "crossed_etau_over_emu": [min(ratios), max(ratios)],
             "tau3eB_over_mu3eR": [(tau["value"]-correction)*min(ratios),
                                    (tau["value"]+correction)*max(ratios)],
             "two_charge_cross_sections_at209GeV_max_fb": upper,
             "ideal_min_common_amplitude_gain_to_OPAL_limits":
                    [math.sqrt(x/y) for x,y in zip(measured.values(), upper)],
             "favorable_angular_response_only_gain": {
                str(transfer): [math.sqrt(x/(transfer*y)) for x,y in zip(measured.values(), upper)]
                for transfer in (1.125,1.5)}}
    # One charge: density relative to uniform. The conjugate charge reverses z.
    def density(z):
        return 3/16*(4+(1-z)**2)
    grid = [-1 + i/12 for i in range(25)]
    assert min(density(z) for z in grid) == .75
    assert max((density(z)+density(-z))/2 for z in grid) == 1.125
    # Finite CP phases in a general produced pure state: neutral projector
    # commutes with a vacuum propagation phase and with mass dephasing.
    e = sum(ans["solar_no"]["electron_probability"])/2
    probabilities = [e,p_mu,1-e-p_mu]
    c, q = .002,.23
    tagged = [complex(math.sqrt(p)*math.cos(a), math.sqrt(p)*math.sin(a))
              for p,a in zip(probabilities, [0,.79,-1.22])]
    scattered = [q*(i==1) + c*amp*tagged[1].conjugate() for i,amp in enumerate(tagged)]
    inclusive = sum(abs(x)**2 for x in scattered)
    assert math.isclose(inclusive,(q+c*p_mu)**2+c*c*p_mu*(1-p_mu),abs_tol=1e-15)
    return {"sectors":ans, "CHARM_inclusive_complex_witness":inclusive, "p_mu":p_mu,
            "comparison":"OPAL 200-209 GeV Table 3; two charge states; angular-only and nominal rate budgets separated."}


def rounded(x):
    if isinstance(x,dict):
        return {k:rounded(v) for k,v in x.items()}
    if isinstance(x,list):
        return [rounded(y) for y in x]
    return float(f"{x:.11g}") if isinstance(x,float) else x


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    opts = parser.parse_args()
    assert not(opts.write and opts.check)
    out = HERE / "independent_results.json"
    result = rounded(make_result())
    if opts.write:
        with out.open("x",encoding="utf8") as f:
            json.dump(result,f,ensure_ascii=False,indent=2)
            f.write("\n")
    else:
        assert json.loads(out.read_text(encoding="utf8"))==result
        print("PASS 1617 independent projector and crossed collider test")
