"""Bound the OPAL off-pole sign test over the entire accepted hard-energy interval."""
import itertools
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "independent_results.json"


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def scalar_interval(center, error):
    return (center - 2 * error, center + 2 * error)


def afb(same, opposite):
    return 0.75 * (same - opposite) / (same + opposite)


def positive_square_bounds(x, y):
    return (0.0 if x <= 0 <= y else min(x*x, y*y), max(x*x, y*y))


def angular_norm_distance(a0, a1):
    # Nearest vector on a ray with the required forward probability; total rate is free.
    p0, p1 = (1 + a0)/2, (1 + a1)/2
    return max(0.0, math.sqrt(p1*(1-p0)) - math.sqrt(p0*(1-p1)))


def calculate():
    weak = read(HERE.parent / "1394/inputs.json")["formula"]["Gmu_GeV_minus2"]
    mz, mz_error = read(HERE.parent / "1475/inputs.json")["Z_mass_GeV"]
    e = [scalar_interval(.26963, .00030), scalar_interval(.23148, .00029)]
    tau = [scalar_interval(.26930, .00058), scalar_interval(.23274, .00062)]
    # Pole chiral couplings in (gV +/- gA)/2 convention; emit both signs in each tau channel.
    s_lo, s_hi = .7225 * 188.5**2, 188.7**2
    m_lo, m_hi = mz - 2*mz_error, mz + 2*mz_error
    c_lo = weak * m_lo**2 * 125 / (2*math.sqrt(2)*math.pi)
    c_hi = weak * m_hi**2 * 130 / (2*math.sqrt(2)*math.pi)
    # Real Breit-Wigner bounded by x/(x-M²) times 1/(1+r²); Im/Re = -r.
    r_max = 2.51*m_hi/(s_lo-m_hi**2)
    re_lo = c_lo*s_hi/(s_hi-m_lo**2)/(1+r_max**2)
    re_hi = c_hi*s_lo/(s_lo-m_hi**2)
    im_max = re_hi*r_max
    data = [.591-2*(.054+.012), .591+2*(.054+.012)]
    cases = {}
    for f_left, f_right in itertools.product((1, -1), repeat=2):
        h = []
        for i, j in itertools.product(range(2), repeat=2):
            low, high = 4*e[i][0]*tau[j][0], 4*e[i][1]*tau[j][1]
            sign = (1 if i == j else -1) * (f_left, f_right)[j]
            re = sorted([1+sign*low*re_lo, 1+sign*high*re_hi])
            a, b = positive_square_bounds(*re)
            h.append((a, b+(high*im_max)**2))
        direct = [h[0][k] + h[3][k] for k in range(2)]
        crossed = [h[1][k] + h[2][k] for k in range(2)]
        window = [afb(direct[0], crossed[1]), afb(direct[1], crossed[0])]
        # Independent evaluation of a physical point, integrating forward and backward analytically.
        s = 188.6**2
        chi = weak*mz**2*128/(2*math.sqrt(2)*math.pi)*s/complex(s-mz**2,mz*2.4952)
        nominal = [abs(1+4*ge*gt*eta*chi)**2 for ge in [-.26963,.23148]
                   for gt, eta in zip([-.26930,.23274], [f_left, f_right])]
        a0 = afb(nominal[0]+nominal[3], nominal[1]+nominal[2])
        assert window[0] <= a0 <= window[1]
        cases[f"{f_left},{f_right}"] = dict(hard_AFB_enclosure=window, nominal_at_188_6=a0,
                    response_needed=max(0., data[0]-window[1]),
                    minimum_extra_amplitude_norm_if_unfolding_error_le_0_03=
                        angular_norm_distance(window[1], data[0]-.03))
    assert cases["1,1"]["response_needed"] == 0
    assert all(v["response_needed"] > .3 for k,v in cases.items() if k != "1,1")
    return {"pole_abs_e_LR": e, "pole_abs_tau_LR": tau, "inverse_alpha_matching_range": [125,130],
            "energy_squared_range_GeV2": [s_lo,s_hi], "propagator_bounds_Re": [re_lo,re_hi],
            "abs_Im_ceiling": im_max, "data_two_marginal_widths": data, "tau_signs_relative_to_e": cases,
            "sources": ["https://arxiv.org/abs/hep-ex/0509008","https://arxiv.org/abs/hep-ex/0309053"]}


if __name__ == "__main__":
    report = json.loads(json.dumps(calculate()))
    if sys.argv[1:] == ["--write-new"]:
        with OUT.open("x", encoding="utf-8") as handle:
            json.dump(report, handle, indent=2)
            handle.write("\n")
    elif sys.argv[1:] in ([], ["--check"]):
        assert report == read(OUT)
        print("1506 independent: tau sign branches and inclusive-amplitude norm bounds verified.")
    else:
        raise SystemExit("Use --write-new or --check.")
