"""Four-dimensional anomalies, index cost and marginal two-bin check."""
import argparse
import itertools
import json
import math
from fractions import Fraction as F
from pathlib import Path

HERE = Path(__file__).resolve().parent
# Colour sign, weak dimension, six times hypercharge.
REPS = {"Q": (1, 2, 1), "U": (-1, 1, -4), "D": (-1, 1, 2),
        "L": (0, 2, -3), "E": (0, 1, 6), "N": (0, 1, 0)}


def anomalies(rows):
    keys = ("333", "33Y", "22Y", "gravY", "YYY", "33X", "22X", "YYX",
            "YXX", "gravX", "XXX", "witten")
    totals = dict.fromkeys(keys, 0)
    index = [F(0), F(0), F(0)]
    for (colour, weak, six_y), degree, copies in rows:
        nc = 3 if colour else 1
        dim = nc * weak
        y = F(six_y, 6)
        strong, left = abs(colour)*weak, (weak-1)*nc
        vector = (colour*weak, strong*y, left*y, dim*y, dim*y**3,
                  strong*degree, left*degree, dim*y*y*degree, dim*y*degree**2,
                  dim*degree, dim*degree**3, left)
        for key, value in zip(keys, vector):
            totals[key] += copies*value
        for j, val in enumerate((F(6,5)*dim*y*y, left, strong)):
            index[j] += copies*val
    totals["witten"] %= 2
    return totals, index


def beta(extra):
    return [v/3 for v in anomalies(extra)[1]]


def run():
    cfg = json.loads((HERE/"inputs.json").read_text(encoding="utf8"))
    d = cfg["degree"]
    signs = list(itertools.product((-1, 1), repeat=6))
    trial = []
    for sign in signs:
        rows = [(rep, d*s, d) for rep,s in zip(REPS.values(),sign)]
        a, _ = anomalies(rows)
        trial.append((sign, a["33X"], a["22X"]))
    same_up = [r for r in trial if r[0][0] == r[0][1]]
    same_both = [r for r in same_up if r[0][0] == r[0][2]]
    floors = [min(abs(r[2]) for r in trial), min(abs(r[1]) for r in same_up),
              min(abs(r[1]) for r in same_both)]
    assert floors == [18,18,36]
    for sign, _, _ in trial:
        q,u,dn = (3*s for s in sign[:3])
        if (q+u) == -(q+dn):
            assert u == dn == -q
    ordinary = [(r,d,d) for r in REPS.values()]
    extra = []
    for colour,weak,y in ((1,1,-2),(0,2,3)):
        for conjugate in (-1,1):
            extra.append(((conjugate*colour,weak,conjugate*y),-d,2*d))
    extra.append((REPS["N"],d,4*d))
    a, _ = anomalies(ordinary+extra)
    assert all(value == 0 for value in a.values())
    a0, _ = anomalies(ordinary)
    ordinary_mixed = [F(6,5)*a0["YYX"],a0["22X"],a0["33X"]]
    scalar = cfg["scalar_charge"]
    assert ordinary_mixed == [36]*3
    assert beta(extra) == [F(2,3)*abs(x)/scalar for x in ordinary_mixed] == [4]*3
    # Opposite down-quark and lepton degree signs saturate the weaker SU(2,3) bounds.
    base2 = [(r,s*d,d) for r,s in zip(REPS.values(),(1,1,-1,-1,1,1))]
    less = []
    for colour,weak,y,weight in ((1,1,4,-d),(0,2,3,-d),(0,1,6,d)):
        for conj in (-1,1):
            less.append(((conj*colour,weak,conj*y),weight,d))
    less.append((REPS["N"],d,2*d))
    assert not any(anomalies(base2+less)[0].values())
    assert beta(less) == [F(34,5),2,2]
    low,high = cfg["cms_momentum_GeV"]
    a1,a2 = cfg["cms_alpha3"]
    w1,w2 = cfg["cms_marginal_widths"]
    delta = math.log(high/low)/(2*math.pi)
    observed = 1/a2-1/a1
    interval = [1/(a2+w2)-1/(a1-w1),1/(a2-w2)-1/(a1+w1)]
    return dict(sign_choices=len(trial),absolute_mixed_index_floors=floors,
                universal_extra_weyl_index=[int(v) for v in anomalies(extra)[1]],
                universal_beta_shift=[int(v) for v in beta(extra)],
                weaker_beta_shift=[str(v) for v in beta(less)],
                anomaly_channels=list(a),universal_anomalies=[int(v) for v in a.values()],
                inverse_alpha_shift_one_decade=-4*math.log(10)/(2*math.pi),
                cms_inverse_change=observed,cms_one_width_marginal_box=interval,
                one_loop_bin_change={str(b):(7-b)*delta for b in (0,2,4)})


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    args = parser.add_mutually_exclusive_group(required=True)
    args.add_argument("--write-new",action="store_true")
    args.add_argument("--check",action="store_true")
    opt = parser.parse_args()
    result = run()
    output = HERE/"results.json"
    if opt.write_new:
        with output.open("x",encoding="utf8") as f:
            json.dump(result,f,indent=2)
            f.write("\n")
    else:
        assert json.loads(output.read_text()) == result
        print("1530 verified: 64 flux signs, both anomaly-free spectra, sharp index and threshold tests.")
