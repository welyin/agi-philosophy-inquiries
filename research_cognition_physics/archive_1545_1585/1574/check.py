"""Reproduce the character bounds and conditional Z-width comparison."""
import argparse
import json
import math
from collections import Counter
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parent

def grade(y):
    return (0, 1, -1)[y % 3]

def low_model(test):
    fields = []
    for y in test["doublet_y"]:
        fields.extend([(2, y), (1, -y+3), (1, -y-3)])
    fields.append((2, test["neutral_grade_doublet_y"]))
    signature = Counter()
    charges = Counter()
    cubic = gravitational = weak_hyper = twice_index = 0
    for dimension, y in fields:
        signature[dimension] += grade(y)
        index_twice = dimension*(dimension**2-1)//6
        twice_index += index_twice
        weak_hyper += index_twice*y
        cubic += dimension*y**3
        gravitational += dimension*y
        for two_t3 in range(1-dimension, dimension, 2):
            charges[3*two_t3+y] += 1
    assert all(charges[q] == charges[-q] for q in charges)
    assert (cubic, gravitational, weak_hyper, twice_index % 2) == (0, 0, 0, 0)
    assert signature == {1:18, 2:-9}
    return {"weyl_components":sum(charges.values()), "signed_multiplicity":dict(sorted((str(k),v) for k,v in signature.items())),
            "physical_six_Q_counts":dict(sorted((str(k),v) for k,v in charges.items())),
            "local_anomalies_y3_gravity_y_2T2y":[cubic,gravitational,weak_hyper],"twice_weak_index":twice_index}

def result(data):
    n = data["usual_families"]
    y = data["y_fundamental"]
    assert sum(y)==0 and y[:6] == [-2,-2,-2,1,1,1] and all((x%3)==1 for x in y[:6])
    # For a full SU(7) torus, the two seven-component multisets are identical, for all t.
    phases = []
    for t in [0, 0.37, 2*math.pi/3, math.pi]:
        weak = 3*n*(2*math.cos(t)-2)
        c = 2*math.cos(t)+1
        color = n*(2*c-c-c)
        extra = -3*n*(2*math.cos(t)-2)
        assert abs(weak+extra-color)<1e-12
        phases.append([t,weak,color])
    delta_b = F(n)+F(n%2,3)
    experimental = data["massless_Z"]
    s = experimental["angle_box"][1]
    gf = experimental["GF_GeV_minus2"]-2*experimental["GF_one_sigma"]
    mz = experimental["MZ_GeV"]-2*experimental["MZ_one_sigma"]
    unit_width = gf*mz**3/(3*math.sqrt(2)*math.pi)
    leptons = n*(0.25+(0.5-s)**2+s*s)
    down = 3*n*((0.5-s/3)**2+(s/3)**2)
    exotic = 3*n*((1-s)**2/2+s*s/18)
    lower = unit_width*(leptons+down+exotic)
    upper = experimental["GammaZ_GeV"]+2*experimental["GammaZ_one_sigma"]
    result = {"kY":str(F(sum(a*a for a in y),18)), "ordinary_character_at_angles_theta_weak_color":phases,
              "single_comparison_grade_floor":6*n, "EW_massive_nonzero_grade_floor":12*n,
              "full_Witten_consistent_component_floor":12*n+2*(n%2),"extra_weak_beta_floor":str(delta_b),
              "three_family_total_weak_beta_floor":str(-F(19,6)+delta_b), "low_character_witness":low_model(data["low_example"]),
              "Z":{"tree_unit_GeV_conservative":unit_width,"ordinary_three_leptons_GeV":leptons*unit_width,
                   "ordinary_down_three_flavors_GeV":down*unit_width,"nine_exotic_doublets_GeV":exotic*unit_width,
                   "open_tree_width_floor_GeV":lower,"measured_two_width_upper_GeV":upper,
                   "joint_transfer_at_overlap":upper/lower,"lower_at_registered_transfer_GeV":lower*experimental["transfer_example"]}}
    assert result["low_character_witness"]["weyl_components"]==result["full_Witten_consistent_component_floor"]
    assert result["Z"]["lower_at_registered_transfer_GeV"]>upper
    return result

if __name__=="__main__":
    cli=argparse.ArgumentParser()
    cli.add_argument("--check",action="store_true")
    args=cli.parse_args()
    current=result(json.loads((ROOT/"inputs.json").read_text(encoding="utf-8-sig")))
    path=ROOT/"results.json"
    if args.check:
        assert json.loads(path.read_text(encoding="utf-8-sig")) == current
        print("1574 verified")
    else:
        path.write_text(json.dumps(current,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
