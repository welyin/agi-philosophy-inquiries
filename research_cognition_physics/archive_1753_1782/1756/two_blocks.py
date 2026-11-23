"""Universal screening floor with both up- and down-type compensators sharing one doublet."""
import argparse
import json
import math
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent


def solve():
    data=json.loads((HERE/"results.json").read_bytes())
    floor=data["min_mass_GeV"]
    low=data["low_reference_GeV"]
    strong_scale=data["U_GeV"]
    B=data["higgs_entry_cap_GeV"]
    # At least one of two mixed charge sectors supplies half the required negative amplitude.
    k=data["lighter_state_negative_Higgs_log_slope_required"]/2
    pair_cap=B**2*(1/k**2+2/k)
    doublet_cap=pair_cap/floor
    other_light_cap=math.hypot(doublet_cap,B)
    mean_cap=(B*B*pair_cap*other_light_cap)**.2
    clipped_cap=mean_cap*max(1.,low/floor)
    best_inverse=data["inverse_color_at_old_U"]-5*math.log(strong_scale/clipped_cap)/(3*math.pi)
    assert best_inverse<0 and other_light_cap<strong_scale
    # Check the second sector can decouple one singlet without removing its light pole.
    samples={}
    mq=400.
    for ms in (1000.,1e5,1e8,1e12):
        other=np.array([[mq,-225.],[190.,ms]])
        poles=np.linalg.svd(other,compute_uv=False)[::-1]
        assert poles[0]<=math.hypot(mq,225.)
        samples[str(ms)]=poles.tolist()
    return dict(minimum_negative_slope_one_of_two_sectors=k,compensating_pair_product_upper_GeV2=pair_cap,
                shared_doublet_hard_mass_upper_GeV=doublet_cap,other_charge_light_pole_upper_GeV=other_light_cap,
                upper_five_guaranteed_mass_geometric_mean_GeV=mean_cap,upper_effective_five_mean_GeV=clipped_cap,
                best_formal_inverse_alpha3_at_old_U=best_inverse,decoupling_examples_GeV=samples)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    args=p.parse_args()
    out=json.loads(json.dumps(solve()),parse_float=lambda x:float(f"{float(x):.10g}"))
    path=HERE/"two_block_results.json"
    if args.check:
        assert out==json.loads(path.read_bytes())
        print("1756 both-charge theorem: five compulsory low color thresholds verified")
    else:
        with path.open("x",encoding="utf8") as f:
            json.dump(out,f,indent=2,allow_nan=False)
            f.write("\n")
