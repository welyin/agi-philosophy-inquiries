"""Use the common bottom diagonal Higgs vertex as a separate readout."""
import argparse
import json
import math
from pathlib import Path

HERE=Path(__file__).resolve().parent


def load(path):
    return json.loads((HERE/path).read_text(encoding="utf-8-sig"))


def outer(row):
    return row[0]-2*sum(row[1:]),row[0]+2*sum(row[1:])


def run():
    masses=load("../1386/inputs.json")["yukawa"]
    blo,bhi=outer(masses["down"][2])
    shi=outer(masses["down"][1])[1]
    tlo=outer(masses["up"][2])[0]
    chi=outer(masses["up"][1])[1]
    cb=outer(load("../1412/inputs.json")["BelleII_Vcb_BGL"])[1]
    ubcb=outer(load("../1381/inputs.json")["ub_over_cb"])[1]
    k=cb**2*(1+ubcb**2)
    xbb_floor=tlo*(1-k)-chi*k
    xu_ct_cap=bhi*cb+shi*math.sqrt(k)
    coeff_bottom_min=xbb_floor-bhi
    assert coeff_bottom_min>0
    mtlo,mthi=170.,175.
    mhlo=124.7
    charm_kin_max=1.8
    bottom_kin_max=5.
    gamma_t_floor=.6
    beta_bb=math.sqrt(1-4*bottom_kin_max**2/mhlo**2)
    # Define Y by -h*Y/sqrt(2) at a single chiral vertex.
    # Positive bound using the limiting scalar diagonal rate even for arbitrary admixtures.
    direct_yhi=math.sqrt(.330*16*math.pi/(.8*3*mhlo*beta_bb**3))
    z_lineshape=(direct_yhi+bhi)/coeff_bottom_min
    meas=load("../1435/inputs.json")["atlas_table7"]["bb_ZZ"]
    ratio=(meas["observed"][0]+2*meas["observed"][2])/(meas["SM"][0]-2*meas["SM"][1])
    z_atlas=bhi*(1+math.sqrt(ratio/.8))/coeff_bottom_min
    rh_min=mhlo/mthi
    rc_max=charm_kin_max/mtlo
    kt=1.25*mthi/(32*math.pi)*(1-rh_min**2)*((1+rc_max)**2-rh_min**2)
    bound=lambda z: kt*(z*xu_ct_cap)**2/gamma_t_floor
    direct=load("../1535/inputs.json")["atlas_combined_tHc_BR_95"]
    return {"min_abs_bottom_opposite_doublet_coefficient":coeff_bottom_min,
            "max_abs_Xu_ct":xu_ct_cap,"width_ratio_b_to_Z_SM_normalized_upper":ratio,
            "max_abs_z_from_on_shell_width_and_bb_positivity":z_lineshape,
            "max_abs_z_over_radial_from_ATLAS_bbZZ":z_atlas,
            "max_B_t_to_ch_from_direct_Higgs_width_plus_diagonal":bound(z_lineshape),
            "max_B_t_to_ch_from_ATLAS_ratio_plus_diagonal":bound(z_atlas),
            "minimum_extra_rms_vertex_to_reach_direct_ATLAS_sensitivity":
                math.sqrt(direct*gamma_t_floor/kt)-z_atlas*xu_ct_cap}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--check",action="store_true")
    group.add_argument("--write-new",action="store_true")
    args=parser.parse_args()
    result=run()
    target=HERE/"independent_results.json"
    if args.write_new:
        with target.open("x",encoding="utf-8") as f:
            json.dump(result,f,indent=2)
            f.write("\n")
    else:
        assert result==json.loads(target.read_text())
    print(json.dumps(result))
