"""Invert the charged-rate geometry with an independent rare neutral-current cap."""
import argparse
import json
from math import hypot, sqrt
from pathlib import Path

HERE=Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE/name).read_text(encoding="utf-8-sig"))


def edges(record):
    mean=float(record["mean"])
    error=float(record["sigma"]) if "sigma" in record else hypot(*map(float,record["sigma_components"]))
    half_width=error+float(record["rounding_half_unit"])
    return mean-half_width, mean+half_width


def maximum_light_row(x,y,r,P):
    # t=|Ld-Rd|. Beyond the larger unperturbed value, the measured weak cap
    # requires (t-x)(rt-y)<=4P, even with unrelated complex phases.
    t=(r*x+y+sqrt((r*x-y)**2+16*r*P))/(2*r)
    row_d=(x+t)/2
    row_s=(y+r*t)/2
    assert abs((row_d-x)*(row_s-y)-P)<2.e-16
    assert abs((row_s+(row_s-y))/(row_d+(row_d-x))-r)<1.e-12
    return row_d**2+row_s**2


def run():
    p=read("inputs.json")
    inputs=read(p["charged"])
    work=read("results.json")
    cap=work["joint_max_abs_RdRs"]
    x,y,rlo,rhi=edges(inputs["x"])[1],edges(inputs["y"])[1],*edges(inputs["r"])
    max_S=max(maximum_light_row(x,y,r,cap) for r in [rlo,rhi])
    ub=work["charged_windows"][1]["Vub_left_upper"]
    missing=sqrt(max(0.,1-max_S-ub**2))
    Bmax=work["charged_windows"][1]["VuB_left_upper"]
    assert missing>Bmax
    w=work["full_UV_example"]
    masses=p["full_witness"]
    # Projection of the common bare column on the light-u row is fixed directly by the physical masses and charges.
    left=w["V_left_u_row"]
    a=w["right_doublet_components_up"]
    b=w["right_doublet_components_down"]
    Bs=w["T_B_masses_GeV"][1]
    prediction=(masses["up_light_masses"][0]*a[0]-sum(q*i*j for q,i,j in
                zip(masses["down_light_masses"],left[:3],b[:3])))/(Bs*b[3])
    assert abs(left[3]-prediction)<1e-15
    return {"one_width_largest_light_ds_squared_at_rare_cap":max_S,
            "one_width_heavy_left_row_amplitude_required_at_least":missing,
            "heavy_row_author_mass_upper":Bmax,
            "stronger_than_allowed_by_factor":missing/Bmax,
            "central_uv_left_heavy_entry_from_common_column":prediction}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--write-new",action="store_true")
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=json.loads(json.dumps(run()),parse_float=lambda x:float(f"{float(x):.11g}"))
    if args.write_new:
        with (HERE/"independent_results.json").open("x",encoding="utf-8") as fp:
            json.dump(result,fp,indent=2)
            fp.write("\n")
    else:
        assert result==read("independent_results.json")
        print("1460 inverse charged-row/neutral bound checked")
