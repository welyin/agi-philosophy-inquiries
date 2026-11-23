"""Test the pure-ten proton-decay operator on the actual passing lepton mass rows."""
import argparse
import json
from pathlib import Path
import runpy
import numpy as np

HERE=Path(__file__).resolve().parent
PORT=runpy.run_path(str(HERE.parent/"1748"/"check.py"))
GAUGE=runpy.run_path(str(HERE.parent/"1744"/"independent_check.py"))


def solve():
    b=PORT["OLD"]["mother"](.125)
    _,f,_,_=PORT["terms"](b,0,(1,-1,1))
    gain,_=PORT["daughter"](b,f)
    split=GAUGE["gauge_and_wilson"](b["Mweak"])
    u=split["cross_GeV"]
    cu=4*np.pi/(split["a_at_cross"][0]*u**2)
    # The scalar calibration applies only to the charged source above the gate.
    calibrated=dict(b)
    calibrated["yW"]=b["yW"].copy()
    calibrated["yW"][2]*=gain
    low=GAUGE["at_low"](calibrated,f)
    pt=b["pt"]
    left=b["left_quark_to_GJ"]
    hi=b["high"]
    re=np.linalg.solve(hi[2,:,:3],pt["E"].T)
    ru=np.linalg.solve(hi[0,:,:3],left*np.linalg.svd(hi[0,:,:3],compute_uv=False)[::-1])
    v=np.array([1.,0.,0.])
    d0=np.linalg.svd(pt["D"])[0][:,-1]
    e0=np.linalg.svd(pt["E"])[0][:,-1]
    dl=np.linalg.svd(left.conj().T@low[1,:,:3])[0][:,-1]
    el=np.linalg.svd((low[2,:,:3]@re).T)[0][:,-1]
    ul,_,ur=np.linalg.svd(left.conj().T@low[0,:,:3]@ru)
    chords=[GAUGE["projective_chord"](x,y) for x,y in ((d0,dl),(e0,el),(v,ul[:,-1]),(v,ur[-1,:]))]
    coeff=abs(e0*(d0+v*np.vdot(v,d0)))
    floor=2*max(coeff)-sum(coeff)-2*sum(chords)
    # A complementary direct coefficient test before taking the uniform phase envelope.
    up,uc=ul[:,-1],ur[-1,:]
    t=np.vdot(uc,up)*dl+np.vdot(uc,dl)*up
    direct=abs(el.conj()*t)
    direct_floor=2*max(direct)-sum(direct)
    assert direct_floor>=floor>1.9
    hadron=runpy.run_path(str(HERE.parent/"1725"/"independent_check.py"))
    out={}
    for label,side in (("central",0),("hadronic_rectangular_floor",1)):
        limit=hadron["coefficients"]("pi_e",side)
        absolute=limit/(split["down_to_2GeV_gain"]*floor)
        ratio=cu/absolute
        out[label]=dict(max_matching_coefficient_GeV_m2=absolute,
                        rate_to_SuperK_single_mode_limit_at_least=ratio**2,
                        direct_coefficient_fraction_max=1/ratio,
                        minimal_amplitude_correction_in_units_of_original_floor=1-1/ratio)
    return dict(calibrated_tau_gain=gain,U_GeV=u,direct_matching_GeV_m2=cu,
                wilson_gauge_gain_to_2GeV=split["down_to_2GeV_gain"],quark_lepton_external_chords=chords,
                triangle_flavour_lower=floor,explicit_rephasing_envelope=direct_floor,limits=out)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    args=p.parse_args()
    obj=json.loads(json.dumps(solve()),parse_float=lambda x:float(f"{float(x):.10g}"))
    path=HERE/"proton_results.json"
    if args.check:
        assert obj==json.loads(path.read_bytes())
        print("1750 successful port: proton matching coefficient and uniform finite-electron rate bound verified")
    else:
        with path.open("x",encoding="utf8") as f:
            json.dump(obj,f,indent=2,allow_nan=False)
            f.write("\n")
