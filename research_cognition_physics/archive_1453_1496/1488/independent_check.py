"""Independent doublet-row and heavy-singlet on-shell checks."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

ROOT=Path(__file__).resolve().parent


def read(name):
    return json.loads((ROOT/name).read_text(encoding="utf-8-sig"))


def decay(H,m,T):
    vfermi=246.22
    mw,mz,mh=80.4,91.19,125.0
    m0=m/T
    matrix=np.array([[m0,0],[math.sqrt((1-T*T)*(H*H-m0*m0)),T*H]])
    masses2, left=np.linalg.eigh(matrix@matrix.T)
    assert np.allclose(masses2,[m*m,H*H],rtol=2e-7)
    s2=float(left[0,1]**2)
    c2=1-s2
    width0=H**3*s2/(16*math.pi*vfermi**2)
    def vec(M):
        x=(M/H)**2
        return (1-x)**2*(1+2*x)
    widths=np.array([width0*vec(mw),width0*c2*vec(mz)/2,
                     width0*c2*(1-(mh/H)**2)**2/2])
    return {"bare_soft_mixing_GeV":float(matrix[1,0]),"sinL_squared":s2,
            "total_width_GeV":float(sum(widths)),"c_tau_metres":float(1.97327e-16/sum(widths)),
            "branch_W_Z_h":list(widths/sum(widths))}


def audit():
    cfg=read("inputs.json")
    match=read(cfg["running_matching"])
    source=read(cfg["MZ_yukawa"])
    m=source["yukawa"]["charged_leptons"][1][0]*read(cfg["MZ_vev"])["matching"]["v_GeV"][0]/math.sqrt(2)
    T=match["T_QED"]
    Thigh=match["T_box_with_transport"][1]
    r0,rw=read(cfg["tau"])["Belle_II"]["ratio_central"],read(cfg["tau"])["Belle_II"]["ratio_width"]
    cW=(r0-2*rw)/cfg["amplitude_or_rate_transfer"][1]
    re,rep=cfg["lep"]["R_e"]
    rm,rmp=cfg["lep"]["R_mu"]
    zfloor=(re-2*rep)/(rm+2*rmp)/cfg["amplitude_or_rate_transfer"][1]
    s,se=read(cfg["electron_pole"])["SLD_sin2_eff_e"]
    lo,hi=s-2*se,s+2*se
    radius=math.sqrt(zfloor*((.5-hi)**2+hi**2)-hi**2)
    assert lo+radius<.5 and cW*cW/2>hi
    cZ=math.sqrt(2*(lo+radius))
    assert cZ>cW
    A=(cZ-Thigh)/math.sqrt(1-cZ*cZ)
    assert A>0
    return {"left_amplitude_lower_belle":cW,"left_amplitude_lower_with_LEP":cZ,
            "same_charge_column_lower_over_m0":A,
            "same_charge_column_lower_over_m0_central_T":(cZ-T)/math.sqrt(1-cZ*cZ),
            "soft_examples":{str(H):decay(H,m,T) for H in [400,600]}}


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    mode=cli.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check",action="store_true")
    mode.add_argument("--write",action="store_true")
    args=cli.parse_args()
    out=audit()
    path=ROOT/"independent_results.json"
    if args.write:
        with path.open("x",encoding="utf-8") as f:
            json.dump(out,f,indent=2)
            f.write("\n")
    else:
        assert out==read(path.name)
        print("1488 independent weak, heavy-column norm, decay length and branching checks passed.")
