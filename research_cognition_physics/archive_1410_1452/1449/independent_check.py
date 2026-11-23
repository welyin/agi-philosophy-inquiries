"""Check the octet overlap and an independent second physical pseudoscalar pole."""

import argparse
import json
from math import hypot, pi, sqrt
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent


def results():
    cfg=json.loads((HERE/"inputs.json").read_text(encoding="utf8"))
    masses=cfg["mass_GeV"]
    # The extrema of I_s-I_d=2(m_K+^2-m_pi+^2) require only two inputs.
    xs_minus_xd=masses["kc"][0]**2-masses["pic"][1]**2-cfg["max_residual_I_GeV2"]
    assert xs_minus_xd>0
    # For any fixed u-d, increasing s-d adds a positive 88 projection. Test signs and extremes.
    v=np.array([1.,1/sqrt(3)])
    for diff in np.r_[-np.geomspace(.0001,1e7,50),0.,np.geomspace(.0001,1e7,50)]:
        n=diff*np.outer(v,v)+np.diag([0,4*xs_minus_xd/3])
        _,vectors=np.linalg.eigh(n)
        assert abs(vectors[1,-1])>=0.5-1e-12
    t=cfg["LO_transfer"]
    amp0=t["mpi_min_GeV"]*t["meta_min_GeV"]*0.5/t["F_max_GeV"]
    upper_width=cfg["eta_gammagamma_width_upper_GeV"]/cfg["eta_gammagamma_BR_lower"]
    def cap(br, width, parent):
        ml,mh=parent
        momentum=sqrt(1-(2*masses["pic"][1]/ml)**2)
        return sqrt(16*pi*mh*width*br/momentum)
    eta=cap(cfg["eta_to_pipi_BR_90"],upper_width,masses["eta"])
    # COSY-11: 0.226 +/- .017 +/- .014 MeV. Sum two widths conservatively.
    etap_width=(.226+2*(.017+.014))*.001
    etap=cap(5.5e-5,etap_width,[.95,.97])
    joint=hypot(eta,etap)
    # One channel can be exactly cancelled, even when the octet source is present.
    rotation=np.array([[.8,.6],[-.6,.8]])
    signal=np.array([amp0,-.8*amp0/.6])
    observed=rotation@signal
    assert abs(observed[0])<1e-14 and np.isclose(np.linalg.norm(observed),np.linalg.norm(signal))
    return {"lower_Xs_minus_Xd_GeV2":round(xs_minus_xd,10),
            "LO_octet_component_lower":0.5,"source_lower_GeV":round(amp0,10),
            "eta_amplitude_cap_GeV":float(format(eta,".11g")),
            "eta_prime_width_upper_GeV":etap_width,"eta_prime_amplitude_cap_GeV":round(etap,11),
            "both_poles_amplitude_norm_cap_GeV":round(joint,11),
            "critical_two_pole_retention":round(joint/amp0,10),
            "one_pole_cancellation_other_amplitude_GeV":round(abs(observed[1]),8)}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    flags=parser.add_mutually_exclusive_group()
    flags.add_argument("--write",action="store_true")
    flags.add_argument("--check",action="store_true")
    args=parser.parse_args()
    value=results()
    target=HERE/"independent_results.json"
    if args.write:
        with target.open("x",encoding="utf8") as stream:
            json.dump(value,stream,indent=2)
            stream.write("\n")
    else:
        assert value==json.loads(target.read_text(encoding="utf8"))
    print(json.dumps(value))
