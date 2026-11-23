"""Sharp overlap bound for three actual charged partners observed through one neutral window."""
import argparse
import json
from math import asin, hypot, sin, sqrt
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
PHASE = 1.26693268


def overlap_metrics(width_lower, data_lower, span, baseline):
    minimum = sqrt(width_lower-2)
    overlap = (1-minimum)/(1+minimum)
    distance = PHASE*span*baseline
    finite = sin(min(np.pi/2, asin(overlap)+distance))**2
    effective_overlap = sin(max(0., asin(sqrt(data_lower))-distance))
    required_width = 2+((1-effective_overlap)/(1+effective_overlap))**2
    wide_window = sin(asin(overlap)+PHASE*span*10.)**2
    # Equality at L=0 in the rotation with two mass axes at 45 degrees.
    rot = np.array([[1,-1,0],[1,1,0],[0,0,sqrt(2)]])/sqrt(2)
    n = rot @ np.diag([1.,sqrt(minimum),1.])
    gram = n@n.T
    assert np.isclose(np.trace((n.T@n)@(n.T@n)),width_lower)
    assert np.isclose(abs(gram[0,1])/sqrt(gram[0,0]*gram[1,1]),overlap)
    return {"width_lower":width_lower,"zero_distance_upper":overlap**2, "finite_window_upper":finite,
            "LSND_transfer_factor_needed":data_lower/finite,"width_max_to_allow_excess":required_width,
            "extra_invisible_units_or_width_bias_needed":width_lower-required_width,
            "critical_L_over_E":(asin(sqrt(data_lower))-asin(overlap))/(PHASE*span),
            "accepted_tail_fraction_needed_if_bounded_by_10":(data_lower-finite)/(wide_window-finite),
            "five_percent_tail_10_m_per_MeV_upper":.95*finite+.05*wide_window}


def run():
    c=json.loads((HERE/"inputs.json").read_text(encoding="utf-8"))
    nu=json.loads((HERE/c["sources"]["splittings"]).read_text(encoding="utf-8"))
    span=c["span_eVsquared"]
    assert all(nu[order]["atm_abs"][2]+nu[order]["solar"][2]<span for order in ("NO","IO"))
    w=c["matching_widths"]
    raw, stat, sys = c["lsnd_percent"]
    tail=(raw-w*hypot(stat,sys))/100
    boxes={name:overlap_metrics(val-w*err,tail,span,c["example_max_L_over_E_m_per_MeV"])
           for name,(val,err) in [(k,c[k]) for k in ("Z_reanalysis","Z_historical")]}
    assert all(row["finite_window_upper"]<tail for row in boxes.values())
    # One extra Dirac mass is below Z kinematics but above stopped muon decay.
    target=raw/100
    defect=2*sqrt(target)/(1+sqrt(target))
    x=(c["heavy_window_example_GeV"]/c["MZ_GeV"])**2
    mixed=(1-x)**2*(1+x/2)
    both=sqrt(1-4*x)*(1-x)
    light=3-2*defect+defect**2
    inv=light+2*defect*(1-defect)*mixed+defect**2*both
    electron_mass=(c["me_GeV"]/c["Mpi_GeV"])**2
    heavy=(c["heavy_window_example_GeV"]/c["Mpi_GeV"])**2
    phase_space=sqrt((1-electron_mass-heavy)**2-4*electron_mass*heavy)
    enhancement=phase_space*(electron_mass+heavy-(electron_mass-heavy)**2)/(electron_mass*(1-electron_mass)**2)
    assert np.isclose((defect/(2-defect))**2,target) and inv>boxes["Z_reanalysis"]["width_lower"]
    return {"round":1422,"LSND_two_width_lower":tail,"precision_comparison":boxes,"threshold_counterexample":{
        "heavy_active_fraction":defect,"electron_heavy_fraction":defect/2,"low_window_strength_squared_trace":light,
        "Z_width_if_all_escape":inv,"pion_e_heavy_over_e_light":defect/(2-defect)*enhancement,
        "pion_rest_positron_total_MeV":1000*c["Mpi_GeV"]*(1+electron_mass-heavy)/2},
        "independent_SM_parameters_calculated":0}


def near(a,b):
    if isinstance(a,dict):
        return a.keys()==b.keys() and all(near(a[k],b[k]) for k in a)
    return np.isclose(a,b,atol=2e-14,rtol=1e-9)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=run()
    if args.check:
        assert near(result,json.loads((HERE/"results.json").read_text(encoding="utf-8")))
        print("1422 common-window sharp overlap, finite propagation and independent pion line verified")
    else:
        print(json.dumps(result,indent=2))
