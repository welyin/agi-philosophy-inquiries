"""Independent check of the threshold displacement and its observable cost."""
import argparse
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
stage = HERE.parent
weak = json.loads((stage/"1621"/"inputs.json").read_text(encoding="utf8"))
poles = json.loads((stage/"1619"/"inputs.json").read_text(encoding="utf8"))
mtau = json.loads((stage/"1619"/"results.json").read_text(encoding="utf8"))["roots_MeV"][1]


def crossing(g_y, g_2):
    return weak["mu0"]*math.exp(8*math.pi**2*(1/g_y**2-3/g_2**2)/(49/3))


def inverse_alpha_difference(m, panels):
    qlo, qhi = 12.25, 3434.
    terms = []
    for j in range(panels+1):
        x = j / panels
        k = x*(1-x)
        integrand = k*math.log1p((qhi-qlo)*k/(m*m+qlo*k))
        terms.append((1 if j in (0,panels) else 4 if j%2 else 2)*integrand)
    return 2/math.pi*math.fsum(terms)/(3*panels)


def photon_pb(mass):
    s = 206.6**2
    fine = 1/poles["alpha_inverse"]
    beta = math.sqrt(1-4*mass*mass/s)
    return 4*math.pi*fine*fine/(3*s)*beta*(1+2*mass*mass/s)*3.8937937e8


def calculate():
    m_e, m_mu = [poles["pole_MeV"][k][0] for k in ("electron", "muon")]
    weights = [math.sqrt(m/mtau) for m in (m_e,m_mu,mtau)]
    g_y, g_2 = weak["box"]["gY"], weak["box"]["g2"]
    scale_lo, scale_hi = crossing(g_y[1],g_2[0]), crossing(g_y[0],g_2[1])
    at_upper = [scale_hi*w for w in weights]
    qed = 2*sum(inverse_alpha_difference(m,512) for m in at_upper)
    assert abs(qed-2*sum(inverse_alpha_difference(m,1024) for m in at_upper))<1e-10
    zfactor = 1/(1+16*(.22-.25)**2)
    floor = zfactor*photon_pb(at_upper[0])
    beta_doublet, beta_singlet = 4/3*(2*(1/2)**2-3*(1/2)), 4/3
    assert abs(beta_singlet+beta_doublet)<1e-12
    return {"SU2_instanton_winding_adjacent": [-2,1,4],
            "delta_bY_minus_3b2_D_E":[beta_doublet,beta_singlet],
            "SM_crossing_GeV":[scale_lo,crossing(weak["central"]["gY"],weak["central"]["g2"]),scale_hi],
            "tau_to_e_flavon":1/weights[0], "paired_max_mass_vector_GeV":at_upper,
            "absolute_min_mass_ceiling_GeV":at_upper[0],
            "min_two_tower_extra_inverse_alpha_difference":qed,
            "fraction_of_L3_error":qed/1.29,
            "OPAL_qed_unit_charge_pb":photon_pb(at_upper[0]),"possible_Z_factor_min":zfactor,
            "OPAL_signal_min_pb":floor,"minimum_transfer_to_0p028":.028/floor,
            "nonuniversal_yukawa_condition_number_needed_for_102_GeV_paired":102/at_upper[0],
            "equal_500GeV_counterexample_delta_inverse_alpha":6*inverse_alpha_difference(500,512),
            "exception_120GeV_lightest_heaviest_GeV":120/weights[0]}


def normal(x):
    if isinstance(x, dict):
        return {k:normal(v) for k,v in x.items()}
    if isinstance(x, list):
        return [normal(v) for v in x]
    return float(f"{x:.11g}") if isinstance(x,float) else x


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    options=parser.parse_args()
    assert not (options.write and options.check)
    data=normal(calculate())
    path=HERE/"independent_results.json"
    if options.write:
        with path.open("x",encoding="utf8") as output:
            json.dump(data,output,indent=2)
            output.write("\n")
    else:
        assert data==json.loads(path.read_text(encoding="utf8"))
        print("PASS 1622 independent completion and charged spectrum")
