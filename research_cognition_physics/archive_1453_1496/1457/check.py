"""Inclusive SU(5) 15 branching fractions and two missing-flavour projections."""
import argparse
from itertools import product
import json
from math import sqrt
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE / name).read_text(encoding="utf-8-sig"))


def box(record):
    return record[0] - 2 * record[1], record[0] + 2 * record[-1]


def hierarchy(order, solar, reactor, atmospheric, gap21, gap31):
    # At m_lightest=0 each ratio is linear-fractional in each input separately.
    solar0 = []
    for a, b, s, t in product(solar, reactor, gap21, gap31):
        if order == "NO":
            solar0.append((s * a * (1-b) + t * b) / (s+t))
        else:
            solar0.append((t * (1-b) + s * a * (1-b)) / (2*t+s))
    if order == "IO":
        tau3 = (1-atmospheric[1]) * (1-reactor[1])
        lower = tau3 / (2 + gap21[1] / gap31[0])
    else:
        tau1 = max(0, sqrt(solar[0]*atmospheric[0])
                   - sqrt((1-solar[0])*(1-atmospheric[0])*reactor[1]))**2
        ratio = gap21[0] / gap31[1]
        lower = ratio * tau1 / (1 + ratio)**2
    lower = min(lower, 1/3)
    return {"electron_at_zero_lightest": [min(solar0), max(solar0)],
            "electron_all_mass_closure": [min(1/3, min(solar0)), max(1/3, max(solar0))],
            "light_jets_e_plus_mu_lower": lower,
            "same_flavour_pair_on_two_light_jets_lower": lower**2/2}


def explicit_example(order, seed):
    # Majorana phases affect individual tags. The unrestricted quark rotation can
    # align the omitted b jet exactly with the first or last output singular vector.
    rng = np.random.default_rng(seed)
    angles = (.306, .022, .515)
    s12, s13, s23 = np.sqrt(angles)
    c12, c13, c23 = np.sqrt(1-np.array(angles))
    phase = np.exp(.74j)
    mixing = np.array([[c12*c13,s12*c13,s13/phase],
                       [-s12*c23-c12*s23*s13*phase,c12*c23-s12*s23*s13*phase,s23*c13],
                       [s12*s23-c12*c23*s13*phase,-c12*s23-s12*c23*s13*phase,c23*c13]])
    inverted = order == "IO"
    squared = .0001 + np.array([.00248,.002554,0] if inverted else [0,.000074,.00251])
    mixing = mixing @ np.diag(np.exp(1j*np.array([0.,.6,.97])))
    tensor = mixing.conj() @ np.diag(np.sqrt(squared)) @ mixing.conj().T
    total = np.vdot(tensor,tensor).real
    qrotation = np.linalg.qr(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))[0]
    yields = abs(qrotation@tensor)**2/total
    double = np.triu(abs(tensor)**2*(2-np.eye(3))/total)
    delta_participation = (np.sum(double,axis=0)+np.sum(double,axis=1))/2
    assert np.allclose(delta_participation,yields.sum(axis=0))
    left, strengths, _ = np.linalg.svd(tensor[:, :2],full_matrices=True)
    actual = []
    for rows in ([1,2,0],[0,1,2]):
        selected = left.conj().T[rows] @ tensor
        actual.append(np.linalg.norm(selected[:2,:2])**2/total)
    assert np.allclose(actual,[strengths[1]**2/total, np.sum(strengths**2)/total])
    return {"inclusive_single_lepton": yields.sum(axis=0).tolist(),
            "sharp_light_jet_sum_range": actual}


def solve():
    cfg=read("inputs.json")
    osc=read(cfg["oscillation"])
    gaps=read(cfg["solar_and_mass_gaps"])["juno_207_days"]
    solar_no=read(cfg["solar_no"])["NO_solar_sin2"]
    reactor=tuple((1-sqrt(1-v))/2 for v in box(osc["DayaBay_amplitude"]))
    result={}
    for order in ["NO","IO"]:
        solar=box(solar_no if order=="NO" else osc[order]["solar_sine_squared"])
        result[order]=hierarchy(order,solar,reactor,box(osc[order]["atmospheric_sine_squared"]),
                                box(gaps["solar_21"]),
                                box(gaps["normal_atmospheric_31" if order=="NO" else "inverted_abs_31"]))
        result[order]["attainable_example"]=explicit_example(order,1457)
    assert result["IO"]["electron_at_zero_lightest"][0]>1/3
    assert result["NO"]["electron_at_zero_lightest"][1]<1/3
    atlas=cfg["ATLAS_139fb"]
    events={}
    for mass, light, bottom_by_tag in zip(atlas["reference_mass_TeV"],atlas["light_e_pretag"],
                                           atlas["bottom_e_012tag"]):
        bottom=sum(bottom_by_tag)
        response=1/(1/sqrt(light)+1/sqrt(bottom))**2
        events[str(mass)]={"at_least_one_template_mean_at_IO_all_mass_floor": response/9,
                            "at_IO_zero_lightest_e_lower":response*result["IO"]["electron_at_zero_lightest"][0]**2,
                            "equal_signal_light_fraction_within_e":sqrt(bottom)/(sqrt(light)+sqrt(bottom))}
    result["IO_e_ATLAS_nominal_expected_yield_floor"]=events
    return json.loads(json.dumps(result),parse_float=lambda x:round(float(x),9))


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    parser.add_argument("--write-new",action="store_true")
    args=parser.parse_args()
    result=solve()
    if args.write_new:
        with (HERE/"results.json").open("x",encoding="utf-8") as out:
            json.dump(result,out,ensure_ascii=False,indent=2)
            out.write("\n")
    else:
        assert result==read("results.json"),"saved results differ"
        print("1457 inclusive, partial visibility and ATLAS template checks passed")
