"""Independent two-body trace and the e/mu centroid of the unresolved pair."""
import argparse
import json
from math import sqrt
from pathlib import Path
ROOT = Path(__file__).resolve().parent

def read(path):
    return json.loads((ROOT/path).read_text(encoding="utf-8-sig"))

def interval(row):
    return row[0]-2*row[1], row[0]+2*row[-1]

def run():
    data = read("inputs.json")
    oscill = read(data["oscillation_reactor_and_accelerator"])
    juno = read(data["oscillation_JUNO"])["juno_207_days"]
    M = data["reference_bare_heavy_mass_MeV"]*1e6
    m2min,m2max = map(sqrt,interval(juno["solar_21"]))
    m3min,m3max = map(sqrt,interval(juno["normal_atmospheric_31"]))
    a,A = interval(data["JUNO_sin2_theta12_207d"])
    b,B = [(1-sqrt(1-x))/2 for x in interval(oscill["DayaBay_amplitude"])]
    c,C = interval(oscill["NO"]["atmospheric_sine_squared"])
    weak_ratio = (sqrt(2)-1)**2
    xl = weak_ratio*(1-4*m3max/M)
    xh = weak_ratio*(1+4*m3max/M)
    # Ratios of efficiencies times two-body weights for individual unresolved states.
    # Each is assumed within +/-10% of the channel's common reference.
    low_ratio, high_ratio = .9/1.1, 1.1/.9
    emin, emax = a*(1-B), A*(1-b)
    e3max = 1/(1+xl*emin/B*low_ratio)
    e3min = 1/(1+xh*emax/b*high_ratio)
    mu2max = (sqrt((1-a)*(1-c)) + sqrt(A*C*B))**2
    mu2min = max(0., sqrt((1-A)*(1-C))-sqrt(A*C*B))**2
    mu3min = 1/(1+xh*mu2max/(c*(1-B))*high_ratio)
    mu3max = 1/(1+xl*mu2min/(C*(1-b))*low_ratio)
    dl = m3min-m2max
    dh = m3max-m2min
    shifts = [(mu3min-e3max)*dl,(mu3max-e3min)*dh]
    assert 0 < shifts[0] < shifts[1] < dh
    raw = []
    mk = data["kaon_mass_MeV_reference"]
    for lepton in ["electron","muon"]:
        m = data[lepton+"_mass_MeV_reference"]
        def width(neutral):
            p = sqrt((mk**2-(m+neutral)**2)*(mk**2-(m-neutral)**2))/(2*mk)
            # Spin contraction after the pseudoscalar current is replaced by its momentum.
            amp = (m*m+neutral**2)*mk**2-(m*m-neutral**2)**2
            return p*amp
        raw.append(width(M/1e6)/width(0))
    results = {"two_body_rate_enhancement_from_rest_momenta": raw,
            "heavy_state3_fraction_e":[e3min,e3max],
            "heavy_state3_fraction_mu":[mu3min,mu3max],
            "mean_missing_mass_mu_minus_e_eV":shifts,
            "actual_e_bin_excess":236745-234678,
            "source_local_significance_approx":(236745-234678)/sqrt(236745+314**2)}
    return {k:[round(t,10) for t in v] if isinstance(v,list) else round(v,10) for k,v in results.items()}

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    group=parser.add_mutually_exclusive_group()
    group.add_argument("--write-new",action="store_true")
    group.add_argument("--check",action="store_true")
    opt=parser.parse_args()
    result=run()
    if opt.write_new:
        with (ROOT/"independent_results.json").open("x",encoding="utf8") as f:
            json.dump(result,f,indent=2)
            f.write("\n")
    elif opt.check:
        assert result==read("independent_results.json")
        print("1448 independent kaon kinematics and unresolved centroid verified")
    else:
        print(json.dumps(result,indent=2))
