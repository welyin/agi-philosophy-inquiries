"""Independent zero-temperature topology check at fixed renormalized masses."""
import json
import math
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "independent_results.json"

def energy_over_sigma(masses, theta):
    """Minimum of sum m_i (1-cos phi_i) on the small-theta positive branch."""
    lo, hi = 0., min(masses)
    for _ in range(65):
        trial = (lo+hi)/2
        if sum(math.asin(trial/m) for m in masses) > theta:
            hi = trial
        else:
            lo = trial
    return sum(2*m*math.sin(math.asin((lo+hi)/(2*m))/2)**2 for m in masses)

def calculate():
    lab = json.loads((ROOT.parent/"1480/inputs.json").read_text(encoding="utf-8"))["lattice"]
    masses = (2.130, lab["md"][0], lab["ms"][0])  # MILC 2 GeV, MeV
    inv = [1/m for m in masses]
    reduced = 1/sum(inv)
    b2 = -sum(i**3 for i in inv)*reduced**3/12
    def finite_b2(theta):
        return (2*energy_over_sigma(masses,theta)/(reduced*theta**2)-1)/theta**2
    richardson = (4*finite_b2(.04)-finite_b2(.08))/3
    assert abs(b2-richardson) < 1e-8
    fourth_root = 75.6  # Borsanyi et al., calibrated physical-point continuum estimate (MeV)
    quoted_stats, quoted_sys = 1.8,.9
    two_width = 2*(quoted_stats+quoted_sys)
    return dict(masses_MeV=masses, reduced_mass_MeV=reduced,
                b2_LO_fixed_positive_QCD_masses=b2, b2_checked_finite_theta=richardson,
                inferred_condensate_cube_root_MeV_at_LO=(fourth_root**4/reduced)**(1/3),
                physical_lattice_chi_fourth_root_MeV=fourth_root,
                positive_two_width_rectangle_chi_MeV4=[(fourth_root-two_width)**4,(fourth_root+two_width)**4],
                topology_when_one_exact_zero={"chi":0,"b2":"undefined (ratio of vanishing cumulants)"},
                assumptions="theta derivatives at fixed QCD mass source; LO condensate inferred, not independent. Independent four-flavour MS masses used as light ratios with common leading decoupling and scheme for the example.",
                sources=["https://arxiv.org/abs/1606.07494","https://arxiv.org/abs/1802.04248","https://arxiv.org/abs/1512.06746"])
if __name__ == "__main__":
    data = calculate()
    if sys.argv[1:]==["--write-new"]:
        with OUT.open("x",encoding="utf-8") as stream:
            json.dump(data,stream,indent=2)
            stream.write("\n")
    elif sys.argv[1:] in ([],["--check"]):
        assert json.loads(json.dumps(data)) == json.loads(OUT.read_text(encoding="utf-8"))
        print("1500: fixed-source susceptibility and independent finite-theta response verified.")
    else:
        raise SystemExit("Expected --check or --write-new.")
