"""Finite electroweak current tests of the seven-doublet running threshold."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

HERE = Path(__file__).resolve().parent
POINTS, WEIGHTS = leggauss(64)
X, WEIGHTS = (POINTS + 1)/2, WEIGHTS/2
U = X*(1-X)


def read(name):
    return json.loads((HERE/name).read_text(encoding="utf-8-sig"))


def secant(mass, q1, q2, alpha):
    # log1p preserves accuracy for the heavy singlet as well as the light poles.
    increment = np.log1p((q2-q1)*U/(mass**2+q1*U))
    return float(alpha/(4*math.pi)*np.dot(WEIGHTS, (1-2*X)**2*increment)/(q2-q1))


def pair_test(cap, source, cfg):
    energy = cfg["partonic_energy_GeV"]
    t3, charge = source
    def factor(s2):
        left, right = t3-charge*s2, -charge*s2
        return (left-right)**2/(2*(left**2+right**2))
    # Reciprocal is convex in s^2, hence its maximum is on an endpoint.
    k = min(factor(s) for s in cfg["sin2theta_scheme_range"])
    beta = math.sqrt(max(0., 1-(2*cap/energy)**2))
    base = math.pi*charge**2/(9*cfg["illustrative_inverse_alpha_hard"]**2*energy**2)
    return {"helicity_floor_fraction":k,"beta_floor":beta,
            "partonic_cross_section_floor_fb":cfg["GeV_minus2_to_fb"]*base*k*beta**3}


def mixed_matrix_check():
    masses2 = np.array([[10.,1.,3.], [1.,15.,-4.], [3.,-4.,8.]])
    masses, basis = np.linalg.eigh(masses2)
    principal = np.linalg.eigvalsh(masses2[:2,:2])
    assert min(masses)>0
    answer = {}
    for label, charge in {"T3":[.5,.5,0.],"hypercharge":[.5,.5,1.],"photon":[1.,1.,1.]}.items():
        generator = basis.T@np.diag(charge)@basis
        primitive = 0.
        spectral = 0.
        for i, a in enumerate(masses):
            for j, b in enumerate(masses):
                coefficient = abs(generator[i,j])**2
                kernel = np.dot(WEIGHTS,U**2/(X*a+(1-X)*b))
                primitive += coefficient*kernel
                edge = (math.sqrt(a)+math.sqrt(b))**2
                asymmetry = (math.sqrt(a)-math.sqrt(b))**2/edge
                # z=1-(s_threshold/s), followed by z=t^2 for a smooth endpoint.
                z = X**2
                dispersion = np.dot(WEIGHTS,2*X*z**1.5*(1-asymmetry*(1-z))**1.5)/(3*edge)
                spectral += coefficient*dispersion
        compressed = charge[0]**2*sum(1/principal)/30
        assert math.isclose(primitive,spectral,rel_tol=1e-12) and primitive>=compressed
        answer[label] = {"loop":primitive,"dispersion":spectral,"compressed_floor":compressed}
    return answer


def calculate():
    cfg = read("inputs.json")
    prior, gauges = read(cfg["prior"]), read(cfg["gauge"])
    experiment = cfg["small_angle_L3"]
    n, mw, alpha = cfg["n"], cfg["MW_GeV"], 1/read(cfg["em"])["electromagnetic_normalization"]["inverse_alpha_Thomson"]
    q1, q2 = experiment["transfer_Q2_GeV2"]
    assert math.sqrt(q1*q2)/4 < gauges["mu_ref_GeV"]**2
    sigma = math.hypot(*experiment["errors"])
    central = prior["central"]["seven_doublet_geomean_cap_GeV"]
    largest = prior["two_width_working_box"]["seven_doublet_geomean_cap_GeV"][1]
    benchmark = prior["seven_positive_witness"]["physical_charged_scalar_masses_GeV"][0]
    rows = {}
    for name, mass, widths, stretch in [("central",central,0,1.),
                ("outer",largest,2,1.),("central_200_GeV",benchmark,0,1.),
                ("outer_10pct_squared_mass_uplift",largest,2,cfg["test"]["physical_doublet_block_factor"])]:
        g2 = gauges["g2"][0]-widths*gauges["g2"][1]
        gy = math.sqrt(3/5)*(gauges["g1_GUT"][0]-widths*gauges["g1_GUT"][1])
        cap = math.sqrt(stretch)*mass
        norm = n*mw**2/(16*math.pi**2*60*cap**2)
        response = n*secant(cap,q1,q2,alpha)
        slope = n*alpha/(120*math.pi*cap**2)
        rows[name] = {"mass_scale_GeV":cap,"intact_W":norm*g2**2,"intact_Y_doublets":norm*gy**2,
                      "compressed_charged_W_floor":norm*g2**2/2,"compressed_charged_Y_floor":norm*gy**2/2,
                      "physical_alpha_slope_at_zero":slope,"finite_spacelike_secant":response,
                      "absolute_change_Delta_alpha":response*(q2-q1),"fraction_of_published_L3_error":response/sigma,
                      "max_first_pair_threshold_GeV":2*cap,
                      "up_pair":pair_test(cap,(.5,2/3),cfg["test"]),
                      "down_pair":pair_test(cap,(-.5,-1/3),cfg["test"])}
        assert response<=slope and response>=slope/(1+q2/(4*cap**2))
    return {"L3_positive_slope_and_hypot_error":[-experiment["S_GeV_minus2"],sigma],
            "L3_positive_slope_working_upper_at_two_widths":-experiment["S_GeV_minus2"]+2*sigma,
            "rows":rows,"mixed_normalized_toy":mixed_matrix_check()}


def matching(first,second):
    if isinstance(first,dict):
        return first.keys()==second.keys() and all(matching(v,second[k]) for k,v in first.items())
    if isinstance(first,list):
        return len(first)==len(second) and all(matching(i,j) for i,j in zip(first,second))
    return math.isclose(first,second,rel_tol=3e-11,abs_tol=2e-23)


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check",action="store_true")
    args = parser.parse_args()
    result = calculate()
    if args.check:
        if not matching(result,read("results.json")):
            raise SystemExit("stored current bounds differ")
        print("1593: weak and hypercharge moments, mixed spectral identity, L3 window and pair bounds verified")
    else:
        (HERE/"results.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
