"""Two final-state mass derivatives, physical two-Higgs annihilation, and data slices."""
import argparse
import json
from math import pi, sqrt
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

HERE = Path(__file__).resolve().parent


def load(name):
    return json.loads((HERE / name).read_text(encoding="utf-8-sig"))


def derivatives(powers, weights):
    terms = np.array(weights, dtype=float)
    exponents = np.array(powers)
    output = [float(terms.sum())]
    for n in range(3):
        terms = terms * (exponents - n)
        output.append(float(terms.sum()))
    return output


def gamma_matrices():
    identity = np.eye(2)
    zero = np.zeros((2, 2))
    pauli = (np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.diag([1, -1]))
    return [np.diag([1, 1, -1, -1])] + [np.block([[zero, s], [-s, zero]]) for s in pauli]


GAMMA = gamma_matrices()
ID = np.eye(4)


def slash(p):
    return p[0] * GAMMA[0] - sum(p[j] * GAMMA[j] for j in (1, 2, 3))


def hh_rate(energy, mass, higgs, vev, kappa, eta, conversion, order):
    """Unpolarized colourless f fbar -> hh, exact Born s,t,u and contact."""
    energy_beam = energy / 2
    momentum = sqrt(energy_beam**2 - mass**2)
    outgoing = sqrt(energy_beam**2 - higgs**2)
    p1 = np.array([energy_beam, 0., 0., momentum])
    p2 = np.array([energy_beam, 0., 0., -momentum])
    coeff = mass / vev**2 * (eta + kappa * 3 * higgs**2 / (energy**2 - higgs**2))
    cosines, weights = leggauss(order)
    integral = 0.
    for cosine, weight in zip(cosines, weights):
        direction = outgoing * np.array([sqrt(1 - cosine**2), 0., cosine])
        vertex = coeff * ID
        for sign in (1, -1):
            internal = p1 - np.r_[energy_beam, sign * direction]
            denom = internal[0]**2 - internal[1:] @ internal[1:] - mass**2
            vertex = vertex + (kappa * mass / vev)**2 * (slash(internal) + mass * ID) / denom
        barred = GAMMA[0] @ vertex.conj().T @ GAMMA[0]
        spin_sum = np.trace((slash(p2) - mass * ID) @ vertex @ (slash(p1) + mass * ID) @ barred)
        assert abs(spin_sum.imag) < 1e-15 and spin_sum.real >= 0
        integral += weight * spin_sum.real / 4
    return float(conversion * integral * outgoing / momentum / (64 * pi * energy**2))


def run():
    data = load("inputs.json")
    shared = data["shared"]
    mu_cal = load(shared["muon_mass_and_higgs_mass"])
    muon = mu_cal["masses"]["muon_MeV"] / 1000
    higgs = mu_cal["cms"]["higgs_mass_GeV"]
    top = load(shared["top_alternative_pole_mass"])["alternative_top"][0]
    vev = (sqrt(2) * load(shared["fermi_constant"])["GF"])**-.5
    checks = {}
    for powers, weights, label in (([1],[1],"linear"), ([3],[1],"cubic"),
                                  ([3,5],[2,-1],"two_aligned"), ([3,5,7],[3,-3,1],"three_aligned")):
        checks[label] = derivatives(powers, weights)
    assert checks["two_aligned"] == [1,1,-8,-48]
    assert checks["three_aligned"] == [1,1,0,48]
    for p in range(3,23,2):
        for q in range(p+2,31,2):
            slope = p+q-1
            gap = (p-1)*(q-1)
            delta = slope * gap / (3+slope**2)
            weights = [(q-1-delta)/(q-p), (1+delta-p)/(q-p)]
            _,kappa,eta,_ = derivatives([p,q],weights)
            assert np.isclose(3*(kappa-1)**2+eta**2, 3*gap**2/(3+slope**2))
            assert 3*gap**2/(3+slope**2) >= 48/13-1e-13
    fixed = float(data["CMS_HIG_19_008_kt_observed_95_percent"][1][1])
    top_eta_bound = 7*(fixed-1)-8
    atlas = data["ATLAS_bbtautau_table4_ctthh_observed_95_percent"]
    assert -4 < atlas[0] <= 0 <= atlas[1]
    rates = {}
    for energy in data["muon_collider_prospective_energies_GeV"]:
        rates[str(energy)] = {}
        for label, (_,kappa,eta,_) in checks.items():
            args = (energy,muon,higgs,vev,kappa,eta,data["conversion_fb_per_GeV_minus2"])
            rate = hh_rate(*args, 192)
            assert np.isclose(rate, hh_rate(*args,384), rtol=2e-7, atol=1e-15)
            rates[str(energy)][label] = rate
    return {"vev_GeV":vev, "muon_mass_GeV":muon, "Higgs_mass_GeV":higgs, "top_example_GeV":top,
            "M_and_first_three_normalized_derivatives":checks,"sharp_joint_strength_squared":48/13,
            "at_sharp_minimum":{"kappa":27/13,"eta":-6/13},
            "two_powers_aligned_max_ctthh":-4, "ATLAS_fixed_slice":atlas,
            "max_top_eta_if_CMS_range_transferred":top_eta_bound,
            "aligned_leading_J0_half_crossing_GeV":{"muon":8*pi*vev**2/(muon*8),
                                                  "top":8*pi*vev**2/(sqrt(3)*top*8)},
            "unpolarized_muon_direct_hh_born_fb":rates}


def same(actual,saved):
    if isinstance(actual,dict):
        return actual.keys()==saved.keys() and all(same(actual[k],saved[k]) for k in actual)
    return np.allclose(actual,saved,rtol=2e-10,atol=1e-15)


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument("--write-new",action="store_true")
    mode.add_argument("--check",action="store_true")
    options=parser.parse_args()
    result=run()
    if options.write_new:
        with (HERE/"results.json").open("x",encoding="utf-8") as stream:
            json.dump(result,stream,indent=2)
            stream.write("\n")
    elif options.check:
        assert same(result,load("results.json"))
        print("1436 two-power envelope, Higgs contact slice and finite scattering verified")
    else:
        print(json.dumps(result,indent=2))
