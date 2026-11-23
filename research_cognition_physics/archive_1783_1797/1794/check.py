"""Neutral-Higgs flavour completeness, spacelike response, and finite mass windows."""

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def gramians():
    """Independent weak-basis construction with four doublets and three real singlets."""
    rng = np.random.default_rng(1794)
    smallest, maximum_defect = math.inf, 0.0
    for _ in range(48):
        # Three ordinary full complex flavour matrices; combine two selected LFV ports.
        flavour = (rng.normal(size=(3, 3, 3)) + 1j*rng.normal(size=(3, 3, 3))) / 4
        neutral_basis = np.zeros((10, 9), dtype=complex)
        neutral_basis[0] = np.diag([.001, .01, .05]).ravel()
        neutral_basis[1:4] = flavour.reshape(3,9)
        neutral_basis[4:7] = 1j*flavour.reshape(3,9)
        mixing = np.linalg.qr(rng.normal(size=(10,10)))[0]
        g = mixing @ neutral_basis
        c = mixing[:,0]
        mass = np.r_[125., rng.uniform(150,1000,size=9)]
        m = neutral_basis[0]
        residual = m-c[0]*g[0]
        np.testing.assert_allclose(c@g,m,atol=2e-14)
        np.testing.assert_allclose(g.T@g,np.outer(m,m),atol=2e-14)
        for q in (0,100,500):
            denominator = mass[1:]**2+q**2
            # Columns in channel space use conjugate conjugation on their second index.
            response = (g[1:].T/denominator) @ g[1:].conjugate()
            moment = c[1:]**2 @ denominator
            bound = np.outer(residual,np.conj(residual))/moment
            smallest = min(smallest, np.linalg.eigvalsh(moment*(response-bound)).min())
        for j in (1, 5):
            phase = g[0,j].conjugate()/abs(g[0,j])
            x, y = np.real(g[:,j]*phase), np.imag(g[:,j]*phase)
            normx, normy = x[1:] @ x[1:], y[1:] @ y[1:]
            assert normx >= c[0]**2*abs(g[0,j])**2/(1-c[0]**2)-1e-13
            assert normy >= abs(g[0,j])**2/(1-c[0]**2)-1e-13
            maximum_defect = max(maximum_defect,abs(x[1:]@y[1:]),abs(normy-normx-x[0]**2))
    assert smallest > -3e-12 and maximum_defect < 2e-13
    return {"random_mixtures":48,"min_eigenvalue_of_scaled_gram_difference":round(float(smallest),12),
            "max_orthogonal_partner_identity_defect":round(float(maximum_defect),12)}


def example(inputs, branching):
    v, mh, c = (inputs[k] for k in ("v_GeV","h_mass_GeV","example_kappa_V"))
    a = 1-c*c
    # The massless two-charge LFV width is mh/(8pi) times this sum; finite-mass
    # corrections for arbitrary chiral phases only decrease it at fixed squared norm.
    signal = 8*math.pi*inputs["assumed_h_total_width_GeV"]*branching/mh
    budget = 4*math.pi*v**2 - c**2*mh**2
    residual = c*math.sqrt(signal)
    qpoints = inputs["spacelike_momenta_GeV"]
    extras = [residual**2/(budget+a*q*q) for q in qpoints]
    floors = [e+signal/(mh**2+q*q) for q,e in zip(qpoints,extras)]
    # If there is a uniform bound on the complete positive high-energy chiral
    # scattering matrix, one coherent flavour direction accounts for both hands.
    trial_norms = inputs["finite_window_portal_norm_examples"]
    windows = {}
    for cutoff2, label in [(8*math.pi, "scalar_Yukawa_elastic_8pi"),
                           *((y*y, "example_total_norm_"+str(y)) for y in trial_norms)]:
        if cutoff2 < 2*signal/a:
            windows[label] = {"incompatible_with_doublet_Gram_identity":True}
        else:
            real_budget = cutoff2/2 - signal
            crossing = math.sqrt(budget*real_budget)/residual
            finite_windows = {}
            for edge in (3000.,5000.,10000.):
                needed = max(0.,residual-math.sqrt(budget*real_budget)/edge)
                finite_windows[str(int(edge))] = {
                    "min_fraction_of_signed_residual_in_window":needed/residual,
                    "min_squared_flavour_coupling_in_window":needed**2/a}
            windows[label] = {"impossible_strictly_all_above_GeV":crossing,
                              "finite_windows_GeV":finite_windows}
    return {"branching_fraction_scenario":branching, "h_ordered_LFV_strength_massless_or_lower":signal,
            "gauge_weight_remaining":a, "gauge_mass_moment_remaining_upper_GeV2":budget,
            "two_polarizations_strength_min_per_h_squared":{"real":c*c/a,"imag":1/a},
            "spacelike_q_GeV":qpoints,"other_states_exchange_floor_GeV_minus2":extras,
            "including_h_exchange_floor_GeV_minus2":floors,"mass_windows":windows}


def cancellation(inputs, signal):
    v, mh, mH, c = (inputs[k] for k in ("v_GeV","h_mass_GeV","heavier_CP_even_mass_GeV","example_kappa_V"))
    a = 1-c*c
    mA = 1/math.sqrt(a/mh**2 + c*c/mH**2)
    masses = np.array([mh,mH,mA])
    g = math.sqrt(signal/a)*np.array([math.sqrt(a), c, 1j])
    gauge = np.array([c, -math.sqrt(a), 0])
    assert abs(g@g) < 1e-15 and abs(gauge@g) < 1e-15
    neutral_moment = np.sum((gauge*masses)**2)
    assert neutral_moment < 4*math.pi*v**2
    rows = []
    for q in inputs["spacelike_momenta_GeV"]:
        positive = sum(abs(g)**2/(masses**2+q*q))
        holomorphic = sum(g**2/(masses**2+q*q))
        rows.append({"q_GeV":q,"positive_GeV_minus2":positive, "same_flavour_change_GeV_minus2":holomorphic.real})
    assert abs(rows[0]["same_flavour_change_GeV_minus2"]) < 1e-24
    return {"masses_h_H_A_GeV":masses.tolist(), "gauge_mass_moment_GeV2":neutral_moment,
            "Yukawa_real_imag":[[z.real,z.imag] for z in g],"responses":rows}


def run():
    inputs = json.loads((HERE/"inputs.json").read_bytes())
    emu = example(inputs,inputs["emu_branching_at_CMS_observed_95pct_limit"])
    mutau = example(inputs,inputs["mutau_branching_at_ATLAS_two_marginal_width_lower_edge"])
    return {"gram_checks":gramians(),"emu_conditional_sensitivity":emu,"mutau_two_width_condition":mutau,
            "one_chiral_emu_countermodel":cancellation(inputs,emu["h_ordered_LFV_strength_massless_or_lower"])}


def normalize(value):
    if isinstance(value,dict):
        return {k:normalize(v) for k,v in value.items()}
    if isinstance(value,list):
        return [normalize(v) for v in value]
    if isinstance(value,(float,np.floating)):
        return float(f"{value:.10g}")
    return value


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    output=json.dumps(normalize(run()),ensure_ascii=False,indent=2,allow_nan=False)+"\n"
    file=HERE/"results.json"
    if args.check:
        assert json.loads(file.read_bytes())==json.loads(output)
    else:
        file.write_text(output,encoding="utf-8")
    print("1794 neutral flavour partners",hashlib.sha256(output.encode()).hexdigest()[:12])
