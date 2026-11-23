"""Round 1418. Run --check for a read-only replay or --write-new once."""
import argparse
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def vacuum_matrix(spectrum, s13sq, delta):
    s12, s13, s23 = np.sqrt([spectrum["s12sq"], s13sq, spectrum["s23sq"]])
    c12, c13, c23 = np.sqrt(1 - np.square([s12, s13, s23]))
    phase = np.exp(1j * delta)
    mixing = np.array([
        [c12*c13, s12*c13, s13*phase.conjugate()],
        [-s12*c23-c12*s23*s13*phase, c12*c23-s12*s23*s13*phase, s23*c13],
        [s12*s23-c12*c23*s13*phase, -c12*s23-s12*c23*s13*phase, c23*c13]
    ])
    assert np.allclose(mixing.conj().T @ mixing, np.eye(3))
    return (mixing * [0, spectrum["dm21"], spectrum["dm31"]]) @ mixing.conj().T


def probabilities(energies, mass, phase_distance, matter, factor, anti):
    spectral = mass.conjugate() if anti else mass
    matrix = phase_distance * spectral / energies[:, None, None]
    matrix[:, 0, 0] += phase_distance * matter * factor * (-1 if anti else 1)
    roots, vectors = np.linalg.eigh(matrix)
    amplitude = np.sum(vectors[:, 0, :] * np.exp(-1j*roots)
                       * np.conjugate(vectors[:, 1, :]), axis=1)
    return np.abs(amplitude)**2


def run(d):
    gp, gn = 0.5-2*d["sW2"], -0.5
    charges = {name: z*gp+n*gn for name, (z,n) in d["targets"].items()}
    scattering = {}
    for label in ["CsI_x10minus40_cm2", "Ar_x10minus39_cm2"]:
        scattering[label] = d[label]["measurement"][0] / d[label]["sm"][0]

    nu = d["NuTeV"]
    left, right = nu["G_L"], nu["G_R"]
    parity = left[0]-right[0]
    lower = parity-2*(left[1]+right[1])
    sm = nu["SM"][0]-nu["SM"][1]
    prop = d["propagation"]
    ye = prop["electron_fraction"]
    yn = (1-ye)/ye
    nuclear = gp+yn*gn
    factors = dict(SM=1, e_flip=1-2*nuclear, mutau_flip=1+2*nuclear)
    # L[km]/(2 E[GeV]) for an eV^2 mass squared; hbar*c is in GeV cm.
    phase_distance = prop["baseline_km"] * 1e-13/(2*prop["hbarc_GeV_cm"])
    matter = (2e18*np.sqrt(2)*prop["GF_GeV_minus2"]*prop["NA_mol_inverse"]*
              prop["mol_inverse_grams"]*prop["hbarc_GeV_cm"]**3*
              prop["density_g_cm3"]*ye)
    energies = np.linspace(*prop["energy_window_GeV"], prop["grid_intervals"]+1)
    middle = len(energies)//2
    s13sq = (1-np.sqrt(1-prop["theta13_double_angle_amplitude"]))/2
    alpha = prop["delta_q_rad"]
    phases = {"minus_a":-alpha, "plus_a":alpha, "pi_minus_a":np.pi-alpha,
              "pi_plus_a":np.pi+alpha}
    central = {}
    certificates = {}
    for ordering in ["NO", "IO"]:
        central[ordering] = {}
        spectrum = prop[ordering]
        span = np.ptp([0, spectrum["dm21"], spectrum["dm31"]])
        # For P1-P2 the 4 derivatives of amplitudes each cost L*span/(2*Emin²).
        lip = 2*phase_distance*span/energies[0]**2
        remainder = lip*(energies[-1]-energies[0])/(2*prop["grid_intervals"])
        for name, delta in phases.items():
            mass = vacuum_matrix(spectrum, s13sq, delta)
            central[ordering][name] = {}
            for anti, tag in [(False,"nu"), (True,"anti")]:
                pp = {k:probabilities(energies, mass, phase_distance, matter, val, anti)
                      for k,val in factors.items()}
                central[ordering][name][tag] = {k:arr[middle] for k,arr in pp.items()}
                if name == "minus_a" and ordering == "NO":
                    gap = pp["e_flip"]-pp["SM"]
                    certificates[tag] = dict(min=np.min(gap)-remainder,
                        max=np.max(gap)+remainder, grid_margin=remainder)
    assert certificates["nu"]["min"]>0 and certificates["anti"]["max"]<0
    return dict(tree_vector_charges=charges, measured_over_SM=scattering,
                nutev=dict(signed_difference=parity,two_marginal_widths_lower=lower,
                           flipped_SM_value=-sm, correction_to_reach_box=lower+sm),
                relative_potentials=factors, a_over_E_eV2_per_GeV=matter,
                central_probability_at_2GeV=central,
                NO_minus_alpha_signed_probability_gap_all_energies=certificates)


def verify(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:
            verify(a[k],b[k])
    else:
        assert np.isclose(a,b,rtol=1e-10,atol=1e-13), (a,b)


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check",action="store_true")
    mode.add_argument("--write-new",action="store_true")
    args=parser.parse_args()
    output=run(json.loads((HERE/"inputs.json").read_text(encoding="utf-8-sig")))
    target=HERE/"results.json"
    if args.write_new:
        with target.open("x",encoding="utf-8") as stream:
            json.dump(output,stream,ensure_ascii=False,indent=2)
            stream.write("\n")
    else:
        verify(output,json.loads(target.read_text(encoding="utf-8-sig")))
        print("1418 verified")
