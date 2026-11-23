"""Electrical and magnetic charge lattices in units of e and 2*pi/e."""
import argparse
import json
from fractions import Fraction
from math import gcd
from pathlib import Path

HERE = Path(__file__).resolve().parent


def read_json(name):
    return json.loads((HERE / name).read_text(encoding="utf-8-sig"))


def distance(value, step):
    remainder = value % step
    return min(remainder, step - remainder)


def build():
    data = read_json("inputs.json")
    forms = data["global_forms"]
    k = data["positive_test"]["electric_k"]
    n = data["positive_test"]["magnetic_n"]
    q = Fraction(k, 6)
    margins = data["positive_test"]["worked_error_box"]
    electric_error = Fraction(margins["q"])
    magnetic_error = Fraction(margins["n"])
    competition = {}
    for p in forms:
        electric_step = Fraction(p, 6)
        magnetic_step = 6 // p
        pure_magnetic_step = magnetic_step * 3 // gcd(magnetic_step, 3)
        dq = distance(q, electric_step)
        dn = distance(n, magnetic_step)
        competition[str(p)] = {
            "colorless_electric_step": str(electric_step),
            "magnetic_step": magnetic_step,
            "zero_color_center_magnetic_step": pure_magnetic_step,
            "distance_to_electric_lattice": str(dq),
            "distance_to_magnetic_lattice": dn,
            "within_demo_box": dq <= electric_error and dn <= magnetic_error,
        }
        for weight in data["weights"]:
            assert (weight["k"] + 2 * weight["triality"]) % p == 0

    # Independently enumerate actual centre endpoints rather than using divisibility.
    for p in forms:
        allowed = {(j % 3, j % 2, j % 6) for j in range(0, 6, 6 // p)}
        for flux in range(-24, 25):
            endpoint = (flux % 3, flux % 2, flux % 6)
            assert (endpoint in allowed) == ((p * flux) % 6 == 0)
    for sample_k in range(-12, 13):
        for sample_n in range(-12, 13):
            lower = 6 // gcd(6, sample_n)
            upper = gcd(6, sample_k)
            expected = [p for p in forms if p % lower == 0 and upper % p == 0]
            enumerated = [p for p in forms if sample_k % p == 0 and sample_n % (6 // p) == 0]
            assert expected == enumerated

    electric_forms = [p for p in forms if k % p == 0]
    magnetic_forms = [p for p in forms if (n * p) % 6 == 0]
    selected = sorted(set(electric_forms).intersection(magnetic_forms))
    assert selected == [2]
    hidden = data["hidden_test"]
    induced_charge = Fraction(hidden["charge"])
    hidden_examples = {}
    for winding in hidden["raw_hidden_fluxes"]:
        dark_flux = winding - n * induced_charge
        pairing = n * induced_charge + dark_flux
        assert pairing == winding
        hidden_examples[str(winding)] = {"visible_flux": n, "hidden_flux": str(dark_flux),
                                         "full_pairing": str(pairing)}
    cms = data["experiments"]["CMS"]["electric_mass_intervals_GeV"][str(k)]
    moedal = data["experiments"]["MoEDAL"]["table_I_mass_limits_GeV"]
    benchmark = data["benchmarks_GeV"]
    return {
        "global_forms_from_electric_alone": electric_forms,
        "global_forms_from_magnetic_alone": magnetic_forms,
        "both_positive_readings": selected,
        "per_form": competition,
        "hidden_U1_at_visible_p6": hidden_examples,
        "extra_Coulomb_Gram_readout_in_demo": str(Fraction(hidden["g_hidden_over_e"]) ** 2),
        "finite_benchmark_in_published_mass_ranges_if_matching": {
            "electric": cms[0] <= benchmark["electric"] <= cms[1],
            "magnetic": {process: benchmark["magnetic"] < boundary
                         for process, boundary in moedal.items()},
        },
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--check", action="store_true")
    group.add_argument("--write-new", action="store_true")
    args = parser.parse_args()
    result = build()
    if args.write_new:
        with (HERE / "results.json").open("x", encoding="utf-8") as output:
            json.dump(result, output, ensure_ascii=False, indent=2)
            output.write("\n")
    else:
        assert result == read_json("results.json")
        print("1471: centre endpoints, finite charge margins, experimental windows and hidden lattice OK")
