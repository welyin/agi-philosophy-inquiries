"""Unit/sign check of the published low-z propagation contract; no data fit.

Default execution is read-only. --save creates the initial receipt exclusively.
This is not a likelihood, a cosmological propagation solver, or new science.
"""
from decimal import Decimal, localcontext
from pathlib import Path
import argparse
import json


def calculate():
    with localcontext() as ctx:
        ctx.prec = 65
        D = Decimal
        pi = D("3.1415926535897932384626433832795028841971693993751058209749445923")
        au = D(149597870700)
        c = D(299792458)
        pc = D(648000) * au / pi
        distance = D(26) * D(1000000) * pc
        factor = c / distance
        rows = []
        for source in (D(0), D(10)):
            lag = D("1.74")
            r = factor * (lag - source)
            delta = r / (1-r)
            reconstructed = source + (distance/c) * delta/(1+delta)
            assert abs(reconstructed-lag) < D("1e-60")
            assert abs(delta-r) < D("1e-28")
            assert (delta > 0) == (source == 0)
            rows.append({"source_delay_s": str(source), "r_linear": str(r),
                         "delta_constant_speed_exact": str(delta),
                         "absolute_linearization_error": str(abs(delta-r))})
        assert D("6.4e-16") < D(rows[0]["r_linear"]) < D("6.6e-16")
        assert D("-3.2e-15") < D(rows[1]["r_linear"]) < D("-3.0e-15")
        # Same observed lag also permits delta=0 and source lag=1.74 s.
        # This merely checks the dictionary; it is not a new degeneracy theorem.
        source_same_cone = D("1.74")
        assert factor * (D("1.74") - source_same_cone) == 0
        return {"all_checks_passed": True, "new_science_groups": 0,
                "new_cognitive_axioms": 0, "raw_data_fit": False,
                "new_confidence_interval": False,
                "distance_input_Mpc": 26,
                "distance_role": "published GW 90% credible lower bound, not exact distance",
                "parsec_unit_convention": "648000 AU / pi; AU=149597870700 m",
                "c_over_distance_per_s": str(factor), "central_value_checks": rows,
                "same_cone_source_delay_s": str(source_same_cone),
                "published_rounded_interval": ["-3e-15", "7e-16"],
                "scope": "historical low-redshift effective-distance approximation"}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--save", action="store_true")
    args = parser.parse_args()
    result = calculate()
    path = Path(__file__).with_name("results.json")
    if args.save:
        with path.open("x", encoding="utf-8") as out:
            json.dump(result, out, ensure_ascii=False, indent=2)
            out.write("\n")
    else:
        assert result == json.loads(path.read_text(encoding="utf-8"))
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
