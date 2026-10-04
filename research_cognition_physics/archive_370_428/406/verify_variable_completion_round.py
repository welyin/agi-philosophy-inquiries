"""Verify round 406 against frozen round 405; no image rendering."""
import argparse
import json
import verify_scattering_dictionary_round as previous

CONFIG = {406: ("variable_completion_time_audit", 7)}
BASES = {406: 405}
core = previous.core
core.CONFIG.update(CONFIG)
loaded = previous
while loaded is not None:
    if hasattr(loaded, "BASES"):
        loaded.BASES.update(BASES)
    loaded = getattr(loaded, "previous", None)


def verify(number, pending=False):
    result = previous.verify(number, pending)
    result.update(date="2026-09-23", parallel_batch=[406],
                  execution_mode="variable completion time exact implementation audit",
                  scientific_base_through_round=405,
                  additional_frozen_dependency_rounds=[], batch_scientific_dependencies=[],
                  final_science_review="primary-agent source, analytic-set proof and executable review; no independent agent final review")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("round", type=int, choices=CONFIG)
    parser.add_argument("--write-checks", action="store_true")
    args = parser.parse_args()
    result = verify(args.round, args.write_checks)
    if args.write_checks:
        with (core.HERE/f"research_round_{args.round}_checks.json").open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2)+"\n")
    print(json.dumps(result, ensure_ascii=False, indent=2))

