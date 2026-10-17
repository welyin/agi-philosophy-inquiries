"""Read-only finite arithmetic for a joint-menu admission screen, not a research round."""
from fractions import Fraction as F
from pathlib import Path
import argparse
import json


def calculate():
    ell, width, eta = F(1, 16), F(1, 4), F(1, 5)
    A_wedge_upper = 1 / (2 * ell)
    # beta/T < eta, beta=2*pi/A, 2T<=width.
    coefficient_pi = 4 / (width * eta)
    A_thermal_lower_using_pi_gt_3 = 3 * coefficient_pi
    assert A_wedge_upper == 8
    assert coefficient_pi == 80
    assert A_thermal_lower_using_pi_gt_3 == 240 > A_wedge_upper
    # r=AT>2*pi/eta; whole center trajectory needs 2*sinh(r)/A<=width.
    assert 2 / width == 8
    assert 2 / eta == 10
    return {
        "status": "PASS_LIMITED_MENU_SCREEN",
        "new_science_rounds": 0,
        "new_adopted_axioms": 0,
        "ell": str(ell), "coordinate_window": str(width), "log_ratio_tolerance": str(eta),
        "compact_wedge_requires_A_strictly_less_than": str(A_wedge_upper),
        "old_certificate_requires_A_greater_than_pi_times": str(coefficient_pi),
        "strict_thermal_lower_using_pi_greater_than_3": str(A_thermal_lower_using_pi_gt_3),
        "two_existing_numerical_menus_compatible": False,
        "whole_flat_center_worldtube_also_requires": "A > 8*sinh(10*pi)",
        "not_claimed": [
            "This flat comparison worldline is an implemented FLRW trajectory.",
            "All finite thermometry is impossible.",
            "A different finite smear scale, switching proof, or calibrated window is excluded.",
            "The observation is a new uncertainty theorem or a full GR conclusion."
        ]
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = calculate()
    path = Path(__file__).with_name("results.json")
    if args.write:
        with path.open("x", encoding="utf-8") as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
        action = "FIRST_WRITE"
    else:
        assert json.loads(path.read_text("utf-8")) == result
        action = "READ_ONLY_MATCH"
    print(json.dumps({"status": result["status"], "action": action, "science_rounds": 0}))


if __name__ == "__main__":
    main()
