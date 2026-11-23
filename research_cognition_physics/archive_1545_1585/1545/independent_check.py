"""Checks of the Z width for arbitrary mass/mixing and a conditional QCD threshold budget."""
import argparse
import json
from math import exp, log, pi, sqrt
from pathlib import Path

import numpy as np


MZ = 91.1876
LEP = 2.9963 + 2 * .0074
LATTICE = .11823
LATTICE_ERROR = .00084
R32 = [(474., .0936, .0041), (664., .0894, .0031), (896., .0889, .0034)]


def axial_weight(m):
    return max(1 - (2 * m / MZ) ** 2, 0) ** 1.5


def ordered_kernel(a, b, c):
    """Returns the ordered-sum convention for a massive Majorana pair."""
    if a + b >= 1:
        return 0.
    x, y = a*a, b*b
    momentum = sqrt((1-(a+b)**2) * (1-(a-b)**2))
    return momentum * ((1-(x+y)/2-(x-y)**2/2) * abs(c)**2
                       - 3*a*b*(c*c).real)


def invisible(projector, masses):
    return sum(ordered_kernel(a/MZ, b/MZ, projector[i, j])
               for i, a in enumerate(masses) for j, b in enumerate(masses))


def minimum_heavy_weight(n, r, allowed):
    if allowed < n-r:
        raise ValueError("incompatible heavy rank")
    return r * (1-sqrt(max(0., 1-max(n-allowed, 0.)/r)))


def projection_witness(n, r, heavy_weight):
    # First n masses are accessible; last r exceed MZ so all their pairs close.
    v = np.zeros((n, n+r))
    for i in range(n):
        if i < r:
            fraction = heavy_weight/r
            v[i, i], v[i, n+i] = sqrt(1-fraction), sqrt(fraction)
        else:
            v[i, i] = 1
    return v.T @ v


def calculate_z():
    n, widths = 6, {}
    rng = np.random.default_rng(1545)
    for r in (3, 4, 5, 6):
        weight = minimum_heavy_weight(n, r, LEP)
        gram = projection_witness(n, r, weight)
        assert np.allclose(gram @ gram, gram)
        exact = invisible(gram, [0.]*n+[100.]*r)
        assert abs(exact - LEP) < 2e-14
        widths[str(r)] = {"minimum_hidden_active_weight": weight,
                          "saturating_width": round(exact, 12)}
    for _ in range(30):
        a = rng.normal(size=(10, n))+1j*rng.normal(size=(10, n))
        u = np.linalg.qr(a)[0]
        gram = u @ u.conj().T
        masses = [*rng.uniform(0, 20, 7), 35., 60., 130.]
        compressed = gram[:7, :7]
        width_floor = axial_weight(20.)*float(np.vdot(compressed, compressed).real)
        assert invisible(gram, masses) >= width_floor - 1.e-11
    # A Dirac state is two degenerate Majorana states of opposite relative phase.
    gram_d = np.array([[1, 1j], [-1j, 1]])/2
    for m in range(0, 45):
        x = (m/MZ)**2
        assert np.isclose(invisible(gram_d, [m, m]), (1-x)*sqrt(1-4*x))
        assert np.isclose(invisible(np.eye(1), [m]), axial_weight(m))
    scales = {str(k): MZ/2*sqrt(1-(LEP/k)**(2/3)) for k in (4, 5, 6, 9)}
    return {"two_width_upper": LEP, "masses_GeV_for_k_unhidden_full_directions": scales,
            "N6_light_massless_other_directions": widths,
            "arbitrary_complex_projectors_checked": 30}


def calculate_qcd():
    mz, mt = MZ, 172.5
    upper_anchor = 1/(LATTICE-2*LATTICE_ERROR)
    data = {}
    for q, alpha, error in R32:
        no_new = upper_anchor + ((23/3)*log(mt/mz)+7*log(q/mt))/(2*pi)
        end = 1/(alpha+2*error)
        bounds = {}
        for response_error in (0., .2):
            ln_budget = 3*pi*(no_new - end + response_error)
            bounds[str(response_error)] = {
                "sum_positive_logs_upper": ln_budget,
                "six_clipped_masses_geometric_mean_lower_GeV":
                    max(mz, q*exp(-ln_budget/6))}
        data[str(int(q))] = bounds
    # Once all are crossed, the measured endpoint depends on their product only.
    split = [250., 250., 400., 400., 640., 640.]
    q = 896.
    assert abs(sum(log(q/m) for m in split)-6*log(q/400.)) < 1.e-12
    limit = data["664"]["0.0"]["sum_positive_logs_upper"]
    missing = (6*log(664./250.) - limit)/(3*pi)
    return {"two_width_working_box": data, "equal_product_above_threshold_GeV": split,
            "minimum_664_inverse_coupling_transfer_six_250_GeV": round(missing, 11)}


def main():
    parser = argparse.ArgumentParser()
    opts = parser.add_mutually_exclusive_group(required=True)
    opts.add_argument("--check", action="store_true")
    opts.add_argument("--write-new", action="store_true")
    args = parser.parse_args()
    data = {"invisible_Z": calculate_z(), "conditional_CMS_transfer": calculate_qcd()}
    target = Path(__file__).with_name("independent_results.json")
    if args.write_new:
        with target.open("x", encoding="utf-8") as out:
            json.dump(data, out, ensure_ascii=False, indent=2)
            out.write("\n")
    else:
        assert data == json.loads(target.read_text(encoding="utf-8"))
        print("PASS 1545 independent neutral-current width and split coloured thresholds")


if __name__ == "__main__":
    main()
