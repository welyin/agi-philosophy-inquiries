"""Count0 algebra checks: a zero direct mass entry is not a zero LNV kernel."""
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import numpy as np


def polynomial_add(a, b):
    result = [F(0)] * max(len(a), len(b))
    for i, value in enumerate(a):
        result[i] += value
    for i, value in enumerate(b):
        result[i] += value
    return result


def calculate():
    # All masses/momenta below are in one arbitrary common unit mu.
    masses = (F(1), F(3))
    residues = (F(-3, 4), F(3, 4))
    assert sum(residues) == 0
    numerator = polynomial_add(
        [residues[0]*masses[1]**2, residues[0]],
        [residues[1]*masses[0]**2, residues[1]])
    assert numerator == [F(-6), F(0)]
    assert sum(c*m*m for c, m in zip(residues, masses)) == 6

    rows = []
    for x in (F(1), F(2), F(4)):
        rational = sum(c/(x+m*m) for c, m in zip(residues, masses))
        assert rational == -6/((x+1)*(x+9))
        rows.append({"dimensionless_spacelike_s": str(x),
                     "mu_times_kernel": str(rational)})
    # The interval bound follows from the increasing positive denominator.
    lower, upper = F(6, 65), F(3, 10)
    assert lower == 6/((F(4)+1)*(F(4)+9))
    assert upper == 6/((F(1)+1)*(F(1)+9))

    root3 = np.sqrt(3.0)
    mass_matrix = np.array([[0.0, root3], [root3, 2.0]])
    U = np.array([[1j*root3/2, 0.5], [-0.5j, root3/2]])
    D = np.diag([1.0, 3.0])
    unitarity_error = float(np.linalg.norm(U.conj().T@U-np.eye(2), ord=2))
    takagi_error = float(np.linalg.norm(U.T@mass_matrix@U-D, ord=2))
    reconstruction_error = float(np.linalg.norm(U.conj()@D@U.conj().T-mass_matrix, ord=2))
    computed_residues = U[0, :]**2 * np.diag(D)
    residue_error = float(np.max(np.abs(computed_residues-np.array([-0.75, 0.75]))))
    kernel_error = 0.0
    for row in rows:
        x = float(F(row["dimensionless_spacelike_s"]))
        spectral = np.sum(computed_residues/(x+np.diag(D)**2))
        # Here M is real symmetric; its same-chirality resolvent is M/(s+M^2).
        direct = np.linalg.solve(x*np.eye(2)+mass_matrix@mass_matrix, mass_matrix)[0, 0]
        exact = float(F(row["mu_times_kernel"]))
        kernel_error = max(kernel_error, float(abs(spectral-exact)), float(abs(direct-exact)))
    assert max(unitarity_error, takagi_error, reconstruction_error,
               residue_error, kernel_error) < 2e-15

    # Degenerate pair: grouped residue vanishes, so the full rational function is zero.
    dirac_residues = (F(-1, 2), F(1, 2))
    assert sum(dirac_residues) == 0
    dirac_polynomial = polynomial_add([dirac_residues[0]], [dirac_residues[1]])
    assert dirac_polynomial == [F(0)]
    UD = np.array([[1j, 1.0], [-1j, 1.0]])/np.sqrt(2.0)
    MD = np.array([[0.0, 1.0], [1.0, 0.0]])
    dirac_takagi_error = float(np.linalg.norm(UD.T@MD@UD-np.eye(2), ord=2))
    assert dirac_takagi_error < 1e-15
    # Negative control: beta weights |U_ei|^2 are not coherent U_ei^2.
    wrong_coherent_sum = F(3, 4)*masses[0]+F(1, 4)*masses[1]
    assert wrong_coherent_sum == F(3, 2)
    return {
        "scope": "count0 exact algebra; no nuclear rates, empirical fit, or off-shell measurement",
        "new_scientific_groups": 0,
        "new_empirical_groups": 0,
        "new_cognitive_principles": 0,
        "common_mass_unit_mu_is_arbitrary_positive": True,
        "toy_masses_over_mu": [str(m) for m in masses],
        "toy_residues_over_mu": [str(c) for c in residues],
        "rational_numerator_ascending": [str(c) for c in numerator],
        "leading_high_s_residue_sum": str(sum(residues)),
        "next_mass_weighted_moment": "6",
        "finite_identity_checks": rows,
        "interval": {"x_min": "1", "x_max": "4", "abs_mu_K_min": str(lower),
                     "abs_mu_K_max": str(upper), "proof": "monotonic positive denominator"},
        "matrix_diagnostics": {"unitarity_error": unitarity_error,
                               "takagi_error": takagi_error,
                               "reconstruction_error": reconstruction_error,
                               "residue_error": residue_error,
                               "direct_vs_spectral_kernel_error": kernel_error,
                               "dirac_takagi_error": dirac_takagi_error},
        "dirac_grouped_residue": str(sum(dirac_residues)),
        "wrong_use_of_beta_weights_leading_sum": str(wrong_coherent_sum),
        "all_assertions_passed": True}


if __name__ == "__main__":
    result = calculate()
    path = Path(__file__).with_name("results.json")
    if sys.argv[1:] == ["--save-exclusive"]:
        with path.open("x", encoding="utf-8") as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
    elif sys.argv[1:]:
        raise SystemExit("Use no arguments for read-only checks, or --save-exclusive once.")
    else:
        assert result == json.loads(path.read_text(encoding="utf-8"))
    print("PASS: zero-block, finite LNV kernel, and Dirac grouped-residue identities; count0.")
