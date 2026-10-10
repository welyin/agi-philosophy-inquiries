"""Count-zero Michel POVM checks; no experimental fit or finite-time QFT claim."""
from fractions import Fraction as F
from pathlib import Path
import argparse
import json
import numpy as np

HERE = Path(__file__).resolve().parent


def integral(coeff):
    return sum((F(c, i + 1) for i, c in enumerate(coeff)), F(0))


def run():
    norm = integral([0, 0, 6, -4])
    anisotropy = integral([0, 0, -2, 4])
    hemisphere = anisotropy / 4
    total_variation = anisotropy / 2
    assert (norm, anisotropy, hemisphere, total_variation) == (
        F(1), F(1, 3), F(1, 12), F(1, 6)
    )
    sx = np.array([[0, 1], [1, 0]], complex)
    sy = np.array([[0, -1j], [1j, 0]], complex)
    sz = np.diag([1, -1]).astype(complex)
    pauli = np.array([sx, sy, sz])
    ident = np.eye(2)
    rng = np.random.default_rng(1055003)
    max_difference = 0.0
    min_eigenvalue = 1.0
    max_trace_error = 0.0
    for _ in range(24):
        # Unknown mixed input on mu (2) times R (3), not a product preparation.
        z = rng.normal(size=(6, 6)) + 1j * rng.normal(size=(6, 6))
        rho = z @ z.conj().T
        rho /= np.trace(rho)
        u = rng.normal(size=3)
        u /= np.linalg.norm(u)
        e = ident / 2 + np.einsum("i,ijk->jk", u, pauli) / 12
        effects = (e, ident - e)
        marginal_sum = np.zeros((3, 3), complex)
        for effect in effects:
            val, vec = np.linalg.eigh(effect)
            root = (vec * np.sqrt(val)) @ vec.conj().T
            k = np.kron(root, np.eye(3))
            full = k @ rho @ k.conj().T
            marginal = np.trace(full.reshape(2, 3, 2, 3), axis1=0, axis2=2)
            left = np.kron(effect, np.eye(3)) @ rho
            direct = np.trace(left.reshape(2, 3, 2, 3), axis1=0, axis2=2)
            max_difference = max(max_difference, float(np.max(np.abs(direct - marginal))))
            min_eigenvalue = min(min_eigenvalue, float(np.linalg.eigvalsh(marginal)[0]))
            marginal_sum += marginal
        original_r = np.trace(rho.reshape(2, 3, 2, 3), axis1=0, axis2=2)
        max_trace_error = max(max_trace_error, float(np.max(np.abs(marginal_sum - original_r))))
    assert max_difference < 2e-14
    assert min_eigenvalue >= -2e-14
    assert max_trace_error < 2e-14
    return {
        "kind": "mature_michel_povm_algebra_only",
        "exact": {"normalization": str(norm), "angular_analyzing_coefficient": str(anisotropy),
                  "hemisphere_spin_coefficient": str(hemisphere), "antipodal_angular_tv": str(total_variation)},
        "reference_diagnostic": {"seed": 1055003, "samples": 24, "input_dimension": 6,
                                 "partial_trace_identity_passed": True,
                                 "reference_marginal_preserved": True, "all_branches_positive": True,
                                 "tolerance": 2e-14},
        "new_science_groups": 0, "new_cognitive_axioms": 0,
        "experimental_data_fit": False, "full_decay_instrument_certified": False,
        "all_checks_passed": True,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    target = HERE / "results.json"
    if args.write:
        with target.open("x", encoding="utf-8") as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
    else:
        assert result == json.loads(target.read_text(encoding="utf-8"))
    print(json.dumps(result, ensure_ascii=False))


if __name__ == "__main__":
    main()
