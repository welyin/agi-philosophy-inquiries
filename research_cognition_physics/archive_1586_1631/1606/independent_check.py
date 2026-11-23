"""Independent bordered-Hermitian completion test; input Yukawas at 10^16 GeV."""
import argparse
import itertools
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def frozen(number):
    return json.loads((HERE.parent / str(number) / "inputs.json").read_text(encoding="utf-8-sig"))


def quartet_product(h):
    return (h[0, 1] * h[1, 2] * h[2, 0]).imag


def rotation(i, j, angle, phase=0):
    matrix = np.eye(3, dtype=complex)
    matrix[i, i] = matrix[j, j] = np.cos(angle)
    matrix[i, j] = np.sin(angle) * np.exp(-1j * phase)
    matrix[j, i] = -matrix[i, j].conjugate()
    return matrix


def bordered_reading(h, charged, m):
    others = [a for a in range(3) if a != m]
    column = h[others, m]
    spread = np.linalg.norm(column)
    mean = h[m, m].real
    weights = []
    for i, eigen in enumerate(charged):
        other = np.delete(charged, i)
        weights.append((spread**2 + np.prod(mean - other)) / np.prod(eigen - other))
    result = {"row_mass_probabilities": weights}
    if min(weights) < 0:
        return result
    vandermonde = np.prod([charged[j] - charged[i] for i,j in itertools.combinations(range(3), 2)])
    size = vandermonde * np.sqrt(np.prod(weights)) / spread**2
    cap = abs(np.prod(column)) * size
    demand = abs(9 * quartet_product(h))
    last_diag = mean + np.prod(mean - charged) / spread**2
    mid_diag = sum(charged) - mean - last_diag
    # Rephase the distinguished column to real, then rotate into its direction.
    n = abs(column) / spread
    axes = np.array([[1,0,0], [0,n[0],-n[1]], [0,n[1],n[0]]])
    fraction = demand / cap
    for angle in [np.pi/2, np.arcsin(fraction)]:
        link = size * np.exp(1j * angle)
        tri = np.array([[mean,spread,0], [spread,mid_diag,link], [0,link.conjugate(),last_diag]])
        witness = axes @ tri @ axes.T
        assert np.allclose(np.linalg.eigvalsh(witness), charged, atol=2e-12)
        assert np.allclose(witness[1:,0], abs(column))
        assert np.isclose(abs(quartet_product(witness)), cap * np.sin(angle), rtol=1e-9)
    result.update(required_over_max=fraction, max_absolute_clebsch=np.sqrt(cap/abs(quartet_product(h))))
    return result


def calculate():
    high = frozen(1397)["high_scale"]
    a,b,c,phase = np.asarray(high["angles_theta12_theta13_theta23_delta"])[:,0]
    mixing = rotation(1,2,c) @ rotation(0,2,b,phase) @ rotation(0,1,a)
    down = np.asarray(high["down"])[:,0]
    charged = np.asarray(frozen(1402)["scales"]["1e16"])[3:,0]
    h = mixing @ np.diag((down/down[-1])**2) @ mixing.conj().T
    rows = {label: bordered_reading(h, (charged/down[-1])**2, m) for m,label in enumerate("uct")}
    assert 0.45 < rows["u"]["required_over_max"] < 0.46
    assert 0.070 < rows["c"]["required_over_max"] < 0.071
    assert min(rows["t"]["row_mass_probabilities"]) < 0
    # Bare real rays and two independent phases. Include all choices of the omitted row.
    rng = np.random.default_rng(1606)
    errors = []
    for omitted in range(3):
        for trial in range(32):
            base = rng.normal(size=(3,3)) + 3*np.eye(3)
            indices = [i for i in range(3) if i != omitted]
            rays = np.linalg.solve(base, np.eye(3)[:,indices])
            rays /= np.linalg.norm(rays, axis=0)
            z = rng.normal(size=2) + 1j*rng.normal(size=2)
            change = (rays*z) @ rays.T
            grams = []
            for factor in [1,-3]:
                mat = base @ (np.eye(3) + factor*change)
                grams.append(mat@mat.conj().T)
            a1, a2 = grams
            assert np.allclose(a1[:,omitted], a2[:,omitted], atol=1e-12)
            errors.append(abs(quartet_product(a2)-9*quartet_product(a1)) /
                          max(1.,abs(quartet_product(a2))))
    assert max(errors) < 2e-10
    return {"scale_GeV": 1e16, "k": 1, "columns": rows, "random_bare_examples": len(errors),
            "max_scaled_invariant_residual_rounded": float('%.2g' % max(errors))}


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--check", action="store_true")
    args = p.parse_args()
    out = HERE / "independent_results.json"
    result = calculate()
    if args.check:
        assert result == json.loads(out.read_text(encoding="utf-8"))
    else:
        with out.open("x", encoding="utf-8") as fp:
            fp.write(json.dumps(result, indent=2, allow_nan=False)+"\n")
    print("1606 independent check passed")
