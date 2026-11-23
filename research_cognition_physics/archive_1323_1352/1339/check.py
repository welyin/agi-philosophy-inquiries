"""Rebuild the finite joint formation isometry and two light charged-current sectors."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


def tensor(factors):
    product = np.array([1.0 + 0j])
    for factor in factors:
        product = np.kron(product, factor)
    return product


def invariant(left, right):
    return float(np.linalg.det(left @ right - right @ left).imag)


def quartet(mixing):
    return float((mixing[0, 0] * mixing[1, 1] * mixing[0, 1].conj()
                  * mixing[1, 0].conj()).imag)


def evaluate(gate, initial, down_weights, reverse_access=False):
    identity = np.eye(3)
    omega = identity.reshape(-1) / math.sqrt(3)
    rays = np.stack([np.linalg.matrix_power(gate, j) @ initial for j in range(3)], axis=1)
    gram = rays.conj().T @ rays
    loop = gram[0, 1] * gram[1, 2] * gram[2, 0]
    norm = np.linalg.norm
    # Contract a physical conjugate-basic program port against incoming data.
    choi = gate / math.sqrt(3)  # output, conjugate input
    assert norm(choi.conj().T @ choi - identity / 3) < 1e-13
    assert norm(choi @ choi.conj().T - identity / 3) < 1e-13
    joint = np.einsum("oi,xy->oixy", choi, identity)
    bell = omega.reshape(3, 3)
    good = np.einsum("ix,oixy->oy", bell, joint)
    bad = joint - np.einsum("ix,oy->oixy", bell, good)
    assert norm(good - gate / 3) < 1e-13
    assert norm(bad.reshape(-1, 3).conj().T @ bad.reshape(-1, 3) - 8 * identity / 9) < 1e-12
    wedge = np.zeros((9, 3), dtype=complex)
    for j in range(3):
        k, l = (j + 1) % 3, (j + 2) % 3
        wedge[3 * k + l, j] = 1 / math.sqrt(2)
        wedge[3 * l + k, j] = -1 / math.sqrt(2)
    assert norm(wedge.conj().T @ np.kron(gate, gate) @ wedge
                - np.linalg.det(gate) * gate.conj()) < 1e-12
    resources = tensor([v for ray in rays.T for v in (ray, ray.conj())])
    free = np.kron(resources.reshape(-1, 1), identity).reshape((3,) * 7 + (3,))
    branches = [free] + [np.swapaxes(free, 2 * j, 6) for j in range(3)]
    formed = np.vstack([v.reshape(-1, 3) for v in branches]) / 2
    assert norm(formed.conj().T @ formed - identity) < 1e-12
    bell_bank = tensor([omega] * 3)
    spectra, light, leading = [], [], []
    matrices = []
    success = []
    up_weights = [1, 1, 2, 0]
    if reverse_access:
        # Physically relabel the first and last controlled swap in both sectors.
        up_weights = [up_weights[i] for i in (0, 3, 2, 1)]
        down_weights = [down_weights[i] for i in (0, 3, 2, 1)]
    for weights, t in ((up_weights, 1 / 3), (down_weights, 1 / 4)):
        weights = np.asarray(weights)
        chi = weights / norm(weights)
        port = np.kron(tensor([chi, bell_bank]).reshape(-1, 1), identity)
        reduced = port.conj().T @ formed
        target = weights[0] * identity + (rays * weights[1:]) @ rays.conj().T
        factor = 1 / (2 * math.sqrt(27) * norm(weights))
        assert norm(reduced - factor * target) < 2e-13
        rejected = formed - port @ reduced
        assert norm(rejected.conj().T @ rejected + reduced.conj().T @ reduced - identity) < 1e-12
        val, vec = np.linalg.eigh(target)
        d = factor * val
        total = 1 + t * t
        z = 2 * (t * d)**2 / (total + np.sqrt(total**2 - 4 * (t * d)**2))
        k = t*t*d / (total - z)
        eig = (formed @ vec - (port @ vec) * k) / np.sqrt(1 + k*k - 2*k*d)
        h_eig = t*t * port @ (port.conj().T @ eig) + eig - formed @ (formed.conj().T @ eig)
        assert norm(h_eig - eig * z) < 2e-13 and norm(eig.conj().T @ eig - identity) < 2e-13
        light.append(eig)
        leading.append(vec)
        spectra.append(np.sqrt(z).tolist())
        success.append((d*d).tolist())
        matrices.append(target)
    physical = light[0].conj().T @ light[1]
    return dict(a=float(abs(gram[0, 1])**2), b=float(abs(gram[0, 2])**2),
                t=float(loop.real), h=float(loop.imag), det_gram=float(np.linalg.det(gram).real),
                imag_det_comm=invariant(*matrices), sector_success_eigenvalues=success,
                leading_masses=[np.linalg.eigvalsh(x).tolist() for x in matrices],
                finite_masses_over_M=spectra, finite_abs_left_current=np.abs(physical).tolist(),
                leading_quartet=quartet(leading[0].conj().T @ leading[1]),
                finite_quartet=quartet(physical))


def run():
    gate = np.diag([1, 1j, -1])
    start = np.sqrt([1/2, 1/3, 1/6])
    base = evaluate(gate, start, [1, 0, 1, 3])
    competitor = evaluate(gate, start, [1, 0, 3, 1])
    reverse = evaluate(gate, start, [1, 0, 1, 3], reverse_access=True)
    conjugate = evaluate(gate.conj(), start, [1, 0, 1, 3])
    assert abs(base["imag_det_comm"] + 32/27) < 1e-12
    assert abs(3*competitor["imag_det_comm"] - base["imag_det_comm"]) < 1e-12
    assert np.allclose(base["finite_masses_over_M"], competitor["finite_masses_over_M"])
    assert np.allclose(base["finite_masses_over_M"], reverse["finite_masses_over_M"])
    assert abs(base["finite_quartet"] + reverse["finite_quartet"]) < 1e-12
    assert abs(base["finite_quartet"] + conjugate["finite_quartet"]) < 1e-12
    assert abs(base["finite_quartet"] - competitor["finite_quartet"]) > .02
    rng = np.random.default_rng(1339)
    for _ in range(10):
        q, _ = np.linalg.qr(rng.normal(size=(3, 3)) + 1j*rng.normal(size=(3, 3)))
        v = rng.normal(size=3) + 1j*rng.normal(size=3)
        entry = evaluate(q, v / np.linalg.norm(v), [1, 0, 1, 3])
        assert abs(entry["imag_det_comm"] + 36*entry["det_gram"]*entry["h"]) < 1e-12
    delta, risk = .001, .01
    # Four balanced blocks of a single registered neutral calibration protocol.
    samples = math.ceil(2*math.log(8/risk) / delta**2)
    return dict(round=1339, reference_dimension=3, complete_left_dimension=4*3**7,
                direct_checks_random_cases=10, neutral_precision=delta, failure_bound=risk,
                accepted_triples_per_setting=samples, sample=base, competitor=competitor,
                reverse_access=reverse, conjugate=conjugate)


def numerically_equal(x, y):
    if isinstance(x, dict):
        return x.keys() == y.keys() and all(numerically_equal(x[k], y[k]) for k in x)
    return np.allclose(x, y, atol=1e-12, rtol=1e-12)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    output = Path(__file__).with_name("results.json")
    result = run()
    if args.write:
        with output.open("x", encoding="utf8") as handle:
            json.dump(result, handle, ensure_ascii=False, indent=2)
            handle.write("\n")
    else:
        assert numerically_equal(result, json.loads(output.read_text(encoding="utf8")))
    print(json.dumps({k: result["sample"][k] for k in ("imag_det_comm", "finite_quartet")}))
