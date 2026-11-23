"""Test if two fixed, tagged Yukawas admit one reciprocal family identification."""
import argparse
import itertools
import json
import runpy
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def links(matrix):
    left, singular, right_h = np.linalg.svd(matrix)
    assert singular[-1] > 0 and np.all(np.diff(singular) < 0)
    # A simultaneous overall minus sign is absorbed into the tested scalar phase.
    for tail in itertools.product((1, -1), repeat=2):
        signs = np.array((1, *tail))
        bridge = (right_h.conj().T * signs) @ left.conj().T
        lifted = matrix @ bridge
        np.testing.assert_allclose(lifted, lifted.conj().T, atol=1e-15)
        yield signs.tolist(), bridge


def compare(down, lepton):
    scale = np.linalg.norm(lepton)
    e_links = list(links(lepton))
    gaps, distances = [], []
    for signs, bridge in links(down):
        x = lepton @ bridge
        phase = np.exp(-.5j * np.angle(np.trace(x @ x)))
        hermitian = (phase*x + (phase*x).conj().T) / 2
        best = hermitian @ bridge.conj().T / phase
        relative_error = np.linalg.norm(lepton - best) / scale
        analytic = (1 - abs(np.trace(x @ x)) / scale**2) / 2
        np.testing.assert_allclose(relative_error**2, analytic, atol=1e-14)
        gaps.append((relative_error, signs))
        for e_signs, e_bridge in e_links:
            overlap = np.trace(bridge.conj().T @ e_bridge)
            shifted = e_bridge * np.exp(-1j*np.angle(overlap))
            distance = np.linalg.norm(bridge - shifted)
            np.testing.assert_allclose(distance**2, 6-2*abs(overlap), atol=2e-14)
            distances.append((distance, signs, e_signs))
    return dict(minimum_relative_lepton_vertex_change=min(gaps)[0],
                down_signs_at_minimum_vertex_change=min(gaps)[1],
                minimum_Frobenius_bridge_distance=min(distances)[0],
                closest_signs_down_and_lepton=min(distances)[1:],
                all_relative_gaps_given_down_signs=[dict(down_signs=b, relative_change=a)
                                                   for a, b in gaps])


def compute():
    source = runpy.run_path(str(HERE / "../1715/check.py"))
    _, rep = source["charged_background"](.2)
    down = source["unpack"](rep["down_texture_10PeV"])
    lepton = source["unpack"](rep["charged_texture_10PeV"]).T
    result = compare(down, lepton)
    rng = np.random.default_rng(123)
    def unitary():
        return np.linalg.qr(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))[0]
    left, right = unitary(), unitary()
    moved = compare(left @ down @ right, np.exp(.67j) * left @ lepton @ right)
    for key in ("minimum_relative_lepton_vertex_change", "minimum_Frobenius_bridge_distance"):
        np.testing.assert_allclose(result[key], moved[key], rtol=1e-10, atol=1e-14)
    # Exact positive control: the same bridge with independent Hermitian matrices.
    reciprocal = np.diag([.03, .4, .8]) @ next(links(down))[1].conj().T
    assert compare(down, reciprocal)["minimum_relative_lepton_vertex_change"] < 2e-14
    return result


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()
    result = json.loads(json.dumps(compute()), parse_float=lambda x: float(f"{float(x):.10g}"))
    file = HERE / "independent_results.json"
    if args.check:
        assert json.loads(file.read_bytes()) == result
        print("1723 verified all simultaneous signed polar factors and the exact nearest vertex.")
    else:
        with file.open("x", encoding="utf8") as stream:
            stream.write(json.dumps(result, indent=2, allow_nan=False) + "\n")
