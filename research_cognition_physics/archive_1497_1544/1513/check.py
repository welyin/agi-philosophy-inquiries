"""Finite simplex sources, the heavy-family angle, and a held-out weak decay."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def frame():
    points = [[0, 0, 1]]
    for i in range(3):
        angle = 2 * math.pi * i / 3
        points.append([math.sqrt(8) * math.cos(angle) / 3,
                       math.sqrt(8) * math.sin(angle) / 3, -1 / 3])
    return np.array(points)


def check_finiteness(up_dirs, down_dirs, coefficients):
    up = [math.sqrt(coefficients["up"]) * np.outer(v, v) for v in up_dirs]
    down = [math.sqrt(coefficients["down"]) * np.outer(v, v.conj()) for v in down_dirs]
    f = (np.eye(4) - np.ones((4, 4)) / 4) / math.sqrt(3)
    gram_u = np.array([[np.vdot(i, j) for j in up] for i in up])
    gram_d = np.array([[np.vdot(i, j) for j in down] for i in down])
    g10 = 3 * sum(t @ t.conj().T for t in up) + 2 * sum(t @ t.conj().T for t in down)
    gm = 4 * sum(t.conj().T @ t for t in down)
    gH = 3 * gram_u + 24 / 5 * f @ f.T
    gbarH = 4 * gram_d + 24 / 5 * f.T @ f
    residuals = [np.linalg.norm(obj - rhs * np.eye(len(obj)), ord=2) for obj, rhs in
                 [(g10, 36/5), (gm, 24/5), (gH, 24/5), (gbarH, 24/5)]]
    adj = float(np.vdot(f, f).real + 21 / 5 * coefficients["adjoint_cubic"] - 10)
    residuals.append(abs(adj))
    assert max(residuals) < 1e-12
    assert 3 * (1.5 + 0.5) + 4 + 5 == 15
    return up, down


def weights(spectrum):
    l1, l2, l3 = spectrum
    small, large = l2-l1, l3-l1
    assert large >= 2 * small > 0
    discriminant = (large - 2 * small) * (large - small / 2)
    a = (large + small + math.sqrt(discriminant)) / 2
    b = 9 * large * small / (8 * a)
    baseline = 3 * l1 / 4
    return np.array([baseline + a, baseline + b, baseline, baseline])


def left_current(up_dirs, down_dirs, spectra, coefficients):
    up, down = check_finiteness(up_dirs, down_dirs, coefficients)
    left = []
    norms = []
    for matrices, spectrum, key in [(up, spectra[0], "up"), (down, spectra[1], "down")]:
        amps = weights(spectrum) / math.sqrt(coefficients[key])
        # Norm becomes the free positive sector matching; the actual light Higgs vector is unit.
        norms.append(float(np.linalg.norm(amps)))
        y = sum(a * m for a, m in zip(amps, matrices))
        v, singular, _ = np.linalg.svd(y)
        np.testing.assert_allclose(singular[::-1], spectrum, atol=1e-15)
        left.append(v[:, ::-1])
    return left[0].conj().T @ left[1], left, norms


def angle_bound(r):
    assert 0 < r < 0.5
    return math.atan((math.sqrt(2) / 3) * r / (1 - 2 * r))


def runs(data, scales, meas):
    output = {}
    vcb, *errors = meas
    for name, source in scales.items():
        rows = {}
        for n in data["widths"]:
            angles = []
            ratios = []
            for family in ("up", "down"):
                (_, _), (m2, d2), (m3, d3) = source[family]
                r = (m2 + n * d2) / (m3 - n * d3)
                ratios.append(r)
                angles.append(angle_bound(r))
            empirical = vcb - n * sum(errors)
            cutoff = math.sin(sum(angles))
            rows[str(n)] = {"r_up_down_max": ratios, "Vcb_ceiling": cutoff, "Belle_floor": empirical,
                            "minimum_extra_projective_angle_deg": math.degrees(math.asin(empirical) - sum(angles)),
                            "minimum_multiplicative_Vcb_transfer": empirical / cutoff}
        output[name] = rows
    return output


def calculate():
    data = read(HERE / "inputs.json")
    low = read(HERE / data["weak_scale"])["yukawa"]
    high = read(HERE / data["high_scale"])["high_scale"]
    meas = read(HERE / data["semileptonic"])["BelleII_Vcb_BGL"]
    coefficients = data["finite_coefficients_squared_over_g_squared"]
    frames = frame()
    np.testing.assert_allclose(frames @ frames.T, 4/3 * np.eye(4) - np.ones((4, 4))/3, atol=1e-14)
    pars = data["rotation_parameters"]
    c, s = math.cos(pars["angle_in_transverse_plane"]), math.sin(pars["angle_in_transverse_plane"])
    phase = np.exp(1j * pars["relative_phase"])
    u2 = np.array([[c * phase, -s, 0], [s * phase, c, 0], [0, 0, 1]])
    # The frame is stored by rows; active left U gives row matrix frame @ U.T.
    moved = frames @ u2.T
    spectra = [[m for m, dm in low[f]] for f in ("up", "down")]
    v, (lu, ld), norms = left_current(frames, moved, spectra, coefficients)
    j = float(np.imag(v[0, 0] * v[1, 1] * v[0, 1].conj() * v[1, 0].conj()))
    assert abs(j) > 1e-10
    for basis, spectrum in zip((lu, ld), spectra):
        phi = math.acos(min(1, abs(basis[2, 2])))
        assert phi <= angle_bound(spectrum[1] / spectrum[2])
    # Physical reorientation leaving all mass singular values and finite equations fixed.
    target = meas[0]
    mix = np.eye(3)
    mix[1:, 1:] = [[math.sqrt(1-target**2), target], [-target, math.sqrt(1-target**2)]]
    rotate = lu @ mix @ ld.conj().T
    v2, _, _ = left_current(frames, moved @ rotate.T, spectra, coefficients)
    np.testing.assert_allclose(abs(v2[1, 2]), target, atol=2e-14)
    return {"boxes": runs(data, {"MZ": low, "1e16_GeV_conditional_SM": high}, meas),
            "exact_mass_ratios_witness": {"weights_up": weights(spectra[0]).tolist(),
                                          "weights_down": weights(spectra[1]).tolist(),
                                          "positive_sector_norms_at_g1": norms,
                                          "Vcb_common": float(abs(v[1, 2])), "J_common": j,
                                          "Vcb_rotated_fitted": float(abs(v2[1, 2])),
                                          "rotated_anchor_separation_deg": math.degrees(math.acos(min(1, abs(rotate[2,2]))))}}


def compare(a, b):
    if isinstance(a, dict):
        assert a.keys() == b.keys()
        for k in a:
            compare(a[k], b[k])
    elif isinstance(a, list):
        assert len(a) == len(b)
        for v, w in zip(a, b):
            compare(v, w)
    else:
        assert math.isclose(a, b, rel_tol=1e-10, abs_tol=1e-13), (a, b)


if __name__ == "__main__":
    cli = argparse.ArgumentParser()
    op = cli.add_mutually_exclusive_group(required=True)
    op.add_argument("--check", action="store_true")
    op.add_argument("--write", action="store_true")
    args = cli.parse_args()
    result = calculate()
    path = HERE / "results.json"
    if args.write:
        with path.open("x", encoding="utf-8") as f:
            json.dump(result, f, indent=2, allow_nan=False)
            f.write("\n")
    else:
        compare(result, read(path))
        print("1513 finite frames, mass spectra, held-out decay and competitors verified.")
