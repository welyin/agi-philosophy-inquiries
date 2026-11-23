"""Replay the reciprocal three-family determinant test. No fit or phase scan."""
import argparse
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def load(name):
    return json.loads((HERE / name).read_text(encoding="utf-8-sig"))


def coefficients(values, t):
    a, b, c, d = values
    A, B, C, D = a*a-d*d, a*b-d*c, b*b-c*c, a*c-d*b
    assert D*D == B*B-A*C
    sine = [(t+t**5)*D+t**3*C, (t**2+t**4)*B, t**3*A]
    difference = [(1-t**6)*A+(t**2-t**4)*C, 2*(t-t**5)*B, 2*(t**2-t**4)*D, 0]
    return sine, difference


def exact_verification():
    count = 0
    for values in itertools.product(range(-2, 3), repeat=4):
        for t in (Fraction(1, 4), Fraction(1), Fraction(3, 2)):
            p = [v*t**i for i, v in enumerate(values)]
            q = [v*t**i for i, v in enumerate(reversed(values))]
            sine, dif = coefficients(values, t)
            for k in range(4):
                corr = sum(p[j+k]*p[j]-q[j+k]*q[j] for j in range(4-k))
                assert dif[k] == (1 if k == 0 else 2)*corr
                if k:
                    assert sine[k-1] == sum(p[j]*q[j+k]-q[j]*p[j+k] for j in range(4-k))
            assert (not any(sine)) == (p == q or p == [-x for x in q])
            count += 1
    return count


def mass_ratio_spectrum(data, widths):
    def ratio(k):
        return math.prod(a+k*s for a, s in data["up"])/math.prod(a-k*s for a, s in data["down"])
    low = ratio(-widths)
    return {"central": ratio(0), "two_widths": [low, ratio(widths)],
            "independent_symmetric_threshold_minimum": (low**(1/3)-1)/(low**(1/3)+1)}


def finite_bound(t, ratio_low):
    a = t**-3
    b = 1/(t**2+t**4)
    w = t**-2+t*t
    d = (w*a+math.hypot(w*a, 2*math.hypot(b, a)))/2
    c = a+w*d
    k = abs(1-t**6)*a+abs(t*t-t**4)*c+2*abs(t-t**5)*b+2*abs(t*t-t**4)*d
    return {"K": k, "max_imag_over_up_volume_squared_minimum": 2*(1-ratio_low**-2)/((1+math.sqrt(2))*k)}


def fixed_vacuum(data, charged, t):
    s12, s13, s23 = [entry[0] for entry in data["quark_sines_and_phase"][:3]]
    c12, c13, c23 = [math.sqrt(1-s*s) for s in (s12, s13, s23)]
    gamma = math.radians(charged["central_degrees"])
    delta = gamma+math.asin(c12*s23*s13*math.sin(gamma)/(s12*c23))
    z = np.exp(1j*delta)
    V = np.array([[c12*c13,s12*c13,s13/z],
                  [-s12*c23-c12*s23*s13*z,c12*c23-s12*s23*s13*z,s23*c13],
                  [s12*s23-c12*c23*s13*z,-c12*s23-s12*c23*s13*z,c23*c13]])
    U = np.diag([v[0] for v in data["yukawa"]["up"]])
    D = V @ np.diag([v[0] for v in data["yukawa"]["down"]]) @ V.conj().T
    A, B = (U-t*D)/(1-t*t), (D-t*U)/(1-t*t)
    assert np.allclose(A, A.conj().T) and np.allclose(B, B.conj().T)
    assert np.max(np.abs(A+t*B-U)) < 1e-14 and np.max(np.abs(B+t*A-D)) < 1e-14
    resolved = np.angle(-V[0,0]*V[0,2].conj()/(V[1,0]*V[1,2].conj()), deg=True)
    determinant = np.linalg.det(A+t*B)*np.linalg.det(B+t*A)
    assert determinant.real > 0 and abs(np.angle(determinant)) < 1e-12
    assert abs(resolved-math.degrees(gamma)) < 1e-10
    return {"gamma_reconstructed_deg": float(resolved), "delta_deg": math.degrees(delta),
            "J": float(np.imag(V[0,0]*V[1,1]*np.conj(V[0,1]*V[1,0]))),
            "strong_phase_zero_verified": True}


def isolated_complex_vacuum():
    discriminant = 924**2+4*337*685
    assert math.isqrt(discriminant)**2 < discriminant
    k = (924+math.sqrt(discriminant))/(2*337)
    A = np.diag([1, 2, 4])

    def second_kernel(x):
        return np.array([[x, 1, 2], [1, 3, 1], [2, 1, 2]])

    B = second_kernel(k)
    U, D = A-0.5j*B, B+0.5j*A
    H, G = U@U.conj().T, D@D.conj().T
    product = np.linalg.det(U)*np.linalg.det(D)
    assert product.real > 470 and abs(product.imag) < 1e-11
    # Exact reduction of 2U and 2D in F5 with k -> 0 and i -> 2.
    prime, x, imaginary_unit = 5, 0, 2
    assert (337*x*x-924*x-685) % prime == (imaginary_unit**2+1) % prime == 0
    B = second_kernel(x)
    H_mod = (2*A-imaginary_unit*B)@(2*A+imaginary_unit*B) % prime
    G_mod = (2*B+imaginary_unit*A)@(2*B-imaginary_unit*A) % prime
    a, b, c = (H_mod@G_mod-G_mod@H_mod) % prime
    residue = int(a@np.cross(b, c)) % prime
    assert residue == 4
    return {"k": k, "positive_determinant_product": float(product.real),
            "up_volume": float(abs(np.linalg.det(U))), "down_volume": float(abs(np.linalg.det(D))),
            "weak_commutator_imag_determinant": float(np.linalg.det(H@G-G@H).imag),
            "nonzero_scaled_commutator_mod5": residue}


def rounded(value):
    if isinstance(value, dict):
        return {key: rounded(item) for key, item in value.items()}
    if isinstance(value, list):
        return [rounded(item) for item in value]
    if isinstance(value, float):
        return float(format(value, ".11g"))
    return value


def certificate():
    config = load("inputs.json")
    low = load(config["low_mass_and_angles"])
    high = load(config["high_mass"])["high_scale"]
    charged = load(config["charged_B"])["chargedB"]
    results = {tag: mass_ratio_spectrum(table, config["mass_box_marginal_widths"])
               for tag, table in (("MZ", low["yukawa"]), ("1e16_GeV_SM", high))}
    response = {str(t): finite_bound(t, results["MZ"]["two_widths"][0])
                for t in config["background_ratio_examples"]}
    holdout = [charged["central_degrees"]-charged["lower_error"]-config["extra_transfer_deg"],
               charged["central_degrees"]+charged["upper_error"]+config["extra_transfer_deg"]]
    return rounded({"exact_coefficient_checks": exact_verification(), "mass_volume_ratio": results,
                    "three_phase_bounds_MZ": response, "charged_B_gamma_work_box_mod180_deg": holdout,
                    "fixed_real_vacuum": fixed_vacuum(low, charged, config["fixed_vacuum_ratio"]),
                    "isolated_complex_vacuum": isolated_complex_vacuum()})


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--write-new", action="store_true")
    group.add_argument("--check", action="store_true")
    args = parser.parse_args()
    output = HERE / "results.json"
    result = certificate()
    if args.write_new:
        with output.open("x", encoding="utf-8") as stream:
            json.dump(result, stream, indent=2, ensure_ascii=False)
            stream.write("\n")
    elif json.loads(output.read_text(encoding="utf-8")) != result:
        raise SystemExit("Saved certificate differs")
    print("1579 verified")
