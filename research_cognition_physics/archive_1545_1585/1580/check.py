"""Finite word and physical-window certificates for round 1580."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE / name).read_text(encoding="utf-8-sig"))


def cyclic(word):
    return min(word[i:]+word[:i] for i in range(len(word)))


def orientation(word):
    key, reverse = cyclic(word), cyclic(word[::-1])
    return (key < reverse)-(key > reverse)


def audit_words():
    weights = {}
    for n in range(5):
        values = {}
        for letters in itertools.product("EF", repeat=n):
            word = "".join(letters)
            routes = [word, word[::-1]]
            total = sum(orientation("EF"+v)-orientation("FE"+v) for v in routes)
            assert total % 2 == 0
            if total:
                values[word] = total//2
        weights[str(n)] = values
    assert all(not weights[str(n)] for n in range(4))
    expected = {"EEFF": -1, "FFEE": -1, "EFEF": -1, "FEFE": -1, "EFFE": 2, "FEEF": 2}
    assert weights["4"] == expected
    return weights


def reduce_box(function, ends):
    values = [function(*x) for x in itertools.product(*ends)]
    return [min(values), max(values)]


def sin_extrema(lo, hi):
    points = [lo, hi]
    points.extend(.5+n for n in range(math.ceil(lo-.5), math.floor(hi-.5)+1))
    values = [math.sin(math.pi*x) for x in points]
    return [min(values), max(values)]


def pmns(q12, q13, q23, delta):
    s, r, v = np.sqrt([q12, q13, q23])
    c, w, u = np.sqrt([1-q12, 1-q13, 1-q23])
    z = np.exp(1j*delta)
    return np.array([[c*w, s*w, r/z],
                     [-s*u-c*v*r*z, c*u-s*v*r*z, v*w],
                     [s*v-c*u*r*z, -c*v-s*u*r*z, u*w]])


def point_shapes(x, r, data, q23, delta, ordering):
    U = pmns(data["s12"][0], data["s13"][0], q23, math.radians(delta))
    eigen = [0, r, 1] if ordering == "NO" else [1-r, 1, 0]
    F = U@np.diag(eigen)@U.conj().T
    E = np.diag([0, x, 1])
    C = E@F-F@E
    D = C.conj().T@C
    J = (U[0,0]*U[1,1]*np.conj(U[0,1]*U[1,0])).imag
    actual = float(np.trace(C@D).imag)
    analytic = -6*x*(1-x)*r*(1-r)*J
    assert abs(actual-analytic) < 1e-13
    eig = np.linalg.eigvalsh(D).clip(min=0)
    assert abs(math.sqrt(eig[2])-math.sqrt(eig[0])-math.sqrt(eig[1])) < 2e-12
    H = -1j*C
    w, basis = np.linalg.eigh(H)
    distance = abs(actual)/sum(abs(w))
    shift = -np.sign(actual)*distance*(basis@np.diag(np.sign(w))@basis.conj().T)
    assert abs(np.trace(C@(D+shift)).imag) < 1e-14
    assert abs(np.trace(C@(E@F+F@E)).imag) < 1e-13
    return {"J": float(J), "source": actual, "dirac_shape_spectrum": [float(e) for e in eig],
            "exact_unrestricted_hermitian_error_to_zero": float(distance)}


def ordering_result(ordering, inputs, data, phases, x_window, x0):
    q12, q13, q23 = data["s12"], data["s13"], inputs["theta23_sin2"][ordering]
    solar, atm = data["solar"], data["atm_abs"]
    r = [solar[1]/atm[2], solar[2]/atm[1]]
    nu_spread = [z*(1-z) for z in r]
    ell_spread = [z*(1-z) for z in x_window]
    solar_factor = [math.sqrt(z*(1-z)) for z in q12[1:]]
    reactor_factor = [(1-z)*math.sqrt(z) for z in q13[1:]]
    atmos_factor = [min(math.sqrt(z*(1-z)) for z in q23[1:]), .5]
    assert q23[1] <= .5 <= q23[2]
    amplitude = [6*math.prod(v[i] for v in [nu_spread, ell_spread, solar_factor,
                                              reactor_factor, atmos_factor]) for i in range(2)]
    delta = phases["phase_credible_arcs_pi"][ordering+"_99.73"]
    pad = phases["joint_published_rounding_guard_pi"]
    sine = sin_extrema(delta[0]-pad, delta[1]+pad)
    source = reduce_box(lambda a, b: -a*b, [amplitude, sine])
    absolute_min = max(source[0], -source[1], 0.0)
    result = {"solar_fraction": r, "coefficient_of_minus_sin_delta": amplitude,
              "delta_sine_window": sine, "signed_source_window": source,
              "all_box_guaranteed_absolute_source": absolute_min,
              "example": point_shapes(x0, solar[0]/atm[0], data, q23[0],
                                      inputs["order_example_delta_deg"][ordering], ordering)}
    if absolute_min:
        result["absolute_theta_transfer_cap_if_no_other_contributions"] = (
            inputs["theta_max_under_neutron_dictionary"]/absolute_min)
    return result


def type_ii_example(inputs):
    d = inputs["typeII_example"]
    # An independently chosen gain and heavy-mass case at local tree matching.
    M0 = d["min_mass_eV"]
    MR_min = d["vR_eff_GeV"]*M0/d["vL_abs_eV"]
    md = .25*d["kappa_GeV"]*d["dirac_gain"]
    type_i_eV = md*md/MR_min*1e9
    max_mass = math.sqrt(M0*M0+.0026)
    delta_H = type_i_eV*(2*max_mass+type_i_eV)
    # Delta of the target reference II spectrum >=.0024 eV^2; normalizing
    # the two spectra separately costs up to four times the squared-mass change.
    F_error = 4*delta_H/.0024
    return {"dirac_mass_norm_max_GeV": md, "heavy_min_GeV": MR_min,
            "typeI_zero_momentum_norm_upper_eV": type_i_eV,
            "mixing_expansion_parameter": md/MR_min, "max_F_change_leading": F_error,
            "max_normalized_quartic_source_change_leading": 2.25*F_error}


def rounded(item):
    if isinstance(item, dict):
        return {key: rounded(value) for key, value in item.items()}
    if isinstance(item, list):
        return [rounded(value) for value in item]
    if isinstance(item, float):
        return float(format(item, ".11g"))
    return item


def certificate():
    inputs = read("inputs.json")
    charged = read(inputs["charged"])["yukawa"]["charged_leptons"]
    ranges = [[a-inputs["charged_widths"]*s, a+inputs["charged_widths"]*s] for a,s in charged]
    relative = lambda e,m,t: (m-e)/(t-e)
    x = reduce_box(relative, ranges)
    x0 = relative(*(a for a,_ in charged))
    assert 0 < x[0] <= x[1] < .5
    oscillation, phase = read(inputs["angles_gaps"]), read(inputs["phase"])
    return rounded({"word_orientation_weights": audit_words(), "mu_relative_charged_position": x,
                    "ordering": {o: ordering_result(o, inputs, oscillation[o], phase, x, x0)
                                 for o in ("NO", "IO")}, "typeII_example": type_ii_example(inputs)})


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--check", action="store_true")
    group.add_argument("--write-new", action="store_true")
    args = parser.parse_args()
    target = HERE / "results.json"
    result = certificate()
    if args.write_new:
        with target.open("x", encoding="utf-8") as stream:
            json.dump(result, stream, indent=2)
            stream.write("\n")
    elif json.loads(target.read_text(encoding="utf-8")) != result:
        raise SystemExit("Saved certificate differs")
    print("1580 verified")
