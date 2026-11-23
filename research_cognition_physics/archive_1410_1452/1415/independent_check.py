"""Independent CKM quartet bound and explicit two-doublet vacuum, using NumPy."""
import argparse
import itertools
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def rotation(i, j, sine, phase=0.0):
    matrix = np.eye(3, dtype=complex)
    matrix[i, i] = matrix[j, j] = np.sqrt(1 - sine**2)
    matrix[i, j] = sine * np.exp(-1j * phase)
    matrix[j, i] = -matrix[i, j].conjugate()
    return matrix


def quartet(matrix, rows, cols):
    i, j = rows
    a, b = cols
    return np.array([matrix[i, a], matrix[j, b],
                     matrix[i, b].conjugate(), matrix[j, a].conjugate()])


def weight(entries):
    return np.sqrt(sum(np.prod(np.delete(np.abs(entries), j))**2
                       for j in range(4)))


def limits(data, width, matrix, v, higgs):
    low, high = data[:, 0] - width * data[:, 1], data[:, 0] + width * data[:, 1]
    s12, s13, s23, delta = high
    c12, c13, c23 = np.sqrt(1 - high[:3]**2)
    assert 0 < low[3] <= delta < np.pi / 2
    j_lower = np.prod(low[:3]) * c12 * c13**2 * c23 * np.sin(low[3])
    if width:
        # The (u,c ; d,b) choice. Bound each modulus over the entire box.
        w = weight([1.0, s23, s13, s12 + s23*s13])
    else:
        pairs = list(itertools.combinations(range(3), 2))
        w = min(weight(quartet(matrix, rows, cols)) for rows in pairs for cols in pairs)
    phase_weighted_chi = j_lower**2 / (2 * w**2)
    gap_upper = 4 * (16 * np.pi) * v**2 / 3 - 3*higgs**2
    return {
        "width": width, "quartet_modulus_weight_upper": w, "abs_J_lower": j_lower,
        "each_chi_times_sin2phi_lower": phase_weighted_chi,
        "author_convention_all_three_each_S_lower": 4 * phase_weighted_chi,
        "each_125_ordered_pairs_floor_at_600_GeV_ceiling":
            phase_weighted_chi * (higgs**2/(600**2-higgs**2))**2,
        "each_125_ordered_pairs_floor_with_original_unitarity":
            phase_weighted_chi * (higgs**2/gap_upper)**2,
    }


def scalar_example(phi, v, higgs):
    theta = phi / 2
    common, mixed, fifth = 0.8, -0.55, 0.9
    kernel = np.array([
        [common+mixed+2*fifth*np.cos(theta)**2, 0, -fifth*np.sin(phi)],
        [0, common-mixed, 0],
        [-fifth*np.sin(phi), 0, 2*fifth*np.sin(theta)**2],
    ])
    eigen, basis = np.linalg.eigh(kernel)
    scale = (higgs/v)**2 / eigen[0]
    masses = v*np.sqrt(scale*eigen)
    l1 = l2 = scale * common
    l5, l4, l3 = scale*fifth, 0.0, scale*(mixed+fifth)
    positive_eigenchannels = [
        l1+l2+np.hypot(l1-l2, 2*l5),
        l1+l2+np.hypot(l1-l2, 2*l4),
        3*(l1+l2)+np.hypot(3*(l1-l2), 2*(2*l3+l4)),
        *(2*abs(value) for value in [
            l3+l4, l3-l4, l3+l5, l3-l5, l3+2*l4+3*l5, l3+2*l4-3*l5]),
    ]
    assert eigen[0] > 0 and l5 > l4 and l1 > abs(scale*mixed)
    assert max(positive_eigenchannels) < 16*np.pi
    assert np.isclose(masses[0], higgs)
    assert abs(kernel[0, 2]) <= (eigen[-1]-eigen[0])*np.sqrt(1-basis[0, 0]**2)
    return {
        "theta": theta, "tan_beta": 1, "quartics_in_paper_order": [l1,l2,l3,l4,l5],
        "m_neutral_GeV": masses.tolist(), "m_charged_GeV": v*np.sqrt(l5-l4),
        "R1h_squared": basis[0,0]**2, "largest_abs_unitarity_channel": max(positive_eigenchannels),
    }


def top_blind(matrix, v, higgs):
    # Local reordered coordinates (t,u,c) and (b,d,s), restored below.
    order = [2,0,1]
    inverse = np.argsort(order)
    r = matrix[np.ix_(order, order)]
    s13 = abs(r[0, 2])
    c13 = np.sqrt(1-s13**2)
    s12, s23 = abs(r[0, 1])/c13, abs(r[1, 2])/c13
    c12, c23 = np.sqrt(1-s12**2), np.sqrt(1-s23**2)
    cos_phase = (abs(r[1,0])**2-(s12*c23)**2-(c12*s23*s13)**2) / (2*s12*c23*c12*s23*s13)
    j = np.prod(quartet(r, (0,1), (0,1))).imag
    sin_phase = j/(s12*c12*s23*c23*s13*c13**2)
    assert np.isclose(np.hypot(cos_phase, sin_phase), 1)
    phi = np.arctan2(sin_phase, cos_phase)
    left = rotation(1,2,s23).real[inverse,:]
    right = (rotation(0,2,s13) @ rotation(0,1,s12)).real[:,inverse]
    reconstructed = left @ np.diag([1,1,np.exp(1j*phi)]) @ right
    mod_error = np.max(abs(abs(reconstructed)-abs(matrix)))
    j_error = abs(np.prod(quartet(reconstructed, (0,1),(0,1))).imag - j)
    assert mod_error < 1e-13 and j_error < 1e-15
    spectra = scalar_example(phi, v, higgs)
    sector = {}
    for name, direction, oriented in [
        ("u",left[:,2],reconstructed), ("d",right[2,:],reconstructed.T),
    ]:
        probabilities = direction**2
        row_min = 0.5*sum(1-abs(np.sum(oriented**2, axis=1)))
        direct = 0.5*sum(1-np.sqrt(1-4*probabilities*(1-probabilities)*np.sin(phi)**2))
        chi = 1-sum(probabilities**2)
        assert abs(row_min-direct)<1e-14 and row_min <= 2*chi*np.sin(phi)**2+1e-14
        sector[name] = {
            "n":direction.tolist(), "chi":chi, "row_or_column_imaginary_norm_squared_min":direct,
            "125_ordered_F": 4*(1-spectra["R1h_squared"])*chi,
        }
    assert sector["u"]["n"][2] == 0
    return {
        "sector":sector, "abs_CKM_max_error": mod_error,
        "signed_J":j, "phase":phi, "potential":spectra,
    }


def run():
    inp = json.loads((HERE.parent/"1386"/"inputs.json").read_text(encoding="utf-8"))
    gf = json.loads((HERE.parent/"1404"/"inputs.json").read_text(encoding="utf-8"))["const_GeV"]["GF"]
    pars = np.array(inp["quark_sines_and_phase"])
    s12,s13,s23,delta = pars[:,0]
    matrix = rotation(1,2,s23) @ rotation(0,2,s13,delta) @ rotation(0,1,s12)
    v, higgs = (2**0.5*gf)**(-0.5), 125.0
    return {
        "lower_bounds":[limits(pars, width, matrix, v, higgs) for width in [0, 2]],
        "explicit_CKM_and_vacuum_example":top_blind(matrix,v,higgs),
        "v_GeV_from_1404_GF":v,
    }


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    action=parser.add_mutually_exclusive_group()
    action.add_argument("--write",action="store_true")
    action.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result = json.dumps(run(), ensure_ascii=False, indent=2, allow_nan=False) + "\n"
    dest = HERE/"independent_results.json"
    if args.write:
        dest.write_text(result,encoding="utf-8")
    elif json.loads(dest.read_text(encoding="utf-8")) != json.loads(result):
        raise SystemExit("independent results mismatch")
    print("1415 independent checks passed.")
