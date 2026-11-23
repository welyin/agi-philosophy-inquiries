"""Round 1598. Run with --check; --write-new refuses to overwrite the saved result."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE / path).read_text(encoding="utf-8-sig"))


def box(record, widths):
    center, *errors = record
    spread = widths * sum(errors)
    return center - spread, center + spread


def mixing(us, cb, ratio, gamma):
    s13 = ratio * cb
    c13 = math.sqrt(1 - s13**2)
    s12, s23 = us/c13, cb/c13
    c12, c23 = math.sqrt(1-s12**2), math.sqrt(1-s23**2)
    delta = gamma + math.asin(c12*s23*s13*math.sin(gamma)/(s12*c23))
    phase = np.exp(1j*delta)
    return np.array([[c12*c13, s12*c13, s13/phase],
                     [-s12*c23-c12*s23*s13*phase, c12*c23-s12*s23*s13*phase, s23*c13],
                     [s12*s23-c12*c23*s13*phase, -c12*s23-s12*c23*s13*phase, c23*c13]])


def distance_box(matrix, cells):
    edges = np.array([matrix[0, 1], matrix[1, 2], matrix[2, 0]])
    a, b, c = sorted(abs(edges))
    angle = abs(math.remainder(np.angle(np.prod(edges)), math.pi))
    grid = np.linspace(0, angle, cells+1)
    s = b*b + c*c
    product = 4*b*b*c*c*np.sin(angle-grid)**2
    second = product / (2*(s + np.sqrt(np.maximum(s*s-product, 0))))
    first = (a*np.sin(grid))**2
    low = np.min(first[:-1]+second[1:])
    high = np.min(first+second)
    guard = 1e-12  # relative padding for evaluation, much smaller than the cell enclosure
    return [float(math.sqrt(max(0, low))*(1-guard)), float(math.sqrt(high)*(1+guard))]


def central(y, us, cb, ratio, gamma, cells):
    V = mixing(us, cb, ratio, gamma)
    h = (V*y**2) @ V.conj().T
    gap = y[2]**2-y[0]**2
    T = (h-np.eye(3)*y[0]**2)/gap
    q = (y[1]**2-y[0]**2)/gap
    J = (V[0, 0]*V[1, 1]*np.conj(V[0, 1]*V[1, 0])).imag
    edge_product = T[0, 1]*T[1, 2]*T[2, 0]
    assert np.allclose(V@V.conj().T, np.eye(3), atol=3e-15)
    assert math.isclose(abs(edge_product.imag), abs(J)*q*(1-q), rel_tol=2e-13)
    reconstructed = np.angle(-V[0, 0]*V[0, 2].conj()/(V[1, 0]*V[1, 2].conj()))
    assert abs(math.remainder(reconstructed-gamma, 2*math.pi)) < 1e-12
    lower, upper = distance_box(T, cells)
    tags = {}
    for j, label in enumerate(["u", "c", "t"]):
        r = V[j].conj()/y
        intensity = 1/float(np.vdot(r, r).real)
        r *= math.sqrt(intensity)
        column = V@(y*r)
        assert np.allclose(column, math.sqrt(intensity)*np.eye(3)[j], atol=2e-16)
        remaining = h-np.outer(column, column.conj())
        eig, vec = np.linalg.eigh(remaining)
        assert abs(eig[0])/gap < 1e-12 and eig[1] > 0
        factory = np.column_stack([column, vec[:, 1:]*np.sqrt(eig[1:])])
        assert np.allclose(factory@factory.conj().T, h, rtol=1e-10, atol=1e-19)
        tags[label] = {"left_column_strength":intensity, "down_flavor_fractions":(abs(r)**2).tolist(),
                       "second_eigenvalue_over_gap":float(eig[1]/gap)}
    return {"abs_J":abs(float(J)), "q":q, "phase_mod_pi_deg":abs(math.degrees(math.remainder(np.angle(edge_product), math.pi))),
            "signed_distance_over_gap_interval":[lower,upper], "positive_extra_eigenvalue_necessary_over_gap":2*lower,
            "positive_extra_eigenvalue_necessary_Y_squared":2*lower*gap, "tags":tags}


def certificate():
    inputs = read("inputs.json")
    y_data = read(inputs["yukawa_file"])["yukawa"]["down"]
    trees = read(inputs["trees_file"])
    cb_data = read(inputs["belle_file"])["BelleII_Vcb_BGL"]
    gamma_data = read(inputs["gamma_file"])["chargedB"]
    y = np.array(y_data)[:, 0]
    central_values = (trees["K_product"][0]/trees["K_form"][0], cb_data[0], trees["ub_over_cb"][0])
    branches = {str(offset):central(y, *central_values, math.radians(gamma_data["central_degrees"]+offset), inputs["grid"]) for offset in [0,180]}
    n = inputs["n_width"]
    dl, dh = box(y_data[0], n)
    sl, sh = box(y_data[1], n)
    bl, bh = box(y_data[2], n)
    qlo, qhi = (sl**2-dh**2)/(bh**2-dl**2), (sh**2-dl**2)/(bl**2-dh**2)
    kl, kh = box(trees["K_product"], n)
    fl, fh = box(trees["K_form"], n)
    usl, ush = kl/fh, kh/fl
    cbl, cbh = box(cb_data,n)
    rl, rh = box(trees["ub_over_cb"], n)
    zlo, zhi = cbl*rl, cbh*rh
    c13 = math.sqrt(1-zhi*zhi)
    s12hi, s23hi = ush/c13, cbh/c13
    c12lo, c23lo = math.sqrt(1-s12hi**2), math.sqrt(1-s23hi**2)
    delta_pad = math.asin(s23hi*zhi/(usl*c23lo))
    gl = math.radians(gamma_data["central_degrees"]-n*gamma_data["lower_error"]-gamma_data["additional_degrees"])-delta_pad
    gh = math.radians(gamma_data["central_degrees"]+n*gamma_data["upper_error"]+gamma_data["additional_degrees"])+delta_pad
    assert 0 < gl < gh < math.pi/2
    jlo = usl*c12lo*cbl*c23lo*zlo*c13**2*math.sin(gl)
    tshi = s23hi+s12hi*zhi
    aa = zhi*cbh + qhi*ush
    bb = cbh + qhi*tshi
    cc = zhi + qhi*tshi*ush
    distance = jlo*qlo*(1-qhi)/math.sqrt((aa*bb)**2+(bb*cc)**2+(cc*aa)**2)
    # Absolute values bounding the d, s and b columns for each row. Triangle
    # inequalities avoid any phase/ordering correlations; they also hold at CP=0.
    dlow = [c12lo*c13, usl*c23lo-s23hi*zhi, usl*cbl-zhi]
    shi = [ush, 1, tshi]
    bhi = [zhi, cbh, 1]
    fractions_lo = [1/(1+(smax*dh/(dmin*sl))**2+(bmax*dh/(dmin*bl))**2)
                    for dmin,smax,bmax in zip(dlow,shi,bhi)]
    assert min(dlow)>0
    ycoh = np.diag([1.,2.,3.]) + np.outer(np.ones(3), np.conj([1,1j,1+1j]))
    hcoh = ycoh@ycoh.conj().T
    witness_cp = (hcoh[0,1]*hcoh[1,2]*hcoh[2,0]).imag
    assert abs(witness_cp)>1
    return {"central_gamma_offsets_degrees":branches, "two_width_box":{
        "Vus":[usl,ush], "Vcb":[cbl,cbh], "ratio_ub_cb":[rl,rh], "q":[qlo,qhi],
        "absolute_J_lower":jlo, "distance_over_gap_lower":distance,
        "additional_positive_eigenvalue_over_gap_lower":2*distance,
        "additional_positive_eigenvalue_Y_squared_lower":2*distance*(bl**2-dh**2),
        "down_tag_fraction_lower_u_c_t":fractions_lo},
        "coherent_rank_one_Y_addition_imaginary_triangle":float(witness_cp)}


def equivalent(left, right):
    if isinstance(left, dict):
        return left.keys() == right.keys() and all(equivalent(left[k], right[k]) for k in left)
    if isinstance(left, list):
        return len(left) == len(right) and all(equivalent(a,b) for a,b in zip(left,right))
    return math.isclose(left, right, rel_tol=2e-10, abs_tol=3e-21)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--check", action="store_true")
    group.add_argument("--write-new", action="store_true")
    args = parser.parse_args()
    result = certificate()
    if args.write_new:
        with (HERE/"results.json").open("x", encoding="utf-8") as stream:
            json.dump(result, stream, indent=2, ensure_ascii=False, allow_nan=False)
            stream.write("\n")
    else:
        if not equivalent(result, read("results.json")):
            raise SystemExit("1598 saved certificate differs")
        print("1598 verified: two tree branches, spectral bounds, six right-end constructions, and coherent competitor.")
