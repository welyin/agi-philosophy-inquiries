"""Cross-flavour simultaneous-plane necessary conditions; analytic box enclosures."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
CFG = json.loads((ROOT/"inputs.json").read_text(encoding="utf-8-sig"))


def box(entry, n):
    return np.array([entry[0]-n*entry[1], entry[0]+n*entry[-1]])


def normalized_middle(values, n):
    small, middle, large = (box(entry, n)**2 for entry in values)
    return np.array([(middle[0]-small[1])/(large[1]-small[0]),
                     (middle[1]-small[0])/(large[0]-small[1])])


def parameters(label, n):
    if label == "quark":
        x = normalized_middle(CFG["yukawa"]["up"], n)
        y = normalized_middle(CFG["yukawa"]["down"], n)
        angles = [box(val, n) for val in CFG["quark_sines_and_phase"][:3]]
        return x, y, 1, angles
    spec = CFG[label]
    solar = box(spec["solar_gap"], n)
    atm = box(spec["atmospheric_gap"], n)
    if label == "NO":
        y = np.array([solar[0]/(atm[1]+solar[0]), solar[1]/(atm[0]+solar[1])])
        z = 1
    else:
        y = solar / atm[::-1]
        z = -1
    x = normalized_middle(CFG["yukawa"]["charged_leptons"], n)
    amp = box(CFG["DayaBay_amplitude"], n)
    p13 = amp/(2*(1+np.sqrt(1-amp)))
    angles = [np.sqrt(box(spec["solar_sine_squared"], n)), np.sqrt(p13),
              np.sqrt(box(spec["atmospheric_sine_squared"], n))]
    return x, y, z, angles


def mixing(s12, s13, s23, delta):
    c12, c13, c23 = np.sqrt(1-np.square([s12,s13,s23]))
    e = np.exp(1j*delta)
    return np.array([[c12*c13, s12*c13, s13/e],
                     [-s12*c23-c12*s23*s13*e, c12*c23-s12*s23*s13*e, s23*c13],
                     [s12*s23-c12*c23*s13*e, -c12*s23-s12*c23*s13*e, c23*c13]])


def matrix_invariants(label, delta):
    xv, yv, z, angle_boxes = parameters(label, 0)
    x, y = xv[0], yv[0]
    angles = [s[0] for s in angle_boxes]
    u = mixing(*angles, delta)
    a = np.diag([0, x, 1])
    b = (u*np.array([0,y,z])) @ u.conj().T
    a = a - np.trace(a)*np.eye(3)/3
    b = b - np.trace(b)*np.eye(3)/3
    c = 1j*(a@b - b@a)
    norm2 = np.vdot(c,c).real
    area2 = np.vdot(a,a).real*np.vdot(b,b).real - np.vdot(a,b).real**2
    imag = abs((u[0,0]*u[1,1]*u[0,1].conj()*u[1,0].conj()).imag)
    det = 2*x*(1-x)*y*abs(z-y)*imag
    assert abs(det-abs(np.linalg.det(c))) < 2e-14
    return {"spectral_shape": float(3*math.sqrt(6)*det/norm2**1.5),
            "plane_ratio": float(norm2/area2), "abs_J":float(imag)}


def interval_certificate(label, n):
    x, y, z, s = parameters(label, n)
    c12, c13, c23 = [np.sqrt(1-a[::-1]**2) for a in s]
    s12, s13, s23 = s
    # Off-diagonal B is a signed column-3 projector plus y times column 2.
    off23_low = s23[0]*c23[0]*c13[0]**2-y[1]/2
    assert off23_low > 0
    off12_hi = s13[1]*s23[1] + y[1]/2
    off13_hi = s13[1] + y[1]/2
    off23_hi = s23[1]*c23[1] + y[1]/2
    nlo = 2*((1-x[1])*off23_low)**2
    nhi = 2*((x[1]*off12_hi)**2 + off13_hi**2 + off23_hi**2)
    diag11_hi, diag22_hi = s13[1]**2+y[1], s23[1]**2+y[1]
    diag33_hi = 1+y[1]
    diag_norm_low = math.sqrt((1-x[1])**2 + 1 + x[0]**2)
    diag_norm_hi = math.sqrt((1-x[0])**2 + 1 + x[1]**2)
    perp_diag_hi = (diag11_hi + diag22_hi + x[1]*diag33_hi)/diag_norm_low
    perp_diag_lo = max(0, s23[0]**2*c13[0]**2-y[1]-diag11_hi-x[1]*diag33_hi)/diag_norm_hi
    anorm_lo = 2*(1-x[1]+x[0]**2)/3
    anorm_hi = 2*(1-x[0]+x[1]**2)/3
    bnorm_hi = 2*(1+y[1]**2+(-z)*y[0 if z==1 else 1])/3
    k_lo = (1-x[1])**2/anorm_hi/(1+(perp_diag_hi**2+2*off12_hi**2)/(2*off23_low**2))
    k_hi = (1-perp_diag_lo**2/bnorm_hi)/anorm_lo
    j0 = s12*c12*s13*c13**2*s23*c23
    # Each gap product is positive, the boxes cover both signs of the atmospheric gap.
    last = np.sort(abs(z-y))
    coef = 6*math.sqrt(6)*x*(1-x[::-1])*y*last*j0
    strength = coef / np.array([nhi, nlo])**1.5
    if label == "quark":
        strength *= np.sin(box(CFG["quark_sines_and_phase"][-1], n))
    return {"plane_lo_hi":[float(k_lo),float(k_hi)],
            "shape_or_shape_over_abs_sine":strength.tolist(),
            "norm2_lower":float(nlo)}


def phase_root(label, target, origin):
    # On the allowed CP bands there is one positive offset near each real matrix.
    lo, hi = 0., 0.01
    assert matrix_invariants(label, origin+hi)["spectral_shape"] > target
    for _ in range(48):
        mid = (lo+hi)/2
        if matrix_invariants(label, origin+mid)["spectral_shape"] > target:
            hi = mid
        else:
            lo = mid
    return (lo+hi)/2


def evaluate():
    dq = CFG["quark_sines_and_phase"][-1][0]
    quark = matrix_invariants("quark", dq)
    out = {"central_quarks":quark, "orderings":{}, "boxes":{}}
    for label in ("NO","IO"):
        offsets = [phase_root(label,quark["spectral_shape"],origin) for origin in (0.,math.pi)]
        out["orderings"][label] = {"shape_only_offset_from_zero_rad":offsets[0],
                                    "shape_only_offset_from_pi_rad":offsets[1],
                                    "leptons_at_pi_over_two":matrix_invariants(label,math.pi/2),
                                    "plane_ratio_at_shape_only_root":matrix_invariants(label,offsets[0])["plane_ratio"]}
    for n in CFG["error_multipliers"]:
        cert = {label:interval_certificate(label,n) for label in ("quark","NO","IO")}
        for label in ("NO","IO"):
            q, ell = cert["quark"], cert[label]
            lo, hi = np.array(q["shape_or_shape_over_abs_sine"])/np.array(ell["shape_or_shape_over_abs_sine"])[::-1]
            ell["shape_only_abs_sine_necessary"] = [float(lo),float(hi)]
            ell["forbidden_plane_gap"] = q["plane_lo_hi"][0] - ell["plane_lo_hi"][1]
            assert ell["forbidden_plane_gap"] > 0 and 0 < lo < hi < .001
        out["boxes"][str(n)] = cert
    return out


def same(a,b):
    if isinstance(a, dict):
        return a.keys() == b.keys() and all(same(a[k],b[k]) for k in a)
    if isinstance(a,list):
        return len(a)==len(b) and all(same(x,y) for x,y in zip(a,b))
    return math.isclose(a,b,rel_tol=1e-9,abs_tol=1e-14)


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write-new", action="store_true")
    mode.add_argument("--check", action="store_true")
    args=parser.parse_args()
    result=evaluate()
    path=ROOT/"results.json"
    if args.write_new:
        with path.open("x",encoding="utf8") as out:
            json.dump(result,out,indent=2)
            out.write("\n")
    else:
        assert same(result,json.loads(path.read_text(encoding="utf-8-sig")))
    print("1386 OK")
