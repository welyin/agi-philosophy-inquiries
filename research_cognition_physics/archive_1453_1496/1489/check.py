"""Certify the three-edge loop phase with closed intervals; construct a matched countermodel."""
import argparse
import importlib.util
import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("interval1486", ROOT.parent / "1486/check.py")
interval = importlib.util.module_from_spec(spec)
spec.loader.exec_module(interval)
I = interval.I
SCALE = interval.SCALE


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def box(row, widths):
    central, error = map(I.point, row)
    return interval.hull([central - widths * error, central + widths * error])


def cosine(sine):
    return (1 - sine.square()).sqrt()


def bound_phase(data, inputs, widths):
    yd, ys, yb = [box(row, widths) for row in data["yukawa"]["down"]]
    p, u, t = [box(row, widths) for row in data["quark_sines_and_phase"][:3]]
    a, b, r = (yd/yb).square(), (ys/yb).square(), (yd/ys).square()
    cp, ct = cosine(p), cosine(t)
    ratio = u * (1 - a * cp.square() - b * p.square()) / (p * cp * (b-a))
    alpha = t * ct * (ratio - 1/ratio)
    upper_left = p.square() + r * cp.square()
    lower_right = u.square() * (cp.square() + r * p.square()) + a * (1-u.square())
    beta = t * ct * (upper_left-lower_right) / (u*p*cp*(1-r))
    d = 1-2*t.square()
    bound = I.point(0)
    for z in I(-SCALE,SCALE).parts(inputs["cosine_cells"]):
        sin_delta = (1-z.square()).sqrt()
        # Two factors have positive real parts and the same sign of imaginary part.
        tan1 = sin_delta / (alpha + d*z)
        tan2 = sin_delta / (beta - d*z)
        tangent = (tan1+tan2)/(1-tan1*tan2)
        bound = interval.hull([bound, tangent])
    # The strict comparison is done on the outward rational tangent interval.
    # Degrees are only a readable representation.
    assert bound.square().hi < 3*SCALE
    return {"alpha":alpha.values(), "beta":beta.values(), "upper_tan":bound.values()[1],
            "upper_loop_angle_deg":math.degrees(math.atan(bound.values()[1])),
            "min_real_extra_over_net":(1-bound.values()[1]/math.sqrt(3))/math.hypot(1,bound.values()[1])}


def matrix(sines, delta):
    p,u,t = sines
    cp,cu,ct = np.sqrt(1-np.square(sines))
    z = np.exp(1j*delta)
    return np.array([[cp*cu,p*cu,u/z],
                     [-p*ct-cp*t*u*z,cp*ct-p*t*u*z,t*cu],
                     [p*t-cp*ct*u*z,-cp*t-p*ct*u*z,ct*cu]])


def construction(data):
    d = np.array([row[0] for row in data["yukawa"]["down"]])
    d /= d[-1]
    p,u,t,delta = [row[0] for row in data["quark_sines_and_phase"]]
    V = matrix([p,u,t], delta)
    H = (V*d**2)@V.conj().T
    basis = np.diag([1,np.exp(1j*np.angle(H[0,1])),np.exp(1j*np.angle(H[0,2]))])
    K = basis@H@basis.conj().T
    a = np.sqrt(K[0,0].real)
    c,e = (K[0,1:].real/a)
    f = np.sqrt(K[2,2].real-e**2)
    z = (K[1,2]-c*e)/f
    b = -np.sqrt(K[1,1].real-c**2-abs(z)**2)
    Y = np.array([[0,a,0],[b,c,z],[0,e,f]])
    assert np.linalg.norm(Y@Y.conj().T-K)<1e-13
    det = np.linalg.det(Y)
    assert det.real>0 and abs(det.imag)<1e-15
    invariant = -H[0,1]*H[2,0]*np.linalg.inv(H)[1,2]
    assert abs(np.angle(invariant/z))<1e-10
    j = (V[0,0]*V[1,1]*V[0,1].conj()*V[1,0].conj()).imag
    eta = -2*math.pi/3
    d1 = z.imag/math.sin(eta)
    d0 = z.real-d1*math.cos(eta)
    assert abs(z-d0-d1*np.exp(1j*eta))<1e-14
    gamma = np.angle(-V[0,0]*V[0,2].conj()/(V[1,0]*V[1,2].conj()),deg=True)
    return {"normalized_real_edges_A_B_C_E_F":[a,b,c,e,f], "Z_real_imag":[z.real,z.imag],
            "measured_loop_phase_deg":np.angle(z,deg=True), "tree_det_in_yb_units":det.real,
            "CP_quartet":j, "gamma_deg_at_matched_point":gamma,
            "eta_example_deg":-120,"two_real_coefficients":[d0,d1],"real_extra_over_net":abs(d0/z),
            "C_and_E_removed_would_make_CP_zero":True}


def solve():
    inputs=read(ROOT/"inputs.json")
    data=read(ROOT/inputs["matched"])
    certificates={str(n):bound_phase(data, inputs, n) for n in inputs["widths"]}
    return {"cos_delta_continuous_certificates":certificates,"reconstruction":construction(data)}


if __name__=="__main__":
    p=argparse.ArgumentParser()
    g=p.add_mutually_exclusive_group(required=True)
    g.add_argument("--check",action="store_true")
    g.add_argument("--write-new",action="store_true")
    args=p.parse_args()
    result=solve()
    path=ROOT/"results.json"
    if args.write_new:
        with path.open("x",encoding="utf-8") as f:
            json.dump(result,f,ensure_ascii=False,indent=2)
            f.write("\n")
    else:
        assert result==read(path)
        print("1489: full CKM phase interval, determinant and coherent two-path witness verified.")
