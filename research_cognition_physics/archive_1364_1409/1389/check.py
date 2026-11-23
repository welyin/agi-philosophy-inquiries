"""Round 1389: symmetry-composed reverse operators at a fixed Nelson-Barr background."""
import argparse
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def det_integer(matrix):
    """Exact elimination for a small integer matrix (used for the real source algebra)."""
    a = [[Fraction(int(v)) for v in row] for row in matrix]
    det = Fraction(1)
    for j in range(len(a)):
        p = next((k for k in range(j,len(a)) if a[k][j]), None)
        if p is None:
            return 0
        if p != j:
            a[j], a[p] = a[p], a[j]
            det *= -1
        det *= a[j][j]
        for k in range(j+1,len(a)):
            factor = a[k][j]/a[j][j]
            for i in range(j+1,len(a)):
                a[k][i] -= factor*a[j][i]
    assert det.denominator == 1
    return det.numerator


def mass_matrix(m, b, heavy, reverse):
    return np.block([[m, reverse[:,None]], [b[None,:], np.array([[heavy]])]])


def observables(matrix):
    h = matrix @ matrix.conj().T
    _, u = np.linalg.eigh(h)
    v = u[:3, :3]
    quartet = (v[0,0]*v[1,1]*v[0,1].conj()*v[1,0].conj()).imag
    return {"theta_rad":float(np.angle(np.linalg.det(matrix))), "light_quartet":float(quartet)}


def right_rotation(m, b, heavy, alpha):
    t, s = math.cos(alpha), math.sin(alpha)
    cos_phi = (1+t*t-abs(b[0]/heavy)**2*s*s)/(2*t)
    phi = math.acos(cos_phi)
    endpoint = np.exp(1j*phi)
    p = heavy*(1-t*endpoint)/(b[0]*s)
    assert abs(abs(p)-1) < 1e-11
    rot = np.eye(4,dtype=complex)
    rot[0,0],rot[0,3] = t,s*p
    rot[3,0],rot[3,3] = -s*endpoint*p.conjugate(),t*endpoint
    original = mass_matrix(m,b,heavy,np.zeros(3))
    changed = original @ rot
    assert np.linalg.norm(rot.conj().T@rot-np.eye(4)) < 1e-12
    assert abs(changed[-1,-1]-heavy) < 1e-13
    assert np.max(abs(changed[:3,:3].imag)) == 0.
    assert np.linalg.det(changed[:3,:3]).real > 0
    residual = np.linalg.norm(changed@changed.conj().T-original@original.conj().T)
    assert residual < 1e-12
    # With <S>=1+i any complex vector Z is g S + h S* for two real vectors.
    for vector in (changed[3,:3],changed[:3,3]):
        g=(vector.real+vector.imag)/2
        h=(vector.real-vector.imag)/2
        assert np.max(abs(g*(1+1j)+h*(1-1j)-vector)) < 1e-13
    return {"phase_shift_rad":phi, "actual":observables(changed),
            "left_gram_residual":float(residual), "constant_mass_real":float(changed[-1,-1].real)}


def evaluate():
    inputs = read(ROOT/"inputs.json")
    m = np.array(inputs["ordinary_mass"],dtype=float)
    b = np.array(inputs["mixing_row_real"])+1j*np.array(inputs["mixing_row_imag"])
    h = np.diag(inputs["H_u_diag"])
    heavy = inputs["heavy_mass"]
    scale = inputs["cutoff_in_mass_units"]
    old = read(ROOT/inputs["previous_example"])
    out = {"aligned":[],"up_insertion":[]}
    for key,left in (("aligned",np.eye(3)), ("up_insertion",h)):
        w = np.linalg.solve(m,left@m)
        z0 = (b@w@b.conj())/(scale*heavy)
        out[key+"_z_per_unit_wilson"] = [float(z0.real),float(z0.imag)]
        for k in inputs["real_wilson_tests"]:
            rev = k/scale*left@m@b.conj()
            mat = mass_matrix(m,b,heavy,rev)
            obs = observables(mat)
            schur = np.linalg.det(m)*heavy*(1-k*z0)
            assert abs(np.linalg.det(mat)-schur) < 5e-13
            if k == 0:
                assert abs(obs["light_quartet"]-old["light_weak_quartets"]["01;01"]) < 1e-12
            out[key].append({"wilson":k,**obs})
    # Products of existing dimensionless Yukawa tensors; positive scalar rescalings
    # remove fractions and never affect the real linear span of the nine words.
    ym = m.astype(np.int64)
    x = ym.T@ym
    z = ym.T@np.rint(9*h).astype(np.int64)@ym
    words = [np.linalg.matrix_power(x,i)@z@np.linalg.matrix_power(x,j)
             for i,j in itertools.product(range(3),repeat=2)]
    out["full_real_algebra_exact_basis_determinant"] = str(det_integer(np.stack(words).reshape(9,9)))
    assert int(out["full_real_algebra_exact_basis_determinant"]) != 0
    out["same_left_fixed_positive_heavy"] = right_rotation(m,b,heavy,inputs["right_rotation_alpha"])
    edm = inputs["edm"]
    slope,stat,syst = edm["abs_slope_e_fm"]
    low = (slope-stat-syst)*edm["fm_in_cm"]
    angle_limit = (edm["limit_90_ecm"]+edm["assumed_other_ecm"])/low
    r,i = out["up_insertion_z_per_unit_wilson"]
    tanmax = math.tan(angle_limit)
    kval = inputs["comparison_wilson_small"]
    theta_small = math.atan2(-kval*i,1-kval*r)
    out["conditional_edm"] = {"slope_lower_ecm":low,"theta_limit_rad":angle_limit,
                "real_wilson_interval_near_zero":[-tanmax/(i-r*tanmax), tanmax/(i+r*tanmax)],
                "test_wilson":kval,"theta_test":theta_small,
                "abs_dn_test_central_slope_ecm":abs(theta_small)*slope*edm["fm_in_cm"]}
    return out


def same(x,y):
    if isinstance(x,dict):
        return x.keys()==y.keys() and all(same(x[k],y[k]) for k in x)
    if isinstance(x,list):
        return len(x)==len(y) and all(same(a,b) for a,b in zip(x,y))
    if isinstance(x,float):
        return math.isclose(x,y,rel_tol=1e-10,abs_tol=1e-30 if 0 < abs(y) < 1e-15 else 5e-14)
    return x==y


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write-new",action="store_true")
    mode.add_argument("--check",action="store_true")
    args=parser.parse_args()
    data=evaluate()
    result=ROOT/"results.json"
    if args.write_new:
        with result.open("x",encoding="utf8") as stream:
            json.dump(data,stream,indent=2)
            stream.write("\n")
    else:
        assert same(data,read(result))
    print("1389 OK")
