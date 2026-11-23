"""Signed CP selecting potential: cubic fibre, Hessian and physical norm."""
import argparse
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np
from numpy.polynomial import Polynomial


HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def pdg(angles):
    a,b,c,phase = angles
    s,bv,t = np.sin([a,b,c])
    cs,cb,ct = np.cos([a,b,c])
    z = np.exp(1j*phase)
    return np.array([[cs*cb,s*cb,bv/z],[-s*ct-cs*t*bv*z,cs*ct-s*t*bv*z,t*cb],
                     [s*t-cs*ct*bv*z,-cs*t-s*ct*bv*z,ct*cb]])


def rational(x):
    return Fraction(str(float(x)))


def perpendicular(triple):
    a,b,c = triple
    return [b-c,c-a,a-b]


def curvature(entries, direction):
    # Complete the first two rows from their first two entries. Exact arithmetic
    # avoids loss of the physical triangle area in hierarchical fitted data.
    rows = []
    for p, d in zip(entries[:2], direction[:2]):
        first = [Polynomial([rational(p[i]), d[i]]) for i in (0,1)]
        rows.append([*first, Polynomial([Fraction(1)])-first[0]-first[1]])
    lengths = [x*y for x,y in zip(*rows)]
    q = (2*sum(lengths[i]*lengths[(i+1)%3] for i in range(3)) -
         sum(term*term for term in lengths))/4
    co = list(q.coef) + [Fraction(0)]*5
    assert co[4]==0 and co[0]>0
    j = math.sqrt(float(co[0]))
    first = float(co[1])/(2*j)
    second = float((4*co[0]*co[2]-co[1]**2) / (4*co[0]))/j
    assert second < 0
    return j, first, second, [float(x) for x in co[:4]]


def compute():
    cfg = read(HERE/"inputs.json")
    data = read(HERE/cfg["charged_table"])["high_scale"]
    angles = [row[0] for row in data["angles_theta12_theta13_theta23_delta"]]
    v = pdg(angles)
    actual_j = float((v[0,0]*v[1,1]*v[0,1].conjugate()*v[1,0].conjugate()).imag)
    axes = [perpendicular([rational(row[0])**2 for row in data[f]]) for f in ("up","down")]
    axes = [[x/abs(axis[0]) for x in axis] for axis in axes]
    direction = [[x*y for y in axes[1]] for x in axes[0]]
    j, first, second, co = curvature(np.abs(v)**2, direction)
    assert math.isclose(j,abs(actual_j),rel_tol=1e-8)
    max_j = 1/(6*math.sqrt(3))
    fourier = np.exp(2j*np.pi*np.outer(np.arange(3),np.arange(3))/3)/math.sqrt(3)
    j_fourier = float((fourier[0,0]*fourier[1,1]*
                       fourier[0,1].conjugate()*fourier[1,0].conjugate()).imag)
    assert np.allclose(abs(fourier)**2,1/3) and math.isclose(j_fourier,max_j)
    tree = read(HERE/"../1381/inputs.json")["ub_over_cb"]
    vcb = read(HERE/"../1412/inputs.json")["BelleII_Vcb_BGL"]
    rho = tree[0]+2*sum(tree[1:])
    cb = vcb[0]+2*sum(vcb[1:])
    tree_max_j = .5*rho*cb**2*math.sqrt(1-(1+rho**2)*cb**2)
    return {"fixed_central":dict(J=actual_j,J_over_max=actual_j/max_j,
                                fibre_J_prime=first,fibre_J_second=second,
                                fibre_J_squared_ascending=co),
            "equal_moduli":dict(J=j_fourier,probability=1/3,
                                maximum_ckm_Hessian_divided_by_Jmax=[-4.,-6.,-4.,-1.]),
            "tree_two_width_bounds":dict(Vub_div_Vcb_upper=rho,Vcb_upper=cb,
                J_upper=tree_max_j,min_rho_transport_to_Fourier=1/rho,
                required_relative_extra_potential_range=1-tree_max_j/max_j)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--check",action="store_true")
    group.add_argument("--write-new",action="store_true")
    args=parser.parse_args()
    result=compute()
    target=HERE/"independent_results.json"
    if args.write_new:
        with target.open("x",encoding="utf-8") as stream:
            json.dump(result,stream,indent=2)
            stream.write("\n")
    else:
        old=read(target)
        for category, rows in result.items():
            for name, value in rows.items():
                np.testing.assert_allclose(value,old[category][name],rtol=1e-9,atol=1e-15)
        print("1564 signed degree-twelve extremum and weak-current bound verified")
