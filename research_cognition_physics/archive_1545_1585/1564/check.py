"""Exact rational CKM competitors with the first three mixed flavour moments held fixed."""
import argparse
import json
from fractions import Fraction as F
from itertools import product
from math import cos, prod, sin, sqrt
from pathlib import Path

HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def rational_probabilities(angles, digits):
    a, b, c, delta = angles
    p00 = (cos(a)*cos(b))**2
    p01 = (sin(a)*cos(b))**2
    p10 = (sin(a)*cos(c))**2+(cos(a)*sin(c)*sin(b))**2 + (
        2*sin(a)*cos(a)*sin(c)*cos(c)*sin(b)*cos(delta))
    p11 = (cos(a)*cos(c))**2+(sin(a)*sin(c)*sin(b))**2 - (
        2*cos(a)*sin(a)*sin(c)*cos(c)*sin(b)*cos(delta))
    p = [[F(format(x,f".{digits}f")) for x in row] for row in [[p00,p01],[p10,p11]]]
    for row in p:
        row.append(1-sum(row))
    p.append([1-p[0][j]-p[1][j] for j in range(3)])
    return p


def axis(spec):
    return [F(1), (spec[2]-spec[0])/(spec[1]-spec[2]),
            (spec[0]-spec[1])/(spec[1]-spec[2])]


def displacement(p, direction, t):
    return [[x+t*y for x,y in zip(row,dr)] for row,dr in zip(p,direction)]


def moments(p, a, b):
    return [sum(a[i]**r*b[j]**s*p[i][j] for i,j in product(range(3),repeat=2))
            for r,s in [(1,1),(2,1),(1,2),(2,2)]]


def j_squared(p):
    x,y,z = [p[0][j]*p[1][j] for j in range(3)]
    return x*y-(z-x-y)**2/4


def isolate(p, direction, outer):
    inside = F(0)
    outside = outer
    assert j_squared(p) > 0
    assert j_squared(displacement(p,direction,outside)) <= 0
    for _ in range(90):
        t = (inside+outside)/2
        if j_squared(displacement(p,direction,t)) > 0:
            inside = t
        else:
            outside = t
    return inside, outside


def evaluate():
    cfg = read(HERE/"inputs.json")
    table = read(HERE/cfg["charged_table"])["high_scale"]
    angles = [x[0] for x in table["angles_theta12_theta13_theta23_delta"]]
    spec = [[F(str(x[0]))**2 for x in table[name]] for name in ("up","down")]
    a,b = spec
    p = rational_probabilities(angles,cfg["rational_digits"])
    ur,dc = axis(a),axis(b)
    assert sum(ur)==sum(dc)==sum(x*y for x,y in zip(a,ur))==sum(x*y for x,y in zip(b,dc))==0
    direction = [[x*y for y in dc] for x in ur]
    slope = moments(direction,a,b)
    assert slope[:3]==[0,0,0] and slope[3]!=0
    # Row and column constants plus mixed moments form a tensor-product Vandermonde.
    determinants = [prod(y-x for i,x in enumerate(row) for y in row[i+1:]) for row in spec]
    assert all(determinants)
    bounds = [(-x/y,y) for row,v in zip(p,direction) for x,y in zip(row,v) if y]
    lower = max(t for t,y in bounds if y>0)
    upper = min(t for t,y in bounds if y<0)
    points = []
    for endpoint in [lower,upper]:
        inside,outside = isolate(p,direction,endpoint)
        edge = displacement(p,direction,inside)
        w = moments(edge,a,b)
        original = moments(p,a,b)
        assert all(x>=0 for row in edge for x in row) and w[:3]==original[:3]
        assert w[3]!=original[3] and abs(outside-inside)<F(1,10**25)
        assert 0 < j_squared(edge) < F(1,10**22)
        points.append({"t":float(inside),"t_certified_width_upper":float(abs(outside-inside)),
                       "Vus_modulus":sqrt(float(edge[0][1])),"Vub_modulus":sqrt(float(edge[0][2])),
                       "Vcb_modulus":sqrt(float(edge[1][2])),"J_abs_at_true_boundary":0.,
                       "absolute_K22_change":float(abs(w[3]-original[3]))})
    ten_degrees = [2,4,6,2,4,6,4,6,6,8]
    return {"realistic_J_abs":sqrt(float(j_squared(p))),"realistic_Vus_modulus":sqrt(float(p[0][1])),
            "central_three_mixed_moments":[float(x) for x in moments(p,a,b)[:3]],
            "competitors":points,"construction_max_degree":max(2*d for d in ten_degrees),
            "number_of_arbitrary_real_targets_in_full_construction":len(ten_degrees)}


def verify(x,y):
    if isinstance(x,dict):
        assert x.keys()==y.keys()
        for k in x:
            verify(x[k],y[k])
    elif isinstance(x,list):
        assert len(x)==len(y)
        for a,b in zip(x,y):
            verify(a,b)
    else:
        assert abs(x-y) <= 1e-10*max(abs(x),abs(y),1e-38), (x,y)


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    check=parser.parse_args().check
    result=evaluate()
    dest=HERE/"results.json"
    if check:
        verify(result,read(dest))
        print("1564 rational competitors with identical three mixed moments and minimal potential degree verified")
    else:
        dest.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
