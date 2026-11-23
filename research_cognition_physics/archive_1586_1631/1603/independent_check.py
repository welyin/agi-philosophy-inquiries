"""Exact rank-one perturbations and a finite vectorlike SU(5) threshold."""
import argparse
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "independent_results.json"


def unit(vector):
    return np.array(vector, dtype=complex)/np.linalg.norm(vector)


def threshold(z, uv, heavy, line):
    norm = math.hypot(1, z)
    canonical = np.eye(3) + (1/norm-1)*np.outer(line, line.conj())
    embedding = np.vstack([canonical, -(z/norm)*line.conj()])
    mass_row = np.append(z*line.conj(), 1)
    np.testing.assert_allclose(mass_row @ embedding, 0, atol=1e-14)
    np.testing.assert_allclose(embedding.conj().T @ embedding, np.eye(3), atol=1e-14)
    return np.column_stack((uv,heavy)) @ embedding


def constructed_example():
    d = np.diag([0.0000065, 0.000129, 0.00606])
    b = unit([3, 1j, 0.12])
    v = unit([0.4j, 1, .8j])
    a1 = np.array([0.0001j, .002, .004j])
    a2 = np.array([.0004, .002j, .001])
    k, t = 1.13, .10
    e = k*d + np.outer(a1, b.conj()) + np.outer(a2, v.conj())
    zd, zl = 2.5*t, -2.5*t
    inv_c = math.hypot(1, zd)
    # Solve the inverse light matching for ONE heavy-Higgs column common to both contractions.
    heavy = (a2 + (inv_c-1)*(e-k*d)@v) / (k*zd-zl)
    def undo(light, z):
        return light + np.outer((math.hypot(1,z)-1)*(light@v)+z*heavy, v.conj())
    raw_e, raw_d = undo(e, zl), undo(d, zd)
    np.testing.assert_allclose(raw_e-k*raw_d, np.outer(a1,b.conj()), atol=3e-15)
    np.testing.assert_allclose(threshold(zd,raw_d,heavy,v),d,atol=3e-15)
    np.testing.assert_allclose(threshold(zl,raw_e,heavy,v),e,atol=3e-15)
    return {"equal_mass_over_M":inv_c, "second_source_norm":float(np.linalg.norm(a2)),
            "single_shared_heavy_yukawa_norm":float(np.linalg.norm(heavy))}


def signed_update(base, target):
    base, target = np.array(base,dtype=float), np.array(target,dtype=float)
    signs = np.where(target >= base, 1., -1.)
    x = np.sort(signs*base)
    z = np.sort(signs*target)
    assert x[0] <= z[0] <= x[1] <= z[1]
    squares = np.array([-np.prod(a-z)/(a-x[1-i]) for i,a in enumerate(x)])
    v = np.sqrt(np.maximum(squares,0))
    updated = np.diag(x)+np.outer(v,v)
    np.testing.assert_allclose(np.linalg.svd(updated, compute_uv=False)[::-1],target,atol=3e-14)
    norm = np.linalg.norm(np.outer(v,v),ord=2)
    np.testing.assert_allclose(norm, np.abs(base-target).sum(), atol=3e-14)
    return norm


def invert(f):
    left, right = 0.0, 1.0
    assert f(left) < 0 < f(right)
    for _ in range(80):
        middle = (left+right)/2
        if f(middle) < 0:
            left = middle
        else:
            right = middle
    return (left+right)/2


def compute():
    def inputs(idx):
        return json.loads((HERE.parent/str(idx)/"inputs.json").read_text(encoding="utf8"))
    masses=np.asarray(inputs(1402)["scales"]["1e16"])
    lo,hi=masses[:,0]-2*masses[:,1],masses[:,0]+2*masses[:,1]
    nu, atmospheric = inputs(1370), inputs(1580)["theta23_sin2"]
    no, io = nu["NO"],nu["IO"]
    pe=(1-no["s12"][2])*(1-no["s13"][2])
    ptno=(math.sqrt(no["s12"][2]*atmospheric["NO"][2])+
          math.sqrt((1-no["s12"][1])*(1-atmospheric["NO"][1])*no["s13"][2]))**2
    ptio=(1-atmospheric["IO"][1])*(1-io["s13"][1])
    # Round out the previously checked high-scale full box.
    pb_max=0.03513
    ratio=1/math.sqrt(1-pb_max)
    second_min=lo[2]/ratio
    delta=math.sqrt(pe)*lo[4] - hi[5]*(hi[1]/lo[2])*ratio
    def finite_factors(t):
        return math.hypot(1,3*t)-1, math.hypot(1,2*t)-1
    def source_cap(t):
        a,b=finite_factors(t)
        return a*hi[5]+ratio*b*(hi[5]+delta)
    t_floor=invert(lambda t:source_cap(t)-delta)
    def pole_ratio(t):
        return math.hypot(1,3*t)/math.hypot(1,2*t)
    t_test=math.sqrt((1.02**2-1)/(9-4*1.02**2))
    slope_at_t01=.01*(3+2*(hi[5]+delta)/second_min)
    fixed={}
    for order, pt, min_e in [("NO",ptno,pe),("IO",ptio,io["s13"][1])]:
        r=max(0.,lo[4]*math.sqrt(min_e)-hi[1])+max(0.,lo[5]*math.sqrt(1-pt)-hi[2])
        tf=invert(lambda t:np.dot(finite_factors(t),[hi[5],hi[2]])-r)
        fixed[order]={"minimum_correction":r,"mixing_parameter_lower":tf,
                      "heavy_doublet_over_triplet_lower":pole_ratio(tf)}
    minimum_log=math.log(pole_ratio(t_floor))
    return {"signed_2_by_2_updates": [signed_update([2,5], z) for z in
                           [[3,7],[1,4],[1,7],[3,4],[2,7],[5,7],[2,5]]],
            "generic_free_k_NO":{"light_source_lower":delta,"t_lower":t_floor,
                "heavy_doublet_over_triplet_lower":pole_ratio(t_floor),
                "abs_inverse_alpha_12_threshold_lower_relative_to_common_at_ML":(4/15)*minimum_log/(2*math.pi),
                "abs_inverse_alpha_23_threshold_lower_relative_to_common_at_ML":(2/3)*minimum_log/(2*math.pi),
                "ray_angle_lower_deg_if_heavy_ratio_is_1_02":math.degrees(math.asin(delta/source_cap(t_test))),
                "min_heavy_yukawa_if_t_equals_0_01":(delta-source_cap(.01))/slope_at_t01},
            "fixed_k1":fixed,"bare_mixing_competitor":constructed_example()}


def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:
            compare(a[k],b[k])
    else:
        np.testing.assert_allclose(a,b,rtol=1e-11,atol=1e-14)


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    mode=cli.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check",action="store_true")
    mode.add_argument("--write",action="store_true")
    args=cli.parse_args()
    result=compute()
    if args.write:
        with OUT.open("x",encoding="utf8") as stream:
            json.dump(result,stream,indent=2)
            stream.write("\n")
    else:
        compare(result,json.loads(OUT.read_text(encoding="utf8")))
        print("1603 independent checks passed")
