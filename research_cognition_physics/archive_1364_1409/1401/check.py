"""Best simultaneous left attenuation for two full-rank quark Yukawa maps."""
import argparse
import importlib.util
import json
from math import cos, sqrt
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
INPUT = json.loads((ROOT/"inputs.json").read_text(encoding="utf-8-sig"))
spec = importlib.util.spec_from_file_location("quark_measurement_1397",ROOT/INPUT["shared_python"])
prior = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prior)


def condition(x):
    sv = np.linalg.svd(x,compute_uv=False)
    return float(sv[0]/sv[-1])


def split(yu,yd,p):
    ratio = np.linalg.solve(yu,yd)
    U,sv,Vh = np.linalg.svd(ratio)
    # G maps from the interim singular coordinates to the common weak space.
    G = yu @ (U * sv**p)
    A,ds,Bh = np.linalg.svd(G)
    F = (A * ds) @ A.conj().T
    W = A @ Bh
    Cu = W @ (sv[:,None]**(-p)*U.conj().T)
    Cd = W @ (sv[:,None]**(1-p)*Vh)
    norm = ds[0]
    F,Cu,Cd = F/norm,Cu*norm,Cd*norm
    assert np.allclose(F@Cu,yu,atol=1e-13,rtol=1e-10)
    assert np.allclose(F@Cd,yd,atol=1e-13,rtol=1e-10)
    assert np.min(np.linalg.eigvalsh(F)) > 0
    k = sv[0]/sv[-1]
    assert np.isclose(condition(Cu),k**p,rtol=2e-10)
    assert np.isclose(condition(Cd),k**(1-p),rtol=2e-10)
    return [condition(Cu),condition(Cd)]


def bounds(spectrum_u,spectrum_d,angles,k,radian):
    def endpoints(x):
        center,width = np.array(x).T
        return center-k*width,center+k*width
    ul,uh = endpoints(spectrum_u)
    dl,dh = endpoints(spectrum_d)
    al,ah = endpoints(angles)
    if radian:
        al[:3],ah[:3]=np.sin(al[:3]),np.sin(ah[:3])
    sl,sh=al[:3],ah[:3]
    cl,ch=np.sqrt(1-sh**2),np.sqrt(1-sl**2)
    assert 0<al[3]<np.pi/2 and ah[3]<np.pi/2
    row_u_lower=np.array([cl[0]*cl[1],sl[0]*cl[1],sl[1]])
    # Vtd = a - b exp(i delta): positive real separation plus its phase contribution.
    a0,a1=sl[0]*sl[2],sh[0]*sh[2]
    b0,b1=cl[0]*cl[2]*sl[1],ch[0]*ch[2]*sh[1]
    difference=max(0.,a0-b1,b0-a1)
    td_min=sqrt(difference**2 + 2*a0*b0*(1-cos(al[3])))
    ts_min=cl[0]*sl[2]+sl[0]*cl[2]*sl[1]*cos(ah[3])
    row_t_lower=np.array([td_min,ts_min,cl[1]*cl[2]])
    row_norm=np.linalg.norm(row_u_lower*dl)/uh[0]
    inverse_col_norm=ul[2]*np.linalg.norm(row_t_lower/dh)
    mass_only=sqrt((ul[2]/uh[0])/(dh[2]/dl[0]))
    # Triangle upper bounds for CKM magnitudes, Frobenius norm upper-bounds spectral norm.
    s12,s13,s23=sh
    c12,c13,c23=ch
    vmax=np.minimum(1.,np.array([[c12*c13,s12*c13,s13],
        [s12*c23+c12*s23*s13,c12*c23+s12*s23*s13,s23*c13],
        [s12*s23+c12*c23*s13,c12*s23+s12*c23*s13,c23*c13]]))
    norm_upper=np.linalg.norm(vmax*dh[None,:]/ul[:,None])
    inv_upper=np.linalg.norm(vmax*uh[:,None]/dl[None,:])
    detlo=np.prod(dl)/np.prod(uh)
    dethi=np.prod(dh)/np.prod(ul)
    intervals=[[row_norm,norm_upper],[detlo*inverse_col_norm/norm_upper,dethi*inv_upper/row_norm],
               [1/inv_upper,1/inverse_col_norm]]
    assert intervals[1][1]<intervals[0][0] and intervals[2][1]<intervals[1][0]
    return {"mass_only_minimax_floor":mass_only,"row_norm_lower":row_norm,
            "inverse_top_column_lower":inverse_col_norm,
            "with_mixing_minimax_floor":sqrt(row_norm*inverse_col_norm),
            "min_product_condition":row_norm*inverse_col_norm,"ratio_singular_value_intervals":intervals}


def single_direction(yu, yd):
    left, spectrum, right_h = np.linalg.svd(np.linalg.solve(yu, yd))
    up_ratio = spectrum[0] / spectrum[1]
    down_ratio = spectrum[1] / spectrum[2]
    up = np.diag([1 / up_ratio, 1., 1.]) @ left.conj().T
    down = spectrum[1] * np.diag([1., 1., 1 / down_ratio]) @ right_h
    common = yu @ left @ np.diag([up_ratio, 1., 1.])
    polar_left, gains, polar_right_h = np.linalg.svd(common)
    positive = (polar_left * gains) @ polar_left.conj().T
    basis = polar_left @ polar_right_h
    assert np.min(np.linalg.eigvalsh(positive)) > 0
    assert np.allclose(positive @ basis @ up, yu)
    assert np.allclose(positive @ basis @ down, yd)
    for bare, expected in ((up, up_ratio), (down, down_ratio)):
        singular = np.linalg.svd(bare, compute_uv=False)
        assert np.isclose(singular[0], singular[1])
        assert np.isclose(condition(bare), expected)
    return {"up_minimum_condition": up_ratio, "down_minimum_condition": down_ratio}


def calculate():
    results={}
    for name,up,down,angles,_g,radian in prior.data():
        u,d,a=(np.array(z)[:,0] for z in (up,down,angles))
        s=np.sin(a[:3]) if radian else a[:3]
        V=prior.mix(s,a[3])
        yu,yd=np.diag(u),V*d
        ratio=np.linalg.solve(yu,yd)
        sv=np.linalg.svd(ratio,compute_uv=False)
        value=condition(ratio)
        gain_equivariant=condition(np.linalg.solve(2*yu,0.7*yd))
        assert np.isclose(gain_equivariant,value)
        compare=condition(np.diag(d/u))
        boxes={str(k):bounds(up,down,angles,k,radian) for k in INPUT["widths"]}
        for box in boxes.values():
            assert box["with_mixing_minimax_floor"]<sqrt(value)
            assert all(lo<=x<=hi for x,(lo,hi) in zip(sv,box["ratio_singular_value_intervals"]))
        results[name]={"relative_singular_values":sv.tolist(),"relative_condition":value,"single_direction":single_direction(yu,yd),
            "minimum_larger_department":sqrt(value),"equal_left_eigenframes_minimum":sqrt(compare),
            "Pareto":[[p,*split(yu,yd,p)] for p in INPUT["pareto_p"]],"boxes":boxes}
    eps=INPUT["anarchy_comparison_epsilon"]
    results["bounded_entries_can_be_ill_conditioned"]=condition(np.ones((3,3))+eps*np.eye(3))
    return results


def close(x,y):
    if isinstance(x,dict):
        return x.keys()==y.keys() and all(close(v,y[k]) for k,v in x.items())
    return bool(np.allclose(x,y,rtol=1e-9,atol=1e-12))


if __name__=="__main__":
    flags=argparse.ArgumentParser()
    group=flags.add_mutually_exclusive_group(required=True)
    group.add_argument("--check",action="store_true")
    group.add_argument("--write-new",action="store_true")
    args=flags.parse_args()
    answer=calculate()
    path=ROOT/"results.json"
    if args.write_new:
        with path.open("x",encoding="utf-8") as file:
            json.dump(answer,file,indent=2,allow_nan=False)
            file.write("\n")
    else:
        assert close(answer,json.loads(path.read_text(encoding="utf-8-sig")))
        print("1401 optimal common attenuation, reconstruction and whole-box lower bounds verified")
