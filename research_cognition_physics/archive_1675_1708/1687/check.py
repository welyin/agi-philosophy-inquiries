"""Normalized same-decay Higgs production for a shared top and colored mass source."""
import argparse
import json
import math
from pathlib import Path
import runpy

import numpy as np


HERE = Path(__file__).resolve().parent
LOOPS = runpy.run_path(str(HERE.parent/"1373/check.py"))
DOWN = runpy.run_path(str(HERE.parent/"1653/check.py"))


def read(path):
    return json.loads((HERE/path).read_bytes())


def couplings(cfg):
    power = np.array(read(cfg["down_powers"])["powers_A_B_C_E_F_R_D"])
    real = np.array(read(cfg["down_paths"])["quark_path_amplitudes"])
    paths = real*np.exp(1j*math.pi*cfg["phase_rad_over_pi"]*power)
    mass = DOWN["assemble"](paths)
    vertex = DOWN["assemble"](power*paths)
    left,sv,vh = np.linalg.svd(mass)
    slopes = (np.diag(left.conj().T@vertex@vh.conj().T)/sv)[::-1]
    np.testing.assert_allclose(sum(slopes),4,atol=3e-13)
    return slopes


def form(mh,mass):
    z=(mh/(2*mass))**2
    if z<0.002:
        # Stable converging series for the two very heavy masses. Error <2e-12 here.
        return (1/3+8*z/45+4*z*z/35+128*z**3/1575+128*z**4/2079,
                4/3+14*z/45+8*z*z/63+104*z**3/1575+2048*z**4/51975,
                2+2*z/3+16*z*z/45+8*z**3/35+256*z**4/1575)
    f = np.arcsin(np.sqrt(z))**2 if z<=1 else -(2*np.arccosh(np.sqrt(z))-1j*np.pi)**2/4
    return (f-z)/z**2, 2*(z+(z-1)*f)/z**2, 2*f/z


def coefficients(cfg,slopes,heavy):
    mh=cfg["higgs_GeV"]
    old=LOOPS["read"]("inputs.json")
    _,_,fermions=LOOPS["common_inputs"](old)
    light=[form(mh,mass) for mass,_,_ in fermions[:6]]
    ds=np.r_[cfg["up_exponents_uct"],slopes]
    g=sum(scalar for _,scalar,_ in light)
    even=sum(scalar*d.real for (_,scalar,_),d in zip(light,ds))
    odd=sum(pseudo*d.imag for (_,_,pseudo),d in zip(light,ds))
    chain=read(cfg["chain"])
    shift=0.
    for mass,p in zip(chain["conditional_masses_TeV_descending"],chain["individual_log_mass_radial_slopes"]):
        spin0,spin_half,_=form(mh,1000*mass)
        shift+=p*(spin_half+heavy*spin0)
    gain=(even+shift)/g
    coeff=np.array([1,2*gain.real,abs(gain)**2+abs(odd/g)**2])
    bb=np.array([1,2*slopes[-1].real,abs(slopes[-1])**2])
    return float(shift),coeff,bb,even,odd,g


def windows(atlas,n):
    out={}
    for key in ["top","vector","gluon"]:
        value,minus,plus,sm=atlas[key]
        out[key]=np.array([max(0.,value-n*minus),value+n*plus])/sm
    return out


def roots(poly):
    p=np.trim_zeros(np.asarray(poly,dtype=float),"b")
    return [float(x.real) for x in np.roots(p[::-1]) if abs(x.imag)<1e-9] if len(p)>1 else []


def project(box,gpoly,bb,bounds,br_times_mixing=None):
    # For a fixed z, one common A>0 multiplies all three rates. The unknown
    # small bbH fraction lets the gluon group's rate traverse a segment.
    v,t,g=(box[k] for k in ["vector","top","gluon"])
    top=np.array([1.,4.,4.])
    flat=np.array([1.,0.,0.])
    group=[(1-b)*gpoly+b*bb for b in bounds]
    dividers=roots(group[0]-group[1])+roots(top)
    dividers+=roots(t[0]*flat-v[1]*top)+roots(t[1]*flat-v[0]*top)
    for gp in group:
        dividers+=roots(g[0]*flat-v[1]*gp)+roots(g[1]*flat-v[0]*gp)
        dividers+=roots(t[0]*gp-g[1]*top)+roots(t[1]*gp-g[0]*top)
    if br_times_mixing is not None:
        for measured,poly in [(v,flat),(t,top),*((g,p) for p in group)]:
            dividers+=roots(measured[0]*br_times_mixing-poly)
    def allowed(z):
        powers=np.array([1,z,z*z])
        rt=top@powers
        gg=sorted(float(poly@powers) for poly in group)
        if min(gg[0],rt)<=0:
            return False
        maximum=min(v[1],t[1]/rt,g[1]/gg[0])
        if br_times_mixing is not None:
            maximum=min(maximum,1/(br_times_mixing@powers))
        return max(v[0],t[0]/rt,g[0]/gg[1])<=maximum+1e-12
    points=sorted(set(round(y,13) for y in dividers))
    intervals=[]
    assert not allowed(points[0]-1000) and not allowed(points[-1]+1000)
    for a,b in zip(points,points[1:]):
        if allowed((a+b)/2):
            if intervals and intervals[-1][1]==a:
                intervals[-1][1]=b
            else:
                intervals.append([a,b])
    for a in points:
        if allowed(a) and not any(x-1e-11<=a<=y+1e-11 for x,y in intervals):
            intervals.append([a,a])
    return sorted(intervals)


def compute():
    cfg=read("inputs.json")
    sl=couplings(cfg)
    cap=cfg["bbH_fraction"]
    out={"mass_basis_radial_down_dsb":[[d.real,d.imag] for d in sl]}
    uv=read("../1685/inputs.json")
    eps=uv["scalar_magnitude_over_cutoff"]
    l2,l3=uv["lambda2_lambda3_left_right"][:2]
    kinetic=1+(2*l2*eps)**2+(3*l2*l3*eps**2)**2
    effective_f=math.sqrt(2*kinetic)*eps*uv["cutoff_TeV"]*1000
    v=LOOPS["read"]("inputs.json")["v_MZ_GeV"]
    br=cfg["SM_diphoton_branching_fraction"]
    br_times_mixing=br*np.array([1,0,(effective_f/v)**2])
    out["valley_metric"]=kinetic
    out["radial_scale_TeV"]=effective_f/1000
    out["unit_branching_VBF_bound_on_abs_z"]={
        str(n):v/effective_f*math.sqrt(1/(br*windows(cfg["atlas"],n)["vector"][0])-1)
        for n in (1,2)}
    for label,n_scalars in [("dirac",0),("unsplit_new_chiral_pairs",2)]:
        shift,gg,bb,even,odd,sm=coefficients(cfg,sl,n_scalars)
        cases={}
        fixed={}
        for n in (1,2):
            box=windows(cfg["atlas"],n)
            feasible=project(box,gg,bb,cap)
            assert any(a<0<b for a,b in feasible)
            cases[str(n)]=feasible
            fixed[str(n)]=project(box,gg,bb,cap,br_times_mixing)
        out[label]={"heavy_even_amplitude_derivative":shift,"gg_polynomial":gg.tolist(),
                    "three_rate_common_normalization_windows":cases,
                    "actual_chain_scale_and_branching_at_most_one":fixed}
    out["standard_gluon"]=[sm.real,sm.imag]
    out["ordinary_even_derivative"]=[even.real,even.imag]
    out["ordinary_odd_derivative"]=[odd.real,odd.imag]
    out["bb_polynomial"]=bb.tolist()
    # Independent production-only ratio bounds with no heavy-loop hypothesis.
    bounds={}
    for n in (1,2):
        box=windows(cfg["atlas"],n)
        lo=math.sqrt(box["top"][0]/box["vector"][1])
        hi=math.sqrt(box["top"][1]/box["vector"][0])
        bounds[str(n)]=[[(-hi-1)/2,(-lo-1)/2],[(lo-1)/2,(hi-1)/2]]
    out["top_vector_only_z"]=bounds
    return out


def rounded(x):
    if isinstance(x,dict):
        return {k:rounded(v) for k,v in x.items()}
    if isinstance(x,list):
        return [rounded(v) for v in x]
    return round(float(x),10)


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    group=cli.add_mutually_exclusive_group(required=True)
    group.add_argument("--check",action="store_true")
    group.add_argument("--write",action="store_true")
    args=cli.parse_args()
    out=rounded(compute())
    if args.write:
        (HERE/"results.json").write_text(json.dumps(out,indent=2)+"\n",encoding="utf-8")
    else:
        assert out==read("results.json")
        print("PASS: shared colored spectrum, full complex down couplings and three ATLAS production windows")
