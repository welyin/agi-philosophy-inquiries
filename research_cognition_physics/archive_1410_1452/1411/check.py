"""Signed path budgets for a cyclic three-state Yukawa operation."""
import argparse
import itertools
import json
from math import sqrt
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
OMEGA=np.exp(2j*np.pi/3)
ROOTS=OMEGA**np.arange(3)


def read(name):
    return json.loads((HERE/name).read_text(encoding="utf8"))


def decompose(vals):
    return np.fft.fft(vals)/3


def shape(pairs):
    s,err=np.array(pairs).T
    assert (np.diff(s)>0).all()
    gap=s[2]-s[0]-s[1]
    signed=np.array([-s[0],-s[1],s[2]])
    vals=signed*ROOTS**2
    coeff=decompose(vals)
    assert np.allclose(np.abs(vals),s) and np.isclose(abs(coeff[2]),gap/3)
    assert np.allclose(np.polyval(coeff[::-1],ROOTS),vals)
    roots=np.roots(coeff[::-1])
    recon=coeff[2]*(ROOTS-roots[0])*(ROOTS-roots[1])
    assert np.allclose(recon,vals,rtol=2e-9,atol=1e-14)
    lo=1-(s[0]+s[1]+2*(err[0]+err[1]))/(s[2]-2*err[2])
    assert lo>0
    return dict(missing_path_min=gap/3,minimum_over_largest=gap/(3*s[2]),
                two_width_min_fraction_of_largest=lo/3,
                two_step_max_reconstruction_error=float(max(abs(recon-vals))))


def mixing(angles,radian):
    val=np.array(angles)[:,0]
    s12,s13,s23=np.sin(val[:3]) if radian else val[:3]
    a,b,c=np.sqrt(1-np.square([s12,s13,s23]))
    z=np.exp(1j*val[3])
    return np.array([[a*b,s12*b,s13/z],[-s12*c-a*s23*s13*z,a*c-s12*s23*s13*z,s23*b],
                     [s12*s23-a*c*s13*z,-a*s23-s12*c*s13*z,c*b]])


def common_cube(ds,angles,radian):
    d,derr=np.array(ds).T
    eig=(d/d[2])**2
    V=mixing(angles,radian)
    H=(V*eig)@V.conj().T
    nearest=np.diag(np.diag(H))
    free=float(np.linalg.norm(H-nearest))
    fixed=min(float(np.linalg.norm(H-np.diag(p))) for p in itertools.permutations(eig))
    central=np.array(angles)[:3,0]
    width=np.array(angles)[:3,1]
    smin,smax=central-2*width,central+2*width
    if radian:
        smin,smax=np.sin(smin),np.sin(smax)
    ratios=(d[:2]+2*derr[:2])/(d[-1]-2*derr[-1])
    bound=sqrt(2)*(smin[2]*sqrt(1-smax[2]**2)*(1-smax[1]**2)-.5*(ratios@ratios))
    assert 0<bound<free<=fixed+1e-12
    assert np.isclose(np.linalg.norm(np.diag(ROOTS)@H-H@np.diag(ROOTS)),sqrt(3)*free)
    return dict(distance_unfixed_spectrum=free,distance_fixed_spectrum=fixed,
                rectangular_lower_bound=bound)


def calculate():
    cfg=read("inputs.json")
    old=read(cfg["low_spectra"])
    high=read(cfg["high_quarks_and_angles"])["high_scale"]
    high_lep=read(cfg["high_leptons"])["scales"]["1e16"][3:]
    result={}
    for name,up,down,lep,ang,rad in [
        ("MZ",old["yukawa"]["up"],old["yukawa"]["down"],old["yukawa"]["charged_leptons"],old["quark_sines_and_phase"],False),
        ("high",high["up"],high["down"],high_lep,high["angles_theta12_theta13_theta23_delta"],True)]:
        result[name]={"three_spectra":{k:shape(v) for k,v in [("up",up),("down",down),("leptons",lep)]},
                      "down_if_up_fixes_exact_cycle":common_cube(down,ang,rad)}
    # A single-step admissible triangle, including nondegenerate lengths.
    example=abs(1+(.27+.38j)*ROOTS)
    assert 2*max(example)<sum(example)
    return result


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    result=calculate()
    if p.parse_args().check:
        def compare(a,b):
            if isinstance(a,dict):
                assert a.keys()==b.keys()
                for k in a: compare(a[k],b[k])
            else:
                assert np.isclose(a,b,rtol=1e-9,atol=1e-13)
        compare(result,read("results.json"))
        print("1411 cyclic path and common-quark-basis certificates verified")
    else:
        print(json.dumps(result,indent=2))
