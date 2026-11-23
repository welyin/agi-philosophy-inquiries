"""Conditional weak-CP tests for a common single-projector flavour potential."""
import json
import math
from pathlib import Path
import sys

import numpy as np

BASE = Path(__file__).resolve().parent

def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))

def driver_distance(j_floor, upper_moduli):
    # Monotone four-factor comparison certificate for entrywise distance.
    def fourth_difference(eps):
        return math.fsum(math.prod(entry if bit & (1 << k) else eps
                                   for k,entry in enumerate(upper_moduli))
                         for bit in range(15))
    lo, hi = 0., .3
    for _ in range(55):
        middle = (lo+hi)/2
        if fourth_difference(middle)<j_floor:
            lo = middle
        else:
            hi = middle
    return math.floor(lo*1e9)/1e9

def quartet(v):
    return abs((v[0,1]*v[1,2]*v[0,2].conjugate()*v[1,1].conjugate()).imag)

def calculate():
    cfg=read(BASE/"inputs.json")
    old=read(BASE/cfg["moduli"])
    low,high = (np.array(old["tree_only_moduli_outer_"+k]) for k in ("low","high"))
    phase=read(BASE/cfg["charged_b"])["chargedB"]
    ends=(phase["central_degrees"]-2*phase["lower_error"]-phase["additional_degrees"],
          phase["central_degrees"]+2*phase["upper_error"]+phase["additional_degrees"])
    assert 0 < ends[0] < ends[1] < 90
    sine_floor=math.floor(math.sin(math.radians(ends[0]))*1e10)/1e10
    jmin=math.floor(np.prod(low[[0,0,1,1],[0,2,0,2]])*sine_floor*1e12)/1e12
    jmax=math.ceil(np.prod(high[[0,0,1,1],[0,2,0,2]])*1e12)/1e12
    moduli=high[[0,1,0,1],[1,2,2,1]]
    eps=driver_distance(jmin,moduli)
    sine=eps*math.sqrt(1-eps**2/4)
    driver_floor=math.floor(sine/(1+sine)*1e9)/1e9
    # A generic noncommuting A and rank-one P; the phase convention of B follows X.
    A=np.diag([.01,.23,1.])
    n=np.array([1+2j,-2+1j,1-3j])
    n/=np.linalg.norm(n)
    P=np.outer(n,n.conj())
    x=.4*A+P+.3*(A@P+P@A)-.1*A@P@A
    eig,v=np.linalg.eigh(x)
    assert min(np.diff(eig))>.01 and quartet(v)<2e-14
    bspec=np.array([.03,.2,.92])
    B=(v*bspec)@v.conj().T
    assert np.linalg.norm(B@x-x@B)<1e-13
    F=np.exp(2j*math.pi*np.outer(np.arange(3),np.arange(3))/3)/math.sqrt(3)
    Bfour=(F*bspec)@F.conj().T
    comm=A@Bfour-Bfour@A
    spectral_prod=lambda vals: math.prod(vals[i]-vals[j] for i in range(3) for j in range(i))
    jfour=abs(np.trace(comm@comm@comm).imag/(6*spectral_prod(A.diagonal())*spectral_prod(bspec)))
    assert abs(jfour-quartet(F))<2e-14 and abs(jfour-1/(6*math.sqrt(3)))<2e-14
    # A doubly degenerate real driver allows both a CP-conserving and violating basis,
    # but only if its remaining real line also has nonzero overlaps with all mass lines.
    xflat=np.ones((3,3))/3
    assert np.linalg.norm(Bfour@xflat-xflat@Bfour)<1e-13
    # The unique high-eigenvalue of X must meet the heaviest mass for a global minimum.
    Bflat=(F*np.array([bspec[-1],bspec[0],bspec[1]]))@F.conj().T
    _,real_basis=np.linalg.eigh(xflat)
    Breal=(real_basis*bspec)@real_basis.T
    assert abs(np.trace(Bflat@xflat)-np.trace(Breal@xflat))<1e-14
    assert quartet(real_basis)<1e-14
    return dict(tree_gamma_degrees=list(ends), mirrored_by_180=True,
        tree_abs_J_lower=jmin,tree_abs_J_upper=jmax,
        minimum_complex_CKM_entry_distance=eps,
        necessary_complex_driver_over_minimum_gap=driver_floor,
        determinant_squared_competitor_max_abs_J=jfour,
        competitor_us=round(abs(F[0,1]),12),
        tree_upper_us=high[0,1], extra_flavour_reference="one rank-one projector")

if __name__=="__main__":
    result=calculate()
    path=BASE/"results.json"
    if sys.argv[1:]==["--write-new"]:
        with path.open("x",encoding="utf-8") as stream:
            json.dump(result,stream,indent=2)
            stream.write("\n")
    elif sys.argv[1:] in ([],["--check"]):
        assert result==read(path)
        print("1501: real-driver minima, CKM distance, gap floor and nonlinear CP test passed.")
    else:
        raise SystemExit("Expected --check or --write-new.")
