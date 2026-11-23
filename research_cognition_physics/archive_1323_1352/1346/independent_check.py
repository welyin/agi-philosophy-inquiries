"""Independent set inversion. Mesh cells use integer outward tests, not point sampling."""
from fractions import Fraction as F
import json
import math
from pathlib import Path

L=(F(26,45),F(19,25))
H=(F(29,28),F(91,75))
GRID=10000


def value(n,d,q):
    return n*q*q+d*(1-q)**2


def intersect_cell(n,d,k,box):
    # q lies between k/grid and (k+1)/grid. Bound both squared terms separately.
    loq=0 if k<=0<=k+1 else min(k*k,(k+1)**2)
    hiq=max(k*k,(k+1)**2)
    lod=(GRID-k-1)**2
    hid=(GRID-k)**2
    m=n*loq+d*lod
    M=n*hiq+d*hid
    return (M*box[0].denominator>=box[0].numerator*GRID**2
            and m*box[1].denominator<=box[1].numerator*GRID**2)


def windows(a,d):
    hits=[]
    for k in range(-7800,7800):
        if intersect_cell(2-a,d,k,L) and intersect_cell(2,d,k,H):
            if hits and hits[-1][1]==k: hits[-1][1]+=1
            else: hits.append([k,k+1])
    return [[x/GRID,y/GRID] for x,y in hits]


def finite_roots(a,d):
    bounds={-1.,1.}
    for n,box in ((2-a,L),(2,H)):
        if n+d==0: continue
        for val in box:
            disc=d*d-(n+d)*(d-val)
            if disc>=0:
                for s in (-1,1):
                    bounds.add((d+s*math.sqrt(disc))/(n+d))
    seq=sorted(bounds)
    return [[x,y] for x,y in zip(seq,seq[1:]) if
            L[0]<=value(2-a,d,(x+y)/2)<=L[1] and H[0]<=value(2,d,(x+y)/2)<=H[1]]


def run():
    numer=F(6741,10000)
    err=F(27,10000)
    Umin=(numer-err)/(1-numer+err)/F(11,10)
    Umax=(numer+err)/(1-numer-err)
    assert 1<Umin<2<Umax<3
    # Bound all integer d; M<=91/75 gives |q|<.78
    assert H[1]<2*F(78,100)**2
    assert L[1]/(1-F(78,100))**2<16
    witnesses={(1,0):[F(-77,100),F(77,100)],(1,1):[F(72,100)],
               (1,2):[F(2,3)],(2,2):[F(39,100)],(2,3):[F(55,100)]}
    found={}
    exact_cands=[]
    for a in (1,2):
        for d in range(16):
            seg=windows(a,d)
            if not seg: continue
            assert (a,d) in witnesses
            for q in witnesses[a,d]:
                assert L[0]<=value(2-a,d,q)<=L[1] and H[0]<=value(2,d,q)<=H[1]
            for k in range(-2,3):
                q=F(k,3)
                if L[0]<=value(2-a,d,q)<=L[1] and H[0]<=value(2,d,q)<=H[1]:
                    exact_cands.append((a,d,str(q)))
            true=finite_roots(a,d)
            assert len(true)==len(seg)
            for (lo,hi),(plo,phi) in zip(true,seg):
                assert plo-1e-11<=lo<hi<=phi+1e-11
            found[f"{a},{d}"]=dict(rational_outer=seg,quadratic_boundaries=true,
                                   simultaneous_witnesses=[str(v) for v in witnesses[a,d]])
    assert len(found)==5 and exact_cands==[(1,2,"2/3")]
    q0=(1+math.sqrt(7/45))/2
    eps0=(math.sqrt(7/5)-1)/2
    assert abs(value(1,1,q0)-float(L[0]))<1e-14
    assert H[0]<value(2,1,q0)<H[1] and abs(3*q0-2-eps0)<1e-15
    assert F(45,100)<=F(9,10)/(F(9,10)+F(11,10))<=F(55,100)
    for key in ("1,0","2,2","2,3"):
        for lo,hi in found[key]["rational_outer"]:
            assert all(lo*3-eps0>k or hi*3+eps0<k for k in range(-4,5))
    return dict(row_interval_exact=[str(Umin),str(Umax)],five_classes=found,
                integer_charge_results=exact_cands,nearest_new_charge=q0,
                sharp_error_threshold=eps0,CMS_reported=[.489,.005,.019,.020],
                CMS_two_templates=[0.5,0.5],full_new_group_tag_all_rows=1.0)


if __name__=="__main__":
    data=run()
    target=Path(__file__).with_name("independent_results.json")
    if target.exists():
        expected=json.loads(target.read_text(encoding="utf-8"))
        assert json.loads(json.dumps(data))==expected
    print(json.dumps(data,indent=2))

