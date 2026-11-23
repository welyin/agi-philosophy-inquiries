"""Finite KK threshold and a held-out strong coupling in inverse-alpha units."""

import argparse
import itertools
import json
import math
from pathlib import Path


BASE = Path(__file__).resolve().parent
P = (5, -12, 7)
SM = (41/10, -19/6, -7)
MSSM = (33/5, 1, -3)
EVEN = (6/5, -2, -6)
ODD = (-46/5, -6, -2)
XY = (-10, -6, -4)
D = (2/5, 0, 1)
E = (3/5, 1, 0)
CASIMIR = (0, 2, 3)
TAU = 2*math.pi


def read(name):
    return json.loads((BASE/name).read_text(encoding="utf-8-sig"))


def project(v):
    return sum(w*x for w, x in zip(P, v))


def robin_shift(level, strength):
    """Return log(k_R/k_N); x=kL-x0 solves (x0+x) tan x=strength."""
    if strength == 0:
        return 0.
    x0 = (level-0.5)*math.pi
    lo, hi = 0., math.pi/2
    for _ in range(70):
        mid = (lo+hi)/2
        if (x0+mid)*math.tan(mid) <= strength:
            lo = mid
        else:
            hi = mid
    return math.log1p((lo+hi)/(2*x0))


def k_theory(n, strength, ms, mz, ratio):
    log_wall = sum(math.log(2*j/(2*j-1)) for j in range(1,n+1))
    shifts = sum(robin_shift(j,strength) for j in range(1,n+1))
    return ((9.5*math.log(ms/mz) + 12*log_wall + 9*math.log(ratio) + 6*shifts)/TAU
            + project(CASIMIR)/(12*math.pi))


def norm_ell_over_MD(rows, angle, ratio):
    d=[r[0] for r in rows[:3]]
    e=[r[0] for r in rows[3:]]
    poles=[0.,max(d[0],e[0])*min(d[1],e[1]),max(d[1],e[1])*min(d[2],e[2]),(1.2*max(d[2],e[2]))**2]
    def residues(s):
        return [math.prod(t-y*y for y in s)/math.prod(t-u for u in poles if u!=t) for t in poles]
    wd,we=residues(d),residues(e)
    cs=math.cos(math.radians(angle))
    squares=[]
    for j in (1,2,3):
        w=ratio**2*we[j]/wd[j]
        b=(6+4*w)*cs
        x=2*(1-w)/(b+math.sqrt(b*b-4*(9-4*w)*(1-w)))
        squares.append(wd[j]*x*x/(1+4*x*cs+4*x*x))
    ell=math.sqrt(sum(squares))
    assert 1-ratio < 5*ell < 1+ratio
    return ell


def full_match(n, strength, ratio, ms, mz, observed, ell):
    """Exhibit positive tree coefficients at the retained order with an explicit split threshold."""
    k = k_theory(n,strength,ms,mz,ratio)
    local = [0., (k-project(observed))/12, 0.]
    assert min(local) >= 0
    even = sum(math.log(n/j) for j in range(1,n+1))
    odd = sum(math.log(2*n/(2*j-1)) for j in range(1,n+1))
    shift = sum(robin_shift(j,strength) for j in range(1,n+1))
    constant = [d + c/(12*math.pi)
                +(sm*math.log(ms/mz)-eb*math.log(ratio)
                  +be*even+bo*odd-xy*shift)/TAU
                for d,c,sm,eb,be,bo,xy in zip(local,CASIMIR,SM,E,EVEN,ODD,XY)]
    slopes = [m/TAU for m in MSSM]
    z = (observed[0]-observed[1]-constant[0]+constant[1])/(slopes[0]-slopes[1])
    cutoff = ms*math.exp(z)
    length = math.pi*n/cutoff
    md=math.sqrt(strength)/(5*ell*length)
    me=md*ratio
    full_pair=math.log(cutoff/md)/TAU
    common = observed[0]-constant[0]-slopes[0]*z-full_pair
    restored = [common+full_pair+c+b*z for c,b in zip(constant,slopes)]
    assert cutoff > max(md,me) and min(md,me)>ms and common>0
    assert all(abs(x-y)<1e-10 for x,y in zip(restored,observed))
    assert abs((5*ell*md*length)**2-strength)<1e-12
    return dict(cutoff_GeV=cutoff, Mc_GeV=cutoff/(2*n), MD_GeV=md, ME_GeV=me,
                splitting_norm_over_MD=ell, same_order_parameter_bL=(5*ell*md*length)**2,
                boundary_inverse_alpha=local, universal_inverse_alpha=common,
                recovered_inverse_alphas=restored)


def certificate():
    cfg=read("inputs.json")
    g=read(cfg["matched_couplings"])
    rows=read(cfg["matched_matter"])["scales"][cfg["mass_scale_label"]]
    widths=cfg["widths"]
    raw=[g[k] for k in ("g1_GUT","g2","g3")]
    central=[4*math.pi/x**2 for x,_ in raw]
    ranges=[(4*math.pi/(x+widths*e)**2,4*math.pi/(x-widths*e)**2) for x,e in raw]
    krange=[min(map(project,itertools.product(*ranges))),max(map(project,itertools.product(*ranges)))]
    def ratio(sign):
        return math.prod((x+sign*widths*dx)/(y-sign*widths*dy)
                         for (x,dx),(y,dy) in zip(rows[:3],rows[3:]))
    rlo,r0,rhi=[ratio(s) for s in (-1,0,1)]
    mz=g["mu_ref_GeV"]
    ms=cfg["common_susy_threshold_GeV"]
    assert abs(project(SM)-9.5)<1e-14 and project(MSSM)==0
    assert (project(EVEN),project(ODD),project(XY),project(D),project(E),project(CASIMIR))==(-12,12,-6,9,-9,-3)
    results={}
    for n in cfg["level_pairs"]:
        kmin=k_theory(n,0,ms,mz,rlo)
        k0=k_theory(n,0,ms,mz,r0)
        deficit=kmin-krange[1]
        alpha3hi=7/(kmin-5*ranges[0][1]+12*ranges[1][0])
        cap=mz*math.exp((krange[1]-k_theory(n,0,mz,mz,rlo))*TAU/9.5)
        results[str(n)]=dict(K_prediction_min=kmin, K_prediction_central=k0,
                             alpha3_prediction_max=alpha3hi,
                             alpha3_prediction_central=7/(k0-5*central[0]+12*central[1]),
                             signed_residual_at_most=-deficit, positive_brane_budget_min=deficit/12,
                             allowed_common_MS_GeV_at_most=cap,
                             largest_adj_triplet_to_octet_ratio_if_XY_chiral_at_octet=math.exp(-TAU*deficit/24),
                             robin_additions={str(b):k_theory(n,b,ms,mz,r0)-k0 for b in cfg["boundary_mass_times_length"]})
    measured=sorted(x**-1 for x in ranges[2])
    assert results["1"]["alpha3_prediction_max"]<measured[0] and results["1"]["allowed_common_MS_GeV_at_most"] < 161
    point=cfg["same_spurion_example"]
    ell=norm_ell_over_MD(rows,point["phase_degrees"],r0)
    return dict(readout="P=(5,-12,7); A_i=inverse alpha_i with GUT hypercharge",
                R_lo_central_hi=[rlo,r0,rhi], observed_K=project(central),
                two_width_K=krange, observed_alpha3_two_width=measured, even_odd_pairs=results,
                matched_threshold_example=full_match(point["paired_levels"],point["bL"],r0,ms,mz,central,ell))


def equivalent(x,y):
    if isinstance(x,dict):
        return isinstance(y,dict) and x.keys()==y.keys() and all(equivalent(x[k],y[k]) for k in x)
    if isinstance(x,list):
        return isinstance(y,list) and len(x)==len(y) and all(equivalent(i,j) for i,j in zip(x,y))
    if isinstance(x,float):
        return isinstance(y,(float,int)) and math.isclose(x,y,rel_tol=1e-11,abs_tol=1e-12)
    return x==y


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--check", action="store_true")
    args=parser.parse_args()
    result=certificate()
    if args.write:
        with (BASE/"results.json").open("x",encoding="utf-8") as out:
            json.dump(result,out,indent=2,allow_nan=False)
            out.write("\n")
    elif not equivalent(result,read("results.json")):
        raise SystemExit("stored numerical certificate differs")
    print("1519 finite towers, held-out alpha_s, boundary and matter budgets verified.")
