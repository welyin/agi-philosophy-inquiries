"""Profile the additional real Majorana arm on newly shot charged matter."""
import argparse
import json
import math
import runpy
from pathlib import Path
from types import SimpleNamespace
from fractions import Fraction

import numpy as np
from numpy.polynomial import Polynomial as Poly


HERE = Path(__file__).resolve().parent
ARCH = HERE.parent
DP = .0851
KAM = (7.54e-5, .18e-5, .19e-5)
SK = (2.50e-3, .31e-3, .13e-3)
CB = .04216


def sturm_real_count(coefficients):
    """Count distinct real roots of the stored binary coefficients exactly."""
    first=[Fraction.from_float(float(a)) for a in coefficients]
    polys=[first,[i*a for i,a in enumerate(first) if i]]
    while len(polys[-1])>1:
        rem=polys[-2].copy()
        divisor=polys[-1]
        while len(rem)>=len(divisor):
            ratio=rem[-1]/divisor[-1]
            offset=len(rem)-len(divisor)
            for j,c in enumerate(divisor):
                rem[offset+j]-=ratio*c
            while rem and rem[-1]==0:
                rem.pop()
        assert rem, "A repeated root needs separate isolation."
        polys.append([-a for a in rem])
    def changes(side):
        signs=[(p[-1]>0) == (side>0 or len(p)%2==1) for p in polys]
        return sum(a!=b for a,b in zip(signs,signs[1:]))
    return changes(-1)-changes(1)


def actual(u, step=.25):
    family = runpy.run_path(str(ARCH/"1755"/"joint.py"))
    captured = {}
    def load(path):
        module = runpy.run_path(str(path))
        if Path(path).parent.name == "1748":
            original = module["daughter"]
            def measured(bg, F):
                nu = bg["nu"]
                old = nu["observables"]
                evaluations = []
                def detector(C, b, rho):
                    evaluations.append((C.copy(), b.copy(), rho))
                    return old(C,b,rho)
                nu["observables"] = detector
                try:
                    scale,row = original(bg,F)
                finally:
                    nu["observables"] = old
                C,b,rho = evaluations[-2]  # daughter followed by its mother-scale control
                np.testing.assert_allclose(old(C,b,rho)["angle_readings"],
                                           row["neutrinos"]["angle_readings"], atol=2e-13)
                captured.update(C=C,b=b,rho=rho)
                return scale,row
            module["daughter"] = measured
        return module
    family["point"].__globals__["runpy"] = SimpleNamespace(run_path=load)
    row = family["shoot"](ARCH/"1755",1/19.7,u,CB,step)
    return row,captured


def model(C,b):
    _,s,vh = np.linalg.svd(C)
    assert s[2]/s[1] < 1e-10 and abs(vh[-1].conj() @ b) < 2e-9
    Q = vh.conj().T[:,:2]
    Qe = Q[0]
    C = Q.T@C@Q
    b = Q.T@b
    gain = np.linalg.norm(C,2)/np.vdot(b,b).real
    B = np.outer(b,b)*gain
    n = np.linalg.norm(C,2)
    O, V = C/n, B/n
    H0, H1, H2 = O.conj().T@O, O.conj().T@V+V.conj().T@O, V.conj().T@V
    T = Poly([np.trace(H).real for H in (H0,H1,H2)])
    E = Poly([(Qe@H@Qe.conj()).real for H in (H0,H1,H2)])
    m = np.vdot(Qe,Qe).real
    d = Poly([np.linalg.det(O),O[0,0]*V[1,1]+O[1,1]*V[0,0]-2*O[1,0]*V[0,1]])
    P = Poly((d*Poly(d.coef.conj())).coef.real)
    gap2 = T*T-4*P
    def forward(x):
        _,spec,rh = np.linalg.svd(O+x*V)
        e3 = abs(Qe@rh.conj().T[:,0])**2
        return dict(rho=x*gain, reactor=4*e3*(1-e3), R=spec[0]/spec[1],
                    SK_dm32_using_KamLAND=KAM[0]*((spec[0]/spec[1])**2-1),
                    JUNO_dm31_using_JUNO_solar=7.388e-5*(spec[0]/spec[1])**2)
    def roots(poly, key, target):
        roots = poly.roots()
        assert min(abs(np.subtract.outer(roots,roots)+10*np.eye(len(roots))).flat) > 1e-7
        assert sum(abs(v.imag)<1e-7 for v in roots)==sturm_real_count(poly.coef)
        out = []
        for value in roots:
            if abs(value.imag)>1e-7:
                continue
            x = value.real
            obs = forward(x)
            if abs(obs[key]-target)>1e-7:
                continue
            obs["dA_dlog_abs_rho"] = (forward(x*math.exp(1e-5))["reactor"]-
                                        forward(x*math.exp(-1e-5))["reactor"])/2e-5
            out.append(obs)
        return sorted(out,key=lambda v:v["rho"])
    def reactor(a):
        out = []
        for p in ((1-math.sqrt(1-a))/2, (1+math.sqrt(1-a))/2):
            if p > m:   # No heavy column can have this much electron weight.
                continue
            eq = (2*E-m*T)**2-(2*p-m)**2*gap2
            out.extend(roots(eq,"reactor",a))
        return out
    def clock(r):
        return roots(T*T-(r+1/r)**2*P,"R",r)
    return gain,forward,clock,reactor


def at(u, step=.25):
    row,parts=actual(u,step)
    gain, calc, clock, reactor=model(parts["C"],parts["b"])
    fixed=calc(parts["rho"]/gain)
    np.testing.assert_allclose(fixed["reactor"],row["angle_readings"][0],atol=2e-12)
    h=math.sqrt(1+SK[0]/KAM[0])
    ratio_lo=math.sqrt(1+(SK[0]-SK[1])/(KAM[0]+KAM[2]))
    ratio_hi=math.sqrt(1+(SK[0]+SK[2])/(KAM[0]-KAM[1]))
    observations=reactor(DP)
    for obs in observations:
        obs["passes_marginal_free_clocks"]=ratio_lo<=obs["R"]<=ratio_hi
    intersections={}
    for widths in (1,2):
        a_lo, a_hi=DP-widths*.0024,DP+widths*.0024
        endpoints=sorted(p["rho"] for group in (reactor(a_lo),reactor(a_hi),clock(ratio_lo),clock(ratio_hi))
                         for p in group)
        def allowed(rho):
            reading=calc(rho/gain)
            return a_lo-1e-10<=reading["reactor"]<=a_hi+1e-10 and ratio_lo-1e-9<=reading["R"]<=ratio_hi+1e-9
        intervals=[[left,right] for left,right in zip(endpoints,endpoints[1:]) if allowed((left+right)/2)]
        isolated=[val for val in endpoints if allowed(val) and not any(l<=val<=r for l,r in intervals)]
        # At both infinities the mass ratio diverges because V has rank one.
        intersections[str(widths)]=dict(intervals=intervals,isolated_endpoints=isolated)
    return dict(Vub=u,Vus=row["low_moduli"]["Vus"],Vcb=row["low_moduli"]["Vcb"],
                d_over_s=row["down_over_s"], frozen=fixed, reactor_central_real_roots=observations,
                KamLAND_SK_central_real_roots=clock(h),
                JUNO_central_real_roots=clock(math.sqrt(.002509/.00007388)),
                all_real_rho_joint_boxes=intersections,
                independent_frequency_mass_ratio_window=[ratio_lo,ratio_hi])


def inverse(rows):
    """Reconstruct the root jointly selected by the three detector scalars."""
    target=math.sqrt(1+SK[0]/KAM[0])
    def residual(row):
        positive, = [r for r in row["reactor_central_real_roots"] if r["rho"]>0]
        return positive["R"]-target
    l,r=rows[:2]
    assert residual(l)<0<residual(r)
    for _ in range(7):
        if abs(residual(r))<3e-9: break
        u=r["Vub"]-residual(r)*(r["Vub"]-l["Vub"])/(residual(r)-residual(l))
        l,r=r,at(u)
    assert abs(residual(r))<3e-9
    p, = [t for t in r["reactor_central_real_roots"] if t["rho"]>0]
    belle=[.00393,.00009,.00013,.00019]
    return dict(Vub=r["Vub"],Vus=r["Vus"],rho=p["rho"],reactor=p["reactor"],
                dm32_eV2_using_KamLAND=p["SK_dm32_using_KamLAND"],
                exclusive_Bpi_marginal_offset=(r["Vub"]-belle[0])/math.hypot(*belle[1:]),
                inclusive_BXu_marginal_offset=(r["Vub"]-.00410)/math.hypot(.00009,.00022,.00015))


def clean(o):
    if isinstance(o,dict): return {k:clean(v) for k,v in o.items()}
    if isinstance(o,(tuple,list)): return [clean(v) for v in o]
    if isinstance(o,(bool,np.bool_)): return bool(o)
    if isinstance(o,(float,np.floating)): return float(f"{o:.10g}")
    return o


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    args=p.parse_args()
    rows=[at(u) for u in (.0036,.00385,.0041,.00435)]
    fine = at(.0041,.125)
    coarse=rows[2]
    error=max(abs(p["rho"]-q["rho"]) for p,q in zip(fine["reactor_central_real_roots"],
                                                       coarse["reactor_central_real_roots"],strict=True))
    assert error<2e-8
    result=clean(dict(reactor=DP,independent_KamLAND=KAM,independent_SK=SK,
                      rows=rows,combined_central_solution=inverse(rows),step_halving_max_delta_rho=error))
    dest=HERE/"independent_results.json"
    if args.check:
        assert result==json.loads(dest.read_bytes())
        print("1757 archived port, open Majorana gain, two independent frequency calibrations verified")
    else:
        with dest.open("w",encoding="utf8") as f:
            json.dump(result,f,indent=2,allow_nan=False)
            f.write("\n")
