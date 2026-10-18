#!/usr/bin/env python3
"""1311 independent check. Reuses only the frozen independent 1310 Dirac kernel.
Default/--check are read-only; --write-new creates the results once.
General quantifiers belong to proof.md, not these finite diagnostics.
"""
import argparse
import hashlib
import importlib.util
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
RESULT = HERE / "independent_results.json"
KERNEL_PATH = STAGE / "1310" / "independent_check.py"
KERNEL_SHA = "8c9444e9ea5d4ee38152ee2d6ce4067e4d0aed6bb63cb81a84379c56e3a58a04"

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def load_kernel():
    assert sha(KERNEL_PATH) == KERNEL_SHA, "The independently accepted kernel changed"
    spec = importlib.util.spec_from_file_location("independent1310", KERNEL_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def kernel_measure(kernel, order):
    mu = 2.0**-12
    m, n, energy = mu, 2*mu, .1
    r, s = math.sqrt(energy**2-m**2), math.sqrt(energy**2-n**2)
    nodes, weights = np.polynomial.legendre.leggauss(order)
    cosines, weights = nodes/2, weights/2
    zs, masses = [], []
    for c, wt in zip(cosines, weights):
        point = kernel.sample(m, n, energy, float(c))
        # Independently integrate dOmega instead of the author's z kernel.
        dmu = (1/9)*abs(point["jj"])**2*(s/r)/(128*math.pi*energy**2)
        zs.append(point["z"])
        masses.append(float(wt*dmu))
    zs, masses = np.asarray(zs), np.asarray(masses)
    assert np.all(zs > 0) and np.all(masses > 0)
    z0, z1 = r*r+s*s-r*s, r*r+s*s+r*s
    return zs, masses, z0, z1

def rate(zs, weights, x, w):
    return float(np.sum(weights/(zs*x+w)**2)/4)

def derivatives(zs, weights, x, w):
    d = zs*x+w
    fx = -float(np.sum(weights*zs/d**3))/2
    fw = -float(np.sum(weights/d**3))/2
    fxx = 1.5*float(np.sum(weights*zs*zs/d**4))
    fxw = 1.5*float(np.sum(weights*zs/d**4))
    fww = 1.5*float(np.sum(weights/d**4))
    slope = -fx/fw
    curvature = -(fxx+2*fxw*slope+fww*slope*slope)/fw
    return slope, curvature, fx, fw

def decreasing_root(fn, target, hi):
    assert fn(0) >= target
    assert fn(hi) <= target*(1+1e-13)
    lo = 0.
    for _ in range(85):
        mid = (lo+hi)/2
        if fn(mid) > target:
            lo = mid
        else:
            hi = mid
    return (lo+hi)/2

def run():
    kernel = load_kernel()
    z64, k64, z0, z1 = kernel_measure(kernel, 64)
    zs, ks, _, _ = kernel_measure(kernel, 96)
    x0, w0 = 25., 1/64
    sig0 = rate(zs, ks, x0, w0)
    target = sig0/2
    sig64 = rate(z64, k64, x0, w0)
    relquad = abs(sig64-sig0)/sig0
    assert relquad < 2e-8
    ktot, c2 = float(sum(ks)), float(sum(ks/zs**2))
    xc, wc = math.sqrt(c2/target)/2, math.sqrt(ktot/target)/2
    assert sig0 > 0 and xc > x0 and wc > w0
    assert abs(rate(zs, ks, xc, 0)/target-1) < 5e-15
    assert abs(rate(zs, ks, 0, wc)/target-1) < 5e-15

    def boundary(x):
        assert 0 <= x <= xc
        return decreasing_root(lambda w: rate(zs, ks, x, w), target, wc)

    curves = []
    for fraction in (.1, .3, .5, .7, .9):
        x = fraction*xc
        w = boundary(x)
        slope, curvature, fx, fw = derivatives(zs, ks, x, w)
        h = xc*1e-4
        wm, wp = boundary(x-h), boundary(x+h)
        slope_num = (wp-wm)/(2*h)
        curvature_num = (wp-2*w+wm)/(h*h)
        assert -z1 < slope < -z0 and curvature > 0
        assert abs(rate(zs, ks, x, w)/target-1) < 5e-15
        assert abs(slope_num-slope) < 2e-7*abs(slope)
        assert abs(curvature_num-curvature) < 2e-5*curvature
        curves.append({"x":x,"w":w,"slope":slope,"positive_curvature":curvature,
                       "finite_difference_curvature":curvature_num})
    # A convex boundary's subgraph must not be mislabeled convex.
    xa, xb = .25*xc, .75*xc
    wa, wb = boundary(xa), boundary(xb)
    xmid, wmid = (xa+xb)/2, (wa+wb)/2
    assert wmid > boundary(xmid)
    assert rate(zs, ks, xmid, wmid) < target

    # Exact running-domain reduction, checked independently against a time grid.
    betas = [Fraction(-19,6),Fraction(-11,6),Fraction(-1,2)]
    assert [(-43+8*n)/6 for n in (3,4,5)] == [float(b) for b in betas]
    budget_checks = 0
    domain_equal_points = 0
    positive_contractions = 0
    for factor in (.1,.5,1.,3.,10.):
        sstar = target*factor
        for length in (.1,1.,15.):
            for gmax in (.15,.3,.6):
                xf = gmax**-2
                for x in (xf*.8,xf,xf*1.01,xf*2,xf*8):
                    for w in (0.,.002,.01,.0225,.05):
                        flags=[]
                        for beta in betas+[Fraction(5,6),Fraction(19,6)]:
                            bb = float(beta)
                            bound = xf+max(bb,0)*length/(8*math.pi**2)
                            implicit = x >= bound and w>0 and rate(zs,ks,x,w)>=sstar
                            flow = min(x-bb*t/(8*math.pi**2)
                                       for t in np.linspace(0,length,17)) >= xf
                            direct = flow and w>0 and rate(zs,ks,x,w)>=sstar
                            assert implicit == direct
                            flags.append(implicit)
                            budget_checks+=1
                        assert flags[0] == flags[1] == flags[2]
                        assert not flags[3] or flags[0]
                        assert not flags[4] or flags[3]
                        domain_equal_points += 1
                        positive_contractions += bool(flags[0] and not flags[4])
    assert positive_contractions > 0

    # Closed-box endpoints and no rescue by negative beta when the initial cap fails.
    xf = .3**-2
    wminus, wplus = .2**2/4, .3**2/4
    edge_target = rate(zs,ks,xf,wminus)
    for dx,dw in ((1e-3,0),(0,1e-5),(1.,.001)):
        assert rate(zs,ks,xf+dx,wminus+dw) < edge_target
    empty_target = rate(zs,ks,.8*xf,wminus)
    assert rate(zs,ks,xf,wminus) < empty_target
    assert all(xf+max(float(b),0)/(8*math.pi**2)==xf for b in betas)
    roots = []
    for mult in (.5,1.,1.2):
        st = edge_target*mult
        if rate(zs,ks,0,wminus) > st:
            upper = math.sqrt(c2/st)/2
            xm = decreasing_root(lambda x:rate(zs,ks,x,wminus),st,upper)
            B = 8*math.pi**2*(xm-xf)
            roots.append({"target_factor":mult,"X_minus":xm,"B_minus":B,
                          "nonempty":bool(rate(zs,ks,xf,wminus)>=st)})
    assert roots[0]["B_minus"]>0 and roots[2]["B_minus"]<0
    assert abs(roots[1]["B_minus"])<1e-10

    # Exact integer transport: diagnostic of the inherited arithmetic and endpoints.
    integer_counts=[]
    for B in (Fraction(0),Fraction(5,6),Fraction(2),Fraction(10)):
        kmax = (43+6*B)//24
        formula={(4*k-1,n) for k in range(1,kmax+1)
                 for n in range(3,int((43+6*B)//(8*k))+1)}
        direct={(4*k-1,n) for k in range(1,15) for n in range(3,100)
                if Fraction(-43+8*k*n,6)<=B}
        assert formula==direct
        integer_counts.append({"B":str(B),"count":len(formula),"pairs":[list(pair) for pair in sorted(formula)]})

    # Uniform rectangle: exact rational denominator bound, not just point sampling.
    eta = Fraction(19,20)**4
    assert eta > Fraction(1,2)
    assert Fraction(26,25)**2 < Fraction(20,19)**2
    probability_lower = 1-math.exp(-2*math.log(2)*float(eta))
    assert probability_lower > .5
    rectangle_min = min(rate(zs,ks,g**-2,v*v/4)/sig0
                        for g in (.19,.20,.21) for v in (.24,.25,.26))
    assert rectangle_min >= float(eta)
    assert all(.21**-2-float(b)/(8*math.pi**2) >= .3**-2 for b in betas)

    # Same calibrated flavor columns for all three candidate branches.
    mu=2**-12
    omega=np.exp(2j*math.pi/3)
    U3=np.array([[omega**(a*b) for b in range(3)] for a in range(3)])/math.sqrt(3)
    mixing=[]
    for n in (3,4,5):
        U=np.eye(n,dtype=complex);U[:3,:3]=U3
        for j in range(3,n):
            theta=2e-5/(j+1)
            R=np.eye(n);c,s=math.cos(theta),math.sin(theta)
            R[2,2]=R[j,j]=c;R[2,j]=s;R[j,2]=-s
            U=U@R
        fixed=np.vstack([U3[:,:2],np.zeros((n-3,2))])
        assert np.max(abs(U[:,:2]-fixed))<1e-14
        mass=mu*np.arange(1,n+1)
        H=U@np.diag(mass**2)@U.conj().T
        visited={0}
        while True:
            expanded=visited|{j for i in visited for j in range(n) if abs(H[i,j])>1e-24}
            if expanded==visited:break
            visited=expanded
        assert len(visited)==n
        mcc=U[1,0]*U[0,1].conjugate()
        assert abs(abs(mcc)-1/3)<1e-14
        mixing.append({"generations":n,"abs_mixing":float(abs(mcc)),
                       "first_two_columns_fixed":True,"flavor_graph_connected":True})

    # Common effective range is input; finite numerical values only check its compatibility.
    vmax,vmin,gy,gmax=.3,.2,.1,.3
    yukawa_bound=math.sqrt(2)*.003/vmin
    massW=gmax*vmax/2
    massZ=math.sqrt(gmax*gmax+gy*gy)*vmax/2
    massH=math.sqrt(2*.04)*vmax
    assert yukawa_bound<.022 and max(massW,massZ,massH)<.1
    mh_pair=[math.sqrt(2*l)*.25 for l in (.02,.04)]
    assert mh_pair[0]!=mh_pair[1]
    luminosity=2*math.log(2)/sig0
    assert math.isfinite(luminosity)
    sourcepaths={
        "admission":STAGE/"_admission"/"round1311_admission.md",
        "proof":HERE/"proof.md","author_check":HERE/"check.py",
        "author_results":HERE/"results.json","note":STAGE/"research_note_1311.md",
        "next":HERE/"NEXT.md","1310_proof":STAGE/"1310"/"proof.md",
        "1310_acceptance":STAGE/"1310"/"acceptance.json",
        "1310_independent_kernel":KERNEL_PATH,
        "1309_proof":STAGE/"1309"/"proof.md","self":Path(__file__).resolve()}
    bindings={name:{"path":str(p.relative_to(STAGE.parent)),"sha256":sha(p)}
              for name,p in sourcepaths.items()}
    return {
        "status":"PASS_WITH_SCOPE","bindings":bindings,
        "method":"Independent frozen 1310 Ward-null dOmega kernel; implicit Hessian and root checks; exact rational margins",
        "kernel":{"m":mu,"n":2*mu,"E":.1,"z0":z0,"z1":z1,
                  "K_integral":ktot,"C2":c2,"sigma0":sig0,"relative_quadrature_change":relquad,
                  "S_star":target,"finite_exposure":luminosity},
        "boundary":{"x_c":xc,"w_c":wc,"derivative_checks":curves,
                    "subgraph_midpoint_infeasible":True},
        "domain":{"finite_flow_checks":budget_checks,"three_branch_equal_points":domain_equal_points,
                  "strict_positive_b_contractions":positive_contractions,"closed_box":roots,
                  "integer_transport":integer_counts},
        "common_slice":{"rectangle_exact_lower":str(eta),"rectangle_sample_min_ratio":rectangle_min,
                        "probability_uniform_lower":probability_lower,"mixing":mixing,
                        "yukawa_upper":yukawa_bound,"mass_W_upper":massW,
                        "mass_Z_upper":massZ,"mass_H_upper":massH,
                        "different_Higgs_masses_at_fixed_v":mh_pair},
        "scope":["Analytic proof carries all-budget equality and strict curvature; finite tests do not.",
                 "Equality concerns this calibrated leading task and common effective box, not all predictions.",
                 "Fixed masses require changing Yukawas; each point retains its own complete action and source.",
                 "No computed actual-device or higher-order error; no finite-resource feasibility claim beyond the declared budget."]
    }

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    parser.add_argument("--write-new",action="store_true")
    args=parser.parse_args()
    result=run()
    if args.write_new:
        with RESULT.open("x",encoding="utf-8",newline="\n") as f:
            json.dump(result,f,ensure_ascii=False,indent=2,allow_nan=False);f.write("\n")
    elif RESULT.exists():
        stored=json.loads(RESULT.read_text(encoding="utf-8"))
        assert stored==result,"Frozen results or source versions changed"
    print(json.dumps({"status":result["status"],"bindings":len(result["bindings"]),
                      "flow_checks":result["domain"]["finite_flow_checks"],
                      "sigma0":result["kernel"]["sigma0"],
                      "write_new":args.write_new},ensure_ascii=False))

if __name__=="__main__":
    main()
