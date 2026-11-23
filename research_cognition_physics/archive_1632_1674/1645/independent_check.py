"""Independent CP-even calibration of two finite quark-mass textures."""
import argparse
import cmath
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE / path).read_bytes())


def window(values, widths=2):
    spread = widths * sum(values[1:])
    return values[0] - spread, values[0] + spread


def left_diagonalizer(matrix):
    left, values, _ = np.linalg.svd(matrix)
    order = np.argsort(values)
    return left[:, order], values[order]


def mixing_kernel():
    def equal(x):
        return x * (4 - 11*x + x*x) / (4*(1-x)**2) - 1.5*x**3*math.log(x)/(1-x)**3
    def unequal(x, y):
        def coef(t):
            return .25 + 1.5/(1-t) - .75/(1-t)**2
        return x*y*((coef(x)*math.log(x)-coef(y)*math.log(y))/(x-y)-.75/((1-x)*(1-y)))
    weak = read("../1404/inputs.json")
    masses = [.0022, 1.27, weak["short_distance"]["mtMS"][0]]
    xs = [(m / weak["const_GeV"]["MW"])**2 for m in masses]
    return np.array([[equal(x) if x == y else unequal(x, y) for x in xs] for y in xs])


def phase_readout(V, kernel):
    lambdas = V[:, 0].conj()*V[:, 2]
    assert abs(sum(lambdas)) < 2e-12
    kk = V[1, 2]*V[0, 0].conjugate()*V[1, 1].conjugate()*V[0, 1]
    M = lambdas @ kernel @ lambdas
    top = lambdas[2]**2 * kernel[-1, -1]
    beta2 = cmath.phase(top/lambdas[1]**2)
    phase = cmath.phase(M/kk**2)
    return {"Kud_Kus_Kub_over_Kcb_Kcb": [abs(V[0, 0]), abs(V[0, 1]), abs(V[0,2]/V[1,2]), abs(V[1,2])],
            "sin_2beta": math.sin(beta2), "S_KK_LO": math.sin(phase),
            "phase_2beta_to_LO_KK_rad": cmath.phase(M/kk**2/(top/lambdas[1]**2))}


def rotate12(tan):
    c = 1/math.hypot(1, tan)
    return np.array([[c, -c*tan, 0], [c*tan, c, 0], [0, 0, 1]])


def fx_mass_pair(witness, yukawa, sign):
    u = witness["up_tangent"]
    d = math.sqrt(witness["down_yukawa"]/witness["strange_yukawa"])
    s = witness["Vcb"]*math.hypot(1, u)
    c = math.sqrt(1-s*s)
    U, D = rotate12(u), np.array([[1,0,0],[0,c,s],[0,-s,c]])@rotate12(d)
    up = U@np.diag([-yukawa["up"][0][0], yukawa["up"][1][0], yukawa["up"][2][0]])@U.T
    bare = D@np.diag([-witness["down_yukawa"],witness["strange_yukawa"],yukawa["down"][2][0]])@D.T
    assert abs(bare[0,0]) < 5e-17
    down = np.diag([-sign*1j,1,1])@bare
    return left_diagonalizer(up)[0].conj().T@left_diagonalizer(down)[0]


def fx_range(tree, belle, yukawa, kernel):
    # Monotone branch u<.1. g increases with u,d and decreases with w throughout these input boxes.
    u0,u1 = window(tree["ub_over_cb"])
    w0,w1 = window(belle["BelleII_Vcb_BGL"])
    yd,ys = window(yukawa["down"][0]), window(yukawa["down"][1])
    d0,d1 = math.sqrt(yd[0]/ys[1]),math.sqrt(yd[1]/ys[0])
    kmax = window(tree["K_product"])[1]/window(tree["K_form"])[0]
    def kaon(u,d,w):
        x,y = u*u,d*d
        return math.sqrt((x+y)/(1+x)/(1+y)-w*w*x/(1+y))
    # Positive lower derivatives with respect to u²,d², respectively.
    assert (1-d1*d1)/(1+u1*u1)**2 > w1*w1 and u1<1
    def cap(lo,hi,fun):
        assert fun(lo) < kmax < fun(hi)
        for _ in range(60):
            mid = (lo+hi)/2
            if fun(mid) <= kmax:
                lo = mid
            else:
                hi = mid
        return hi + 2e-12
    umax = cap(u0, u1, lambda u:kaon(u,d0,w1))
    dmax = cap(d0, d1, lambda d:kaon(u0,d,w1))
    cmin = math.sqrt(1-w1*w1*(1+umax*umax))
    cmax = math.sqrt(1-w0*w0*(1+u0*u0))
    # |Vtd Vts| <= (1-c²)d, |Vud| >= sqrt(1-kmax²-ubcap²).
    ubcap = w1*umax
    udmin = math.sqrt(1-kmax*kmax-ubcap*ubcap)
    tagcap = math.asin((1-cmin*cmin)*dmax/(udmin*kaon(u0,d0,w1)))
    cdmin = cmin*d0/math.sqrt((1+umax**2)*(1+dmax**2))
    rcap = umax/cdmin
    coeff = [rcap/(1-rcap), 1/(1-rcap), 1.]
    gimcap = sum(coeff[i]*coeff[j]*kernel[i,j]/kernel[2,2] for i in range(3) for j in range(3) if (i,j)!=(2,2))
    phasepad = 2*tagcap + math.asin(gimcap) + 1e-10
    lower = 2*math.atan(u0/(cmax*dmax))-phasepad
    upper = 2*math.atan(umax/(cmin*d0))+phasepad
    udmax = window(tree["ud_beta"])[1]
    nuclear_floor = math.sqrt(1-udmax**2-kmax**2)
    assert 0 < lower < upper < math.pi/2 and ubcap < nuclear_floor
    return {"upper_tangent_u":umax,"upper_tangent_d":dmax,"phasepad_rad":phasepad,
            "abs_S_KK_LO_closed_outer":[math.sin(lower),math.sin(upper)],
            "Kub_max_before_beta":ubcap,"Kub_min_with_nuclear_beta":nuclear_floor}


def free_texture(delta, kernel, yukawa):
    # Prescribed CP-even semileptonic and first-row quantities only.
    us, cb, ratio = .224999, .04420, .0989
    s13 = cb*ratio
    c13 = math.sqrt(1-s13*s13)
    s12, s23 = us/c13, cb/c13
    c12, c23 = math.sqrt(1-s12*s12), math.sqrt(1-s23*s23)
    z = cmath.exp(1j*delta)
    K = np.array([[c12*c13, s12*c13, s13/z],
                  [-s12*c23-c12*s23*s13*z, c12*c23-s12*s23*s13*z, s23*c13],
                  [s12*s23-c12*c23*s13*z,-c12*s23-s12*c23*s13*z, c23*c13]])
    y = np.array([row[0] for row in yukawa["down"]])
    h = (K*y**2)@K.conj().T
    rephase = np.diag([1,cmath.exp(1j*cmath.phase(h[0,1])),cmath.exp(1j*cmath.phase(h[0,2]))])
    h = rephase@h@rephase.conj().T
    a = math.sqrt(h[0,0].real)
    c,e = h[0,1:].real/a
    f = math.sqrt(h[2,2].real-e*e)
    coefficient = (h[1,2]-c*e)/f
    b = -math.sqrt(h[1,1].real-c*c-abs(coefficient)**2)
    eta = -2*math.pi/3
    d1 = coefficient.imag/math.sin(eta)
    d0 = coefficient.real-d1*math.cos(eta)
    Y = np.array([[0,a,0],[b,c,d0+d1*cmath.exp(1j*eta)],[0,e,f]])
    L, measured = left_diagonalizer(Y)
    np.testing.assert_allclose(measured,y,rtol=2e-11)
    np.testing.assert_allclose(Y@Y.conj().T,h,atol=1e-16)
    assert np.linalg.det(Y).real>0 and min(abs(d0),abs(d1))>0
    row = phase_readout(L,kernel)
    np.testing.assert_allclose(row["Kud_Kus_Kub_over_Kcb_Kcb"],
                               [math.sqrt(1-us**2-s13**2),us,ratio,cb],atol=3e-13)
    return {"delta_deg":math.degrees(delta),"Z_real_paths":[d0,d1],**row}


def offset_pair(tree, belle, yukawa, kernel):
    point = read("inputs.json")["calibrated_competitor"]
    u, w, us, f = (point[k] for k in ("up_tangent","cb","us","fplus"))
    d = math.sqrt(point["yd"]/point["ys"])
    s = w*math.hypot(1,u)
    c = math.sqrt(1-s*s)
    cosphi = (d*d+u*u*c*c-us*us*(1+u*u)*(1+d*d))/(2*u*d*c)
    assert abs(cosphi)<1
    phi = math.acos(cosphi)
    ou = rotate12(u)
    od = np.array([[1,0,0],[0,c,s],[0,-s,c]])@rotate12(d)
    up = ou@np.diag([(-1 if i==0 else 1)*v[0] for i,v in enumerate(yukawa["up"])])@ou.T
    bare = od@np.diag([-point["yd"],point["ys"],yukawa["down"][-1][0]])@od.T
    U = left_diagonalizer(up)[0]
    signed = []
    incl = read("../1640/inputs.json")["inclusive_holdout"]["Vcb"]
    for branch in (-1,1):
        mass = np.diag([cmath.exp(1j*branch*phi),1,1])@bare
        V = U.conj().T@left_diagonalizer(mass)[0]
        row = phase_readout(V,kernel)
        q, k, uread, vread = row["Kud_Kus_Kub_over_Kcb_Kcb"]
        np.testing.assert_allclose([k,uread,vread],[us,u,w],atol=1e-11)
        for x,v in [(q,tree["ud_beta"]),(k*f,tree["K_product"]),(f,tree["K_form"]),
                    (uread,tree["ub_over_cb"]),(vread,belle["BelleII_Vcb_BGL"]),(vread,incl)]:
            assert window(v)[0] < x < window(v)[1]
        signed.append(row)
    assert abs(signed[0]["S_KK_LO"]+signed[1]["S_KK_LO"])<1e-12
    return {"abs_phase_deg":math.degrees(phi),"conjugates":signed}


def solve():
    kernel = mixing_kernel()
    yukawa = read("../1386/inputs.json")["yukawa"]
    tree = read("../1381/inputs.json")
    belle = read("../1412/inputs.json")
    p = read("../1502/inputs.json")["witness"]
    fx = [phase_readout(fx_mass_pair(p,yukawa,sign),kernel) for sign in (1,-1)]
    assert abs(fx[0]["S_KK_LO"]+fx[1]["S_KK_LO"])<1e-12
    interval = fx_range(tree,belle,yukawa,kernel)
    assert interval["abs_S_KK_LO_closed_outer"][0] <= abs(fx[0]["S_KK_LO"]) <= interval["abs_S_KK_LO_closed_outer"][1]
    feasible = [free_texture(math.radians(x),kernel,yukawa) for x in [.001,60.,-60.,90.]]
    ud,us,rho,cb = feasible[0]["Kud_Kus_Kub_over_Kcb_Kcb"]
    for value, limits in [(ud,window(tree["ud_beta"])),(us*.9640001,window(tree["K_product"])),
                          (.9640001,window(tree["K_form"])),(rho,window(tree["ub_over_cb"])),
                          (cb,window(belle["BelleII_Vcb_BGL"]))]:
        assert limits[0] < value < limits[1]
    assert abs(feasible[0]["S_KK_LO"])<.0001 and feasible[1]["S_KK_LO"]>.6
    assert abs(feasible[1]["S_KK_LO"]+feasible[2]["S_KK_LO"])<1e-12
    obj = {"1502_physical_mass_pair":fx,"1502_domain":interval,"1489_fixed_eta_texture":feasible,
           "offset_cosine_competitor":offset_pair(tree,belle,yukawa,kernel)}
    def stabilize(v):
        if isinstance(v,dict): return {k:stabilize(x) for k,x in v.items()}
        if isinstance(v,list): return [stabilize(x) for x in v]
        return float(f"{v:.11g}")
    return stabilize(obj)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args = parser.parse_args()
    result = solve()
    path = HERE/"independent_results.json"
    if args.check:
        assert json.loads(path.read_bytes())==result
        print("1645 independent PASS: physical conjugates, first-row veto, free two-path counterexample")
    else:
        with path.open("x", encoding="utf-8") as out:
            out.write(json.dumps(result,indent=2)+"\n")

