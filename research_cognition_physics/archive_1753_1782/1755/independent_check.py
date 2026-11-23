"""Independent two-angle inversion and a bottom-mass-free branch test (round 1755)."""
import argparse
import cmath
import json
import math
from pathlib import Path


HERE = Path(__file__).resolve().parent
ARCH = HERE.parent


def read(path):
    return json.loads(path.read_bytes())


def sine_range(lo, hi):
    assert 0 <= lo <= hi <= math.pi
    return min(math.sin(lo), math.sin(hi)), (1 if lo <= math.pi/2 <= hi else
                                              max(math.sin(lo), math.sin(hi)))


def angles(sign, n, data):
    centre, *unc = data["S_d"]
    direct, *unc_direct = data["C_d"]
    s_unc, c_unc = math.hypot(*unc), math.hypot(*unc_direct)
    cmax = abs(direct) + n*c_unc
    cmin = max(0., abs(direct)-n*c_unc)
    a = math.asin((centre-n*s_unc)/math.sqrt(1-cmin*cmin))
    b = math.asin((centre+n*s_unc)/math.sqrt(1-cmax*cmax))
    beta = (a/2, b/2) if sign == "acute_2beta" else ((math.pi-b)/2, (math.pi-a)/2)
    gammas = data["time_integrated_gamma_deg"]
    gamma = [gammas["central"]]*2 if n == 0 else gammas["68" if n == 1 else "95"]
    gamma = list(map(math.radians, gamma))
    sl, sh = sine_range(*beta)
    tl, th = sine_range(beta[0]+gamma[0], beta[1]+gamma[1])
    return beta, gamma, (sl/th, sh/tl)


def mass_root(D, B, ub_high):
    def residual(x):
        c = math.sqrt(1-x*x-ub_high*ub_high)
        return x*x - c*D - ub_high**2*B/math.sqrt((1+c)/2)
    lo, hi = 0., .75
    assert residual(lo) < 0 < residual(hi)
    for _ in range(60):
        mid = (lo+hi)/2
        if residual(mid) < 0:
            lo = mid
        else:
            hi = mid
    return (lo+hi)/2


def central(branch, data, cb):
    beta, gamma, _ = angles(branch, 0, data)
    beta, gamma = beta[0], gamma[0]
    apex = math.sin(beta)/math.sin(beta+gamma)*cmath.exp(1j*gamma)
    ub = .005
    # Concrete interior transport of the declared effective box.
    for _ in range(60):
        us = mass_root(1/19.7, .973*52., ub*1.025)
        ud = math.sqrt(1-us*us-ub*ub)
        s23 = cb/math.sqrt(1-ub*ub)
        w = us/ud * math.sqrt(1-s23*s23)*s23*apex/(1-s23*s23*apex)
        ub = abs(w)
    c13 = math.sqrt(1-ub*ub)
    s12, c12 = us/c13, ud/c13
    c23 = math.sqrt(1-s23*s23)
    cd = -s12*c23-c12*s23*w
    td = s12*s23-c12*c23*w
    tb = c23*c13
    z_ckm = -ud*w/(cd*cb)
    beta_ckm = cmath.phase(-cd*cb/(td*tb))
    assert abs(z_ckm-apex) < 3e-15 and abs(beta-beta_ckm) < 3e-15
    assert abs(us-mass_root(1/19.7, .973*52., 1.025*ub)) < 2e-15
    return dict(Vus=us, Vub=ub, Vud=ud, Vcb=cb, delta_deg=math.degrees(cmath.phase(w)),
                beta_deg=math.degrees(beta), gamma_deg=math.degrees(gamma))


def project(branch, n, data, q, tree_ub):
    beta, _, R = angles(branch, n, data)
    qmin, qmax = q
    assert qmax < .05
    # Intentionally enlarged R1755 envelopes. D is d/s just below the color threshold.
    d0, d1 = math.exp(-.0002)/(19.7+n*.9), math.exp(.0002)/(19.7-n*.9)
    umax = .05
    s23max = qmax/math.sqrt(1-umax**2)
    c23min = math.sqrt(1-s23max*s23max)
    # x_H² >= D sqrt(1 - x_H² - u_H²) for every nonnegative bottom mass.
    h = 2*d0*(1-umax**2)/(math.sqrt(d0*d0+4*(1-umax**2))+d0)
    xmin = math.sqrt(h) - .00015
    tanmin = xmin/math.sqrt(1-xmin*xmin)
    lower = tanmin*c23min*qmin*R[0]/(1+s23max*s23max*R[1])
    # Converse, after imposing the direct Vub upper edge.
    t_bound = tree_ub/(tanmin*c23min*qmin - tree_ub*s23max**2)
    assert t_bound < 1
    phase_gap = max(0., 2*(beta[0]-math.asin(t_bound)))
    # A narrower projection with an explicit, unsigned 0 <= (b/s)_H <= 100.
    x1 = math.sqrt(d1+umax**2*100*math.sqrt(2))+.00015
    x0, u0, u1 = xmin, 0., umax
    for _ in range(40):
        c_low = math.sqrt(1-x1*x1-u1*u1)
        c_high = math.sqrt(1-x0*x0-u0*u0)
        u0 = max(u0, x0/c_high*c23min*qmin*R[0]/(1+s23max*s23max*R[1]) - 1e-13)
        u1 = min(u1, x1/c_low*s23max*R[1]/(1-s23max*s23max*R[1]) + 1e-13)
        uH1 = min(umax, 1.06*u1+.00001)
        cH0 = math.sqrt(1-(x1+.00015)**2-uH1*uH1)
        cH1 = math.sqrt(1-(x0-.00015)**2)
        x0 = max(x0, math.sqrt(cH0*d0)-.00015-1e-13)
        x1 = min(x1, math.sqrt(cH1*d1+uH1*uH1*100/math.sqrt((1+cH0)/2))+.00015+1e-13)
    return dict(beta_deg=list(map(math.degrees,beta)), apex_radius=list(R),
                any_positive_bottom=dict(min_Vus=xmin,min_Vub=lower,
                    max_apex_radius_if_Belle_ub=t_bound,
                    min_extra_golden_phase_deg_if_Belle_and_reflected=math.degrees(phase_gap)),
                bottom_below_100=dict(Vus=[x0,x1],Vub=[u0,u1]))


def tuning_table():
    # Last-place half-unit is retained to bound the printed f_pi sensitivity.
    S, r, B = 27.178, .4556, 53.94
    def result(S, r, B):
        return mass_root(2/((1+r)*S), B*.973, .00393*1.025)
    centre = result(S, r, B)
    out = {}
    for label, k in [("total_fpi",1.),("nuclear_component_linear",.03/.13)]:
        span = [result(S+.0015*k, r+.00005*k, B-.055*k),
                result(S-.0015*k, r-.00005*k, B+.055*k)]
        out[label] = dict(range=span, max_shift=max(abs(x-centre) for x in span))
    return dict(central=centre, **out)


def compute():
    b = read(ARCH/"1707"/"inputs.json")
    tree = read(ARCH/"1702"/"inputs.json")["direct_semileptonic_moduli"]
    cb, *cb_errors = tree["Vcb"]
    ub, *ub_errors = tree["Vub"]
    outputs = {}
    for branch in ["acute_2beta","reflected_2beta"]:
        rows = {"central_interior_model": central(branch,b,cb)}
        for n in [1,2]:
            dq, du = n*math.hypot(*cb_errors), n*math.hypot(*ub_errors)
            row = project(branch,n,b,(cb-dq,cb+dq),ub+du)
            bounds = row["bottom_below_100"]
            assert all(bounds[key][0] <= rows["central_interior_model"][key] <= bounds[key][1]
                       for key in ["Vus","Vub"])
            if branch == "reflected_2beta":
                assert row["any_positive_bottom"]["min_Vub"] > 2*(ub+du)
            rows[f"width_{n}"] = row
        outputs[branch] = rows
    return dict(mass_inputs={"QCDSF_s_over_d":[19.7,.9], "conditional_bottom_over_s":[0,100],
                            "illustrative_bottom_at_H":52*.973}, direct_Bpi_2width_upper=ub+2*math.hypot(*ub_errors),
                branches=outputs, MILC_pion_normalisation_mass_side=tuning_table())


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    args = p.parse_args()
    out = compute()
    path = HERE/"independent_results.json"
    if args.check:
        assert read(path) == out
        print("1755 two CP branches, independent QCD mass bound, B-pion holdout and scale audit verified")
    else:
        with path.open("x",encoding="utf8") as f:
            json.dump(out,f,indent=2,allow_nan=False)
            f.write("\n")
