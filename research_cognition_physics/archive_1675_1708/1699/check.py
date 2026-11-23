"""One-loop spectral thresholds and two-Higgs projections of the round 1698 point."""
import argparse
import json
import math
from pathlib import Path
import runpy
from fractions import Fraction

import numpy as np

HERE = Path(__file__).resolve().parent
LOOP = 16 * math.pi**2


def read(path):
    return json.loads((HERE / path).read_bytes())


cfg = read("inputs.json")
spec = read(cfg["shared_spectrum"])["shared_witness"]
weak = read(cfg["gauge_yukawa_input"])
quartic = read(cfg["quartic_input"])
soft = spec["soft_GeV"]
cutoff = cfg["comparison_scale_GeV"]
thresholds = [(m, np.array(b)) for key, b in [
    ("colored_fermions_GeV", [4/15, 0, 2/3]), ("colored_scalars_GeV", [2/15, 0, 1/3]),
    ("weak_fermions_GeV", [2/5, 2/3, 0]), ("weak_scalars_GeV", [1/5, 1/3, 0])]
              for m in spec[key]]


def neutral_port():
    prev = runpy.run_path(str(HERE / cfg["shared_graph"]))
    dc, cc, hc = prev["neutral_graph"](read(cfg["shared_inputs"]))
    b, t = spec["neutral_b"], spec["neutral_t"]
    d, c, h = math.sqrt(b)*dc, b*math.sqrt(t)*cc, b*t*hc
    mass = np.block([[0*c, c], [c.T, h]])
    _, v = np.linalg.eigh(mass.conj().T @ mass)
    u = read("../1492/inputs.json")["v246_GeV"] / math.sqrt(2)
    y = np.column_stack([d, 0*d]) @ v / u
    for lo in (0, 2):
        assert abs(np.sum(abs(y[:, lo:lo+2])**2) -
                   sum(spec["neutral_pole_Yukawa_squared"][lo:lo+2])) < 3e-8
    return y


Y0 = neutral_port()


def beta(z, mu, susy, keep_extra=True, sterile_susy_delay=False):
    # The last four columns form the active heavy-neutrino Yukawa matrix.
    g = z[:3].real
    yt, yb, lam = z[3:6].real
    e = z[6:9]
    n = z[9:].reshape(3, 4)
    p = np.outer(e, e.conj())
    if sterile_susy_delay and susy:
        mask = np.array(spec["neutral_scalar_poles_if_common_shift_GeV"]) <= mu
        n = n * mask
    s = n @ n.conj().T
    t, b, trp, trs = yt*yt, yb*yb, np.trace(p).real, np.trace(s).real
    g1, g2, g3 = g*g
    gy = 3*g1/5
    bg = np.array([33/5, 1, -3.] if susy else [41/10, -19/6, -7.])
    if keep_extra:
        for mass, extra in thresholds:
            if mass < mu:
                bg += extra
    if susy:
        top = 6*t+b+trs-13*g1/15-3*g2-16*g3/3
        bot = 6*b+t+trp-7*g1/15-3*g2-16*g3/3
        de = (3*p + s) @ e + (3*b+trp-9*g1/5-3*g2)*e
        dn = (3*s + p) @ n + (3*t+trs-3*g1/5-3*g2)*n
        dl = 0.
    else:
        x = 3*t+3*b+trp+trs
        top = x+1.5*(t-b)-17*gy/12-9*g2/4-8*g3
        bot = x+1.5*(b-t)-5*gy/12-9*g2/4-8*g3
        de = 1.5*(p-s) @ e + (x-15*gy/4-9*g2/4)*e
        dn = 1.5*(s-p) @ n + (x-3*gy/4-9*g2/4)*n
        dl = (24*lam**2+(4*x-9*g2-3*gy)*lam
              -2*(3*t*t+3*b*b+trp*trp+np.trace(s@s).real)
              +3/8*(2*g2*g2+(g2+gy)**2))
    return np.r_[bg*g**3, yt*top, yb*bot, dl, de, dn.ravel()] / LOOP


def transport(z, low, high, susy, keep_extra=True, delay=False, step=.1):
    # Thresholds split the flow exactly. All step-dependent coefficients use the segment midpoint.
    midpoint = math.sqrt(low*high)
    steps = max(1, math.ceil(math.log(high/low)/step))
    dt = math.log(high/low)/steps
    def rhs(v):
        return beta(v, midpoint, susy, keep_extra, delay)
    for k in range(steps):
        a = rhs(z)
        b = rhs(z+dt*a/2)
        c = rhs(z+dt*b/2)
        d = rhs(z+dt*c)
        z = z+dt*(a+2*b+2*c+d)/6
        if max(abs(z[3:5])) > math.sqrt(cfg["perturbativity_Y_squared_upper"]):
            return z, low*math.exp((k+1)*dt)
    return z, high


def lower(neutrals=True, extra=True, tau_index=2, step=.1, top_z=0., higgs_z=0.):
    v = np.zeros(21, complex)
    v[:5] = [weak[key][0] for key in ("g1_GUT", "g2", "g3", "yt", "yb")]
    v[3] += top_z*weak["yt"][1]
    v[5] = (quartic["lambda_paper"][0]+higgs_z*quartic["lambda_paper"][1]) / quartic["lambda_paper_factor"]
    v[6+tau_index] = weak["ytau"][0]
    mz = weak["mu_ref_GeV"]
    events = sorted(set([mz, soft] + [m for m,_ in thresholds if mz<m<soft]
                        + spec["neutral_poles_GeV"]))
    for a,b in zip(events,events[1:]):
        if neutrals:
            for j,m in enumerate(spec["neutral_poles_GeV"]):
                if m == a:
                    v[9:].reshape(3,4)[:,j] = Y0[:,j]
        v,_ = transport(v,a,b,False,extra,step=step)
    return v


def upper(v, t, delay=False):
    v = v.copy()
    # One-loop MS-to-DR inverse-alpha offset used in the gauge witness.
    v[:3] = 1/np.sqrt(1/v[:3].real**2-np.array([0,2,3])/(48*math.pi**2))
    cs = 1/math.hypot(1,t)
    sn = t*cs
    v[3] /= sn
    v[4] /= cs
    v[6:9] /= cs
    v[9:] /= sn
    if np.max(abs(v[3:5])) >= math.sqrt(cfg["perturbativity_Y_squared_upper"]):
        return {"endpoint_GeV":soft, "perturbative":False, "Y_at_soft":v[3:5].real.tolist()}
    events = sorted(set([soft,cutoff] + [m for m,_ in thresholds if soft<m<cutoff] +
                        list(spec["neutral_scalar_poles_if_common_shift_GeV"])))
    for a,b in zip(events,events[1:]):
        v,end = transport(v,a,b,True,delay=delay)
        if end < b or max(abs(v[3:5])) >= math.sqrt(cfg["perturbativity_Y_squared_upper"]):
            return {"endpoint_GeV":end, "perturbative":False, "Y":v[3:5].real.tolist()}
    return {"endpoint_GeV":cutoff, "perturbative":True, "Y":v[3:5].real.tolist(),
            "tau_norm":float(np.linalg.norm(v[6:9])),
            "ratio_top_bottom":float((v[3]/v[4]).real)}


def scalar_projection(z):
    g1,g2 = z[:2].real
    yt,yb,lam = z[3:6].real
    D = (g2*g2+3*g1*g1/5)/8
    coeff = 3*yt**4/LOOP
    stop_max = 3*coeff
    def angles(shift):
        a = (lam-shift)/D
        if not (0 <= a < 1):
            return []
        f=math.sqrt((1+math.sqrt(a))/(1-math.sqrt(a)))
        return [1/f,f]
    ratio=yt/yb
    equal_needed = lam-D*((1-ratio**2)/(1+ratio**2))**2
    return {"g1_GUT":g1,"g2":g2,"gY":math.sqrt(.6)*g1,"g3":z[2].real,
            "yt":yt,"yb":yb,"lambda":lam,"D":D,"stop_max":stop_max,"R":ratio,
            "tan_at_max_stop":angles(stop_max),"tan_at_zero_stop":angles(0.),
            "tan1_extra_beyond_max_stop":lam-stop_max,
            "equal_at_soft_total_shift":equal_needed,
            "equal_at_soft_Xt_over_M": (math.sqrt(6*(1-math.sqrt(1-equal_needed/stop_max)))
                                       if 0<=equal_needed<=stop_max else None)}


def finite_quartic(p, tan_beta, mix_squared):
    """Central degenerate ordinary-MSSM threshold; Bagnaschi et al. Eqs. 8--13.

    Their lambda multiplies (H^dagger H)^2/2. Electroweak fermion mass ratios
    are +1, stop soft masses are equal, so all the loop functions below are one.
    """
    gy, gw, yt = p["gY"], p["g2"], p["yt"]
    cs = 1/math.hypot(1,tan_beta)
    sn = tan_beta*cs
    c = cs**2-sn**2
    gauge = gy**2+gw**2
    tree = p["D"]*c*c
    reg = -gy**4/4 -gy**2*gw**2/2 -(3/4-c*c/6)*gw**4
    heavy_H = -3/16*gauge**2*(4*sn*cs*c)**2
    stop = (6*yt**4*mix_squared*(1-mix_squared/12)
            +yt**2*mix_squared*gauge*(3*c-c*c)/4)
    d1,u1,d2,u2 = gy*cs,gy*sn,gw*cs,gw*sn
    ewinos = (-7/12*(d1**4+u1**4)-9/4*(d2**4+u2**4)
              -3/2*d1**2*u1**2-7/2*d2**2*u2**2
              -8/3*d1*u1*d2*u2-7/6*((d1*d2)**2+(u1*u2)**2)
              -1/6*((d1*u2)**2+(u1*d2)**2)
              -4/3*(d1*u2+u1*d2)*(d1*d2+u1*u2)
              +2/3*d1*u1*(2*tree-2*gy**2)+2*d2*u2*(2*tree-2*gw**2)
              +2*tree*(gy**2/3+gw**2))
    spectator = 0.
    for family,bY,b2 in [("weak",2/3,2/3),("colored",4/9,0.)]:
        for mf,mb in zip(spec[family+"_fermions_GeV"],spec[family+"_scalars_GeV"]):
            spectator += (bY*gy**4+b2*gw**4)*math.log(mb/max(mf,soft))/(4*LOOP)
    return tree + (reg+heavy_H+stop+ewinos)/(2*LOOP) + c*c*spectator


def bisect(fun, l, h):
    fl=fun(l)
    assert fl*fun(h)<0
    for _ in range(60):
        m=(l+h)/2
        if fun(m)*fl>0:
            l=m
        else:
            h=m
    return (l+h)/2


def root_certificate(p, shift, upper_envelope, positive_c):
    """Exact sign variations for the fixed numerical coefficients, after clearing (1+t^2)^4."""
    F=Fraction
    Y,W,T,D,A,L=map(F,(p["gY"]**2,p["g2"]**2,p["yt"]**2,p["D"],
                       LOOP,p["lambda"]-shift))
    B=F(finite_quartic(p,0.,0.) - p["D"])  # used only for a consistency assertion below
    speccoef=0.
    for fam,yy,ww in [("weak",2/3,2/3),("colored",4/9,0.)]:
        for f,b in zip(spec[fam+"_fermions_GeV"],spec[fam+"_scalars_GeV"]):
            speccoef+=(yy*float(Y)**2+ww*float(W)**2)*math.log(b/max(soft,f))/(4*LOOP)
    def poly_value(t):
        d=1+t*t
        c,v=F(1-t*t,d),F(2*t,d)
        q = 6 if upper_envelope else 0
        if upper_envelope and not positive_c:
            q += (Y+W)*(3*c-c*c)/(4*T)
        reg= -Y*Y/4-Y*W/2-(F(3,4)-c*c/6)*W*W
        h= -F(3,16)*(Y+W)**2*(2*v*c)**2
        st=6*T*T*q*(1-F(q,12)) + T*q*(Y+W)*(3*c-c*c)/4
        mix=(-(F(7,12)*Y**2+F(9,4)*W**2)*(1-v*v/2)
             -(F(3,2)*Y**2+F(7,2)*W**2)*v*v/4
             -F(2,3)*Y*W*v*v-F(7,6)*Y*W*(1-v*v/2)
             -Y*W*v*v/12-F(4,3)*Y*W*v
             +Y*v/3*(2*D*c*c-2*Y)+W*v*(2*D*c*c-2*W)
             +2*D*c*c*(Y/3+W))
        value=(D+F(speccoef))*c*c+(reg+h+st+mix)/(2*A)
        if t==0 and not upper_envelope:
            assert abs(float(value-D-B))<1e-14
        return (value-L)*d**4
    # Integer forward differences give the exact Newton polynomial passing through j=0,...,8.
    vals=[poly_value(k) for k in range(9)]
    coeff=[F(0)]*9
    base=[F(1)]
    for k in range(9):
        for j,x in enumerate(base):
            coeff[j]+=vals[0]*x
        vals=[b-a for a,b in zip(vals,vals[1:])]
        nxt=[F(0)]*(len(base)+1)
        for j,x in enumerate(base):
            nxt[j]-=k*x/(k+1)
            nxt[j+1]+=x/(k+1)
        base=nxt
    assert sum(x*10**j for j,x in enumerate(coeff)) == poly_value(10)
    transformed=[F(0)]*9
    for j,a in enumerate(coeff):
        if positive_c:  # (1+w)^8 P(w/(1+w)), 0<t<1
            for k in range(9-j):
                transformed[j+k]+=a*math.comb(8-j,k)
        else:  # P(1+w), t>1
            for k in range(j+1):
                transformed[k]+=a*math.comb(j,k)
    signs=[v>0 for v in transformed if v]
    return sum(a!=b for a,b in zip(signs,signs[1:]))


def inequalities(p, z):
    # Safe comparison domains for all t in the indicated intervals, no sampling.
    dt=math.log(cutoff/soft)/LOOP
    gdr2=1/(1/z[:3].real**2-np.array([0,2,3])/(48*math.pi**2))
    g12max=1/(1/gdr2[:2]-2*np.array([8.6,3])*dt)
    G=13*g12max[0]/15+3*g12max[1]+16*gdr2[2]/3
    low=.50
    w0=(low/math.hypot(1,low)/p["yt"])**2
    upper_time=LOOP/(2*G)*math.log((6/G-1/(4*math.pi))/(6/G-w0))
    # Prove preservation of the comparison box from a first exit contradiction.
    lo,hi=2.,6.
    a,b,e,n=1.22,.15,.15,.30
    sn=lo/math.hypot(1,lo)
    cs=1/math.hypot(1,hi)
    x0=np.array([p["yt"]/sn,p["yb"]/cs,np.linalg.norm(z[6:9])/cs,
                 np.linalg.norm(z[9:])/sn])
    B=np.array([6*a*a+b*b+n*n,6*b*b+a*a+e*e,
                4*e*e+3*b*b+n*n,4*n*n+e*e+3*a*a])
    upp=x0*np.exp(dt*B)
    assert np.all(upp < [a,b,e,n])
    qmin=p["R"]/hi*math.exp(-dt*(5*b*b+e*e+2*g12max[0]/5))
    return {"all_tan_below":low,"threshold_crossed_by_GeV":soft*math.exp(upper_time),
            "max_gauge_drag":G,"all_high_tan_in":[lo,hi],
            "positive_flow_max_top_bottom_tau_neutrino":upp.tolist(),
            "UV_top_over_bottom_lower":qmin}


def main():
    a = lower()
    p = scalar_projection(a)
    candidates = {"at_soft_equal_sources":p["R"], "tan1":1.}
    for key in ("tan_at_max_stop","tan_at_zero_stop"):
        candidates.update({key+str(i):v for i,v in enumerate(p[key])})
    profiles={k:upper(a,t) for k,t in candidates.items()}
    pmin=lambda t:finite_quartic(p,t,0)
    def best_x(t):
        c=(1-t*t)/(1+t*t)
        return min(6.,max(0.,6+(p["gY"]**2+p["g2"]**2)*(3*c-c*c)/(4*p["yt"]**2)))
    # derivative of 6*y_t^4*u*(1-u/12) with respect to u equals y_t^4*(6-u).
    pmax=lambda t:finite_quartic(p,t,best_x(t))
    improved={}
    for label,fun in [("minimum",pmin),("maximum",pmax)]:
        improved[label+"_at_equal"]=fun(p["R"])
        improved[label+"_at_tan1"]=fun(1)
        for name,l,h in [("low",.03,.55),("high",1.,20.)]:
            improved[label+"_root_"+name]=bisect(lambda t:fun(t)-p["lambda"],l,h)
    eps=cfg["residual_quartic_halfwidth_test"]
    corners={str(r):[bisect(lambda t:fun(t)+r-p["lambda"],lo,hi)
                     for fun in (pmin,pmax) for lo,hi in ((.01,.6),(1.,40.))]
             for r in (-eps,eps)}
    counts={str(shift):[root_certificate(p,shift,edge,side)
                        for edge in (False,True) for side in (True,False)]
            for shift in (-eps,0.,eps)}
    assert all(v == [1,1,1,1] for v in counts.values())
    old=runpy.run_path(str(HERE/"../1678/check.py"))
    old_angle=bisect(lambda angle:old["residual"](angle,1),58.,62.)
    edges,_,states=old["parameters"](old_angle)
    te,normal,_=states[1]
    A,B,C,E,F,Zr,Zi=edges
    Z=complex(Zr,Zi)
    t=te*np.exp(1j*old["phi"])
    source=normal*np.array([[0,A,0],[B,3*C,t*Z],[0,t*E,t*F]])
    cs=1/math.hypot(1,improved["minimum_root_high"])
    sv=np.linalg.svd(source,compute_uv=False)
    proj=np.linalg.svd(cs*source,compute_uv=False)
    gj={"old_input_scale_GeV":1e16,"old_gamma_degrees":old_angle,
        "common_Hd_projection":cs,
        "maximum_relative_change_in_normalized_charged_singular_values":
            float(max(abs(proj/proj[0]/(sv/sv[0])-1))),
        "new_running_Ytau_over_Yb_at_high_zero_stop":
            (lambda u:u["tau_norm"]/u["Y"][1])(upper(a,improved["minimum_root_high"]))}
    data_corners=[]
    for it in (-1.,1.):
        for ih in (-1.,1.):
            q=scalar_projection(lower(top_z=it,higgs_z=ih))
            data_corners.append({"top_widths":it,"quartic_widths":ih,"lambda_at_soft":q["lambda"],
                                 "equal_sources_extra_needed_at_zero_mixing":
                                     q["lambda"]-finite_quartic(q,q["R"],0.)})
    result={"soft_interface":p,"running_examples":profiles,
            "ordinary_finite_matching":improved,"residual_quartic_tests":corners,
            "Descartes_root_bounds_minlow_minhigh_maxlow_maxhigh":counts,
            "GJ_common_down_projection_test":gj,
            "independent_top_and_quartic_marginal_corners":data_corners,
            "analytic_branch_tests":inequalities(p,a),
            "without_neutral":scalar_projection(lower(neutrals=False)),
            "SM_only_below_soft":scalar_projection(lower(neutrals=False,extra=False)),
            "delayed_neutral_SUSY_at_equal":upper(a,p["R"],delay=True),
            "max_endpoint_quartic_change_under_two_tau_permutations":max(
                abs(lower(tau_index=i)[5].real-p["lambda"]) for i in (0,1)),
            "RK4_halfstep_lambda_error":abs(lower(step=.05)[5].real-p["lambda"])}
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=json.dumps(main(),ensure_ascii=False,indent=2)+"\n"
    out=HERE/"results.json"
    if args.check:
        assert json.loads(result)==json.loads(out.read_bytes())
        print("1699 verified")
    else:
        out.write_text(result,encoding="utf-8")
