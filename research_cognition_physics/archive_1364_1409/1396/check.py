"""Top-only one-loop SM, one common matching scale. All energies in GeV."""
import argparse
import importlib.util
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
INPUT = json.loads((HERE/"inputs.json").read_text(encoding="utf-8-sig"))
MATCH = json.loads((HERE/INPUT["matching"]).read_text(encoding="utf-8-sig"))
PARENT = json.loads((HERE/INPUT["top_model"]).read_text(encoding="utf-8-sig"))
B = np.array(PARENT["beta_g"])
C = np.array(PARENT["C_gut_norm"])
MZ = MATCH["mu_ref_GeV"]
TMAX = math.log(INPUT["window_end_GeV"]/MZ)
K = 8 * math.pi**2
NODES, WEIGHTS = np.polynomial.legendre.leggauss(40)
COUPLINGS = np.array([MATCH[key] for key in ["g1_GUT", "g2", "g3", "yt"]])
GAUGE, ERROR = COUPLINGS[:3, 0], COUPLINGS[:3, 1]
Y, DY = MATCH["yt"]
OBS, DOBS = np.array(INPUT["lambda_paper"]) / INPUT["lambda_paper_factor"]


def gauge_sq(t, g):
    return g**2 / (1-np.asarray(t)[..., None] * B*g**2/K)


def top_sq(t, g, y):
    # Same analytic kernel as 1395, batched for all quartic integration nodes.
    t = np.atleast_1d(t)
    def kernel(times):
        return np.exp(np.sum(np.log1p(-times[...,None] * B*g**2/K) * (C/B), axis=-1))
    integral = t/2 * (kernel(t[:,None] * (NODES+1)/2) @ WEIGHTS)
    return kernel(t) / (1/y**2 - PARENT["a"]*integral/K)


def zero(fun, left, right):
    lo_val, hi_val = fun(left), fun(right)
    assert lo_val*hi_val < 0, (lo_val, hi_val)
    for _ in range(42):
        mid = (left+right)/2
        if fun(mid)*lo_val > 0:
            left = mid
        else:
            right = mid
    return (left+right)/2


def crossing(g, y, weak=None):
    compare_g = g if weak is None else weak
    return zero(lambda t:top_sq(t,g,y)[0] - gauge_sq(t,compare_g)[1],0,TMAX)


def terms(t, g, y, envelope):
    if not envelope:
        sq = gauge_sq(t,g)
        q, w, u = 0.6*sq[:,0],sq[:,1],top_sq(t,g,y)
        bias = 12*u-9*w-3*q
        constant = -6*u**2 + 3/8*(2*w**2+(w+q)**2)
        return bias,bias,constant,constant,None,None
    lower,upper = gauge_sq(t,GAUGE-ERROR),gauge_sq(t,GAUGE+ERROR)
    ul,uu = top_sq(t,GAUGE+ERROR,Y-DY),top_sq(t,GAUGE-ERROR,Y+DY)
    ql,qu,wl,wu = 0.6*lower[:,0],0.6*upper[:,0],lower[:,1],upper[:,1]
    return (-9*wu-3*qu,-9*wl-3*ql,3/8*(2*wl**2+(wl+ql)**2),
            3/8*(2*wu**2+(wu+qu)**2),ul,uu)


def transport(start, stop, initial, g, y, bound=None, steps=192):
    grid = np.linspace(start,stop,2*steps+1)
    low_b,high_b,low_c,high_c,low_u,high_u = terms(grid,g,y,bound)
    def beta(value,i):
        if bound == "upper":
            linear = max(low_b[i]*value,high_b[i]*value)
            constant = high_c[i]
        elif bound == "lower":
            linear = min(low_b[i]*value,high_b[i]*value)
            constant = low_c[i]
        else:
            linear,constant = low_b[i]*value,low_c[i]
        if bound:
            # 12*u*lambda - 6*u^2 is concave on the same uncertain top interval.
            def yukawa_term(u):
                return 12*u*value - 6*u*u
            if bound == "upper":
                constant += yukawa_term(min(high_u[i],max(low_u[i],value)))
            else:
                constant += min(yukawa_term(low_u[i]),yukawa_term(high_u[i]))
        return (24*value**2+linear+constant)/(2*K)
    step = (stop-start)/steps
    value = initial
    for i in range(0,len(grid)-1,2):
        a = beta(value,i)
        b = beta(value+step*a/2,i+1)
        c = beta(value+step*b/2,i+1)
        d = beta(value+step*c,i+2)
        value += step*(a+2*b+2*c+d)/6
    return float(value)


def coupled(g,y,steps=192):
    t = crossing(g,y)
    return t,transport(t,0,0.,g,y,steps=steps)


def one_sided_certificate():
    """Outward decimal comparison. Natural-log time endpoints are exact rationals."""
    from decimal import Decimal, localcontext

    with localcontext() as ctx:
        ctx.prec = 48
        D = Decimal
        pad = D("1e-24")
        klo, khi = D("78.95683520871486"), D("78.95683520871488")
        coeff = (D(".85"),D("2.25"),D(8))
        betas = (D("4.1"),-D(19)/D(6),-D(7))
        ranges = []
        for key in ("g1_GUT","g2","g3","yt"):
            central,width = map(lambda x:D(str(x)),MATCH[key])
            ranges.append((central-width,central+width))
        n = 16384
        end = D("16.25")
        step = end/D(n)

        def values(k, j):
            t = step*j
            b = betas[k]
            x,y = (g*g for g in ranges[k])
            kmin,kmax = (khi,klo) if b>0 else (klo,khi)
            return (x/(1-t*b*x/kmin)-pad, y/(1-t*b*y/kmax)+pad)

        def analytic_step(p,rate,intercept):
            growth=(rate*step).exp()
            return p*growth-intercept*(growth-1)/rate

        lower_p,upper_p = 1/ranges[3][1]**2-pad,1/ranges[3][0]**2+pad
        previous = [values(k,0) for k in range(3)]
        cells=[]
        for j in range(1,n+1):
            next_g=[values(k,j) for k in range(3)]
            lower_g=[min(a[0],b[0]) for a,b in zip(previous,next_g)]
            upper_g=[max(a[1],b[1]) for a,b in zip(previous,next_g)]
            amin=sum(a*b for a,b in zip(coeff,lower_g))/khi-pad
            amax=sum(a*b for a,b in zip(coeff,upper_g))/klo+pad
            lp=analytic_step(lower_p,amin,D("4.5")/klo+pad)-pad
            up=analytic_step(upper_p,amax,D("4.5")/khi-pad)+pad
            assert lp>0
            cells.append((1/max(up,upper_p)-pad,1/min(lp,lower_p)+pad,
                          lower_g[1],upper_g[1],D(".6")*lower_g[0],D(".6")*upper_g[0]))
            lower_p,upper_p,previous=lp,up,next_g
        margin=1/upper_p-previous[1][1]-pad
        assert margin>0
        # All true crossings occur later than end. The backward solution initialized
        # at zero there is a lower subsolution because beta_lambda(0) stays negative.
        lower=D(0)
        speed=D(".12")
        for ul,uu,wl,wu,ql,qu in reversed(cells):
            l,r=lower-speed*step,lower+speed*step
            def part(u,x):
                return 12*u*x-6*u*u
            umax=max(ul,min(uu,r))
            high=24*max(l*l,r*r)+part(umax,r)+max(-(9*wl+3*ql)*l,-(9*wu+3*qu)*l)
            high+=D(".375")*(2*wu**2+(wu+qu)**2)+pad
            low=min(part(ul,l),part(uu,l))+min(-(9*wu+3*qu)*r,-(9*wl+3*ql)*r)-pad
            assert max(abs(low),abs(high))/(2*klo)<speed
            lower-=step*max(high/(2*klo),high/(2*khi))+pad
        assert lower>D(".140025")
        return {"steps":n,"fixed_log_time":float(end),"min_top_minus_weak_squared_at_fixed_time":float(margin),
                "all_inputs_lower_lambda":float(lower),"separation_from_matched_max":float(lower-D(".140025"))}


def calculate():
    t,low = coupled(GAUGE,Y)
    scale = MZ*math.exp(t)
    target_uv = transport(0,t,OBS,GAUGE,Y)
    assert abs(transport(t,0,target_uv,GAUGE,Y)-OBS)<1e-9
    assert abs(coupled(GAUGE,Y,384)[1]-low)<1e-9
    # Cross-check the inherited round's implementation at exactly this scale.
    previous = HERE/"../1395/check.py"
    spec = importlib.util.spec_from_file_location("inherited_top_map",previous)
    implementation = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(implementation)
    slope,intercept = implementation.affine(scale,GAUGE)
    assert abs(implementation.infer(Y,slope,intercept)-math.sqrt(gauge_sq(t,GAUGE)[1]))<1e-11

    tl,th = crossing(GAUGE+ERROR,Y-DY,GAUGE+ERROR),crossing(GAUGE-ERROR,Y+DY,GAUGE-ERROR)
    box = [transport(tl,0,0.,GAUGE,Y,"upper"),transport(th,0,0.,GAUGE,Y,"lower")]
    assert box[0]<low<box[1]
    # At crossings 16*pi^2*d ln(y/g2) has a negative numerator on the entire window.
    # The extreme comparison uses external weak running in the denominator conservatively.
    gr0 = gauge_sq(0,GAUGE-ERROR)[2]/gauge_sq(0,GAUGE+ERROR)[1]
    gr1 = gauge_sq(TMAX,GAUGE-ERROR)[2]/gauge_sq(TMAX,GAUGE+ERROR)[1]
    lower_color_ratio = min(gr0,gr1)
    assert 8*lower_color_ratio > 4.5 + 19/6
    # Matrix factors below 0.54 form a downward barrier through the complete matching window.
    g_end_low = gauge_sq(th,GAUGE-ERROR)
    assert 4.5*0.54**2 < 8*g_end_low[2]
    assert top_sq(th,GAUGE+ERROR,Y-DY)[0] > 0.54**2
    qmax = .6*gauge_sq(th,GAUGE+ERROR)[0]
    wmax = (GAUGE[1]+ERROR[1])**2
    beta_zero_max = -6*.54**4 + 3/8*(2*wmax**2 + (wmax+qmax)**2)
    assert beta_zero_max < 0
    variations = {}
    for index,name in enumerate(("g1","g2","g3","yt")):
        data = COUPLINGS[:,0].copy()
        data[index] += COUPLINGS[index,1]
        varied_t,varied_l = coupled(data[:3],data[3])
        variations[name] = {"scale_factor":math.exp(varied_t-t),"low_lambda_change":varied_l-low}
    rivals = {str(boundary):transport(t,0,boundary,GAUGE,Y)
              for boundary in INPUT["independent_matching_tests"]}
    zero_l = zero(lambda time:transport(0,time,OBS,GAUGE,Y),0,TMAX)
    old = INPUT["two_loop_table"]
    margins = []
    for index in range(2):
        top,dt = old["yt"][index]
        w,dw = old["g2"][index]
        margins.append([top-dt-w-dw,top+dt-w+dw])
    assert margins[0][0]>0 and margins[1][1]<0
    return {"one_sided_certificate":one_sided_certificate(),"reference_scale_GeV":MZ,"unique_scale_GeV":scale,
            "coupling_at_crossing":math.sqrt(gauge_sq(t,GAUGE)[1]),
            "predicted_lambda_MZ":low,"observed_matched_lambda_MZ":OBS,
            "abs_low_gap":abs(low-OBS),"needed_lambda_at_matching":target_uv,
            "scale_enclosure_GeV":[MZ*math.exp(tl),MZ*math.exp(th)],
            "low_lambda_continuous_enclosure":box,"target_low_range":[OBS-DOBS,OBS+DOBS],
            "separation_of_enclosures":max(box[0]-(OBS+DOBS),(OBS-DOBS)-box[1],0),
            "one_input_width_changes":variations,"low_lambda_nonzero_boundaries":rivals,
            "standalone_quartic_zero_GeV":MZ*math.exp(zero_l),
            "needed_y_over_g2_at_quartic_zero":math.sqrt(top_sq(zero_l,GAUGE,Y)[0]/gauge_sq(zero_l,GAUGE)[1]),
            "min_color_to_weak_squared_ratio":lower_color_ratio,
            "beta_lambda_zero_upper_times_16pi2":beta_zero_max,"two_loop_table_y_minus_g2_boxes":margins}


def close(x,z):
    if isinstance(x,dict):
        return x.keys()==z.keys() and all(close(v,z[k]) for k,v in x.items())
    return bool(np.allclose(x,z,rtol=2e-8,atol=1e-10))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    flags=parser.add_mutually_exclusive_group(required=True)
    flags.add_argument("--write-new",action="store_true")
    flags.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=calculate()
    path=HERE/"results.json"
    if args.write_new:
        with path.open("x",encoding="utf-8") as handle:
            json.dump(result,handle,indent=2,allow_nan=False)
            handle.write("\n")
    else:
        assert close(result,json.loads(path.read_text(encoding="utf-8-sig")))
        print("1396 simultaneous boundary prediction and scalar comparison envelopes verified")



