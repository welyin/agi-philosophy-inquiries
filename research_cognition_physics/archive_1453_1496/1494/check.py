"""Third-family SM evolution and a conditional supersymmetric quartic boundary."""
import argparse
from decimal import Decimal, localcontext
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


CFG = read(HERE / "inputs.json")
MASS = read(HERE / CFG["gauge_yukawa"])
HIGGS = read(HERE / CFG["quartic"])
MZ = MASS["mu_ref_GeV"]
KEYS = ("yt", "yb", "ytau")
GAUGE = ("g1_GUT", "g2", "g3")
LOOP = 16 * math.pi**2


def beta(time, state):
    t, b, e, quartic = state
    g = np.array([MASS[key][0]**2 for key in GAUGE])
    g1, w, s = g / (1 - time * g * np.array([4.1, -19/6, -7]) * 2 / LOOP)
    q = 0.6 * g1
    trace = 3*t + 3*b + e
    yukawa = [2*t*(trace + 1.5*(t-b) - 17*q/12 - 2.25*w - 8*s),
              2*b*(trace + 1.5*(b-t) - 5*q/12 - 2.25*w - 8*s),
              2*e*(trace + 1.5*e - 15*q/4 - 2.25*w)]
    scalar = (24*quartic**2 + quartic*(4*trace - 9*w - 3*q)
              - 2*(3*t*t+3*b*b+e*e) + 3/8*(2*w*w+(w+q)**2))
    return np.array(yukawa + [scalar]) / LOOP


def central(scale):
    time = math.log(scale/MZ)
    step = time/256
    x = np.array([MASS[key][0]**2 for key in KEYS] +
                 [HIGGS["lambda_paper"][0]/HIGGS["lambda_paper_factor"]])
    for k in range(256):
        t = k*step
        a = beta(t,x)
        b = beta(t+step/2,x+a*step/2)
        c = beta(t+step/2,x+b*step/2)
        d = beta(t+step,x+c*step)
        x += step*(a+2*b+2*c+d)/6
    g = np.array([MASS[key][0]**2 for key in GAUGE[:2]])
    q, weak = g/(1 - time*g*np.array([4.1,-19/6])*2/LOOP)
    gauge_factor = (0.6*q + weak)/8
    ratio = math.sqrt(x[0]/x[1])
    prediction = gauge_factor*((ratio**2-1)/(ratio**2+1))**2
    result = dict(ratio=ratio, lambda_D=prediction, lambda_SM=x[3],
                  additive_threshold=x[3]-prediction, gauge_factor=gauge_factor)
    if 0 < x[3] < gauge_factor:
        shift = math.sqrt(x[3]/gauge_factor)
        large = math.sqrt((1+shift)/(1-shift))
        result["both_tanbeta_for_zero_threshold"] = [1/large,large]
        result["Yt_over_Yb"] = [ratio*large,ratio/large]
    return result


def continuous_certificate(cut_low,cut_high):
    # Two cooperative Yukawa comparison systems; a scalar comparison for lambda.
    with localcontext() as ctx:
        ctx.prec = 44
        D = Decimal
        pad = D("1e-23")
        loop = D("157.91367041742973790135")
        def data_box(pair):
            mid, err = map(lambda x: D(str(x)),pair)
            return mid-err,mid+err
        ybox = [data_box(MASS[k]) for k in KEYS]
        low = [a*a for a,b in ybox]
        high = [b*b for a,b in ybox]
        lam0,lam1 = data_box(HIGGS["lambda_paper"])
        lam0,lam1 = lam0/4,lam1/4
        gbox = [data_box(MASS[k]) for k in GAUGE]
        group = [D("4.1"),-D(19)/6,D(-7)]
        speeds = [D(".25"),D(".00015"),D(".000015")]
        scalar_speed = D(".09")
        count = 6144
        finish = D(str(math.ceil(math.log(cut_high/MZ)*100)/100))
        start = D(str(math.floor(math.log(cut_low/MZ)*100)/100))
        step = finish/count
        range_all = {}
        def update(name,value):
            current = range_all.setdefault(name,[value,value])
            current[0] = min(current[0],value)
            current[1] = max(current[1],value)
        def gauge(time):
            return [(a*a/(1-2*k*a*a*time/loop)-pad,b*b/(1-2*k*b*b*time/loop)+pad)
                    for (a,b),k in zip(gbox,group)]
        prev = gauge(D(0))
        for j in range(count):
            now = gauge((j+1)*step)
            gl = [min(x[0],y[0]) for x,y in zip(prev,now)]
            gh = [max(x[1],y[1]) for x,y in zip(prev,now)]
            def gauge_loss(g):
                q,w,s = D(".6")*g[0],g[1],g[2]
                return [17*q/12+D("2.25")*w+8*s,5*q/12+D("2.25")*w+8*s,
                        D("3.75")*q+D("2.25")*w]
            loss_lo,loss_hi = gauge_loss(gl),gauge_loss(gh)
            def yuk_increments(x,lower):
                u,v = [a-b*step for a,b in zip(x,speeds)],[a+b*step for a,b in zip(x,speeds)]
                assert min(u)>0
                def gain(xx):
                    t,b,e=xx
                    return [D("4.5")*t+D("1.5")*b+e,D("1.5")*t+D("4.5")*b+e,3*t+3*b+D("2.5")*e]
                first = [a-b for a,b in zip(gain(u),loss_hi)]
                last = [a-b for a,b in zip(gain(v),loss_lo)]
                inc=[]
                for lo,hi,fl,fh,speed in zip(u,v,first,last,speeds):
                    choices=[2*xx*ff/loop for xx in (lo,hi) for ff in (fl,fh)]
                    assert max(map(abs,choices))+pad<speed
                    inc.append(step*(min(choices) if lower else max(choices))+(-pad if lower else pad))
                return inc
            dl,dh=yuk_increments(low,True),yuk_increments(high,False)
            lows=[a-b*step for a,b in zip(low,speeds)]
            highs=[a+b*step for a,b in zip(high,speeds)]
            ql,qh=D(".6")*gl[0],D(".6")*gh[0]
            wl,wh=gl[1],gh[1]
            def scalar_increment(anchor,lower):
                l,h=anchor-scalar_speed*step,anchor+scalar_speed*step
                assert l>0
                bottom=D(24)*l*l-(9*wh+3*qh)*h+D(".375")*(2*wl*wl+(wl+ql)**2)
                top=D(24)*h*h-(9*wl+3*ql)*l+D(".375")*(2*wh*wh+(wh+qh)**2)
                for k,yl,yh in zip((3,3,1),lows,highs):
                    def term(y,lv):
                        return k*(4*y*lv-2*y*y)
                    bottom+=min(term(yl,l),term(yh,l))
                    top+=term(min(yh,max(yl,h)),h)
                assert max(abs(bottom),abs(top))/loop+pad<scalar_speed
                return step*(bottom if lower else top)/loop+(-pad if lower else pad)
            sl,sh=scalar_increment(lam0,True),scalar_increment(lam1,False)
            if (j+1)*step>=start:
                rr=(lows[0]/highs[1]).sqrt()-pad
                rr_hi=(highs[0]/lows[1]).sqrt()+pad
                g0,g1=(ql+wl)/8-pad,(qh+wh)/8+pad
                l,h=lam0-scalar_speed*step,lam1+scalar_speed*step
                assert rr>1 and 0<l<h<g0
                f0=g0*((rr*rr-1)/(rr*rr+1))**2-pad
                f1=g1*((rr_hi**2-1)/(rr_hi**2+1))**2+pad
                beta_hi=((1+(h/g0).sqrt())/(1-(h/g0).sqrt())).sqrt()+pad
                for name,left,right in (("ratio",rr,rr_hi),("lambda_SM",l,h),("lambda_D",f0,f1),
                                         ("D_minus_target",f0-h,f1-l),("Yt_over_Yb_lower",rr/beta_hi-pad,rr/beta_hi-pad)):
                    update(name,left)
                    update(name,right)
            low=[a+b for a,b in zip(low,dl)]
            high=[a+b for a,b in zip(high,dh)]
            lam0,lam1=lam0+sl,lam1+sh
            prev=now
        out={key:[float(minimum),float(maximum)] for key,(minimum,maximum) in range_all.items()}
        out["Yt_over_Yb_lower"]=out["Yt_over_Yb_lower"][0]
        out.update(steps=count,log_domain=[float(start),float(finish)])
        return out


def run():
    interval=CFG["scale_in_GeV"]
    inherited=read(HERE/CFG["neutrino_example_window"])["smaller_mass_GeV_outer"]
    assert interval[0]<inherited[0]<inherited[1]<interval[1]
    table=CFG["table2_at_1e7"]
    ratio=table["yt"][0]/table["yb"][0]
    den=(0.6*table["g1"][0]**2+table["g2"][0]**2)/8
    lo,hi=1000.,100000.
    assert central(lo)["additive_threshold"]*central(hi)["additive_threshold"]<0
    for _ in range(38):
        mid=math.sqrt(lo*hi)
        if central(mid)["additive_threshold"]>0:
            lo=mid
        else:
            hi=mid
    cert=continuous_certificate(*interval)
    assert cert["D_minus_target"][0]>0 and cert["Yt_over_Yb_lower"]>10
    return {"certificate":cert,"scale_GeV":interval,"central":{str(x):central(x) for x in [475000,1000000,3000000,10000000,12300000]},
            "central_lower_scale_crossing_GeV":[lo,hi],
            "published_two_loop_1e7":{"ratio":ratio,"lambda_D":den*((ratio**2-1)/(ratio**2+1))**2}}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    modes=parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--write-new",action="store_true")
    modes.add_argument("--check",action="store_true")
    args=parser.parse_args()
    outcome=run()
    path=HERE/"results.json"
    if args.write_new:
        with path.open("x",encoding="utf-8") as out:
            json.dump(outcome,out,indent=2)
            out.write("\n")
    else:
        assert outcome==read(path)
        print("1494: continuous threshold gap, mass ratios and lower-scale crossing verified.")
