"""A texture-specific bound on the neutrino-option quartic and the SM infrared input."""
import argparse
from decimal import Decimal, localcontext
import importlib.util
import json
import math
from pathlib import Path

import numpy as np


ROOT=Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def spectral_inputs(data):
    def edges(row):
        return [row[0]-2*row[1],row[0]+2*row[-1]]
    ds,dt=data["KamLAND_dm21_eV2"],data["T2K_dm32_NO_eV2"]
    def norms(x,y):
        m2,m3=math.sqrt(x),math.sqrt(x+y)
        product=11*m2*m3/3
        total=math.sqrt(m2*m2+m3*m3+10*m2*m3)
        ta=(total+math.sqrt(total*total-4*product))/2
        return m2,m3,total,ta,product/ta
    return [norms(a,b) for a,b in [tuple([ds[0],dt[0]]), *zip(edges(ds),edges(dt))]]


def thresholds(i,j,overlap,ratio,phase_cos):
    f=5/(32*math.pi**2)
    if ratio==1:
        return -f*(i*i+j*j+(6+4*phase_cos)*overlap)
    L=math.log1p(ratio-1)/(ratio-1)
    H=2*ratio/(1+ratio)*L
    return -f*(i*i+j*j+(1-4*L-2*H*phase_cos)*overlap)


def proven_lower_quartic(old, log_end):
    """Forward scalar differential comparison including gauge/top one-width input boxes."""
    with localcontext() as ctx:
        ctx.prec=48
        D=Decimal
        pad=D("1e-24")
        # Upper round-off in log_end extends the certified time domain.
        time=D(str(math.ceil(log_end*100)/100))
        count=8192
        dt=time/count
        kmin,kmax=D("78.95683520871486"),D("78.95683520871488")
        b=[D("4.1"),-D(19)/6,-D(7)]
        coeff=[D(".85"),D("2.25"),D(8)]
        boxes=[]
        for key in ("g1_GUT","g2","g3","yt"):
            value,err=map(lambda v:D(str(v)),old.MATCH[key])
            boxes.append((value-err,value+err))
        def gauge(k,time_index):
            t=dt*time_index
            l,r=[x*x for x in boxes[k]]
            kb,kt=(kmax,kmin) if b[k]>0 else (kmin,kmax)
            return (l/(1-b[k]*t*l/kb)-pad,r/(1-b[k]*t*r/kt)+pad)
        def top_step(value,slope,loss):
            growth=(slope*dt).exp()
            return growth*value-loss*(growth-1)/slope
        pl=1/boxes[3][1]**2-pad
        ph=1/boxes[3][0]**2+pad
        value=(D(str(old.INPUT["lambda_paper"][0]))-D(str(old.INPUT["lambda_paper"][1])))/4
        floor=value
        rate=D(".09")
        prev=[gauge(k,0) for k in range(3)]
        for j in range(1,count+1):
            nxt=[gauge(k,j) for k in range(3)]
            gl=[min(x[0],y[0]) for x,y in zip(prev,nxt)]
            gh=[max(x[1],y[1]) for x,y in zip(prev,nxt)]
            al=sum(x*y for x,y in zip(coeff,gl))/kmax-pad
            ah=sum(x*y for x,y in zip(coeff,gh))/kmin+pad
            nl=top_step(pl,al,D("4.5")/kmin+pad)-pad
            nh=top_step(ph,ah,D("4.5")/kmax-pad)+pad
            assert min(pl,nl)>0
            yl,yh=1/max(ph,nh)-pad,1/min(pl,nl)+pad
            l,h=value-rate*dt,value+rate*dt
            assert l>0
            ql,qh=D(".6")*gl[0],D(".6")*gh[0]
            wl,wh=gl[1],gh[1]
            def top(u,x):
                return 12*u*x-6*u*u
            low=24*l*l+min(top(yl,l),top(yh,l))-(9*wh+3*qh)*h+D(".375")*(2*wl**2+(wl+ql)**2)-pad
            high=24*h*h+top(max(yl,min(yh,h)),h)-(9*wl+3*ql)*l+D(".375")*(2*wh**2+(wh+qh)**2)+pad
            assert max(abs(low),abs(high))/(2*kmin)<rate
            floor=min(floor,l)
            value+=dt*min(low/(2*kmin),low/(2*kmax))-pad
            prev,pl,ph=nxt,nl,nh
        return {"log_endpoint":float(time),"scale_ceiling_GeV":float(old.MZ)*math.exp(float(time)),
                "continuous_lower_bound":float(min(floor,value)),"steps":count}


def solve():
    cfg=read(ROOT/"inputs.json")
    ew=read(ROOT/cfg["electroweak"])
    mh,dmh=ew["atlas"]["mass_GeV"]
    u=cfg["vev_over_sqrt2_GeV"]
    central,lo,hi=spectral_inputs(read(ROOT/cfg["oscillations"]))
    light_lo,light_hi=cfg["light_coefficient_transport"]
    threshold_lo,threshold_hi=cfg["quadratic_threshold_to_higgs_mass_squared_transport"]
    limit=(mh+2*dmh)**2*8*math.pi**2*u*u/threshold_lo
    tmax=hi[2]*1e-9*light_hi
    tmin=lo[2]*1e-9*light_lo
    smaller_max=(limit/tmin)**(1/3)
    smaller_min=((mh-2*dmh)**2*8*math.pi**2*u*u/(threshold_hi*tmax))**(1/3)/10
    sum_norm_max=tmax**(2/3)*limit**(1/3)/(u*u)
    gram=8/11
    # The stronger split bound uses only log(x)/(x-1)<=1 and 2x log(x)/(x²-1)<=1.
    degenerate_abs=(5/(32*math.pi**2))*(1+max(0,10*gram-2)/4)*sum_norm_max**2
    split_abs=(5/(32*math.pi**2))*(1+max(0,5*gram-2)/4)*sum_norm_max**2
    spec=importlib.util.spec_from_file_location("rge1396",ROOT/cfg["rge_source"])
    old=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(old)
    low=proven_lower_quartic(old,math.log(smaller_max/old.MZ))
    assert low["continuous_lower_bound"]>1e-3
    def target(M):
        t=math.log(M/old.MZ)
        return old.transport(0,t,old.OBS,old.GAUGE,old.Y)
    witnesses={}
    for ratio in cfg["examples"]:
        for first in ("a","b"):
            t1,t2=central[3:5] if first=="a" else central[3:5][::-1]
            m1=(8*math.pi**2*u*u*mh**2/((t1+ratio**3*t2)*1e-9))**(1/3)
            m2=ratio*m1
            aa,bb=t1*1e-9*m1/u**2,t2*1e-9*m2/u**2
            dl=thresholds(aa,bb,gram*aa*bb,ratio,-.5)
            assert dl<0 and abs(dl)<max(degenerate_abs,split_abs)
            assert math.isclose((aa*m1**2+bb*m2**2)/(8*math.pi**2),mh**2,rel_tol=6e-14)
            # Separate the author's equal-scale finite match from the scale where lambda0 is quoted.
            witnesses[str(ratio)+first]={"heavy_masses_GeV":[m1,m2],"squared_norms":[aa,bb],
                "delta_lambda":dl,"needed_lambda_at_M1":target(m1)-dl}
    return {"central_t_eV":list(central[3:5]),"raw_sum_t_eV_outer":[lo[2],hi[2]],
        "smaller_mass_GeV_outer":[smaller_min,smaller_max],"sum_squared_norms_upper":sum_norm_max,
        "max_abs_delta_equal":degenerate_abs,"max_abs_delta_split":split_abs,
        "one_loop_SM_target":low,"observed_boundary_quartic":old.OBS,
        "reference_boundary_lower":low["continuous_lower_bound"]-max(degenerate_abs,split_abs),
        "examples":witnesses}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write-new",action="store_true")
    mode.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=solve()
    path=ROOT/"results.json"
    if args.write_new:
        with path.open("x",encoding="utf-8") as f:
            json.dump(result,f,indent=2)
            f.write("\n")
    else:
        assert result==read(path)
        print("1493: same-source quartic size, mass window and continuous one-loop SM floor verified.")
