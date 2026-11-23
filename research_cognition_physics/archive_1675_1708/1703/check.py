"""Transport the protected null column with the full charged matrices and actual pole events."""
import argparse
import json
import math
import runpy
from pathlib import Path

import numpy as np

ROOT=Path(__file__).resolve().parent
LOOP=16*math.pi**2


def read(s):
    return json.loads((ROOT/s).read_bytes())


def unpack(z):
    return np.array(z["real"])+1j*np.array(z["imag"])


def growth_gate(high,charged,ratio,soft,mz,poles,tan):
    # Bound both rays for any permitted singlet rotation/tree threshold, including close pairs.
    # Dropping positive Yukawa terms in the downward norm derivatives bounds their growth.
    high_N=np.linalg.norm(high)**2
    high_E=sum(np.square(charged))
    d1=math.log(1e7/soft)/LOOP
    d2=math.log(soft/min(poles))/LOOP
    d3=math.log(soft/mz)/LOOP
    n1=high_N*math.exp(2*d1*(3*.70**2+.6*.55**2))
    e1=high_E*math.exp(2*d1*(3*.70**2+1.8*.55**2))
    n2_try,e2_try=.10,.002
    e2=e1/(1+tan*tan)*math.exp(2*d3*(9*.70**2/4+9*.55**2/4+1.5*n2_try))
    n2=n1*math.exp(2*d2*(9*.70**2/4+9*.55**2/20+1.5*e2_try))
    assert n2<n2_try and e2<e2_try
    # Common scale terms cancel in singular-value ratios.
    r=ratio*math.exp(d1*(3*e1+n1)+1.5*(d3*e2+d2*n2))
    factor=(1+r*r)/(1-r*r)
    assert r<.10
    ray=.5*(d1*e1+1.5*d3*e2)
    charged_ray=.5*factor*(d1*n1+1.5*d2*n2)
    return ray,charged_ray


def integrate_null(fr, gamma, dz, rep, step):
    frame=unpack(rep["charged_left_MZ"])
    svals=np.array(rep["Y_charged_at_10PeV"])
    y=np.zeros((4,3,4),complex)
    y[0,:,:3]=np.diag(fr["central"]("up"))
    down=fr["central"]("down")
    down[0]+=dz*fr["low"]["yukawa"]["down"][0][1]
    y[1,:,:3]=fr["low_ckm"](gamma,fr["angles"]("direct"))@np.diag(down)
    y[2,:,:3]=frame@np.diag(fr["central"]("charged_leptons"))
    track=[]
    for a,susy,dt,gg in fr["segments"](step):
        if a==fr["MS"]:
            cs=1/math.hypot(1,fr["TAN"])
            y/=np.array([fr["TAN"]*cs,cs,cs,fr["TAN"]*cs])[:,None,None]
        for j,m in enumerate(fr["fields"]["neutral_poles_GeV"]):
            if m==a:
                y[3,:,j]=fr["physics"]["Y0"][:,j]
        for i in range(0,len(gg)-1,2):
            u=fr["beta"](y,gg[i],susy)
            ym=y+dt*u/2
            v=fr["beta"](ym,gg[i+1],susy)
            yn=y+dt*v/2
            w=fr["beta"](yn,gg[i+1],susy)
            x=fr["beta"](y+dt*w,gg[i+2],susy)
            end=y+dt*(u+2*v+2*w+x)/6
            Es=[v[2]@v[2].conj().T for v in [end,yn,ym,y]]
            track.append((dt,susy,Es,end[3],y[3]))
            y=end
    np.testing.assert_allclose(y[3],unpack(rep["neutral_columns_10PeV"]),atol=3e-10)
    uc,sv,_=np.linalg.svd(y[3])
    n=uc[:,-1]
    assert sv[-1]/sv[-2]<1e-11
    ehigh=unpack(rep["charged_left_10PeV"])
    phigh=abs(ehigh.conj().T@n)**2
    rotation=np.linalg.qr(np.array([[1,2j,3],[2j,3,1j],[2,3,4j]],complex))[0]
    moved=np.linalg.svd(rotation.conj().T@y[3],full_matrices=True)[0][:,-1]
    np.testing.assert_allclose(abs(ehigh.conj().T@rotation@moved)**2,phigh,atol=2e-14)
    bound=growth_gate(y[3],svals,max(svals[:-1]/svals[1:]),fr["MS"],fr["MZ"],
                      fr["fields"]["neutral_poles_GeV"],fr["TAN"])
    max_null_error=0
    for dt,susy,Es,end,start in reversed(track):
        max_null_error=max(max_null_error,np.linalg.norm(end.conj().T@n))
        ce=1 if susy else -1.5
        def reverse(v,E):
            return ce*(E@v-v*np.vdot(v,E@v)/np.vdot(v,v))
        u=reverse(n,Es[0])
        v=reverse(n+dt*u/2,Es[1])
        w=reverse(n+dt*v/2,Es[2])
        x=reverse(n+dt*w,Es[3])
        n=n+dt*(u+2*v+2*w+x)/6
        n/=np.linalg.norm(n)
        max_null_error=max(max_null_error,np.linalg.norm(start.conj().T@n))
    plow=abs(frame.conj().T@n)**2
    np.testing.assert_allclose(sum(plow),1,atol=1e-12)
    assert np.max(np.abs(np.arccos(np.sqrt(plow))-np.arccos(np.sqrt(phigh))))<sum(bound)
    assert max_null_error<1e-11
    return phigh,plow,bound,max_null_error


def constraints(probs,bounds_r,windows):
    theta=np.arccos(np.sqrt(probs))
    lo=np.cos(np.minimum(math.pi/2,theta+sum(bounds_r)))**2
    hi=np.cos(np.maximum(0,theta-sum(bounds_r)))**2
    a,solar,atm=windows
    zlo,zhi=(1-np.sqrt(1-a))/2
    # The second reactor branch would require p_e1<=1-z<=zhi.
    assert lo[0]>zhi
    xlo=max(solar[0],1-hi[0]/(1-zhi))
    xhi=min(solar[1],1-lo[0]/(1-zlo))
    result=dict(protected_probability_envelopes=np.array([lo,hi]).T.tolist(),
                solar_reactor_compatible=bool(xlo<=xhi))
    if xlo<=xhi:
        num=lo[1]-xhi*(1-atm[0])-(1-xlo)*atm[1]*zhi
        denmax=2*math.sqrt(xhi*(1-xhi)*.25*zhi)
        denmin=2*math.sqrt(xlo*(1-xlo)*min(atm*(1-atm))*zlo)
        result["necessary_cos_delta_lower"]=num/(denmax if num>0 else denmin)
        result["conditional_s12sq"]=[xlo,xhi]
    return result


def kernel(s):
    q=1/math.sqrt(3)
    v=math.sqrt(7)/5
    x=q*(1-s*q)/(s-q)
    assert x>0
    return v*(1+x)/(x+v*v)


def compute():
    cfg=read("inputs.json")
    fr=runpy.run_path(str(ROOT/cfg["three_family_runner"]))
    past=read(cfg["charged_endpoints"])
    raw=fr["physics"]["Y0"].copy()
    phi=cfg["phase"]
    old=np.exp(1j*phi*np.array([2,0,1]))
    oriented=np.exp(1j*np.array([0,2*phi,3*phi+np.angle(3+2*np.exp(4j*phi))]))
    cases={}
    step=cfg["rk_log_step"]
    windows=np.array(read(cfg["angles"])["joint_windows"])
    for dz,kind in [(0,"central"),(-1,"lower_down")]:
        for k in [1,2]:
            anchor=(past["charged_endpoints"]["direct"][str(k)] if not dz
                    else past["down_one_width_lower"]["direct_k"+str(k)])
            for role in ["supplied_left","holomorphic_raw","unit_third","shared_third"]:
                h=1 if role!="shared_third" else kernel(anchor["large_modulus"])
                for _ in range(8):
                    if role=="supplied_left":
                        columns=raw.copy()
                    elif role=="holomorphic_raw":
                        columns=raw.conj()
                    else:
                        columns=((oriented/old)*np.array([1,1,h])).conj()[:,None]*raw.conj()
                    fr["physics"]["Y0"]=columns
                    rep=fr["solve"](k,"direct",seed=anchor["gamma_MZ_deg"],step=step,d_z=dz)
                    new=kernel(rep["large_modulus"]) if role=="shared_third" else 1
                    if abs(new-h)<3e-10:
                        break
                    h=new
                else:
                    raise RuntimeError("shared transfer did not converge")
                if role=="supplied_left":
                    assert abs(rep["gamma_MZ_deg"]-anchor["gamma_MZ_deg"])<3e-6
                uh,ul,bb,err=integrate_null(fr,rep["gamma_MZ_deg"],dz,rep,step)
                res=dict(gamma_quark_MZ_deg=rep["gamma_MZ_deg"],charged_gain=rep["large_modulus"],
                         third_neutral_modulus=h,p_high=uh.tolist(),p_MZ=ul.tolist(),
                         max_massless_orthogonality_residual=err,
                         projective_angles_bounds_rad=list(bb))
                res.update(constraints(uh,bb,windows))
                if dz==0 and k==2 and role=="shared_third":
                    finer=fr["solve"](k,"direct",seed=rep["gamma_MZ_deg"],step=step/2,d_z=dz)
                    _,almost,_,_=integrate_null(fr,finer["gamma_MZ_deg"],dz,finer,step/2)
                    res["half_logstep_change_in_probabilities"]=max(abs(almost-ul))
                cases[f"{kind}_k{k}_{role}"]=res
    return dict(tan_beta=fr["TAN"],solar_reactor_p_e1_range=[
        (1-windows[1,1])*(1+math.sqrt(1-windows[0,1]))/2,
        (1-windows[1,0])*(1+math.sqrt(1-windows[0,0]))/2],cases=cases)


def neat(y):
    if isinstance(y,dict):
        return {k:neat(v) for k,v in y.items()}
    if isinstance(y,list):
        return [neat(v) for v in y]
    return float(f"{y:.9g}") if isinstance(y,(float,np.floating)) else y


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--check",action="store_true")
    args=ap.parse_args()
    data=neat(compute())
    if args.check:
        assert data==read("results.json")
        print("1703: charged shooting, moving massless column and phase enclosures verified")
    else:
        (ROOT/"results.json").write_text(json.dumps(data,indent=2)+"\n",encoding="utf-8")
