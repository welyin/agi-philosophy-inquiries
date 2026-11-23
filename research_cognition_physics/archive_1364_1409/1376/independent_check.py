"""Find trace intersections by positive-cone linear algebra, then run the flavour matrices."""
import argparse
import json
from itertools import product
from math import exp, log, pi, sqrt
from pathlib import Path
import numpy as np


MZ, END = 91.1876,1e16
LOOP=16*pi*pi
G0=np.array([.461228,.65096,1.2123])
GE=np.array([.000026,.00004,.0046])
Y0=np.array([[7.04e-6,.00356,.967],[1.54e-5,.000306,.01630],
             [2.77713e-6,.000585042,.0099378]])
YE=np.array([[.15e-6,.00006,.004],[.02e-5,.000004,.00009],
             [.00036e-6,.000000075,.0000014]])
SM=np.array([-4.1,19/6,7])*2/LOOP
WEIGHTS=np.arange(-3,4)
# beta of inverse squares, including the GUT hypercharge normalization.
CHANGE=-np.array([.6*sum(2**2 for _ in WEIGHTS),sum(WEIGHTS**2),0])*2/(3*LOOP)
TRACE=np.array([[1.1,.9],[1.5,.5],[2,0]])
TMAX=log(END/MZ)


def crossings(mass,g):
    split=log(mass/MZ)
    initial=1/g**2
    ans=[]
    for side,slopes,intercept in (("below",SM,initial),("above",SM+CHANGE,initial-CHANGE*split)):
        # Unknowns are the two positive coefficients of the traces and log(mu/MZ).
        a,b,t=np.linalg.solve(np.column_stack((TRACE,-slopes)),intercept)
        if 0<=t<=TMAX and (side=="below" and t<=split or side=="above" and t>split) and min(a,b)>0:
            xx=initial+SM*t+CHANGE*max(0.,t-split)
            assert np.allclose(xx,TRACE@[a,b])
            ans.append((side,t,3*a/b))
    return ans


def flav_box():
    # Gluon damping is unchanged by a colour-singlet threshold. The weak
    # inverse square is smaller than without the threshold while perturbative.
    u=(Y0+YE)[:2]
    u0=(u**2).sum()
    l0=((Y0+YE)[2]**2).sum()
    min_g3=G0[2]-GE[2]
    rate=SM[2]*min_g3**2
    eint=7/rate*(1-(1+rate*TMAX)**(-1/7))
    d=1-(9*u0/LOOP)*eint
    uworst=u0/d
    lworst=l0*exp((18+.05)*TMAX/LOOP)
    weak_floor=1/(1/(G0[1]-GE[1])**2+SM[1]*TMAX)
    assert 0<d and uworst<3 and lworst<.01 and .02 < 4.5*weak_floor
    rlo=(Y0[1,1]-YE[1,1])/(Y0[1,0]+YE[1,0])*d**(1/6)
    rhi=(Y0[1,1]+YE[1,1])/(Y0[1,0]-YE[1,0])/d**(1/6)
    elo=(Y0[2,1]-YE[2,1])/(Y0[2,0]+YE[2,0])
    ehi=(Y0[2,1]+YE[2,1])/(Y0[2,0]-YE[2,0])*exp(1.5*lworst*TMAX/LOOP)
    return {"down":[rlo,rhi],"lepton":[elo,ehi],"min_required_weight":
            [wmin(rhi,elo),wmin(rlo,ehi)],"comparison_D":d}


def wmin(r,e):
    return (sqrt(e)-1/sqrt(e))/(sqrt(r)-1/sqrt(r))


def plane_rotation(i,j,sine,phase=0.):
    rot=np.eye(3,dtype=complex)
    rot[i,i]=rot[j,j]=sqrt(1-sine*sine)
    rot[i,j]=sine*np.exp(-1j*phase)
    rot[j,i]=-rot[i,j].conjugate()
    return rot


def central_ratios(tstar,mass,nstep):
    split=log(mass/MZ)
    v=plane_rotation(1,2,.04193) @ plane_rotation(0,2,.00370,1.139) @ plane_rotation(0,1,.2251)
    mats=np.array([np.diag(Y0[0]),v@np.diag(Y0[1]),np.diag(Y0[2])])
    def derivative(t,m):
        g1,g2,g3=1/(1/G0**2+SM*t+CHANGE*max(0.,t-split))
        hu,hd,he=m@m.conj().transpose(0,2,1)
        common=np.trace(3*(hu+hd)+he).real
        scalar=common-np.array([17*g1/20+9*g2/4+8*g3,g1/4+9*g2/4+8*g3,9*(g1+g2)/4])
        return ((1.5*np.array([hu-hd,hd-hu,he])+scalar[:,None,None]*np.eye(3))@m)/LOOP
    # Integrate to the actual root and break the steps at the threshold.
    for a,b in zip([0,min(tstar,split)],[min(tstar,split),tstar]):
        if a==b:
            continue
        h=(b-a)/nstep
        for step in range(nstep):
            t=a+step*h
            k1=derivative(t,mats)
            k2=derivative(t+h/2,mats+h*k1/2)
            k3=derivative(t+h/2,mats+h*k2/2)
            k4=derivative(t+h,mats+h*k3)
            mats+=h*(k1+2*k2+2*k3+k4)/6
    _,dv,lv=np.sort(np.linalg.svd(mats,compute_uv=False),axis=1)
    return dv[1]/dv[0],lv[1]/lv[0]


def profile_at(mass):
    _,t,w=crossings(mass,G0)[-1]
    r,e=central_ratios(t,mass,64)
    return w-wmin(r,e)


def compute():
    box=flav_box()
    cases={}
    for mass in [1e3,1e12,2e13,1e14]:
        center=crossings(mass,G0)
        corners=[crossings(mass,G0+GE*np.array(s)) for s in product([-1,1],repeat=3)]
        assert all([row[0] for row in k]==[row[0] for row in center] for k in corners)
        found=[]
        for j,(branch,t,w) in enumerate(center):
            times,ws=np.array([k[j][1:] for k in corners]).T
            r,e=central_ratios(t,mass,128)
            rr,ee=central_ratios(t,mass,64)
            assert np.allclose([r,e],[rr,ee],rtol=1e-10)
            found.append({"branch":branch,"central_scale_GeV":MZ*exp(t),
                          "log_scale_box":[min(times),max(times)],
                          "central_weight":w,"weight_box":[min(ws),max(ws)],
                          "central_down_lepton_ratio":[r,e],"central_weight_needed":wmin(r,e),
                          "certified_flavour_existence_in_full_box":bool(min(ws)>box["min_required_weight"][1]),
                          "certified_flavour_exclusion_in_full_box":bool(max(ws)<box["min_required_weight"][0])})
        cases[format(mass,"g")]=found
    assert not cases["1000"]
    assert cases["1e+12"][-1]["certified_flavour_exclusion_in_full_box"]
    assert cases["2e+13"][-1]["certified_flavour_existence_in_full_box"]
    assert len(cases["1e+14"])==1
    # Retrospective 1372 check with no new negative contribution.
    radial=Path(__file__).resolve().parents[1]/"1372"/"results.json"
    required=json.loads(radial.read_text(encoding="utf-8-sig"))["bosonic_fourth_root_with_f_equal_v_no_neutral"]
    septet_mass=required["center"]["Higgs_matching_boundary_GeV"]/14**0.25
    assert not crossings(septet_mass,G0)
    assert all(not crossings(required["box_high"]["Higgs_matching_boundary_GeV"]/14**0.25,
                             G0+GE*np.array(sign)) for sign in product([-1,1],repeat=3))
    # Numerical location of one flavour-boundary solution. This is inferred
    # from the target data, and no uniqueness of the coupled non-linear flow is asserted.
    a,b=12.,13.
    assert profile_at(10**a)<0<profile_at(10**b)
    for _ in range(22):
        mid=(a+b)/2
        if profile_at(10**mid)<0:
            a=mid
        else:
            b=mid
    return {"former_inputs":"https://arxiv.org/html/2510.01312v2","crossings":cases,
            "whole_window_flavour_bounds":box,"fitted_threshold_bracket_GeV":[10**a,10**b],
            "retrospective_full_radial_source_septet_mass_GeV":septet_mass,
            "beta_increments_for_g1_g2_g3":(-CHANGE*LOOP/2).tolist()}


def rounded(x):
    if isinstance(x,dict):
        return {k:rounded(v) for k,v in x.items()}
    if isinstance(x,list):
        return [rounded(i) for i in x]
    return float(format(x,".9g")) if isinstance(x,float) else x


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true",help="verify frozen results (default)")
    p.add_argument("--write",action="store_true")
    op=p.parse_args()
    file=Path(__file__).with_name("independent_results.json")
    answer=rounded(compute())
    if op.write:
        with file.open("x",encoding="utf-8") as out:
            json.dump(answer,out,indent=2,ensure_ascii=False)
            out.write("\n")
    else:
        assert answer==json.loads(file.read_text(encoding="utf-8"))
    print("1376 independent intersections, flavour transport and full input-box tests passed.")
