"""One-loop certificates and running Takagi thresholds for the 1556 tensors."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent
LOOP=8*math.pi**2


def read(p):
    return json.loads(p.read_text(encoding="utf-8-sig"))


def starting(p,vc):
    h=np.diag(p["H_diag"])
    tri=np.triu_indices(3)
    f=np.zeros((3,3),complex)
    f[tri]=np.array(p["F_real_upper"])+1j*np.array(p["F_imag_upper"])
    f[tri[::-1]]=f[tri]
    r,s=p["r"],complex(p["s_real"],p["s_imag"])
    return (r*(h+s*f),h+f,h-3*f,r*(h-3*s*f),-(vc*r)**2*1e9*f/p["t_eV"])


def double_ratio(d,e):
    ds=np.linalg.svd(d,compute_uv=False)[::-1]
    es=np.linalg.svd(e,compute_uv=False)[::-1]
    return ds[0]*es[1]/(ds[1]*es[0])


def compare_step(k,g,a,length):
    u=g*length/LOOP
    correction=a*k/g*math.expm1(u)
    return k*math.exp(u)/(1+correction),LOOP/a*math.log1p(correction)


def gauge_squares(t,weak,cfg,T):
    # ell=t runs from MX downward. Boundary data refer to MZ.
    return 1/(1/np.square(weak)-np.array(cfg["gauge_beta"])*(T-t)/LOOP)


def normative_bounds(initial,weak,cfg,T,cut):
    kg=[float(np.linalg.norm(m,ord=2)**2) for m in initial[:4]]
    ku,kd,ke,kn=kg
    major=max(kn,ke)
    sums=np.zeros(4)
    peaks=np.array(kg)
    grid=sorted(set(np.r_[np.linspace(0,T,cfg["comparison_cells"]+1),cut].tolist()))
    for start,stop in zip(grid,grid[1:]):
        step=stop-start
        g2=np.maximum(gauge_squares(start,weak,cfg,T),gauge_squares(stop,weak,cfg,T))
        gu,gd,ge,gn=[float(g2@cfg["gauge_terms"][tag])*(1+1e-13) for tag in ("u","d","e","n")]
        ku,iu=compare_step(ku,gu,4.5,step)
        kd,idown=compare_step(kd,gd,4.5,step)
        new_major,_=compare_step(major,max(ge,gn),2,step)
        major_bound=max(major,new_major)
        new_e,ie=compare_step(ke,ge+.5*major_bound,2.5,step)
        if start<cut:
            kn,inu=compare_step(kn,gn+.5*max(ke,new_e),2.5,step)
        else:
            kn,inu=0,0
        ke,major=new_e,new_major
        sums+=np.array([iu,idown,ie,inu])
        peaks=np.maximum(peaks,[ku,kd,ke,kn])
    return sums,peaks


def derive(state,t,weak,cfg,T):
    u,d,e,n,m=state
    ku,kd,ke,kn=[v.conj().T@v for v in (u,d,e,n)]
    trace=float((3*np.trace(ku+kd)+np.trace(ke+kn)).real)
    g=gauge_squares(t,weak,cfg,T)
    def beta(v,k,other,sector):
        gauge=g@cfg["gauge_terms"][sector]
        return ((gauge-trace)*v-1.5*v@(k-other))/(2*LOOP)
    kr=n@n.conj().T
    return [beta(u,ku,kd,"u"),beta(d,kd,ku,"d"),beta(e,ke,kn,"e"),beta(n,kn,ke,"n"),
            -(kr@m+m@kr.T)/(2*LOOP)]


def advance(state,t,dt,weak,cfg,T):
    def shifted(prev,factor):
        return [x+dt*factor*y for x,y in zip(state,prev)]
    a=derive(state,t,weak,cfg,T)
    b=derive(shifted(a,.5),t+dt/2,weak,cfg,T)
    c=derive(shifted(b,.5),t+dt/2,weak,cfg,T)
    d=derive(shifted(c,1),t+dt,weak,cfg,T)
    return [v+dt*(k1+2*k2+2*k3+k4)/6 for v,k1,k2,k3,k4 in zip(state,a,b,c,d)]


def drop_heavy(state):
    _,sing,vh=np.linalg.svd(state[-1])
    rot=vh.conj().T[:,::-1]
    phases=np.angle(np.diag(rot.T@state[-1]@rot))
    rot=rot*np.exp(-.5j*phases)
    if not np.allclose(rot.T@state[-1]@rot,np.diag(sing[::-1]),atol=sing[0]*1e-11,rtol=1e-10):
        raise RuntimeError("Takagi threshold failed")
    return [*state[:3],(rot.T@state[3])[:-1],np.diag(sing[::-1][:-1]).astype(complex)]


def running(initial,weak,cfg,T,steps):
    state=list(initial)
    t=0.
    thresholds=[]
    for target in np.linspace(0,T,steps+1)[1:]:
        while target-t>2e-13:
            trial=advance(state,t,target-t,weak,cfg,T)
            if not len(state[-1]) or np.linalg.norm(trial[-1],2)<cfg["scale_GeV"]*math.exp(-target):
                state,t=trial,target
                break
            lo,hi=0.,target-t
            for _ in range(40):
                dt=(lo+hi)/2
                m=advance(state,t,dt,weak,cfg,T)[-1]
                if np.linalg.norm(m,2)>=cfg["scale_GeV"]*math.exp(-t-dt):
                    hi=dt
                else:
                    lo=dt
            state=drop_heavy(advance(state,t,hi,weak,cfg,T))
            t+=hi
            thresholds.append(cfg["scale_GeV"]*math.exp(-t))
    assert len(thresholds)==3
    return dict(double_ratio_MZ=double_ratio(state[1],state[2]),running_thresholds_desc_GeV=thresholds,
                Yt_MZ=float(np.linalg.norm(state[0],2)),Ytau_MZ=float(np.linalg.norm(state[2],2)))


def evaluate():
    cfg=read(HERE/"inputs.json")
    source=read(HERE/cfg["tensors"])["candidates"]
    gauge=read(HERE/cfg["weak_gauge"])
    T=math.log(cfg["scale_GeV"]/gauge["mu_ref_GeV"])
    centers,errors=np.array([gauge[k] for k in ("g1_GUT","g2","g3")]).T
    widths=cfg["weak_table_marginal_widths"]
    gauge_upper=centers+widths*errors
    weak=read(HERE/cfg["weak_charged"])["yukawa"]
    d,e=[np.array(weak[k]) for k in ("down","charged_leptons")]
    dl,el=[v[:,0]-widths*v[:,1] for v in (d,e)]
    du,eu=[v[:,0]+widths*v[:,1] for v in (d,e)]
    low=dl[0]*el[1]/(du[1]*eu[0])
    ret={"weak_double_ratio_lower":low,"log_interval":T,"branches":{}}
    for name in cfg["candidates"]:
        initial=starting(source[name]["tensors"],cfg["vev_GeV"])
        high=double_ratio(initial[1],initial[2])
        total_need=4*LOOP/3*math.log(low/high)
        all_sums,_=normative_bounds(initial,gauge_upper,cfg,T,T)
        minimum=np.linalg.svd(initial[-1],compute_uv=False)[-1]
        floor=minimum*math.exp(-all_sums[3]/LOOP)
        cut=min(T,max(0.,math.log(cfg["scale_GeV"]/floor)))
        sums,peaks=normative_bounds(initial,gauge_upper,cfg,T,cut)
        total=sums[0]+sums[1]+sums[3]
        gain=math.exp(3*total/(4*LOOP))
        ratio=low/(high*gain)
        finite_component=(math.sqrt(ratio)-1)/(math.sqrt(ratio)+1)
        assert total < total_need and all_sums[0]+all_sums[1]+all_sums[3]<total_need
        numerical=[running(initial,centers,cfg,T,k) for k in cfg["rk4_steps"]]
        coarse,fine=numerical
        assert abs(coarse["double_ratio_MZ"]-fine["double_ratio_MZ"])<8e-8
        assert max(abs(np.array(coarse["running_thresholds_desc_GeV"])/fine["running_thresholds_desc_GeV"]-1))<1e-9
        assert min(fine["running_thresholds_desc_GeV"])>=floor and fine["double_ratio_MZ"]<high*gain
        ret["branches"][name]=dict(D_high=high,needed_integral_min=total_need,
            initial_Ynu_singular_values=np.linalg.svd(initial[3],compute_uv=False).tolist(),
            integral_all_on_u_d_e_n=all_sums.tolist(),integral_with_threshold_u_d_e_n=sums.tolist(),
            arbitrary_CKM_total_area_upper=total,maximum_gain_from_smooth_flow=gain,
            peak_squared_upper_u_d_e_n=peaks.tolist(),minimum_heavy_running_mass_GeV=floor,
            max_active_log_length=cut,low_double_ratio_upper=high*gain,
            minimum_extra_double_ratio_gain=ratio,necessary_equal_relative_four_factor=finite_component,
            numerical_central=fine,coarse_minus_fine_double_ratio=coarse["double_ratio_MZ"]-fine["double_ratio_MZ"])
    return ret


def check_equal(now,expected):
    if isinstance(now,dict):
        assert now.keys()==expected.keys()
        for k,value in now.items():
            check_equal(value,expected[k])
    else:
        np.testing.assert_allclose(now,expected,rtol=5e-9,atol=5e-12)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    flags=parser.add_mutually_exclusive_group(required=True)
    flags.add_argument("--check",action="store_true")
    flags.add_argument("--write-new",action="store_true")
    args=parser.parse_args()
    result=evaluate()
    path=HERE/"results.json"
    if args.check:
        check_equal(result,read(path))
        print("1557: all-scale inequalities and three evolving mass thresholds verified")
    else:
        with path.open("x",encoding="utf-8") as out:
            json.dump(result,out,indent=2)
            out.write("\n")
