"""Match only the measured candidate's additional dyad to a hypercharged scalar triplet."""
import argparse
import json
import math
from pathlib import Path
import runpy

import numpy as np


HERE = Path(__file__).resolve().parent
MASS = 2000.
VEV = 246.22
ALPHA = 1/137.036
GF = 1/(math.sqrt(2)*VEV**2)
MAX3 = 1e-12
MAXG = 1.5e-13


def read(path):
    return json.loads((HERE/path).read_bytes())


def module(path):
    return runpy.run_path(str(HERE/path))


def below_both_clusters(flow, rep, pole):
    # Keep v=246.22 GeV as the common coefficient normalization. The SM
    # Weinberg scalar and flavour coefficients follow the frozen 1706 convention.
    y=np.zeros((4,3,4),complex)
    y[0,:,:3]=np.diag(flow["central"]("up"))
    y[1,:,:3]=flow["low_ckm"](rep["gamma_MZ_deg"],flow["angles"]("direct"))@np.diag(flow["central"]("down"))
    y[2,:,:3]=np.diag(flow["central"]("charged_leptons"))
    u=np.zeros(5)
    u[0]=flow["physics"]["quartic"]["lambda_paper"][0]/4
    weights=np.array([3.,3.,1.,1.])
    def rhs(a, v, g):
        la=v[0]
        h=a@a.conj().transpose(0,2,1)
        T=weights@np.trace(h,axis1=1,axis2=2).real
        T2=weights@np.trace(h@h,axis1=1,axis2=2).real
        gy,gw=.6*g[0],g[1]
        dl=24*la**2+(4*T-9*gw-3*gy)*la-2*T2+3/8*(2*gw**2+(gy+gw)**2)
        return flow["beta"](a,g,False),np.r_[dl,4*la+2*T-3*gw,-1.5*np.diag(h[2]).real]
    for start,susy,dt,gs in flow["segments"](.06):
        assert not susy
        for i in range(0,len(gs)-1,2):
            t0=math.log(start)+i/2*dt*flow["LOOP"]
            frac=min(1.,(math.log(pole)-t0)/(dt*flow["LOOP"]))
            if frac<=0:
                return np.exp(u[1]),np.exp(u[2:])
            h=dt*frac
            interp=lambda f:1/((1-f)/gs[i]+f/gs[i+2])
            k1,l1=rhs(y,u,gs[i])
            k2,l2=rhs(y+h*k1/2,u+h*l1/2,interp(frac/2))
            k3,l3=rhs(y+h*k2/2,u+h*l2/2,interp(frac/2))
            k4,l4=rhs(y+h*k3,u+h*l3,interp(frac))
            y+=h*(k1+2*k2+2*k3+k4)/6
            u+=h*(l1+2*l2+2*l3+l4)/6
            if frac<1:
                return np.exp(u[1]),np.exp(u[2:])
    raise ValueError("threshold unavailable")


def heavy_W(model, flow, rep, mw):
    poles=np.array(flow["fields"]["neutral_poles_GeV"])
    me=np.array([(poles[0]+poles[1])/2,(poles[2]+poles[3])/2])
    frame=module("../1706/independent_check.py")["charged_at_threshold"](flow,rep,poles[2])
    d=frame.T@model["at_cluster"](flow,rep,.1).conj()*VEV/math.sqrt(2)
    empty=np.zeros((2,2))
    total=np.block([[np.zeros((3,3)),d,np.zeros((3,2))],
                    [d.T,empty,np.diag(me)],[np.zeros((2,3)),np.diag(me),empty]])
    _,s,vh=np.linalg.svd(total)
    active=vh.conj().T[:3,:4]
    loop=module("../1709/independent_check.py")["finite_loop"]((s[:4]/mw)**2)
    amp=(active[0]*active[1].conj())@loop
    loss=np.sum(abs(active)**2,axis=1)
    return 3*ALPHA/(2*np.pi)*abs(amp)**2/((1-loss[0])*(1-loss[1]))


def main():
    old=module("../1716/check.py")
    null,_,_,_,_,(W,A,arms,_),_=old["background"](.1)
    flow,rep=old["FLOW"]["charged_background"](.1)
    common, leg=below_both_clusters(flow,rep,MASS)
    bl=arms[:,1]
    b=leg*bl
    vnorm=np.linalg.norm(b)
    p=abs(b/vnorm)**2
    pair=np.outer(p,p)
    channels={str(i)+str(j):pair[i,j]*(1+(i!=j)) for i in range(3) for j in range(i,3)}
    assert abs(sum(channels.values())-1)<1e-12
    M_W=read("../1631/inputs.json")["masses_GeV"]["W"]
    hbarc=read("../1631/inputs.json")["hbar_c_GeV_mm"]
    w=(M_W/MASS)**2
    phase_space=math.sqrt(1-4*w)*(1-4*w+12*w*w)
    gamma_over_three=27*ALPHA/(16*np.pi*p[0]**2)
    Wbr=heavy_W(old,flow,rep,M_W)
    scalar_max_gamma=(math.sqrt(MAXG)+math.sqrt(Wbr))**2
    max3=min(MAX3,scalar_max_gamma/gamma_over_three)
    three=lambda y:y**4*p[0]**3*p[1]/(4*GF**2*MASS**4)
    ymax=(max3/three(1.))**.25
    assert max3<MAX3
    results=[]
    for point in read("../1716/results.json")["real_b_port_solutions"]:
        z=np.exp(1j*np.r_[0.,point["two_inverse_row_phases_rad"]])
        C=W.T@np.diag(1/(z*A))@W
        base=C-np.outer(C@null,C@null)/(null@C@null)
        rho=point["relative_real_coefficient_on_high_lifted_b"]
        dyad=rho*np.outer(bl,bl)
        m=np.linalg.svd(base+dyad,compute_uv=False)
        schur_scale=math.sqrt(old["CFG"]["solar_normalization_eV2"])/m[1]
        low_norm=schur_scale*np.linalg.norm(dyad,2)
        pole_eV=low_norm*common*(vnorm/np.linalg.norm(bl))**2
        pole_GeV=pole_eV*1e-9
        A_l=MASS*pole_GeV**2/(16*np.pi)
        A_w=MASS**3*phase_space/(4*np.pi*VEV**4)
        balanced=(A_l/A_w)**.25
        ll=A_l/balanced**2
        ww=A_w*balanced**2
        prod=(MASS**2*pole_GeV/(8*np.pi*VEV**2))**2*phase_space
        np.testing.assert_allclose([ll,ww],math.sqrt(prod),rtol=2e-12)
        # Continuum of L-conserving measurements identical at fixed h for different m_b.
        vdelta_floor_eV=pole_eV/(math.sqrt(2)*ymax)
        results.append({"signed_high_rho":rho, "weak_rank_one_mass_eV":low_norm,
                        "base_weak_operator_norm_eV":schur_scale*np.linalg.norm(base,2),
                        "pole_rank_one_mass_eV":pole_eV,"branch_pattern":channels,
                        "rank_one_mbb_magnitude_eV":low_norm*abs(bl[0]/np.linalg.norm(bl))**2,
                        "total_mbb_eV":point["full_forward"]["mbb_eV"],
                        "Gamma_ll_times_Gamma_WW_GeV2":prod,"VEV_eV_at_equal_widths":balanced*1e9,
                        "max_proper_decay_length_micron":hbarc*1000/(ll+ww),
                        "necessary_VEV_eV_with_two_LFV_bounds":vdelta_floor_eV})
    return {"triplet_mass_example_GeV":MASS,"rank_one_flavor_weights_at_pole":p.tolist(),
            "flavor_universal_RGE_factor":common,"one_leg_RGE_factors":leg.tolist(),
            "degenerate_triplet_Br_gamma_over_Br_3e":gamma_over_three,
            "existing_Dirac_cluster_W_only_mu_e_gamma":Wbr,
            "scalar_only_Br_3e_max_from_MEGII":MAXG/gamma_over_three,
            "scalar_Br_3e_max_allowing_any_Dirac_W_phase":max3,
            "maximum_rank_one_Yukawa_norm_with_two_LFV_bounds":ymax,
            "two_1716_solutions":results}


def tidy(x):
    if isinstance(x,dict): return {k:tidy(v) for k,v in x.items()}
    if isinstance(x,list): return [tidy(v) for v in x]
    if isinstance(x,(float,np.floating)): return float(f"{x:.11g}")
    return x


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--check",action="store_true")
    value=tidy(main())
    target=HERE/"independent_results.json"
    if ap.parse_args().check:
        assert json.loads(target.read_bytes())==value
        print("1718 verified the rank-one tagged triplet, its two-width invariant and charged limits")
    else:
        target.write_text(json.dumps(value,indent=2)+"\n",encoding="utf8")
