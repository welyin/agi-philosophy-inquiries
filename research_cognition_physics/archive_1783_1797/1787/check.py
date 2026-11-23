"""Quantify loss of the labelled mass zero and construct its two-insertion repair."""
import argparse
import hashlib
import json
import math
import runpy
from pathlib import Path

import numpy as np


HERE=Path(__file__).resolve().parent
CFG=json.loads((HERE/"inputs.json").read_bytes())


def read(key):
    return json.loads((HERE/CFG[key]).read_bytes())


def floors():
    obs=read("oscillation_windows")
    rect=obs["experimental_rectangles"]
    old=read("charged_and_neutral_point")
    d=old["class_bounds"]["first_axis_tilt_sine"]+1.e-10
    c=math.sqrt(1-d*d)
    weight_min,weight_max=rect["isolated_weight_roots"][0]
    share_min,share_max=rect["solar_share"]
    dlo,dhi=rect["solar_frequency_eV2"]

    def tilt(amplitude,norm2):
        return max(0.,c*amplitude-d*math.sqrt(max(0.,norm2-amplitude**2)))

    atm=obs["assignments"]["NO"]["absolute_atmospheric_frequency_eV2"]
    m2hi,m3hi=math.sqrt(dhi),math.sqrt(atm[1])
    m2lo,m3lo=math.sqrt(dlo),math.sqrt(atm[0])
    pmax=m2hi*(1-weight_max)*share_max+m3hi*weight_max
    pmin=m2lo*(1-weight_min)*share_min+m3lo*weight_min
    bounds={}
    for val in CFG["possible_future_light_exchange_meV"]:
        target=val/1000
        if target>pmax:
            raise ValueError("target is outside all rank-two oscillation predictions")
        # Fixed-target smaller root decreases with P; larger root increases with P.
        small=target/math.sqrt(pmax+math.sqrt(pmax*pmax-target*target))
        p=max(pmin,target)
        large=math.sqrt(p+math.sqrt(p*p-target*target))
        bounds[str(val)]=1000*tilt(small,m2hi+m3hi)*tilt(large,m2hi+m3hi)
    atm=obs["assignments"]["IO"]["absolute_atmospheric_frequency_eV2"]
    # Archive's IO gap |m3^2-m1^2| with m3=0 gives first mass, then add d21.
    # Use conservative extremes valid also if the quoted fast gap is assigned to m2.
    low=math.sqrt(atm[0]-dhi)
    high=math.sqrt(atm[1]+dhi)
    dominant=math.sqrt(low*(1-weight_max)*(1-share_max))
    sub=math.sqrt(high*(1-weight_min)*share_max)
    small=max(dominant-sub,0.)
    large=math.sqrt(low*(1-weight_max))
    io=1000*tilt(small,2*high)*tilt(large,2*high)
    return {"hypothetical_NO_lower_meV":bounds,"IO_lower_either_close_pair_assignment_meV":io}


def phase_minimum(C,b,target):
    # Rephase just the physical low-energy Takagi contribution of mass m2.
    _,sing,vh=np.linalg.svd(C)
    axis=vh.conj().T[:,1]
    mode=(axis.T@C@axis)*np.outer(axis.conj(),axis.conj())
    rest=C-mode
    a,h=mode[0,0],rest[0,0]
    cosine=(target**2-abs(a)**2-abs(h)**2)/(2*abs(a*h))
    assert -1<=cosine<=1
    etas=[t*math.acos(cosine)-np.angle(a/h) for t in [-1,1]]
    values=[]
    for angle in etas:
        new=np.exp(1j*angle)*mode+rest
        assert np.max(abs(new.conj().T@new-C.conj().T@C))<1e-14
        assert abs(abs(new[0,0])-target)<1e-12
        values.append(abs(b.T@new@b))
    return min(values)*1000


def weak_values(C):
    _,sv,vh=np.linalg.svd(C)
    u=vh.conj().T[:,::-1]
    x=abs(u[0,2])**2
    return {"solar_mass2":sv[1]**2,"large_mass2":sv[0]**2,"s12":abs(u[0,1])**2/(1-x),
            "s23":abs(u[1,2])**2/(1-x),"A13":4*x*(1-x),"light_exchange_meV":1000*abs(C[0,0])}


def concrete(step):
    core=runpy.run_path(str(HERE/CFG["original_material"]))
    source=runpy.run_path(str(HERE/CFG["charged_transfer"]))
    flow,rep=source["charged_background"](step)
    W,_=source["congruence"](flow,rep,step)
    sample=core["witness"](W)
    pars=sample["high_angles_ratio"]
    m,H,R,_,_=core["observable"](W,pars)[2]
    r=pars[2]
    scale=sample["common_mass_scale_eV"]
    C0=scale*W.T@R@W
    axis=np.linalg.solve(W,np.eye(3)[:,0])
    norm=np.linalg.norm(axis)
    b=axis/norm
    assert abs(b.T@C0@b)<1e-15
    exact={str(v):phase_minimum(C0,b,v/1000) for v in CFG["possible_future_light_exchange_meV"]}
    # Exactly realize a real high eigenline with variable complex Majorana phase.
    a=(math.sqrt(r)*m+H)/math.sqrt(2)
    v=(math.sqrt(r)*m-H)/math.sqrt(2)
    D=np.column_stack([(a-v)/math.sqrt(2),math.sqrt(2)*v])
    assert abs(D[0,1])<1e-14
    xi=np.exp(2j*(-math.pi/6))
    phases=np.diag([1,xi])
    bmev=CFG["heavy_scale_PeV"]
    osc=read("oscillation_windows")["experimental_rectangles"]
    atm=read("oscillation_windows")["assignments"]["NO"]["absolute_atmospheric_frequency_eV2"]
    # Find an actual two-insertion point for each specified experimental target.
    e_mid=scale*r*(W.T@m)[0]**2
    e_high=-scale*(W.T@H)[0]**2
    angles={"single_insertion":0.,"positive_terminal_masses":math.pi/3,
            "opposite_terminal_parities":-2*math.pi/3,"quadrature":math.pi/2,"flipped":math.pi}
    for val in CFG["possible_future_light_exchange_meV"]:
        cosine=((val/1000)**2-abs(e_mid)**2-abs(e_high)**2)/(2*abs(e_mid*e_high))
        angles[f"target_{val}_meV"]=math.acos(cosine)-np.angle(e_mid/e_high)
    sweep=[]
    for label,angle in angles.items():
        z=np.exp(-1j*angle)-1
        M=np.array([[-1,1],[1,z]])
        high=D@np.linalg.solve(M,D.T)
        np.testing.assert_allclose(high,r*np.exp(1j*angle)*np.outer(m,m)-np.outer(H,H),atol=2e-15)
        C=scale*W.T@high@W
        out=weak_values(C)
        assert osc["solar_frequency_eV2"][0]<=out["solar_mass2"]<=osc["solar_frequency_eV2"][1]
        assert atm[0]<=out["large_mass2"]<=atm[1]
        assert osc["solar_share"][0]<=out["s12"]<=osc["solar_share"][1] and .49<=out["s23"]<=.63
        assert osc["isolated_fast_amplitude"][0]<=out["A13"]<=osc["isolated_fast_amplitude"][1]
        poles=np.linalg.svd(M,compute_uv=False)
        # Match raw M22 with real coefficients multiplying X^2 and A2^2, no new phase.
        mat=np.array([[.25,(xi*xi).real/16],[0,(xi*xi).imag/16]])
        coef=np.linalg.solve(mat,[np.real(xi*xi*z),np.imag(xi*xi*z)])
        assert abs(coef[0]/4+coef[1]*xi**2/16-xi**2*z)<1e-14
        assert np.linalg.norm((D@phases)@np.linalg.solve(phases@M@phases,(D@phases).T)-high)<1e-14
        resid=abs(b.T@C@b)
        assert abs(resid*norm**2-scale*abs(D[0,0])**2*abs(z)/abs(np.linalg.det(M)))<1e-14
        if "terminal" in label:
            sign=1 if label=="positive_terminal_masses" else -1
            links=np.array([[1,0],[-xi,1]],dtype=complex)
            np.testing.assert_allclose(-links@np.diag([1,sign])@links.T,phases@M@phases,atol=1e-14)
        if label.startswith("target_"):
            assert abs(out["light_exchange_meV"]-float(label.split("_")[1]))<1e-10
        sweep.append(dict(scenario=label, high_phase_deg=math.degrees(angle),observables=out,residual_meV=1000*resid,
                         normalized_heavy_poles=poles.tolist(),poles_PeV=(poles*bmev).tolist(),
                         registered_real_X2_A22_coeff=coef.tolist(),abs_M22_over_Mstar=abs(z)))
    norm_meV=1000*scale*abs(D[0,0])**2/norm**2
    return dict(charged_gamma_deg=rep["gamma_MZ_deg"],inverse_axis_norm=norm,source_row_factor_meV=norm_meV,
                maximum_residual_if_each_two_insert_coefficient_at_most_one_meV=norm_meV*(3/8)/(1-3/8),
                at_fixed_full_oscillation_tensor_exact_minimum_meV=exact,high_phase_witnesses=sweep)


def tidy(obj):
    if isinstance(obj,dict): return {k:tidy(v) for k,v in obj.items()}
    if isinstance(obj,(list,tuple)): return [tidy(v) for v in obj]
    if isinstance(obj,(float,np.floating)): return float(f"{obj:.10g}")
    return obj


if __name__=="__main__":
    arg=argparse.ArgumentParser()
    arg.add_argument("--check",action="store_true")
    check=arg.parse_args().check
    coarse,fine=[concrete(x) for x in CFG["half_step"]]
    drift=max(abs(p["observables"]["light_exchange_meV"]-q["observables"]["light_exchange_meV"])
              for p,q in zip(coarse["high_phase_witnesses"],fine["high_phase_witnesses"]))
    assert drift<1e-7
    record=tidy({"round":1787,"data_only_meV_bounds":floors(),"two_insertions":fine,"step_change_meV":drift})
    payload=json.dumps(record,ensure_ascii=False,indent=2)+"\n"
    dest=HERE/"results.json"
    if check:
        assert json.loads(dest.read_bytes())==record
    else:
        dest.write_text(payload,encoding="utf-8")
    print("1787 transported zero and physical heavy self-energy",hashlib.sha256(payload.encode()).hexdigest()[:12])
