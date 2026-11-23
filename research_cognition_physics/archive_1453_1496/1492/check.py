"""A positive heavy threshold, an invariant sum of two effective masses, and a CP-source cap."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


HERE=Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def gap_masses(solar,atm):
    return math.sqrt(solar),math.sqrt(solar+atm)


def source_sum(m2,m3):
    return math.sqrt(m2*m2+m3*m3+10*m2*m3)


def loop(x):
    # log1p is stable in the finite numerical stress tests used here.
    return math.sqrt(x)*(1/(1-x)+1-(1+x)*math.log1p(1/x)),1/(1-x)


def flavour_eps(Y,masses):
    h=Y.conj().T@Y
    rows=[]
    for i,j in [(0,1),(1,0)]:
        f,g=loop((masses[j]/masses[i])**2)
        rows.append(np.imag(Y[:,i].conj()*Y[:,j]*(f*h[i,j]+g*h[j,i]))/(8*math.pi*h[i,i].real))
    return np.array(rows)


def compute():
    cfg=read(HERE/"inputs.json")
    osc=read(HERE/cfg["oscillation"])
    sol=osc["KamLAND_dm21_eV2"]; atm=osc["T2K_dm32_NO_eV2"]
    m2,m3=gap_masses(sol[0],atm[0])
    hi2,hi3=gap_masses(sol[0]+2*sol[-1],atm[0]+2*atm[-1])
    R=m3/m2
    poly=6*(R+1/R)+16
    source_r=8/(poly+math.sqrt(poly**2-1936))
    ma=math.sqrt(m2*m3/(6*source_r)); mb=source_r*ma
    a=np.array([0.,1.,1.]); b=np.array([1.,3.,1.])
    eta=2*math.pi/3
    u=cfg["v246_GeV"]/math.sqrt(2)
    higgs=read(HERE/cfg["higgs_observable"])["atlas"]["mass_GeV"]
    mh=higgs[0]; high=mh+2*higgs[1]
    Q=cfg["heavy_to_light_ratio_minimum"]
    aligned=[ma*2,mb*11]
    assert math.isclose(sum(aligned),source_sum(m2,m3),rel_tol=1e-13)
    thermal_yield=135*cfg["zeta3"]/(4*math.pi**4*cfg["thermal_g_star_min"])
    num,den=cfg["baryon_to_BminusL"]
    response=thermal_yield*(num/den)*cfg["late_entropy_per_photon"]
    urot=read(HERE/cfg["charged_basis"])["textures"]["small_left"]["real_minus_witness"]["rotation_degrees"]
    angle=math.radians(urot)
    charged=np.array([[math.cos(angle),-math.sin(angle),0.],[math.sin(angle),math.cos(angle),0.],[0,0,1.]])
    witnesses={}
    for key,ratios in [("a_lighter",[1,Q]),("b_lighter",[Q,1])]:
        factor=(8*math.pi**2*u**2*mh**2/sum(t*1e-9*z**3 for t,z in zip(aligned,ratios)))**(1/3)
        ms=np.array(ratios)*factor
        Y=charged @ np.column_stack((a*math.sqrt(ma*1e-9*ms[0])/u,
                    np.exp(.5j*eta)*b*math.sqrt(mb*1e-9*ms[1])/u))
        light=-(Y/ms)@Y.T*u**2*1e9
        spectrum=np.linalg.svd(light,compute_uv=False)
        np.testing.assert_allclose(spectrum,[m3,m2,0.],atol=2e-14)
        h=Y.conj().T@Y
        threshold=(h.diagonal().real*ms**2).sum()/(8*math.pi**2)
        assert math.isclose(threshold,mh**2,rel_tol=1e-12)
        eps=flavour_eps(Y,ms)
        witnesses[key]={"masses_a_b_GeV":ms.tolist(),"squared_column_norms":h.diagonal().real.tolist(),
                        "Higgs_quadratic_GeV2":threshold,"sum_abs_flavoured_eps":float(abs(eps).sum()),
                        "optimistic_eta_bound_at_witness":float(abs(eps).sum()*response)}
    gap=cfg["heavy_to_light_ratio_minimum"]
    kernel_H=gap/(gap**2-1)+1/(2*gap)+1/(gap**2-1)
    kernel_L=gap/(gap**2-1)+1+1/(1-1/gap**2)
    rho=4/math.sqrt(22)
    def cap(positive_budget,eff_mass_sum):
        h_sum=(8*math.pi**2*positive_budget*(eff_mass_sum*1e-9)**2/u**4)**(1/3)
        return rho*max(kernel_H,kernel_L)*h_sum/(8*math.pi)*response
    budget=high**2/cfg["positive_Higgs_threshold_transport"][0]
    t_hi=source_sum(hi2,hi3)*cfg["source_transfer_each_upper"]
    eta_cap=cap(budget,t_hi)
    needed=cfg["eta_B_comparison_lower"]/eta_cap
    assert needed>1
    return {"central_m2_m3_eV":[m2,m3],"central_ma_mb_eV":[ma,mb],"central_ta_tb_eV":aligned,
            "sum_t_at_upper_gap_edges_eV":source_sum(hi2,hi3),"phase_overlap":rho,
            "one_relativistic_Majorana_yield":thermal_yield,"eta_per_unit_summed_epsilon":response,
            "hierarchy_loop_kernel_upper_heavy_in_loop":kernel_H,"hierarchy_loop_kernel_upper_light_in_loop":kernel_L,
            "matched_witnesses":witnesses,"eta_central_unity_transport_absolute_bound":cap(mh**2,sum(aligned)),
            "eta_full_work_box_absolute_bound":eta_cap,"minimum_source_enhancement_at_5p9e_minus10":needed,
            "positive_Higgs_budget_factor_needed_at_fixed_other_inputs":needed**3}


def equal(a,b):
    if isinstance(a,dict):
        return a.keys()==b.keys() and all(equal(a[k],b[k]) for k in a)
    return np.allclose(a,b,rtol=1e-10,atol=1e-18)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument("--write",action="store_true");mode.add_argument("--check",action="store_true")
    result=compute(); target=HERE/"results.json"
    if parser.parse_args().write:
        with target.open("x",encoding="utf-8") as f:
            json.dump(result,f,indent=2,allow_nan=False);f.write("\n")
    else:
        assert equal(result,read(target))
    print(json.dumps(result,indent=2))
