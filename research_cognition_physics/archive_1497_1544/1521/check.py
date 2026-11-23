"""UV flavor witnesses with the same three traces, triplet thresholds and protected mass data."""
import argparse
import importlib.util
import json
import math
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent


def load(path):
    return json.loads((HERE/path).read_text(encoding="utf8"))


def family(cfg):
    base=(HERE/cfg["flavor_source"]).resolve()
    s=importlib.util.spec_from_file_location("phase1516",base/"check.py")
    mod=importlib.util.module_from_spec(s)
    s.loader.exec_module(mod)
    core=mod.parent_core(mod.read("inputs.json"))
    poles,d,e,hdabs,_,_=core
    degrees=cfg["selected_phase_deg"]
    v,roots,_=mod.mixing_core(core,degrees)
    selected=mod.read("results.json")["solutions"]["direct"][str(degrees)]["one_passing_frame"]
    O=np.array(selected["O"])
    phase=np.exp(1j*math.radians(degrees))
    bare=hdabs[1:]/abs(1+2*phase*roots)
    hd=np.r_[hdabs[0],bare*(1+2*phase*roots)]
    he=np.r_[hdabs[0],bare*(1-3*phase*roots)]
    hd,he=hd/np.linalg.norm(hd),he/np.linalg.norm(he)
    P=np.eye(4)-np.outer(hd,hd.conj())
    Q=np.eye(4)-np.outer(he,he.conj())
    Y=np.column_stack([np.zeros(3),np.diag(poles)])
    R=P@Y.T@v/d
    np.testing.assert_allclose(R.conj().T@R,np.eye(3),atol=3e-10)
    return O,Y,R,Q,v


def get_rates(cases,O,R,Q,v,mass,nuc):
    data=json.loads((HERE/"../1520"/nuc["previous_nucleons"]).read_text(encoding="utf8"))
    mpi=nuc["pion_mass_GeV"]
    def exposure(m,t):
        p=data["mass_GeV"]["p"]
        return p/(32*math.pi)*(1-(m/p)**2)**2*t*nuc["year_seconds"]/nuc["hbar_GeV_s"]
    e_pi=exposure(mpi,nuc["pion_lifetime_years"])
    e_K=exposure(data["mass_GeV"]["K+"],data["proton"]["partial_lifetime_lower_yr"])
    pion,ka,kb=np.zeros((3,4),dtype=complex)
    for D,m in zip(cases,mass):
        diq=O[0]@D@R
        lep=v.conj().T@D@Q
        # Final missing flavors represented in the original four-dimensional
        # parent space with the actual light charged-lepton projector.
        pion+=diq[0].conj()*lep[0]/m**2
        ka+=diq[1].conj()*lep[0]/m**2
        kb+=diq[0].conj()*lep[1]/m**2
    p=nuc["pion_W0_GeV2"][0]*pion
    a=data["W0_GeV2"]["us_d"][0]*ka
    b=data["W0_GeV2"]["ud_s"][0]*kb
    pi_rate=e_pi*np.vdot(p,p).real
    k_rate=e_K*np.vdot(a+b,a+b).real
    base=e_K*(np.vdot(a,a).real+np.vdot(b,b).real)
    extent=2*e_K*abs(np.vdot(a,b))
    # Exact minimum over two independent lattice boxes and all relative phases.
    def magnitude_box(key):
        val,*errors=data["W0_GeV2"][key]
        return (abs(val)-2*sum(errors),abs(val)+2*sum(errors))
    xs,ys=magnitude_box("us_d"),magnitude_box("ud_s")
    aa=e_K*np.vdot(ka,ka).real
    bb=e_K*np.vdot(kb,kb).real
    cc=e_K*abs(np.vdot(ka,kb))
    trial=[(xx,float(np.clip(cc*xx/bb,*ys))) for xx in xs]
    trial.extend((float(np.clip(cc*yy/aa,*xs)),yy) for yy in ys)
    lowest=min(aa*xx**2+bb*yy**2-2*cc*xx*yy for xx,yy in trial)
    return {"Gamma_pi_over_limit":float(pi_rate),"Gamma_K_over_limit":float(k_rate),
            "K_if_relative_phase_free":[float(base-extent),float(base+extent)],
            "K_two_width_lattice_box_lowest_any_phase":float(lowest),
            "max_common_amplitude_transfer_from_both_at_registered_sign":float(1/math.sqrt(max(pi_rate,k_rate)))}


def compute():
    c=load("inputs.json")
    old=load(c["gauge_and_triplets"])["witness"]
    m=np.array(old["even_triplets_GeV"])
    g=old["gauge_unified"]
    O,Y,R,Q,v=family(c)
    beta=math.atan(c["tan_beta"])
    Dl=Y/math.cos(beta)
    yu=np.array(load(c["masses"])["high_scale"]["up"])[:,0]
    Ul=O.T@np.diag(yu)@O/math.sin(beta)
    tree_budget=3*(6/5)*g*g
    norm2=float(np.linalg.norm(Dl)**2)
    max_c=math.sqrt(tree_budget/norm2-1)
    assert max_c>1
    result={}
    for tag,coeff in [("shared_tensor_spends_full_trace",max_c),("extra_tensor_hidden_from_up_quark",c["c_shared_comparison"])]:
        U=[Ul,-m[1]/(m[0]*coeff)*Ul,np.zeros_like(Ul)]
        fill=0 if tag.startswith("shared") else math.sqrt(tree_budget-(1+coeff**2)*norm2)
        D=[Dl,coeff*Dl,fill*np.outer(O[2],np.array([0,1,0,0]))]
        np.testing.assert_allclose(sum(np.linalg.norm(x)**2 for x in D),tree_budget)
        np.testing.assert_allclose(O[0]@D[2],0,atol=2e-15)
        full=sum(np.outer(u.ravel(),d.ravel())/t for u,d,t in zip(U,D,m))
        defect=np.linalg.norm(full)*m[0]/(np.linalg.norm(Ul)*np.linalg.norm(Dl))
        assert defect<1e-13
        src_u=sum(np.linalg.norm(u)**2/t**2 for u,t in zip(U,m))
        src_d=sum(np.linalg.norm(d)**2/t**2 for d,t in zip(D,m))
        result[tag]={"c_heavy_to_light":coeff,"remaining_tensor_norm":fill,
                     "dimension_five_full_relative_residual":float(defect),
                     "trace_Ku_GeV_minus2":float(src_u),"trace_Kd_GeV_minus2":float(src_d),
                     "scalar_exchange_only_central":get_rates(D,O,R,Q,v,m,load(c["nucleons"]))}
    r=m[0]/max(m[1:])
    # The equal-cost solution uses one heavy channel with equal relative amplitudes.
    ratio=math.sqrt(max(m[1:])/m[0])
    np.testing.assert_allclose(1/m[0]-ratio**2/max(m[1:]),0.,atol=1e-28)
    limits=result["shared_tensor_spends_full_trace"]["scalar_exchange_only_central"]
    def max_visible_fraction(rate):
        cap=(m[1]/m[0])**2*((1+(max_c*m[0]/m[1])**2)/math.sqrt(rate)-1)
        return min(1.,(1+cap)*norm2/tree_budget)
    bounds={key:max_visible_fraction(limits[key]) for key in
            ("Gamma_K_over_limit","K_two_width_lattice_box_lowest_any_phase")}
    gauge_scale=5*old["Sigma_vev_GeV"]
    previous=load("../1520/results.json")["previous_flavor_example_at_central_matching"]
    x_size=math.sqrt(max(previous["kaon_rate_at_ref_arbitrary_relative_phase"]))*(1e15/gauge_scale)**2
    return {"triplet_masses_GeV":m.tolist(),"total_down_Frobenius_squared":tree_budget,
            "C5_zero_max_of_two_elastic_responses_at_least_light_reference_times":float(1+r),
            "maximum_down_norm_squared_fraction_on_shared_tensor":bounds,
            "central_gauge_X_amplitude_over_experimental_limit_if_unit_matching":x_size,
            "minimal_extra_amplitude_to_rescue_full_alignment_central":math.sqrt(limits["Gamma_K_over_limit"])-1,
            "trace_Kd_floor_GeV_minus2":float(tree_budget/max(m)**2),"witnesses":result}


def check(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:
            if k=="dimension_five_full_relative_residual":
                assert max(abs(a[k]),abs(b[k]))<1e-13
            else:
                check(a[k],b[k])
    else:
        np.testing.assert_allclose(a,b,rtol=5e-7,atol=1e-75)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check",action="store_true")
    mode.add_argument("--write-new",action="store_true")
    opt=parser.parse_args()
    result=compute()
    if opt.write_new:
        with (HERE/"results.json").open("x",encoding="utf8") as f:
            json.dump(result,f,indent=2)
            f.write("\n")
    else:
        check(result,load("results.json"))
    print("Verified matched three-triplet exchange, trace budgets and physical baryon flavor channels.")
