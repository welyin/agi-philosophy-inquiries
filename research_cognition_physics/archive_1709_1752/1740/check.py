"""Unordered ray selection, spurion competitors and tagged two-source pole readout."""
import argparse
import json
import math
from pathlib import Path
import runpy

import numpy as np

HERE=Path(__file__).resolve().parent
def read(name):
    return json.loads((HERE/name).read_bytes())
CFG=read("inputs.json")
PARENT=runpy.run_path(str(HERE/CFG["parent_model"]))
I=np.eye(3,dtype=complex)


def outer(x):
    return np.outer(x,x.conj())


def threshold(step):
    args=PARENT["CFG"]["nearby_measured_inputs_in_marginal_widths"][CFG["fit"]]
    _,pt,nu=PARENT["charged_endpoint"](2,step,args["weak"],args["up"],args["g3"],args["top"],args["charm"])
    flow=pt["flow"]
    y=PARENT["initial_state"](pt)
    poles=flow["fields"]["neutral_poles_GeV"]
    seg=[entry for entry in flow["segments"](step) if entry[0]<poles[2]]
    total=I.copy()
    for start,susy,dt,g in seg:
        assert not susy
        for j,m in enumerate(poles):
            if start==m:
                y[3,:,j]=flow["physics"]["Y0"][:,j]
        y,total=PARENT["advance"](y,total,dt,g,flow)
    return flow,y,seg,(pt,nu,total)


def hadron(fb,flow,y,seg,transport=None,fd=None):
    # Only the upper/right mass row changes. Integrate the daughter's actual matrices.
    U,s,Vh=np.linalg.svd(y[1,:,:3])
    u=flow["basis"](y[0])[0][:,0]
    c=abs(u.conj()@U[:,-1])
    if fd is None:
        fd=c
    lows=y.copy()
    lows[1,:,:3]=(U*(s*np.array([fb,1,fd])))@Vh
    sb=I.copy()
    for start,_,dt,g in reversed(seg):
        lows,sb=PARENT["advance"](lows,sb,-dt,g[::-1],flow)
        for j,m in enumerate(flow["fields"]["neutral_poles_GeV"]):
            if start==m:
                lows[3,:,j]=0.
    q,_,gamma,L,up=flow["extraction"](lows)
    mods=dict(zip(("Vus","Vub","Vcb"),np.sin(q[3:])*[np.cos(q[4]),1,np.cos(q[4])]))
    ref=PARENT["DIRECT"]
    offsets={key[1:]:(v-ref[key][0])/math.hypot(*ref[key][1:]) for key,v in mods.items()}
    m=PARENT["MATCH"]
    ks,_=m["thresholds"](dict(weak_offsets=offsets,moduli_MZ=mods,gamma=gamma),m["source_thresholds"]())
    shape=np.array([q[2]/q[1]/ks[2],1.,q[0]/q[1]/ks[0],up[0]/q[1]/(ks[1]*ks[0])])
    profile=PARENT["project_shared"](shape)
    absolute=PARENT["absolute_table"](shape,up[0]/up[1],profile)
    result=dict(MZ_bottom=q[2], bottom_marginal_width=(q[2]-.0163)/.00009, four_absolute_radius=absolute["radius_at_Kuc_one"],
                four_shape_radius=profile["radius"],low_bs=shape[0],gamma_deg=gamma)
    if transport:
        pt,nu,sv=transport
        E=pt["low_E"]
        phases=np.diag(L.conj().T@E)
        L=L@np.diag(phases/abs(phases))
        H=E.conj().T@np.linalg.solve(sv,np.linalg.solve(sb,L))
        neut=nu["observables"](H.T@pt["C"]@H,H.T@pt["b"],nu["FROZEN"][0])
        cp=PARENT["CP"]["test"](mods,gamma)
        kaon=mods["Vus"]/math.sqrt(1-mods["Vus"]**2-mods["Vub"]**2)
        result.update(CP_B_and_epsilon=[cp["sin2beta"],cp["epsilon"]], quark_moduli_widths=offsets,
                      RM123_and_RBC_kaon_widths=[(kaon-center)/err for tag,(center,err) in PARENT["JOINT"]["KAON"].items()
                                                   if tag in ("RM123_axial","RBC_axial")],
                      neutrino_angles=neut["angle_readings"],neutrino_angles_two_widths=neut["angles_in_box"],
                      frequencies_eV2=[neut["dm31_eV2"],neut["dm32_eV2"]],frequencies_one_width=neut["dm31_and_32_in_box_at_central_solar"],
                      FC68_90=neut["SK_T2K_FC68_90"])
    return result


def pole_analysis(r,u,d):
    a,b=r,1.
    c=abs(u.conj()@d)
    Pu,Pd=outer(u),outer(d)
    K=a*Pu+b*Pd
    eig,U=np.linalg.eigh(K)
    plus=U[:,-1]
    pu,pd=abs(u.conj()@plus)**2,abs(d.conj()@plus)**2
    gap=math.sqrt((a-b)**2+4*a*b*c*c)
    D=(a-b)*(1-c*c)/gap
    assert abs(pu-pd-D)<1e-13
    assert abs(pu-.5*(1+(a-b+2*b*c*c)/gap))<1e-13
    spectrum=CFG["pole_mass_common_squared_units"]+eig
    readings=[]
    for q2 in CFG["Euclidean_q2_in_same_units"]:
        resolvent=np.linalg.inv((q2+CFG["pole_mass_common_squared_units"])*I+K)
        measured=(np.vdot(u,resolvent@u)-np.vdot(d,resolvent@d)).real
        expect=-(a-b)*(1-c*c)/((q2+spectrum[1])*(q2+spectrum[2]))
        assert abs(measured-expect)<2e-14
        readings.append(measured)
    return dict(r=r,common_pole_prob_u=pu,common_pole_prob_d=pd,fb=math.sqrt(pu),
                difference=D,pole_squared_masses_in_units=spectrum,normalized_Gu_minus_Gd=readings)


def ratio_from_difference(D,c):
    k=1-c*c
    t=D*c/math.sqrt(k*(k-D*D))
    return (1+t)/(1-t)


def compute(step, detail=True):
    flow,y,seg,transport=threshold(step)
    u,yu=flow["basis"](y[0])
    d,yd=flow["basis"](y[1])
    V=u.conj().T@d
    J=(V[0,0]*V[1,1]*np.conj(V[0,1]*V[1,0])).imag
    u,d=u[:,0],d[:,0]
    c=abs(np.vdot(u,d))
    d=d*np.exp(-1j*np.angle(np.vdot(u,d)))
    A,B=outer(u),outer(d)
    plus=(u+d)/np.linalg.norm(u+d)
    minus=(u-d)/np.linalg.norm(u-d)
    normal=I-outer(plus)-outer(minus)
    swap=outer(plus)-outer(minus)+normal
    assert np.linalg.norm(swap@A@swap.conj().T-B)<1e-12
    assert np.linalg.norm(A@normal)+np.linalg.norm(B@normal)<1e-12
    x,ytr=yu@yu,yd@yd
    cases={}
    for l in CFG["real_symmetric_trace_multipliers"]:
        a,b=1+l*x,1+l*ytr
        assert min(a,b)>0
        obs=pole_analysis(a/b,u,d)
        obs["hadron"]=hadron(obs["fb"],flow,y,seg,transport)
        cases[str(l)]=obs
    resolutions=CFG["example_absolute_error_on_normalized_pole_weight_difference"]
    window=[ratio_from_difference(D,c) for D in (-resolutions,resolutions)]
    assert max(abs(pole_analysis(r,u,d)["difference"]-D) for r,D in zip(window,(-resolutions,resolutions)))<1e-13
    def radius(r):
        return hadron(pole_analysis(r,u,d)["fb"],flow,y,seg)["four_absolute_radius"]
    crossings=[]
    if detail:
        for lo,hi in ((.05,1.),(1.,10.)):
            sign=radius(lo)>2
            assert (radius(hi)>2)!=sign
            for _ in range(37):
                mid=(lo+hi)/2
                if (radius(mid)>2)==sign: lo=mid
                else: hi=mid
            # Round outwards and recheck signs so the printed endpoint pair remains useful.
            edges=[math.floor(lo*1e8)/1e8, math.ceil(hi*1e8)/1e8]
            residual=[radius(q)-2 for q in edges]
            assert residual[0]*residual[1]<0
            crossings.append(dict(r_bracket=edges,radius_minus_two_at_ends=residual,
                                  pole_difference=pole_analysis((lo+hi)/2,u,d)["difference"]))
    clone_bound={}
    for ctest in (.3,.1,.01):
        required=math.ceil(math.log(ctest)/math.log1p(-2*ctest**2))
        def distance(n): return math.sqrt(-math.expm1(2*n*math.log1p(-2*ctest**2)))
        assert distance(required-1)<math.sqrt(1-ctest*ctest)<=distance(required)
        clone_bound[str(ctest)]=required
    return dict(c=c,singular_up=yu,singular_down=yd,traces=[x,ytr],min_abs_CKM=np.min(abs(V)),
                J_threshold=J,Z12_unbroken_powers=[i for i in range(12) if i*read(CFG["flavor_sources"])["flavor_charge"]["S"]%12==0],
                examples=cases,ratio_for_pole_difference_bounded_by=window,absolute_table_two_crossings=crossings,
                independent_source_derivative_at_equal=(1-c*c)/(2*c),minimum_copies_necessary_for_unknown_input_test=clone_bound,
                exchange_b_and_d_physical_routes=hadron(c,flow,y,seg,transport,fd=math.sqrt((1+c)/2)),
                switch_to_mixed_output_purity=(1+c*c)/2,
                unequal_fixed_ratio_examples={str(r):pole_analysis(r,u,d) for r in CFG["fixed_ratio_competitors"]})


def run():
    result=compute(.25)
    fine=compute(.125,False)
    key="four_absolute_radius"
    err=max(abs(v["hadron"][key]-fine["examples"][k]["hadron"][key]) for k,v in result["examples"].items())
    assert err<1e-5 and abs(result["examples"]["0"]["hadron"][key]-1.0175314092)<5e-8
    assert result["min_abs_CKM"]>.001 and abs(result["J_threshold"])>1e-5
    comp=result["examples"]["-0.2"]["hadron"]
    assert result["exchange_b_and_d_physical_routes"]["four_absolute_radius"]>2
    assert comp["four_absolute_radius"]<2 and abs(comp["bottom_marginal_width"])<2
    assert all(comp["neutrino_angles_two_widths"]) and all(comp["frequencies_one_width"]) and comp["FC68_90"]==[False,True]
    result["half_step_max_absolute_radius_change"]=err
    def tidy(value):
        if isinstance(value,dict): return {k:tidy(v) for k,v in value.items()}
        if isinstance(value,(np.ndarray,tuple,list)): return [tidy(v) for v in value]
        if isinstance(value,(float,np.floating)): return float(f"{value:.11g}")
        return value
    return tidy(result)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    answer=run()
    if args.check:
        assert answer==read("results.json")
        print("1740 verified ray exchange, trace competitors, actual daughter and pole response")
    else:
        (HERE/"results.json").write_text(json.dumps(answer,indent=2,allow_nan=False)+"\n",encoding="utf8")

