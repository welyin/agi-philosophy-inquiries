"""Single-adjoint two-pole matching and a charged-spectrum enclosure on the fixed GJ axes."""
import argparse
import json
import math
from pathlib import Path
import runpy
import numpy as np

HERE=Path(__file__).resolve().parent
I=np.eye(3,dtype=complex)
def read(path):
    return json.loads((HERE/path).read_bytes())
CFG=read("inputs.json")
MODEL=runpy.run_path(str(HERE/CFG["split"]))


def summary(row):
    nu=row["neutrinos"]
    return dict(e_mu_over_tau_residual=row["MZ_relative_e_mu_after_tau_calibration"],
                plane_gain=row["threshold_source_plane_gains"], dm31_eV2=nu["dm31_eV2"],
                dm32_eV2=nu["dm32_eV2"], signed_mother_solar_dm21=nu["dm21_with_mother_solar_gain"],
                angles=nu["angle_readings"], FC68_90=nu["SK_T2K_FC68_90"],
                clocks_at_solar=nu["dm31_and_32_in_box_at_central_solar"],
                absolute_four_quark_radius=row["absolute_four_ratio_radius"], color_check=row["color_scale_c"],
                residual_on_active_columns=row["matching_residual_neutral_columns"])


def common_rows(C,k):
    n=np.sum(C,axis=0)
    n/=np.linalg.norm(n)
    h=C@n
    assert abs(h[0]-h[1])<1e-13
    delta=(C[0]-C[1])/2
    s=np.linalg.norm(delta)
    out={}
    for sign,tag in ((1,"small"),(-1,"reflected")):
        shift=(sign*math.sqrt(k*k-s*s)-h[0])*n
        W=C+shift
        np.testing.assert_allclose(np.linalg.norm(W,axis=1),k)
        np.testing.assert_allclose(W[0]-W[1],C[0]-C[1])
        out[tag]=(W,shift)
    return out


def spectral_box(b,ref,minimum):
    # Ratio f_e/f_tau, f_mu/f_tau refers here to charged singular values and may
    # be recalibrated by a single common scalar. The GJ strange direction is unsuppressed.
    U=b["tagged_light_columns"]
    tagged=U.conj().T@b["yW"][2,:,:3]
    scales=np.linalg.svd(tagged,compute_uv=False)[::-1]
    known=b["model"]["JOINT"]["LEPTONS"]
    widths=2*known[:,1]/known[:,0]
    low_ref=1+ref["MZ_relative_e_mu_after_tau_calibration"]
    mz_box_lo=(1-widths[:2])/(1+widths[2])/low_ref
    mz_box_hi=(1+widths[:2])/(1-widths[2])/low_ref
    flow=b["pt"]["flow"]
    end=b["Mweak"]
    lightest=min(flow["fields"]["neutral_poles_GeV"])
    length=math.log(end/flow["MZ"])/(16*math.pi**2)
    active=math.log(end/lightest)/(16*math.pi**2)
    # Bootstrap physical contractions: downward universal traces are nonpositive,
    # and clipping inactive columns cannot increase their spectral norm.
    cap_e,cap_n=.02,.25
    mx=max(max(np.max(seg[3][:,:2],axis=0)) for seg in b["segments"])
    assert mx<.6
    assert np.linalg.norm(tagged,2)*math.exp((4.5*.6+1.5*cap_n**2)*length)<cap_e
    assert np.linalg.norm(b["yW"][3],2)*math.exp((2.7*.6+1.5*cap_e**2)*active)<cap_n
    # Bound each trajectory's shift in any charged mass ratio by the instantaneous
    # spread of its Hermitian nonscalar anomalous dimension; then compare two trajectories.
    charged_margin=3*(cap_e**2*length+cap_n**2*active)
    lower=mz_box_lo*np.exp(-charged_margin)
    upper=mz_box_hi*np.exp(charged_margin)
    # The unfiltered row forces s_tau'/s_tau >= ||that row||/s_tau.
    # If the two normalized ratios are in their boxes this gives det(F)=fd*fb a floor.
    floor=(np.linalg.norm(tagged[1])/scales[2])**3*np.prod(lower)
    floor=np.nextafter(floor*(1-1e-12),0.)
    lo0=max(minimum,floor)
    cache={}
    def singular(a,bottom):
        key=(a,bottom)
        if key not in cache:
            Y=np.array([a,1.,bottom])[:,None]*tagged
            left,s,rh=np.linalg.svd(Y)
            # Polar orthonormalization of the two approximate singular frames costs
            # at most their Gram errors when those errors are smaller than 1/2.
            errors=np.linalg.norm(left.conj().T@left-I)+np.linalg.norm(rh@rh.conj().T-I)
            res=np.linalg.norm(Y-(left*s)@rh)
            assert res+2*errors*scales[2]<1e-13*scales[2]
            cache[key]=s[::-1]/scales
        return cache[key]
    # Weyl outward allowance > 1000 eps times ||Y|| in absolute mass units,
    # much greater than all checked 3x3 residuals. Integral transport bounds above dominate it.
    pad=1e-11*scales[2]/scales
    queue=[(lo0,1.,lo0,1.)]
    accepted=[]
    rejected=0
    while queue:
        d0,d1,b0,b1=queue.pop()
        bot=np.maximum(0.,singular(d0,b0)-pad)
        top=singular(d1,b1)+pad
        if np.any(bot[:2]>upper*top[2]) or np.any(top[:2]<lower*bot[2]):
            rejected+=1
            continue
        if max(d1-d0,b1-b0)<=CFG["charged_filter_cell_max_side"]:
            accepted.append((d0,d1,b0,b1))
        elif d1-d0>b1-b0:
            half=(d0+d1)/2
            queue.extend(((d0,half,b0,b1),(half,d1,b0,b1)))
        else:
            half=(b0+b1)/2
            queue.extend(((d0,d1,b0,half),(d0,d1,half,b1)))
    assert accepted
    z=np.array(accepted)
    Q=np.linalg.qr(b["N"])[0][:,:2]
    def plane2(d,bb):
        x=U.conj().T@Q
        return np.linalg.eigvalsh(x.conj().T@np.diag([d*d,1.,bb*bb])@x)
    islands={}
    # Mark the ordinary (both factors >0.8) and any far/possibly permuted branch.
    for r in z:
        tag="near_identity" if min(r[0],r[2])>CFG["positive_ancestor_domain_minimum"] else "large_attenuation"
        er=islands.setdefault(tag,dict(cells=0, d=[1.,0.], b=[1.,0.],plane_ratio_upper=1.,plane_ratio_lower=1e100))
        er["cells"]+=1
        for key,idx in (("d",0),("b",2)):
            er[key]=[min(er[key][0],r[idx]),max(er[key][1],r[idx+1])]
        Pmin=plane2(r[0],r[2])
        Pmax=plane2(r[1],r[3])
        er["plane_ratio_upper"]=max(er["plane_ratio_upper"],math.sqrt((Pmax[1]+1e-12)/(Pmin[0]-1e-12)))
        er["plane_ratio_lower"]=min(er["plane_ratio_lower"],math.sqrt((Pmin[1]-1e-12)/(Pmax[0]+1e-12)))
    return dict(det_filter_floor=floor,rejected_cells=rejected,kept_cells=len(accepted),
                charged_log_ratio_RG_allowance=charged_margin,allowed_charged_ratios_at_MW=[lower,upper],possible_components=islands,
                allowed_lepton_norms=[cap_e,cap_n],integration_lengths=[length,active]),accepted


def compute(step,detailed=True):
    b=MODEL["mother"](step)
    run=lambda F:MODEL["daughter"](b,F)
    ref=run(I)
    c=ref["color_scale_c"]
    fc=np.array([c,math.sqrt((1+c)/2)])
    C=np.column_stack((np.sqrt(1-fc*fc),fc))
    k=b["Mweak"]/b["Mcolor"]
    U=b["tagged_light_columns"]
    def weak_route(W):
        lifting=W[:,1]/k
        F=(U*np.array([lifting[0],1.,lifting[1]]))@U.conj().T
        returned=summary(run(F))
        return dict(weak_rows_in_MD=W, signed_weak_gains=lifting,
                    observable_absolute_gains=abs(lifting),observed=returned)
    rows={}
    for name,(W,shift) in common_rows(C,k).items():
        row=weak_route(W)
        row.update(common_shift_in_MD=shift,adjoint_V_mixing_GeV=-shift*b["Mcolor"]/5,
                   singlet_mass_rows_GeV=(C+2*shift/5)*b["Mcolor"])
        rows[name]=row
    models={}
    for tag,W in [("proportional",k*C),("family_bare_only",np.column_stack((C[:,0],np.sqrt(k*k-C[:,0]**2)))),
                  ("family_mix_only",np.column_stack((np.sqrt(k*k-C[:,1]**2),C[:,1])))]:
        models[tag]=weak_route(W)
        models[tag]["adjoint_row_difference_GeV"]=b["Mcolor"]*np.linalg.norm((W-C)[0]-(W-C)[1])/5
    one_slot={}
    for ind,tag in ((0,"mix_only"),(1,"bare_only")):
        a=C[:,ind];avg=np.mean(a);shift=-avg+math.sqrt(avg*avg+k*k-1)
        residual=abs(shift*(a[0]-a[1]))
        one_slot[tag]=dict(best_small_shift_in_MD=shift,minimax_squared_mass_error_MD2=residual,
              two_poles_GeV=b["Mcolor"]*np.sqrt(1+shift*shift+2*shift*a))
    result=dict(c_at_Mc=c,c_at_Mw=b["cW"],color_legs=fc,color_rows=C,mass_ratio=k,
                shared_adjoint=rows,alternative_rows=models,best_one_slot_same_coefficient=one_slot,
                identity=summary(ref))
    if detailed:
        window,accepted=spectral_box(b,ref,0.)
        result["arbitrary_positive_weak_rows_tagged_GJ"]=window
        # Conservative low-energy mass-squared ratio bound using the denominator
        # of the transported mother and the extrema of a congruence restricted to its two-source plane.
        f=b["nu"]["FREQ"]
        sol=f["juno_NO_dm21"]; atm=f["juno_NO_dm31"]
        length,active=window["integration_lengths"]
        logcond=1.5*.02**2*length+.5*.25**2*active
        baseline=ref["neutrinos"]["dm31_eV2"]/sol[0]
        leak=ref["effective_mass_outside_two_Dirac_plane"]
        assert leak<1e-12
        ratio=[((atm[0]-i*atm[1])/((sol[0]+i*sol[2])*baseline*math.exp(8*logcond)*(1+1e-9)))**.25 for i in (1,2)]
        window["required_source_plane_ratio_one_two_widths"]=ratio
        near=window["possible_components"]["near_identity"]
        assert near["plane_ratio_upper"]<min(ratio)
        far=window["possible_components"]["large_attenuation"]
        # Rank-two area from all nine two-by-two minors. In the tagged basis their
        # moduli carry only the four corresponding diagonal attenuation factors.
        pt=b["pt"]
        carrier=np.linalg.inv(b["S"]@pt["low_E"])@U
        mass=carrier.T@(pt["C"]+b["nu"]["FROZEN"][0]*np.outer(pt["b"],pt["b"]))@carrier
        high=np.array([far["d"][1],1.,far["b"][1]])
        pairs=((0,1),(0,2),(1,2))
        # Use explicitly scaled absolute roundoff pads. The expressions are just products
        # and differences; cancellation does not compromise an upper bound on their moduli.
        masspad=1e-12*np.linalg.norm(mass,2)
        total=0.
        for i,j in pairs:
            for u,v in pairs:
                minor=mass[i,u]*mass[j,v]-mass[i,v]*mass[j,u]
                pad=5*masspad*np.linalg.norm(mass,2)
                total+=(high[i]*high[j]*high[u]*high[v]*(abs(minor)+pad))**2
        area=math.sqrt(total)
        ratio_W=(abs(mass[1,1])-masspad)**2/area
        far["all_filter_upper_Takagi_area"]=area
        far["unfiltered_s_s_mass_modulus"]=abs(mass[1,1])
        far["all_filter_lower_m3_over_m2_at_Mw"]=ratio_W
        far["all_filter_lower_m3_over_m2_at_MZ"]=ratio_W*math.exp(-4*logcond)
        far["max_data_m3_over_m2_two_width"]=math.sqrt((atm[0]+2*atm[2])/(sol[0]-2*sol[1]))
        assert far["all_filter_lower_m3_over_m2_at_MZ"]>4000
        assert far["max_data_m3_over_m2_two_width"]<6.
    return result


def tidy(o):
    if isinstance(o,dict): return {k:tidy(v) for k,v in o.items()}
    if isinstance(o,(np.ndarray,list,tuple)): return [tidy(v) for v in o]
    if isinstance(o,(float,np.floating)): return float(f"{o:.12g}")
    if isinstance(o,np.integer): return int(o)
    if isinstance(o,np.bool_): return bool(o)
    return o


def run():
    full=compute(CFG["steps"][1])
    coarse=compute(CFG["steps"][0])
    diffs=[abs(coarse["shared_adjoint"][k]["observed"]["dm31_eV2"]-row["observed"]["dm31_eV2"])
           for k,row in full["shared_adjoint"].items()]
    assert max(diffs)<2e-9
    full["step_refinement_dm31_eV2"]=max(diffs)
    old=coarse["arbitrary_positive_weak_rows_tagged_GJ"]["possible_components"]
    new=full["arbitrary_positive_weak_rows_tagged_GJ"]["possible_components"]
    full["step_refinement_component_gain_upper"]=abs(old["near_identity"]["plane_ratio_upper"]-new["near_identity"]["plane_ratio_upper"])
    full["step_refinement_exotic_mass_ratio_lower"]=abs(old["large_attenuation"]["all_filter_lower_m3_over_m2_at_MZ"]-
                                                        new["large_attenuation"]["all_filter_lower_m3_over_m2_at_MZ"])
    return tidy(full)


if __name__=="__main__":
    args=argparse.ArgumentParser()
    args.add_argument("--check",action="store_true")
    result=run()
    if args.parse_args().check:
        assert result==read("results.json")
        print("1745 verified two equal-pole mass rows, marked light spectra and joint charged/neutral bounds")
    else:
        (HERE/"results.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf8")

