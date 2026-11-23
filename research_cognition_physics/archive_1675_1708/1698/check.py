"""Join the common charged soft thresholds to a four-neutral-pole EFT."""
from pathlib import Path
from itertools import product
import argparse
import json
import math
import runpy

import numpy as np

HERE = Path(__file__).resolve().parent


def read(p):
    return json.loads((HERE/p).read_bytes())


def crossing(fun, low, high):
    assert fun(low)>0 and fun(high)<0
    for _ in range(85):
        mid = math.sqrt(low*high)
        if fun(mid)>0:
            low = mid
        else:
            high = mid
    return math.sqrt(low*high)


def neutral_graph(conf):
    old = read(conf["neutral"])
    reference = read(conf["neutral_reference"])
    specs = read(old["old_model"])
    r = read(old["old_result"])["rank_two_witness"]["small_root_relative_amplitude"]
    parent = read(old["flavor_fields"])
    rho = parent["scalar_magnitude_over_cutoff"]
    s = rho*np.exp(1j*np.pi*parent["scalar_angle_over_pi"])
    x = parent["real_X_over_cutoff"]
    m, p, q = (old[key] for key in ("m_N0H0","coupling_to_s2","coupling_to_x"))
    a,b = np.array(specs["source_a"]),np.array(specs["source_b"])
    rows = (s/rho)**np.array(specs["su5_generation_charges_5bar"])
    one = math.sqrt(r)*b
    two = s**2*(a*x/rho**2+p/m*one)
    unit = old["sterile_unit_GeV"]
    d = rows[:,None]*np.column_stack([one,two])*unit*reference["active_D_prefactor_in_sterile_units"]
    c = unit*np.array([[m,0],[p*s**2,q*x]])
    h = unit*np.diag(old["bare_H0_H6_masses_in_sterile_units"])
    return d,c,h


def compute():
    cfg = read("inputs.json")
    gauge = runpy.run_path(str(HERE/cfg["gauge_round"]))
    g = read(cfg["gauge_inputs"])
    mz,soft = g["mu_ref_GeV"],cfg["ordinary_soft_scale_GeV"]
    sr = read("../1696/results.json")
    col = np.array(sr["color_masses_GeV"])
    tau = 2*math.pi
    central = 4*math.pi/np.array([g[key][0] for key in ["g1_GUT","g2","g3"]])**2
    limits = [(g[k][0]-2*g[k][1],g[k][0]+2*g[k][1]) for k in ["g1_GUT","g2","g3"]]
    corners = [4*math.pi/np.array(row)**2 for row in product(*limits)]

    def need(a, ms):
        return gauge["invert"](a,gauge["branch"](True,ms,mz))
    def impact(ml,ms):
        fermions = np.log(col/ml).sum()
        scalars = .5*np.log((col**2+ms**2)/(ml**2+ms**2)).sum()
        return 2*fermions/3+scalars/3
    def ceiling(a, lower_mass=mz):
        return crossing(lambda ms:impact(lower_mass,ms)-need(a,ms)[1],1000.,1e7)
    ceilings = [ceiling(a) for a in corners]
    length, target = need(central,soft)
    upper = mz*math.exp(length)
    ell = crossing(lambda m:impact(m,soft)-target,mz,1e5)
    lf = np.array([ell,ell])
    sc_col = np.hypot(col,soft)
    sc_ell = np.hypot(lf,soft)
    sc_offset = -gauge["ADJ"]/(12*math.pi)
    events = [(soft,gauge["MSSM"]-gauge["SM"],sc_offset)]
    for masses, coeff in [(col,2*gauge["COLOR"]/3),(lf,2*gauge["WEAK"]/3),
                          (sc_col,gauge["COLOR"]/3),(sc_ell,gauge["WEAK"]/3)]:
        events.extend((m,coeff,0) for m in masses)
    a,scale,b=central.copy(),mz,gauge["SM"].copy()
    for mass,change,shift in sorted(events,key=lambda row:row[0]):
        a=a-b*math.log(mass/scale)/tau+shift
        b+=change
        scale=mass
    a-=b*math.log(upper/scale)/tau
    np.testing.assert_allclose(a,a[0],rtol=0,atol=2e-10)
    old=read(cfg["neutral_threshold"])
    d0,c0,h0=neutral_graph(cfg)
    u=read(cfg["neutral_matching"])["v246_GeV"]/math.sqrt(2)
    htarget=read(cfg["higgs"])["atlas"]["mass_GeV"][0]**2
    moment=lambda d,c:np.linalg.norm(d@c.conj())**2/(8*math.pi**2*u**2)
    b0=(htarget/moment(d0,c0))**(1/3)
    k=cfg["neutral_shape_k"]
    t=k**6
    d,c,h=math.sqrt(b0)/k*d0, b0*k*c0, b0*k**4*h0
    heavy=np.block([[np.zeros_like(c),c],[c.T,h]])
    port=np.column_stack([d,0*d])
    pole2,modes=np.linalg.eigh(heavy.conj().T@heavy)
    mN=np.sqrt(pole2)
    y2=np.sum(np.abs(port@modes/u)**2,axis=0)
    assert max(mN)<soft and min(mN)>mz
    np.testing.assert_allclose(np.dot(pole2,y2)/(8*math.pi**2),htarget,rtol=2e-13)
    # Schur evaluation via the nearly lepton-number conserving block avoids cancellations.
    def active(dd,cc,hh):
        z=np.linalg.solve(cc,dd.T).T
        return z@hh@z.T
    np.testing.assert_allclose(active(d,c,h),active(d0,c0,h0),rtol=1e-13,atol=1e-25)
    F=np.linalg.solve(heavy.T,port.T).T
    W=F@F.conj().T
    tr=np.trace(W).real
    eta_trace=tr/2
    klo,khi=5.5e-4,6.1e-4
    # The two nonnegative trace terms can be bounded independently on this interval.
    c_inv=np.linalg.inv(c0)
    direct=np.linalg.norm(d0@c_inv.T@h0@c_inv)**2
    excess=np.linalg.norm(d0@c_inv.T)**2
    trace_upper=(direct*khi**2+excess/klo**4)/b0
    pole_interval_upper=b0*(khi*np.linalg.norm(c0,2)+khi**4*np.linalg.norm(h0,2))
    interval_yukawa=math.sqrt(b0)/klo*np.linalg.norm(d0,2)/u
    # PSD gives |eta_ab| <= Tr eta / 2 for a != b, in any physical charged basis.
    gf=cfg["general_eta_table6"]
    assert eta_trace<min(gf["lfc_diagonal_upper_95"]) and eta_trace/2<gf["lfv_emu_limit_95"]
    assert eta_trace<gf["preferred_etaee_lower_95"]
    assert trace_upper/4 < gf["lfv_emu_limit_95"] and pole_interval_upper<soft/4 and interval_yukawa<1
    split_local=np.sum(y2*((pole2+soft**2)*np.log1p(soft**2/pole2)-soft**2))/(8*math.pi**2)
    cutoff=cfg["common_flavor_cutoff_GeV"]
    flavor=read("../1685/inputs.json")
    eps=flavor["scalar_magnitude_over_cutoff"]
    vev=eps*cutoff
    cY=np.array([abs(c[1,0])/(eps**2*cutoff),abs(c[1,1])/vev])
    hY=np.diag(h)/vev
    col_yuk=12500/(eps**3*cutoff)
    assert max(cY)<1 and max(hY)<1 and np.linalg.norm(d,2)/u<1
    p_hi=max(mN)
    return dict(soft_ceiling_central_GeV=ceiling(central),soft_ceiling_across_two_widths_GeV=[min(ceilings),max(ceilings)],
        central_soft_ceilings_if_both_weak_poles_GeV_at_least={str(m):ceiling(central,m) for m in [200,300,500]},
        open_neutral_family=dict(k_range=[klo,khi],max_pole_upper_GeV=pole_interval_upper,
            weak_trace_upper=trace_upper,offdiag_eta_upper=trace_upper/4,yukawa_upper=interval_yukawa),
        conservative_fixed_shape_largest_pole_min_GeV=9243.6720155*old["scaled_unit_central_and_window_PeV"][1]*1e3,
        published_changed_shape_largest_pole_GeV=1e6*max(old["fixed_cutoff_changed_heavy_shape_witness"]["poles_PeV"]),
        shared_witness=dict(soft_GeV=soft,required_effective_split=target,
            colored_fermions_GeV=col.tolist(),colored_scalars_GeV=sc_col.tolist(),
            weak_fermions_GeV=lf.tolist(),weak_scalars_GeV=sc_ell.tolist(),
            U_GeV=upper,inverse_alpha_U=a[0],all_gauge_equal_residual=max(a)-min(a),
            neutral_t=t,neutral_b=b0/k**2,sterile_offdiagonal_C_scale_GeV=abs(c).tolist(),
            sterile_H_diagonal_GeV=np.diag(h).tolist(),neutral_poles_GeV=mN.tolist(),
            neutral_pole_Yukawa_squared=y2.tolist(),
            neutral_scalar_poles_if_common_shift_GeV=np.hypot(mN,soft).tolist(),
            neutral_only_Higgs_mass_squared_GeV2=moment(d,c),
            active_singular_masses_eV=(np.linalg.svd(active(d,c,h),compute_uv=False)*1e9).tolist(),
            total_weak_current_weight=tr,upper_eta_diagonal=eta_trace,upper_eta_offdiagonal=eta_trace/2,
            active_largest_Yukawa=np.linalg.norm(d,2)/u,
            sterile_to_soft_mass_ratio=p_hi/soft,
            largest_omitted_relative_Euclidean_propagator_term_bound=(p_hi/soft)**2,
            common_soft_local_CW_abs_H2_GeV2=split_local,
            local_abs_H2_over_target_abs_H2=split_local/(htarget/2),
            common_flavor_cutoff_GeV=cutoff,shared_S_X_B_vev_modulus_GeV=vev,
            color_A3_mixing_abs_yukawas=[col_yuk,col_yuk],
            neutral_A2_X_Yukawa=cY.tolist(),neutral_B_diagonal_Yukawa=hY.tolist(),
            weak_offdiagonal_relative_to_color=ell/12500,weak_diagonal_to_color=0))


def tidy(obj):
    if isinstance(obj,dict):
        return {key:tidy(value) for key,value in obj.items()}
    if isinstance(obj,list):
        return [tidy(value) for value in obj]
    return float(f"{obj:.10g}") if isinstance(obj,(float,np.floating)) else obj


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    cli.add_argument("--check",action="store_true")
    args=cli.parse_args()
    data=tidy(compute())
    if args.check:
        assert data==read("results.json")
        print("1698: soft ceiling, actual charged thresholds, four neutral poles and common vacuum verified")
    else:
        (HERE/"results.json").write_text(json.dumps(data,indent=2)+"\n",encoding="utf-8")
