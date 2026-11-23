"""Unified gauge matching with one adjoint fermion, two split fives, and a neutrino test."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import runpy
import numpy as np

HERE = Path(__file__).resolve().parent
def read(name):
    return json.loads((HERE/name).read_bytes())

CFG=read("inputs.json")
G=runpy.run_path(str(HERE/CFG["physical_matching"]))
LOW=runpy.run_path(str(HERE/CFG["gauge_boxes"]))["inputs"]
PI=math.pi
Z=G["MZ"]
C=G["C"]/(12*PI)
B,W=G["SM"],G["W"]
F3=np.array([0.,4/3,0.])
F8=np.array([0.,0.,2.])
FQ=np.array([10/3,2.,4/3])
X_TEST=CFG["physical_X_test_GeV"]
V=CFG["v_GeV"]
BOUND=read(CFG["light_mass_bounds"])["global_matching_boxes"]["0"]["oscillation"]


def gauge_bounds(data, error, x):
    l=math.log(x/Z)
    lo,hi=data["al"]+C-error,data["ah"]+C+error
    l0=2*PI*(5*hi[0]-3*lo[1]-2*lo[2])/44
    log_g=(44*(l-l0)-2*l)/4
    log_q_over_3=.75*((B[0]-B[1]-8/15-1/3)*l-2*PI*(hi[0]-lo[1]))
    assert log_g>0 and log_q_over_3>0
    g=math.exp(log_g)
    k=(math.sqrt(g)+math.sqrt(g-1))**2
    log_q_over_8=math.log1p(1/k)-math.log(2)
    d5=(2*PI*(hi[0]-lo[2])-10.5*l+2*log_q_over_8)/.4
    sr=(2*k-3)/5
    # The direction defect is defined at a common tree matching with zero
    # additional scalar rescaling, before a possible common multiplicative run.
    mass_budget=Z*((2/3)*k-2)
    defect={order:math.sqrt(2*row["distance_from_rank_one_at_least_eV"]*1e-9*mass_budget)/V
            for order,row in BOUND.items()}
    return dict(log_Q_over_3_lower=log_q_over_3, projected_mass_ratio_G_lower=g,
                mass_8_over_3_lower=k, mass_0_over_3_lower=sr,
                log_product_colored_over_weak_fives_upper=d5,
                equivalent_geometric_mean_fives_mass_ratio_upper=math.exp(d5/2),
                absolute_effective_vector_defect_lower=defect,
                fixed_high_X_triplet_upper_GeV=x/k)


def mass_upper(l0,triplet_floor,light_scalars):
    # Once the 12 constraint selects the light-triplet branch, M8/M3 <= X/floor.
    def excess(l):
        log_ratio=l-math.log(triplet_floor/Z)
        log_g=log_ratio-2*math.log(2)+2*math.log1p(math.exp(-log_ratio))
        return 44*(l-l0)-light_scalars*l-4*log_g
    low,high=20.,45.
    assert excess(low)<0<excess(high)
    for _ in range(90):
        mid=(low+high)/2
        if excess(mid)>0: high=mid
        else: low=mid
    return Z*math.exp(high)


def nuclear_lower(data, error, gate):
    # For positive matter contributions to each one-loop beta-function the
    # actual a_i(mu) is at most its SM reference evolved from the data's upper edge.
    mx=Z
    hi=data["ah"]
    short0=1.
    for _ in range(60):
        l=math.log(mx/Z)
        evolved=hi-B*l/(2*PI)
        logs=np.log(hi/evolved)/B
        short0=math.exp(np.array([23/20,9/4,2])@logs)
        g_min=math.sqrt(4*PI/(evolved[0]+C[0]+error))
        new=gate*g_min*math.sqrt(CFG["common_long_distance_qque"]*short0)
        assert new>=mx*(1-1e-13)
        mx=new
    return mx,short0,g_min


def envelopes(label,error):
    data=LOW(label)
    fit=read(CFG["matter_and_weak_gate"])[label]
    gate=fit["Lambda_geV_lower"]["one_width_shape_box"]
    nlimit,short,gmin=nuclear_lower(data,error,gate)
    # Down-round the converged necessary mass to avoid overstating the iteration.
    nlimit*=1-1e-12
    bound=gauge_bounds(data,error,X_TEST)
    at_gate=gauge_bounds(data,error,nlimit)
    lo,hi=data["al"]+C-error,data["ah"]+C+error
    l0=2*PI*(5*hi[0]-3*lo[1]-2*lo[2])/44
    d,t=fit["D_T_poles_GeV"]
    return dict(at_test=bound,old_log_product_colored_over_weak_fives=2*math.log(t/d),
                original_fives_fixed_excluded_at_and_above_test=bool(2*math.log(t/d)>bound["log_product_colored_over_weak_fives_upper"]),
                X_cap_all_scalars_at_X_above_test_or_unconditional_with_MZ_floor=mass_upper(l0,Z,0),
                X_cap_all_scalars_free_above_test_if_triplet_template=mass_upper(l0,CFG["triplet_direct_template_GeV"],2),
                X_necessary_from_old_flavour_gate_with_gauge_running=nlimit,
                pure_SM_short_gain_floor_at_nuclear_bound=short,g5_lower_from_hypercharge_at_nuclear_bound=gmin,
                triplet_cap_if_old_flavour_gate=at_gate["fixed_high_X_triplet_upper_GeV"],
                min_abs_effective_source_defect_if_gate=at_gate["absolute_effective_vector_defect_lower"])


def example(label, finite=False):
    data=LOW(label)
    a=data["a"]
    cfg=CFG["witness"]
    extra=CFG["finite_threshold_witness"]
    mx=extra["X_GeV"] if finite else X_TEST
    h=np.array(extra["additional_inverse_coupling_thresholds"]) if finite else np.zeros(3)
    matched=a+C-h
    m3=cfg["triplet_GeV"]
    m8=mx*(extra["octet_over_X"] if finite else cfg["octet_over_X"])
    mq=(m3+m8)/2
    m0=(3*m3+2*m8)/5
    hc=mx*cfg["Higgs_colour_over_X"]
    def log(m): return math.log(mx/m)
    l=log(Z)
    l0=2*PI*(W@matched)/44
    sig_log=(44*(l-l0)-4*math.log(mq*mq/(m8*m3)))/2
    sigma=mx*math.exp(-sig_log)
    d5=(.5*sig_log+.1*log(hc)-2*PI*(matched[2]-matched[0])-11.1*l+2*math.log(mq/m8))/.4
    md,mt=(mx*math.exp(-l/2+d5*i/4) for i in (-1,1))
    thresholds=[(m3,F3),(m8,F8),(mq,FQ),(sigma,G["S3"]),(sigma,G["S8"]),(hc,G["HC"]),
                (md,2*G["D"]),(mt,2*G["T"])]
    assert all(Z<=m<=mx for m,_ in thresholds)
    couplings,short=G["run_gauge"](a,mx,thresholds)
    couplings -= h
    np.testing.assert_allclose(couplings,couplings[0],atol=3e-12,rtol=0)
    g=math.sqrt(4*PI/couplings[0])
    total=CFG["common_long_distance_qque"]*short
    reach=mx/g/math.sqrt(total)
    fit=read(CFG["matter_and_weak_gate"])[label]
    gate=fit["Lambda_geV_lower"]["one_width_shape_box"]
    # Independent physical directions, one for each massive eigenstate.
    _,m2,mh=cfg["normal_masses_eV"]
    c=math.sqrt(3/5)
    u0=math.sqrt(2*m2*1e-9*m0)/V
    u3=math.sqrt(2*mh*1e-9*m3)/V
    # Solve with two orthogonal columns of a unitary at a single scale.
    CKM=runpy.run_path(str(HERE.parent/"1805"/"check.py"))["mixing"]
    U=CKM(np.arcsin(np.sqrt([.305,.022,.55])), -1.0)
    low=V*V/2*(np.outer(U[:,1].conj(),U[:,1].conj())*u0*u0/m0+
               np.outer(U[:,2].conj(),U[:,2].conj())*u3*u3/m3)
    np.testing.assert_allclose(np.linalg.svd(low,compute_uv=False)*1e9,[mh,m2,0],atol=1e-15)
    return dict(additional_finite_inverse_coupling_terms=h.tolist(),
                common_nuclear_necessary_bound_satisfied=bool(reach>=gate),
                masses_GeV=dict(X=mx,adjoint_triplet=m3,adjoint_octet=m8,adjoint_leptoquark=mq,
                   adjoint_singlet=m0,Higgs_colour=hc,scalar_triplet_and_octet=sigma,
                   each_weak_five=md,each_color_five=mt),log_product_colored_over_weak_fives=d5,
                inverse_gauge_couplings_at_X=couplings.tolist(),g5=g,short_qque=short,
                full_qque_if_long_gain_1p25=total,effective_qque_scale_GeV=reach,
                required_old_light_flavour_strength_ratio_upper=(reach/gate)**4,
                actual_portal_Y3_if_assigned_upper_light_mass=u3,actual_Y0_if_assigned_solar_mass=u0,
                required_abs_vector_defect_lower_NO=math.sqrt(2*BOUND["NO"]["distance_from_rank_one_at_least_eV"]*1e-9*(m0/(c*c)+m3))/V,
                independent_witness_abs_defect=math.hypot(u0/c,u3),
                parent_norm_needed_if_noncommon_transfer_norm_at_most_1em4=math.sqrt(2*BOUND["NO"]["distance_from_rank_one_at_least_eV"]*1e-9*(m0/(c*c)+m3))/V/1e-4)


def compute():
    return {label:dict(witness=example(label),finite_threshold_witness=example(label,True),
            bounds={str(h):envelopes(label,h)
            for h in CFG["further_inverse_coupling_error"]}) for label in ("common_MZ","independent_lattice")}


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    check=p.parse_args().check
    out=json.dumps(compute(),indent=2,allow_nan=False)+"\n"
    dest=HERE/"results.json"
    if check: assert json.loads(dest.read_bytes())==json.loads(out)
    else: dest.write_text(out,encoding="utf-8")
    print("1814",hashlib.sha256(out.encode()).hexdigest()[:12])
