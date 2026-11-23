"""Heavy triplet threshold and minimum mixed-adjoint flavour coefficient on one physical spectrum."""
import argparse
from itertools import product
import json
import math
from pathlib import Path
import runpy

import numpy as np

HERE = Path(__file__).resolve().parent
CFG = json.loads((HERE / "inputs.json").read_bytes())
G = runpy.run_path(str(HERE / CFG["previous_gauge"]))
PI = math.pi
TRIPLET = np.array([-1., 3., -2.])
GAUGE = np.array([5., -3., -2.])
BS = G["SUSY"] + 4
Z = G["MZ"]
SOFT = G["SOFT"]


def mass_spectrum(original):
    if original:
        return G["floor_spectrum"](G["OLD"]["weak_fermions_GeV"][0], G["NEW_WEAK"])
    return G["floor_spectrum"](G["CFG"]["step_mass_floor_GeV"], G["MC"])


def matching(weak, displacements=(0, 0, 0, 0)):
    s1, s2, s3, st = displacements
    a5, da5 = G["CFG"]["strong_lattice_5f"]["ALPHA_2026"]
    mt, dmt = G["CFG"]["top_MS_matching_GeV"]
    a3 = 1/(a5+s3*da5) + math.log((mt+st*dmt)/Z)/(3*PI)
    low = np.r_[4*PI/G["electroweak"](s1,s2)**2, a3]
    offset = (G["ADJ"]/(12*PI)+(G["SM"]-G["SUSY"])*math.log(SOFT/Z)/(2*PI)
              -(G["C"]*sum(G["logmass"](G["COLOR"],SOFT))
                +G["L"]*sum(G["logmass"](weak,SOFT)))/(2*PI))
    reduced = low-offset
    q = Z*math.exp(2*PI*(GAUGE@reduced)/(GAUGE@BS))
    mtc = Z*math.exp(2*PI*(TRIPLET@reduced)/(TRIPLET@BS))
    inverse_at_q = reduced-BS*math.log(q/Z)/(2*PI)
    # In this branch the X and adjoint threshold logarithms vanish at q.
    a_unified = inverse_at_q-G["C"]*math.log(q/mtc)/(2*PI)
    np.testing.assert_allclose(a_unified, [a_unified[1]]*3, rtol=2.e-14)
    return q, mtc, math.sqrt(4*PI/a_unified[1]), a3, inverse_at_q


def gauge_box(weak, widths):
    vals = [matching(weak, shifts)[:3] for shifts in product(*[(-widths,widths)]*4)]
    vals = np.array(vals)
    return {"MX_squared_times_MSigma_cube_root_GeV": [vals[:,0].min(),vals[:,0].max()],
            "effective_triplet_mass_GeV": [vals[:,1].min(),vals[:,1].max()],
            "g5_on_equal_heavy_branch": [vals[:,2].min(),vals[:,2].max()]}


def charged_point(full_a3):
    archived = runpy.run_path(str(HERE / CFG["charged_background"]))
    bg = archived["mother"](CFG["step_low"])
    model = bg["model"]
    old = model["CFG"]["nearby_measured_inputs_in_marginal_widths"][archived["CFG"]["parent_point"]]
    g3central, g3error = G["LOW"]["g3"]
    g3 = math.sqrt(4*PI/full_a3)
    new_g3_shift = (g3-g3central)/g3error
    _, pt, _ = model["charged_endpoint"](2,CFG["step_low"],old["weak"],old["up"],new_g3_shift,
                                         old["top"],old["charm"])
    y = np.zeros((4,3,4),dtype=complex)
    y[0,:,:3], y[1,:,:3], y[2,:,:3] = np.diag(pt["A"]), pt["D"], pt["E"].T
    y[3] = pt["high_N"]
    return y, pt, {"g3_full_MZ":g3, "g3_shift_in_original_widths":new_g3_shift,
                   "tau_calibrated_tan_beta":pt["row"]["positive_zero_stop_tan_beta"],
                   "gamma_MZ_deg":pt["rep"]["gamma_MZ_deg"]}


def run_yukawa(y, start, q, inverse_at_q, beta, step):
    length = math.log(q/start)
    count = math.ceil(length/step)
    dl = length/count
    dt = dl/(16*PI**2)
    def gauge_at(u):
        return 4*PI/(inverse_at_q+BS*(length-u)/(2*PI))
    for j in range(count):
        a, b, c = (gauge_at((j+k)*dl) for k in (0.,.5,1.))
        k1 = beta(y,a,True)
        k2 = beta(y+.5*dt*k1,b,True)
        k3 = beta(y+.5*dt*k2,b,True)
        k4 = beta(y+dt*k3,c,True)
        y = y+dt*(k1+2*k2+2*k3+k4)/6
    return y


def insertions(d, e, v_over_cutoff):
    delta = d-e
    # With a mixed insertion H, the four matrix entries on (a,b)=(+2,-3)^2
    # define a bilinear polynomial. Multiplying the charged matrices by (2-b)/5
    # keeps both mass Yukawas and zeros both triplet Yukawas exactly.
    f = delta/5
    one = e+3*f
    h_mixed = -f/5
    coeff = [2*one/5, 2*f/5, -one/5, h_mixed]
    def coupling(a,b):
        y0,yf,yh,yfh = coeff
        return y0+a*yf+b*yh+a*b*yfh
    for a,b,target in ((2,-3,d),(-3,-3,e),(2,2,0*d),(-3,2,0*e)):
        np.testing.assert_allclose(coupling(a,b),target,atol=2e-13)
    sd,se = [np.linalg.svd(t,compute_uv=False) for t in (d,e)]
    return dict(quark_yukawa_singular_values=sd,lepton_yukawa_singular_values=se,
                norm_delta_2=np.linalg.norm(delta,2),
                minimal_delta_2_if_both_relative_unitary_frames_free=max(abs(sd-se)),
                exact_dimensionless_mixed_coefficient_norm=np.linalg.norm(h_mixed,2)/v_over_cutoff**2,
                dimensionless_lower_bound_all_relative_frames=max(abs(sd-se))/(25*v_over_cutoff**2),
                minimum_unweighted_sum_norms_squared_single_insertion=np.linalg.norm(delta)**2/2)


def compute():
    rows = {}
    for name,orig in [("archived_weak_masses",True),("1769_floor_endpoint",False)]:
        weak = mass_spectrum(orig)
        q, mtc, g5, _, _ = matching(weak)
        rows[name] = dict(weak_mass_GeV=weak, equal_MX_MSigma_GeV=q, MT_GeV=mtc,
                          g5=g5, abs_lambdaH=g5*mtc/q,
                          gauge_boxes={str(n):gauge_box(weak,n) for n in (1,2)})
    q, mtc, g5, a3, aq = matching(mass_spectrum(True))
    y, pt, charged = charged_point(a3)
    first = CFG["flavour_scale_GeV"]
    # Require the low-scale shooting and the gauge matching to use exactly the same material.
    g2_first = pt["flow"]["segments"](CFG["step_low"])[-1][-1][-1]
    np.testing.assert_allclose(g2_first,4*PI/(aq+BS*math.log(q/first)/(2*PI)),rtol=2e-11)
    high = run_yukawa(y,first,q,aq,pt["flow"]["beta"],CFG["step_high"])
    fine = run_yukawa(y,first,q,aq,pt["flow"]["beta"],CFG["step_high"]/2)
    np.testing.assert_allclose(high,fine,atol=3e-12,rtol=1e-10)
    vev = q/(5*g5)
    eps = vev/CFG["cutoff_for_adjoint_expansion_GeV"]
    info = insertions(high[1,:,:3],high[2,:,:3].T,eps)
    budget = CFG["aggregate_quadratic_coefficient_budget"]
    info["largest_cutoff_GeV_at_budget_tagged"] = 5*vev*math.sqrt(budget/info["norm_delta_2"])
    info["largest_cutoff_GeV_at_budget_free_relative_frames"] = 5*vev*math.sqrt(
        budget/info["minimal_delta_2_if_both_relative_unitary_frames_free"])
    info["minimum_max_triplet_yukawa_norm_at_Planck_budget_tagged"] = max(
        0.,(info["norm_delta_2"]-25*eps**2*budget)/2)
    info["minimum_max_triplet_yukawa_norm_at_Planck_budget_free_frames"] = max(
        0.,(info["minimal_delta_2_if_both_relative_unitary_frames_free"]-25*eps**2*budget)/2)
    # For every mixed monomial multiply its independently bounded coefficient by
    # the product of the two eigenvalue differences. Their absolute sums factorise.
    response = 5*eps/((1-3*eps)*(1+2*eps))
    np.testing.assert_allclose(response,sum(abs(2**n-(-3)**n)*eps**n for n in range(1,16)),
                               rtol=2.e-15)
    for label,delta in (("tagged",info["norm_delta_2"]),
                        ("free_frames",info["minimal_delta_2_if_both_relative_unitary_frames_free"])):
        target = math.sqrt(delta/budget)
        eps_required = 2*target/(5+target+math.sqrt((5+target)**2+24*target**2))
        info["all_orders_"+label] = dict(
            necessary_uniform_mixed_coefficient_budget=delta/response**2,
            largest_cutoff_GeV_at_budget=vev/eps_required,
            minimum_max_triplet_yukawa_norm=max(0.,(delta-budget*response**2)/2))
    return dict(gauge=rows, charged_input_refit=charged,
                common_frame_yukawa_matching=dict(Q_GeV=q,M_T_GeV=mtc,adjoint_vev_GeV=vev,
                    V_over_cutoff=eps,low_norm_delta=np.linalg.norm(y[1,:,:3]-y[2,:,:3].T,2),
                    up_singular_values_Q=np.linalg.svd(high[0,:,:3],compute_uv=False),**info))


def serial(x):
    if isinstance(x,dict):
        return {k:serial(v) for k,v in x.items()}
    if isinstance(x,(list,tuple,np.ndarray)):
        return [serial(v) for v in x]
    return float(f"{x:.9g}")


if __name__ == "__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    check=p.parse_args().check
    result=serial(compute())
    path=HERE/"results.json"
    if check:
        assert result==json.loads(path.read_bytes())
        print("1771 OK: heavy gauge matching, full 10 PeV-to-GUT Yukawa transport, mixed-adjoint certificate")
    else:
        path.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf-8")

