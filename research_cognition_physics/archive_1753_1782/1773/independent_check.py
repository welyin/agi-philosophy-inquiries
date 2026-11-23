"""Cross-test the same CKM invariant in kaons and direct charm CP violation."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
CFG = json.loads((HERE / "independent_inputs.json").read_bytes())


def read(key):
    return json.loads((HERE / CFG[key]).read_bytes())


def unitary_from_tree(row):
    """Recover the standard unitary parametrization from three moduli and true gamma."""
    us, ub, cb = (row["moduli"][key] for key in ("Vus", "Vub", "Vcb"))
    c13 = math.sqrt(1 - ub**2)
    s12, s23 = us / c13, cb / c13
    c12, c23 = math.sqrt(1-s12*s12), math.sqrt(1-s23*s23)
    gamma = math.radians(row["gamma"])
    delta = gamma + math.asin(c12*s23*ub*math.sin(gamma)/(s12*c23))
    e = np.exp(1j*delta)
    v = np.array([[c12*c13, s12*c13, ub*e.conjugate()],
                  [-s12*c23-c12*s23*ub*e, c12*c23-s12*s23*ub*e, s23*c13],
                  [s12*s23-c12*c23*ub*e, -c12*s23-s12*c23*ub*e, c23*c13]])
    np.testing.assert_allclose(v.conj().T@v, np.eye(3), atol=1.e-14)
    assert abs(np.angle(-v[0,0]*v[0,2].conj()/(v[1,0]*v[1,2].conj()))-gamma)<1.e-12
    return v


def charm_response(v):
    # A_f = lambda_q T_f + lambda_b P_f, q=s (KK) or d (pipi).
    lam = v[1].conjugate()*v[0]
    q = lam[[1, 0]]
    weak = lam[2] / q
    j = (q*lam[2].conjugate()).imag
    np.testing.assert_allclose(j[0], -j[1], rtol=2.e-12)
    return abs(j[0]), weak


def asymmetry_limit(budget, weak):
    # Under an upper bound on |P/T|, r can be tuned up to one.
    r = np.minimum(budget*abs(weak), 1.0)
    sin_weak = np.abs(weak.imag) / abs(weak)
    return 2*r*sin_weak / np.hypot(1-r*r, 2*r*sin_weak)


def min_cost_point(asym, weak):
    y = abs(asym)/np.sqrt(1-asym**2)
    sin_weak = abs(weak.imag) / abs(weak)
    r = y / (np.hypot(y, sin_weak) + sin_weak)
    cos_delta = -2*r*(weak.real/abs(weak))/(1+r*r)
    phase = -np.sign(asym*weak.imag)*np.arccos(np.clip(cos_delta, -1, 1))
    strong = r/abs(weak)*np.exp(1j*phase)
    direct = abs(1+weak*strong)**2
    conjugate = abs(1+weak.conjugate()*strong)**2
    np.testing.assert_allclose((direct-conjugate)/(direct+conjugate), asym, atol=2.e-15)
    return np.abs(strong)


def rectangle_distance(cap, mean, cov):
    """Minimum squared Mahalanobis distance over the full two-mode attainable rectangle."""
    if np.all(np.abs(mean) <= cap):
        return 0.0
    candidates = []
    for i in range(2):
        k = 1-i
        for sign in (-1, 1):
            x = np.empty(2)
            x[i] = sign*cap[i]
            x[k] = np.clip(mean[k]+cov[k,i]/cov[i,i]*(x[i]-mean[i]), -cap[k], cap[k])
            diff = x-mean
            candidates.append(diff @ np.linalg.solve(cov, diff))
    return min(candidates)


def charm_budget(v, mean, cov):
    j, weak = charm_response(v)
    point = min_cost_point(mean, weak)
    assert np.max(abs(asymmetry_limit(point,weak)-abs(mean))) < 1.e-14
    limits = {}
    for cutoff in (1., 2.30, 5.99):
        lo,hi = 0.,max(point)
        assert rectangle_distance(asymmetry_limit(lo,weak),mean,cov)>cutoff
        for _ in range(60):
            mid = (lo+hi)/2
            if rectangle_distance(asymmetry_limit(mid,weak),mean,cov)>cutoff:
                lo=mid
            else:
                hi=mid
        limits[str(cutoff)]=hi
    return {"J_abs":j,"individual_min_abs_P_over_T":point.tolist(),
            "joint_min_max_abs_P_over_T_for_delta_chi2":limits}


def rare_rate(z_t, ud, transfer):
    ref = transfer["lambda_reference"]
    x = transfer["Xtop"][0]
    p = ref**4*transfer["Pcharm"][0]
    z_c = -1-z_t
    scale = (ud/ref**4)**2
    charged = (transfer["kplus_e11"][0]*(1+transfer["delta_EM"])*scale
               *((x*z_t.imag)**2+(x*z_t.real+p*z_c.real)**2))
    neutral = transfer["kL_e11"][0]*scale*(x*z_t.imag)**2
    return charged,neutral


def full_inverse(ud, beta, beta_s, rho, rate, transfer):
    # With the d column real V_is= lambda_i^* / V_id.
    w = rho*np.exp(-1j*(beta+beta_s))
    zc,zt = -1/(1+w), -w/(1+w)
    norm = 1-ud**2
    const = abs(zc*zt)*math.sin(beta_s)/math.sin(beta)
    # x = |Vtd|^2, x(norm-x) from the product of two J-invariant formulae.
    poly = [-1/ud**2, norm/ud**2+abs(zc)**2-abs(zt)**2, abs(zt)**2*norm-const]
    out = []
    for root in np.roots(poly):
        if abs(root.imag)>1.e-12 or not 0<root.real<norm:
            continue
        t = root.real
        d = np.sqrt([ud**2,norm-t,t])
        u = math.sqrt((norm-t)*t / const)
        s = u*np.conjugate([1.,zc,zt])/d
        b = np.cross(d, s.conj())
        v = np.column_stack((d,s,b))
        np.testing.assert_allclose(v.conj().T@v,np.eye(3),atol=1.e-10)
        rebuilt_b = np.angle(-v[1,0]*v[1,2].conj()/(v[2,0]*v[2,2].conj()))
        rebuilt_s = np.angle(-v[2,1]*v[2,2].conj()/(v[1,1]*v[1,2].conj()))
        if not np.allclose([rebuilt_b,rebuilt_s],[beta,beta_s],atol=1.e-8,rtol=0):
            continue
        charged,neutral = rare_rate(zt,ud,transfer)
        assert abs(charged-rate)<2.e-9
        j, _ = charm_response(v)
        assert abs(u*u*zt.imag/j-1)<1.e-8
        out.append((v,neutral))
    return out


def invert_kaon(ud,beta,beta_s,charged,transfer):
    theta = beta+beta_s
    def evaluate(rho):
        w = rho*np.exp(-1j*theta)
        return rare_rate(-w/(1+w),ud,transfer)[0]
    lo,hi = 0.,.02
    assert evaluate(lo)<charged<evaluate(hi)
    for _ in range(60):
        mid=(lo+hi)/2
        if evaluate(mid)<charged:
            lo=mid
        else:
            hi=mid
    return full_inverse(ud,beta,beta_s,(lo+hi)/2,charged,transfer)


def run():
    mean = np.array(CFG["LHCb_direct_central"])
    err = np.array(CFG["LHCb_error"])
    corr = CFG["LHCb_correlation"]
    cov = np.outer(err,err)*np.array([[1,corr],[corr,1]])
    transfer = read("kaon_transfer")
    variants = {}
    for key, row in read("archived_fit")["refits"].items():
        v = unitary_from_tree(row)
        jboard = charm_budget(v,mean,cov)
        sd = v[:,1].conj()*v[:,0]
        charged,neutral = rare_rate(sd[2]/sd[0],abs(v[0,0]),transfer)
        assert abs(jboard["J_abs"]-abs((sd[2]*sd[0].conj()).imag))<1.e-15
        lam = v[1].conj()*v[0]
        q = abs(lam[[1,0]])
        b = abs(lam[2])
        # Illustrative independent cap |P/T| <= 1, fixed weak moduli.
        j_required_unit_cap = abs(mean)/np.sqrt(1-mean**2)*(q**2-b**2)/2
        minimum_neutral = neutral*(j_required_unit_cap/jboard["J_abs"])**2
        floor_any_weak_angle = (abs(mean)/(1+np.sqrt(1-mean**2)))*q/b
        variants[key] = {**jboard,"KL_e11":neutral,"Kplus_e11":charged,
                         "min_abs_P_over_T_even_at_optimal_weak_angle": floor_any_weak_angle.tolist(),
                         "necessary_KL_e11_at_unit_hadronic_cap_per_mode":minimum_neutral.tolist()}
    angles = read("B_phases")["data"]
    beta = math.asin(angles["golden_S"][0])/2
    ud = read("nuclear_beta")["ud_beta"][0]
    charged = read("kaon_rate")["B_rare"]["central"]*1.e11
    inverse = {}
    phases = {k:v[0] for k,v in read("comparison_Bs")["phi_s_rad"].items()}
    phases["archived_average"] = angles["phi_s_ccs"][0]
    for key,phi_s in [*phases.items(),("main_inputs",angles["phi_s_ccs"][0])]:
        used_ud = read("main_nuclear_beta")["holdout_nuclear"][0] if key=="main_inputs" else ud
        inverse[key] = {}
        for tag,angle in (("small",beta),("reflected",math.pi/2-beta)):
            branches = invert_kaon(used_ud,angle,-phi_s/2,charged,transfer)
            assert len(branches)==1, (key, tag, len(branches))
            v,neutral = branches[0]
            inverse[key][tag] = {**charm_budget(v,mean,cov),"KL_e11":neutral,
                                 "Vus":abs(v[0,1]),"Vub":abs(v[0,2]),"Vcb":abs(v[1,2])}
    m=CFG["two_body_mass_GeV"]
    phase_ratio=math.sqrt((m["D0"]**2-4*m["charged_pi"]**2)/(m["D0"]**2-4*m["charged_K"]**2))
    ratio, stat, syst = CFG["FOCUS_width_ratio_KK_over_pipi"]
    weight = ratio*phase_ratio
    null_combination = weight*mean[0]+mean[1]
    weights=np.array([weight,1.0])
    error=math.sqrt(weights@cov@weights + (mean[0]*phase_ratio)**2*(stat**2+syst**2))
    return {"six_frozen_material_branches":variants,"kaon_angle_inversions":inverse,
            "charm_common_interference_test":{"phase_space_p_pi_over_p_K":phase_ratio,
            "central_Im_TPstar_KK_over_pipi":-weight*mean[0]/mean[1],
            "common_reduced_interference_null_weighted_asymmetry":null_combination,
            "propagated_error":error,"conditional_standardized_residual":null_combination/error}}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    check=parser.parse_args().check
    # Round only when serializing platform-dependent trigonometric results.
    compact=json.loads(json.dumps(run()),parse_float=lambda s:float(f"{float(s):.10g}"))
    dest=HERE/"independent_results.json"
    if check:
        assert compact==json.loads(dest.read_bytes())
    else:
        with dest.open("x", encoding="utf8") as f:
            json.dump(compact,f,indent=2,allow_nan=False)
            f.write("\n")
    print("1773 independent: unitary completions, invariant kaon-charm test and attained hadronic bounds")
