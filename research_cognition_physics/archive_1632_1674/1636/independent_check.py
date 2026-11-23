"""Independent amplitude-space check of the two fixed round-1635 weak matrices."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
BASE = HERE.parent


def read(round_number, filename="inputs.json"):
    return json.loads((BASE / str(round_number) / filename).read_text(encoding="utf-8-sig"))


def interval(measurement, n):
    radius = n * sum(measurement[1:])
    return [measurement[0] - radius, measurement[0] + radius]


def limits(values):
    return [min(values), max(values)]


def covered_distance(a, b):
    return max(0., a[0] - b[1], b[0] - a[1])


def calculate(matrix_rows, md):
    old, ref, hp = read(1468), read(1404), read(1407)
    fit, match, rates = read(1635)["witness"], read(1464)["HPQCD_tableI_eq15"], read(1464, "results.json")
    k = read(1406)["sm_coefficient_ratio"]["center"]
    raw = np.asarray(matrix_rows)
    mass = raw[:, :, 0] + 1j * raw[:, :, 1]
    u, singular, vh = np.linalg.svd(mass)
    order = np.argsort(singular)
    u, v, singular = u[:, order], vh.conj().T[:, order], singular[order]
    s, b = 1, 2
    ms, mb = singular[s:b+1]
    vev, mh = fit["Higgs_vev_GeV"], fit["h_mass_GeV"]
    weak_y = np.diag([1, 1, 1, 0, 1]) @ mass
    weak_y[4, 4] -= md
    physical_y = u.conj().T @ weak_y @ v / vev

    def product(row):
        return row[s].conjugate() * row[b]

    lt, lc = product(u[2]), product(u[1])
    fl, fr = -product(u[3]) - 2*product(u[4]), -product(v[4])
    t = (fl/fr).real
    tratio = 2*ms*mb/md**2
    assert abs(fl/fr-tratio) < 2e-10
    assert abs(physical_y[s, b]*vev/fr-ms) < 2e-7
    assert abs(physical_y[b, s].conjugate()*vev/fr-mb) < 2e-7
    const = ref["const_GeV"]
    gf, mw, meson = const["GF"], const["MW"], const["M_Bs"]
    norm = 2*meson/(3*hp["hbar_GeV_ps"])
    sc = ref["short_distance"]
    bag = hp["HPQCD_inputs"]["fBs_sqrtB_GeV"]
    base = gf**2 * mw**2 * lt**2/(4*math.pi**2)
    tree = gf*fr**2/math.sqrt(2)
    s0, eb = sc["S0_published"], sc["etaB"]
    ez = old["alpha_s_ref"]**(6/23)
    a = base*s0[0]*eb[0]
    coef = a+ez*(1+t*t)*tree
    delta = np.angle(coef/a)
    phi = np.angle(coef/lc**2)
    gf_v = 1/(math.sqrt(2)*gf*vev**2)
    corrections = {"V_LR": 2*t, "S_LL": -gf_v*mb**2/mh**2,
                   "S_RR": -gf_v*ms**2/mh**2, "S_LR": -gf_v*2*mb*ms/mh**2}
    leverage = sum(map(abs, corrections.values()))
    lepton = -lt*match["Xt_ref"]/(math.pi*match["sine_squared_ref"]*match["alpha_inverse_ref"])
    p = k*lepton+fl-fr
    weights = {key: rates[key] for key in ("calculated_light_response", "run3_central_light_response")}
    h = rates["Hmu_heavy"]*1e9*abs(p)**2
    decay = {"P_modulus":abs(p), "pure_heavy":h}
    for key, weight in weights.items():
        d = np.angle(p/lt)-delta/2
        corrected = h*(1-(1-weight)*math.sin(d)**2)
        before = h*(1-(1-weight)*math.sin(np.angle(p/lt))**2)
        decay[key] = {"corrected_1e9":corrected, "before_1e9":before}
    out = {"M_D_GeV":md, "abs_FL":abs(fl), "abs_FR":abs(fr), "t":t,
           "sum_scalar_LR_leverage":leverage, "tree_coefficient_weights":corrections,
           "h_sb_times_v_over_FR":(physical_y[s,b]*vev/fr).real,
           "h_bs_conj_times_v_over_FR":(physical_y[b,s].conjugate()*vev/fr).real,
           "phase_SM":np.angle(a/lc**2), "central_phase":phi, "central_delta":delta,
           "central_tree_over_top_re_im":[(ez*tree/a).real,(ez*tree/a).imag],
           "central_clock": norm*bag[0]**2*abs(coef),
           "central_clock_required_fsqrtB":math.sqrt(old["omega_ps_inverse"][0]/(norm*abs(coef))),
           "decay_feedback":decay}
    for n in (1, 2):
        obs = old["phase_inputs_rad"]
        kk = [interval(obs[name],n) for name in ("LHCb_KK_2015_18","CMS_KK_13TeV_final")]
        control = old["extra_decay_phase_cap_rad"]
        target = [max(x[0] for x in kk)-control,min(x[1] for x in kk)+control]
        strengths = [x*y for x,y in itertools.product(interval(s0,n),interval(eb,n))]
        coeffs = [base*x+tree*ez0*(1+t*t) for x,ez0 in itertools.product(limits(strengths),old["eta_Z_RGI_range"])]
        phis = limits([np.angle(z/lc**2) for z in coeffs])
        # The norm extrema occur at the corners here: both radial derivatives stay positive.
        radial = lambda z,q: (z.conjugate()*q).real
        assert all(radial(z,base)>0 for z in coeffs)
        increasing = [radial(z,tree)>0 for z in coeffs]
        assert all(increasing) or not any(increasing)
        clocks = limits([norm*f**2*abs(z) for f,z in itertools.product(interval(bag,n),coeffs)])
        eta_feasible = []
        for target_phi, strength in itertools.product(target,limits(strengths)):
            rotation = np.exp(-1j*target_phi)/lc**2
            e = -(rotation*base*strength).imag/(rotation*tree).imag
            assert (rotation*(base*strength+e*tree)).real>0
            eta_feasible.append(e)
        affine_eta = limits(eta_feasible)
        default = [z*(1+t*t) for z in old["eta_Z_RGI_range"]]
        gap = covered_distance(affine_eta,default)
        phase_gap = covered_distance(phis,target)
        clock_measured=interval(old["omega_ps_inverse"],n)
        out[f"width_{n}"] = {"untagged_common_KK_phase_target_rad":target, "leading_phase":phis,
            "phase_gap_rad":phase_gap,"extra_arbitrary_complex_min_relative_modulus":math.sin(phase_gap),
            "leading_clock":clocks,"clock_data":clock_measured,
            "clock_range_gap":covered_distance(clocks,clock_measured),
            "aligned_eta_needed_for_phase":affine_eta,"eta_gap":gap,
            "LR_and_h_only_minimum_max_abs_transport":gap/leverage}
    # Table 5: the (3,2,-5/6) Y^4 contraction is +1/4 in the RR basis.
    # Use two explicit light-flavour identifications before comparing approximation errors.
    assert np.linalg.norm(mass[3, :3])+np.linalg.norm(mass[:4, 4]) < 1e-9
    u0, w0, vh0 = np.linalg.svd(mass[:3, :3])
    low_order = np.argsort(w0)
    u0, v0 = u0[:, low_order], vh0.conj().T[:, low_order]
    overlaps = np.diag(u0.conj().T @ u[:3, :3])
    # Same individual s,b phase conventions in the two calculations.
    align = np.diag(overlaps/abs(overlaps))
    lam_unbroken = (math.sqrt(2)/vev) * mass[4, :3] @ v0 @ align
    lam_physical = (math.sqrt(2)/vev) * mass[4, :3] @ v[:3, :3]
    l_power = product(lam_physical)
    assert abs(l_power /(-fr*2*md**2/vev**2)
               - (1-ms**2/md**2)*(1-mb**2/md**2))<2e-8
    a5 = old["alpha_s_ref"]
    mt = sc["mtMS"][0]
    mz = hp["running"]["MZ"]
    a6 = a5/(1+a5*(23/3)/(2*math.pi)*math.log(mt/mz))
    am = a6/(1+a6*7/(2*math.pi)*math.log(md/mt))
    qcd_rr = (am/a6)**(2/7)*(a6/a5)**(6/23)
    y_norm = np.linalg.svd(weak_y*(math.sqrt(2)/vev),compute_uv=False)[0]
    heavy = {"largest_yukawa":y_norm, "sqrt_4pi":math.sqrt(4*math.pi),
             "Ymax_squared_over_16pi_squared": y_norm**2/(16*math.pi**2),
             "qcd_RR_heavy_over_MZ":qcd_rr,
             "light_state_left_phase_alignment_moduli":abs(overlaps).tolist()}
    for name, lambdas in (("dim6_unbroken_flavour",lam_unbroken),("exact_external_projection",lam_physical)):
        lam = product(lambdas)
        br = lam**2/(64*math.pi**2*md**2)
        ratio_to_tree = br/tree
        central = coef + ez*qcd_rr*br
        gr = {"abs_Lambda_sb":abs(lam),"bare_box_over_tree":[ratio_to_tree.real,ratio_to_tree.imag],
              "central_clock": norm*bag[0]**2*abs(central),
              "central_phase":np.angle(central/lc**2),"central_delta":np.angle(central/a),
              "f_sqrtB_to_match_clock": math.sqrt(old["omega_ps_inverse"][0]/(norm*abs(central)))}
        for key, weight in weights.items():
            angle = np.angle(p/lt)-np.angle(central/a)/2
            gr[key+"_1e9"] = h*(1-(1-weight)*math.sin(angle)**2)
        for n in (1,2):
            # Independently vary an illustrative multiplicative self-box transfer of 20%.
            # This interval does not include the separate Yukawa-induced LR operator by definition.
            obs = old["phase_inputs_rad"]
            kk = [interval(obs[key],n) for key in ("LHCb_KK_2015_18","CMS_KK_13TeV_final")]
            target = [max(x[0] for x in kk)-old["extra_decay_phase_cap_rad"],
                      min(x[1] for x in kk)+old["extra_decay_phase_cap_rad"]]
            points = []
            for b1,b2,z,scale in itertools.product(interval(s0,n),interval(eb,n),
                                                     old["eta_Z_RGI_range"],[0.8,1.2]):
                val = base*b1*b2+z*((1+t*t)*tree+qcd_rr*br*scale)
                assert (val.conjugate()*base).real>0 and (val.conjugate()*br).real>0
                assert (val.conjugate()*(tree*(1+t*t)+qcd_rr*br*scale)).real>0
                points.append(val)
            phases = limits([np.angle(z/lc**2) for z in points])
            c = limits([norm*f**2*abs(z) for f,z in itertools.product(interval(bag,n),points)])
            # Complex phase is monotone on each segment of the positive-coefficient hull.
            phase_gap = covered_distance(phases,target)
            gr[f"width_{n}"] = {"phase_rad":phases,"target_rad":target,"clock":c,
                               "phase_gap":phase_gap,
                               "extra_complex_min_relative":math.sin(phase_gap),
                               "clock_gap":covered_distance(c,interval(old["omega_ps_inverse"],n))}
        # Published mixed-chirality normalization, table 2 and eqs. 35, 37, 61
        # of 1703.04753. Refer its LR amplitude to its own SM top amplitude.
        xt = sc["xt"][0]
        H1_MW = -(xt-7)/(4*(xt-1)) - (xt**2-2*xt+4)*math.log(xt)/(2*(xt-1)**2)
        log_combo = math.log(md/mw)+H1_MW
        lt0 = product((u0@align)[2]) if name=="dim6_unbroken_flavour" else lt
        dC = -vev**2*lam*xt*log_combo/(2*md**2*lt0)
        f1, f4, f5 = 0.0498, 0.0534, 0.0493
        lr_by_top = -dC*(1.18*f4+1.42*f5)/(1.95*f1)
        # The Wilson coefficient is proportional to the counterterm in the
        # unbroken basis; convert back to the finite top convention explicitly.
        lr_by_top *= (lt0/lt)**2
        combined = central+a*lr_by_top
        gr["published_top_Yukawa_LR_central_diagnostic"] = {
            "H1_MW":H1_MW, "Delta_C_LR1":[dC.real,dC.imag],
            "LR_over_finite_top":[lr_by_top.real,lr_by_top.imag],
            "phase_rad":np.angle(combined/lc**2), "clock":norm*bag[0]**2*abs(combined),
            "f_sqrtB_to_match_clock":math.sqrt(old["omega_ps_inverse"][0]/(norm*abs(combined)))}

        # Widen the three published hadronic combinations independently. This
        # outer rectangle deliberately relaxes their correlation with f sqrt(B).
        def lr_bag_range(n):
            return limits([(1.18*v4+1.42*v5)/v1/((1.18*f4+1.42*f5)/f1)
                           for v1,v4,v5 in itertools.product(interval([f1,.0032],n),
                               interval([f4,.0032],n),interval([f5,.0037],n))])
        for n in (1,2):
            xs = interval(sc["xt"],n)
            log_ratios = []
            for x in xs:
                h1 = -(x-7)/(4*(x-1))-(x*x-2*x+4)*math.log(x)/(2*(x-1)**2)
                log_ratios.append(x*(math.log(md/mw)+h1)/(xt*log_combo))
            # Positive derivative of x [ln(M/MW)+H1(x)] on these short intervals:
            # dH1/dx = -1.5/(x-1)^2 +3 ln(x)/(x-1)^3 -(1+3/(x-1)^2)/(2x).
            assert math.log(md/mw)-0.75-xs[1]*(1.5/(xs[0]-1)**2+
                      (1+3/(xs[0]-1)**2)/(2*xs[0])) > 0
            lr_range = limits([g*j for g,j in itertools.product(lr_bag_range(n),log_ratios)])
            vectors = []
            for b1,b2,z,mult,r in itertools.product(interval(s0,n),interval(eb,n),
                                          old["eta_Z_RGI_range"],[.8,1.2],lr_range):
                val = base*b1*b2 + a*lr_by_top*r + z*(tree*(1+t*t)+qcd_rr*br*mult)
                vectors.append(val/lc**2)
            red = limits([w.real for w in vectors])
            imd = limits([w.imag for w in vectors])
            assert red[0]>0 and imd[0]>0
            pr = limits([math.atan2(y,x) for x,y in itertools.product(red,imd)])
            cr = limits([norm*ff**2*abs(lc)**2*math.hypot(x,y) for ff,x,y
                         in itertools.product(interval(bag,n),red,imd)])
            gr["published_top_Yukawa_LR_central_diagnostic"][f"width_{n}"] = {
                "complex_rectangle_phase_rad_outer":pr, "complex_rectangle_clock_outer":cr,
                "clock_gap":covered_distance(cr,interval(old["omega_ps_inverse"],n)),
                "phase_gap":covered_distance(pr,gr[f"width_{n}"]["target_rad"])}

        # Absolute lattice calibration of 1636; equation (61) supplies the LR
        # term directly, without multiplication by the borrowed SM ratio.
        raw_lr = -(base/lt**2)*lt0**2*dC*(1.18*f4+1.42*f5)
        absolute = central + raw_lr/bag[0]**2
        gr["absolute_LR_lattice_calibration"] = {
            "phase":np.angle(absolute/lc**2),"clock":norm*bag[0]**2*abs(absolute),
            "f2hat_SM_over_own_2017_SM_ratio":bag[0]**2*eb[0]*s0[0]/(1.95*f1)}

        heavy[name]=gr
    out["positive_RR_self_box"] = heavy

    return out


def rounded(obj):
    if isinstance(obj, dict):
        return {key: rounded(value) for key,value in obj.items()}
    if isinstance(obj,list):
        return [rounded(value) for value in obj]
    if isinstance(obj,(float,np.floating)):
        return float(f"{obj:.9g}")
    return obj


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--check",action="store_true")
    args=p.parse_args()
    source = read(1635,"results.json")
    data = {name:calculate(source[name]["weak_basis_mass_GeV"],md)
            for name,md in (("common_witness",1600),("high_mass_witness",2600))}
    data=rounded(data)
    target=HERE/"independent_results.json"
    if args.check:
        assert data==json.loads(target.read_text(encoding="utf-8"))
        print("1636 independent matrices, currents, mixing, phase-ray bounds and decay response: PASS")
    else:
        target.write_text(json.dumps(data,indent=2,allow_nan=False)+"\n",encoding="utf-8")
        print(json.dumps(data["common_witness"],ensure_ascii=False))


if __name__=="__main__":
    main()
