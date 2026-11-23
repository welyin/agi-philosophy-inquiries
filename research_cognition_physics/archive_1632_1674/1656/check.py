"""Two mixed CP-even poles in coherent mu-e conversion. All masses in GeV."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE / path).read_bytes())


def run():
    cfg = read("inputs.json")
    old = read(cfg["texture"])
    scalar = read(cfg["portal"])
    mh = read(cfg["texture_inputs"])["mh_GeV"]
    v, f, s = (scalar[k] for k in ("v_GeV", "f_GeV", "sin_theta"))
    c = math.sqrt(1 - s*s)
    mu = read(cfg["muon"])["mu_mass_GeV"]
    me = mu / read("../1480/inputs.json")["muonium"]["pole_ratio"][0]
    mass = np.array(cfg["nucleon_masses_GeV"])
    names = cfg["target_order"]
    targets = [cfg["targets"][n] for n in names]
    overlaps = np.array([t["overlaps"] for t in targets])
    caps = np.array([t["limit90"] for t in targets])
    capture = np.array([t["capture_per_s"] for t in targets]) * cfg["hbar_GeV_s"]
    # b is the cap on the norm of (m_h^2 sum_N S_N C_{chi,N}) in each target.
    cap_amp = mh**2 * np.sqrt(caps * capture / (4 * mu**5))
    sin2R = old["sin_2theta_down_right"]
    q, z = old["quark_q_s"]
    b = math.hypot(q*math.sqrt(1-z*z), z)
    a, cd = q/b, (1-q*q)*z*math.sqrt(1-z*z)/b
    _, _, vh = np.linalg.svd([[0, b], [a, cd]])
    nd, ns = 2 * vh[::-1, 0]**2
    np.testing.assert_allclose([nd+ns, nd*ns], [2, sin2R**2])
    vals = cfg["quark_fractions"]
    means = np.array([row[0] for row in vals.values()])
    widths = np.array([2*math.hypot(*row[1:]) for row in vals.values()])

    def sources(form):
        up, down = form[[0, 2]], form[[1, 3]]
        light = up + down + form[4]
        heavy_one = 2*(1-light)/27
        ordinary = mass/v * (light+3*heavy_one)
        radial_mass = mass * (nd*down+ns*form[4])
        return ordinary, radial_mass, mass*up, mass*heavy_one

    D, T, U, B = sources(means)
    points = {"degenerate_limit": 1.0, "decoupling_limit": 0.0}
    points.update({str(m): (mh/m)**2 for m in cfg["mH_tests_GeV"]})
    factors = np.sqrt(mu**2+me**2)
    response = {}
    bounds = {}
    cms = old["CMS_reference_norm_from_exact_BR_formula"]
    for label, x in points.items():
        weight = s*s+c*c*x
        coeff = s*c*(1-x)
        proton_neutron = coeff*D + weight*T/f
        support = overlaps @ proton_neutron
        h = overlaps @ (s*(c*D+s*T/f))
        H = overlaps @ (x*c*(-s*D+c*T/f))
        np.testing.assert_allclose(h+H, support, rtol=1e-14)
        ycap = s*cap_amp/support
        rates, fmin = {}, {}
        for key, branch in old["branches"].items():
            se = branch["sin_2theta_e"]
            lf = se*factors/f
            np.testing.assert_allclose(s*lf, np.linalg.norm(branch["h_e_mu_and_mu_e_abs_right"]))
            rates[key] = (caps*(lf*support/cap_amp)**2).tolist()
            # D/f + T/f^2, solve the positive root in f.
            aa = se*factors*coeff*(overlaps @ D)/cap_amp
            bb = se*factors*weight*(overlaps @ T)/cap_amp
            fmin[key] = ((aa+np.sqrt(aa**2+4*bb))/2).tolist()
            for j, scale in enumerate(fmin[key]):
                lhs = se*factors*(coeff*(overlaps@D)[j]/scale+weight*(overlaps@T)[j]/scale**2)
                assert math.isclose(lhs, cap_amp[j], rel_tol=1e-12)
        response[label] = {"x":x,"h_contribution":h.tolist(),"H_contribution":H.tolist(),
                           "Yh_norm_caps":ycap.tolist(),"CMS_to_caps":(cms/ycap).tolist(),
                           "rates_at_ref":rates,"conditional_min_f_GeV":fmin}
        # One open heavy response n_c+n_b+n_t; eliminate it by a two-target determinant.
        extra_q = weight/f*(overlaps@B)
        det = support[0]*extra_q[1]-support[1]*extra_q[0]
        ellmax = (cap_amp[0]*abs(extra_q[1])+cap_amp[1]*abs(extra_q[0]))/abs(det)
        # There is a common endpoint of the two target intervals at this norm.
        sols = [tuple(sorted(((-bd/ellmax-A)/bq, (bd/ellmax-A)/bq)))
                for bd, A, bq in zip(cap_amp, support, extra_q)]
        assert math.isclose(max(lo for lo, _ in sols), min(hi for _, hi in sols), abs_tol=1e-11)
        both = np.column_stack((U,B))*weight/f
        simultaneous = np.linalg.solve(both,-proton_neutron)
        np.testing.assert_allclose(both@simultaneous+proton_neutron, [0,0], atol=2e-18)
        bounds[label]={"profile_heavy_only_yh_cap":s*ellmax,
                       "profile_saturating_sum_nQ":max(lo for lo, _ in sols),
                       "both_zero_nu_sum_nQ":simultaneous.tolist()}
    # Endpoints in x and in all linear form factors suffice because all products are multilinear.
    box = []
    rate_max = {}
    stress = cfg["extra_amplitude_relative_stress"]
    for signs in itertools.product((-1,1), repeat=5):
        d,t,_,_ = sources(means+widths*np.array(signs))
        for x in (0.,1.):
            support = overlaps @ (s*c*(1-x)*d+(s*s+c*c*x)*t/f)
            assert np.all(support>0)
            box.extend([(s*cap_amp/(support*k)).tolist() for k in (1-stress,1+stress)])
            for key, branch in old["branches"].items():
                ell = max(branch["sin_2theta_e_interval"])*factors/f
                rates = caps*(ell*support*(1+stress)/cap_amp)**2
                rate_max[key] = np.maximum(rate_max.get(key, np.zeros(2)), rates)
    angle = math.degrees(math.acos(np.dot(*overlaps)/np.prod(np.linalg.norm(overlaps,axis=1))))
    ratio = (np.array(targets[1]["KKO_comparison"])@D)/(overlaps[1]@D)
    return {"order":names,"n_down_strange":[nd,ns],"sources_h0":D.tolist(),
            "sources_phi":(T/f).tolist(),"cap_amp":cap_amp.tolist(),"target_angle_deg":angle,
            "approx_q2_over_mh2":(mu/mh)**2,"Ti_old_to_final_Ycap":math.sqrt(targets[1]["older_limit90"]/caps[1]),
            "Ti_KKO_to_2025_ordinary_amplitude_ratio":ratio,"reference":response,
            "Yh_norm_caps_box_lo_hi":[np.min(box,axis=0).tolist(),np.max(box,axis=0).tolist()],
            "box_max_rates":{k:r.tolist() for k,r in rate_max.items()},
            "profiles_central":bounds}


def tidy(node):
    if isinstance(node, dict):
        return {k:tidy(v) for k,v in node.items()}
    if isinstance(node, list):
        return list(map(tidy,node))
    if isinstance(node,(float,np.floating)):
        return float(format(node,".11g"))
    return node


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    flags=parser.add_mutually_exclusive_group()
    flags.add_argument("--check",action="store_true")
    flags.add_argument("--write-new",action="store_true")
    args=parser.parse_args()
    out=tidy(run())
    if args.check:
        assert out == read("results.json")
    elif args.write_new:
        with (HERE/"results.json").open("x",encoding="utf8") as handle:
            json.dump(out,handle,indent=2)
            handle.write("\n")
    else:
        print(json.dumps(out,indent=2))
    print("1656 PASS: coherent scalar poles, per-target normalization and profiled cancellations")
