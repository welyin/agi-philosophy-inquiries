"""Gauge-referenced common phase and determinant matching on the same copies."""
import argparse
import json
import math
from pathlib import Path
import runpy

import numpy as np


HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE/path).read_bytes())


CFG = read("inputs.json")
VAC = read(CFG["vacuum"])
LOCK = VAC["real_X_over_cutoff"]*VAC["scalar_magnitude_over_cutoff"]**4


def pair(z):
    return [float(z.real), float(z.imag)]


def deg(z):
    return float(np.angle(z, deg=True))


def gram(z, K):
    rho, ct = abs(z), z.real/abs(z)
    discriminant = math.sqrt((1+rho*rho)**2-4*(rho*ct)**2)
    kmin = math.sqrt((1+rho*rho+discriminant)/2)
    pbest = 2*rho*ct/(1+rho*rho+discriminant)
    def floor(p):
        return math.sqrt((1+rho*rho-2*rho*p*ct)/(1-p*p))
    np.testing.assert_allclose(floor(pbest), kmin, rtol=1e-12)
    # Given t_color = 1/K (the worst allowed), the normalized gauge overlap
    # lies on or inside a disk with positive center in the source-aligned gauge.
    radius = math.sqrt((K*K-1)*(K*K-rho*rho))/(K*K)
    center = rho/(K*K)
    assert radius < center
    return {"gain":rho, "wanted_phase_deg":deg(z), "original_source_budget_over_yb":K,
            "max_phase_deg_at_original_budget":math.degrees(math.asin(radius/center)),
            "gauge_modulus_if_original_budget_and_fixed_gain":[center-radius,center+radius],
            "least_source_budget_over_yb_any_nonzero_gauge":kmin,
            "gauge_modulus_at_least_budget":pbest,
            "source_floors_for_named_gauge_amplitudes":[floor(p) for p in CFG["gauge_amplitude_comparisons"]]}


def reference(x,x6_weight):
    return LOCK if x6_weight is None else x6_weight*x**6


def rows(q,x,s,weight=None):
    a = reference(x,weight)
    current = x*s**4
    return np.array([-q*a,a],complex),np.array([current,q*a],complex)


def state(row):
    ray = np.array([row[1],-row[0]])/np.linalg.norm(row)
    return ray*(ray[0].conj()/abs(ray[0]))


def vacuum_readout(q,x,s,weight=None):
    rc,rw = rows(q,x,s,weight)
    u,v = state(rc),state(rw)
    g = np.vdot(u,v)
    w = v[0]/u[0]
    z = w*g.conj()/abs(g)
    # Change each light phase independently and test the measured closed product.
    ur,vr = u*np.exp(.62j),v*np.exp(-.87j)
    gg = np.vdot(ur,vr)
    np.testing.assert_allclose(vr[0]/ur[0]*gg.conj()/abs(gg), z, atol=2e-15)
    return z, g, np.linalg.norm(rw)/np.linalg.norm(rc)


def experiment(q,x,s,weight=None):
    z,g,poles = vacuum_readout(q,x,s,weight)
    k = abs(x*s**4)/reference(x,weight)
    phase = -x*s**4/abs(x*s**4)
    predicted_phase = -np.angle(1+k*phase)
    np.testing.assert_allclose(np.angle(z),predicted_phase,atol=1e-14)
    np.testing.assert_allclose(abs(z)*poles,q,atol=1e-14)
    eps = 1e-5
    readings = [vacuum_readout(q,x,s*np.exp(t),weight) for t in (-eps,eps)]
    slope = np.angle(readings[1][0]/readings[0][0])/(2*eps)
    mass_slope = math.log(readings[1][2]/readings[0][2])/(2*eps)
    denominator = abs(1+k*phase)**2
    np.testing.assert_allclose(slope,-4*k*phase.imag/denominator,rtol=1e-8)
    np.testing.assert_allclose(mass_slope,4*k*k/(q*q+k*k),rtol=1e-8)
    xread = [vacuum_readout(q,x*np.exp(t),s,weight)[0] for t in (-eps,eps)]
    xslope = np.angle(xread[1]/xread[0])/(2*eps)
    desired = (-1 if weight is None else 5)*k*phase.imag/denominator
    np.testing.assert_allclose(xslope,desired,rtol=1e-8)
    return {"z_physical":pair(z),"phase_deg":deg(z),"gauge_conversion":pair(g),
            "light_to_light_gauge_probability":abs(g)**2,"heavy_weak_over_color_mass":poles,
            "color_heavy_over_light_source_gain":q,"source_norm_over_yb":math.hypot(1,q),
            "d_phase_d_log_abs_s_rad":slope,"d_log_heavy_ratio_d_log_abs_s":mass_slope,
            "phase_slope_over_heavy_ratio_slope":slope/mass_slope,
            "d_phase_d_log_x_rad":xslope}


def color_matching(q,x,s):
    old = runpy.run_path(str(HERE/CFG["charged_reconstruction"]))
    _,rep = old["charged_background"](.1)
    d = old["unpack"](rep["down_texture_10PeV"])
    phi = np.angle(s)
    p10 = np.diag(np.exp(1j*phi*np.array([0,0,1])))
    pf = np.diag(np.exp(1j*phi*np.array([2,0,1])))
    bare = math.hypot(1,q)*p10@d@pf
    up = p10@np.diag(rep["up_texture_diag_10PeV"])@p10
    current = x*s**4
    up = up@np.diag([current/abs(current),(current/abs(current))**2,1.])
    old_triplets = np.array([[0,s**3],[-s**3,1.]])
    rc,rw = rows(q,x,s)
    eye = np.eye(3)
    vh = CFG["vev_over_parent_mass_for_direct_test"]
    def master(mix,second):
        row = rc*(1-mix)+rw*mix
        return np.block([[vh*bare,vh*second*bare],[row[0]*eye,row[1]*eye]])
    def theta(mix,second):
        matrix = master(mix,second)
        return np.angle(np.linalg.slogdet(matrix)[0]*np.linalg.slogdet(up)[0]*np.linalg.slogdet(old_triplets)[0])
    base = master(0.,0.)
    colors = state(rc)
    rotation = np.kron(np.column_stack([colors,rc.conj()/np.linalg.norm(rc)]),eye)
    rotated = base@rotation
    np.testing.assert_allclose(rotated[:3,:3],vh*p10@d@pf,atol=1e-12)
    np.testing.assert_allclose(rotated[3:,:3],0,atol=1e-14)
    np.testing.assert_allclose(rotated[3:,3:],np.linalg.norm(rc)*eye,atol=1e-14)
    np.testing.assert_allclose(np.linalg.det(rotation),1,atol=1e-14)
    # The three light colors already existed. The determinant includes each Dirac
    # flavor once for the QCD anomaly, not a separate extra power for each color.
    e,l = (CFG["pollution_example"][t] for t in ("color_mixture","second_Yukawa"))
    a = reference(x,None)
    cycle = (1-e)*(1+l*q)+e*q-e*l*current/a
    np.testing.assert_allclose(theta(e,l),3*np.angle(cycle),atol=2e-13)
    for mix,second in ((0,0),(e,0),(0,l)):
        assert abs(theta(mix,second))<1e-13
    budget = CFG["theta_working_limit_rad_rounded"]
    return {"old_light_plus_new_triangular_colors_phase":theta(0,0),
            "color_leak_only_phase":theta(e,0),"second_Yukawa_only_phase":theta(0,l),
            "combined_phase_rad":theta(e,l),"analytic_phase_rad":3*np.angle(cycle),
            "allowed_abs_product_over_real_cycle_for_single_theta_source":
                math.tan(budget/3)/abs((current/a).imag),
            "new_threshold_Dirac_masses_in_symmetric_phase_per_family":np.linalg.norm(rc),
            "down_and_old_compensators_at_reference_reconstruct":True}


def calculate():
    targets = read(CFG["charged_targets"])
    zt = complex(*targets["exact_and_conditional_implementations"]["second_complex_root"]["z"])
    optimum = complex(*read(CFG["continuous_target"])["tagged"]["complex"]["z"])
    old = read(CFG["vacuum"])
    x = old["real_X_over_cutoff"]
    s = old["scalar_magnitude_over_cutoff"]*np.exp(1j*np.pi*old["scalar_angle_over_pi"])
    radius = targets["five_path_tree_match"]["R"]
    K = 1/math.sqrt(1-radius**-2)
    q = abs(zt)
    cases = {"same_real_entries_q1":experiment(1,x,s),
             "q_calibrated_to_two_path_root":experiment(q,x,s),
             "different_radial_same_vacuum_2x6":experiment(q,x,s,2.0),
             "alternative_unit_coefficient_in_x6":experiment(q,x,s,1.0)}
    np.testing.assert_allclose(cases["q_calibrated_to_two_path_root"]["z_physical"],pair(zt),atol=2e-11)
    return {"targets":{"two_path_root":gram(zt,K),"continuously_optimized":gram(optimum,K)},
            "witness":cases,"color":color_matching(q,x,s),
            "phase_sourced_by_existing_I":pair(x*s**4),"real_mass_reference_a":LOCK}


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args = parser.parse_args()
    result = calculate()
    if args.check:
        assert read("results.json")==result
        print("1728 gauge phase bound, three-family singlet-adjoint matching and colored determinants verified")
    else:
        (HERE/"results.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf-8")
