"""Test all singly-inserted Z12 two-source sectors on the labelled charged background."""
from pathlib import Path
import argparse
import json
import math
import runpy

import numpy as np

HERE = Path(__file__).resolve().parent
CFG = json.loads((HERE/"inputs.json").read_bytes())


def read(filename):
    return json.loads((HERE/filename).read_bytes())


def masks():
    old = runpy.run_path(str(HERE/CFG["pair_enumeration"]))
    pairs = old["check"]()["pairs_and_texture"]
    allowed = lambda h: old["allowed"](h,1)
    return {str(tuple(row["pair"])):{
        "D": [[bool(allowed(l+n)) for n in row["pair"]] for l in CFG["family_charges"]],
        "M": [[bool(allowed(n+z)) for n in row["pair"]] for z in row["pair"]]}
        for row in pairs}


def bounds(W):
    data = read(CFG["data"])
    r, R = data["assignments"]["NO"]["middle_over_heaviest_mass_required"]
    rI, RI = data["assignments"]["IO"]["middle_over_heaviest_mass_required"]
    x, X = data["experimental_rectangles"]["isolated_weight_roots"][0]
    s, S = data["experimental_rectangles"]["solar_share"]
    a, A = read(CFG["atmospheric"])["joint_windows"][2]
    left, v, right = np.linalg.svd(W)
    Q = left@right
    # Compare the normalized inverse positive polar factor with the identity.
    gain = 2/(1/min(v)+1/max(v))
    g = max(abs(gain/v-1))
    eps = g*(2+g)
    distance = lambda gap: 2*math.sin(.5*math.asin(eps/(gap-eps)))
    dm, dh = distance(r), distance(1-R)
    rho, RH = (r-eps)/(1+eps),(R+eps)/(1-eps)
    solar_amp = math.sqrt((1-X)*s)-dm
    atm_amp = math.sqrt(X)+dh

    tagged = np.linalg.inv(W)[:,0]
    tagged /= np.linalg.norm(tagged)
    e, mu, tau = abs(tagged)
    # IO: the null isotropic condition is a union of two complex hyperplanes.
    # Independently of which member of the close pair owns the larger electron weight.
    io_weight_gap = min(math.sqrt(rI)*math.sqrt((1-X)*(1-S))-math.sqrt((1-x)*S),
                        math.sqrt((1-X)*(1-S))-math.sqrt(RI)*math.sqrt((1-x)*S))
    io_sine_min = io_weight_gap/math.sqrt(1+RI)
    physical_tilt = math.hypot(mu,tau)
    assert io_sine_min > physical_tilt

    zero_prob_if_first_absent = (np.linalg.norm(np.linalg.inv(W)[0,1:])*max(v))**2
    assert zero_prob_if_first_absent < x

    one_row = abs(Q[1,0])*math.sqrt(1+RH)
    one_row_needed = math.sqrt(rho)*solar_amp-atm_amp
    assert one_row_needed > one_row
    P = np.exp(1j*math.pi*CFG["phi_over_pi"]*np.array(CFG["family_charges"]))
    leg = P*Q[:,0]
    leg *= np.exp(-1j*np.angle(leg[0]))
    phase_allowed = math.sqrt(1+RH)*np.linalg.norm(leg[1:].real)
    phase_needed = math.sqrt(rho)*math.sqrt(solar_amp**2-np.linalg.norm(leg.imag)**2)-atm_amp
    assert phase_needed > phase_allowed

    # Bound cos(delta) from the exact transported isotropic line, irrespective of Wilson phases.
    hee = e*e*(r*s*(1-X)-X)
    amax = R*math.sqrt(S*(1-S)*(1-a))
    bmax = math.sqrt(X*A)*(1+R*S)
    cross_tau = R*math.sqrt(S*(1-S)*A)+math.sqrt(X*(1-a))*(1+R*S)
    transverse = (mu*math.sqrt(A*(1-x))+tau*math.sqrt((1-a)*(1-x)))**2
    wanted = (hee-2*e*tau*cross_tau-transverse)/(2*e*mu)
    num = amax**2+bmax**2-wanted**2/(1-x)
    assert wanted>0 and num<0
    cos_max = num/(2*amax*bmax)

    # After unobservable column phases the reference real eigenaxes acquire
    # only the imaginary parts of the two measured charged columns of Q.
    q = Q.copy()
    for k in (0,1):
        q[:,k] *= np.exp(-1j*np.angle(q[np.argmax(abs(q[:,k])),k]))
    ei, ui = (np.linalg.norm(q[:,i].imag) for i in (0,1))
    # Quartet e3,mu3,e2,mu2, using real parts of its four physical entries as reference.
    mods = [math.sqrt(X),math.sqrt(A*(1-x)),math.sqrt(S*(1-x)),
            min(1., math.sqrt((1-s)*(1-a))+math.sqrt(S*A*X))]
    im_max = [ei+dh,ui+dh,ei+dm,ui+dm]
    jmax = sum(im_max[i]*math.prod(mods[:i]+mods[i+1:]) for i in range(4))
    pref = (1-X)*math.sqrt(x*s*(1-s)*min(a*(1-a),A*(1-A)))
    sin_max = jmax/pref
    assert sin_max < 1 and cos_max<0
    return {"polar_singulars":v.tolist(),"phase_insensitive_inverse_congruence_error":eps,
            "middle_heavy_vector_distances":[dm,dh],"transported_first_axis_moduli":abs(tagged).tolist(),
            "IO_required_tilt_sine_min":io_sine_min,"first_axis_tilt_sine":physical_tilt,
            "massless_e_weight_cap_for_two_sparse_pairs":zero_prob_if_first_absent,
            "zero_e2_block_08_needed_and_allowed":[one_row_needed,one_row],
            "physical_extra_P_comparison_needed_and_allowed":[phase_needed,phase_allowed],
            "all_phases_cos_delta_upper":cos_max,"common_real_all_coefficients_abs_J_upper":jmax,
            "common_real_max_distance_from_pi_degrees":math.degrees(math.asin(sin_max))}


def observable(W, q):
    theta, psi, ratio = q
    c, h = math.cos(theta),math.sin(theta)
    cs,sn = math.cos(psi),math.sin(psi)
    a, b = np.array([h,-c*cs,-c*sn]),np.array([0.,-sn,cs])
    m = (a+math.sqrt(ratio)*b)/math.sqrt(1+ratio)
    H = (math.sqrt(ratio)*a-b)/math.sqrt(1+ratio)
    R = ratio*np.outer(m,m)-np.outer(H,H)
    C = W.T@R@W
    _, sv, vh = np.linalg.svd(C)
    U = vh.conj().T[:,::-1]
    x = abs(U[0,2])**2
    s, t = abs(U[0,1])**2/(1-x),abs(U[1,2])**2/(1-x)
    j = np.imag(U[0,0]*U[1,1]*U[0,1].conj()*U[1,0].conj())
    st = j/((1-x)*math.sqrt(s*(1-s)*t*(1-t)*x))
    ct = (abs(U[1,0])**2-s*(1-t)-(1-s)*t*x)/(2*math.sqrt(s*(1-s)*t*(1-t)*x))
    return np.array([s,t,sv[1]/sv[0]]), np.array([4*x*(1-x),math.degrees(math.atan2(st,ct)),j]), (m,H,R,C,sv)


def witness(W):
    cfg = CFG["witness"]
    target = np.array([cfg["sin2solar"],cfg["sin2atm"], cfg["middle_over_high"]])
    q = np.array(cfg["seed"])
    for _ in range(12):
        y = observable(W,q)[0]
        jac = np.column_stack([(observable(W,q+step)[0]-observable(W,q-step)[0])/2e-5
                               for step in 1e-5*np.eye(3)])
        dx = np.linalg.solve(jac,y-target)
        q -= dx
        if max(abs(dx))<1e-12:
            break
    y, out, (m,H,R,C,sv) = observable(W,q)
    np.testing.assert_allclose(y,target,atol=1e-10)
    scale = math.sqrt(cfg["solar_delta_eV2"])/sv[1]
    assert .0803<=out[0]<=.0899
    d31 = (scale*sv[0])**2
    assert .002459<=d31<=.002563 and .49<=y[1]<=.63 and .29048<=y[0]<=.31672

    # Off-diagonal heavy matrix plus an allowed unit (0,0) mass.
    u = (math.sqrt(q[2])*m+H)/math.sqrt(2)
    v = (math.sqrt(q[2])*m-H)/math.sqrt(2)
    u += .5*v
    assert abs(v[0])<1e-15
    ell = np.array(CFG["family_charges"])
    phase = np.exp(1j*math.pi*CFG["phi_over_pi"])
    p = np.diag(phase**ell)
    D = p@np.column_stack((u,phase**2*v))
    heavy = np.array([[1,phase**2],[phase**2,0]])
    recovered = np.linalg.solve(p,D)@np.linalg.solve(heavy,np.linalg.solve(p,D).T)
    np.testing.assert_allclose(R,recovered,atol=1e-13)
    return {"high_angles_ratio":q.tolist(), "calibration_s12_s23_mass_ratio":y.tolist(),
            "delta21_eV2":cfg["solar_delta_eV2"],"delta31_eV2":d31,
            "reactor_fast_amplitude":out[0],"dirac_delta_degrees":out[1],"jarlskog":out[2],
            "light_exchange_mbb_meV":1000*scale*abs(C[0,0]),"common_mass_scale_eV":scale,
            "raw_real_D_coefficients_before_common_scale":[(u/CFG["svev"]**ell).tolist(),
                  (v/CFG["svev"]**(ell+2)).tolist()],
            "raw_M00_real_M02_A2_real": [1.,1/CFG["svev"]**2],
            "dimensionless_heavy_singular_values":np.linalg.svd(heavy,compute_uv=False).tolist()}


def calculate(step):
    reader = runpy.run_path(str(HERE/CFG["charged_frame"]))
    flow, rep = reader["charged_background"](step)
    down = reader["unpack"](rep["down_texture_10PeV"])
    phi = math.pi*CFG["phi_over_pi"]
    p10 = np.diag(np.exp(1j*phi*np.array([0,0,1])))
    p5 = np.diag(np.exp(1j*phi*np.array(CFG["family_charges"])))
    # Lift the five-real-edge charged texture back to the underlying CP basis.
    base = down.copy()
    long_edge = down[1,2].imag/math.sin(4*phi)
    base[1,2] -= long_edge*np.exp(4j*phi)
    assert np.max(abs(base.imag))<1e-14
    degrees = np.array([0,0,1])[:,None]+np.array(CFG["family_charges"])
    raw = base.real*np.exp(1j*phi*degrees)
    raw[1,2] += long_edge*np.exp(1j*phi*(degrees[1,2]+4))
    np.testing.assert_allclose(raw,p10 @ down @ p5,atol=1e-14)
    W, tr = reader["congruence"](flow,rep,step)
    bound = bounds(W)
    atlas = read(CFG["phase_literature"])
    rad = math.radians(bound["common_real_max_distance_from_pi_degrees"])
    arc = [-math.pi-rad,-math.pi+rad]
    mid,minus,plus = atlas["joint_2025"]["hpd1_center_minus_plus_in_pi"]
    joint = [(mid-minus)*math.pi,(mid+plus)*math.pi]
    t2k, minus, plus = atlas["t2k_updated"]["fc1_center_minus_plus_rad"]
    assert joint[0]<arc[0]<arc[1]<joint[1] and arc[1]<t2k-minus
    return {"transport":tr,"charged_gamma_deg":rep["gamma_MZ_deg"],"class_bounds":bound,
            "CP_holdouts":{"arc_rad_unwrapped":arc,"inside_joint_NO_HPD68":True,
                "distance_to_update_T2K_FC68_rad":t2k-minus-arc[1]},
            "two_source_witness_02":witness(W)}


def clean(x):
    if isinstance(x,dict): return {k:clean(v) for k,v in x.items()}
    if isinstance(x,list): return [clean(v) for v in x]
    if isinstance(x,(float,np.floating)): return float(f"{x:.10g}")
    return x


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--check",action="store_true")
    coarse,fine = [calculate(h) for h in CFG["steps"]]
    ch = abs(coarse["two_source_witness_02"]["dirac_delta_degrees"]-fine["two_source_witness_02"]["dirac_delta_degrees"])
    assert ch<2e-7
    result = clean({**fine,"half_step_delta_change_below_deg":2e-7,"allowed_fields":masks()})
    path = HERE/"results.json"
    if ap.parse_args().check:
        assert result==json.loads(path.read_bytes())
        print("1784 verified: six singlet charge pairs, transported axes, normal branch, CP holdout")
    else:
        path.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
