"""Bound the product of three pencil roots using traces and physical CKM lines."""
import argparse
import cmath
import json
import math
import runpy
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def integrate(flow, report, step):
    y = np.zeros((4, 3, 4), complex)
    y[0, :, :3] = np.diag(flow["central"]("up"))
    v0 = flow["low_ckm"](report["gamma_MZ_deg"], flow["angles"]("direct"))
    y[1, :, :3] = v0 @ np.diag(flow["central"]("down"))
    y[2, :, :3] = unpack(report["charged_left_MZ"]) @ np.diag(flow["central"]("charged_leptons"))
    value = 0.
    checks = {}
    def snapshots():
        singular = [np.linalg.svd(a[:,:3], compute_uv=False) for a in y[:3]]
        ratio = np.prod(singular[2])/np.prod(singular[1])
        u = np.linalg.eigh(y[0] @ y[0].conj().T)[1][:,0]
        d = np.linalg.eigh(y[1] @ y[1].conj().T)[1][:,0]
        return {"mean_root_modulus": ratio**(1/3), "abs_Vud": abs(u.conj() @ d)}
    def stage(a, g, susy):
        u,d,e,n = np.sum(np.abs(a)**2, axis=(1,2))
        g1, _, g3 = g
        if susy:
            s = 3*(e-d)+n-u+16*g3-4*g1
        else:
            s = 1.5*(e-n-d+u)+24*g3-6*g1
        return flow["beta"](a,g,susy), s
    for lo, susy, dt, g in flow["segments"](step):
        if lo==flow["MZ"]:
            checks["MZ"] = snapshots()
        if lo==flow["MS"]:
            checks["soft_below"] = snapshots()
            tan = flow["TAN"]
            y /= np.array([tan,1,1,tan])[:,None,None]/math.hypot(tan,1)
            checks["soft_above"] = snapshots()
        for n, pole in enumerate(flow["fields"]["neutral_poles_GeV"]):
            if lo == pole:
                y[3,:,n] = flow["physics"]["Y0"][:,n]
        for j in range(0,len(g)-1,2):
            k1,v1 = stage(y, g[j], susy)
            k2,v2 = stage(y+dt*k1/2, g[j+1], susy)
            k3,v3 = stage(y+dt*k2/2, g[j+1], susy)
            k4,v4 = stage(y+dt*k3, g[j+2], susy)
            y += dt*(k1+2*k2+2*k3+k4)/6
            value += dt*(v1+2*v2+2*v3+v4)/6
    checks["10PeV"] = snapshots()
    np.testing.assert_allclose(value, 3*math.log(checks["10PeV"]["mean_root_modulus"] /
                                                 checks["MZ"]["mean_root_modulus"]), atol=2e-9)
    return checks, value


def unpack(a):
    return np.array(a["real"]) + 1j*np.array(a["imag"])


def certify(flow):
    """Frobenius-norm barriers, then Fubini-Study lengths of u and d lines."""
    log_unit = 16*math.pi**2
    T = math.log(flow["MS"]/flow["MZ"])/log_unit
    S = math.log(flow["END"]/flow["MS"])/log_unit
    # g columns are squared gauge couplings. For each threshold-free segment
    # they are monotone analytically, so endpoints suffice for these bounds.
    endpoints = np.vstack([group[-1][[0,-1]] for group in flow["segments"](.1)])
    g1min,g2min,g3min = endpoints.min(axis=0)
    g1max = max(endpoints[:,0])
    # First exit argument, starting below all the indicated barriers.
    un,dn,en = (np.linalg.norm(flow["central"](key)) for key in ("up","down","charged_leptons"))
    nn = np.linalg.norm(flow["physics"]["Y0"])
    U,D,E,N = 1., .017, .025, .24
    trace = 3*U**2+3*D**2+E**2+N**2
    assert un<U and dn<D
    assert trace+1.5*U**2-(.85*g1min+2.25*g2min+8*g3min)<0
    assert trace+1.5*D**2-(.25*g1min+2.25*g2min+8*g3min)<0
    e_soft = en*math.exp(T*(trace+1.5*E**2))
    n_soft = nn*math.exp(T*(trace+1.5*N**2))
    assert e_soft<E and n_soft<N
    su,sd,se,sn = 1.6,.17,.15,.35
    tan = flow["TAN"]
    c, s = 1/math.hypot(tan,1), tan/math.hypot(tan,1)
    u0 = U/s
    a = sd**2+sn**2
    # Comparison equation Y'=6Y^3+aY. Its squared norm has an elementary integral.
    discr = 1-(6/a)*u0**2*math.expm1(2*a*S)
    u_bound = u0*math.exp(a*S)/math.sqrt(discr)
    u2int = -math.log(discr)/12
    d_bound = D/c*math.exp(S*(6*sd**2+su**2+se**2))
    e_bound = e_soft/c*math.exp(S*(4*se**2+sn**2+3*sd**2))
    n_bound = n_soft/s*math.exp(S*(4*sn**2+se**2+3*su**2))
    assert max(u_bound/su, d_bound/sd, e_bound/se, n_bound/sn)<1
    low_rate = 24*g3min-6*g1max-1.5*(N**2+D**2)
    high_rate = 16*g3min-4*g1max-3*sd**2-su**2
    assert min(low_rate,high_rate)>0
    # D12 and U12 are mass (not squared-mass) ratios. The logarithmic
    # drift of each is bounded by the spectral spread of its left generator.
    u12, d12 = (flow["central"](key)[0]/flow["central"](key)[1] for key in ("up","down"))
    low_ratio = 1.5*T*(U**2+D**2)
    u12 *= math.exp(low_ratio+3*u2int+S*sd**2)
    d12 *= math.exp(low_ratio+u2int+3*S*sd**2)
    ku,kd = ((1+r*r)/(1-r*r) for r in (u12,d12))
    delta = .5*(1.5*T*(ku*D**2+kd*U**2) + ku*S*sd**2+kd*u2int)
    return {"sm_and_susy_log_intervals_over_16pi2": [T,S],
            "upper_yukawa_norms_u_d_e_n_susy": [u_bound,d_bound,e_bound,n_bound],
            "integral_susy_top_Frobenius_squared_upper": u2int,
            "positive_log_determinant_ratio_derivative_floors_sm_susy": [low_rate,high_rate],
            "mass_ratio_bounds_u_over_c_and_d_over_s": [u12, d12],
            "angle_speed_upper_sm_susy": [.75*(ku*D**2+kd*U**2), .5*(ku*sd**2+kd*su**2)],
            "physical_u_d_Fubini_Study_angle_change_upper_rad":delta}


def run():
    old = runpy.run_path(str(HERE/"../1715/check.py"))
    flow, rep = old["charged_background"](.1)
    checks, val = integrate(flow, rep, .10)
    finer, vf = integrate(flow, rep, .05)
    for key in checks:
        np.testing.assert_allclose(list(checks[key].values()),list(finer[key].values()),atol=2e-10)
    cert = certify(flow)
    d,e = (unpack(rep[k]) for k in ("down_texture_10PeV","charged_texture_10PeV"))
    e = e.T
    prod = np.linalg.det(e)/np.linalg.det(d)
    mu = abs(prod)**(1/3)
    np.testing.assert_allclose(mu, checks["10PeV"]["mean_root_modulus"],atol=2e-9)
    vmin = math.cos(math.acos(checks["MZ"]["abs_Vud"])
                    +cert["physical_u_d_Fubini_Study_angle_change_upper_rad"])
    # The ratio V/mean and even the difference V-mean decrease upwards
    # throughout the frozen flow. This locates the true worst scale, not
    # just an envelope assembled from values at unrelated scales.
    log_v_max = math.sqrt(1-vmin*vmin)/vmin*np.array(cert["angle_speed_upper_sm_susy"])
    log_mean_min = np.array(cert["positive_log_determinant_ratio_derivative_floors_sm_susy"])/3
    assert all(log_v_max < log_mean_min)
    assert all(log_v_max < checks["MZ"]["mean_root_modulus"]*log_mean_min)
    gain = checks["10PeV"]["abs_Vud"]/mu
    assert gain>1
    return {"snapshots":finer,"quadrature_log_det_gain":vf,"certificate":cert,
            "physical_Vud_lower_at_every_scale":vmin,
            "absolute_log_Vud_derivative_upper_sm_susy":log_v_max.tolist(),
            "one_root_must_lie_below_physical_Vud_by_at_least": checks["10PeV"]["abs_Vud"]-mu,
            "all_three_moduli_equal_Vud_necessary_net_finite_det_gain":gain**3,
            "each_of_six_mass_factors_symmetric_relative_width_needed":(gain-1)/(gain+1),
            "anchored_product_phase_deg":math.degrees(cmath.phase(prod)),
            "one_root_angle_separation_from_axis_zero_floor_deg":abs(math.degrees(cmath.phase(prod)))/3,
            "one_root_angle_separation_from_axis_minus30_floor_deg":
                abs(math.degrees(cmath.phase(prod*cmath.exp(.5j*math.pi))))/3}


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args=parser.parse_args()
    result=json.loads(json.dumps(run()),parse_float=lambda v:float(f"{float(v):.11g}"))
    path=HERE/"independent_results.json"
    if args.check:
        assert result==json.loads(path.read_bytes())
        print("1729 certified the determinant and CKM-line bounds")
    else:
        path.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf-8")
