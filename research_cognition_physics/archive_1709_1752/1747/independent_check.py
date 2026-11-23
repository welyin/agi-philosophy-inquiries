"""Inverse spectral moments give all-orientation bounds for the physical dark line."""
import argparse
import itertools
import json
import math
from pathlib import Path
import runpy

import numpy as np


HERE = Path(__file__).resolve().parent
BASE = runpy.run_path(str(HERE / "../1744/check.py"))
OLD = json.loads((HERE / "../1746/independent_results.json").read_bytes())


def upper_root(fun):
    """A monotone upper envelope starts excluded and eventually meets the data."""
    assert fun(0.) < 0 and fun(1.) > 0
    lo, hi = 0., 1.
    for _ in range(60):
        mid = (lo + hi) / 2
        if fun(mid) < 0:
            lo = mid
        else:
            hi = mid
    return max(0., lo - 1e-9)


def min_mixing(gains, bounds):
    # Variance of an eigenvalue distribution is concave on the simplex.
    vertices = []
    caps = 1 - np.square(bounds)
    for i, j in itertools.permutations(range(3), 2):
        k = 3-i-j
        for x in [0., caps[i]]:
            for y in [0., caps[j]]:
                p = np.zeros(3)
                p[[i,j,k]] = [x,y,1-x-y]
                if all(p >= -1e-13) and all(p <= caps+1e-13):
                    v = sum(p[a]*p[b]*(gains[a]-gains[b])**2
                            for a,b in itertools.combinations(range(3),2))
                    vertices.append((v,p))
    return min(vertices, key=lambda row:row[0])


def solve(step):
    b = BASE["mother"](step)
    Y = b["yW"][2,:,:3]
    lc, d, _ = np.linalg.svd(Y)
    d = d[::-1]
    electron = lc[:,-1]
    pt = b["pt"]
    S = np.linalg.inv(b["S"] @ pt["low_E"])
    M = S.T @ (pt["C"] + b["nu"]["FROZEN"][0]*np.outer(pt["b"],pt["b"])) @ S
    source_basis, _, _ = np.linalg.svd(b["N"])
    q = source_basis[:,-1]
    mu = np.linalg.svd(M,compute_uv=False)
    source_residual = np.linalg.norm(M@q)/mu[1]
    assert source_residual < 1e-11
    qmass = np.linalg.svd(M)[2].conj().T[:,-1]
    assert np.linalg.norm(b["N"].conj().T @ q) < 1e-13
    moment = sum((d[0]/d[i])**2 * abs(np.vdot(lc[:,2-i],q))**2 for i in range(3))
    p0 = abs(np.vdot(electron,q))**2
    spectrum = np.array([b["cW"],b["fW"],1.])
    inverse2 = spectrum**-2
    vol = np.prod(spectrum)**(1/3)
    # Require both control and daughter to stay in the same registered effective range.
    # A weak-scale common charged factor up to 1.05 includes the tau calibrations.
    caps = [.02,.25]
    field = pt["flow"]
    t = math.log(b["Mweak"]/field["MZ"])/(16*math.pi**2)
    ta = math.log(b["Mweak"]/min(field["fields"]["neutral_poles_GeV"]))/(16*math.pi**2)
    assert max(np.max(g[:,:2]) for _,_,_,g in b["segments"]) < .6
    assert 1.05*d[-1]*math.exp((4.5*.6+1.5*caps[1]**2)*t)<caps[0]
    assert np.linalg.norm(b["yW"][3],2)*math.exp((2.7*.6+1.5*caps[0]**2)*ta)<caps[1]
    e2, n2 = caps[0]**2*t, caps[1]**2*ta
    charged_log = 3*(e2+n2)
    charged_limits = b["model"]["JOINT"]["LEPTONS"]
    exp_frac = 2*charged_limits[:,1]/charged_limits[:,0]
    control = BASE["daughter"](b,np.eye(3))
    ratios = 1+control["MZ_relative_e_mu_after_tau_calibration"]
    lo = (1-exp_frac[:2])/(1+exp_frac[2])/ratios*math.exp(-charged_log)
    hi = (1+exp_frac[:2])/(1-exp_frac[2])/ratios*math.exp(charged_log)
    e_lo, e_hi = vol*(lo[0]**2/hi[1])**(1/3), vol*(hi[0]**2/lo[1])**(1/3)
    e_mu_square = (d[0]/(spectrum.min()*d[1]))**2*math.exp(charged_log)
    electron_drift = .75*(1+e_mu_square)/(1-e_mu_square)*n2
    neutral_drift = .5*(1.5*e2 + .5*n2)
    # Outward pad is >100 times the measured mass/source line residual.
    drift = electron_drift + neutral_drift + 1e-8
    ang = b["nu"]["ANG"]
    reactor_zmax = (1-math.sqrt(1-ang[0,1]))/2
    observed_pmin = (1-ang[1,1])*(1-reactor_zmax)
    # S = e^† F^-2 e. The leading singular value of Y^-1 F^-1 lies between
    # sqrt(S)/d_e and the full Frobenius norm. q_h may rotate away from q.
    pad = 1e-10
    remainder = sum((d[0]/d[1:])**2)*inverse2.max()+pad
    def weight_lower(s):
        return max(p0-s-pad,0.)
    def electron_gain_upper(s):
        # h=c, nearest competitors' inverse squares >= min(inv2).
        return (inverse2.min()+(inverse2[0]-inverse2.min())*weight_lower(s)-pad)**-.5
    def electron_gain_lower(s):
        # h=1, nearest competitors' inverse squares <= max(inv2).
        return (inverse2.max()+(1-inverse2.max())*weight_lower(s)+remainder)**-.5
    H = inverse2[1]
    delta = max(abs(inverse2-H))
    def null_e_upper(s):
        # F^-2 q differs from H q by <= delta*s, and its expectation by <= delta*s^2.
        return e_hi**2*(H*math.sqrt(moment)+delta*s+pad)**2/(H-delta*s*s-pad) + drift
    excluded_sin = np.array([upper_root(lambda s:electron_gain_upper(s)-e_lo),
                             upper_root(lambda s:null_e_upper(s)-observed_pmin),
                             upper_root(lambda s:e_hi-electron_gain_lower(s))])
    var, at = min_mixing(spectrum,excluded_sin)
    W = OLD["fine"]
    F = np.array(W["F_real"])+1j*np.array(W["F_imag"])
    f,axes = np.linalg.eigh(F)
    np.testing.assert_allclose(f,spectrum,atol=2e-10)
    ws = np.sqrt(np.maximum(0,1-abs(q.conj()@axes)**2))
    assert np.all(ws>excluded_sin)
    coupling = np.linalg.norm((np.eye(3)-np.outer(q,q.conj()))@F@q)
    assert coupling**2>var
    # Check the measured electron-column identity with a general noncommuting F.
    qd = np.linalg.solve(F,qmass)
    qd /= np.linalg.norm(qd)
    L,ds,_ = np.linalg.svd(F@Y)
    observed_sum = sum((d[0]/ds[k])**2*abs(np.vdot(L[:,k],qd))**2 for k in range(3))
    qq = np.linalg.solve(F,np.linalg.solve(F,qmass))
    expected = np.linalg.norm(d[0]*np.linalg.solve(Y,qq))**2/np.linalg.norm(np.linalg.solve(F,qmass))**2
    assert abs(observed_sum-expected)<1e-11
    return dict(weak_scale_GeV=b["Mweak"],dark_spectrum=spectrum,p_e_mother=p0,full_inverse_moment=moment,
                source_to_mass_null_sin_upper=source_residual,geometric_mean=vol,
                charged_gain_ratio_intervals=[lo,hi],electron_gain_required=[e_lo,e_hi],
                exact_dark_c_electron_gain_upper=electron_gain_upper(0),
                exact_dark_one_electron_gain_lower=electron_gain_lower(0),
                angle_drift_bound=dict(charged=electron_drift,neutral=neutral_drift,plus_pad=drift),
                dark_half_angle_p_e_MZ_upper=null_e_upper(0),experimental_p_e1_low=observed_pmin,
                all_three_eigenline_minimum_sines=excluded_sin,
                all_three_eigenline_minimum_angles_deg=np.degrees(np.arcsin(excluded_sin)),
                minimum_norm_source_projection_F_dark=math.sqrt(var)*(1-1e-10),
                relaxing_probability_vertex=at,
                competing_1746_actual_sines=ws,competing_1746_projected_amplitude=coupling,
                competing_1746_p_e1=W["actual_mass_zero_column"][0],competing_1746_clocks=[W["dm31"],W["dm32"]],
                noncommuting_inverse_moment_check=abs(observed_sum-expected))


def clean(o):
    if isinstance(o,dict):
        return {k:clean(v) for k,v in o.items()}
    if isinstance(o,np.ndarray):
        return clean(o.tolist())
    if isinstance(o,list):
        return [clean(v) for v in o]
    if isinstance(o,(float,np.floating)):
        return float(f"{o:.10g}")
    return o


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    fine,coarse=solve(.125),solve(.25)
    assert max(abs(fine["all_three_eigenline_minimum_sines"]-coarse["all_three_eigenline_minimum_sines"]))<3e-5
    result=clean(dict(fine=fine,step_halving_max_sine_shift=max(abs(fine["all_three_eigenline_minimum_sines"]-coarse["all_three_eigenline_minimum_sines"]))))
    path=HERE/"independent_results.json"
    if args.check:
        assert result==json.loads(path.read_bytes())
        print("1747 three dark-line possibilities and quantitative leakage bound verified")
    else:
        with path.open("x",encoding="utf-8") as handle:
            json.dump(result,handle,indent=2,allow_nan=False)
            handle.write("\n")

