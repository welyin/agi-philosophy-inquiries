"""H1 published DIS form-factor covariances and a pole-matched spectator loop."""
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
CFG = json.loads((HERE / "inputs.json").read_bytes())
W, Z = np.square(CFG["MW_MZ_GeV"])
K = CFG["alpha_Thomson_inverse"]/(4*math.pi)
C = W/Z
S = 1-C
V0 = W*S*K
nodes, weights = np.polynomial.legendre.leggauss(64)
u = (1-nodes**2)/4


def covariance(sample):
    cfg = CFG[sample]
    corr = np.eye(len(cfg["values"]))
    for i, row in enumerate(cfg["correlation_lower"]):
        corr[i, :i+1] = row
        corr[:i+1, i] = row
    assert np.linalg.eigvalsh(corr)[0] > 0
    return corr * np.outer(cfg["sigmas"], cfg["sigmas"])


def form(t, mass):
    kernel = np.log1p(-np.asarray(t)[..., None]*u/mass**2)
    return -3/(8*math.pi**2)*(kernel @ (weights*u/2))


def response(q2, mass, x, y):
    """Fractions of Z amplitude, effective neutral angle, W amplitude, photon amplitude."""
    hW, hZ = form(np.array([W,Z]), mass)
    h = form(-q2, mass)
    dk = ((C*C*x+S*S*y)*hZ - C*x*hW)/S
    dV = W*(dk+x*hW)
    k2 = dk+x*h
    kY = -dk+y*h
    n = -(dV+q2*(C*C*k2+S*S*kY))/(V0+q2*K*C*S)
    kappa = (k2-S*(k2+kY))/(K*S)
    cc = -(dV+q2*k2)/(V0+q2*K*S)
    photon = -(k2+kY)/K
    return np.array([n,kappa,cc,photon])


def exact_central_check():
    """Independent nonlinear calibration and propagator inversion as a derivative test."""
    mass = 100.
    q2 = 12345.
    x, y = 2/3, -2/3
    hW, hZ = form(np.array([W,Z]), mass)
    def normalized(e):
        left, right = K*.1, K*.4
        def poles(k):
            k2z, kYz = k+e*x*hZ, K-k+e*y*hZ
            return Z*k2z*kYz/(k2z+kYz) - W*(k+e*x*hW)
        assert poles(left) > 0 > poles(right)
        for _ in range(70):
            mid = (left+right)/2
            if poles(mid) > 0:
                left = mid
            else:
                right = mid
        k = (left+right)/2
        v = W*(k+e*x*hW)
        h = form(-q2,mass)
        a, b = k+e*x*h, K-k+e*y*h
        full = np.linalg.inv(q2*np.diag([a,b])+v*np.array([[1.,-1.],[-1.,1.]]))
        # In the (T3,Q) basis Y=Q-T3. Read Z, angle and massless photon residues.
        t = np.array([[1,0],[-1,1]])
        decomposed = t.T @ full @ t
        nc = decomposed[0,0]
        ang = -decomposed[0,1]/nc
        photon = decomposed[1,1]-nc*ang**2
        return np.array([nc*(V0+q2*K*C*S),ang/S,(V0+q2*K*S)/(v+q2*a),photon*q2*K])
    eps = .005
    np.testing.assert_allclose((normalized(eps)-normalized(-eps))/(2*eps),
                               response(q2,mass,x,y),atol=2e-12,rtol=4e-9)


def pull(v, cov):
    return math.sqrt(max(0., v @ np.linalg.solve(cov,v)))


def row_vector(sample, q2, mass, x, y):
    rows = response(q2,mass,x,y)
    return np.r_[rows[0],rows[1]] if sample == "NC" else rows[2]


def fit(sample, design):
    v = np.array(CFG[sample]["values"])-1
    Vinv = np.linalg.inv(covariance(sample))
    cov = np.linalg.inv(design.T @ Vinv @ design)
    mean = cov @ design.T @ Vinv @ v
    residual = v-design@mean
    return dict(mean=mean, sigma=np.sqrt(np.diag(cov)), correlation=cov/np.sqrt(
        np.outer(np.diag(cov),np.diag(cov))), chi2=float(residual@Vinv@residual))


def empirical(sample):
    edge = np.array(CFG[sample]["Q2_edges_GeV2"])
    logq = np.log(np.sqrt(np.sqrt(edge[1:]*edge[:-1]))/60.)
    X = np.column_stack([np.ones(4),logq])
    design = np.kron(np.eye(2),X) if sample == "NC" else X
    return fit(sample,design)


def full_activity_scale(nc, cc):
    # Full-activity high-Q normalization only. The actual H1 interval crosses
    # the weak scale and the finite kernel is used for its physical examples.
    s, alpha = .2312, 1/128
    c = 1-s
    A = alpha/(2*math.pi)*np.array([[c/s,s/c],[-c/s,1]])
    chosen = [1,3]
    sig = nc["sigma"][chosen]
    cov = nc["correlation"][np.ix_(chosen,chosen)]*np.outer(sig,sig)
    inv = np.linalg.inv(A)
    budget = inv @ cov @ inv.T
    return dict(ideal_NC_sigma_b2_bY=np.sqrt(np.diag(budget)),
                ideal_NC_corr=budget[0,1]/math.sqrt(np.prod(np.diag(budget))),
                ideal_CC_sigma_b2=cc["sigma"][1]*2*math.pi*s/alpha,
                angle_only_null_direction=[1.,c/s], examples_signal_per_logQ=A@np.array([2/3,-2/3]))


def spectator(sample, mass):
    x, y = np.subtract(*CFG["UV_beta_shifts_b2_bY"].values())
    edge = np.array(CFG[sample]["Q2_edges_GeV2"])
    pivots = np.sqrt(edge[:-1]*edge[1:])
    cov = covariance(sample)
    centers = np.array(CFG[sample]["values"])
    at = row_vector(sample,pivots,mass,x,y)
    ends = np.array([row_vector(sample,q2,mass,x,y) for q2 in (edge[:-1],edge[1:])])
    # Monotonic response proved in proof.md. Independent positive-bin means may
    # occupy the enclosing box, whose maximal quadratic distance is at a vertex.
    corners = itertools.product(*ends.T)
    distance_max = max(pull(np.array(v),cov) for v in corners)
    fits = {}
    for key, amounts in CFG["UV_beta_shifts_b2_bY"].items():
        shift = row_vector(sample,pivots,mass,*amounts)
        fits[key] = dict(response_centers=shift, conditional_chi2_pivot=pull(centers-1-shift,cov)**2)
    physical_design = np.column_stack([row_vector(sample,pivots,mass,*pair) for pair in [(1.,0.),(0.,1.)]])
    # Local errors for a two-coefficient oblique plane, with independently
    # removed common photon response in NC. Not a large-beta global fit.
    fisher = physical_design.T @ np.linalg.solve(cov,physical_design)
    local_cov = np.linalg.inv(fisher)
    return dict(shift_difference=at, distance=pull(at,cov), distance_positive_bin_upper=distance_max,
                max_abs_difference_in_bins=np.max(np.abs(ends),axis=0),
                SM_chi2=pull(centers-1,cov)**2, named_hypotheses=fits,
                local_sigma_xy=np.sqrt(np.diag(local_cov)),
                local_corr_xy=local_cov[0,1]/np.sqrt(np.prod(np.diag(local_cov))))


def run():
    exact_central_check()
    assert all(0.5*math.sqrt(Z)<m for m in CFG["comparison_masses_GeV"])
    assert abs(form(0,100.))==0
    # UV matching of the common Euclidean threshold to the beta convention.
    hi = abs(form(-1e14,100.)-form(-1e12,100.))
    assert abs(hi-math.log(10)/(8*math.pi**2))<1e-8
    named = CFG["UV_beta_shifts_b2_bY"]
    assert abs(sum(named["doublet"])-sum(named["charged_plus_neutral_singlets"]))<1e-14
    for m in CFG["comparison_masses_GeV"]:
        np.testing.assert_allclose(response(np.array([561,6000,77000]),m,2/3,-2/3)[3],0,atol=1e-16)
    slopes = {key:empirical(key) for key in ("NC","CC")}
    return dict(empirical_per_logQ=slopes,
                hypothetical_full_activity=full_activity_scale(slopes["NC"],slopes["CC"]),
                finite_thresholds={str(m):{key:spectator(key,m) for key in ("NC","CC")}
                                   for m in CFG["comparison_masses_GeV"]})


def tidy(o):
    if isinstance(o, dict):
        return {k:tidy(v) for k,v in o.items()}
    if isinstance(o,np.ndarray):
        return tidy(o.tolist())
    if isinstance(o,list):
        return [tidy(v) for v in o]
    if isinstance(o,float):
        return float(f"{o:.10g}")
    return o


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--check",action="store_true")
    result = tidy(run())
    path = HERE/"results.json"
    payload = json.dumps(result,indent=2)+"\n"
    if ap.parse_args().check:
        assert json.loads(path.read_bytes())==result
    else:
        path.write_text(payload,encoding="utf8")
    print("1783 H1 correlated NC and CC pole matching",hashlib.sha256(payload.encode()).hexdigest()[:12])
