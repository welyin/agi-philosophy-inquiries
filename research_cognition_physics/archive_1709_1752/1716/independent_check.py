"""A phase-free singular-value bound and a finite exclusion cone for the named ray."""
import argparse
import json
import math
from pathlib import Path
import runpy

import numpy as np


HERE = Path(__file__).resolve().parent
PAIRS = [(0, 1), (0, 2), (1, 2)]


def read(path):
    return json.loads((HERE / path).read_bytes())


def certificate(X, v):
    # Each column is projected onto the Hermitian complement of the holomorphic null line.
    W = np.linalg.svd(v.conj()[None, :], full_matrices=True)[2][1:].conj().T
    F = W.conj().T @ X
    pieces = [np.outer(x, x.conj()) for x in X.T]
    candidates = []
    for ii in range(3):
        L = 2*pieces[ii]-sum(pieces)
        candidates.append((np.linalg.eigvalsh(W.conj().T @ L @ W)[-1], L))
    pick = int(np.argmax([t[0] for t in candidates]))
    ell, L = candidates[pick]
    area = np.array([abs(np.linalg.det(F[:, [i,j]]))**2 for i,j in PAIRS])
    exterior = np.column_stack([np.cross(X[:,i], X[:,j]).conj() for i,j in PAIRS])
    Q = exterior @ exterior.conj().T
    den = float(np.vdot(v, Q@v).real)
    np.testing.assert_allclose(den, area.sum(), rtol=2e-13)
    return W, F, {"strongest_row_index":pick, "projected_singular_value_floor":ell,
                  "absolute_minor_terms":area.tolist(), "upper_determinant":den,
                  "lower_determinant":2*max(area)-sum(area), "mass_ratio_floor":ell**2/den,
                  "support_spread":np.ptp(np.linalg.eigvalsh(L)),
                  "exterior_spread":np.ptp(np.linalg.eigvalsh(Q))}


def tolerance(bound, target, condition):
    ell = bound["projected_singular_value_floor"]
    den = bound["upper_determinant"]
    dl, dd = bound["support_spread"], bound["exterior_spread"]
    ceiling = target*condition**2
    residual = lambda s:max(0., ell-dl*s)**2-ceiling*(den+dd*s)
    assert residual(0) > 0 and residual(1) < 0
    lo, hi = 0., 1.
    for _ in range(64):
        mid = (lo+hi)/2
        if residual(mid)>0:
            lo=mid
        else:
            hi=mid
    floor = bound["mass_ratio_floor"]/condition**2
    return {"condition_number_budget":condition, "low_mass_ratio_floor":floor,
            "necessary_ray_angle_deg":math.degrees(math.asin(lo)),
            "necessary_fractional_low_mass_operator_correction":(1-target/floor)/(1+target)}


def main():
    archive = runpy.run_path(str(HERE / "../1715/check.py"))
    flow, rep = archive["charged_background"](.10)
    G, _ = archive["congruence"](flow,rep,.10)
    D = archive["unpack"](rep["down_texture_10PeV"])
    N = archive["unpack"](rep["neutral_columns_10PeV"])
    up = np.array(rep["up_texture_diag_10PeV"])
    U, singular, _ = np.linalg.svd(N)
    v = U[:,-1]
    assert np.linalg.norm(N.conj().T@v)<1e-14 and singular[1]>.03
    X = np.linalg.solve(D,np.diag(np.sqrt(up)))
    W,F,cert = certificate(X,v)
    assert cert["lower_determinant"]>0
    # Direct inverse, flavour congruence and determinant lemma, all at a complex point.
    phases = np.exp(1j*np.array([.13,-.72,1.88]))
    C = D.T @ np.diag(1/(up*phases)) @ D
    w = C@v
    scalar = v@w
    Cv = C-np.outer(w,w)/scalar
    K = F@np.diag(phases)@F.T
    test = W.T @ Cv @ W
    np.testing.assert_allclose(test@K,np.eye(2),atol=1e-12)
    s = np.linalg.svd(Cv,compute_uv=False)
    assert s[0]/s[1] > cert["mass_ratio_floor"]
    # A self mass of the added singlet directly modifies the cancelling coefficient.
    mu = .18j*scalar
    Cmu = C+np.outer(w,w)/(mu-scalar)
    np.testing.assert_allclose(np.linalg.det(Cmu)/np.linalg.det(C),mu/(mu-scalar),rtol=1e-13)
    exp = read("../1668/inputs.json")
    max_ratio = math.sqrt((exp["juno_NO_dm31"][0]+exp["juno_NO_dm31"][2])/
                          (exp["juno_NO_dm21"][0]-exp["juno_NO_dm21"][1]))
    eta = read("../1715/independent_results.json")["all_interval_conditional_log_spread"]
    return {"high_scale_GeV":flow["END"], "named_port_norm_residual":np.linalg.norm(N.conj().T@v),
            "named_port_electron_probability_high":abs(np.vdot(archive["unpack"](rep["charged_left_10PeV"])[:,0],v))**2,
            "high_scale_phase_independent_bound":cert, "JUNO_normal_one_width_maximum_ratio":max_ratio,
            "minimum_congruence_condition_to_reach_data":math.sqrt(cert["mass_ratio_floor"]/max_ratio),
            "passive_Weinberg_competitor":tolerance(cert,max_ratio,math.exp(eta)),
            "registered_active_flavour_transport":tolerance(cert,max_ratio,np.linalg.cond(G)),
            "arbitrary_transport_within_one_percent":tolerance(cert,max_ratio,1.01),
            "complex_point_inverse_compression_error":np.linalg.norm(test@K-np.eye(2),2)}


def tidy(x):
    if isinstance(x,dict): return {k:tidy(v) for k,v in x.items()}
    if isinstance(x,list): return [tidy(v) for v in x]
    if isinstance(x,(float,np.floating)): return float(f"{x:.10g}")
    return x


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    val=tidy(main())
    target=HERE/"independent_results.json"
    if parser.parse_args().check:
        assert json.loads(target.read_bytes())==val
        print("1716 verified the Schur phase bound, ray perturbation and contact restriction")
    else:
        target.write_text(json.dumps(val,indent=2)+"\n",encoding="utf8")
