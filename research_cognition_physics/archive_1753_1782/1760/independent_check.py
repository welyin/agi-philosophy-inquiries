"""Seven-state spacelike propagator, its controlled expansion and an electron-rate ambiguity."""
import argparse
import json
import math
import runpy
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
EV = 1.e-9


def read(path):
    return json.loads(path.read_bytes())


def material():
    prior = read(HERE.parent/"1757"/"independent_results.json")["combined_central_solution"]
    audit = runpy.run_path(str(HERE.parent/"1758"/"independent_check.py"))
    _, masses, L, kernels, _ = audit["physical"](prior["Vus"], prior["Vub"])
    unscaled = kernels[0]+prior["rho"]*kernels[1]
    kappa = math.sqrt(audit["SOLAR"])/np.linalg.svd(L@unscaled@L.T,compute_uv=False)[1]
    fine = (L@np.diag(masses), np.diag(masses), kappa*unscaled*EV)
    # Use exactly the main study's coarse-step endpoint in its fractional bound.
    archive = runpy.run_path(str(HERE.parent/"1758"/"check.py"))
    _, p = archive["endpoint"](prior["Vub"],.25)
    K = p["C0"]+prior["rho"]*np.outer(p["b0"],p["b0"])
    scale = math.sqrt(audit["SOLAR"])/np.linalg.svd(K,compute_uv=False)[1]
    inv = np.linalg.pinv(p["L"])
    D = np.diag(p["M"])
    frozen = (p["L"]@D,D,inv@K@inv.T*(EV*scale))
    return frozen,fine


def mass_matrix(md,D,mu):
    H = np.zeros((5,5), complex)
    H[3:,3:] = mu
    K = np.vstack([md,D])
    return np.block([[H,K],[K.T,np.zeros((2,2))]])


def full(q,md,D,mu):
    M = mass_matrix(md,D,mu)
    G = q*q*np.linalg.solve(q*q*np.eye(7)+M@M.conj().T,M)
    return G[:3,:3]


def resummed(q,md,D,mu):
    K = np.vstack([md,D])
    P = np.eye(5)-K@np.linalg.solve(q*q*np.eye(2)+K.conj().T@K,K.conj().T)
    a = P[:3,3:]
    p = P[3:,3:]
    c = np.eye(2)+p.T@mu.conj().T@p@mu/(q*q)
    return a@mu@np.linalg.solve(c,a.T), a@mu@a.T, a


def remainder(q,port,norm):
    return np.linalg.norm(port[0])**2*norm**3/(q*q-norm**2)


def as_complex(a):
    return {"real":a.real.tolist(),"imag":a.imag.tolist()}


def ambiguity(md,D,qs):
    # Scale out M1 in the quadratic in z=Q^2/M1^2.
    m = md/D[0,0]
    d = D/D[0,0]
    h = d@d+m.conj().T@m
    coeff = np.array([(m@(np.trace(h)*np.eye(2)-h)@d)[0],(m@d)[0]])
    diag = np.diag([5.e-6,3.2e-6*np.exp(.85j)])
    poly = coeff@diag@coeff.T
    c = np.array([poly[1,1],2*poly[0,1],poly[0,0]])
    roots = np.roots(c)
    roots = roots[np.argsort(roots.real)]
    flipped = roots.copy()
    flipped[0] = flipped[0].conjugate()
    t = c[0]*np.poly(flipped)
    target = np.array([[t[2],t[1]/2],[t[1]/2,t[0]]])
    inverse = np.linalg.inv(coeff)
    rival = inverse@target@inverse.T
    adj = (rival+rival.T)/2
    assert np.linalg.norm(adj-diag)>1e-6 and abs(adj[0,1])>1e-6
    records = []
    for q in qs:
        f0, g0, a = resummed(q,md,D,diag)
        f1, g1, _ = resummed(q,md,D,adj)
        kernel = abs(g0[0,0])
        np.testing.assert_allclose(abs(g1[0,0]),kernel,rtol=2e-13,atol=0)
        for model,test in [(diag,f0),(adj,f1)]:
            np.testing.assert_allclose(full(q,md,D,model),test,atol=1e-22,rtol=1e-12)
        r = (remainder(q,a,np.linalg.norm(diag,2))+remainder(q,a,np.linalg.norm(adj,2)))/kernel
        one = remainder(q,a,np.linalg.norm(diag,2))/kernel
        limit = r/(1-one)
        records.append(dict(q_GeV=q,ee_rate_analytic_relative_upper=2*limit+limit**2,
                            em_rate_ratio=abs(f1[0,1]/f0[0,1])**2,
                            muon_rate_ratio=abs(f1[1,1]/f0[1,1])**2))
    return dict(diagonal_mu_keV=as_complex(diag*1e6),reflected_mu_keV=as_complex(adj*1e6),
                roots_dimensionless=as_complex(roots),checks=records)


def evaluate():
    frozen, fine = material()
    md,D,mu = frozen
    inp = read(HERE/"inputs.json")
    qs = inp["spacelike_momenta_GeV"]
    bound = read(HERE/"results.json")["certified_minimum_max_fractional_rate_error"]
    # The cutoff specified here concerns candidate Majorana masses, not momenta.
    cap = .010  # GeV
    eta = np.linalg.norm(mu,2)
    rows = []
    amounts = []
    for q in qs:
        original = full(q,md,D,mu)
        fine_ee = full(q,*fine)[0,0]
        exact, first, a = resummed(q,md,D,mu)
        np.testing.assert_allclose(original,exact,rtol=3e-12,atol=1.e-21)
        high = (md@np.diag(np.diag(D)/(q*q+np.diag(D)**2)))[0]
        scale = abs(first[0,0])
        btarget = remainder(q,a,eta)/scale
        bcomp = remainder(q,a,cap)/scale
        amounts.append(bcomp+math.sqrt(1.01)*btarget)
        rows.append(dict(q_GeV=q,full_ee_meV=abs(original[0,0])*1e12,
                         exact_ee_meV=abs(exact[0,0])*1e12,
                         independent_half_step_ee_rate_shift=abs(fine_ee/original[0,0])**2-1,
                         mu_linearization_relative_amp_upper=btarget,
                         diagonal_10MeV_absolute_error_over_target=bcomp,
                         no_Dirac_dressing_relative_rate=abs(high@mu@high.T/first[0,0])**2-1))
    loss = max(amounts)
    all_order_floor = bound-2*math.sqrt(1.01)*loss-loss**2
    assert all_order_floor>1e-4
    # Use an intentionally large Majorana block to test the finite-order formula
    # where the cubic correction can be numerically resolved in double precision.
    large = mu/eta*80.
    check, first, a = resummed(qs[0],md,D,large)
    np.testing.assert_allclose(full(qs[0],md,D,large),check,rtol=1e-11,atol=1.e-20)
    assert 1e-9<abs(check[0,0]/first[0,0]-1)<.01
    rephase = np.diag(np.exp(1j*np.array([1.22,-.73])))
    np.testing.assert_allclose(full(qs[0],md@rephase.conj(),D,rephase@large@rephase),
                               check,rtol=1e-11,atol=1.e-20)
    return dict(physical_majorana_norm_keV=eta*1e6,checkpoints=rows,
                conditional_all_order_diagonal_mu_budget_MeV=cap*1e3,
                conditional_all_order_fractional_rate_floor=all_order_floor,
                large_80GeV_control_relative_cubic=abs(check[0,0]/first[0,0]-1),
                electron_rate_ambiguity=ambiguity(md,D,qs))


if __name__=="__main__":
    cli = argparse.ArgumentParser()
    cli.add_argument("--check",action="store_true")
    args = cli.parse_args()
    result = json.loads(json.dumps(evaluate()),parse_float=lambda v:float(f"{float(v):.9g}"))
    dest = HERE/"independent_results.json"
    if args.check:
        assert result==read(dest)
        print("1760 seven-state resolvent, bounded all-order diagonal test and root-reflected rivals verified")
    else:
        dest.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf-8")
