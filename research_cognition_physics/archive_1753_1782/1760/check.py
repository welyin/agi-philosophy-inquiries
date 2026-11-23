"""Four central-angle rates test all diagonal complex LNV sources in the same two-port completion."""
import argparse
import json
import math
from pathlib import Path
import runpy

import numpy as np


HERE = Path(__file__).resolve().parent
CFG = json.loads((HERE/"inputs.json").read_bytes())


def load_model():
    source = json.loads((HERE/CFG["endpoint_from"]).read_bytes())["combined_central_solution"]
    stage = runpy.run_path(str(HERE/CFG["masses_and_ports"]))
    row, data = stage["endpoint"](source["Vub"], .25)
    assert abs(row["low_moduli"]["Vus"] - source["Vus"]) < 1e-7
    rho = source["rho"]
    # The low-energy amplitude is fixed by the independent clocks at this point.
    low = data["C0"] + rho*np.outer(data["b0"], data["b0"])
    scale = math.sqrt(CFG["solar_dm2_eV2"]) / np.linalg.svd(low, compute_uv=False)[1]
    theta = data["L"]
    inverse = np.linalg.pinv(theta)
    mu = scale*inverse@low@inverse.T  # eV
    np.testing.assert_allclose(theta@mu@theta.T, scale*low, rtol=1e-11)
    return data["M"], theta, mu


def central_rate(q, amplitude_eV):
    mw = CFG["mW_GeV"]
    s = 2*(q*q + mw*mw)
    t = -q*q
    r = mw**2/s
    beta = math.sqrt(1-4*r)
    bt = (1-4*r)*t**2 - 4*r*(1-2*r)*s*t + 4*(1-r)*r**2*s**2
    btu = (1-4*r)*t*t + 4*r**3*s*s
    dsigma = CFG["GF_GeVminus2"]**2*beta*(2*bt+2*btu)/(32*math.pi)
    return math.sqrt(s), dsigma*abs(amplitude_eV*1e-9/t)**2*CFG["GeVminus2_to_fb"]


def solve():
    masses, theta, mu = load_model()
    D = np.diag(masses)
    md = theta@D
    eta = np.linalg.norm(theta,2)**2
    h = D@D + md.conj().T@md

    def legs(q):
        f = masses**2/(q*q + masses**2)
        # Mixing retained to all orders, light masses negligible compared with q.
        exact = (md@np.linalg.solve(q*q*np.eye(2)+h,D))[0]
        return exact, theta[0]*f, f

    momenta = np.array(CFG["spacelike_momenta_GeV"])
    calculated = [legs(q) for q in momenta]
    a = np.array([v[0] for v in calculated])
    f = np.array([v[2] for v in calculated])
    amp = np.einsum("bi,ij,bj->b",a,mu,a)
    contact = theta[0]@mu@theta[0]
    assert abs(contact)>1e-5
    # Normalize by the first mass propagator and the reference zeroth-order m_ee.
    factor = f[:,0]**2*abs(contact)
    y = abs(amp/factor)**2
    x = (f[:,1]/f[:,0])**2
    # All leading diagonal intensities are polynomials of degree at most two in x.
    w = np.array([1/np.prod(t-np.delete(x,i)) for i,t in enumerate(x)])
    w /= sum(abs(w))
    for k in range(3):
        assert abs(w@x**k)<1e-12
    # Restore exact Dirac mixing. For a diagonal competitor define its complex
    # coefficients in units theta_ei^2 / |m_ee|. Each row below maps its two
    # coefficients to the normalized complex amplitude.
    H = (a/theta[0])**2/f[:,0,None]**2
    Q = H.conj().T@(w[:,None]*H)
    transfer = np.linalg.inv(H[[0,-1],:])
    intensity = abs(w@y)
    # If every normalized rate is within relative error eps of the target,
    # the outer two rates bound ||coefs||^2 <= (1+eps)*coef_cap.
    coef_cap = np.linalg.norm(transfer,2)**2*(y[0]+y[-1])
    delta = np.linalg.norm(Q,2)*coef_cap
    signal_scale = abs(w)@y
    relative_floor = max(0.,(intensity-delta)/(signal_scale+delta))
    assert relative_floor>1e-5
    energies,rates = zip(*(central_rate(q,m) for q,m in zip(momenta,amp)))
    nuclear_q = CFG["probe_for_nuclear_kernel_only_GeV"]
    n0 = legs(0)[0]
    nq = legs(nuclear_q)[0]
    quasi_zero = n0@mu@n0
    # The q=0 expression here means the mass-window limit m_light << q << M.
    return dict(M_GeV=masses.tolist(),theta_norm_squared=eta,
                weak_m_ee_meV=1000*abs(contact),finite_mixing_low_m_ee_meV=1000*abs(quasi_zero),
                fractional_kernel_change_at_200MeV=abs((nq@mu@nq-quasi_zero)/quasi_zero),
                energy_GeV=list(energies),q_GeV=momenta.tolist(),central_dsigma_dcos_fb=list(rates),
                four_reduced_rates=y.tolist(),diagonal_test_weights=w.tolist(),four_weight_residual=float(w@y),
                finite_mixing_allowance_at_exact_outer_rates=delta,weighted_target_scale=signal_scale,
                certified_minimum_max_fractional_rate_error=relative_floor,
                comparison_contact_dsigma_dcos_fb=[central_rate(q,contact)[1] for q in momenta],
                max_relative_mu_heavy_mass=np.linalg.norm(mu,2)*1e-9/masses.min())


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=solve()
    out=json.loads(json.dumps(result),parse_float=lambda v:float(f"{float(v):.10g}"))
    path=HERE/"results.json"
    if args.check:
        assert out==json.loads(path.read_bytes())
        print("1760: two-pole momentum test, finite Dirac mixing and central cross sections verified")
    else:
        with path.open("x",encoding="utf-8") as handle:
            json.dump(out,handle,indent=2,allow_nan=False)
            handle.write("\n")
