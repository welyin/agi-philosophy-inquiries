"""Exact rank-one identities, finite perturbations and MZ data, using NumPy only."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

BASE = Path(__file__).resolve().parent


def read(name):
    return json.loads((BASE/name).read_text(encoding="utf-8-sig"))


def mixing_measurement(unitary):
    s13 = abs(unitary[0, 2])
    c13 = math.sqrt(1-s13**2)
    s12, s23 = abs(unitary[0, 1])/c13, abs(unitary[1, 2])/c13
    c12, c23 = math.sqrt(1-s12**2), math.sqrt(1-s23**2)
    j = (unitary[0,0]*unitary[1,1]*unitary[0,1].conj()*unitary[1,0].conj()).imag
    sin_delta = j/(s12*c12*s23*c23*s13*c13**2)
    cos_delta = (s12**2*s23**2+c12**2*c23**2*s13**2-abs(unitary[2,0])**2)/(2*s12*s23*c12*c23*s13)
    return [s12, s13, s23, math.atan2(sin_delta,cos_delta)], float(j)


def lower_invariant(weak, width):
    d, s, b = weak["yukawa"]["down"]
    params = weak["quark_sines_and_phase"]
    lo = [x-width*e for x,e in params]
    hi = [x+width*e for x,e in params]
    j = lo[0]*lo[1]*lo[2]*math.sin(lo[3])*(1-hi[1]**2)*math.sqrt((1-hi[0]**2)*(1-hi[2]**2))
    rd_hi = (d[0]+width*d[1])/(b[0]-width*b[1])
    rs_lo = (s[0]-width*s[1])/(b[0]+width*b[1])
    rs_hi = (s[0]+width*s[1])/(b[0]-width*b[1])
    return j*(rs_lo**2-rd_hi**2)*(1-rd_hi**2)*(1-rs_hi**2)


def verify_texture(d, u, v, c):
    y = np.diag(d)+c*np.outer(u,v)
    w = d*v
    h = y@y.conj().T
    inv = float((h[0,1]*h[1,2]*h[2,0]).imag)
    area = [(u[i]*w[j]-u[j]*w[i]) for i,j in ((0,1),(1,2),(0,2))]
    np.testing.assert_allclose(inv, abs(c)**2*c.imag*np.prod(area), rtol=1e-9, atol=2e-14)
    det_ratio = np.linalg.det(y)/np.prod(d)
    np.testing.assert_allclose(det_ratio, 1+c*np.dot(u/d,v), rtol=1e-11, atol=2e-13)
    q0 = np.abs(c.imag*u*v/d)
    assert abs(inv) <= abs(c)**2*max(d)**4*max(q0)+5e-14
    return y, q0


def certificate():
    data = read("inputs.json")
    weak = read(data["weak_matching_inputs"])
    edm = read(data["neutron_EDM_input"])["edm"]
    theta_lim = {}
    for name, values in [("historical",edm["abs_slope_e_fm"]),("alternative",data["second_lattice"]["slope_magnitude_e_fm"])]:
        center, e1, e2 = values
        theta_lim[name] = [edm["limit_90_ecm"]/(edm["fm_in_cm"]*(center-n*(e1+e2))) for n in range(3)]
    inv = [lower_invariant(weak,n) for n in range(3)]
    budgets = data["channel_budgets_over_physical_yb"]
    guard = {}
    for a in budgets["rank_one_max"]:
        q = [t/(a**2*budgets["direct_max"]**4) for t in inv]
        guard[str(a)] = dict(q_minimum=q, two_width_direct_change_limit={
            basis: math.tan(bounds[2])/q[2] for basis,bounds in theta_lim.items()})
    witness = data["witness_in_units_of_yb"]
    d = np.array(witness["direct"])
    u = np.array(witness["u"])
    u /= np.linalg.norm(u)
    k = u/d
    v = np.array(witness["v_seed"])
    v -= k*(np.dot(k,v)/np.dot(k,k))
    v /= np.linalg.norm(v)
    assert witness["phase_in_units_of_pi"] == -.5
    c = -1j*witness["A"]
    y, q = verify_texture(d,u,v,c)
    columns, eig, _ = np.linalg.svd(y)
    mixing, j = mixing_measurement(columns[:, ::-1])
    masses = eig[::-1]*weak["yukawa"]["down"][-1][0]
    target_m = [x[0] for x in weak["yukawa"]["down"]]
    np.testing.assert_allclose(masses, target_m, rtol=2e-10)
    np.testing.assert_allclose(mixing, [x[0] for x in weak["quark_sines_and_phase"]], rtol=2e-10)
    assert abs(np.angle(np.linalg.det(y))) < 2e-14
    eps = data["inverse_direct_independent_fraction"]
    axis = int(np.argmax(q))
    shifted = []
    for sign in (-1,1):
        shifted_d = d.copy()
        shifted_d[axis] /= 1+sign*eps
        z = np.diag(shifted_d)+c*np.outer(u,v)
        shifted.append(float(np.angle(np.linalg.det(z))))
    np.testing.assert_allclose(np.abs(shifted), math.atan(eps*max(q)),atol=1e-13)
    # Independent real/imag directions in 32 general rank-one matrices.
    gen = np.random.default_rng(1414)
    for _ in range(32):
        p,r = gen.normal(size=(2,3))
        verify_texture(gen.uniform(.05,2.,3), p/np.linalg.norm(p),r/np.linalg.norm(r),
                       gen.uniform(.1,2)*np.exp(1j*gen.uniform(-np.pi,np.pi)))
    return dict(weak_invariant_over_yb6_lower=inv, theta_target=theta_lim, guards=guard,
                witness=dict(masses_MZ=masses.tolist(), CKM_s12_s13_s23_phase=mixing,
                             J=j, inverse_axis_sensitivities=q.tolist(),
                             five_percent_signed_theta=shifted,
                             historical_two_width_inverse_change=math.tan(theta_lim["historical"][2])/max(q)),
                identity_checks=32)


def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:
            compare(a[k],b[k])
    else:
        np.testing.assert_allclose(a,b,rtol=5e-10,atol=0)


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    cli.add_argument("--write",action="store_true")
    result=certificate()
    if cli.parse_args().write:
        with (BASE/"results.json").open("x",encoding="utf-8",newline="\n") as f:
            json.dump(result,f,indent=2,allow_nan=False)
            f.write("\n")
    else:
        compare(result,read("results.json"))
    print(json.dumps(result,allow_nan=False))
