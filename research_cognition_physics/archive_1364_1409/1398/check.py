"""Direct semileptonic amplitudes bound weak accessibility of the Yukawa remainder."""
import argparse
import json
from math import sqrt
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE/path).read_text(encoding="utf-8"))


def geometry(up, down, width):
    def extrema(rows):
        return np.array([v-width*e for v,e in rows]), np.array([v+width*e for v,e in rows])
    ul, uh = extrema(up)
    dl, dh = extrema(down)
    q_u, q_c = (uh[:2]/ul[2])**2
    q_d, q_s = (dh[:2]/dl[2])**2
    return float(1+q_u+q_c+q_s), [float(q_u+q_d),float(q_c+q_s)]


def bound(ud, cs, up, down, n):
    a = ud[0]-n*sum(ud[1:])
    b = cs[0]-n*sum(cs[1:])
    visible, diagonal_caps = geometry(up, down, n)
    observed = a*a+b*b
    assert 0 < visible < observed < 2 and visible/2 < ((a+b)/2)**2
    active_block = sqrt((observed-visible)/(2-visible))
    v, mean = visible/2, (a+b)/2
    receiving = sqrt(1-v)*mean-sqrt(v)*sqrt(1-mean**2)
    assert 0 < receiving < active_block < 1
    return {"amplitude_lower_ud_cs":[a,b],"two_measured_squares_lower":observed,
            "visible_two_row_strength_upper":visible, "direct_excess_strength_lower":observed-visible,
            "minimum_common_apparent_amplitude_scale_for_inert":sqrt(observed/visible),
            "visible_ordinary_W_max_ud_cs":diagonal_caps,
            "disjoint_weak_blocks_eta_lower":active_block,
            "arbitrary_coherent_receiving_norm_eta_lower":receiving}


def counterexample(ud, cs, up, down, witness):
    # No input CKM fit. Two measured moduli determine these two trigonometric parameters.
    s13 = witness["s13"]
    phase = witness["delta"]
    assert abs(np.cos(phase)) < 1e-14
    c13 = sqrt(1-s13*s13)
    c12 = ud/c13
    s12 = sqrt(1-c12*c12)
    c23 = sqrt((cs**2-s12**2*s13**2)/(c12**2-s12**2*s13**2))
    s23 = sqrt(1-c23*c23)
    r12 = np.array([[c12,s12,0],[-s12,c12,0],[0,0,1]], dtype=complex)
    r23 = np.array([[1,0,0],[0,c23,s23],[0,-s23,c23]])
    z = np.exp(1j*phase)
    u13 = np.array([[c13,0,s13/z],[0,1,0],[-s13*z,0,c13]])
    v = r23 @ u13 @ r12
    yu, yd = np.array(up)[:,0], np.array(down)[:,0]
    # Pick an independently specified illustrative strength large enough for both ports.
    k = 2.
    visible_part = np.vstack([np.diag(yu)/k,np.diag(yd) @ v.conj().T/k])
    eig, basis = np.linalg.eigh(np.eye(3)-visible_part.conj().T @ visible_part)
    assert min(eig) > .5
    t = np.vstack([visible_part,np.diag(np.sqrt(eig)) @ basis.conj().T])
    right = t @ v
    weak = t.conj().T @ right
    np.testing.assert_allclose(weak,v,atol=1e-14)
    np.testing.assert_allclose([abs(weak[0,0]),abs(weak[1,1])],[ud,cs],atol=1e-13)
    weak0 = visible_part.conj().T @ visible_part @ v
    expected = (yu[:,None]**2+yd[None,:]**2)*v/k**2
    np.testing.assert_allclose(weak0,expected,atol=1e-14)
    j = float((weak[0,0]*weak[1,1]*weak[0,1].conjugate()*weak[1,0].conjugate()).imag)
    assert abs(j) > 1e-6
    return {"k":k, "s12_s13_s23_delta":[s12,s13,s23,phase],
            "abs_ordinary_all_active_ud_cs":[float(abs(weak[0,0])),float(abs(weak[1,1]))],
            "imag_quartet":j, "D_eigenvalues":eig.tolist(),
            "same_spectrum_only_visible_abs_ud_cs":[float(abs(weak0[0,0])),float(abs(weak0[1,1]))]}


def run():
    cfg = read("inputs.json")
    measured = read(cfg["charged"])
    masses = read(cfg["light_mass_inputs"])["yukawa"]
    ud,cs = measured["ud_beta"], measured["charm_s"]
    u,d = masses["up"],masses["down"]
    return {"widths":{str(n):bound(ud,cs,u,d,n) for n in cfg["width_multipliers"]},
            "all_active_counterexample":counterexample(ud[0],cs[0],u,d,cfg["all_active_example"])}


def verify(a,b):
    if isinstance(a, dict):
        assert a.keys()==b.keys()
        for k in a:
            verify(a[k],b[k])
    elif isinstance(a, list):
        assert len(a)==len(b)
        for x,y in zip(a,b):
            verify(x,y)
    else:
        np.testing.assert_allclose(a,b,rtol=1e-9,atol=3e-15)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    if args.write:
        (HERE/"results.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    elif args.check:
        verify(result,read("results.json"))
        print("1398 weak-current capacity bounds and fully active counterexample verified")
    else:
        print(json.dumps(result,indent=2))
