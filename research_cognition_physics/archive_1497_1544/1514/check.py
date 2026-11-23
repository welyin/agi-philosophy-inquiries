"""One versus two antifive flavour directions in an adjoint mass splitting."""

import argparse
import itertools
import json
import math
from pathlib import Path
import numpy as np


HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def one_direction(d, e):
    ratios = e/d
    # Arrays use ascending singular values.
    grow = [float(max(e[:-1]/d[1:])), float(min(ratios))]
    shrink = [float(max(ratios)), float(min(e[1:]/d[:-1]))]
    assert grow[0] < grow[1] < 1 < shrink[0] < shrink[1]
    volume = float(np.prod(ratios))
    return {"mass_ratios": ratios.tolist(), "k_for_growing": grow, "k_for_shrinking": shrink,
            "minimum_strain_growing": volume / grow[1]**3,
            "minimum_strain_shrinking": shrink[0]**3 / volume}


def downdate_weights(parent, target):
    """Squared overlaps inside the three nonzero singular directions of Y."""
    c, t = parent**2, target**2
    return np.array([np.prod(c[j]-t)/(c[j]*np.prod(c[j]-np.delete(c,j)))
                     for j in range(len(c))])


def unit_heavy(parent, target):
    weights = downdate_weights(parent, target)
    assert np.min(weights) > -1e-13 and sum(weights)<1+1e-12
    return np.sqrt(np.append(np.maximum(weights,0), max(0.0,1-sum(weights))))


def light_spectrum(Y, n):
    # Vh[1:] spans the orthogonal complement of the heavy mass column.
    _, _, vh = np.linalg.svd(n.conj()[None,:], full_matrices=True)
    physical = np.linalg.svd(Y @ vh[1:].conj().T, compute_uv=False)[::-1]
    return physical, vh[1:].conj().T


def aligned_witness(d,e,k):
    scaled = e/k
    rank_delta = np.array([np.prod(d[j]**2-scaled**2) /
                            np.prod(d[j]**2-np.delete(d**2,j))
                            for j in range(3)])
    rank_delta = -rank_delta  # H_target = diag(d^2) + signed vv*
    sign = np.sign(rank_delta[np.argmax(np.abs(rank_delta))])
    assert min(sign*rank_delta) > -1e-13
    a = np.sqrt(np.maximum(sign*rank_delta,0))/d
    v = a/np.linalg.norm(a)
    rho = float(np.prod(scaled/d))
    assert math.isclose(rho*rho, 1+sign*np.vdot(a,a).real, rel_tol=3e-10)
    ad = min(1.0,1/rho)/2
    ae = ad*rho
    bare = np.diag(d) @ (np.eye(3)+(1/ad-1)*np.outer(v,v))
    Y = np.column_stack((bare,np.zeros(3)))
    nd = np.append(math.sqrt(1-ad*ad)*v,ad)
    ne = np.append(-math.sqrt(1-ae*ae)*v,ae)
    le = 1.5*math.sqrt((1-ad*ad)/(1-ae*ae))
    m = (3*nd+2*le*ne)/5
    vlam = (nd-le*ne)/5
    assert np.linalg.norm(Y@m) < 1e-12
    for col, norm, scale, wanted in [(m+2*vlam,1,1,d),(m-3*vlam,le,k,e)]:
        assert math.isclose(np.linalg.norm(col),norm)
        matched, _ = light_spectrum(Y,col/norm)
        assert np.allclose(scale*matched,wanted,atol=3e-13,rtol=3e-9)
    return {"k": k, "relative_filter_e_over_d": rho, "alpha_d":ad,"alpha_e":ae}


def certificate():
    a = read(HERE/"inputs.json")
    data = np.array(read(HERE/a["same_scale_inputs"])["scales"][a["scale_key"]])
    d,e = np.split(data[:,0],2)
    central = one_direction(d,e)
    corners = []
    for indices in itertools.product((-1,1),repeat=6):
        ds,es = np.split(data[:,0]+a["halfwidths"]*data[:,1]*indices,2)
        corners.append(one_direction(ds,es))
    outer = {key:[min(v[key][0] for v in corners),max(v[key][1] for v in corners)]
             for key in ("k_for_growing","k_for_shrinking")}
    cost = {key:min(v[key] for v in corners)
            for key in ("minimum_strain_growing","minimum_strain_shrinking")}
    assert min(cost.values()) > 24

    parent = np.maximum(d,e)
    Y = np.column_stack((np.diag(parent),np.zeros(3)))
    nd,ne = unit_heavy(parent,d),unit_heavy(parent,e)
    ld,le = a["heavy_norm_d_and_e_for_witness"]
    m, vl = (3*ld*nd+2*le*ne)/5,(ld*nd-le*ne)/5
    assert np.allclose(m+2*vl,ld*nd) and np.allclose(m-3*vl,le*ne)
    sd,Ud = light_spectrum(Y,nd)
    se,Ue = light_spectrum(Y,ne)
    assert np.allclose(sd,d,rtol=5e-12) and np.allclose(se,e,rtol=5e-12)
    overlap = abs(np.vdot(nd,ne))
    gauge_sigmas = np.linalg.svd(Ud.conj().T @ Ue,compute_uv=False)
    assert np.allclose(gauge_sigmas,[1,1,overlap],atol=4e-13)
    assert abs(np.vdot(nd[:3],ne[:3])) < 1e-12
    e_upper = float(data[3,0]+2*data[3,1])
    d_lower = float(data[0,0]-2*data[0,1])
    s_upper = float(data[1,0]+2*data[1,1])
    mu_lower = float(data[4,0]-2*data[4,1])
    return {"central":central, "outer_k_intervals_two_widths":outer,
            "minimum_filter_strain_across_two_widths":cost,
            "test_without_using_electron": {"mu_lower_over_s_upper":mu_lower/s_upper,
                    "electron_predicted_lower":d_lower,"observed_electron_upper":e_upper,
                    "missing_yukawa":d_lower-e_upper,"gap_factor":d_lower/e_upper},
            "pure_offdiagonal_adjoint_max_filter_strain":1.5,
            "collinear_examples_with_extra_k":[aligned_witness(d,e,k) for k in (.3,5.)],
            "two_direction_witness": {"parent":parent.tolist(),"h_d_unit":nd.tolist(),
                    "h_e_unit":ne.tolist(),"m":m.tolist(),"v_lambda":vl.tolist(),
                    "light_X_singular_values":[1.0,1.0,float(overlap)],
                    "loss_in_third_squared_coupling":float(1-overlap**2)}}


def compare(x,y):
    if isinstance(x,dict):
        assert x.keys()==y.keys()
        for key in x:
            compare(x[key],y[key])
    elif isinstance(x,list):
        assert len(x)==len(y)
        for a,b in zip(x,y):
            compare(a,b)
    else:
        assert math.isclose(x,y,abs_tol=1e-14,rel_tol=5e-10),(x,y)


if __name__ == "__main__":
    cli = argparse.ArgumentParser()
    mode = cli.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write",action="store_true")
    mode.add_argument("--check",action="store_true")
    args = cli.parse_args()
    result = certificate()
    path = HERE/"results.json"
    if args.write:
        with path.open("x",encoding="utf-8") as out:
            json.dump(result,out,indent=2,allow_nan=False)
            out.write("\n")
    else:
        compare(result,read(path))
        print("Round 1514 spectra, matching and gauge overlaps passed.")
