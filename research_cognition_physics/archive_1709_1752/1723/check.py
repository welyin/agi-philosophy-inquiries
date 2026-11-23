"""Rank cost of an adjoint-selected heavy-matter port and its radial flavour rates."""
import argparse
import json
import math
import runpy
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent
CFG=json.loads((HERE/"inputs.json").read_bytes())


def code(key):
    return runpy.run_path(str(HERE/CFG[key]))


def sv(a):
    return np.linalg.svd(a,compute_uv=False)


def one_path_floor(d,e):
    # Eckart–Young gives sigma_2 exactly at a fixed complex z. Cover the
    # finite square using Weyl Lipschitz bounds; tail is controlled separately.
    radius=CFG["pencil_tail_radius"]
    lower=CFG["certified_one_path_error"]
    guard=2e-12
    ds,es=sv(d),sv(e)
    assert radius*ds[1]-es[0]>lower
    centers=np.zeros(1,dtype=complex)
    centers_examined=0
    width=float(radius)
    for depth in range(35):
        centers_examined+=len(centers)
        s2=sv(e-centers[:,None,None]*d)[:,1]
        centers=centers[s2-(math.sqrt(2)*width*ds[0]+guard)<=lower]
        if not len(centers):
            return {"lower_bound":lower, "cells":centers_examined,"max_depth":depth,
                    "tail_lower_bound":radius*ds[1]-es[0]}
        width/=2
        centers=(centers[:,None]+width*np.array([1+1j,1-1j,-1+1j,-1-1j])).ravel()
    raise ValueError("certificate is too strong for the actual flavour point")


def main():
    frame=code("charged_frames")
    _,rep=frame["charged_background"](.10)
    d=frame["unpack"](rep["down_texture_10PeV"])
    e=frame["unpack"](rep["charged_texture_10PeV"]).T
    z0=e[0,1]/d[0,1]
    zt=e[2,2]/d[2,2]
    re=np.array([[0,z0,0],[z0,3*z0,zt],[0,zt,zt]])
    np.testing.assert_allclose(e,re*d,atol=1e-15)
    # With the first 1,2 rows and columns the minor is -D12 D21(z0-z)^2.
    # At z0 the independent lower-right 2x2 minor survives.
    minor=np.linalg.det((e-z0*d)[1:,1:])
    assert abs(minor)>5e-5
    # A pencil rank drop has only these determinant roots, including the
    # double root z0. Neither cancels an arbitrary two-dimensional subspace.
    for z in (0,1,1j,2+.2j,z0,zt):
        np.testing.assert_allclose(np.linalg.det(e-z*d),np.linalg.det(d)*(z0-z)**2*(zt-z),
                                   atol=2e-20,rtol=1e-10)

    color=code("five_path_exchange")
    R=json.loads((HERE/CFG["register"]).read_bytes())["central_conditional_combined_point"]["R"]
    fraction=math.sqrt(1-R**-2)
    l,ds,rh=np.linalg.svd(d)
    u=(rh.conj().T*[-1,1,-1])@l.conj().T
    a,w=color["dilation"](d@u,ds[0])
    embedding=a@w@u.conj().T
    common=ds[0]/fraction*w.conj().T
    np.testing.assert_allclose(common@embedding*fraction,d,atol=1e-16)
    le,es,reh=np.linalg.svd(e)
    def off(a):
        return a-np.diag(np.diag(a))
    ed=off(l.conj().T@e@rh.conj().T)
    de=off(le.conj().T@d@reh.conj().T)
    pair_d=sum(abs(ed[i,j])**2 for i,j in ((0,1),(1,0)))
    pair_e=sum(abs(de[i,j])**2 for i,j in ((0,1),(1,0)))
    cases={}
    for name,z in (("common",1.),("free_real",z0),("second_complex_root",zt)):
        delta=e/z-d
        # Minimum operator-norm lift to the five heavy directions. The full light
        # vector is (I/R, fraction*embedding), and the new vertex has zero old block.
        add=delta@embedding.conj().T/fraction
        sl=np.linalg.svd(add,compute_uv=False)
        rank=np.linalg.matrix_rank(delta,tol=1e-13)
        assert rank==(3 if name=="common" else 2)
        np.testing.assert_allclose(z*(common+add)@embedding*fraction,e,atol=1e-15)
        np.testing.assert_allclose(sl[:3],sv(delta)/fraction,atol=1e-15)
        def read_pair(x):
            return ((common+(2-2*(1+x))/5*add)@embedding*fraction,
                    z*(common+(2+3*(1+x))/5*add)@embedding*fraction)
        dplus,eplus=read_pair(.01)
        dminus,eminus=read_pair(-.01)
        gd=l.conj().T@(dplus-dminus)@rh.conj().T/.02
        ge=le.conj().T@(eplus-eminus)@reh.conj().T/.02
        np.testing.assert_allclose(off(gd),-2/(5*z)*ed,atol=1e-14)
        np.testing.assert_allclose(off(ge),-3*z/5*de,atol=1e-14)
        cases[name]={"z":[z.real,z.imag],"rank":int(rank),
                     "lift_singular_values":sl,"best_one_factor_absolute_error_on_E":sv(e-z*d)[1],
                     "fraction_tagged_s_b_over_tau_mu":3*4*pair_d/(9*abs(z)**4*pair_e)}
    cert=one_path_floor(d,e)
    return {"original_10PeV_masses_down_tau_desc":[ds,es],"nonzero_minor_at_double_root":abs(minor),
            "five_path_tree_match":{"R":R,"required_common_vertex_norm_at_this_boundary":ds[0]/fraction,
                "signed_up_return_over_yd_here":(d@u)[0,0].real/ds[-1]},
            "uniform_one_factor_error":cert,"exact_and_conditional_implementations":cases,
            "radial_flavour_test":{
                "fraction_bs_among_all_flavor_changing_down":pair_d/np.linalg.norm(ed)**2,
                "fraction_tau_mu_among_all_flavor_changing_charged_leptons":pair_e/np.linalg.norm(de)**2,
                "geometric_mean_of_dimensionless_offdiag_norms":math.sqrt(6*np.linalg.norm(ed)*np.linalg.norm(de))/5}}


def tidy(x):
    if isinstance(x,dict): return {k:tidy(v) for k,v in x.items()}
    if isinstance(x,(list,tuple,np.ndarray)): return [tidy(v) for v in x]
    if isinstance(x,(np.floating,float)): return float(f"{x:.11g}")
    return x


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    check=p.parse_args().check
    output=tidy(main())
    target=HERE/"results.json"
    if check:
        assert json.loads(target.read_bytes())==output
        print("1723 verified the full complex pencil, five-path matching and radial flavor rates")
    else:
        with target.open("x",encoding="utf8") as out:
            json.dump(output,out,indent=2,allow_nan=False)
            out.write("\n")
