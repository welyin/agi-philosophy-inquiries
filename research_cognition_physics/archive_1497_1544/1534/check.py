"""Radial strong angle of an affine Yukawa pencil and two-channel Higgs holdout."""
import argparse
import cmath
import json
import math
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE/name).read_text(encoding="utf-8-sig"))


def vertex(z):
    return 3-2/(1+z)


def pole_response(z, x):
    return np.prod([1+x*k for k in z])


def multiple_flavours():
    # General non-normal matrices with positive undeformed determinant.
    a=np.array([[.07+.19j,.15,-.11j],[0,.07-.19j,.16],[0,0,-.08]])
    transform=np.array([[1,.5j,-.2],[.2,1,0],[.1,.3,1]])
    a=transform@a@np.linalg.inv(transform)
    y=np.array([[.3,.1j,0],[0,.7,-.2],[0,0,1.3]])
    for radial in (.2,.7,1.,1.2):
        x=radial**2
        mass=radial*y@(np.eye(3)+x*a)
        deriv=y@(np.eye(3)+3*x*a)
        left,sv,vh=np.linalg.svd(mass)
        slopes=np.diag(left.conj().T@deriv@vh.conj().T)*radial/sv
        assert abs(np.linalg.det(mass).imag)<2e-14
        assert abs(slopes.imag.sum())<5e-15
        np.testing.assert_allclose(slopes.sum(),radial*np.trace(np.linalg.solve(mass,deriv)),
                                   rtol=1e-13,atol=1e-13)
    return np.poly(a).real.tolist()


def solve():
    cfg=read("inputs.json")
    edm=read(cfg["neutron"])["edm"]
    lattices=[edm["abs_slope_e_fm"],read(cfg["second_lattice"])["lattice_2026_abs_slope_e_fm"]]
    limits={}
    for year,(center,stat,systematic) in zip(("2023","2026"),lattices):
        theta=edm["limit_90_ecm"]/(edm["fm_in_cm"]*(center-2*(stat+systematic)))
        limits[year]={"local_theta_box":theta,"odd_top_if_even_one":2*math.tan(theta)}
    t=cfg["working_example_t"]
    k=cfg["finite_angle_competitor_k"]
    local=[1j*t,k*(1-1j*t)-1]
    pair=[1j*t,-1j*t]
    for z in (local,pair):
        np.testing.assert_allclose(pole_response(z,1).imag,0,atol=1e-15)
    d=sum(vertex(z).imag for z in local)
    step=1e-6
    approx=cmath.phase(pole_response(local,(1+step)**2)/pole_response(local,(1-step)**2))/(2*step)
    np.testing.assert_allclose(d,approx,rtol=1e-8)
    np.testing.assert_allclose([pole_response(pair,x).imag for x in (0,.1,1,2)],0)
    assert abs(pole_response(local,.5).imag)>1e-3
    a,b=.25+.14j,-.19-.06j
    p=(1+a)*(1+b)
    theta=cmath.phase(p)
    ratio=a/(1+a)+b/(1+b)
    im_p=abs(p)*math.sin(theta)
    im_der=abs(p)*(ratio.real*math.sin(theta)+ratio.imag*math.cos(theta))
    np.testing.assert_allclose([2*im_p-im_der,im_der-im_p],[(a+b).imag,(a*b).imag])
    mass=cfg["higgs_mass_GeV"]
    tau=cfg["tau_mass_GeV"]
    beta_tau_sq=1-(2*tau/mass)**2
    finite_max=beta_tau_sq/(1-(2*max(cfg["b_mass_work_range_GeV"])/mass)**2)
    atlas=read(cfg["atlas_rates"])["atlas_table7"]
    def rate_interval(name):
        (m,d,u),(reference,error)=atlas[name]["observed"],atlas[name]["SM"]
        return m/reference,(m-2*d)/(reference+2*error),(m+2*u)/(reference-2*error)
    rb,lb,ub=rate_interval("bb_ZZ")
    rl,ll,ul=rate_interval("tautau_ZZ")
    return dict(theta_transfer=limits,protected_ckm_independent_polynomial=multiple_flavours(),
                protected_example_vertex=[vertex(pair[0]).real,vertex(pair[0]).imag],
                one_vacuum_competitor_normalized_radial_angle_derivative=d,
                one_vacuum_competitor_theta_at_half_x=cmath.phase(pole_response(local,.5)),
                finite_mass_bottom_over_tau_normalized=[1.,finite_max],
                with_rate_transfer=[cfg["width_ratio_transfer"][0],finite_max*cfg["width_ratio_transfer"][1]],
                atlas_central_normalized_bb_over_tautau=rb/rl,
                atlas_two_width_outer=[lb/ul,ub/ll])


if __name__=="__main__":
    p=argparse.ArgumentParser()
    g=p.add_mutually_exclusive_group(required=True)
    g.add_argument("--check",action="store_true")
    g.add_argument("--write-new",action="store_true")
    arg=p.parse_args()
    result=solve()
    if arg.write_new:
        with (HERE/"results.json").open("x",encoding="utf8") as out:
            json.dump(result,out,indent=2)
            out.write("\n")
    else:
        def compare(a,b):
            if isinstance(a,dict):
                assert a.keys()==b.keys()
                for key in a:
                    compare(a[key],b[key])
            else:
                np.testing.assert_allclose(a,b,rtol=1e-10,atol=1e-16)
        compare(result,read("results.json"))
        print("1534: determinant phases, soft derivative, neutron and branching-ratio checks passed.")
