"""Flavour-changing top width required by a single protected c,t mass deformation."""
import argparse
import json
import math
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE/name).read_text(encoding="utf-8-sig"))


def width(top,charm,h,left_to_right,right_to_left):
    lam=math.sqrt((top**2-(h+charm)**2)*(top**2-(h-charm)**2))
    norm=abs(left_to_right)**2+abs(right_to_left)**2
    interference=(left_to_right*right_to_left).real
    return lam*((top**2+charm**2-h**2)*norm+4*top*charm*interference)/(32*math.pi*top**3)


def exact_floor(top,charm,h,normalized_product):
    return (math.sqrt((top**2-(h+charm)**2)*(top**2-(h-charm)**2))*
            ((top+charm)**2-h**2)*normalized_product/(16*math.pi*top**3))


def solve():
    cfg=read("inputs.json")
    masses=read(cfg["mass_calibration"])["one_TeV_DR"]["up"]
    (yc,ec),(yt,et)=masses[1:]
    v=cfg["v_GeV"]
    h=cfg["mH_GeV"]
    low,top,high=cfg["kinetic_top_mass_GeV"]
    _,charm,chigh=cfg["kinetic_charm_mass_GeV"]
    prod=yc*yt/2
    expected=exact_floor(top,charm,h,prod)
    gs=math.sqrt(prod)
    np.testing.assert_allclose(expected,width(top,charm,h,gs,gs),rtol=1e-15)
    # Explicit Wilson term at the physical vacuum, with rank 1 and zero mass-trace.
    mf=v*np.array([yc,yt])/math.sqrt(2)
    m=np.diag(mf)
    x=np.array([[-1j*mf[0],gs*v],[gs*v,1j*mf[1]]])
    assert np.linalg.matrix_rank(x,tol=1e-12)==1
    assert abs(np.trace(np.linalg.solve(m,x)))<2e-14
    for r in (.3,.78,1.,1.15):
        field=r*(m+(r**2-1)*x/2)
        np.testing.assert_allclose(np.linalg.det(field),r**2*np.linalg.det(m),rtol=2e-13)
    # Test both interference signs, chiral ratios and CP-even diagonal shifts.
    for a in (-3,-.5,0,.5,3):
        for b in (-3,-.5,0,.5,3):
            gp=prod*(a-1j)*(b+1j)
            for ratio in (.25,1.,4):
                g=math.sqrt(abs(gp))*ratio
                assert width(top,charm,h,g,gp/g)>=expected*(1-2e-15)
    # A protected rank-two diagonal deformation permits pseudoscalar slopes without FCNC.
    z=.1j
    slopes=np.array([3-2/(1+z),3-2/(1+z.conjugate())])
    assert abs(sum(slopes.imag))<1e-15 and abs(slopes[0].imag)>.19
    outer_prod=(yc-2*ec)*(yt-2*et)/2
    kinematic_floor=(low**2-(h+chigh)**2)*(low**2-h**2)/(16*math.pi*high**3)
    outer=kinematic_floor*outer_prod*cfg["matching_and_rate_response_floor"]
    reference=cfg["total_width_GeV"]["reference"]
    wmax=cfg["total_width_GeV"]["outer"]
    # Allow a measured soft-slope mismatch d and a second normalized flavour port.
    # These toy data check the inequality without assuming exact factorization.
    for a,b in ((0.,0.),(1.4,-.6),(-.7,-.5)):
        eta,slack=1.2,.05
        diagonal_product=(a-1j*(eta-slack))*(b+1j*eta)
        for determinant in (.1, .03+.09j, -.07j):
            gp=prod*(diagonal_product-determinant)
            g=math.sqrt(abs(gp))
            safe=expected*max(0,eta*(eta-slack)-abs(determinant))
            assert width(top,charm,h,g,gp/g)>=safe*(1-2e-15)
    return dict(central_GeV_per_odd_squared=expected,central_BR_per_odd_squared=expected/reference,
                conservative_visible_GeV_per_odd_squared=outer,reference_width_GeV=reference,
                conservative_visible_BR_per_odd_squared=outer/wmax,
                conditional_odd_limits={key:{"reference":math.sqrt(cfg[key]*reference/expected),
                     "outer":math.sqrt(cfg[key]*wmax/outer)} for key in
                     ("atlas_combined_tHc_BR_95","cms_combined_tHc_BR_95")},
                saturated_abs_chiral_coupling_at_unit_odd=gs,
                normalized_determinant_required_for_eta_1p2_soft_slack_0p05=max(0.,
                    1.2*(1.2-.05)-cfg["atlas_combined_tHc_BR_95"]*wmax/outer),
                protected_rank_two_odd_slopes=slopes.imag.tolist())


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--check",action="store_true")
    group.add_argument("--write-new",action="store_true")
    opt=parser.parse_args()
    result=solve()
    if opt.write_new:
        with (HERE/"results.json").open("x",encoding="utf8") as f:
            json.dump(result,f,indent=2)
            f.write("\n")
    else:
        def compare(a,b):
            if isinstance(a,dict):
                assert a.keys()==b.keys()
                for key in a:
                    compare(a[key],b[key])
            else:
                np.testing.assert_allclose(a,b,rtol=1e-12,atol=1e-15)
        compare(result,read("results.json"))
        print("1535: sharp chiral interference, rank-one saturation, rate and measurement transfer verified.")
