"""Independent local-potential and finite-alignment checks (NumPy only)."""
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np
from numpy.polynomial import Polynomial


HERE = Path(__file__).resolve().parent
V = 246.22
MH = 125.25
HEAVY = 350.


def quartics(c, k3, k4, b_shift=0., singlet_fraction=1.):
    """Return coefficients of a globally nonnegative renormalizable potential.

    singlet_fraction selects what share of the free k4 margin stays in eps*S^4.
    The sign of b_shift selects the heavy-square amplitude branch.
    """
    s = math.sqrt((1-c)*(1+c))
    a = (k3-c**3)/(2*V*s*s)
    b = c*c/(2*V*s)
    excess = k4-k3*k3
    if excess < 0:
        raise ValueError("potential would fall below the specified vacuum")
    if b_shift:
        hhh_heavy_row = math.copysign(MH/HEAVY/V*math.sqrt((1-singlet_fraction)*excess), b_shift)
        b += hhh_heavy_row/(2*s*s)
    eps = (MH/V)**2 * excess*singlet_fraction/(8*s**4)
    return a,b,eps


def ray(c, coef, direction):
    s = math.sqrt((1-c)*(1+c))
    t = Polynomial([0,1])
    radial, singlet = direction
    x = radial*t + (radial*t)**2/(2*V)
    y = singlet*t
    a,b,eps = coef
    light_row = c*x+s*y+a*y*y
    heavy_row = -s*x+c*y+b*y*y
    return MH**2/2*light_row**2 + HEAVY**2/2*heavy_row**2 + eps*y**4


def example(c,k3,k4,sign=0):
    frac = 1. if sign == 0 else 0.
    coef = quartics(c,k3,k4,sign,frac)
    s = math.sqrt((1-c)*(1+c))
    h, p = np.array([c,s]),np.array([-s,c])
    hs, ps = ray(c,coef,h),ray(c,coef,p)
    plus,minus = ray(c,coef,h+p),ray(c,coef,h-p)
    gphh = (plus.deriv(3)(0)-minus.deriv(3)(0)-2*ps.deriv(3)(0))/6
    expected = 2*MH**2*c*(k3-c)/(V*s)
    expected += sign*MH*HEAVY/V*math.sqrt(max(0,k4-k3**2))
    assert abs(gphh-expected) < 2e-7
    g3, g4 = hs.deriv(3)(0), hs.deriv(4)(0)
    assert np.allclose([g3*V/(3*MH**2),g4*V**2/(3*MH**2)],[k3,k4])
    mp = math.sqrt(ps.deriv(2)(0))
    mh = math.sqrt(hs.deriv(2)(0))
    mix_mass = (plus.deriv(2)(0)-minus.deriv(2)(0))/4
    assert abs(mp-HEAVY)<1e-8 and abs(mh-MH)<1e-8 and abs(mix_mass)<1e-8
    assert coef[2]>=0
    singlet_quartic = ray(c,coef,[0.,1.]).deriv(4)(0)
    width = gphh**2/(32*math.pi*HEAVY)*math.sqrt(1-(2*MH/HEAVY)**2)
    return {"c":c,"k3":k3,"k4":k4,"heavy_square_sign":sign,"A_inverse_GeV":coef[0],
            "B_inverse_GeV":coef[1],"epsilon":coef[2],"G_SSSS":singlet_quartic,
            "G_Hhh_GeV":gphh,"tree_partial_width_H_to_hh_GeV":width}


def capped(c, cap):
    """Rigorous all-doublet/singlet bound and sharper singlet-only bound."""
    s = math.sqrt((1-c)*(1+c))
    max_mu = (cap*V**2/3-c*c*MH**2)/(s*s)
    if max_mu < 0:
        raise ValueError("radial mass sum already exceeds the cap")
    cubic_bound = math.sqrt(3*cap*max(MH**2,max_mu))
    # 0<=Q_4<=cap evaluated at 2*n +/- d and n +/- d gives |Q_nnnd|<=11cap/16.
    # |T_nuu| <= (4sqrt(2)+2) K/6 using polarization of n +/- u.
    qmixed = 11*cap/16
    tmixed = (4*math.sqrt(2)+2)*cubic_bound/6
    all_err = (2*c*c*s*V*qmixed + 3*c*s*s*tmixed + s**3*cubic_bound)*V/(3*MH**2)
    singlet_m2 = s*s*MH**2+c*c*max_mu
    singlet_err = (3*c*s*s*V*(2*cap/3)+s**3*math.sqrt(3*cap*singlet_m2))*V/(3*MH**2)
    return {"c":c,"directional_quartic_cap":cap,"all_k3_range":[c**3-all_err,c**3+all_err],
            "one_doublet_any_singlets_k3_range":[c**3-singlet_err,c**3+singlet_err]}


def exact_mixed_test():
    """Use r-polynomials directly for one transverse doublet and a real singlet."""
    doublet_mix, singlet_mix = 320., -920.
    r = Polynomial([V,1.])
    source_d = doublet_mix/(2*V**2) * r*(r*r-V**2)
    source_s = singlet_mix/(2*V)*(r*r-V**2)
    np.testing.assert_allclose([source_d(0),source_s(0)],[0.,0.],atol=1e-9)
    np.testing.assert_allclose([source_d.deriv()(0),source_s.deriv()(0)],[doublet_mix,singlet_mix])
    return [V*source_d.deriv(2)(0)/doublet_mix, V*source_s.deriv(2)(0)/singlet_mix]


def heavy_threshold(c):
    s = math.sqrt((1-c)*(1+c))
    delta = MH**2/4
    mass2 = MH**2+delta/s**2
    cap = 3*(MH**2+delta)/V**2
    k4 = c**4*(1+delta/MH**2)
    heavy_vertex = -s*c**2*(mass2+2*MH**2)/V
    low = k4 - V**2*heavy_vertex**2/(MH**2*mass2)
    np.testing.assert_allclose(low,c**4*(1-4*s*s*(1+MH**2/mass2)-s*s),atol=1e-11)
    # The exact three heavy propagators for hh->hh at physical threshold.
    threshold = k4 - V**2*heavy_vertex**2/(3*MH**2)*(1/(mass2-4*MH**2)+2/mass2)
    return {"c":c,"M_GeV":math.sqrt(mass2),"all_direction_quartic_cap":cap,"k3":c**3,
            "k4_full_contact":k4,"k4_low_local":low,"k4_plus_heavy_exchange_at_hh_threshold":threshold,
            "G_Hhh_GeV":heavy_vertex}


def invisible():
    mass_s = 40.
    visible_width = 4.1e-3
    phase = math.sqrt(1-4*mass_s**2/MH**2)
    gamma_per_portal2 = V**2/(32*math.pi*MH)*phase
    bound = .107
    return {"mS_GeV":mass_s,"assumed_visible_width_GeV":visible_width,"branching_upper":bound,
            "portal_abs_max":math.sqrt(bound/(1-bound)*visible_width/gamma_per_portal2),
            "at_portal_0p01_BR":1/(1+visible_width/(.01**2*gamma_per_portal2))}


def round_tree(x):
    if isinstance(x,dict):
        return {k:round_tree(v) for k,v in x.items()}
    if isinstance(x,list):
        return [round_tree(v) for v in x]
    return float(f"{x:.12g}") if isinstance(x,(float,np.floating)) else x


if __name__=="__main__":
    out=round_tree({"round":1790,"v_GeV":V,"light_mass_GeV":MH,"heavy_mass_GeV":HEAVY,
                    "mixed_Ward_factors":exact_mixed_test(),
                    "stable_witnesses":[example(c,k3,k4,sgn) for c,k3,k4,sgn in
                                        [(0.99,1,1,0),(0.99,2,4,0),(0.9999,2,4,0),
                                         (0.99,1,2,-1),(0.99,1,2,1),(0.99,1,2,0)]],
                    "finite_ranges":[capped(c,cap) for cap in (1.,4*math.pi) for c in (.99,.999,.9999)],
                    "threshold_with_bounded_quartics":[heavy_threshold(c) for c in (.99,.999,.9999)],
                    "exact_alignment_invisible_example":invisible(),
                    "experimental_slices":{"ATLAS_CMS_HH_k3_95":[-.71,6.1],
                                           "ATLAS_HHH_k4_at_k3_1_95":[-230.,240.],
                                           "ATLAS_CMS_HH":"https://cms-results.web.cern.ch/cms-results/public-results/publications/HIG-25-014/",
                                           "ATLAS_HHH":"https://arxiv.org/abs/2411.02040"}})
    body=json.dumps(out,indent=2,ensure_ascii=False)+"\n"
    target=HERE/"independent_results.json"
    if "--check" in sys.argv:
        assert json.loads(target.read_text(encoding="utf-8")) == out
    else:
        target.write_text(body,encoding="utf-8")
    print("1790 independent", hashlib.sha256(body.encode()).hexdigest()[:12])
