"""Spectral-norm certificates and explicit Majorana competitors for the equal-edge seed."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
# Independent extraction of the registered two marginal widths.
RATIO_MAX = (7.54+2*.19)*1e-5/((2.51-2*.10)*1e-3)
AMP_MAX = .0851+2*.0024
POLE_R = (.1056583755**2-.00051099895**2)/(1.77686**2-.00051099895**2)


def seed(r):
    return np.array([[0,r,-1],[r,0,1-r],[-1,1-r,0]],dtype=float)


def certificates(r):
    lam, vec = np.linalg.eigh(seed(r))
    i,k,j = np.argsort(abs(lam))
    alpha,beta = abs(lam[[i,k]]/lam[j])
    p0 = vec[0,i]**2
    no = (beta**2-alpha**2-RATIO_MAX*(1-beta**2))/(2*(beta+alpha)+2*RATIO_MAX*(1+beta))
    shifts = np.array([alpha+no,beta-no,1+no])
    np.testing.assert_allclose((shifts[1]**2-shifts[0]**2)/(shifts[2]**2-shifts[1]**2), RATIO_MAX)
    corners = [(1+sgn*math.sqrt(1-AMP_MAX))/2 for sgn in [-1,1]]
    bounds = []
    for p in corners:
        cos_angle = math.sqrt(p0*p)+math.sqrt((1-p0)*(1-p))
        s2 = 1-cos_angle**2
        # The coefficient of sin²(chi) is positive below the cut.
        cut = (beta**2-alpha**2)/(2*(1-alpha))-alpha
        interior = math.sqrt(s2*((beta**2-1)+(1-alpha)**2*(1-s2)))
        bounds.append(min(cut,interior))
    return {"small_singular_ratio": alpha, "middle_singular_ratio": beta,
            "electron_weight_of_seed_light_state": p0, "NO_sharp_spectral_only": no,
            "IO_bounds_two_reactor_branches": bounds,
            "combined_lower": min(no,*bounds)}


def pmns():
    p = (1-math.sqrt(1-.0851))/2
    a, b, c = [math.asin(math.sqrt(v)) for v in [.307,p,.514]]
    ca,cb,cc,sa,sb,sc = map(float,[math.cos(a),math.cos(b),math.cos(c),
                                     math.sin(a),math.sin(b),math.sin(c)])
    return np.array([[ca*cb,sa*cb,sb],
                     [-sa*cc-ca*sc*sb,ca*cc-sa*sc*sb,sc*cb],
                     [sa*sc-ca*cc*sb,-ca*sc-sa*cc*sb,cc*cb]])


def witnesses():
    # A finite grid yields examples, never a lower bound. Phases +- signs of real symmetric M
    # correspond to admissible relative Majorana CP parities.
    mixing = pmns()
    source=seed(POLE_R)
    source=source/np.linalg.norm(source,2)
    b2, a2 = 7.54e-5, 2.51e-3
    out={}
    weights=np.linspace(.2,1.4,61)
    for ordering in ["NO","IO"]:
        best={"upper": 9.}
        for light in [0.,.01,.03,.04,.055,.075,.10,.2]:
            masses=np.sqrt((np.array([0,b2,b2+a2]) if ordering=="NO" else
                            np.array([2.48e-3,2.48e-3+b2,0]))+light**2)
            for signs in itertools.product([1,-1],repeat=3):
                for rows2 in itertools.product([1,-1],repeat=2):
                    rows=(1,*rows2)
                    mat=mixing*np.array(rows)[:,None]
                    target=(mat*(masses*np.array(signs)/max(masses)))@mat.T
                    err=np.max(abs(np.linalg.eigvalsh(weights[:,None,None]*target-source)),axis=1)
                    idx=np.argmin(err)
                    if err[idx] < best["upper"]:
                        best={"upper":float(err[idx]), "higgs_independent_light_mass_eV":light,
                              "target_largest_mass_over_seed_norm":float(weights[idx]),
                              "mass_signs":signs,"row_signs":rows}
        assert .45 < best["upper"] < 1
        out[ordering]=best
    return out


def run():
    ref=certificates(POLE_R)
    zero=certificates(0.)
    np.testing.assert_allclose(zero["NO_sharp_spectral_only"],1/(2+4*RATIO_MAX))
    np.testing.assert_allclose(zero["IO_bounds_two_reactor_branches"],math.sqrt(1-AMP_MAX)/2)
    rmax=.004
    # ||seed(r)-seed(0)|| / ||seed(0)|| = |r|.
    np.testing.assert_allclose(np.linalg.norm(seed(rmax)-seed(0),2),math.sqrt(2)*rmax)
    return {"pole_r":POLE_R,"ratio_upper":RATIO_MAX,"reactor_amp_upper":AMP_MAX,
            "finite_r_certificates":ref,"rank_two_limit":zero,
            "conditional_rmax":rmax,
            "uniform_joint_lower_for_r_at_most_rmax":(zero["combined_lower"]-rmax)/(1+rmax),
            "compatible_mass_matrix_upper_witnesses":witnesses()}


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    mode=cli.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write-new",action="store_true")
    mode.add_argument("--check",action="store_true")
    opt=cli.parse_args()
    path=HERE/"independent_results.json"
    result=run()
    if opt.write_new:
        with path.open("x",encoding="utf-8") as f:
            json.dump(result,f,indent=2)
            f.write("\n")
    else:
        def check(a,b):
            if isinstance(a,dict):
                assert set(a)==set(b)
                for name in a:
                    check(a[name],b[name])
            else:
                np.testing.assert_allclose(a,b,rtol=1e-10,atol=1e-12)
        check(result,json.loads(path.read_text(encoding="utf-8")))
    print("1537 independent",result["finite_r_certificates"],result["compatible_mass_matrix_upper_witnesses"])
