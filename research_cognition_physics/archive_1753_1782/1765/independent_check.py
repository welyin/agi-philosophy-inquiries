"""Reconstruct the db triangle without any bottom semileptonic decay or flavour texture."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from numpy.polynomial import Polynomial as Poly


HERE = Path(__file__).resolve().parent
STAGE = HERE.parent


def polynomial_test(xi, us, angle, clock, fixed_beta=False):
    md, ms = clock["delta_md"][0], clock["delta_ms"][0]
    # Standard common top box, absolute coefficient only for the strange meson.
    b2 = ms / 10180.
    a2 = b2 * xi**2 * (md/ms) * clock["mass_s"][0]/clock["mass_d"][0]
    discriminant = 1-4*(a2+b2)
    small = 2*(a2+b2)/(1+math.sqrt(discriminant))
    p = us**2
    cosgamma = math.cos(math.radians(angle))
    branches = []
    for u in (small, 1-small):
        w = 1-u
        td2 = a2/u
        x = Poly([0., w])  # rescaling keeps small physical roots well-conditioned
        left = x*(1-p-x)
        middle = (p+x-td2)*(w-x)
        total = left+middle-a2
        beta_side = middle+a2-left
        equation = (beta_side**2-4*cosgamma**2*a2*middle if fixed_beta
                    else total**2-4*cosgamma**2*left*middle)
        rows = []
        for root in equation.roots():
            if abs(root.imag) > 1e-8:
                continue
            xx = w*root.real
            if min(xx, p+xx-td2, 1-p-xx, w-xx) <= 0:
                continue
            signed = beta_side(root.real) if fixed_beta else total(root.real)
            if signed*cosgamma <= 0:
                continue
            cos_gamma_here = total(root.real)/(2*math.sqrt(left(root.real)*middle(root.real)))
            if abs(cos_gamma_here) > 1+1e-11:
                continue
            gamma_here = math.acos(np.clip(cos_gamma_here,-1.,1.)) if fixed_beta else math.radians(angle)
            d = np.sqrt([1-p-xx,p+xx-td2,td2])
            bu = math.sqrt(xx)*np.exp(-1j*gamma_here)
            bc = -math.sqrt(w-xx)
            bt = -(d[0]*bu+d[1]*bc)/d[2]
            b = np.array([bu,bc,bt])
            v = np.column_stack([d,np.conj(np.cross(d,b)),b])
            np.testing.assert_allclose(v.conj().T@v,np.eye(3),atol=5e-10,rtol=0)
            np.testing.assert_allclose(abs(v[0,1]),us,atol=2e-10)
            np.testing.assert_allclose([abs(v[2,2]*v[2,0])**2,abs(v[2,2]*v[2,1])**2],
                                       [a2,b2],atol=2e-10,rtol=0)
            beta = np.angle(-v[1,0]*v[1,2].conjugate()/(v[2,0]*v[2,2].conjugate()))
            if fixed_beta:
                assert abs(math.degrees(beta)-angle) < 1e-6
            rows.append(dict(ub=abs(bu),cb=abs(bc),beta_deg=math.degrees(beta),
                             golden_S=math.sin(2*beta),gamma_deg=math.degrees(gamma_here)))
        # Exact stationary points of cos(gamma) when all three triangle sides are nonzero.
        deriv = 2*total.deriv()*left*middle-total*(left*middle).deriv()
        gamma_max = None
        if u > .5:
            vals = []
            for y in deriv.roots():
                if abs(y.imag) < 1e-8 and 0 < y.real < 1 and left(y.real)*middle(y.real)>0:
                    val = total(y.real)/(2*math.sqrt(left(y.real)*middle(y.real)))
                    if -1 <= val <= 1:
                        vals.append(math.degrees(math.acos(val)))
            # For this large-u branch L or M vanish at either end with total>0.
            assert total(0) > 0 and total(1) > 0 and len(vals) == 1
            gamma_max = vals[0]
        branches.append(dict(tb2=u, allowed_by_CMS_top=math.sqrt(u)>.975,
                             physical_solutions=rows, max_gamma_deg=gamma_max))
    return branches


def run():
    clock = json.loads((STAGE/"1736/independent_results.json").read_bytes())["data"]
    us = json.loads((STAGE/"1755/inputs.json").read_bytes())["holdout_Kl3"][0]
    points, golden = {}, {}
    beta = math.degrees(math.asin(clock["golden_S"][0]))/2
    for source in ("HPQCD","FNAL_MILC"):
        xi = clock[source][0]
        points[source] = {str(g): polynomial_test(xi, us, g, clock)
                          for g in (60.6, 63.6, 66.6)}
        golden[source] = {}
        for label, value in (("central",xi),("low",xi-math.hypot(*clock[source][1:])),
                             ("high",xi+math.hypot(*clock[source][1:]))):
            # Show sensitivity to one lattice input, with all other data held fixed.
            golden[source][label] = {str(b):polynomial_test(value,us,b,clock,True)[1]["physical_solutions"]
                                     for b in (beta,90-beta)}
    return dict(gamma_central_and_68_deg=[63.6,60.6,66.6], kaon_us=us,
                theta_calibration="No bottom semileptonic. Additional observed Δmd and K_l3 used here.",
                branches=points, golden_calibration=golden)


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=json.loads(json.dumps(run()),parse_float=lambda x:float(f"{float(x):.9g}"))
    file=HERE/"independent_results.json"
    if args.check:
        assert result==json.loads(file.read_bytes())
        print("1765 independent: db triangle inverse branches and two B clocks verified")
    else:
        with file.open("w",encoding="utf-8") as output:
            json.dump(result,output,indent=2,allow_nan=False)
            output.write("\n")
