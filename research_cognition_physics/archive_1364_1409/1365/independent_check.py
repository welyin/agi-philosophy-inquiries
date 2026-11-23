"""Round 1365: cut triangle quadrature, invariant scalar couplings and an EFT alias."""
import argparse
import json
from functools import lru_cache
from pathlib import Path

import numpy as np


@lru_cache(maxsize=None)
def quadrature(order):
    points, weights = np.polynomial.legendre.leggauss(order)
    return (points+1)/2, weights/2


def triangle(mass, q, order):
    """Integrate 3 (1-4xy)/(1-4*tau*xy-i0) over the Feynman triangle."""
    tau = (q/(2*mass))**2
    u, weights = quadrature(order)
    # Cubic endpoint mapping regularises a logarithmic singularity on either end.
    a, b = u**3, (1-u)**3
    shape, reverse = a/(a+b), b/(a+b)
    slope = 3*u**2*(1-u)**2/(a+b)**2
    endpoints = [0., 1.]
    if tau > 1:
        beta = np.sqrt(1 - 1/tau)
        lower = 1/(2*tau*(1+beta))
        endpoints = [0., lower, 1-lower, 1.]
    integral = 0j
    for lo, hi in zip(endpoints,endpoints[1:]):
        span = hi-lo
        x = lo+span*shape
        txy = 4*tau*x*(1-x)
        if tau < 1:
            real_log = np.log1p(-txy)
        else:
            real_log = np.full_like(x,np.log(4*tau))
            for root in endpoints[1:3]:
                if root == lo:
                    distance = span*shape
                elif root == hi:
                    distance = span*reverse
                else:
                    distance = abs(x-root)
                real_log += np.log(distance)
            # Near x=0 the small difference must be taken before dividing by x.
            small = txy < .5
            real_log[small] = np.log1p(-txy[small])
        cut = (np.pi if tau>1 and lo == endpoints[1] else 0.)
        integral += np.dot(weights,span*slope*(-real_log+1j*cut)/x)
    ans = 1.5/tau + 3*(1 - 1/tau)*integral/(4*tau)
    return ans


def virtual_rates(up_masses, down_masses, kappas, q, order):
    up = sum(triangle(m,q,order) for m in up_masses)
    down = np.array([triangle(m,q,order) for m in down_masses])
    reference = up+sum(down[:3])
    return (up+down@kappas)/reference


def run():
    # Both spectra are certified independently in the previous round; the two 8-decimal
    # rational endpoints suffice to bracket changes in the second branch below the recorded accuracy.
    mass = np.array([[1,0,0,0],[-1,2,0,0],[1,-1,3,0],
                     [-1+2j,-2+1j,-2+2j,4]],complex)
    squares, left = np.linalg.eigh(mass@mass.conj().T)
    masses = np.sqrt(squares)
    right = (mass.conj().T@left)/masses
    assert np.allclose(right.conj().T@right,np.eye(4))
    vm = np.diag([1,1,1,0])@mass
    k = np.diag(left.conj().T@vm@right)/masses
    assert np.allclose(k.imag,0)
    k = k.real
    assert np.isclose(sum(k),3)
    values = []
    up = [.5,2.,9.]
    other_mass = 0.
    for heavy in [squares[-1],42.81815218,42.81815219]:
        others = masses.copy()
        others[-1] = np.sqrt(heavy)
        amp64 = virtual_rates(up,others,k,8,64)
        amp128 = virtual_rates(up,others,k,8,128)
        assert abs(amp128-amp64) < 2.e-9
        values.append(float(abs(amp128)**2))
        other_mass = others[-1]
    assert np.ptp(values[1:]) < 2.e-11
    assert np.allclose(values[:2],[.9772179933832474,.9766236468033801],rtol=0,atol=2.e-10)
    # First alternative in the admission uses only the 44 element of the derivative.
    changed = vm.copy()
    changed[3,3] += 1
    shifted_k = np.diag(left.conj().T@changed@right)/masses
    assert abs(shifted_k.imag).max() < 2.e-12
    simple_gg = abs(virtual_rates(up,masses,shifted_k.real,8,128))**2
    # Broader gauge-allowed lower-row derivatives. Design the row in the gauge basis
    # and diagonalise it back; off-diagonal scalar vertices survive but do not enter this loop.
    delta_k = k[-1]*(triangle(other_mass,8,128)/triangle(masses[-1],8,128)-1)
    assert abs(delta_k.imag) < 1.e-13
    lower = delta_k.real*masses[-1]*right[:,-1].conj()/left[-1,-1].conj()
    assert np.isclose(lower[-1],4*delta_k.real)
    broadened = vm.copy()
    broadened[-1] += lower
    broadened_k = np.diag(left.conj().T@broadened@right)/masses
    assert np.allclose(broadened_k-k,[0,0,0,delta_k])
    alias_amp = virtual_rates(up,masses,broadened_k.real,8,128)
    target = masses.copy()
    target[-1] = other_mass
    target_amp = virtual_rates(up,target,k,8,128)
    assert abs(alias_amp-target_amp)<2.e-13
    # A single fixed retuning is tied to the measured scalar invariant mass.
    other_invariant = abs(virtual_rates(up,masses,broadened_k.real,6,128))**2
    target_other = abs(virtual_rates(up,target,k,6,128))**2
    assert abs(other_invariant-target_other)>1.e-5
    # Added after inspecting the registered example: two charge conjugate flavour-changing
    # two-body final states. Compare against the same synthetic SM b-bbar partial width.
    outside = []
    light, bottom = masses[1], masses[2]
    coupling_square = (1-k[1])*(1-k[-1])
    for heavy in [masses[-1],other_mass]:
        q = 8.
        threshold = light+heavy
        ratio = 0.
        if threshold < q:
            a, b = (light/q)**2, (heavy/q)**2
            momentum_factor = np.sqrt((1-a-b)**2-4*a*b)
            spin_trace = (light**2+heavy**2)*(1-a-b)-4*(light*heavy/q)**2
            ratio = momentum_factor*coupling_square*spin_trace/(bottom**2*(1-4*(bottom/q)**2)**1.5)
        outside.append({"threshold":float(threshold),"partial_width_over_reference_bb":float(ratio)})
    assert np.isclose(outside[0]["partial_width_over_reference_bb"],.0204960089679)
    assert outside[1]["partial_width_over_reference_bb"]==0.
    return {"kinematic_reference":"same synthetic six masses and canonical diagonal doublet slopes",
            "bb_width_ratio":float(k[2]**2), "gg_width_ratio":values[:2],
            "exclusive_equal_halfwidth":float(abs(values[0]-values[1])/2),
            "single_mass_slope_competitor":{"bb":float(abs(shifted_k[2])**2),"gg":float(simple_gg)},
            "posthoc_flavour_changing_tree_threshold":outside,
            "matched_both_on_shell":{"heavy_kappa_shift":float(delta_k.real),
                                    "bb":float(abs(broadened_k[2])**2),"gg":float(abs(alias_amp)**2),
                                    "lower_row_real":lower.real.tolist(),"lower_row_imag":lower.imag.tolist(),
                                    "q6_alias":float(other_invariant),"q6_target":float(target_other)}}


def equal(data,ref):
    if isinstance(data,dict):
        assert data.keys() == ref.keys()
        for key,val in data.items():
            equal(val,ref[key])
    elif isinstance(data,list) and data and isinstance(data[0],dict):
        assert len(data) == len(ref)
        for val, target in zip(data,ref):
            equal(val,target)
    elif isinstance(data,str):
        assert data==ref
    else:
        assert np.allclose(data,ref,atol=2.e-10,rtol=0)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--write",action="store_true")
    args=parser.parse_args()
    output=run()
    filename=Path(__file__).with_name("independent_results.json")
    if args.write:
        with filename.open("x",encoding="utf-8") as stream:
            json.dump(output,stream,allow_nan=False,indent=2)
            stream.write("\n")
    else:
        equal(output,json.loads(filename.read_text(encoding="utf-8")))
    print("PASS",output)
