"""Charged-current normalization and independently accessible ordered rare-B channels."""
import argparse
import json
from math import hypot, pi, sqrt
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent


def unitary(rng):
    a=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
    return np.linalg.qr(a)[0]


def run():
    rng=np.random.default_rng(1434)
    eye=np.eye(3)
    for i in range(35):
        V,W,U=[unitary(rng) for _ in range(3)]
        j,b=i%3,(i//3)%3
        k=.027
        A=V[j,b]*U+k*np.outer(W[b].conj(),(V@W@U)[j])
        total=np.vdot(A,A).real
        assert np.isclose(total,(3+2*k)*abs(V[j,b])**2+k*k)
        # A neutral threshold across experiments changes the inclusive transfer.
        cut=eye[:,:2]
        assert np.linalg.norm(A@cut)**2 <= total+1e-13
        x,y=W[b],W[(b+1)%3]
        raw=np.outer(x.conj(),y)
        assert np.linalg.norm(raw-np.diag(np.diag(raw)))**2 >= .5-1e-13
        x=(.3+.1j)*x
        y=(.1+.4j)*x+.5*y
        product=np.outer(x.conj(),y)
        lower=(np.vdot(x,x).real*np.vdot(y,y).real-abs(np.vdot(x,y))**2)/2
        assert np.linalg.norm(product-np.diag(np.diag(product)))**2 >= lower-1e-13
    x=np.array([1,1,0])/sqrt(2)
    y=np.array([1,-1,0])/sqrt(2)
    pair=abs(np.outer(x,y))**2
    assert np.isclose(pair.sum()-np.trace(pair),.5)
    # HPQCD correction to Table IV for B+, vector/axial coefficients of Eq. (10).
    # LHCb 1909.01010 and Belle 2212.04128, individual PHSP and scalar/pseudoscalar
    # comparison templates, respectively. Their efficiency ratios remain explicit.
    channels=[("electron_muon",19.2,1.3,19.2,1.3,7.0e-9,6.4e-9),
              ("electron_tau",11.99,.75,11.99,.75,1.71e-5,2.08e-5),
              ("muon_tau",11.73,.75,12.18,.75,.65e-5,2.97e-5)]
    inv_alpha,reference_ckm=132.32,.04185
    norm=(2*pi*inv_alpha/reference_ckm)**2*1e-9
    coeffs={group:(a+b)*norm for group,a,ae,b,be,p,q in channels}
    def bound(mult,ratio):
        # Denominator fixed to table central calibration, optional rectangular
        # variation of it for conservative unnormalizing of the form-factor tables.
        lesser=(2*pi*(inv_alpha-mult*.05)/(reference_ckm+mult*.00093))**2*1e-9
        total=0.
        for group,a,ae,b,be,p,q in channels:
            bfloor=(a+b-mult*(ae+be))*lesser
            total+=(p+q)/(bfloor*ratio)
        return sqrt(2*total)
    return {"charged_neutral_flavour_checks":35,"orthogonal_saturation_squared_LFV":.5,
            "BplusKplus_branch_per_unit_ordered_LL_coefficient_squared":coeffs,
            "conditional_lambda_times_G0_over_Gmu_upper_if_efficiency_ratios_at_least_one":bound(0,1),
            "conditional_outer_two_width_lambda_if_efficiency_ratios_at_least_one":bound(2,1),
            "conditional_outer_two_width_lambda_if_efficiency_ratios_at_least_half":bound(2,.5)}


def compare(a,b):
    if isinstance(a,dict):
        return a.keys()==b.keys() and all(compare(a[k],b[k]) for k in a)
    return np.allclose(a,b,rtol=1e-10,atol=1e-14)


def crosscheck_charged_b():
    """Reconstruct two alternatives for the same tau row from primary inputs."""
    report=json.loads((HERE/"results.json").read_text(encoding="utf-8"))
    muon_mass=.10565838
    electron_mass=.000510999
    vector_floor=lambda mass: (1-mass**2/16)**2*(1+mass**2/32)
    high_e=(.380+2*(hypot(.024,.022)+hypot(.024,.074)))*1e-4
    high_mu=(.414+2*(hypot(.021,.019)+hypot(.023,.045)))*1e-4
    high_norm=(1.499-2*hypot(.013,.008))*(2.07-2*.57)*.98
    light=(high_e/vector_floor(electron_mass)+high_mu/vector_floor(muon_mass))/high_norm
    semi_tau=.000250/(high_norm*vector_floor(1.77686))
    assert np.isclose(sqrt(light+semi_tau),report["experimental_work_boxes"][2]["lambda_upper"])

    mass=5.27938-2*hypot(.00011,.00033)
    momentum=(mass**2-1.776**2)/(2*mass)
    lifetime=(1.639-2*hypot(.009,.009))/6.582119569e-13
    unit=(1.1663787e-5*(.1894-2*.0014)*1.776*momentum)**2*lifetime/(2*pi*mass)
    tau=(1.24+2*hypot(.41,.19))*1e-4/(unit*.98)
    assert np.isclose(unit,6.04983710194,rtol=2e-11)
    assert np.isclose(sqrt(light+tau),report["charged_B_tau_replacement"]["lambda_upper"])
    assert np.isclose(246.22/sqrt(2*sqrt(light+tau)),1650.53515081,rtol=2e-11)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument("--write",action="store_true")
    mode.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=run()
    path=HERE/"independent_results.json"
    if args.write:
        with path.open("x",encoding="utf-8") as f:
            json.dump(result,f,indent=2)
            f.write("\n")
    else:
        assert compare(result,json.loads(path.read_text(encoding="utf-8")))
        crosscheck_charged_b()
        print("1434 charged trace, restricted Gram determinant and ordered rare-B matching verified")
