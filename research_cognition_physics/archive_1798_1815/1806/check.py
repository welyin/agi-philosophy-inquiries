"""A common Majorana scalar in the same finite heat-kernel trace."""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE/name).read_bytes())


def tr(a):
    return float(np.trace(a).real)


def couplings(A,G,S,omega,Bch):
    bnu=omega*tr(A@A)
    return ((Bch+bnu)/S, 4*omega*tr(A@G)/S, 2*omega*tr(G@G)/S)


def finite_matrix_check(S,omega,D,Bch):
    rng=np.random.default_rng(1806)
    Z=lambda: rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
    Y,R=Z(),Z()
    R=(R+R.T)/2
    Y*=math.sqrt(D/(omega*tr(Y.conj().T@Y)))
    R*=math.sqrt(S/(omega*tr(R.conj().T@R)))
    A,G=Y.conj().T@Y,R.conj().T@R
    lr,lp,ls=couplings(A,G,S,omega,Bch)
    def spectra(h,s):
        m=np.block([[np.zeros((3,3)), h*Y/math.sqrt(2)],
                    [h*Y.T/math.sqrt(2), s*R]])
        gram=m.conj().T@m
        return 2*omega*tr(gram),2*omega*tr(gram@gram)
    largest_error=0.
    for h,s in [(0.,1.),(1.,0.),(.3,-.9),(-.7,.8),(1.2,1.6)]:
        t2,t4=spectra(h,s)
        expected2=2*D*h*h+2*S*s*s
        expected4=(lr-Bch/S)*S*h**4+lp*S*h*h*s*s+ls*S*s**4
        largest_error=max(largest_error,abs(t2-expected2),abs(t4-expected4))
    assert largest_error<1e-12
    return dict(finite_complex_trace_test_max_absolute_error=largest_error,
                normalized_example=dict(lambda_H=lr,lambda_portal=lp,lambda_sigma=ls,
                    shift=lp**2/(4*ls),ratio_quadratic=lp/(2*ls)))


def calculate():
    cfg=read("inputs.json")
    old=read(cfg["source"])
    c,o=old["central"],old["outer"]
    S,om,D,B=c["S"],c["omega"],c["DeltaA"],c["B_ch"]
    data=c["lambda_data"]
    xmax=(1+math.sqrt(3))/2
    # Sharp suprema over positive/semidefinite spectra with at most three singlet lines.
    rmax=xmax*D/S
    max_shift=2*D**2/(om*S)
    low=(B-D**2/om)/S
    eq_shift=max_shift/3
    eqrange=[c["critical_lambda_rank"][str(k)]-eq_shift for k in (3,1)]
    top_norm=D/om
    # Diagonal witness saturates the maximal quadratic ratio; scale fixes the sigma kinetic term.
    traceRshape=np.array([1.,1/(1+math.sqrt(3)),1/(1+math.sqrt(3))])
    A=np.diag([top_norm,0.,0.])
    G=np.diag(traceRshape*S/(om*traceRshape.sum()))
    _,lp,ls=couplings(A,G,S,om,B)
    assert abs(lp/(2*ls)-rmax)<1e-14
    physical_mass_floor=2/math.sqrt(xmax)
    assert abs(2*math.sqrt(tr(G@G)/(tr(G)*max(np.diag(G))))-physical_mass_floor)<1e-14
    # Equal M: vary the physical A spectrum. Both extrema give the same shift.
    eqG=np.eye(3)*S/(3*om)
    for mat,target in [(np.eye(3)*top_norm/3,eqrange[0]),(A,eqrange[1])]:
        lh,lp,ls=couplings(mat,eqG,S,om,B)
        assert abs(lh-lp**2/(4*ls)-target)<1e-14
    gauge=read(cfg["gauge"])
    om_min=gauge["one_loop_corner_bounds"]["wl_wq"][0]-1e-9
    smin,smax=o["S"]
    dmax=smax-o["A_ch"][0]
    fraction=1-o["A_ch"][0]/smax
    bound_r=fraction*xmax
    weak_input=cfg["lower_matching_weak_mass_coefficient_GeV2"]
    hz=gauge["one_loop_corner_bounds"]["mu_GeV"]
    log_span=math.log(hz[1]/cfg["weak_mass_bound_reference_scale_GeV"])
    trace_upper=cfg["singlet_RGE_test_T_upper"]
    quartic_upper=cfg["SM_lambda_upper"]
    scale_factor=math.exp((12*quartic_upper+2*trace_upper)*log_span/(16*math.pi**2))
    weak_upper=weak_input*scale_factor
    max_sigma=math.sqrt(2*weak_upper/(1-bound_r))
    # The inequalities below apply only with the neutral trace weights commensurate.
    n_min=lambda f: math.ceil((2/f-1)**2)
    mass_boxes={}
    for k in cfg["mass_window_ratio_examples"]:
        factor=max(0.,min(1.,2/(1+2/k**4)-1/3))
        leff=o["B_ch_lower"]/smax-factor*dmax**2/(om_min*smin)
        mass_boxes[str(k)]=dict(lambda_eff_lower=leff,
                               extra_lambda_eff_minus_data_lower=leff-o["lambda_data"][1])
    return dict(central=dict(any_majorana_spectrum_lower_lambda_eff=low,matched_low_lambda=data,
                        necessary_additional_displacement=low-data,max_tree_downward_shift=max_shift,
                        equal_three_masses_shift=eq_shift,equal_three_masses_possible_lambda_eff=eqrange,
                        max_ratio_shared_quadratic=rmax,
                        minimal_fraction_extra_quadratic_while_H_light=1-rmax,
                        minimal_equal_trace_right_lines_to_reach_full_quadratic_cancellation=n_min(D/S),
                        m_sigma_over_M_max_interval=[physical_mass_floor,2.]),
                complex_check=finite_matrix_check(S,om,D,B),
                full_box=dict(max_missing_fraction=fraction,max_quadratic_ratio=bound_r,
                    min_relative_quadratic_gap=1-bound_r,
                    integrated_SM_mass_growth_upper=scale_factor,weak_quadratic_upper_GeV2=weak_upper,
                    max_symmetric_origin_sigma_mass_GeV=max_sigma,
                    common_heavy_scale_min_GeV=hz[0],
                    max_relative_allowed_residue_for_sigma_at_least_common_scale=2*weak_upper/hz[0]**2,
                    min_same_weight_R_lines_for_exact_quadratic_cancellation=n_min(fraction),
                    three_nearby_masses_by_max_over_min=mass_boxes))


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    val=calculate()
    raw=json.dumps(val,indent=2,allow_nan=False)+"\n"
    p=HERE/"results.json"
    if args.check:
        assert json.loads(p.read_bytes())==val
    else:
        p.write_bytes(raw.encode())
    print("1806",hashlib.sha256(raw.encode()).hexdigest()[:12])

