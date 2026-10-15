"""1231: exact conservative 1+1 asymptotic scattering and unknown internal state.
Fraction identities/intervals carry all certificates. Finite matrices diagnose
only the logical reference channel. No local collision or Unruh implementation.
Default and --verify are read-only; --write-results creates a receipt once.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
TARGET=HERE/"results.json"
L=F(10**14)
SIGMA=F(1,100)
S=F(1,10**7)
S_BAD=F(1,10**16)
A_UP=F(3)
MU_MAX=F(4)
SOURCES=[
    "05_认知动力学公理与狭义相对论基础动力学_阶段论文.md",
    "archive_1106_1119/1118/proof.md",
    "archive_1106_1119/1119/proof.md",
    "archive_1219_/1219/candidate_proof.md",
    "archive_1219_/_admission/after1230_finite_motion_thermal_admission.md",
    "archive_1219_/1231/proof.md",
    "archive_1219_/1231/check.py",
]


def exp_interval(x,n=24):
    assert 0<=x<1
    terms=[x**j/math.factorial(j) for j in range(n+1)]
    low=sum(terms,F(0))
    first=x**(n+1)/math.factorial(n+1)
    high=low+first/(1-x/F(n+2))
    return low,high


def poly_integral(n,extra=0):
    return sum((F(2*(-1)**j*math.comb(n,j),2*j+extra+1)
                for j in range(n+1)),F(0))


def collision(mu,q,v):
    light_plus=mu*q+L*v
    light_minus=mu/q+L/v
    cmexp2=light_plus/light_minus
    return cmexp2/q,cmexp2/v


def ep(mu,q):
    return mu*(q+1/q)/2,mu*(q-1/q)/2


def exact_shell_and_measure():
    cases=[]
    for mu in (F(3),F(4)):
        for q,v in ((F(1),F(3,2)),(F(1,2),F(2)),(F(2),F(3,2))):
            qo,vo=collision(mu,q,v)
            ei,pi=ep(mu,q); es,ps=ep(L,v)
            eo,po=ep(mu,qo); eb,pb=ep(L,vo)
            assert ei+es==eo+eb and pi+ps==po+pb
            assert eo*eo-po*po==mu*mu and eb*eb-pb*pb==L*L
            assert collision(mu,qo,vo)==(q,v)
            u=q/v
            derivative=mu*u/(L+mu*u)+(mu/u)/(L+mu/u)
            det=(-1+derivative)*(1-derivative)-(2-derivative)*derivative
            assert det==-1
            boost=F(7,5)
            assert collision(mu,boost*q,boost*v)==(boost*qo,boost*vo)
            # Exact differences avoid cancellation against the huge source energy.
            assert eb-es==-(eo-ei) and pb-ps==-(po-pi)
            cases.append(dict(mu=int(mu),exp_eta=str(q),exp_xi=str(v),
                target_delta_E=float(eo-ei),source_delta_E=float(eb-es),
                target_delta_P=float(po-pi),source_delta_P=float(pb-ps),
                jacobian_exact=str(det),conservation_exact=True))
    return dict(status="PASS",mass_shell_measure="deta dxi = dp_A/E_A dp_B/E_B",
                map="(d,xi) -> (-d,xi+R_mu(d))",jacobian=-1,
                full_four_momentum_conserved=True,transverse_momenta_zero=True,
                involution=True,common_boost_covariance=True,cases=cases,
                exact_identities_not_a_grid_proof=True)


def joint_packet_bounds():
    i0=poly_integral(16)
    d1=256*poly_integral(14,2)/i0
    var=poly_integral(16,2)/i0
    assert d1==F(44,5)<A_UP*A_UP and var==F(1,35)
    dmax=F(1,2)+SIGMA+S
    assert dmax<F(52,100)
    _,exp_upper=exp_interval(F(52,100))
    assert exp_upper<2  # sinh D < e^D/2 < 1; cosh D < (2+1)/2 < 2.
    rmax=8/L
    rprime=16/L
    source_velocity=A_UP*2*S/(5*SIGMA)  # sqrt(35)>5.
    recoil_target=rmax*A_UP/SIGMA
    recoil_source=rmax*A_UP/S
    error=source_velocity+recoil_target+recoil_source
    assert error<F(15,10**6)
    energy=2*L+5
    assert F(1,2)+S<F(52,100)
    overlap_distance=2*A_UP/L*(1/SIGMA+1/S)
    coherence_loss=overlap_distance**2/2
    assert coherence_loss<F(181,10**15)
    return dict(status="PASS",target_masses=[3,4],source_mass=str(L),
                eta_center=0,xi_center=.5,target_output_center=1,
                target_rapidity_halfwidth=str(SIGMA),source_rapidity_halfwidth=str(S),
                envelope_order=8,regularity="C7 and H8, not C-infinity",
                I0_exact=str(i0),derivative_ratio_exact=str(d1),variance_ratio_exact=str(var),
                relative_rapidity_support_upper=float(dmax),R_upper=float(rmax),
                R_derivative_upper=float(rprime),
                three_error_terms=dict(source_velocity_spread=float(source_velocity),
                    target_packet_recoil=float(recoil_target),source_packet_recoil=float(recoil_source)),
                all_output_isometry_norm_error=float(error),exact_error_upper=str(error),
                half_diamond_error_upper=float(error),arbitrary_passive_reference=True,
                full_source_retained=True,source_preparation_independent_of_internal_input=True,
                initial_energy_support_upper_exact=str(energy),
                joint_motional_overlap_loss_upper=float(coherence_loss),
                exact_overlap_loss_upper=str(coherence_loss),
                reduced_source_not_assumed_pure=True,
                internal_coherence_coefficient="Gamma >= 1 - stated loss",
                arbitrary_unknown_motion_boost=False),error,rmax,rprime,coherence_loss


def narrow_source_and_reference(loss):
    xmin=F(489,1000)
    xmax=F(511,1000)
    sinh_upper=xmax/(1-xmax*xmax/6)
    assert sinh_upper<F(535,1000)
    _,exp_upper=exp_interval(xmax)
    assert exp_upper<2
    assert F(1,2)-SIGMA-S_BAD>xmin
    assert F(1,2)+SIGMA+S_BAD<xmax
    gap=(8*xmin/(1+8/L)-6*F(535,1000))/L-2*S_BAD
    assert gap>F(6,10**15)
    # Diagnostic channel on an entangled logical reference. The environment
    # coefficient is a certified Gram bound; no QFT/source wavefunction is replaced.
    g=1-float(loss)
    bell=np.array([1.,0.,0.,1.])/math.sqrt(2)
    rho=np.outer(bell,bell)
    def dephase(gamma):
        out=rho.astype(complex).copy()
        out[0,3]*=gamma
        out[3,0]*=gamma
        return out
    good=dephase(g)
    bad=dephase(0.)
    phase_flipped=dephase(-g)
    assert min(np.linalg.eigvalsh(good))>-1e-14
    assert abs(np.trace(good)-1)<1e-14
    assert np.allclose(np.diag(good),np.diag(bad))
    assert np.allclose(np.diag(good),np.diag(phase_flipped))
    bell_good=float(np.real(bell@good@bell))
    bell_bad=float(np.real(bell@bad@bell))
    assert abs(bell_bad-.5)<1e-14
    assert bell_good>1-1e-12
    return dict(status="PASS",bad_source_halfwidth=str(S_BAD),
                source_branch_support_gap_lower=float(gap),exact_support_gap=str(gap),
                bad_joint_environment_overlap_exact=0,
                internal_populations_unchanged=True,all_four_momentum_still_exact=True,
                full_global_information_not_erased=True,
                logical_reference_diagnostic_only=True,
                good_Bell_return_diagnostic=bell_good,bad_Bell_return_diagnostic=bell_bad,
                constant_relative_pi_phase_also_preserves_fourP=True,
                relative_scattering_phase_not_fixed_by_conservation=True)


def clock_position_phase(error,rmax,rprime):
    width=SIGMA+2*S+rmax
    assert width<F(10001,10**6)
    clock_error=2*width*width/(1-width*width)
    combined=error+clock_error
    assert combined<F(215,10**6)
    # Alternating cos(1) partial with last negative term is a strict lower bound.
    cos_lower=1-F(1,2)+F(1,24)-F(1,720)
    phase_record_gap=cos_lower-2*combined
    assert phase_record_gap>F(539,1000)
    xA=((1+rprime)*A_UP/SIGMA+rprime*A_UP/S+F(1,2))/3
    xB=((2+rprime)*A_UP/SIGMA+(1+rprime)*A_UP/S+F(1,2))/L
    assert xA<101 and xB<F(31,10**8)
    # An illustrative finite freely propagated tail; no finite collision is inferred.
    lab_wait=F(10)
    radius=F(10**6)
    position_tail=((xA+lab_wait*width)/radius)**2
    assert position_tail<F(11,10**9)
    return dict(status="PASS",output_rapidity_halfwidth_upper=float(width),
                central_comoving_proper_parameter=1,
                free_internal_phase_norm_error=float(clock_error),
                scattering_plus_free_phase_norm_error=float(combined),
                actual_contact_proper_time_not_constructed=True,
                comoving_generator="u_* . P_A = mu cosh(eta-1)",
                optional_X_and_rapidity_record=dict(
                    actual_readout_permission_not_constructed=True,
                    ideal_zero_phase_probability=(1+math.cos(1))/2,
                    ideal_pi_phase_probability=(1-math.cos(1))/2,
                    certified_probability_gap_if_same_readout_is_licensed=float(phase_record_gap),
                    rapidity_window_halfwidth=.02),
                target_asymptotic_X_norm_upper=float(xA),
                source_asymptotic_X_norm_upper=float(xB),
                physical_CO_position_identification_is_additional_contract=True,
                optional_free_position_tail=dict(lab_time=float(lab_wait),radius=float(radius),
                    probability_upper=float(position_tail),finite_collision_not_inferred=True),
                phase_zero_is_selected_model_input=True,
                arbitrary_phase_error_must_be_added=True,
                source_free_wait_retained=True)


def main():
    parser=argparse.ArgumentParser()
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument("--write-results",action="store_true")
    mode.add_argument("--verify",action="store_true")
    args=parser.parse_args()
    shell=exact_shell_and_measure()
    packets,error,rmax,rprime,loss=joint_packet_bounds()
    deletion=narrow_source_and_reference(loss)
    clock=clock_position_phase(error,rmax,rprime)
    result=dict(round=1231,status="PASS",check_groups=4,
                groups=dict(exact_shell_measure_and_recoil=shell,
                    one_normal_source_and_full_isometry=packets,
                    coherence_resource_and_phase_freedom=deletion,
                    free_clock_position_and_scope=clock),
                scope=dict(effective_asymptotic_scattering_boundary=True,
                    coherent_mass_input_permission_is_declared=True,
                    same_source_before_unknown_input=True,
                    source_and_target_identity_retained=True,
                    unknown_internal_state_and_arbitrary_reference=True,
                    arbitrary_unknown_motional_boost=False,
                    phase_zero_is_selected_scattering_law=True,
                    local_continuous_collision_dynamics=False,
                    finite_contact_window_or_worldtube=False,
                    actual_internal_X_readout_constructed=False,
                    original_288_dim_clock_not_replaced=True,
                    Dirac_field_coupling=False,Unruh_thermometer=False,
                    empirical_certification=False,stage_complete=False),
                source_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCES})
    if args.write_results:
        with TARGET.open("x",encoding="utf8",newline="\n") as out:
            json.dump(result,out,ensure_ascii=False,indent=2,allow_nan=False)
            out.write("\n")
    else:
        saved=json.loads(TARGET.read_text(encoding="utf8"))
        if saved!=result:
            raise SystemExit("Saved receipt differs; no files changed.")
    print(json.dumps(dict(round=1231,status="PASS",check_groups=4,
        all_output_norm_error=packets["all_output_isometry_norm_error"],
        coherence_loss_upper=packets["joint_motional_overlap_loss_upper"],
        narrow_source_support_gap=deletion["source_branch_support_gap_lower"],
        read_only=not args.write_results),ensure_ascii=False))


if __name__=="__main__":
    main()
