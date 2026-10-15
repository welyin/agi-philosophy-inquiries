"""1228: one Dirac vacuum, two full-record menus, same-source/clock budget.
Default is read-only. --write-results creates a new receipt, never overwrites.
Exact Fraction guards certify analytic margins; finite CAR NumPy is a diagnostic.
No full QFT spectrum, trajectory, thermal response, or experimental calibration.
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
K=10**6
B=F(21,4)
C=16*K
A0=2**28
T=F(2)
DELTA=F(1,10**24)
SIGMA=F(1,10**120)
M=F(10**260)
SOURCES=[
    "archive_1120_1158/1126/proof.md",
    "archive_1120_1158/1128/proof.md",
    "archive_1201_1218/1203/proof.md",
    "archive_1219_/1227/proof.md",
    "archive_1219_/1227/check.py",
    "archive_1219_/1227/results.json",
    "archive_1219_/_admission/after1227_reference_flow_priority.md",
    "archive_1219_/1228/proof.md",
    "archive_1219_/1228/check.py",
]


def alternating_arctan(q, n):
    partial=sum(((-1)**j*F(q)**(2*j+1)/F(2*j+1)
                 for j in range(n)), F(0))
    next_term=(-1)**n*F(q)**(2*n+1)/F(2*n+1)
    return min(partial,partial+next_term),max(partial,partial+next_term)


def pi_interval():
    a,b=alternating_arctan(F(1,5),60)
    c,d=alternating_arctan(F(1,239),20)
    return 16*a-4*d,16*b-4*c


def r_interval(pilo,pihi):
    # r=sum (-1)^j*pi^(2j)/(j!*(35/2)_j); decreasing terms from j=0.
    lo=hi=F(0)
    den=F(1)
    n=30
    for j in range(n):
        if j:
            den*=j*F(35+2*(j-1),2)
        low=pilo**(2*j)/den
        high=pihi**(2*j)/den
        if j%2:
            lo-=high
            hi-=low
        else:
            lo+=low
            hi+=high
    den*=n*F(35+2*(n-1),2)
    # n even: the omitted alternating tail is positive and <= its first term.
    hi+=pihi**(2*n)/den
    assert pihi*pihi/F(35,2)<1
    return lo,hi


def analytic_flow():
    pilo,pihi=pi_interval()
    assert F(314,100)<pilo<pihi<F(22,7)
    assert pihi-pilo<F(1,10**80)
    d1=F(8*(4*8+1),2*(2*8-1))
    assert d1==F(44,5) and 1+3*d1==F(137,5)<B*B
    rlo,rhi=r_interval(pilo,pihi)
    assert F(43,100)<rlo<rhi<F(59,100)
    # Independent broad moment route, not needed for the tighter Fraction sum.
    kappa=2*pihi
    rlow_moment=1-kappa*kappa/70
    rhigh_moment=1-(2*pilo)**2/70+(2*pihi)**4/F(8*35*37)
    assert rlow_moment>F(43,100) and rhigh_moment<F(59,100)
    eps_t=2*kappa*B/K
    eps_W=kappa*(3*B+F(1,2)+3*kappa)/K
    covariance=B/K
    gap_broad=(1-F(59,100)**2)/4-covariance-2*(eps_W+eps_t)
    gap_tight=(1-rhi*rhi)/4-covariance-2*(eps_W+eps_t)
    assert gap_broad>F(1623,10000)
    assert 4*kappa/K<F(1,10000)
    return dict(status="PASS",K=K,mass=1,
                envelope_order=8,envelope_regularity="C7 and H8, not C-infinity",
                d1_rational=str(d1),B_squared="137/5",B_upper=str(B),
                pi_lower_display=float(pilo),pi_upper_display=float(pihi),
                pi_width_upper="1e-80",
                r_lower_display=float(rlo),r_upper_display=float(rhi),
                r_certified_broad_interval=[.43,.59],
                eta_upper=float(kappa/K),tau_upper=float(2*kappa/K),
                physical_center_separation_upper=float(4*kappa/K),
                high_carrier_time_error_upper=float(eps_t),
                high_carrier_boost_error_upper=float(eps_W),
                actual_vacuum_covariance_error_upper=float(covariance),
                ideal_event_gap_lower_broad=float(gap_broad),
                ideal_event_gap_lower_tighter=float(gap_tight),
                certificate="exact rational bounds, no floating quadrature"),gap_broad


def annihilator(j,n):
    a=np.zeros((2**n,2**n),complex)
    for s in range(2**n):
        if (s>>j)&1:
            lower=(s&((1<<j)-1)).bit_count()
            a[s^(1<<j),s]=(-1)**lower
    return a


def full_car_pointers(r):
    # Four-mode realization of only the finite comparator algebra.
    # Modes 0,1 empty, 2,3 filled. This is NOT the continuum physical vacuum.
    aa=[annihilator(j,4) for j in range(4)]
    I=np.eye(16)
    f=(aa[0]+aa[2])/math.sqrt(2)
    g=(r*aa[0]+math.sqrt(1-r*r)*aa[1]+
       r*aa[2]+math.sqrt(1-r*r)*aa[3])/math.sqrt(2)
    nf=f.conj().T@f
    ng=g.conj().T@g
    omega=np.zeros(16,complex)
    omega[12]=1
    assert np.linalg.norm(nf@nf-nf)<1e-14
    assert np.linalg.norm(ng@ng-ng)<1e-14
    # Full PE1 PE2 space: each PE retains 4 states, so 16*4*4=256 amplitudes.
    out=np.zeros((16,4,4),complex)
    probs={}
    for x in (0,1):
        for y in (0,1):
            branch=(ng if y else I-ng)@(nf if x else I-nf)@omega
            out[:,3*x,3*y]=branch
            probs[f"{x}{y}"]=float(np.vdot(branch,branch).real)
    p=float(np.vdot(nf@omega,ng@nf@omega).real)
    covariance=np.vdot(omega,f.conj().T@g@omega)
    formula=.25-abs(covariance)**2+.5*r*r
    residual=max(abs(p-formula),abs(p-(1+r*r)/4),
                 abs(np.vdot(out,out).real-1))
    assert residual<1e-13
    # Exact equal-mode endpoint must reduce to p11=1/2.
    same=float(np.vdot(nf@omega,nf@nf@omega).real)
    assert abs(same-.5)<1e-14
    # PE pointers are not erased; tracing E produces classical P records.
    return dict(status="PASS",comparison_only=True,
                full_comparator_field_dimension=16,full_PE_dimension=16,
                full_joint_dimension=256,comparison_r=r,
                comparison_probabilities=probs,
                same_mode_p11=same,max_float_residual=residual,
                finite_comparator_not_actual_QFT_vacuum=True,
                full_probability_formula="Cff*Cgg-|Cfg|^2+|<f,g>|^2*Cff")


def source_and_clocks(gap):
    assert A0>8*2**8*16**4
    assert (F(15,16)**8)/2>F(1,4)  # integral b over [-1/4,1/4].
    assert 16**2*3<28**2
    assert 2*F(22,7)*(2*K+14)<C
    assert DELTA<F(1,10**10)
    L=16/(DELTA*DELTA)
    bT=SIGMA+T/(2*M*SIGMA)
    epsC=L*T*bT
    epsp=2*C*DELTA
    RH=F(2*C)
    RQ=F(A0+2*C)
    RHA=2*T*C/DELTA
    RQA=A0+RHA
    errQ=(RQA+RQ)*epsC+2*RQ*epsp
    errH=(RHA+RH)*epsC+2*RH*epsp
    gap_actual=gap-4*(epsC+epsp)
    assert errQ<F(1,10**6) and errH<F(1,10**6)
    assert gap_actual>F(16,100)
    group=dict(status="PASS",stress_component="vacuum-normal-ordered T00",
               stress_Q="rho(t)=b8(t-4)/integral(b8); q=b8(x)b8(y)b8(z-2)",
               stress_role="future free-field source functional of the same terminal state",
               stress_not_actual_driven_future_history=True,
               stress_vacuum_norm_upper=A0,
               per_gate_H_and_Q_commutator_upper=C,
               classical_and_ideal_H_norm_upper=float(RH),
               classical_and_ideal_Q_norm_upper=float(RQ),
               actual_H_norm_upper=float(RHA),
               actual_Q_norm_upper=float(RQA),
               clock_to_classical_vector_error_upper=float(epsC),
               classical_to_ideal_vector_error_upper=float(epsp),
               per_branch_Q_absolute_error_upper=float(errQ),
               full_field_energy_absolute_error_upper=float(errH),
               actual_two_menu_event_gap_lower=float(gap_actual),
               operator_Schatten_scope="HS controls T(Q) Omega only; not bounded Fock operator",
               source_values_not_numerically_integrated=True)
    endpoint=L*(3*SIGMA+L*T*T/(2*M))
    n=(SIGMA*L*T+L*T*T/(4*M*SIGMA))**2
    EC=2*(M/2+1/(8*M*SIGMA*SIGMA)+1)
    assert endpoint<F(1,10**60)
    assert n<F(1,10**140)
    # For n<1e-140, monotone n(2-ln n), ln10<3 proves sum S<1e-136.
    Sbound=2*F(1,10**140)*(2+140*3)
    assert Sbound<F(1,10**136)
    resources=dict(status="PASS",clock_count=2,delta="1e-24",
                   gaussian_sigma="1e-120",each_clock_M="1e260",T=2,
                   force_upper=float(L),initial_two_clock_energy_display=float(EC),
                   initial_pointer_energy=4,initial_field_energy=0,
                   initial_and_terminal_interaction_sum_upper=float(endpoint),
                   field_energy_gain_upper=float(RH+errH),
                   energy_reaction="|Delta E_C1+C2 + Delta E_field| <= endpoint bound",
                   each_clock_ref_free_occupation_upper=float(n),
                   both_clock_entropy_certified_upper="1e-136",
                   apparatus_baseline_energy=6,
                   finite_not_practically_feasible=True,
                   clock_has_no_physical_boost_or_accelerated_proper_time_identity=True,
                   source_model_permission="finite-region CAR effective kernel, not pointwise causal implementation")
    return group,resources


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    flow,gap=analytic_flow()
    r=(flow["r_lower_display"]+flow["r_upper_display"])/2
    pointer=full_car_pointers(r)
    source,resources=source_and_clocks(gap)
    result=dict(round=1228,status="PASS",check_groups=4,
                groups=dict(analytic_reference_flow=flow,
                            full_CAR_and_entity_record=pointer,
                            same_terminal_source_transport=source,
                            actual_clock_energy_entropy=resources),
                scope=dict(massive_free_Dirac=True,flavors=2,
                           both_energy_signs_retained=True,
                           continuum_not_truncated=True,
                           accelerated_thermometer=False,
                           microscopic_local_implementation=False,
                           all_apparatus_boost_conservation=False,
                           finite_region_effective_gate_permission_is_input=True,
                           empirical_certification=False,
                           stage_complete=False),
                source_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
                               for p in SOURCES})
    if args.write_results:
        with TARGET.open("x",encoding="utf8",newline="\n") as out:
            json.dump(result,out,ensure_ascii=False,indent=2,allow_nan=False)
            out.write("\n")
    else:
        saved=json.loads(TARGET.read_text("utf8"))
        if result!=saved:
            raise SystemExit("Saved receipt differs; no files changed.")
    print(json.dumps(dict(status="PASS",round=1228,check_groups=4,
                          event_gap_lower=source["actual_two_menu_event_gap_lower"],
                          branch_stress_error_upper=source["per_branch_Q_absolute_error_upper"],
                          read_only=not args.write_results),ensure_ascii=False))


if __name__=="__main__":
    main()

