"""1229: finite same-field information, actual wedge/auxiliary modular domain.
Four groups. Exact Fraction guards certify bounds; finite CAR matrices only
check universal instrument identities, not a truncation of the physical QFT.
Default read-only; --write-results creates a new receipt exclusively.
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
B2=F(137,5)
BUP=F(21,4)
C=16*K
CB=128*K
DELTA=F(1,10**24)
SIGMA=F(1,10**120)
M=F(10**260)
L=16/DELTA**2
PIUP=F(22,7)
SOURCES=[
    "archive_1120_1158/1126/proof.md",
    "archive_1120_1158/1128/proof.md",
    "archive_1201_1218/1203/proof.md",
    "archive_1219_/1228/proof.md",
    "archive_1219_/1228/check.py",
    "archive_1219_/1228/results.json",
    "archive_1219_/1228/NEXT.md",
    "archive_1219_/_admission/after1228_geometric_information_admission.md",
    "archive_1219_/1229/proof.md",
    "archive_1219_/1229/check.py",
]


def d(T):
    return SIGMA+F(T)/(2*M*SIGMA)


def epsA(T):
    return L*T*d(T)


def nbar(T):
    return (SIGMA*L*T+L*T*T/(4*M*SIGMA))**2


def region_and_aux_domain():
    # Both physical kernels for all t in [0,4].
    a=8
    phys_zmin=-1  # second physical kernel is e^(-2ih)u.
    max_time_distance=3
    margin=phys_zmin+a-max_time_distance
    assert margin==4
    cauchy_zmin=1-3
    assert cauchy_zmin+a==6
    assert 14*(K+BUP)+F(1,2)<15*K
    assert 2*PIUP*15*K<CB
    aa=32*d(4)/DELTA**2
    bb=64*d(4)**2/DELTA**3
    aux_m2=(1+2*(aa+bb))**4
    assert aux_m2<4
    assert aa+bb<F(1,10**69)
    assert nbar(4)<F(1,10**140)
    # ln 10 < 3, g(n)<=n*(2-ln n) and n<1e-140.
    clock_entropy=2*F(1,10**140)*(2+140*3)
    assert clock_entropy<F(1,10**136)
    return dict(status="PASS",wedge="z+8 > |t-1|",
                total_time_window=[0,4],preparation_cut=2,gate_centers=[1,3],
                physical_spacetime_support_margin=margin,
                reference_slice_Cauchy_support_margin=6,
                first_slice_one_particle_generator="k_z+8*h",
                initial_global_slice_generator="k_z+8*h+p_z",
                final_Schrodinger_generator="k_z+8*h-(T-1)*p_z",
                field_boost_commutator_upper=CB,
                uniform_aux_comm_a=float(aa),uniform_aux_comm_b=float(bb),
                uniform_aux_Nplus1_second_moment_certified_upper=4,
                uniform_aux_N_tail_first_moment="<=4/R",
                each_clock_occupation_upper=float(nbar(4)),
                both_clock_entropy_upper="1e-136",
                finite_pointer_entropy_upper="4*ln(2)",
                actual_aux_entropy_upper=2.8,
                faithful_aux_reference="exp(-N)/Z tensor I_PE/16; Z=(1-exp(-1))^-2",
                initial_aux_pure_not_equal_reference=True,
                relative_entropy_proof="finite auxiliary identity; fixed-reference lsc; then chain rule",
                relative_entropy_equality_not_claimed=True)


def ann(j):
    a=np.zeros((4,4),complex)
    for s in range(4):
        if (s>>j)&1:
            a[s^(1<<j),s]=(-1)**((s&((1<<j)-1)).bit_count())
    return a


def car_channel():
    a,b=ann(0),ann(1)
    I=np.eye(4)
    na=a.conj().T@a
    nb=b.conj().T@b
    ff=(a+b)/math.sqrt(2)
    nf=ff.conj().T@ff
    lhs=nf@na@nf+(I-nf)@na@(I-nf)
    identity_error=float(np.linalg.norm(lhs-(na+nb)/2))
    assert identity_error<1e-14
    # Full two-PE isometry, keeping both 4-state PE spaces.
    V=np.zeros((4,4,4,4),complex)
    for x in (0,1):
        for y in (0,1):
            V[:,3*x,3*y,:]=(na if y else I-na)@(nf if x else I-nf)
    VV=V.reshape(64,4)
    iso=float(np.linalg.norm(VV.conj().T@VV-I))
    assert iso<1e-14
    rows=[]
    for s in (F(1,10),F(1,1000),F(1,2)):
        sf=float(s)
        rho=np.diag([(1-sf)*sf,sf*sf,(1-sf)**2,sf*(1-sf)])
        prepared=nf@rho@nf+(I-nf)@rho@(I-nf)
        p=float(np.trace(prepared@na).real)
        out=VV@rho@VV.conj().T
        q2=sum(float(out[(z*4+x)*4+3,(z*4+x)*4+3].real)
               for z in range(4) for x in range(4))
        assert abs(p-.5)<1e-14 and abs(q2-.5)<1e-14
        assert abs(np.trace(rho@na).real-sf)<1e-14
        rows.append(dict(s_rational=str(s),ideal_after_first=p,
                         final_second_pointer=q2))
    return dict(status="PASS",full_finite_identity_domain_dimension=4,
                full_PE_dimension=16,full_isometry_output_dimension=64,
                identity_error=identity_error,isometry_error=iso,rows=rows,
                same_field_effect="n_(f_plus in flavor1)",
                physical_s_not_numerically_evaluated=True,
                actual_continued_readout_not_exact_field_only_channel=True)


def information():
    smin=1/(4*(K*K+B2))
    smax=B2/(K*K)
    qstar=F(1,5*10**12)
    alpha=2*(epsA(2)+C*DELTA)
    gamma=2*(epsA(4)+2*C*DELTA)
    eta=alpha+gamma
    zeta=2*(epsA(4)+C*DELTA)
    assert smin>qstar and smax<1-qstar
    assert eta<F(1,10**15) and zeta<F(1,10**16)
    assert F(1,2)-alpha>F(49,100)
    assert F(1,2)-gamma>F(49,100)
    # e < 2.72, proved by positive series and geometric tail.
    e10=sum((F(1,math.factorial(j)) for j in range(11)),F(0))
    e_upper=e10+F(1,math.factorial(11))*F(12,11)
    assert e_upper<F(68,25)
    assert F(68,25)**24<3*10**10
    assert smax+zeta<F(1,3*10**10)  # < exp(-24).
    # exp(.7)>2, hence ln 2<.7.
    exp07lower=sum((F(7,10)**j/math.factorial(j) for j in range(6)),F(0))
    assert exp07lower>2
    KLlower=F(49,100)*24-F(7,10)
    assert KLlower>11
    # h2(eta)<=eta*(1-ln eta), ln10<3; ln(1/qstar)<38.
    entropy_error=F(1,10**15)*(1+15*3)
    cross_error=2*F(1,10**15)*38
    ref_error=F(1,10**16)/(qstar-F(1,10**16))
    transport=entropy_error+cross_error+ref_error
    assert transport<F(1,1000)
    return dict(status="PASS",vacuum_s_lower=float(smin),vacuum_s_upper=float(smax),
                vacuum_s_lower_formula="1/[4(K^2+137/5)]",
                vacuum_s_upper_formula="(137/5)/K^2",
                common_probability_floor=float(qstar),
                actual_prep_field_probability_error=float(alpha),
                actual_full_second_pointer_error=float(gamma),
                actual_pointer_to_same_field_channel_TV_upper=float(eta),
                vacuum_actual_readout_error=float(zeta),
                finite_record_KL_transport_upper=float(transport),
                fixed_field_KL_lower=float(KLlower),
                actual_record_KL_lower=float(KLlower),
                exact_KL_values_not_evaluated=True,
                actual_record_DPI_with_correction_only=True,
                correction_claim="Drecord - 1e-3 + I(W:E) <= DeltaK + S(E)")


def source_geometry_resources():
    eC=epsA(2)
    ep=C*DELTA
    RH=F(2*C)
    RB=F(2*CB)
    RHA=F(4*C)/DELTA
    RBA=F(4*CB)/DELTA
    errH=(RHA+RH)*eC+2*RH*ep
    errB=(RBA+RB)*eC+2*RB*ep
    errK=2*PIUP*errB
    assert max(errH,errB,errK)<F(1,10**6)
    Ebound=RH+errH
    Bbound=RB+errB
    Kbound=2*PIUP*Bbound
    geo=F(1,10**12)
    geoerr=geo*Kbound+(1+geo)*2*PIUP*(2*geo)*Ebound
    assert geoerr<F(1,100)
    endpoint=L*(3*SIGMA+L*16/(2*M))
    assert endpoint<F(1,10**60)
    free_initial=M+1/(4*M*SIGMA*SIGMA)+6
    # Positive information and finite actual auxiliary entropy imply true boost >8.2.
    assert F(11)-F(28,10)==F(82,10)
    return dict(status="PASS",field_energy_transport_error=float(errH),
                field_boost_expectation_transport_error=float(errB),
                modular_energy_transport_error=float(errK),
                actual_field_energy_upper=float(Ebound),
                actual_field_boost_absolute_upper=float(Bbound),
                actual_modular_energy_absolute_upper=float(Kbound),
                independent_uniform_geometry_error={
                    "longitudinal_origin_absolute":1e-12,
                    "time_origin_absolute":1e-12,
                    "area_common_relative_scale":1e-12},
                geometry_transport_error_upper=float(geoerr),
                geometry_calibration_is_conditional_input_not_new_instrument=True,
                source_role="same-slice renormalized expectation distribution / generator form",
                no_arbitrary_sharp_time_stress_operator_claim=True,
                free_initial_total_energy_display=float(free_initial),
                exact_total_initial_energy="free initial energy plus <V(0)>",
                initial_plus_terminal_interaction_bound=float(endpoint),
                actual_modular_energy_strict_lower=8.2,
                clock_energy_boost_identity_not_claimed=True,
                no_area_law_or_gravity_claim=True)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--write-results",action="store_true")
    args=ap.parse_args()
    groups=dict(actual_wedge_and_auxiliary_domain=region_and_aux_domain(),
                full_field_channel_identity=car_channel(),
                finite_same_field_information=information(),
                source_geometry_and_resource_transport=source_geometry_resources())
    result=dict(round=1229,status="PASS",check_groups=4,groups=groups,
                source_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
                               for p in SOURCES},
                scope=dict(actual_untruncated_field_and_clocks=True,
                           full_auxiliary_kept=True,
                           strict_finite_window_wedge_localization=True,
                           finite_relative_entropy_proven_before_chain_rule=True,
                           actual_entropy_energy_inequality=True,
                           entropy_energy_equality=False,
                           actual_readout_exact_field_channel=False,
                           finite_readout_error_included=True,
                           independent_geometry_input=True,
                           local_microscopic_gate_realization=False,
                           accelerated_thermometer=False,area_law=False,
                           stage_complete=False))
    if args.write_results:
        with TARGET.open("x",encoding="utf8",newline="\n") as out:
            json.dump(result,out,ensure_ascii=False,indent=2,allow_nan=False)
            out.write("\n")
    elif result!=json.loads(TARGET.read_text("utf8")):
        raise SystemExit("Saved receipt differs; no files changed.")
    print(json.dumps(dict(round=1229,status="PASS",check_groups=4,
                          record_KL_lower=groups["finite_same_field_information"]["actual_record_KL_lower"],
                          KL_transport_upper=groups["finite_same_field_information"]["finite_record_KL_transport_upper"],
                          read_only=not args.write_results),ensure_ascii=False))


if __name__=="__main__":
    main()

