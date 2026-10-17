"""1287 exact finite causal-record/source certificate.
One author group. Standard-library Fraction only; no field discretization,
propagator scan, Bessel library or quadrature. Default is read-only replay;
--write exclusively creates results.json. Continuous claims are in proof.md.
"""
import argparse
import hashlib
import json
import math
from fractions import Fraction as F
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in HERE.parents if (p/'research_cognition_physics').is_dir() and (p/'scripts').is_dir())
RESEARCH = ROOT/'research_cognition_physics'
METRIC = (-1, 1, 1, 1)
MASS = [[F(5,2), F(-3,2)], [F(-3,2), F(5,2)]]


def fs(x):
    x = F(x)
    return f'{x.numerator}/{x.denominator}'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def mm(a, b):
    return [[sum(a[i][k]*b[k][j] for k in range(len(b))) for j in range(len(b[0]))] for i in range(len(a))]


def mv(a, v):
    return [sum(x*y for x,y in zip(row,v)) for row in a]


def ms(a):
    return [[fs(x) for x in row] for row in a]


def exp_upper(x, n):
    # Positive Taylor series; after term n the successive ratio is <=x/(n+2).
    partial = sum((x**j/F(math.factorial(j)) for j in range(n+1)), F(0))
    tail = x**(n+1)/F(math.factorial(n+1))/(1-x/F(n+2))
    return partial+tail


def bessel_bracket(x, order):
    # [2 J1(2 sqrt(x))-J1(sqrt(x))]/sqrt(x); odd order lower, even upper.
    return sum(((-1)**n*(F(2)-F(1,2*4**n))*x**n/F(math.factorial(n)*math.factorial(n+1))
                for n in range(order+1)), F(0))


def stress(phi, first, mass):
    q = sum(METRIC[a]*sum(first[i][a]**2 for i in range(2)) for a in range(4))
    q += sum(x*y for x,y in zip(phi,mv(mass,phi)))
    return [[sum(first[i][a]*first[i][b] for i in range(2))
             -(METRIC[a]*q/2 if a == b else 0) for b in range(4)] for a in range(4)]


def divergence(phi, first, second, mass_in_stress):
    # Differentiate the full Hilbert tensor componentwise, not a precomputed Ward residual.
    mphi = mv(mass_in_stress,phi)
    answer=[]
    for b in range(4):
        kinetic = sum(METRIC[a]*sum(second[i][a][a]*first[i][b]
                      +first[i][a]*second[i][a][b] for i in range(2)) for a in range(4))
        dq = 2*sum(METRIC[a]*sum(first[i][a]*second[i][b][a] for i in range(2)) for a in range(4))
        dq += 2*sum(mphi[i]*first[i][b] for i in range(2))
        answer.append(kinetic-dq/2)
    return answer


def source_checks():
    max_residual=F(0)
    deleted=[]
    diag=[[MASS[0][0],F(0)],[F(0),MASS[1][1]]]
    for j in range(1,25):
        phi=[F(2*j-11,7),F(3*j-13,11)]
        first=[[F(((j+2)*(i+1)*(a+2))%19-9,5) for a in range(4)] for i in range(2)]
        second=[]
        for i in range(2):
            h=[[F((j*(i+2)*(min(a,b)+1)+max(a,b))%17-8,9) for b in range(4)] for a in range(4)]
            second.append(h)
        mphi=mv(MASS,phi)
        for i in range(2):
            second[i][0][0]=sum(second[i][a][a] for a in (1,2,3))-mphi[i]
        div=divergence(phi,first,second,MASS)
        assert div == [0]*4
        max_residual=max(max_residual,*(abs(x) for x in div))
        t=stress(phi,first,MASS)
        assert t == stress([-x for x in phi],[[-x for x in row] for row in first],MASS)
        assert t[0][0]>=0 and all(abs(x)<=t[0][0] for row in t for x in row)
        deleted_div=divergence(phi,first,second,diag)
        if any(deleted_div) and not deleted:
            deleted=[fs(x) for x in deleted_div]
    assert deleted
    # Flavor singlet stress commutes with both Pauli generators, including offdiagonal one.
    for scalar in (F(2,3),F(-7,5)):
        tf=[[scalar,F(0)],[F(0),scalar]]
        for gen in ([[0,1],[1,0]],[[1,0],[0,-1]]):
            assert mm(tf,gen)==mm(gen,tf)
    return {'finite_jet_diagnostics':24, 'complete_divergence_max':fs(max_residual),
            'omitted_mixing_divergence_witness':deleted,
            'opposite_coherent_sign_same_full_classical_tensor':True,
            'component_absolute_bound_by_classical_energy_checked':True,
            'flavor_source_offdiagonal_zero':True,
            'continuous_identity_proved_analytically_not_by_jet_sampling':True}


def calculate():
    identity=[[F(1),F(0)],[F(0),F(1)]]
    pp=[[F(1,2),F(1,2)],[F(1,2),F(1,2)]]
    pm=[[F(1,2),F(-1,2)],[F(-1,2),F(1,2)]]
    assert mm(pp,pp)==pp and mm(pm,pm)==pm
    assert mm(pp,pm)==[[0,0],[0,0]]
    assert [[pp[i][j]+pm[i][j] for j in range(2)] for i in range(2)]==identity
    assert mm(MASS,pp)==pp
    assert mm(MASS,pm)==[[4*x for x in row] for row in pm]
    comp=[[MASS[i][j]-identity[i][j] for j in range(2)] for i in range(2)]
    assert mv(comp,[F(1),F(1)])==[0,0]
    assert comp[0][0]*comp[1][1]-comp[0][1]**2==0
    # PSD Gram factor: B=(3/2) vv^T, v=(1,-1).
    assert comp==[[F(3,2),F(-3,2)],[F(-3,2),F(3,2)]]

    a,L,t0,t1,zeta=F(1,32),F(1,4),F(7,16),F(1,2),F(1,256)
    rmin,rmax=L-2*a,L+2*a
    s2min,s2max=t0*t0-rmax*rmax,t1*t1-rmin*rmin
    assert 0<rmin<rmax<t0 and (s2min,s2max)==(F(3,32),F(55,256))
    xupper=F(1,4)
    assert s2max<xupper
    # Infinite-series ratio bounds: x/((n+1)(n+2)) <= 1/8;
    # for the combined difference use numerator ratio <=4/3, hence <=1/6.
    ratio_mass2=xupper/2
    ratio_mass1=xupper/8
    ratio_combined=F(4,3)*xupper/2
    assert ratio_mass2==F(1,8) and ratio_mass1==F(1,32) and ratio_combined==F(1,6)
    kernel_num_lower=2-xupper-F(1,2)
    assert kernel_num_lower==F(5,4)
    assert F(1,2)-xupper/16>0
    samples=[]
    for x in (s2min,(s2min+s2max)/2,s2max):
        low,high=bessel_bracket(x,9),bessel_bracket(x,10)
        assert F(5,4)<low<high<2
        samples.append({'s_squared':fs(x),'numerator_lower':fs(low),'numerator_upper':fs(high)})

    # b<=e^-1<1/2; on inner half ball b>e^-4/3>1/4.
    e_upper=exp_upper(F(1),8)
    assert F(2)<e_upper<F(11,4)
    assert F(11,4)**4<4**3
    l2upper=4/a**3
    variance=F(3,8)*l2upper
    exponent=zeta*zeta*variance/2
    assert (l2upper,variance,exponent)==(F(131072),F(49152),F(3,8))
    assert 1/(1-exponent)<2
    # pi>3 gives zeta*m < zeta/12<1. sin(x)>=x-x^3/6>x/2.
    phase_upper=zeta/12
    assert phase_upper<1 and 1-phase_upper**2/6>F(1,2)
    # m>=5/(32 pi); replacing pi by upper4 gives strict final bound.
    gap=zeta*F(5,32*4)/4
    assert gap==F(5,131072)
    robust_gap=gap-2*F(1,10**6)
    assert robust_gap>0
    ekg=l2upper/2
    etag=3+F(3,2)/6
    etotal=ekg+etag
    assert etag==F(13,4) and etotal==F(262157,4)
    tag_derivative_square=9+F(3,2)
    assert tag_derivative_square<16
    source_constant=ekg+4
    assert source_constant==65540

    sources=[HERE/'proof.md',Path(__file__).resolve(),
        HERE.parent/'_shared/source_scope_and_round1287_admission.md',
        HERE.parent/'1282/proof.md',
        RESEARCH/'archive_935_955/research_note_942.md',
        RESEARCH/'archive_1120_1158/_shared/1120_1156/after1124_execution_audit.md',
        RESEARCH/'archive_1159_1160/_shared/field_bridge_after_1160.md']
    return {
        'round':1287,'status':'PASS_WITH_SCOPE','author_groups':1,'new_axioms':0,
        'method':'Fraction mass projectors, globally bounded alternating Bessel series, causal support, finite-source component differentiation; no quadrature or propagator truncation',
        'source_bindings':{p.relative_to(ROOT).as_posix():sha(p) for p in sources},
        'parameters':{'mass_squared_matrix':ms(MASS),'eigen_masses':['1/1','2/1'],
            'tag_mass':'3/1','support_radius':fs(a),'support_center_separation':fs(L),
            'time_window':[fs(t0),fs(t1)],'weyl_effect_strength':fs(zeta),
            'input':'zero configuration; first KG velocity +/- normalized C_infinity compact bump',
            'tag':'two degenerate neutral Dirac flavors; same normal positive-energy Gaussian momentum packet'},
        'causal_certificate':{'source_receiver_distance_interval':[fs(rmin),fs(rmax)],
            'proper_interior_s_squared_interval':[fs(s2min),fs(s2max)],
            'no_signal_before':fs(rmin),'lightcone_delta_cancellation':True,
            'mean_interval':['5/(32*pi)','1/(4*pi)'],
            'global_alternating_ratio_bounds':[fs(ratio_mass2),fs(ratio_mass1),fs(ratio_combined)],
            'tail_numerator_uniform_lower':fs(kernel_num_lower),
            'finite_series_diagnostics':samples},
        'normal_resource_certificate':{'e_upper_rational':fs(e_upper),
            'normalization_lower':'pi*a^3/24','bump_squared_L2_strict_upper':fs(l2upper),
            'receiver_variance_strict_upper':fs(variance),'gaussian_exponent_strict_upper':fs(exponent),
            'coherent_energy_strict_upper':fs(ekg),'tag_energy_upper':fs(etag),
            'total_energy_strict_upper':fs(etotal),'tag_derivative_L2_squared':fs(tag_derivative_square),
            'full_source_smearing_constant_strict_upper':fs(source_constant)},
        'record_certificate':{'effect':'(I+sin(zeta*phi_2(t,h_B)))/2',
            'difference':'p_plus-p_minus','whole_window_strict_lower':fs(gap),
            'decimal_diagnostic':float(gap),'nonnegative_normalized_time_average_same_gap':True,
            'optional_each_event_error':'1/1000000','optional_gap_lower':fs(robust_gap)},
        'source_diagnostics':source_checks(),
        'comparison_role_reuse':{'independent_bare_mass_squared':'1/1','positive_comparison_coefficient':'3/2',
            'comparison_matrix':ms(comp),'common_direction_kernel':True,
            'old_PSD_calibration_reused_not_new_classification':True},
        'scope':{'full_normal_unknown_Q_R_isometry':True,
            'initial_specified_receiver_algebra_and_R_product':True,
            'initial_all_fields_in_B_blind':False,
            'same_coherent_covariance_from_coupled_vacuum':True,
            'whole_receiver_Fock_tensor_factor_assumed':False,
            'actual_autonomous_local_causal_field_evolution':True,
            'finite_spatial_record_in_second_field':True,
            'full_Hilbert_source_including_mixing_tag_reference':True,
            'source_covers_unread_process':True,
            'same_mean_source_for_every_unknown_code':True,
            'same_noise_or_conditional_source_certified':False,
            'external_switch_used':False,
            'initial_preparation_factory_constructed':False,
            'actual_terminal_instrument_constructed':False,
            'global_quantum_information_cloned_to_receiver':False,
            'background_with_nonzero_source_self_gravitating':False,
            'new_Einstein_solution_or_response_error_claimed':False,
            'deletion_is_same_state_mass_quench':False,
            'new_axiom_or_stage_completion':False}}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    result=calculate()
    path=HERE/'results.json'
    if args.write:
        with path.open('x',encoding='utf-8',newline='\n') as out:
            json.dump(result,out,ensure_ascii=False,indent=2)
            out.write('\n')
        print('PASS_WITH_SCOPE: EXCLUSIVE_FIRST_WRITE')
    else:
        assert result==json.loads(path.read_text(encoding='utf-8')), 'Saved result mismatch'
        print('PASS_WITH_SCOPE: READ_ONLY_MATCH')
    print('whole_window_gap > '+result['record_certificate']['whole_window_strict_lower'])
    print('total_energy < '+result['normal_resource_certificate']['total_energy_strict_upper'])
    print('results_sha256 = '+sha(path))


if __name__=='__main__':
    main()
