"""Round 462: bilateral activation from fixed raw pair exchanges and node roles.

No per-edge controller or state-dependent classical Hamiltonian is used.
Roles, possible contacts, initial preparations and readout permissions remain
explicit model inputs. The gate is eligibility, not a physical adjacency.
"""
import argparse
from fractions import Fraction
import io
import itertools
import json
import math
from pathlib import Path
import platform
import unittest
import numpy as np
import exchange_relation_audit as old
from recursive_exchange_interface_audit import permutation

TARGET=Path(__file__).with_name('bilateral_role_exchange_audit_results.json')
OBS={}


def edges_for(n):
    internal=[(3*a,3*a+1) for a in range(n)]+[(3*a+1,3*a+2) for a in range(n)]
    external=[(3*a+i,3*b+j) for a in range(n) for b in range(a+1,n)
              for i in (0,1) for j in (0,1)]
    return internal,external


def raw_h(n):
    internal,external=edges_for(n)
    return sum(old.swap(3*n,*edge) for edge in internal+external)


def triplet(n,a):
    return (np.eye(2**(3*n))+old.swap(3*n,3*a,3*a+1))/2


def local_response(n,a):
    return sum(old.swap(3*n,3*a+i,3*a+j) for i,j in [(0,1),(0,2),(1,2)])/3


def cross(n,a,b):
    return sum(old.swap(3*n,3*a+i,3*b+j) for i in (0,1) for j in (0,1))-2*np.eye(2**(3*n))


def integer_columns():
    """Initial real code columns before forming source Y+/Y-.

    First eight columns: source L0, B L1, C L0; norm squared 24.
    Last eight columns: source L1, B L1, C L0; norm squared 72.
    All eight source/receiver G choices are retained.
    """
    v=old.encoding().real
    groups=[[],[]]
    for ga,gb,gc in itertools.product(range(2),repeat=3):
        b=np.rint(math.sqrt(6)*v[:,2*gb+1]).astype(np.int64)
        c=np.rint(math.sqrt(2)*v[:,2*gc]).astype(np.int64)
        for l,d in [(0,2),(1,6)]:
            a=np.rint(math.sqrt(d)*v[:,2*ga+l]).astype(np.int64)
            groups[l].append(np.kron(np.kron(a,b),c))
    return np.column_stack(groups[0]+groups[1]).astype(object)


def integer_signal_certificate():
    degree=40
    denominator=2**degree*math.factorial(degree)
    internal,external=edges_for(3)
    permutations=[permutation(9,*e) for e in internal+external]
    initial=integer_columns()
    real=np.zeros(initial.shape,dtype=object)
    imag=np.zeros(initial.shape,dtype=object)
    power=initial.copy()
    small=[]
    for k in range(degree+1):
        if k<8:
            small.append(power.copy())
        coefficient=denominator//(2**k*math.factorial(k))
        if k%4==0:
            real+=coefficient*power
        elif k%4==1:
            imag-=coefficient*power
        elif k%4==2:
            real-=coefficient*power
        else:
            imag+=coefficient*power
        if k<degree:
            power=sum(power[p] for p in permutations)
    epsilon=Fraction(3**9*9**41,math.factorial(41))
    probability_difference_error=4*epsilon+2*epsilon**2
    sqrt_low=Fraction(17320508075688772,10**16)
    sqrt_high=Fraction(17320508075688773,10**16)
    assert sqrt_low**2<3<sqrt_high**2
    rows=[]
    for node in (1,2):
        qs=[permutation(9,3*node+a,3*node+b) for a,b in [(0,1),(0,2),(1,2)]]
        qreal=sum(real[p,8:] for p in qs)
        qimag=sum(imag[p,8:] for p in qs)
        numerator=int(np.sum(real[:,:8]*qimag-imag[:,:8]*qreal))
        coefficient=Fraction(numerator,288*denominator**2)
        assert coefficient>0
        lower=coefficient/sqrt_high-probability_difference_error
        upper=coefficient/sqrt_low+probability_difference_error
        derivatives=[]
        for k in (1,3,5,7):
            integer=0
            for r in range(k+1):
                s=k-r
                sign=[0,1,0,-1][(r-s)%4]
                integer+=math.comb(k,r)*sign*int(np.sum(
                    small[r][:,:8]*sum(small[s][p,8:] for p in qs)))
            derivatives.append(integer)
        rows.append(dict(node=node, polynomial_numerator=str(numerator),
            polynomial_denominator=str(288*denominator**2),
            derivative_numerators_k_1_3_5_7=derivatives,
            lower=lower,upper=upper))
    return rows,epsilon,probability_difference_error


class Audit(unittest.TestCase):
    def close(self,a,b,tolerance=8e-12):
        self.assertLess(float(np.linalg.norm(a-b)),tolerance)

    def test_01_common_rule_and_raw_exchange_resource_count(self):
        h=raw_h(3)
        self.assertEqual(len(edges_for(3)[0]),6)
        self.assertEqual(len(edges_for(3)[1]),12)
        # Every node permutation preserves the one fixed raw Hamiltonian.
        for a,b in [(0,1),(1,2)]:
            idx=np.arange(512)
            for port in range(3):
                idx=idx[permutation(9,3*a+port,3*b+port)]
            self.close(h[idx,:],h[:,idx])
        # Uniformity is at the role-defined node-pair level, not all raw pairs.
        self.assertLess(len(edges_for(3)[0])+len(edges_for(3)[1]),math.comb(9,2))
        OBS['rule_and_inputs']=dict(raw_qubits_per_node=3,
            internal_pairs_per_node=2,raw_pairs_per_potential_node_pair=4,
            node_pair_specific_weights=False,independent_edge_registers=0,
            full_potential_contact_count='4 binomial(N,2)',
            total_raw_pairs='2N + 4 binomial(N,2)',
            node_permutation_covariant=True,
            grouping_roles_contacts_preparation_and_readout_are_inputs=True,
            full_raw_permutation_invariance_claimed=False)

    def test_02_exact_bilateral_gate_and_internal_update(self):
        n=2
        ta,tb=triplet(n,0),triplet(n,1)
        k=cross(n,0,1)
        self.close(k,ta@tb@k@ta@tb)
        self.close(k@ta,ta@k)
        self.close(k@tb,tb@k)
        h=raw_h(n)
        local=old.swap(6,0,1)+old.swap(6,1,2)
        self.close(h@ta-ta@h,local@ta-ta@local)
        self.assertGreater(np.linalg.norm(h@ta-ta@h),1)
        v=old.encoding()
        t=(np.eye(8)+old.swap(3,0,1))/2
        loc=old.swap(3,0,1)+old.swap(3,1,2)
        self.close(v.conj().T@(1j*(loc@t-t@loc))@v,
                   -math.sqrt(3)*np.kron(np.eye(2),old.PAULI[1])/2)
        # Activation dynamics are affected at the next derivative, despite
        # direct commutation of cross interactions with the record.
        mixed=k@(local@ta-ta@local)-(local@ta-ta@local)@k
        self.assertGreater(np.linalg.norm(mixed),1)
        OBS['bilateral_activation']=dict(
            centered_cross_generator='2 J_port_A dot J_port_B',
            eligibility='T_A T_B', both_triplets_required=True,
            cross_terms_do_not_flip_local_triplet_population=True,
            same_internal_exchange_updates_population=True,
            compressed_initial_record_derivative='-sqrt(3) Y_L / 2',
            feedback_backaction_at_higher_derivative=True,
            expectation_value_substituted_into_H=False)

    def test_03_all_state_actual_response_current(self):
        # A has a two-spin port and its third constituent; B contributes its
        # two-spin port. The spectator third B spin is unnecessary here.
        swaps={(a,b):old.swap(5,a,b).real.astype(np.int64)
               for a in range(5) for b in range(a+1,5)}
        k=sum(swaps[a,b] for a in (0,1) for b in (3,4))-2*np.eye(32,dtype=np.int64)
        q=swaps[0,1]+swaps[0,2]+swaps[1,2]  # 3 * spin-3/2 response of A
        d=k@q-q@k  # current = i D / 3
        d2=-d@d
        self.assertTrue(np.array_equal(d2@(d2-8*np.eye(32,dtype=np.int64))@
                         (d2-20*np.eye(32,dtype=np.int64)),np.zeros((32,32),dtype=np.int64)))
        self.assertEqual(int(np.trace(d2)),192)
        self.assertEqual(int(np.trace(d2@d2)),3456)
        self.assertEqual(int(np.trace(d2@(d2-8*np.eye(32,dtype=np.int64)))),1920)
        gate=(np.eye(32)+swaps[0,1])@(np.eye(32)+swaps[3,4])/4
        self.close(gate@d@gate,d)
        self.assertAlmostEqual(np.linalg.norm(1j*d/3,2),2*math.sqrt(5)/3)
        OBS['actual_response_current']=dict(
            local_effect='total j=3/2 projector of the three actual node qubits',
            exact_integer_polynomial='D2 (D2-8I) (D2-20I)=0, D2=-D^2',
            exact_trace_D2=192,exact_trace_D2_squared=3456,
            largest_eigenvalue_present_integer_witness=1920,
            exact_current_norm='2 sqrt(5)/3',gate_support='T_A T_B',
            arbitrary_joint_state_speed_bound='|d<E_A>/dt| <= (2 sqrt(5)/3) sum_(B!=A) |g_AB| <T_A T_B>',
            local_internal_terms_commute_with_response=True,
            bound_is_necessary_not_sufficient_for_response=True,
            no_external_axis_needed_for_effect=True)

    def test_04_exact_preferential_response_certificate(self):
        rows,epsilon,error=integer_signal_certificate()
        self.assertGreater(rows[0]['lower'],Fraction(1,8))
        self.assertGreater(rows[1]['lower'],0)
        self.assertLess(rows[1]['upper'],Fraction(1,200))
        self.assertEqual(rows[0]['derivative_numerators_k_1_3_5_7'],[0,4608,-138240,1741824])
        self.assertEqual(rows[1]['derivative_numerators_k_1_3_5_7'],[0,0,23040,-4161024])
        self.assertGreater(Fraction(1,8)-6*Fraction(1,2048),Fraction(3,25))
        self.assertLess(Fraction(1,200)+6*Fraction(1,2048),Fraction(1,100))
        for row in rows:
            row['certified_lower']=float(row.pop('lower'))
            row['certified_upper']=float(row.pop('upper'))
        OBS['integer_response_certificate']=dict(raw_qubits=9,time='1/2',
            series_degree=40,hamiltonian_norm_upper=18,
            norm_tail_bound=str(epsilon),contrast_error_bound=str(error),
            rows=rows,source='Y+ versus Y- in original A private L',
            receiver_B='L1 (triplet port initially)',receiver_C='L0 (singlet port initially)',
            all_three_G_inputs_in_this_witness='independent I/2',
            initial_full_states_have_no_inter_node_correlations=True,
            short_time_B='8 t^3/(3 sqrt(3)) + O(t^5)',
            short_time_C='2 t^5/(3 sqrt(3)) + O(t^7)',
            window_half_width='1/2048',B_contrast_lower_on_window='3/25',
            absolute_C_contrast_upper_on_window='1/100',
            old_three_spin_codes_claimed_preserved=False)

    def test_05_independent_raw_evolution_and_unknown_reference(self):
        h=raw_h(3)
        v=old.encoding()
        rows=[]
        u=old.evolve(h,.5)
        for sign in (1,-1):
            columns=[]
            for ga,gb,gc in itertools.product(range(2),repeat=3):
                a=(v[:,2*ga]+sign*1j*v[:,2*ga+1])/math.sqrt(2)
                columns.append(np.kron(np.kron(a,v[:,2*gb+1]),v[:,2*gc]))
            out=u@np.column_stack(columns)
            rows.append([float(f'{np.vdot(out,local_response(3,a)@out).real/8:.12g}')
                         for a in (1,2)])
        self.assertGreater(rows[1][0]-rows[0][0],1/8)
        self.assertLess(rows[1][1]-rows[0][1],1/200)
        # Whole unknown source/receiver states and references are retained;
        # the above special witness is not claimed G-independent in general.
        embedding=np.kron(np.kron(v,v),v)
        rng=np.random.default_rng(46205)
        psi=rng.normal(size=(64,3))+1j*rng.normal(size=(64,3))
        psi/=np.linalg.norm(psi)
        initial=embedding@psi
        final=u@initial
        self.close(u.conj().T@final,initial,tolerance=1e-10)
        self.close(final.conj().T@final,initial.conj().T@initial,tolerance=1e-10)
        OBS['raw_reproduction']=dict(
            receiver_probabilities_Y_plus=rows[0],receiver_probabilities_Y_minus=rows[1],
            arbitrary_64_dimensional_old_input_and_reference_preserved_globally=True,
            numerical_reference_dimension=3,
            source_or_receiver_reduced_quantum_states_unchanged=False,
            response_witness_not_general_gauge_independent_channel=True,
            effective_interface_forced_to_single_total_j=False)

    def test_06_configuration_limits_correlations_and_size_budget(self):
        for bits in itertools.product((0,1),repeat=3):
            a,b,c=bits
            self.assertLessEqual(a*b+b*c-b,a*c)
        # Per-record configuration the eligibility graph is one active clique.
        for bits in itertools.product((0,1),repeat=5):
            edges=[(a,b) for a in range(5) for b in range(a+1,5) if bits[a]*bits[b]]
            active=sum(bits)
            self.assertEqual(len(edges),math.comb(active,2))
        # Correlations can mix K2 branches into a non-clique marginal support.
        configurations=[(1,1,0),(0,1,1)]
        p=[sum(bits[a]*bits[b] for bits in configurations)/2 for a,b in [(0,1),(1,2),(0,2)]]
        self.assertEqual(p,[.5,.5,0])
        # Each port constituent touches its sibling, optionally local memory,
        # and 2(N-1) external constituents; the worst row load is 2N.
        sizes=[]
        for n in (2,3,4,8):
            internal,external=edges_for(n)
            degree=[0]*(3*n)
            for a,b in internal+external:
                degree[a]+=1
                degree[b]+=1
            self.assertEqual(max(degree),2*n)
            sizes.append(dict(nodes=n,raw_qubits=3*n,raw_exchange_contacts=len(internal+external),
                maximum_unit_coefficient_load=2*n,
                maximum_common_strength_at_kappa_one=str(Fraction(1,2*n))))
        OBS['configuration_and_resource_boundaries']=dict(
            every_joint_record_configuration_is_one_active_clique=True,
            arbitrary_correlated_marginal_support_claimed_clique=False,
            correlated_path_gate_probabilities=['1/2','1/2','0'],
            exact_all_state_inequality='p_AB+p_BC-p_B <= p_AC',
            deterministic_two_edge_path_without_third_edge_not_supported=True,
            independent_partner_choices_implemented=False,
            no_edge_register_does_not_remove_quadratic_contacts=True,
            finite_local_absolute_coefficient_budget_is_extra_input=True,
            example_sizes=sizes,
            coupling_budget_identified_with_thermodynamic_work=False,
            cognitive_axioms_or_physical_space_derived=False)


def run():
    OBS.clear()
    output=io.StringIO()
    result=unittest.TextTestRunner(stream=output).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise RuntimeError(output.getvalue())
    return dict(round=462,baseline_round=461,tests_run=result.testsRun,
        failures=len(result.failures),errors=len(result.errors),
        python=platform.python_version(),numpy=np.__version__,observations=OBS,
        scope=dict(raw_two_body_bilateral_activation_with_node_roles=True,
            common_rule_without_edge_registers=True,
            actual_preferential_response_certified=True,
            arbitrary_unknown_joint_information_retained=True,
            pair_specific_partner_choice_or_contacts_generated=False,
            node_role_partition_derived_from_cognitive_principles=False,
            full_recursive_cognition_proved=False,
            physical_space_dimension_derived=False,
            full_GR_goal_completed=False,phase_closure_triggered=False))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dry-run',action='store_true')
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    result=run()
    if args.check:
        assert result==json.loads(TARGET.read_text(encoding='utf-8'))
    elif not args.dry_run:
        with TARGET.open('x',encoding='utf-8',newline='\n') as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False,indent=2))
