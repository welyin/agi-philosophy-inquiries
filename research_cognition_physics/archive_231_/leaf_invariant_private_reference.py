"""Round 483: leaf-invariant private inputs preserve the actual tree reference.

Baseline 475. The ququart model is H4=F+sum S_ab(A) S_ab(B) n_ab.
This is a reference-moment selection rule, not a decoupled private subsystem
or an encoding of arbitrary 64-dimensional input into a 20-dimensional code.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import io
import itertools
import json
import math
from pathlib import Path
import platform
import unittest
import numpy as np
import rigid_leaf_reference_audit as frame
import sequential_tree_distance_readout_audit as reader

TARGET = Path(__file__).with_name('leaf_invariant_private_reference_results.json')
OBS = {}
DIM = 64*64*6


def short(x):
    return float(f'{float(x):.11g}')


def integer_norm_squared(x):
    assert frame.gauss_integer(x)
    return int(np.sum(x.real*x.real+x.imag*x.imag))


@lru_cache(None)
def swap_maps():
    out = {}
    for a,b in itertools.combinations(range(6),2):
        mapping=[]
        for v in range(64):
            bit_a=(v>>(5-a))&1
            bit_b=(v>>(5-b))&1
            mapping.append(v^((bit_a^bit_b)<< (5-a))^((bit_a^bit_b)<<(5-b)))
        out[a,b]=np.array(mapping,dtype=np.int64)
    return out


@lru_cache(None)
def neighbours():
    trees,_,flip=frame.system()
    indices=np.arange(DIM).reshape(64,64,6)
    result=np.empty((9,DIM),dtype=np.int64)
    for g,tree in enumerate(trees):
        src=indices[:,:,g].ravel()
        for k,target in enumerate(np.flatnonzero(flip[g])):
            result[k,src]=indices[:,:,target].ravel()
        for k,edge in enumerate(sorted(tree),4):
            p=swap_maps()[edge]
            result[k,src]=indices[p[:,None],p[None,:],g].ravel()
    return result


def h_apply(columns):
    out=np.zeros_like(columns)
    for mapping in neighbours():
        out+=columns[mapping]
    return out


def evolve(columns,time=.05,degree=18):
    term=columns.copy()
    out=term.copy()
    for k in range(1,degree+1):
        term=(-1j*time/k)*h_apply(term)
        out+=term
    return out


@lru_cache(None)
def code_columns():
    # Unnormalised Dicke columns; Gram diagonal is binomial(4,k).
    columns=[]
    weights=[]
    for inside in itertools.product(range(2),repeat=2):
        for k in range(5):
            v=np.zeros(64,dtype=np.int64)
            for leaves in itertools.product(range(2),repeat=4):
                if sum(leaves)!=k:
                    continue
                bits=[0]*6
                bits[1],bits[2]=inside
                for a,b in zip(frame.LEAVES,leaves):
                    bits[a]=b
                v[sum(b<<(5-a) for a,b in enumerate(bits))]=1
            columns.append(v)
            weights.append(math.comb(4,k))
    return np.column_stack(columns),np.array(weights,dtype=np.int64)


def source_columns(with_reference=False):
    # A=eta(0), G=|s><s|. B either |excitation at 1>, or an entangled
    # (|excitation at 1>|0_R>+|excitation at 2>|1_R>)/sqrt(2).
    out=np.zeros((64,64,6,8,2 if with_reference else 1),dtype=complex)
    for a,inside in enumerate(itertools.product(range(2),repeat=3)):
        for v3,v4 in itertools.product(range(2),repeat=2):
            bits=list(inside)+[v3,v4,0]
            ai=sum(v<<(5-j) for j,v in enumerate(bits))
            out[ai,16,:,a,0]=1j**v4
            if with_reference:
                out[ai,8,:,a,1]=1j**v4
    return out.reshape(DIM,-1),384 if with_reference else 192


@lru_cache(None)
def evolved_private_reference():
    c,den=source_columns(True)
    return evolve(c),den


def reference_moment(columns,den,observable):
    y=columns.reshape(DIM,8,2)
    weights=np.tile(observable,64*64)
    return np.einsum('bar,bas,b->rs',y,y.conj(),weights,optimize=True)/den


def graph_density(columns,den):
    y=columns.reshape(64,64,6,-1)
    return np.einsum('abgc,abhc->gh',y,y.conj(),optimize=True)/den


def local_swap(d):
    out=np.zeros((d*d,d*d),dtype=np.int64)
    for a,b in itertools.product(range(d),repeat=2):
        out[b*d+a,a*d+b]=1
    return out


class Audit(unittest.TestCase):
    def test_01_exact_leaf_characters_and_multiplicities(self):
        classes=[(1,1,1,1),(1,1,2),(2,2),(1,3),(4,)]
        sizes=[1,6,3,8,6]
        chars={'[4]':[1,1,1,1,1],'[31]':[3,1,-1,0,-1],
               '[22]':[2,0,2,-1,0],'[211]':[3,-1,-1,0,1],
               '[1111]':[1,-1,1,1,-1]}
        leaf_char=[2**len(c) for c in classes]
        multiplicities={name:sum(s*a*b for s,a,b in zip(sizes,leaf_char,char))//24
                        for name,char in chars.items()}
        self.assertEqual(multiplicities,{'[4]':5,'[31]':3,'[22]':1,'[211]':0,'[1111]':0})
        self.assertEqual(sum(multiplicities[n]*chars[n][0] for n in chars),16)
        sym_sum=np.zeros((64,64),dtype=np.int64)
        sign_sum=np.zeros_like(sym_sum)
        for p,dm,_,_ in frame.permutations():
            u=np.eye(64,dtype=np.int64)[dm]
            sym_sum+=u
            sign_sum+=frame.sign(p)*u
            self.assertEqual(np.trace(u),4*2**len(frame.cycles(p)))
        self.assertTrue(np.array_equal(sym_sum@sym_sum,24*sym_sum))
        self.assertEqual(int(np.trace(sym_sum)),24*20)
        self.assertFalse(np.any(sign_sum))
        OBS['leaf_representation'] = dict(
            conjugacy_class_sizes=sizes, four_leaf_Hilbert_character=leaf_char,
            four_leaf_irrep_multiplicities=multiplicities,
            full_B_multiplicities={'[4]':20,'[31]':12,'[22]':4},
            full_B_sector_dimensions={'[4]':20,'[31]':36,'[22]':8},
            full_B_dimension=64, sign_isotypic_projector_exactly_zero=True,
            pure_isometric_invariant_code_maximum_dimension=20,
            maximum_statement_not_extended_to_general_mixed_encoding=True)

    def test_02_explicit_code_unknown_reference_and_marginal_boundary(self):
        c,weights=code_columns()
        self.assertEqual(c.shape,(64,20))
        self.assertTrue(np.array_equal(c.T@c,np.diag(weights)))
        for _,dm,_,_ in frame.permutations():
            self.assertTrue(np.array_equal(c[dm],c))
        v=c/np.sqrt(weights)
        rng=np.random.default_rng(483)
        unknown=rng.normal(size=(20,3))+1j*rng.normal(size=(20,3))
        unknown/=np.linalg.norm(unknown)
        encoded=v@unknown
        self.assertLess(np.linalg.norm(v.T@encoded-unknown),1e-14)
        for _,dm,_,_ in frame.permutations():
            self.assertLess(np.linalg.norm(encoded[dm]-encoded),1e-14)
        b=np.eye(64)[:,16]
        self.assertLess(np.linalg.norm(v@(v.T@b)-b),1e-14)
        self.assertEqual(F(1,6),F(int(sum(b[x] for x in range(64) if x.bit_count()==1))**2,6))
        trans=next(dm for p,dm,_,_ in frame.permutations() if frame.cycles(p)==(1,1,2))
        overlap=F(int(np.sum(trans==np.arange(64))),64)**2
        self.assertEqual(overlap,F(1,4))
        OBS['safe_private_input_contract'] = dict(
            isometry='C^2_(B1) x C^2_(B2) x Sym^4(C^2)_leaves -> (C^2)^6',
            input_dimension=20, entire_unknown_private_input_and_R_preserved_by_isometry=True,
            source_A_and_graph_must_still_be_independent_of_BR=True,
            joint_BR_leaf_conjugation_invariance_required=True,
            B_I64_over64_without_R_is_allowed_and_not_restricted_to_pure_code=True,
            maximally_entangled_64_by64_BR_has_same_marginal_but_transposition_fidelity=str(overlap),
            invariant_B_marginal_does_not_imply_joint_BR_hypothesis=True,
            no_claim_that_every_noninvariant_BR_has_failed_reference=True,
            internal_one_excitation_codeword_full_Sym6_weight='1/6',
            arbitrary_unknown_64_dimensional_input_not_compressed_to_20=True,
            code_is_initial_reference_safe_support_not_a_dynamically_noiseless_subsystem=True)

    def test_03_exact_H4_symmetry_and_selection_transfer(self):
        near=neighbours()
        indices=np.arange(DIM).reshape(64,64,6)
        for _,dm,gm,_ in frame.permutations():
            perm=indices[dm[:,None,None],dm[None,:,None],gm[None,None,:]].ravel()
            self.assertTrue(np.array_equal(np.sort(perm[near],axis=0),np.sort(near[:,perm],axis=0)))
        self.assertEqual(near.shape,(9,DIM))
        s4=local_swap(4)
        for pauli in frame.P:
            local=np.kron(pauli,np.eye(2))
            total=np.kron(local,np.eye(4))+np.kron(np.eye(4),local)
            self.assertTrue(np.array_equal(s4@total,total@s4))
        n=frame.twirled_numerators()
        for numerator in n:
            projected=sum(char*numerator[np.ix_(dm,dm)] for _,dm,_,char in frame.permutations())
            self.assertFalse(np.any(projected))
        for a,b in itertools.combinations(frame.LEAVES,2):
            q=3*frame.distance(a,b)-8
            self.assertTrue(np.array_equal(sum(char*q[gm] for _,_,gm,char in frame.permutations()),12*q))
        OBS['all_time_reference_selection_rule'] = dict(
            actual_ABG_dimension=DIM, nonzero_H_row_terms_with_multiplicity=9,
            full_basis_diagonal_leaf_S4_covariance_verified_for_24_elements=True,
            independent_collective_A_SU2_commutator_verified_locally=True,
            H4_norm_upper_at_unit_couplings=9,
            general_H4_norm_upper='4*abs(kappa)+5*abs(J)',
            initial_state='eta_A(r) x sigma_BR x gamma_G',
            sigma_BR_joint_leaf_invariant=True, gamma_G_independent_leaf_invariant=True,
            E22_of_A_twirl_initial_state_exactly_zero=True,
            all_time_conditional_leaf_distance_moment='(8/3)*rho_R',
            all_time_matching_class_reference_moment='rho_R/3',
            arbitrary_real_kappa_J_time_covered_by_proof=True,
            no_requirement_that_BR_remain_independent_during_evolution=True)

    def test_04_private_state_and_code_really_change(self):
        c,den=source_columns(False)
        self.assertEqual(integer_norm_squared(c),den)
        hc=h_apply(c)
        x=hc.reshape(64,64,6,8)
        away=x.copy()
        away[:,16]=0
        coefficient_leave=F(integer_norm_squared(away),den)
        sym=sum(x[:,dm] for _,dm,_,_ in frame.permutations())
        leak_numerator=24*x-sym
        # Prior exact-integer upper bound, covering all Gaussian arithmetic.
        self.assertLess(2*DIM*8*(24*18)**2,2**53)
        coefficient_leak=F(integer_norm_squared(leak_numerator),den*24**2)
        self.assertEqual(coefficient_leave,F(3))
        self.assertEqual(coefficient_leak,F(11,8))
        t=F(1,2000)
        tail_over_t2=972*t/(1-F(18,4)*t)
        self.assertLess(tail_over_t2,F(1,2))
        y=evolve(c)
        moved=y.reshape(64,64,6,8).copy()
        moved[:,16]=0
        moved_probability=np.linalg.norm(moved)**2/den
        sym=sum(y.reshape(64,64,6,8)[:,dm] for _,dm,_,_ in frame.permutations())/24
        leak_probability=np.linalg.norm(y.reshape(64,64,6,8)-sym)**2/den
        self.assertGreater(moved_probability,.001)
        self.assertGreater(leak_probability,.001)
        OBS['nontrivial_private_dynamics'] = dict(
            witness_A_source='eta(0)', witness_B='one excitation at internal site 1; all other B=0',
            graph='uniform pure six-tree state',
            leave_original_B_vector_t_squared_coefficient=str(coefficient_leave),
            leave_initial_20_dimensional_code_t_squared_coefficient=str(coefficient_leak),
            remainder_upper='(18*t)^3/(6*(1-18*t/4))',
            strictly_positive_interval='0 < t <= 1/2000',
            leave_B_probability_lower='(5/2)*t^2', code_leakage_probability_lower='(7/8)*t^2',
            interval_tail_over_t_squared_upper=str(tail_over_t2),
            numerical_time='1/20', numerical_leave_B_probability=short(moved_probability),
            numerical_code_leakage_probability=short(leak_probability),
            no_B_freezing_or_factorisation_of_H4_claim=True)

    def test_05_true_H4_evolution_with_entangled_private_reference(self):
        y,den=evolved_private_reference()
        ref=reference_moment(y,den,np.ones(6))
        errors=[]
        for a,b in itertools.combinations(frame.LEAVES,2):
            errors.append(np.max(abs(reference_moment(y,den,frame.distance(a,b))-(8/3)*ref)))
        matching_errors=[]
        for b in (3,4,5):
            matching_errors.append(np.max(abs(reference_moment(y,den,3-frame.distance(0,b))-ref/3)))
        self.assertLess(np.max(abs(ref-np.eye(2)/2)),2e-12)
        self.assertLess(max(errors+matching_errors),2e-12)
        c,_=source_columns(True)
        inverse=evolve(y,time=-.05)
        self.assertLess(np.linalg.norm(inverse-c)/math.sqrt(den),2e-12)
        tail=F(9,20)**19/math.factorial(19)/(1-F(9,400))
        self.assertLess(tail,F(1,10**22))
        OBS['conditional_reference_crosscheck'] = dict(
            B_R_input='(|B excitation 1>|0_R> + |B excitation 2>|1_R>)/sqrt(2)',
            independent_A_mixed_source_rank=8,
            full_active_vector_dimension=DIM, reference_dimension=2,
            time='1/20', vector_Taylor_degree=18,
            operator_Taylor_tail_bound=str(tail),
            finite_Taylor_bound_does_not_include_floating_roundoff=True,
            maximum_leaf_reference_matrix_residual=short(max(errors)),
            maximum_matching_reference_matrix_residual=short(max(matching_errors)),
            unknown_global_information_preserved_by_unitarity_not_by_private_marginal_freezing=True,
            inverse_only_used_as_a_numerical_unitarity_check=True)

    def test_06_ququart_probe_intertwines_actual_data_instrument(self):
        v=np.array([[1,0],[0,0],[0,1],[0,0]],dtype=np.int64)
        pair=np.kron(v,v)
        self.assertTrue(np.array_equal(local_swap(4)@pair,pair@local_swap(2)))
        z4=np.kron(np.diag([1,-1]),np.eye(2,dtype=np.int64))
        self.assertTrue(np.array_equal(z4@v,v@np.diag([1,-1])))
        # All 384 embedded basis columns have the same nine transitions as H2.
        embedding=np.array([(a*64)*6+g for a in range(64) for g in range(6)])
        h2=frame.system()[1]
        for j,raw in enumerate(embedding):
            counts=np.bincount(neighbours()[:,raw],minlength=DIM)
            self.assertEqual(int(np.sum(counts)),9)
            self.assertTrue(np.array_equal(counts[embedding],h2[:,j]))
        y,den=evolved_private_reference()
        graph=graph_density(y,den)
        probe=F(1,65536)
        records=[]
        for edge in ((1,0),(1,3),(1,4),(1,5)):
            minus,plus,_=reader.instrument(edge,probe)
            p=[float(np.trace(np.einsum('ghij,ij->gh',op,graph)).real) for op in (minus,plus)]
            self.assertLess(abs(sum(p)-1),2e-12)
            self.assertGreaterEqual(min(p),-2e-12)
            records.append(dict(edge=list(edge),probabilities=[short(q) for q in p]))
        path_values=[]
        for path in ((0,1,3),(0,2,3)):
            state=graph.copy()
            for edge in zip(path[:-1],path[1:]):
                minus,plus,_=reader.instrument(edge,probe)
                state=np.einsum('ghij,ij->gh',(plus-minus)/float(probe),state)
            path_values.append(float(np.trace(state).real))
        measured_leaf=3-sum(path_values)
        single_bias=(18*probe)**2/(2*(1-6*probe)*probe)
        path_bias=2*single_bias+single_bias*single_bias
        self.assertLess(abs(measured_leaf-8/3),float(2*path_bias)+1e-5)
        leaf_error=F(1,1000)
        path_error=leaf_error/2
        leaf_probe=path_error/(32*81)
        leaf_gamma=path_error*leaf_probe**2/8
        d=(18*leaf_probe)**2/(2*(1-6*leaf_probe)*leaf_probe)
        self.assertLess(2*d+d*d+leaf_gamma/leaf_probe**2,path_error/2)
        leaf_copies=math.ceil(80/(leaf_probe**4*path_error**2))
        self.assertGreaterEqual(leaf_copies*leaf_probe**4*path_error**2/8,10)
        a=F(1,1000)
        h=a/(16*81)
        gamma=a*h/4
        x=18*h
        self.assertLessEqual((x*x/(2*(1-x/3))+gamma)/h,a/2)
        n=math.ceil(80/(h*h*a*a))
        self.assertGreaterEqual(n*h*h*a*a/8,10)
        OBS['physical_d4_readout'] = dict(
            local_embedding='|a> -> |a>_A |0>_B', fresh_B_probe_vacuum_dimension=1,
            true_ququart_SWAP_and_binary_port_effect_intertwining=True,
            exact_full_384_column_generator_intertwining=True,
            frozen_instrument_round=472, actual_CP_outcome_records=records,
            two_actual_signed_paths=[[0,1,3],[0,2,3]],
            two_path_signed_estimates=[short(q) for q in path_values],
            actual_leaf_D03_estimate=short(measured_leaf),
            actual_current_leaf_mean_estimation_error=str(leaf_error), leaf_path_probe_time=str(leaf_probe),
            ideal_8_over3_comparison_error_bound='a + delta_source/2',
            source_preparation_and_natural_evolution_error_is_a_separate_complete_process_budget=True,
            whole_two_step_path_diamond_budget=str(leaf_gamma),
            independent_copies_per_path=str(leaf_copies),
            joint_two_path_failure_bound='4*exp(-10) < 1/100',
            all_old_AB_saved_in_internal_noninteracting_S_graph_G_stays_active=True,
            S_G_R_correlations_retained_final_probe_can_disturb_G=True,
            fresh_independent_A_probe_and_B_vacuum_are_resources=True,
            operator_valued_instrument_identity_valid_for_arbitrary_stored_reference=True,
            illustrative_each_adjacency_mean_error=str(a), probe_time=str(h),
            complete_storage_handoff_and_probe_diamond_budget=str(gamma),
            independent_copies_per_edge=str(n), joint_four_edge_failure_bound='8*exp(-10) < 1/100',
            no_free_state_reset_or_cloning_or_autonomous_reader_generation=True)


def run():
    OBS.clear()
    stream=io.StringIO()
    result=unittest.TextTestRunner(stream=stream).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise RuntimeError(stream.getvalue())
    return dict(round=483,baseline_round=475,tests_run=result.testsRun,
        failures=len(result.failures),errors=len(result.errors),
        python=platform.python_version(),numpy=np.__version__,observations=OBS,
        scope=dict(all_time_reference_safe_joint_leaf_invariant_private_states=True,
            twenty_dimensional_pure_isometric_code_is_maximal_for_this_invariance_contract=True,
            arbitrary_mixed_encodings_not_classified=True,
            private_subsystem_may_evolve_entangle_and_leave_initial_code=True,
            invariant_private_marginal_not_substitute_for_joint_unknown_reference_condition=True,
            source_and_graph_independence_preparation_and_probe_permissions_remain_inputs=True,
            no_full_d4_state_quotient_or_physical_dimension_selection=True,
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
