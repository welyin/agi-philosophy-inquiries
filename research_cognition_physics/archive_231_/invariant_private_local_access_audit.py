"""Round 484: private storage capacity versus single-leaf closed-code access.

Scientific inputs are the frozen 475 reference selection rule and 471 interface
scope. Permutation-invariant quantum codes are established constructions.
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

TARGET=Path(__file__).with_name('invariant_private_local_access_audit_results.json')
OBS={}
LEAVES=(0,3,4,5)


def cycles(p):
    seen=set();count=0
    for k in range(len(p)):
        if k in seen:continue
        count+=1;j=k
        while j not in seen:seen.add(j);j=p[j]
    return count

def sign(p):
    return (-1)**sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p)))

@lru_cache(None)
def leaf_columns(q):
    multisets=list(itertools.combinations_with_replacement(range(q),4))
    column={value:j for j,value in enumerate(multisets)}
    c=np.zeros((q**4,len(multisets)),dtype=np.int64)
    for row,word in enumerate(itertools.product(range(q),repeat=4)):
        c[row,column[tuple(sorted(word))]]=1
    weights=c.sum(axis=0)
    assert np.array_equal(c.T@c,np.diag(weights))
    return c,weights,multisets

@lru_cache(None)
def full_columns(q):
    leaf,weights,multisets=leaf_columns(q)
    result=np.zeros((q**6,q*q*len(multisets)),dtype=np.int64)
    for a,b in itertools.product(range(q),repeat=2):
        for row,word in enumerate(itertools.product(range(q),repeat=4)):
            full=(word[0],a,b,word[1],word[2],word[3])
            index=sum(v*q**(5-j) for j,v in enumerate(full))
            col=(a*q+b)*len(multisets)+int(np.flatnonzero(leaf[row])[0])
            result[index,col]=1
    return result,np.tile(weights,q*q)

def image(q,n,sites,p):
    target=[]
    for word in itertools.product(range(q),repeat=n):
        out=list(word)
        for j,site in enumerate(sites):out[sites[p[j]]]=word[site]
        target.append(sum(v*q**(n-1-j) for j,v in enumerate(out)))
    return np.array(target)

def local(v,a,site,q,n):
    # Supports all columns at once, with integer operations when inputs are int.
    shape=(q,)*n+(v.shape[1],)
    x=np.moveaxis(v.reshape(shape),site,0)
    y=np.tensordot(a,x,axes=(1,0))
    return np.moveaxis(y,0,site).reshape(v.shape)

def matrix_units(q):
    out=[]
    for a,b in itertools.product(range(q),repeat=2):
        e=np.zeros((q,q),dtype=np.int64);e[a,b]=1;out.append(e)
    return out

def rational_rank(rows):
    pivots={}
    for raw in rows:
        row=[F(int(v)) for v in raw]
        for k in sorted(pivots):
            if row[k]:
                value=row[k];row=[a-value*b for a,b in zip(row,pivots[k])]
        nonzero=next((j for j,v in enumerate(row) if v),None)
        if nonzero is not None:
            value=row[nonzero];pivots[nonzero]=[a/value for a in row]
    return len(pivots)

def swap(d):
    out=np.zeros((d*d,d*d),dtype=np.int64)
    for a,b in itertools.product(range(d),repeat=2):out[b*d+a,a*d+b]=1
    return out

def short(x):
    return float(f'{float(x):.13g}')

class Audit(unittest.TestCase):
    def test_01_capacity_characters_and_true_probe_embedding(self):
        permutations=list(itertools.permutations(range(4)))
        capacities=[]
        for q in (2,3,4):
            trivial=F(sum(q**cycles(p) for p in permutations),24)
            alternating=F(sum(sign(p)*q**cycles(p) for p in permutations),24)
            self.assertEqual(trivial,math.comb(q+3,4))
            self.assertEqual(alternating,math.comb(q,4) if q>=4 else 0)
            capacities.append(dict(private_local_dimension=q,symmetric_leaf_dimension=int(trivial),
                symmetric_pure_code_dimension=int(q*q*trivial),
                alternating_pure_code_dimension=int(q*q*alternating)))
        self.assertEqual(capacities[1]['symmetric_pure_code_dimension'],135)
        self.assertGreaterEqual(135,64)
        probe=np.zeros((6,2),dtype=np.int64);probe[0,0]=probe[3,1]=1
        self.assertTrue(np.array_equal(swap(6)@np.kron(probe,probe),np.kron(probe,probe)@swap(2)))
        OBS['capacity_and_reference_contract']=dict(records=capacities,old_private_input_dimension=64,
            six_private_qutrit_hardware_dimension=729,full_local_A_times_B_dimension=6,
            pure_pointwise_invariant_contract_only=True,
            sign_sector_not_silently_discarded=True,
            not_maximum_of_all_mixed_or_reference_safe_encodings=True,
            reference_selection_rule_extends_with_joint_BR_leaf_invariance=True,
            A_BR_G_initial_independence_required=True,
            known_new_B_zero_probe_intertwines_complete_472_instrument=True,
            old_unknown_private_data_first_moved_to_internal_storage=True,
            global_encoder_and_qutrit_hardware_are_extra_inputs=True)

    def test_02_explicit_64_dimensional_isometry_and_unknown_R(self):
        c,weights=full_columns(3)
        self.assertEqual(c.shape,(729,135))
        self.assertTrue(np.array_equal(c.T@c,np.diag(weights)))
        for p in itertools.permutations(range(4)):
            self.assertTrue(np.array_equal(c[image(3,6,LEAVES,p)],c))
        chosen=c[:,:64];w=weights[:64]
        v=chosen/np.sqrt(w)[None,:]
        self.assertLess(np.max(abs(v.T@v-np.eye(64))),3e-15)
        self.assertEqual(chosen[0,0],1)
        self.assertEqual(int(w[0]),1)
        rng=np.random.default_rng(484)
        unknown=rng.normal(size=(64,3))+1j*rng.normal(size=(64,3));unknown/=np.linalg.norm(unknown)
        encoded=v@unknown;decoded=v.T@encoded
        reference_error=np.linalg.norm(encoded.conj().T@encoded-unknown.conj().T@unknown)
        self.assertLess(reference_error,2e-14)
        self.assertLess(np.linalg.norm(decoded-unknown),2e-14)
        # |i>source|0>target and |0>source V|i>target are orthogonal for i>0;
        # 63 pairwise swaps plus a fixed i=0 vector extend transfer to a unitary.
        self.assertTrue(np.array_equal(chosen[0,1:],np.zeros(63,dtype=np.int64)))
        OBS['explicit_encoding']=dict(integer_column_shape=list(c.shape),
            selected_columns=list(range(64)),selected_integer_Gram_diagonal=[int(x) for x in w],
            exact_all_24_leaf_permutations_fix_every_column=True,
            exact_isometry_definition='V=C_selected diag(Gram)^(-1/2)',
            numerical_unknown_R_dimension=3,reference_error=short(reference_error),
            decode_error=short(np.linalg.norm(decoded-unknown)),
            coherent_transfer_unitary_exists_by_63_orthogonal_pair_swaps=True,
            input_register_dimension=64,blank_target_dimension=729,
            transfer_not_derived_from_existing_Hamiltonian=True,
            no_unknown_input_discarded_or_cloned=True)

    def test_03_same_leaf_effects_and_direct_identity_error(self):
        c,w,_=leaf_columns(3)
        checks=0
        for e in matrix_units(3):
            first=c.T@local(c,e,0,3,4)
            for site in (1,2,3):
                self.assertTrue(np.array_equal(first,c.T@local(c,e,site,3,4)))
                checks+=1
        # The equal reductions theorem is for all logical inputs and references.
        # Two different old leaf states cannot both be the same new marginal.
        sigma=np.diag([.5,.5,0]);zero=np.diag([1,0,0]);one=np.diag([0,1,0])
        distances=[np.linalg.svd(sigma-target,compute_uv=False).sum()/2 for target in (zero,one)]
        self.assertEqual(distances,[.5,.5])
        OBS['direct_leaf_readout_boundary']=dict(exact_matrix_unit_comparisons=checks,
            all_four_single_leaf_compressions_equal=True,
            equality_extends_to_operator_valued_R_marginals=True,
            incompatible_input_pair=['|0><0|','|1><1|'],
            universal_direct_marginal_worst_case_lower_bound='1/2',
            bound_is_for_same_physical_identification_not_arbitrary_decoders=True,
            equality_does_not_imply_all_compressed_observables_commute=True)

    def test_04_closed_single_leaf_algebra(self):
        ranks=[]
        for q in (2,3):
            c,w,_=leaf_columns(q);terms=[]
            for e in matrix_units(q):
                at0=local(c,e,0,q,4)
                terms.append(np.concatenate([(local(c,e,j,q,4)-at0).reshape(-1) for j in (1,2,3)]))
            rows=np.stack(terms,axis=1)
            rank=rational_rank(rows)
            self.assertEqual(rank,q*q-1)
            ranks.append(dict(q=q,exact_constraint_rank=rank,scalar_solution_dimension=1))
        ghz=np.zeros((16,2),dtype=np.int64);ghz[0,0]=ghz[15,1]=1
        z=np.diag([1,-1]);x=np.array([[0,1],[1,0]])
        for j in range(4):self.assertTrue(np.array_equal(local(ghz,z,j,2,4),ghz@z))
        self.assertTrue(np.array_equal(ghz.T@local(ghz,x,0,2,4),np.zeros((2,2),int)))
        OBS['closed_single_leaf_algebra']=dict(
            general_proof='relocate second code-preserving operator to a disjoint leaf by permutation',
            all_logical_single_leaf_code_endomorphisms_commute=True,
            physical_and_logical_observable_calibrations_may_differ=True,
            full_symmetric_code_local_solutions=ranks,
            proper_GHZ_subcode_has_nontrivial_single_leaf_logical_Z=True,
            nonAbelian_complete_local_write_algebra_impossible_under_this_contract=True,
            requirement_is_static_single_leaf_and_no_code_exit=True)

    def test_05_leakage_and_noncommuting_compressions(self):
        c,w,_=leaf_columns(2)
        v=c/np.sqrt(w)[None,:]
        x=np.array([[0,1],[1,0]],dtype=np.int64);z=np.diag([1,-1])
        ax=local(v,x,0,2,4);bz=local(v,z,1,2,4)
        a=v.T@ax;b=v.T@bz
        ea=ax-v@a;eb=bz-v@b
        comm=a@b-b@a
        self.assertLess(np.max(abs(comm-(-ea.T@eb+eb.T@ea))),3e-15)
        norm=np.linalg.norm(comm,2);leaka=np.linalg.norm(ea,2);leakb=np.linalg.norm(eb,2)
        self.assertAlmostEqual(norm,.5,places=13)
        self.assertLessEqual(norm,2*leaka*leakb+3e-15)
        # Unnormalized rational-basis certificate, with multiset order 0000,
        # 0001,0011,0111,1111: [compressed X,compressed Z]_(1,0)=1/8.
        grama=c.T@local(c,x,0,2,4);gramb=c.T@local(c,z,1,2,4)
        aa=[[F(int(grama[i,j]),int(w[i])) for j in range(5)] for i in range(5)]
        bb=[[F(int(gramb[i,j]),int(w[i])) for j in range(5)] for i in range(5)]
        entry=sum(aa[1][k]*bb[k][0]-bb[1][k]*aa[k][0] for k in range(5))
        self.assertEqual(entry,F(1,8))
        OBS['leakage_boundary']=dict(
            Hermitian_contraction_identity='[a,b]=-L_A^* L_B+L_B^* L_A after disjoint relocation',
            bound='norm([a,b]) <= 2 epsilon_A epsilon_B',
            full_symmetric_qubit_code_compressed_commutator_norm=short(norm),
            rational_nonzero_commutator_entry=str(entry),
            leakage_norms=[short(leaka),short(leakb)],
            naive_claim_all_single_leaf_compressions_commute_is_false=True,
            partial_or_disturbing_readout_not_excluded=True)

    def test_06_collective_and_dynamic_scope_controls(self):
        c=np.zeros((16,2),dtype=np.int64);c[0,0]=c[15,1]=1
        x=np.array([[0,1],[1,0]],dtype=np.int64);z=np.diag([1,-1])
        current=c.copy();leaks=[]
        for site in range(4):
            current=local(current,x,site,2,4)
            leakage=current-c@(c.T@current)
            leaks.append(int(np.trace(leakage.T@leakage)//2))
        self.assertEqual(leaks,[1,1,1,0])
        self.assertTrue(np.array_equal(current,c@x))
        logical_z=c.T@local(c,z,0,2,4);logical_x=c.T@current
        self.assertTrue(np.array_equal(logical_x@logical_z-logical_z@logical_x,x@z-z@x))
        # Internal nodes 1 and 2 are not leaf permuted, and keep full local algebras
        # on the complete 135-dimensional code, independently of the leaf no-go.
        full,w=full_columns(3)
        x3=np.array([[0,1,0],[1,0,0],[0,0,0]],dtype=np.int64)
        z3=np.diag([1,-1,0])
        expected=local(local(full,z3,1,3,6),x3,1,3,6)-local(local(full,x3,1,3,6),z3,1,3,6)
        self.assertGreater(int(np.sum(abs(expected))),0)
        for p in itertools.permutations(range(4)):
            self.assertTrue(np.array_equal(expected[image(3,6,LEAVES,p)],expected))
        OBS['scope_controls']=dict(GHZ_logical_X_uses_all_four_leaves=True,
            sequential_local_X_outside_code_probabilities=leaks,
            noncommuting_logical_X_and_Z_available_with_relaxed_support=True,
            these_gates_are_extra_permissions_not_derived_from_H=True,
            code_exit_not_unknown_information_loss=True,
            internal_unpermuted_nodes_keep_nonAbelian_local_algebras=True,
            no_claim_to_exclude_483_full_dynamical_reference_protection=True,
            no_requirement_that_geometry_summary_predict_all_future=True)

def run():
    OBS.clear();stream=io.StringIO()
    result=unittest.TextTestRunner(stream=stream).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():raise RuntimeError(stream.getvalue())
    output=dict(round=484,scientific_baseline_round=475,additional_frozen_baseline_round=471,
        tests_run=result.testsRun,failures=len(result.failures),errors=len(result.errors),
        python=platform.python_version(),numpy=np.__version__,observations=OBS,
        scope=dict(capacity_increase_does_not_by_itself_restore_independent_leaf_access=True,
            single_leaf_closed_code_condition_is_extra_not_cognitive_axiom=True,
            general_dynamic_private_control_remains_open=True,
            general_mixed_reference_safe_encoding_not_classified=True,
            no_spatial_dimension_or_general_cognitive_no_go=True,
            full_GR_goal_completed=False,phase_closure_triggered=False))
    assert json.loads(json.dumps(output))==output
    return output

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dry-run',action='store_true');parser.add_argument('--check',action='store_true')
    args=parser.parse_args();result=run()
    if args.check:assert result==json.loads(TARGET.read_text(encoding='utf-8'))
    elif not args.dry_run:
        with TARGET.open('x',encoding='utf-8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False,indent=2))
