"""Round 470: short-time two-step tree relation readout from one data port.

Baseline 469. Same continuous graph-flip plus graph-controlled SWAP H.
The two preparation settings, all data probes, common X/Z reference, time and
independent repeatable graph-state source are explicit operation inputs.
"""
import argparse
from fractions import Fraction as Q
from functools import lru_cache
import io
import json
import math
from pathlib import Path
import platform
import unittest
import numpy as np
import branching_tree_distance_audit as old
import unit_tree_endpoint_metric_audit as endpoint

TARGET=Path(__file__).with_name('short_time_tree_relation_readout_results.json')
OBS={}


def short(x):
    return float(f'{float(x):.12g}')


def trace_norm(x):
    assert np.linalg.norm(x-x.conj().T)<1e-10
    return float(np.abs(np.linalg.eigvalsh((x+x.conj().T)/2)).sum())


def probe_vector(n,a,b,sign):
    v=np.zeros(2**n,dtype=np.int64)
    for abit in (0,1):
        for bbit in (0,1):
            index=(abit<<(n-1-a))|(bbit<<(n-1-b))
            v[index]=sign**bbit
    return v  # Divide by two for its normalized ket.


def commutator(h,b):
    return h@b-b@h


@lru_cache(None)
def system():
    trees,_,h,f,_=old.six_vertex_sector()
    hd=h-np.kron(np.eye(64,dtype=np.int64),f)
    hf=np.kron(np.eye(64,dtype=np.int64),f)
    za=np.array([1 if (x>>5)==0 else -1 for x in range(64)])
    b=np.diag(np.repeat(za,6))
    vp=np.kron(probe_vector(6,0,5,1)[:,None],np.eye(6,dtype=np.int64))
    vm=np.kron(probe_vector(6,0,5,-1)[:,None],np.eye(6,dtype=np.int64))
    q=[]
    for tree in trees:
        adj=endpoint.adjacent(tree,6)
        q.append(len(adj[0]&adj[5]))
    return trees,h,hd,hf,f,b,vp,vm,np.diag(q),za


def compress(v,x):
    return v.T@x@v  # Normalization factor 1/4 is kept explicit.


def evolution(h,t):
    e,v=np.linalg.eigh(h)
    return (v*np.exp(-1j*t*e))@v.conj().T


@lru_cache(None)
def reference_input():
    raw=np.array([[1,1j,2],[2,-1j,1],[1+1j,2,0],
                  [-1,1,2j],[2,1j,1-1j],[1,-2,1j]],dtype=complex)
    norm=int(np.vdot(raw,raw).real)
    return raw/math.sqrt(norm),raw,norm


def signed_reference(out,za):
    weights=np.repeat(za,6)
    return np.einsum('ir,i,is->rs',out,weights,out.conj())


def record_blocks(out,za):
    x=out.reshape(64,18)
    return [x[za==sign].T@x[za==sign].conj() for sign in (1,-1)]


class Audit(unittest.TestCase):
    def close(self,a,b,tol=2e-11):
        self.assertLess(np.linalg.norm(a-b),tol)

    def test_01_exact_local_sign_and_degree_correction(self):
        n,a,b=3,0,1
        vp,vm=(probe_vector(n,a,b,sign) for sign in (1,-1))
        z=np.diag([1 if x<4 else -1 for x in range(8)])
        swaps=[old.core.swap(3,0,1).real.astype(np.int64),
               old.core.swap(3,0,2).real.astype(np.int64),
               old.core.swap(3,1,2).real.astype(np.int64)]
        cases=[(1,0,0),(0,1,0),(0,0,1),(1,1,0),(1,0,1),(0,1,1),(1,1,1)]
        rows=[]
        for x,y,w in cases:
            h=x*swaps[0]+y*swaps[1]+w*swaps[2]
            c2=commutator(h,commutator(h,z))
            plus=Q(-int(vp@c2@vp),4)
            minus=Q(-int(vm@c2@vm),4)
            difference=plus-minus
            self.assertEqual(difference,2*(y*w+x*(y-w)))
            rows.append(dict(ab=x,av=y,bv=w,
                plus_second=str(plus),minus_second=str(minus),difference=str(difference)))
        path=[r for r in rows if (r['ab'],r['av'],r['bv'])==(0,1,1)][0]
        self.assertEqual((path['plus_second'],path['minus_second']),('3','1'))
        OBS['exact_local_coefficient'] = dict(
            sign_convention='g = mean(Z_a | +X_b) - mean(Z_a | -X_b)',
            rational_cases=rows,
            general_graph_identity='g_second(0)=2*J^2*E[sum_v n_av*n_bv+n_ab*(D_a-D_b)]',
            fixed_leaves_N_ge_3_remove_direct_edge_term=True,
            fixed_equal_degrees_also_sufficient_for_correction_to_vanish=True)

    def test_02_full_quantum_graph_operator_certificate(self):
        trees,h,hd,hf,f,b,vp,vm,q,za=system()
        first=commutator(hd,b)
        second=commutator(hd,first)
        mixed=commutator(hf,first)
        self.assertTrue(np.array_equal(compress(vp,b)-compress(vm,b),np.zeros((6,6))))
        self.assertTrue(np.array_equal(compress(vp,first),np.zeros((6,6))))
        self.assertTrue(np.array_equal(compress(vm,first),np.zeros((6,6))))
        self.assertTrue(np.array_equal(-(compress(vp,second)-compress(vm,second)),8*q))
        self.assertTrue(np.array_equal(commutator(hf,b),np.zeros_like(b)))
        self.assertTrue(np.array_equal(compress(vp,mixed),np.zeros((6,6))))
        self.assertTrue(np.array_equal(compress(vm,mixed),np.zeros((6,6))))
        v=hd-5*np.eye(384,dtype=np.int64)
        for prep in (vp,vm):
            self.assertTrue(np.array_equal(compress(prep,v),-4*np.eye(6,dtype=np.int64)))
        indices=[6*x+g for x in range(64) if x.bit_count()<=2 for g in range(6)]
        complement=sorted(set(range(384))-set(indices))
        self.assertFalse(np.any(v[np.ix_(complement,indices)]))
        row_bound=int(np.abs(v[np.ix_(indices,indices)]).sum(axis=1).max())
        self.assertLessEqual(row_bound,12)
        OBS['operator_certificate'] = dict(
            full_dimension=384, graph_dimension=6,
            distance_two_projector_diagonal=q.diagonal().tolist(),
            exact_zero_orders=[0,1], exact_second_derivative='2*J^2*Q_ab',
            mixed_kappa_J_coefficient_zero_in_each_preparation=True,
            arbitrary_real_J_and_kappa_covered_symbolically=True,
            reference_extension_is_operator_identity=True,
            probe_mean_interaction='Tr_D(eta_pm * V) = -J*I_G',
            at_most_two_excitation_sector_dimension=len(indices),
            exact_interaction_absolute_row_bound_at_J1=row_bound,
            universal_interaction_bound='v=4*C*abs(J)')

    def test_03_unknown_graph_reference_finite_time_readout(self):
        _,h,_,_,_,_,vp,vm,q,za=system()
        psi,raw,norm=reference_input()
        t=Q(1,64)
        u=evolution(h,float(t))
        plus,minus=(u@(v@psi/2) for v in (vp,vm))
        delta=signed_reference(plus,za)-signed_reference(minus,za)
        qr=psi.T@(q@psi.conj())
        rho_r=psi.T@psi.conj()
        self.assertGreaterEqual(np.linalg.eigvalsh(qr).min(),-1e-13)
        self.assertGreaterEqual(np.linalg.eigvalsh(rho_r-qr).min(),-1e-13)
        residual=trace_norm(delta-float(t*t)*qr)
        # This finite example has ||H||<=9. Its factorial Taylor remainder,
        # independent of any vacuum assumption on the graph/reference,
        # is tighter than the deliberately loose all-size K from round 466.
        x=18*t
        example_bound=2*x**3/Q(6)/(1-x/4)
        self.assertLess(residual,float(example_bound))
        q_exact=Q(int(np.vdot(raw,q@raw).real),norm)
        g=float(np.trace(delta).real)
        graph=psi@psi.conj().T
        self.assertGreater(np.linalg.norm(graph-np.diag(np.diag(graph))),.1)
        OBS['reference_valued_readout'] = dict(
            matrix_time=str(t), unknown_graph_reference_dimensions=[6,3],
            graph_population_distance_two_exact=str(q_exact),
            graph_off_diagonal_Frobenius_norm=short(np.linalg.norm(graph-np.diag(np.diag(graph)))),
            signed_reference_moment_trace=short(g),
            raw_finite_time_q_estimate=short(g/float(t*t)),
            signed_reference_operator_error=short(residual),
            finite_example_remainder_bound=str(example_bound),
            finite_example_bound_uses_H_norm_9_not_all_size_K=True,
            QR_positive_and_bounded_by_reference_marginal=True,
            Delta_A_R_is_a_signed_operator_on_R_not_a_joint_state=True,
            graph_coherence_and_reference_not_measured=True)

    def test_04_actual_disturbance_and_classical_record(self):
        _,h,_,_,f,_,vp,vm,_,za=system()
        psi,_,_=reference_input()
        t=Q(1,64)
        u=evolution(h,float(t))
        free=evolution(f,float(t))@psi
        baseline=np.outer(free.reshape(-1),free.reshape(-1).conj())
        rho_r=psi.T@psi.conj()
        v=12
        marginal_bound=.5*(math.expm1(2*v*float(t))-2*v*float(t))
        joint_record_bound=min(1,v*float(t))
        rows=[]
        for label,prep in [('+X_b',vp),('-X_b',vm)]:
            out=u@(prep@psi/2)
            blocks=record_blocks(out,za)
            marginal=sum(blocks)
            marginal_error=.5*trace_norm(marginal-baseline)
            record_error=.5*sum(trace_norm(block-.5*baseline) for block in blocks)
            self.assertGreater(marginal_error,1e-8)
            self.assertLessEqual(marginal_error,marginal_bound)
            self.assertLessEqual(record_error,joint_record_bound)
            self.close(marginal.reshape(6,3,6,3).trace(axis1=0,axis2=2),rho_r)
            rows.append(dict(preparation=label,
                graph_reference_trace_distance=short(marginal_error),
                record_graph_reference_trace_distance=short(record_error)))
        OBS['disturbance'] = dict(
            matrix_time=str(t), comparison='unprobed vacuum data: graph evolution exp(-i*kappa*F*t)',
            observed=rows, reference_marginal_exactly_preserved=True,
            nonselective_graph_reference_bound='min(1,(exp(2*v*t)-1-2*v*t)/2)',
            nonselective_bound_at_example=short(marginal_bound),
            coarse_classical_record_joint_bound='min(1,v*t)',
            record_joint_bound_at_example=short(joint_record_bound),
            coarse_bound_not_claimed_optimal=True,
            unknown_graph_reference_not_assumed_undisturbed=True,
            conditioning_on_outcome_not_claimed_nondisturbing=True)

    def test_05_all_size_bias_and_statistical_resource_certificate(self):
        c,j,kappa=3,1,1
        lam=c*abs(j)+2*c*(c-1)**2*abs(kappa)
        k=96*lam*c**6
        t=Q(1,100*k**3)
        remainder=2*(k*t)**3/(1-k*t)
        bias=remainder/(j*j*t*t)
        self.assertLess(k*t,Q(1,2))
        self.assertLess(bias,Q(1,20))
        eta=Q(1,20)
        # exp(6)>400 follows already from a finite positive Taylor sum.
        self.assertGreater(sum((Q(6)**m/Q(math.factorial(m)) for m in range(16)),Q(0)),400)
        per_branch=19200*(100*k**3)**4
        self.assertEqual(Q(per_branch),48/(eta**2*t**4))
        self.assertLess(bias+eta,Q(1,10))
        OBS['finite_resource_certificate'] = dict(
            C=c,J=j,kappa=kappa,K=k,
            time_exact=str(t),time_approx=short(float(t)),
            reference_remainder_exact=str(remainder),
            q_bias_exact=str(bias),q_bias_approx=short(float(bias)),
            statistical_q_tolerance=str(eta),failure_probability_at_most='1/100',
            per_preparation_independent_copies=str(per_branch),
            per_preparation_copies_scientific=f'{per_branch:.6e}',
            exp6_greater_than_400_exactly_certified=True,
            total_q_error_less_than='1/10',
            general_copies_formula='M >= 8*log(4/delta)/(eta^2*J^4*t^4) per preparation',
            all_N_data_probes_prepared_for_every_copy=True,
            copies_are_an_extra_source_contract_not_cloned=True,
            no_actual_huge_sampling_run_performed=True,
            inherited_K_is_deliberately_loose_not_practical_protocol_optimization=True)

    def test_06_distance_scope_and_same_department_counterexample(self):
        trees,_,_,_,_,_,_,_,q,_=system()
        ds=[int(endpoint.distances(g,6,6)[0,5]) for g in trees]
        self.assertEqual(ds,[3-int(x) for x in q.diagonal()])
        longer=[]
        for shore in ({1},{3}):
            edges,total,degrees,_,_=endpoint.tree_data(4,shore,4)
            matrix=endpoint.distances(edges,total,4)
            adj=endpoint.adjacent(edges,total)
            longer.append(dict(N=total,degree_vector=degrees,distance=int(matrix[0,3]),
                               common_neighbors=len(adj[0]&adj[3])))
        self.assertEqual(longer[0]['N'],longer[1]['N'])
        self.assertEqual(longer[0]['degree_vector'],longer[1]['degree_vector'])
        self.assertEqual([row['common_neighbors'] for row in longer],[0,0])
        self.assertEqual([row['distance'] for row in longer],[3,5])
        OBS['scope_boundary'] = dict(
            six_tree_only_mean_distance='E[D_ab]=3-q',
            same_unit_edge_fixed_degree_department_examples=longer,
            same_q_does_not_determine_general_mean_distance=True,
            not_claiming_all_time_data_channels_of_examples_equal=True,
            original_data_unknown_states_not_preserved_by_probe_preparation=True,
            operational_three_dimensional_geometry_not_derived=True)


def run():
    OBS.clear()
    stream=io.StringIO()
    result=unittest.TextTestRunner(stream=stream).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise RuntimeError(stream.getvalue())
    return dict(round=470,baseline_round=469,tests_run=result.testsRun,
        failures=len(result.failures),errors=len(result.errors),
        python=platform.python_version(),numpy=np.__version__,observations=OBS,
        scope=dict(same_continuous_H_graph_F_plus_controlled_SWAP=True,
            all_size_fixed_leaf_distance_two_population=True,
            one_data_port_Z_measurement_only=True,
            graph_reference_arbitrary_under_independent_probe_preparation=True,
            complete_reference_signed_moment_bound=True,
            nonzero_graph_backaction_explicitly_allowed=True,
            all_probe_preparations_shared_X_Z_reference_and_clock_are_inputs=True,
            repeated_identical_graph_source_is_an_input=True,
            no_graph_measurement_or_distance_oracle_used=True,
            autonomous_device_derived=False,
            arbitrary_original_data_preserved=False,
            general_full_path_distance_reconstructed=False,
            physical_spatial_dimension_derived=False,
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
