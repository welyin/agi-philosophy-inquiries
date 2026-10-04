"""Round 471: coherent graph dynamics under recursive qubit block interfaces.

Baseline 470. This audits the scalar energy of coherently controlled raw
exchange, not a spatial dimension or autonomous preparation mechanism.
"""
import argparse
from fractions import Fraction as Q
from functools import lru_cache
import io
import itertools
import json
import math
from pathlib import Path
import platform
import unittest
import numpy as np

import recursive_exchange_interface_audit as lift
import branching_tree_distance_audit as graph
import short_time_tree_relation_readout as readout

TARGET=Path(__file__).with_name('dynamic_block_interface_inheritance_audit_results.json')
OBS={}


def short(x):
    return float(f'{float(x):.12g}')


def constant_profile(m, sizes):
    internal=list(sizes[:m])
    leaves=list(sizes[m:])
    assert len(leaves)==m+2 and all(x>0 and x%2 for x in sizes)
    if m==1 or len(set(internal))==1:
        return True
    if len(set(leaves))>1:
        return False
    return m==2 or sum(x!=leaves[0] for x in internal)<=1


def binary_trees(m):
    n=2*m+2
    words=set(itertools.permutations(tuple(i for i in range(m) for _ in range(2))))
    return sorted({graph.prufer_tree(n,word) for word in words},key=lambda t:tuple(sorted(t)))


def edge_sum(tree,sizes):
    return sum(sizes[a]*sizes[b] for a,b in tree)


def effective_data(tree,n):
    d=2**n
    return sum(np.eye(d,dtype=np.int64)[lift.permutation(n,a,b)] for a,b in tree)


def embed(w,state):
    # state indexed exposed G, graph, private M, reference R.
    d,g,m,r=state.shape
    return (w@state.transpose(0,2,1,3).reshape(d*m,g*r)).reshape(w.shape[0],g,r)


@lru_cache(None)
def compatible_system():
    trees,h,_,_,f,_,_,_,_,za=readout.system()
    sizes=(3,1,1,1,1,1)
    w,ms=lift.grouped_embedding(sizes)
    m=math.prod(ms)
    private=np.diag([-1.,1.])
    edges=[]
    for tree in trees:
        raw,_=lift.lifted_edges(sizes,[(a,b,1) for a,b in tree])
        edges.append(raw)
    zraw=np.array([1 if (x>>5).bit_count()<=1 else -1 for x in range(256)])
    return trees,h,f,sizes,w,m,private,edges,zraw,za


def raw_action(state,alpha=0.,beta=1.,internal=True):
    trees,_,f,sizes,_,_,_,edges,_,_=compatible_system()
    out=np.einsum('gh,xhr->xgr',f,state)
    for g,contact in enumerate(edges):
        out[:,g,:]+=alpha*len(contact)*state[:,g,:]
        out[:,g,:]+=beta*lift.action(sum(sizes),contact,state[:,g,:])
    if internal:
        out+=state[lift.permutation(sum(sizes),0,1)]
    return out


def taylor_action(action,state,t,terms=28):
    out=state.astype(complex).copy()
    term=out.copy()
    for k in range(1,terms+1):
        term=(-1j*t/k)*action(term)
        out+=term
    return out


def private_reference(state):
    # Pure joint state with axes G,graph,M,R.
    _,_,m,r=state.shape
    z=state.reshape(-1,m*r)
    return z.T@z.conj()


class Audit(unittest.TestCase):
    def close(self,a,b,tol=6e-11):
        self.assertLess(float(np.linalg.norm(a-b)),tol)

    def test_01_local_primitive_and_genuine_raw_port(self):
        w,ms=lift.grouped_embedding((3,3))
        sw=np.eye(4)[lift.permutation(2,0,1)]
        raw_edges=[(a,b,1) for a in range(3) for b in range(3,6)]
        k=lift.action(6,raw_edges,np.eye(64))
        self.close(k@w,w@np.kron(sw+4*np.eye(4),np.eye(4)))
        primitive=[]
        for alpha,beta in ((Q(0),Q(1)),(Q(-1,2),Q(1)),(Q(2),Q(3))):
            lhs=(float(alpha)*9*np.eye(64)+float(beta)*k)@w
            scalar=alpha*9+beta*4
            rhs=w@np.kron(float(beta)*sw+float(scalar)*np.eye(4),np.eye(4))
            self.close(lhs,rhs)
            primitive.append(dict(alpha=str(alpha),beta=str(beta),
                                  collective_scalar=str(scalar),gamma=str(alpha+beta/2)))
        v=lift.half_encoding(3)
        sign=np.diag([1 if x.bit_count()<=1 else -1 for x in range(8)])
        self.close(sign@v,v@np.kron(np.diag([1.,-1.]),np.eye(2)))
        self.assertEqual(set(np.diag(sign)),{-1,1})
        self.close(sign@sign,np.eye(8))
        for outcome in (-1,1):
            effect=(np.eye(8)+outcome*sign)/2
            self.close(effect@effect,effect)
        OBS['local_interface']=dict(
            raw_blocks=[3,3],full_private_dimension=math.prod(ms),
            collective_identity='sum raw S = S_G +(n_i*n_j-1)/2 * I',
            primitive_checks=primitive,
            raw_binary_port='sgn(J_block^z), eigenvalues exactly +/-1 on full raw space',
            collective_2Jz_not_used_as_full_space_binary_effect=True)

    def test_02_dynamic_generator_and_unknown_private_reference(self):
        trees,h,f,sizes,w,m,hm,edges,zraw,za=compatible_system()
        scalar=[]
        for g,tree in enumerate(trees):
            k=lift.action(8,edges[g],w)
            c=Q(edge_sum(tree,sizes)-5,2)
            self.close(k,w@np.kron(effective_data(tree,6)+float(c)*np.eye(64),np.eye(m)))
            scalar.append(str(c))
        self.assertEqual(set(scalar),{'1'})
        # One whole raw graph action on an arbitrary complex G-graph-M-R ket.
        index=np.arange(64*6*m*3).reshape(64,6,m,3)
        state=((index*7%23)-11)+1j*((index*11%19)-9)
        state=state/np.linalg.norm(state)
        def logical_action(x):
            y=(h@x.reshape(384,-1)).reshape(x.shape)+x
            return y+np.einsum('mn,dgnr->dgmr',hm,x)
        raw=embed(w,state)
        self.close(raw_action(raw),embed(w,logical_action(state)))
        t=Q(1,8)
        raw_out=taylor_action(raw_action,raw,float(t))
        logical_out=taylor_action(logical_action,state,float(t))
        self.close(raw_out,embed(w,logical_out))
        um=readout.evolution(hm,float(t))
        umr=np.kron(um,np.eye(3))
        self.close(private_reference(logical_out),
                   umr@private_reference(state)@umr.conj().T)
        self.close(zraw[:,None,None]*raw_out,
                   embed(w,za[:,None,None,None]*logical_out))
        # ||H_raw|| <= 4 + 7 + 1=12; omitted state-vector Taylor tail.
        x=Q(12)*t
        tail=Q(3)**math.ceil(x)*x**29/Q(math.factorial(29))
        self.assertLess(tail,Q(1,10**23))
        OBS['full_dynamic_intertwining']=dict(
            block_sizes=list(sizes),graph_dimension=6,raw_dimension=1536,
            complete_code_dimension=768,private_dimension=m,reference_dimension=3,
            graph_potential_values=scalar,time=str(t),
            arbitrary_coherent_G_graph_private_reference_state=True,
            private_reference_evolves_only_by_own_hM=True,
            full_raw_action_not_replaced_by_compressed_action=True,
            Taylor_state_tail_bound=str(tail),roundoff_tolerance='6e-11',
            graph_register_and_F_retained_as_input_hardware=True)

    def test_03_necessary_sufficient_full_binary_tree_classification(self):
        profiles={
            1:[(3,1,5,7)],
            2:[(3,3,1,3,5,7),(3,5,1,1,1,1),(3,5,1,1,1,3)],
            3:[(3,3,3,1,3,5,7,9),(5,1,1,1,1,1,1,1),
               (3,3,1,1,1,1,1,1),(3,5,7,1,1,1,1,1),
               (3,5,1,1,1,1,1,3)]}
        rows=[]
        delta_checks=0
        for m,ps in profiles.items():
            ts=binary_trees(m)
            self.assertEqual(len(ts),math.factorial(2*m)//(2**m))
            for sizes in ps:
                vals={edge_sum(t,sizes) for t in ts}
                self.assertEqual(len(vals)==1,constant_profile(m,sizes))
                rows.append(dict(internal_count=m,trees=len(ts),sizes=list(sizes),
                                 constant=len(vals)==1,distinct_edge_sums=sorted(vals)))
                for tree in ts:
                    adj=graph.adjacency(tree,2*m+2)
                    for b,c in tree:
                        for a in adj[b]-{c}:
                            for d in adj[c]-{b}:
                                new=tree.symmetric_difference({
                                    graph.edge(a,b),graph.edge(c,d),
                                    graph.edge(a,c),graph.edge(b,d)})
                                self.assertEqual(edge_sum(new,sizes)-edge_sum(tree,sizes),
                                                 (sizes[a]-sizes[d])*(sizes[c]-sizes[b]))
                                delta_checks+=1
        OBS['classification']=dict(
            finite_certificate_profiles=rows,exact_NNI_delta_checks=delta_checks,
            all_m_ge3='all internal sizes equal OR all leaves and all but at most one internal size equal',
            m2='internal sizes equal OR all leaf sizes equal',
            m1='single labelled star: any sizes',
            gamma_alpha_plus_beta_over2_zero_removes_constraint=True,
            default_alpha_zero_and_J_zero_special_case=True,
            analytic_general_proof_not_extrapolated_from_profiles=True)

    def test_04_scalar_primitive_choice_changes_coherent_graph_energy(self):
        trees,h,_,_,f,_,_,_,_,_=readout.system()
        sizes=(3,3,1,1,1,1)
        w=np.array([edge_sum(t,sizes)-5 for t in trees],dtype=np.int64)
        self.assertEqual(w.tolist(),[8,8,12,8,12,12])
        self.assertGreater(np.linalg.norm(np.diag(w)@f-f@np.diag(w)),0)
        gamma=Q(-1,2)+Q(1,2)
        self.assertEqual(gamma,0)
        self.assertEqual([gamma*int(x) for x in w],[Q(0)]*6)
        # Independent direct raw centered action in the heterogeneous positive
        # example; target is centered logical primitive too.
        _,hh,_,_,ww,mm,_,_,_,_=compatible_system()
        index=np.arange(64*6*mm).reshape(64,6,mm,1)
        x=((index%13)-6)+1j*((index%17)-8)
        x=x/np.linalg.norm(x)
        target=((hh-2.5*np.eye(384))@x.reshape(384,-1)).reshape(x.shape)
        self.close(raw_action(embed(ww,x),alpha=-.5,internal=False),embed(ww,target))
        OBS['scalar_energy_audit']=dict(
            nonconstant_example_sizes=list(sizes),
            sum_edges_ninj_minus1=w.tolist(),
            delta_potential_same_logical_primitive='(alpha+beta/2)*sum_edges(n_i*n_j-1)',
            centered_alpha='-beta/2',centered_delta_potential_exactly_zero=True,
            fixed_graph_phase_does_not_fix_coherently_controlled_phase=True,
            primitive_scalar_energy_is_unselected_model_input=True,
            no_conclusion_all_subjects_must_have_equal_size=True)

    def test_05_actual_470_port_readout_distinguishes_size_potential(self):
        trees,h,_,_,_,b,vp,vm,_,_=readout.system()
        n01=np.array([int((0,1) in t) for t in trees],dtype=np.int64)
        self.assertEqual(n01.tolist(),[0,0,1,0,1,1])
        tilted=h+2*np.kron(np.eye(64,dtype=np.int64),np.diag(n01))
        w=np.array([1,0,1,0,0,0],dtype=np.int64)
        plus,minus=vp@w,vm@w
        self.assertEqual(int(plus@plus),8)
        a=b.copy();c=b.copy();derivatives=[]
        for n in range(6):
            z=(plus@(c-a)@plus-minus@(c-a)@minus)
            if n%2:
                self.assertEqual(int(z),0)
                derivative=Q(0)
            else:
                derivative=Q(int(z)*((-1)**(n//2)),8)
            derivatives.append(str(derivative))
            a=h@a-a@h
            c=tilted@c-c@tilted
        self.assertEqual(derivatives,['0','0','0','0','12','0'])
        t=Q(1,2000)
        lead=t**4/2
        xx=22*t
        tail=4*xx**6/Q(math.factorial(6))/(1-xx**2/56)
        lower=lead-tail
        self.assertGreater(lower,0)
        values=[]
        sample_time=.1
        for hh in (h,tilted):
            u=readout.evolution(hh,sample_time)
            pp=u@(plus/math.sqrt(8));pm=u@(minus/math.sqrt(8))
            values.append(float((np.vdot(pp,b@pp)-np.vdot(pm,b@pm)).real))
        self.assertGreater(abs(values[1]-values[0]),1e-7)
        OBS['actual_port_consequence']=dict(
            effective_tilt_after_common_phase='2*n_01',
            coherent_graph_state='(T0+T2)/sqrt(2)',
            all_G_probe_states='round 470 eta_plus and eta_minus',
            output='g=mean(Z_G0|eta_plus)-mean(Z_G0|eta_minus)',
            exact_derivative_difference_orders_0_to_5=derivatives,
            analytic_leading_difference='t^4/2',
            certificate_time=str(t),leading=str(lead),tail_bound=str(tail),
            positive_readout_difference_lower_bound=str(lower),
            sample_time=sample_time,sample_base_g=short(values[0]),
            sample_lifted_g=short(values[1]),sample_difference=short(values[1]-values[0]),
            not_merely_spectral_difference=True,
            original_470_second_coefficient_still_equal=True,
            no_claim_all_nonconstant_potentials_detected_by_this_one_protocol=True)

    def test_06_actual_raw_470_readout_and_retained_private_information(self):
        trees,h,f,sizes,w,m,hm,_,zraw,za=compatible_system()
        ix=np.arange(6*m*3).reshape(6,m,3)
        source=((ix*7%17)-8)+1j*((ix*5%13)-6)
        source=source/np.linalg.norm(source)
        t=1/8
        u=readout.evolution(h,t)
        um=readout.evolution(hm,t)
        results=[]
        for sign in (1,-1):
            probe=readout.probe_vector(6,0,5,sign)/2
            psi=probe[:,None,None,None]*source[None,:,:,:]
            raw=embed(w,psi)
            actual=taylor_action(raw_action,raw,t)
            expected=(u@psi.reshape(384,-1)).reshape(psi.shape)
            expected=np.einsum('mn,dgnr->dgmr',um,expected)*np.exp(-1j*t)
            self.close(actual,embed(w,expected))
            actual_z=float(np.sum(actual.conj()*zraw[:,None,None]*actual).real)
            expected_z=float(np.sum(expected.conj()*za[:,None,None,None]*expected).real)
            self.assertAlmostEqual(actual_z,expected_z,places=11)
            initial_mr=private_reference(psi)
            updated=private_reference(expected)
            umr=np.kron(um,np.eye(3))
            self.close(updated,umr@initial_mr@umr.conj().T)
            results.append(dict(sign=sign,raw_mean=short(actual_z),
                                inherited_G_mean=short(expected_z)))
        OBS['actual_recursive_relation_readout']=dict(
            graph_private_reference_dimensions=[6,m,3],time='1/8',readouts=results,
            all_G_probes_prepared_independently_of_graph_private_reference=True,
            unknown_private_graph_reference_joint_state_retained=True,
            private_internal_motion_accounted_in_co_rotating_reference=True,
            raw_Z_is_legal_collective_sign_effect=True,
            arbitrary_old_G_state_not_preserved_by_probe_preparation=True,
            original_graph_carrier_F_clock_preparation_and_reference_are_inputs=True,
            primitive_contacts_per_present_edge='n_i*n_j',
            per_raw_unit_exchange_only_coefficient_budget='sum_{j adjacent to i} n_j*abs(beta)',
            per_raw_unit_full_controlled_primitive_sufficient_budget='sum_{j adjacent to i} n_j*(abs(alpha)+abs(beta))',
            full_spatial_geometry_or_complete_cognitive_subject_not_derived=True)


def run():
    OBS.clear()
    stream=io.StringIO()
    result=unittest.TextTestRunner(stream=stream).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise RuntimeError(stream.getvalue())
    return dict(round=471,baseline_round=470,tests_run=result.testsRun,
                failures=len(result.failures),errors=len(result.errors),
                python=platform.python_version(),numpy=np.__version__,
                observations=OBS,scope=dict(
                    exact_dynamic_block_intertwining=True,
                    arbitrary_graph_private_reference_states_in_selected_codes=True,
                    all_fixed_degree_binary_tree_size_classification=True,
                    controlled_scalar_energy_choice_explicit=True,
                    actual_data_port_consequence_and_inheritance=True,
                    graph_register_and_F_recursively_generated=False,
                    arbitrary_raw_state_or_G_probe_preparation_free=False,
                    full_cognitive_subject_derived=False,
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

