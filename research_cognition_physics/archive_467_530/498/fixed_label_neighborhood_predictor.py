"""Round 498: a fixed-label coherent neighborhood predictor.

Baseline 496; independent of 497. The actual network always evolves with H.
H_A is an explicit CPTP mathematical predictor, not an external switch in H.
"""
import argparse
from fractions import Fraction as Q
from functools import lru_cache
import hashlib, io, itertools, json, math, platform, unittest
from pathlib import Path
import numpy as np
import branching_tree_distance_audit as tree
import unknown_background_propagation_audit as prior
import quantum_neighborhood_signal_audit as neighborhood

ROOT=Path(__file__).resolve().parent
TARGET=ROOT/'fixed_label_neighborhood_predictor_results.json'
OBS={}
K=1889568
# Exact sparse Gaussian-integer phase calculations use at most two H factors.
# At most 15 terms per row, initial entries <=1: intermediate <=15^2 and
# inner sums <=23040*15^2*2 <2^53. No dense 23040-square matrix is built.
assert 23040*15**2*2 < 2**53
# The independent integer intertwining vector has entries <23040.
assert 15*23040 < 2**63


def short(x):return float(f'{float(x):.12g}')


def swap_indices(n,a,b):
    z=np.arange(2**n)
    different=((z>>(n-1-a))&1)!=((z>>(n-1-b))&1)
    return z^(different*((1<<(n-1-a))|(1<<(n-1-b))))


@lru_cache(None)
def model():
    trees=sorted({tree.prufer_tree(8,p) for p in set(itertools.permutations((1,1,2,2,3,3)))},key=lambda g:sorted(g))
    atoms=prior.atoms(trees,8)
    local=frozenset(range(6))
    chosen=[x for x in atoms if x[1]<=local]
    g=frozenset({(0,1),(1,2),(2,5),(2,4),(1,3),(3,6),(3,7)})
    h=frozenset({(0,2),(1,2),(1,5),(2,4),(1,3),(3,6),(3,7)})
    ids=[trees.index(g),trees.index(h)]
    gl=sorted({tuple(e for e in sorted(t) if set(e)<=local) for t in trees})
    ge=sorted({tuple(e for e in sorted(t) if not set(e)<=local) for t in trees})
    il=np.array([gl.index(tuple(e for e in sorted(t) if set(e)<=local)) for t in trees])
    ie=np.array([ge.index(tuple(e for e in sorted(t) if not set(e)<=local)) for t in trees])
    la=[]
    for kind,support,targets,parameters in chosen:
        lt=np.full(len(gl),-1,dtype=np.int64)
        for j,edges in enumerate(gl):
            es=set(edges)
            if kind=='D':
                if parameters in es:lt[j]=j
            else:
                a,b,c,d=parameters
                left={tree.edge(a,b),tree.edge(c,d)}
                right={tree.edge(a,c),tree.edge(b,d)}
                if tree.edge(b,c) in es and ((left<=es and not right&es) or (right<=es and not left&es)):
                    lt[j]=gl.index(tuple(sorted(es.symmetric_difference(left|right))))
        la.append((kind,support,lt,parameters))
    return trees,atoms,chosen,ids,gl,ge,il,ie,la


def apply_h(x,atoms,n):
    y=np.zeros_like(x)
    data=np.arange(2**n)
    for kind,support,targets,parameters in atoms:
        src=np.flatnonzero(targets>=0);dst=targets[src]
        perm=swap_indices(n,*parameters) if kind=='D' else data
        y[perm[:,None],dst[None,:]]+=x[data[:,None],src[None,:]]
    return y


def pack(x):
    *_,gl,ge,il,ie,la=model()
    ans=np.zeros((64,len(gl),4*len(ge))+x.shape[2:],dtype=x.dtype)
    for z in range(256):
        ans[z>>2,il,(z&3)*len(ge)+ie]=x[z]
    return ans


def receiver_full(x):
    # Data0 is the high data bit. All other data, graph and environment traced.
    r=x.shape[-1]
    y=x.reshape(2,128,90,r).transpose(0,3,1,2).reshape(2*r,-1)
    return y@y.conj().T


def receiver_local(x):
    r=x.shape[-1]
    y=x.reshape(2,32,54,36,r).transpose(0,4,1,2,3).reshape(2*r,-1)
    return y@y.conj().T


def taylor_action(x,atoms,n,t,degree=16):
    out=x.astype(complex).copy();term=out.copy()
    for j in range(1,degree+1):
        term=(-1j*float(t)/j)*apply_h(term,atoms,n)
        out+=term
    return out


def tail_action(t,degree=16,norm=15):
    t=abs(Q(t));first=(norm*t)**(degree+1)/math.factorial(degree+1)
    return first/(1-norm*t/(degree+2))


def phase_sources():
    trees,atoms,chosen,ids,*_=model()
    ans=[]
    for sign in (1,-1):
        x=np.zeros((256,90),complex)
        for z,c in zip((0,128,64,192),(1,1,1j,1j)):
            x[z,ids[0]]=c;x[z,ids[1]]=sign*1j*c
        ans.append(x)
    return ans # norm squared exactly8


def active_atoms(g,n):
    for e in sorted(g):yield ('D',frozenset(e),g)
    adjacent=tree.adjacency(g,n)
    for b,c in sorted(g):
        for a in sorted(adjacent[b]-{c}):
            for d in sorted(adjacent[c]-{b}):
                gp=g.symmetric_difference({tree.edge(a,b),tree.edge(c,d),tree.edge(a,c),tree.edge(b,d)})
                yield ('F',frozenset((a,b,c,d)),frozenset(gp))


def two_arms(length):
    assert length>=2
    edges={(0,1),(1,2),(2,3)}
    arms=[list(range(4,4+length)),list(range(4+length,4+2*length))]
    nxt=4+2*length
    for parent,arm in zip((1,2),arms):
        edges.add(tree.edge(parent,arm[0]))
        edges.update(tree.edge(a,b) for a,b in zip(arm[:-1],arm[1:]))
        for j,v in enumerate(arm):
            for _ in range(2 if j==length-1 else 1):
                edges.add(tree.edge(v,nxt));nxt+=1
    other=set(edges).symmetric_difference({(0,1),(2,3),(0,2),(1,3)})
    return frozenset(edges),frozenset(other),nxt


def outside_record_signature(g,n,radius):
    region=neighborhood.ball(g,0,radius,n)
    outside=tuple(v for v in range(n) if v not in region)
    # All absent edge positions are determined by this outside vertex mask.
    occupied=tuple(e for e in sorted(g) if e[0] in outside and e[1] in outside)
    return outside,occupied


class Audit(unittest.TestCase):
    def test_01_fixed_factor_CPTP_split_and_exact_atomic_intertwining(self):
        trees,atoms,chosen,ids,gl,ge,il,ie,la=model()
        self.assertEqual(len(trees),90)
        self.assertEqual(len(set(zip(il,ie))),90)
        checks=0
        for original,local in zip(chosen,la):
            kind,support,targets,parameters=original
            self.assertTrue(support<=set(range(6)))
            lt=local[2]
            for j,dest in enumerate(targets):
                self.assertEqual(dest>=0,lt[il[j]]>=0)
                if dest>=0:
                    self.assertEqual(ie[dest],ie[j])
                    self.assertEqual(il[dest],lt[il[j]])
                    self.assertEqual(targets[dest],j)
                checks+=1
            if kind=='D':
                glob=swap_indices(8,*parameters)
                loc=swap_indices(6,*parameters)
                for z in range(256):
                    self.assertEqual(glob[z]>>2,loc[z>>2]);self.assertEqual(glob[z]&3,z&3)
        # Full basis incidence identity, no coherent neighborhood selection.
        x=np.arange(256*90,dtype=np.int64).reshape(256,90)
        self.assertTrue(np.array_equal(pack(apply_h(x,chosen,8)),apply_h(pack(x),la,6)))
        self.assertEqual(sum(len(tree.tree_flips(g,8)) for g in trees),90*8)
        OBS['fixed_factor_predictor']=dict(
            full_data_dimension=256,tree_dimension=90,full_allowed_dimension=23040,
            original_atoms=len(atoms),retained_atoms=len(chosen),all_tree_graph_targets_checked=checks,
            local_subjects=list(range(6)),excluded_subjects=[6,7],
            compressed_local_graph_patterns=len(gl),compressed_external_graph_patterns=len(ge),
            local_input_factor_qubits=6+math.comb(6,2),
            local_pattern_compression_only_a_numeric_representation=True,
            exact_intertwining_V_HA_equals_hA_tensor_I_V=True,
            HA_retains_full_nij_SWAP_atoms_not_one_excitation_minus_L=True,
            every_chosen_atom_preserves_the_full_tree_sector=True,
            zero_edge_bits_and_physical_factor_identities_retained=True,
            ordinary_fixed_partial_trace_and_unitary_decoder_are_CPTP=True)

    def test_02_initial_envelope_not_invariant_and_low_order_connected_words(self):
        trees,atoms,chosen,ids,*_=model();local=set(range(6))
        balls=[sorted(neighborhood.ball(trees[j],0,2,8)) for j in ids]
        self.assertTrue(all(set(ball)<=local for ball in balls))
        self.assertTrue(any(targets[j]>=0 and targets[j] not in ids for _,_,targets,_ in atoms for j in ids))
        # Longer non-global envelope to check the nontrivial k=1,2 word statement.
        g,h,n=two_arms(8);sources=(g,h);r=6
        envelope=set.union(*(set(neighborhood.ball(t,0,r,n)) for t in sources))
        self.assertLess(len(envelope),n)
        count=outside_connected=0
        for initial in sources:
            stack=[(initial,())]
            while stack:
                current,supports=stack.pop()
                if supports and prior.rooted_atom_family(supports,0):
                    support=frozenset({0}).union(*supports)
                    self.assertLessEqual(len(support),3*len(supports)+1)
                    self.assertEqual(len(prior.induced_components(initial,support)),1)
                    count+=1
                    if not support<=envelope:outside_connected+=1
                if len(supports)<2:
                    for kind,supp,gp in active_atoms(current,n):
                        stack.append((gp,supports+(supp,)))
        self.assertEqual(outside_connected,0);self.assertGreater(count,0)
        OBS['initial_support_envelope']=dict(
            example_root=0,example_radius=2,example_initial_graph_indices=ids,example_balls=balls,
            initial_two_graph_support_not_H_invariant=True,
            separate_word_check_vertices=n,word_check_radius=r,word_check_envelope_vertices=len(envelope),
            nonzero_connected_words_checked=count,maximum_word_order=2,outside_envelope_words=outside_connected,
            general_derivative_match='P0 ad_H^k(B) P0=P0 ad_HA^k(B) P0 when 3*k<=r',
            matching_uses_initial_graph_support_not_future_containment=True)

    def test_03_uniform_analytic_CPTP_prediction_error(self):
        c=3;lam=c+2*c*(c-1)**2
        self.assertEqual(96*lam*c**6,K)
        # Same majorant applies after removing atoms, not from ||H_A||<=||H||.
        m=math.floor(29/3)+1;self.assertEqual(m,10)
        error=2*Q(11,70)**m
        self.assertLess(error,Q(1,50000000))
        summary_error=Q(1,100000000);decoder_error=Q(1,100000000)
        whole=error+(summary_error+decoder_error)/2
        self.assertLess(whole,Q(1,25000000))
        OBS['analytic_prediction_contract']=dict(
            all_n_partial_atom_majorant_uses_same_K=True,K=K,
            zero_order_multiplicity='floor(r/3)+1',
            half_trace_norm_error='min(1,2*tanh(pi*K*abs(t)/2)^(floor(r/3)+1))',
            arbitrary_data_background_and_passive_reference=True,
            input_support_restricted_complete_channel_norm=True,
            unrestricted_all_graph_diamond_claim=False,
            radius=29,K_times_abs_t='1/10',rational_prediction_bound=str(error),
            summary_ordinary_diamond_error=str(summary_error),decoder_ordinary_diamond_error=str(decoder_error),
            combined_half_trace_error_bound=str(whole),
            H_A_is_only_the_mathematical_predictor_H_is_not_physically_truncated=True)

    def test_04_actual_H_vs_explicit_local_decoder_with_reference(self):
        trees,atoms,chosen,ids,gl,ge,il,ie,la=model()
        z=np.arange(256)[:,None];ref=np.arange(2)[None,:]
        x=np.zeros((256,90,2),complex)
        x[:,ids[0],:]=(z+1)*(ref+1)+1j*((3*z+ref)%17)
        x[:,ids[1],:]=((5*z+ref)%23)+1j*(z+3)*(2-ref)
        x/=np.linalg.norm(x)
        t=Q(1,64)
        actual=taylor_action(x,atoms,8,t)
        truncated=taylor_action(x,chosen,8,t)
        local=taylor_action(pack(x),la,6,t)
        inter=np.max(np.abs(pack(truncated)-local))
        self.assertLess(inter,1e-13)
        from_full=receiver_full(truncated);from_local=receiver_local(local)
        self.assertLess(np.max(np.abs(from_full-from_local)),1e-12)
        dist=prior.distance(receiver_full(actual),from_local)
        tail=tail_action(t)
        self.assertLess(tail,Q(1,10**20))
        self.assertLess(abs(np.linalg.norm(actual)-1),2*float(tail)+1e-12)
        self.assertGreater(np.linalg.norm(actual-truncated),1e-7)
        # The initial background is not restricted to a single excitation.
        probability_high=sum(np.sum(abs(x[a])**2) for a in range(256) if a.bit_count()>=2)
        self.assertGreater(probability_high,.5)
        OBS['full_background_decoder_diagnostic']=dict(
            actual_input_reference_dimension=2,nonzero_input_graphs=2,time=str(t),Taylor_degree=16,
            operator_action_tail_bound=str(tail),whole_state_intertwining_error=short(inter),
            actual_H_versus_local_prediction_half_trace_error=short(dist),
            probability_at_least_two_data_excitations=short(probability_high),
            actual_evolution_not_equal_to_truncated_global_evolution=True,
            formula_bound_at_this_diagnostic_is_the_trivial_cap_one=True,
            numeric_example_not_the_all_radius_or_all_reference_proof=True)

    def test_05_near_graph_phase_survives_fixed_summary_and_affects_true_readout(self):
        trees,atoms,chosen,ids,*_=model();vectors=phase_sources()
        b=np.array([(z>>7)&1 for z in range(256)])[:,None]
        values=[];zeros=[]
        for x in vectors:
            self.assertEqual(np.vdot(x,x).real,8)
            hx=apply_h(x,atoms,8);h2x=apply_h(hx,atoms,8)
            first=(1j*(np.vdot(hx,b*x)-np.vdot(x,b*hx))).real/8
            second=(2*np.vdot(x,b*h2x).real-2*np.vdot(hx,b*hx).real)/8
            zeros.append(first);values.append(second)
        self.assertEqual(zeros[0],zeros[1])
        self.assertEqual(values,[0.,1.])
        coefficient=Q(int(-(values[0]-values[1])),2);self.assertEqual(coefficient,Q(1,2))
        _,_,_,_,gl,ge,il,ie,la=model()
        self.assertEqual(ie[ids[0]],ie[ids[1]])
        self.assertNotEqual(il[ids[0]],il[ids[1]])
        # Fixed trace retains an orthogonal phase pair; the exterior is the same.
        self.assertAlmostEqual(np.vdot(vectors[0],vectors[1]).real,0)
        g,h=(trees[j] for j in ids)
        self.assertNotEqual(outside_record_signature(g,8,2),outside_record_signature(h,8,2))
        t=Q(1,100000)
        remainder=2*(30*t)**3/math.factorial(3)/(1-30*t/4)
        self.assertLess(remainder,t*t/4)
        OBS['actual_phase_interface']=dict(
            data_source='plus X at qubit0, plus Y at qubit1, other six qubits zero',
            graph_sources='(|g>+i|h>)/sqrt(2), (|g>-i|h>)/sqrt(2)',
            same_passive_reference_may_be_appended=True,
            actual_observable='excitation projector at fixed subject0',
            exact_second_commutator_expectations=[0,1],
            exact_probability_difference_taylor_coefficient=str(coefficient),
            certified_time=str(t),third_and_higher_remainder=str(remainder),
            true_probability_gap_greater_than=str(t*t/4),
            fixed_summary_preserves_graph_phase_pair_trace_distance_one=True,
            round443_radius_two_record_discards_this_graph_cross_term=True,
            separate_source_preparations_required_for_statistical_comparison=True)

    def test_06_random_global_radius_is_not_an_automatic_coherence_repair(self):
        examples=[]
        for length in (3,6):
            g,h,n=two_arms(length);maxradius=length+1
            adj=tree.adjacency(g,n)
            self.assertTrue(all(len(a) in (1,3) for a in adj))
            self.assertEqual(tree.tree_flips(g,n).get(h,0),1)
            sg=[outside_record_signature(g,n,r) for r in range(maxradius+1)]
            sh=[outside_record_signature(h,n,r) for r in range(maxradius+1)]
            self.assertTrue(all(x!=y for x in sg for y in sh))
            sizes_g=[len(neighborhood.ball(g,0,r,n)) for r in range(maxradius+1)]
            sizes_h=[len(neighborhood.ball(h,0,r,n)) for r in range(maxradius+1)]
            self.assertEqual(sizes_g,sizes_h)
            self.assertTrue(all(a<b for a,b in zip(sizes_g,sizes_g[1:])))
            examples.append(dict(vertices=n,arm_internal_length=length,maximum_radius=maxradius,
                                 all_cross_radius_external_pairs_orthogonal=len(sg)*len(sh)))
        t=Q(1,16*K**3)
        uniform_tail=2*(K*t)**3/(1-K*t)
        self.assertLess(uniform_tail,t*t/4)
        OBS['random_radius_boundary_and_resources']=dict(
            examples=examples,
            specific_protocol='round443 split, with one classical or coherent global radius register',
            any_amplitudes_supported_on_declared_finite_radius_window_lose_the_near_cross_term=True,
            not_a_no_go_for_all_soft_or_adaptive_localization_channels=True,
            uniform_phase_visibility_time=str(t),uniform_probability_gap_above=str(t*t/4),
            local_factor_count='|A|+|A|*(|A|-1)/2 data/edge qubits in original carrier',
            local_label_addresses_and_routing_are_inputs=True,
            envelope_can_equal_entire_graph_for_unrestricted_sources=True,
            summary_extraction_and_predictor_compilation_not_autonomously_supplied=True,
            no_cloning_unknown_summary_for_simultaneous_real_and_predicted_trial=True,
            full_information_retained_in_unread_exterior_when_trace_is_used=True,
            any_physical_realization_requires_clock_control_storage_and_full_channel_error_budget=True)


def run():
    OBS.clear();stream=io.StringIO()
    result=unittest.TextTestRunner(stream=stream).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():raise AssertionError(stream.getvalue())
    answer=dict(round=498,scientific_baseline_round=496,reused_frozen_rounds=[221,343,347,399,400,440,443,465,466],
        tests_run=result.testsRun,failures=len(result.failures),errors=len(result.errors),
        python=platform.python_version(),numpy=np.__version__,
        dependency_results_sha256={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in
            ('unknown_background_propagation_audit_results.json','quantum_neighborhood_signal_audit_results.json')},
        observations=OBS,scope=dict(
            explicit_fixed_physical_factor_CPTP_neighborhood_predictor=True,
            fixed_interface_CPTP_predictor_itself_reuses_round399=True,
            new_uniform_initial_support_envelope_bound_for_original_NNI=True,
            original_H_not_modified=True,initial_support_need_not_be_H_invariant=True,
            arbitrary_background_graph_coherence_and_passive_reference_preserved_in_contract=True,
            source_envelope_and_readout_permissions_are_extra_inputs=True,
            not_a_general_dynamic_graph_ball_predictor=True,
            not_an_abstract_minimal_predictive_algebra_rederivation=True,
            no_spatial_dimension_or_actual_coordinate_chart_derived=True,
            full_GR_goal_completed=False,phase_closure_triggered=False))
    assert json.loads(json.dumps(answer))==answer
    return answer


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dry-run',action='store_true');parser.add_argument('--check',action='store_true')
    args=parser.parse_args();answer=run()
    if args.check:assert answer==json.loads(TARGET.read_text(encoding='utf-8'))
    elif not args.dry_run:
        with TARGET.open('x',encoding='utf-8',newline='\n') as f:f.write(json.dumps(answer,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(answer,ensure_ascii=False,indent=2))

