"""Round 497: a hard graph-ball channel loses dynamically visible local phase.

General statements are analytic; finite examples check the original NNI rule,
exact record channel, Gaussian-integer commutators, and finite error budgets.
No result is overwritten: first saving uses x; --dry-run writes nothing.
"""
import argparse
from collections import defaultdict
from fractions import Fraction as F
import hashlib
import io
import itertools
import json
import math
from pathlib import Path
import platform
import unittest

import numpy as np
import branching_tree_distance_audit as tree_model
import quantum_neighborhood_signal_audit as old_channel

TARGET = Path(__file__).with_name('neighborhood_phase_boundary_obstruction_results.json')
OBS = {}


def comb_pair(radius):
    assert isinstance(radius, int) and radius >= 2
    n = 2*radius+2
    branch = [3]+list(range(6, radius+4))
    side = list(range(radius+4, n))
    edges = {(0,1),(1,2),(2,5),(2,4),(1,3)}
    for j, leaf in enumerate(side):
        edges.add(tree_model.edge(branch[j], branch[j+1]))
        edges.add(tree_model.edge(branch[j], leaf))
    g = frozenset(edges)
    h = g.symmetric_difference({(0,1),(2,5),(0,2),(1,5)})
    return n, g, h, branch[-1]


def split_record(tree, bits, radius):
    n = len(bits)
    region = old_channel.ball(tree, 0, radius, n)
    data = [('D',v,bits[v]) for v in range(n)]
    edges = [('E',a,b,int((a,b) in tree))
             for a,b in itertools.combinations(range(n),2)]
    local = tuple(x for x in data if x[1] in region)+tuple(
        x for x in edges if x[1] in region or x[2] in region)
    external = tuple(x for x in data if x[1] not in region)+tuple(
        x for x in edges if x[1] not in region and x[2] not in region)
    return local, external


def gm(a,b):
    return a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0]


def gc(a):
    return a[0],-a[1]


def state_records(radius, sign):
    n,g,h,_ = comb_pair(radius)
    out=[]
    for tree,phase in ((g,(1,0)),(h,(0,sign))):
        for x,y in itertools.product(range(2),repeat=2):
            bits=(x,y)+(0,)*(n-2)
            out.append((*split_record(tree,bits,radius),gm((1,0) if y==0 else (0,1),phase)))
    return out  # squared norm 8, all amplitudes Gaussian integers


def reduced_integer_matrix(records):
    out=defaultdict(lambda:(0,0))
    for l,e,a in records:
        for m,f,b in records:
            if e==f:
                q=gm(a,gc(b)); old=out[l,m]
                out[l,m]=(old[0]+q[0],old[1]+q[1])
    return {key:value for key,value in out.items() if value!=(0,0)}


def swap_matrix(n,a,b):
    out=np.zeros((2**n,2**n),dtype=np.int64)
    for z in range(2**n):
        mask=((z>>(n-1-a))^(z>>(n-1-b)))&1
        w=z^((mask<<(n-1-a))|(mask<<(n-1-b)))
        out[w,z]=1
    return out


def six_model():
    trees,_,legacy,flip,_=tree_model.six_vertex_sector()
    assert np.array_equal(legacy,legacy.real) and np.array_equal(legacy.real,np.rint(legacy.real))
    h=np.kron(np.eye(64,dtype=np.int64),np.asarray(flip,dtype=np.int64))
    for j,g in enumerate(trees):
        block=sum((swap_matrix(6,a,b) for a,b in g),np.zeros((64,64),dtype=np.int64))
        ids=np.arange(64)*6+j
        h[np.ix_(ids,ids)]+=block
    assert np.array_equal(h,legacy)
    return trees,h


def six_sources(trees):
    _,g,h,_=comb_pair(2)
    out=[]
    for sign in (1,-1):
        real=np.zeros(384,dtype=np.int64); imag=real.copy()
        for tree,phase in ((g,(1,0)),(h,(0,sign))):
            for x,y in itertools.product(range(2),repeat=2):
                z=(x<<5)|(y<<4)
                value=gm((1,0) if y==0 else (0,1),phase)
                index=6*z+trees.index(tree)
                real[index],imag[index]=value
        out.append((real,imag))
    return out


def quadratic_integer(a, pair):
    real,imag=pair
    return int(real@(a@real)+imag@(a@imag)), int(real@(a@imag)-imag@(a@real))


def exact_coefficients(h,sources,degree):
    # Every row/column of H has absolute sum <=9. Before multiplication,
    # |ad_H^j B|_max <=18^j, and each quadratic has at most 8x8 terms.
    # This prior absolute-sum bound controls all int64 intermediate sums.
    assert int(np.max(np.sum(np.abs(h),axis=0)))<=9
    assert int(np.max(np.sum(np.abs(h),axis=1)))<=9
    assert 256*18**degree < 2**63
    effect=np.repeat([int((z>>5)&1) for z in range(64)],6)
    a=np.diag(np.asarray(effect,dtype=np.int64))
    coefficients=[]
    for k in range(degree+1):
        plus=quadratic_integer(a,sources[0]); minus=quadratic_integer(a,sources[1])
        q=(plus[0]-minus[0],plus[1]-minus[1])
        rotated=gm(q,((1,0),(0,1),(-1,0),(0,-1))[k%4])
        assert rotated[1]==0
        coefficients.append(F(rotated[0],8*math.factorial(k)))
        if k<degree:a=h@a-a@h
    return coefficients


def uniform_certificate():
    c=3; lam=c+2*c*(c-1)**2; k=96*lam*c**6
    time=F(1,8*k**3)
    tail=(k*time)**3/(1-k*time)
    gap=time*time/4
    assert tail<time*time/4
    return dict(K=k,time=time,tail=tail,gap=gap,worst_prediction_event_error=gap/2)


def ceil_fraction(x):
    return -(-x.numerator//x.denominator)


class Audit(unittest.TestCase):
    def test_01_all_radius_geometry_and_original_nni(self):
        rows=[]
        for r in range(2,17):
            n,g,h,z=comb_pair(r)
            deg=[len(x) for x in tree_model.adjacency(g,n)]
            self.assertEqual(sorted(deg),[1]*(r+2)+[3]*r)
            self.assertEqual(deg,[len(x) for x in tree_model.adjacency(h,n)])
            self.assertEqual(len(g),n-1)
            dg=len(tree_model.path_between(g,n,0,z))-1
            dh=len(tree_model.path_between(h,n,0,z))-1
            self.assertEqual((dg,dh),(r,r+1))
            self.assertEqual(tree_model.tree_flips(g,n)[h],1)
            self.assertEqual(tree_model.tree_flips(h,n)[g],1)
            rows.append(dict(radius=r,N=n,terminal=z,distances=[dg,dh],NNI_amplitude=1))
        OBS['all_radius_construction']=dict(finite_crosschecks=rows,
            analytic_family='N=2r+2, r>=2; one near-root NNI swaps 01,25 to 02,15',
            full_fixed_degree_tree_sector_is_individually_NNI_invariant=True,
            proof_for_all_r_not_inferred_from_finite_list=True)

    def test_02_exact_equal_neighborhood_inputs(self):
        for r in range(2,13):
            n,g,h,z=comb_pair(r)
            gp=state_records(r,1); gm_=state_records(r,-1)
            self.assertEqual(reduced_integer_matrix(gp),reduced_integer_matrix(gm_))
            # Basis-factor identity, not its value, distinguishes environments.
            eg=split_record(g,(0,)*n,r)[1]; eh=split_record(h,(0,)*n,r)[1]
            self.assertNotIn(('D',z,0),eg); self.assertIn(('D',z,0),eh)
            self.assertTrue(any(x[0]=='E' and x[-1]==0 for x in eg+eh))
            diag=sum(v[0] for (a,b),v in reduced_integer_matrix(gp).items() if a==b)
            self.assertEqual(diag,8)
        n,g,h,_=comb_pair(2)
        for tree in (g,h):
            for z in range(64):
                self.assertEqual(split_record(tree,old_channel.bits(z),2),
                                 old_channel.split_record(tree,z,center=0,radius=2))
        # When the SAME finite system is wholly included, no phase is lost.
        def full_records(sign):
            return [(tuple(sorted(l+e)),(),amp) for l,e,amp in state_records(2,sign)]
        self.assertNotEqual(reduced_integer_matrix(full_records(1)),
                            reduced_integer_matrix(full_records(-1)))
        OBS['legal_channel_kernel_witness']=dict(radii_checked=list(range(2,13)),
            exact_integer_reduced_matrices_equal=True,
            equality_holds_for_all_data_cross_terms_by_external_D_factor_identity=True,
            legacy_443_basis_records_compared=128,
            source_squared_normalization=8,
            full_finite_network_channel_still_recovers_phase=True,
            CPTP_legality_and_time_zero_bare_readout_not_refuted=True)

    def test_03_local_symbolic_and_full_model_second_order(self):
        s01=swap_matrix(3,0,1);s02=swap_matrix(3,0,2)
        z=np.diag([1 if k<4 else -1 for k in range(8)]).astype(np.int64)
        vreal=np.zeros(8,dtype=np.int64);vimag=vreal.copy()
        vreal[0]=vreal[4]=1;vimag[2]=vimag[6]=1
        q=quadratic_integer((s01-s02)@z-z@(s01-s02),(vreal,vimag))
        self.assertEqual(q,(0,4))  # /4 = i
        trees,h=six_model();sources=six_sources(trees)
        coeff=exact_coefficients(h,sources,2)
        self.assertEqual(coeff,[F(0),F(0),F(1,2)])
        OBS['actual_receiver_second_order']=dict(
            local_commutator_expectation='i',
            data_state='|+X>_0 |+Y>_1 |0>_others',
            graph_states='(|G> +/- i|H>)/sqrt(2)',
            measured_effect='|1><1| on physical D_0, identity elsewhere',
            exact_probability_difference_coefficients=[str(x) for x in coeff],
            full_unchanged_H_dimension=384,
            symbolic_formula='ad_H^2(B)_{GH}=[D_H-D_G,B]',
            second_derivative_not_in_Phi_r_dual_image=True)

    def test_04_uniform_finite_time_and_factorization_obstruction(self):
        cert=uniform_certificate()
        self.assertEqual(cert['K'],1889568)
        self.assertLess(cert['tail']/cert['time']**2,F(1,7))
        self.assertGreater(cert['time']**2/2-cert['tail'],cert['gap'])
        OBS['uniform_analytic_certificate']={key:str(value) for key,value in cert.items()}
        OBS['uniform_analytic_certificate'].update(
            applies_to_all_r_ge_2=True,
            tail_formula='(K*t)^3/(1-K*t), for 0<t<=1/(8*K^3)',
            ordinary_diamond_error_lower_bound=str(cert['gap']),
            uniformity_is_across_N_equal_2r_plus_2_not_fixed_N=True,
            scalar_event_error_lower_bound_needs_no_reference=True,
            round466_norm_bound_used_but_remote_data_encoding_premise_not_assumed=True)

    def test_05_finite_six_subject_readout_certificate(self):
        trees,h=six_model();sources=six_sources(trees)
        coefficients=exact_coefficients(h,sources,8)
        t=F(1,64)
        polynomial=sum(c*t**j for j,c in enumerate(coefficients))
        # Centering B removes scalar. ||rho+ -rho-||_1=2 cancels ||B-I/2||=1/2.
        tail=(18*t)**9/(math.factorial(9)*(1-18*t/10))
        low=polynomial-tail;high=polynomial+tail
        self.assertGreater(low,F(1,10000))
        evals,evecs=np.linalg.eigh(h.astype(float))
        u=(evecs*np.exp(-1j*float(t)*evals))@evecs.T
        mask=np.repeat([bool((z>>5)&1) for z in range(64)],6)
        probabilities=[]
        for real,imag in sources:
            psi=u@(real+1j*imag)/math.sqrt(8)
            probabilities.append(float(np.vdot(psi[mask],psi[mask]).real))
            # The actual CP instrument has two Kraus projectors on D_0.
            self.assertAlmostEqual(float(np.vdot(psi,psi).real),1,places=12)
        self.assertGreater(probabilities[0]-probabilities[1],float(low)-1e-12)
        self.assertLess(probabilities[0]-probabilities[1],float(high)+1e-12)
        OBS['finite_actual_readout']=dict(time=str(t),
            exact_degree_eight_coefficients=[str(x) for x in coefficients],
            polynomial=str(polynomial),remainder=str(tail),
            rigorous_difference_interval=[str(low),str(high)],
            diagnostic_probabilities=[float(f'{p:.13g}') for p in probabilities],
            diagnostic_probability_difference=float(f'{probabilities[0]-probabilities[1]:.13g}'),
            exact_int64_prior_absolute_sum_bound=256*18**8,
            actual_H_and_bare_two_outcome_port_instrument=True,
            floating_subtraction_not_used_at_uniform_tiny_time=True)

    def test_06_source_phase_and_complete_finite_ledger(self):
        # Explicit finite source permission, NOT autonomous control derived from H.
        # C=(|0> +/- i|1>)/sqrt2; four CNOTs complement the changed edge bits.
        # A CNOT from edge01 to C, then X_C, returns C to zero on both branches.
        def prep_permutation(c,edge_bits):
            x=list(edge_bits)
            if c:x=[1-v for v in x]
            c^=x[0];c^=1
            return c,tuple(x)
        all_out={prep_permutation(c,x) for c in range(2)
                 for x in itertools.product(range(2),repeat=4)}
        self.assertEqual(len(all_out),32)
        self.assertEqual(prep_permutation(0,(1,1,0,0)),(0,(1,1,0,0)))
        self.assertEqual(prep_permutation(1,(1,1,0,0)),(0,(0,0,1,1)))
        cert=uniform_certificate();gap=cert['gap'];k=cert['K']
        source=gap/16;dynamics=gap/16;clock=gap/(16*k);reader=gap/32
        per_setting=source/2+dynamics/2+k*clock/2+reader
        self.assertEqual(per_setting,gap/8)
        accuracy=gap/8;samples=ceil_fraction(160/(gap*gap))
        self.assertGreaterEqual(2*samples*accuracy**2,5)
        self.assertGreater(sum(F(5)**j/math.factorial(j) for j in range(6)),80)
        self.assertEqual(gap-2*per_setting-2*accuracy,gap/2)
        rows=[]
        for radius in (2,3,10):
            n=2*radius+2;e=n*(n-1)//2;gates=n+10
            self.assertEqual(gap/32+gates*(gap/(32*gates)),source)
            rows.append(dict(radius=radius,N=n,listed_carriers_per_trial=n+e+3,
                source_gate_count=gates,per_gate_diamond_budget=str(gap/(32*gates)),
                total_listed_carriers_for_both_banks=str(2*samples*(n+e+3))))
        OBS['finite_source_and_measurement_contract']=dict(
            source_permutation_dimension=32,
            blank_source_preparation_trace_budget=str(gap/32),
            full_source_trace_budget_including_phase_reference=str(source),
            phase_gate_accuracy_included_in_each_gate_diamond_budget=True,
            coherent_preparation_gates_and_isolation_are_additional_permissions=True,
            intermediate_raw_edge_register_states_need_not_be_trees=True,
            endpoint_source_is_a_legal_tree_superposition=True,
            no_free_preparation_by_original_H_claim=True,
            dynamics_ordinary_diamond_budget=str(dynamics),
            absolute_clock_tolerance=str(clock),
            reader_effect_operator_error=str(reader),
            complete_event_error_per_setting=str(per_setting),
            samples_per_setting=str(samples),empirical_mean_accuracy_each=str(accuracy),
            statistical_failure_bound='1/20',
            surviving_empirical_gap=str(gap/2),
            independent_same_actual_protocol_trials_are_an_extra_input=True,
            all_source_controllers_old_trials_and_reading_records_are_retained=True,
            listed_resource_examples=rows,
            gate_hardware_clock_phase_reference_and_routing_not_derived_or_costed_as_free=True,
            total_original_H_wait_for_both_banks=str(2*samples*cert['time']),
            huge_sampling_and_hardware_construction_not_executed=True)


def run():
    OBS.clear()
    stream=io.StringIO()
    result=unittest.TextTestRunner(stream=stream,verbosity=0).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():raise AssertionError(stream.getvalue())
    names=['graph_blind_scaling_limit_results.json',
           'quantum_neighborhood_signal_audit_results.json',
           'unknown_background_propagation_audit_results.json']
    dependencies={name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                  for name in names}
    return dict(round=497,scientific_baseline_round=496,reused_frozen_rounds=[442,443,466],
        tests_run=result.testsRun,failures=len(result.failures),errors=len(result.errors),
        python=platform.python_version(),numpy=np.__version__,dependency_results_sha256=dependencies,
        scope=dict(specified_443_hard_ball_channel_only=True,
            arbitrary_unknown_graph_coherence_is_in_candidate_domain=True,
            uniform_in_N_prediction_sufficiency_refuted=True,
            legal_CPTP_channel_and_time_zero_compatibility_preserved=True,
            no_superluminal_or_universal_spatial_no_go_claim=True,
            no_remote_data_encoding_counterexample_claim=True,
            full_GR_goal_completed=False,phase_closure_triggered=False),observations=OBS)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--dry-run',action='store_true')
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args();data=run()
    if args.check:assert json.loads(TARGET.read_text(encoding='utf-8'))==data
    elif not args.dry_run:
        with TARGET.open('x',encoding='utf-8') as handle:
            json.dump(data,handle,ensure_ascii=False,indent=2);handle.write('\n')
    print(json.dumps(data,ensure_ascii=False,indent=2))
