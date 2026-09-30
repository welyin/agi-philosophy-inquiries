"""Round 489: the six-point mean-tree ruler has no smooth local tangent limit.

Scientific baseline 487; reuses 467/468/469/472/478. Fixed smooth point and
all six retained labels are explicit hypotheses, not cognitive axioms.
The result is about the original mean path LENGTHS, not the statistics TV metric.
"""
import argparse
from fractions import Fraction as F
import io
import itertools
import json
import math
from pathlib import Path
import platform
import unittest
import numpy as np
import averaged_tree_euclidean_audit as old468
import tree_metric_spatial_limit_audit as old467
import unit_tree_endpoint_metric_audit as old469
import sequential_tree_distance_readout_audit as reader
import stable_reference_coordinate_quotient as old478

TARGET=Path(__file__).with_name('local_metric_tangent_obstruction_results.json')
OBS={}
LABEL_MAP=(0,3,4,5,1,2)
WITNESS=np.array([1,-2,-2,1,1,1],dtype=np.int64)
LOWER=F(40,3)
SQUARED_GAP=F(5,24)
PAIRS=list(itertools.combinations(range(6),2))

def short(x):
    return float(f'{float(x):.13g}')

def frozen(module,number):
    r=json.loads(module.TARGET.read_text(encoding='utf-8'))
    assert (r['round'],r['tests_run'],r['failures'],r['errors'])==(number,6,0,0)
    return r

def distance_tables():
    trees=reader.system()[0]
    return [np.array([[len(reader.old.path_between(g,6,a,b))-1 for b in range(6)]
                      for a in range(6)],dtype=np.int64) for g in trees]

def mean_matrix(weights):
    return sum((F(w)*d.astype(object) for w,d in zip(weights,distance_tables())),
               np.full((6,6),F(0),dtype=object))

def admissible_paths():
    trees=reader.system()[0];out={}
    for a,b in PAIRS:
        paths=[]
        for path in reader.all_paths(6,a,b):
            edges=reader.path_edges(path)
            support=np.array([int(all(e in g for e in edges)) for g in trees],dtype=np.int64)
            if np.any(support):paths.append((path,support))
        out[(a,b)]=paths
    return out

def ceil_fraction(v):
    return (v.numerator+v.denominator-1)//v.denominator

class Audit(unittest.TestCase):
    def test_01_frozen_witness_and_current_label_polynomial(self):
        prior=frozen(old468,468)
        self.assertEqual(prior['observations']['all_six_vertices_diagnostic']['universal_conservative_lower_bound'],'40/3')
        old_trees=old468.binary_six_trees()
        old_tables=[old468.ordinary_distances(g,6) for _,g in old_trees]
        current=distance_tables()
        for old in old_tables:
            mapped=np.empty((6,6),dtype=np.int64)
            for i,j in itertools.product(range(6),repeat=2):mapped[LABEL_MAP[i],LABEL_MAP[j]]=old[i,j]
            self.assertTrue(any(np.array_equal(mapped,d) for d in current))
        self.assertEqual(int(WITNESS.sum()),0)
        self.assertEqual(int(sum(abs(WITNESS)))**2,64)
        self.assertEqual(LOWER/64,SQUARED_GAP)
        trees=reader.system()[0]
        incidence=np.array([[int(reader.old.edge(1,j) in g) for g in trees] for j in (0,3,4,5)])
        matching=np.array([[int((frozenset(j for j in (0,3,4,5) if reader.old.edge(1,j) in g)
                                      in (frozenset((0,b)),frozenset(set((0,3,4,5))-{0,b}))))
                            for g in trees] for b in (3,4,5)],dtype=np.int64)
        quadratic=np.array([[WITNESS@(a*b)@WITNESS for b in current] for a in current])
        expected=28*np.ones((6,6),dtype=np.int64)+4*matching.T@matching-8*incidence.T@incidence
        self.assertTrue(np.array_equal(quadratic,expected))
        self.assertTrue(np.array_equal(matching.sum(axis=0),np.ones(6,dtype=np.int64)))
        self.assertTrue(np.array_equal(incidence.sum(axis=0),2*np.ones(6,dtype=np.int64)))
        OBS['frozen_witness_bridge']=dict(source_round=468,
            original_to_current_labels=list(LABEL_MAP),witness=WITNESS.tolist(),
            exact_population_quadratic_matrix=quadratic.tolist(),
            exact_formula='28+4 sum matching_probability^2-8 sum source_leaf_adjacency^2',
            all_population_quantifier_from_frozen_analytic_bounds_not_samples=True,
            lower=str(LOWER),l1_squared=64,squared_tangent_error_lower=str(SQUARED_GAP),
            lower_bound_not_claimed_optimal=True,
            original_mean_path_lengths_not_mean_squared_lengths=True)

    def test_02_finite_error_and_changing_population_certificates(self):
        # Scale-dependent populations have no assumed convergence.
        records=[]
        for n in (4,5,8,9):
            scale=F(1,2**n);p=[F(0)]*6
            p[n%6]=1-scale;p[(n+1)%6]=scale
            m=mean_matrix(p)
            witness=WITNESS@(m*m)@WITNESS
            scaled=scale*scale*witness
            self.assertGreaterEqual(scaled,LOWER*scale*scale)
            records.append(dict(n=n,a=str(scale),p=[str(v) for v in p],
                                normalized_witness=str(witness),scaled_witness=str(scaled)))
        epsilon=F(1,100);radius=F(2);omega=F(1,1000)
        tangent=4*radius*radius*omega
        length_to_square=epsilon*(6+epsilon)
        slack=SQUARED_GAP-length_to_square-tangent
        self.assertGreater(slack,0)
        measurement=F(1,1000)
        measured_witness_lower=LOWER-64*measurement*(6+measurement)
        self.assertGreater(measured_witness_lower,0)
        OBS['finite_error_certificate']=dict(
            generic_necessary_inequality='epsilon*(6+epsilon)+eta >= 5/24',
            epsilon_is_max_length_error_divided_by_scale=True,
            eta_is_max_tangent_squared_distance_error_divided_by_scale_squared=True,
            metric_matrix_comparison='(1-omega)I <= g <= (1+omega)I on normal ball of radius 3*R*a',
            sufficient_tangent_eta='4*R^2*omega',example_R=str(radius),example_omega=str(omega),
            example_length_epsilon=str(epsilon),example_tangent_eta=str(tangent),
            strict_rejection_slack=str(slack),nonconvergent_population_examples=records,
            measurement_error_per_unscaled_distance=str(measurement),
            observed_witness_lower=str(measured_witness_lower),
            source_and_reader_error_cannot_be_omitted=True)

    def test_03_fixed_smooth_metric_positive_control(self):
        # Unit two-sphere, exact exp at north pole; numerical diagnostics only.
        vectors=np.array([[0,0],[1,0],[-1,0],[0,1],[0,-1],[1,1]],dtype=float)
        tangent=np.linalg.norm(vectors[:,None,:]-vectors[None,:,:],axis=-1)**2
        records=[]
        previous=None
        for n in (3,4,5,6):
            a=2.0**(-n);points=[]
            for v in vectors:
                r=a*np.linalg.norm(v)
                points.append(np.r_[a*v*(math.sin(r)/r if r else 1),math.cos(r)])
            points=np.asarray(points)
            distances=np.zeros((6,6))
            for i,j in PAIRS:
                distances[i,j]=distances[j,i]=math.atan2(np.linalg.norm(np.cross(points[i],points[j])),
                                                       float(points[i]@points[j]))
            err=float(np.max(abs(distances**2-a*a*tangent)))
            normalized=err/(a*a)
            if previous is not None:self.assertLess(normalized,previous/3)
            previous=normalized
            self.assertLess(normalized,float(SQUARED_GAP))
            records.append(dict(a=short(a),maximum_squared_tangent_error=short(err),
                                normalized_error=short(normalized)))
        OBS['smooth_positive_control']=dict(
            manifold='fixed unit two-sphere at north pole',tangent_vectors=vectors.tolist(),
            max_normal_coordinate_norm_bound_R=2,records=records,
            numerical_decay_not_used_to_prove_uniform_limit=True,
            general_proof_is_normal_neighborhood_metric_comparison=True,
            applies_to_any_fixed_finite_Riemannian_dimension=True,
            does_not_identify_these_sphere_distances_with_the_six_tree_means=True)

    def test_04_degenerating_and_projected_boundary_contracts(self):
        a467=frozen(old467,467);a469=frozen(old469,469);a478=frozen(old478,478)
        self.assertTrue(a469['scope']['full_GR_goal_completed'] is False)
        rows=[]
        # Frozen 468 short-third-leaf example: dimensionless edge ratio tends to zero.
        # sqrt(1+e^2) lies between 1 and 1+e^2/2; no floating root needed.
        for n in (4,8,12):
            a=F(1,2**n);e=a
            minimum=1+e*e/2
            self.assertGreaterEqual(minimum*minimum,1+e*e)
            self.assertGreater(e-e*e/2,0)
            error_upper=a*e
            self.assertEqual(error_upper/(a*a),1)
            rows.append(dict(a=str(a),relative_short_edge=str(e),absolute_short_edge=str(a*e),
                             max_length_error_upper=str(error_upper),relative_error_upper=str(e)))
        # Four retained uniform leaves alone are the old regular tetrahedron case.
        length=F(8,3)
        gram=(length*length/2)*(np.eye(4,dtype=int).astype(object)-np.full((4,4),F(1,4),dtype=object))
        projection=4*np.eye(4,dtype=np.int64)-np.ones((4,4),dtype=np.int64)
        self.assertTrue(np.array_equal(projection@projection,4*projection))
        self.assertEqual(int(np.trace(projection)),12)
        self.assertTrue(all(gram[i,i]+gram[j,j]-2*gram[i,j]==length*length for i,j in itertools.combinations(range(4),2)))
        OBS['scope_boundaries']=dict(
            frozen_boundaries_reused=[467,468,469,478],
            short_edge_examples=rows,short_edge_changes_dimensionless_shape_so_not_fixed_unit_six_tree_contract=True,
            fixed_six_label_uniform_leaves_alone_have_regular_tetrahedron_Gram=True,
            excluding_two_internal_labels_changes_the_observation_contract=True,
            round469_increasing_N_selected_leaf_approximation_not_excluded=True,
            same_all_label_obstruction_not_transferred_to_endpoint_projection=True,
            round478_statistics_TV_metric_not_a_two_subject_length=True,
            moving_basepoint_requires_uniform_tangent_approximation=True,
            bounded_curvature_alone_not_substitute_for_normal_neighborhood_control=True,
            no_curvature_divergence_claim=True,
            no_exclusion_of_finite_curved_six_point_realization=True)

    def test_05_actual_path_reader_and_finite_statistics_ledger(self):
        paths=admissible_paths();tables=distance_tables()
        total=0;weightmax=0;maxlen=0
        for pair,records in paths.items():
            a,b=pair;total+=len(records)
            recovered=np.zeros(6,dtype=np.int64)
            supports=np.array([s for _,s in records])
            self.assertTrue(np.array_equal(supports.sum(axis=0),np.ones(6,dtype=np.int64)))
            weightsum=0
            for path,support in records:
                length=len(path)-1;maxlen=max(maxlen,length);weightsum+=length
                recovered+=length*support
            self.assertTrue(np.array_equal(recovered,np.array([d[a,b] for d in tables])))
            weightmax=max(weightmax,weightsum)
        self.assertEqual((total,weightmax,maxlen),(41,10,3))
        # One actual longest-path instrument with all signed classical records,
        # on the original graph/reference preparation; not a fresh model.
        probe=F(1,512);path=(0,1,2,3)
        estimates=reader.signed_reference_for_path(path,probe)
        initial=reader.system()[-1]
        projection=np.eye(6,dtype=np.int64)
        for edge in reader.path_edges(path):projection=projection@reader.graph_projection(edge)
        product=np.kron(projection,np.eye(3))
        ideal=reader.reference_marginal(product@initial@product)
        x=18*probe;r=x*x/(2*(1-x/3));delta=r/probe
        error=reader.trace_norm(estimates-ideal)
        diagnostic_bound=(1+delta)**3-1
        self.assertLess(error,float(diagnostic_bound)+2e-9)
        # Common source before every path choice permits source error <=3*gamma_src.
        epsilon=F(1,1000);source=epsilon/6;zeta=epsilon/(2*weightmax)
        h=zeta/(48*maxlen*81);gamma=zeta*h/(8*maxlen)
        x=18*h;r=x*x/(2*(1-x/3));delta=(r+gamma)/h
        self.assertLessEqual((1+delta)**maxlen-1,zeta/2)
        copies=ceil_fraction(80/(zeta*zeta*h**(2*maxlen)))
        self.assertGreaterEqual(copies*zeta*zeta*h**(2*maxlen)/8,10)
        self.assertGreater(sum(F(10**k,math.factorial(k)) for k in range(21)),8200)
        self.assertEqual(3*source+weightmax*zeta,epsilon)
        OBS['actual_reader_and_resource_contract']=dict(
            source_instrument_round=472,all_pairs=15,nonzero_path_programs=total,
            maximum_path_weight_sum=weightmax,maximum_path_length=maxlen,
            exact_path_diagonal_decomposition_in_original_six_graphs=True,
            diagnostic_path=list(path),diagnostic_probe_wait=str(probe),
            diagnostic_reference_matrix_error=short(error),
            diagnostic_error_upper=str(diagnostic_bound),
            diagnostic_role_only_not_tangent_limit_measurement=True,
            unscaled_distance_error=str(epsilon),scaled_distance_error='a/1000',
            common_source_diamond_error=str(source),each_path_product_error=str(zeta),
            certified_probe_wait=str(h),complete_each_stage_instrument_diamond_error=str(gamma),
            independent_copies_per_path=str(copies),total_known_history_preparations=str(total*copies),
            joint_failure_probability_less_than='1/100',
            common_actual_source_for_all_path_settings_required=True,
            old_complete_data_moved_to_internal_storage_and_graph_remains_active=True,
            old_storage_graph_reference_correlations_retained=True,
            every_stage_uses_fresh_probe_and_preserves_old_probe_information_in_storage=True,
            no_old_random_record_feedback_or_unknown_state_cloning=True,
            no_sampling_experiment_executed=True,
            graph_readout_disturbance_allowed_not_reversible_measurement_claim=True,
            preparation_clock_storage_control_and_repetition_are_extra_inputs=True)

def run():
    OBS.clear();stream=io.StringIO()
    result=unittest.TextTestRunner(stream=stream).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():raise RuntimeError(stream.getvalue())
    output=dict(round=489,baseline_round=487,scientific_baselines=[467,468,469,472,478,487],
        tests_run=result.testsRun,failures=len(result.failures),errors=len(result.errors),
        python=platform.python_version(),numpy=np.__version__,observations=OBS,
        scope=dict(
            all_six_original_mean_path_distance_readouts_and_uniform_scale_required=True,
            populations_may_vary_arbitrarily_with_scale=True,
            fixed_smooth_Riemannian_point_and_O_scale_normal_neighborhood_required=True,
            additive_length_error_little_o_scale_cannot_remove_obstruction=True,
            moving_centers_only_under_explicit_uniform_tangent_approximation=True,
            no_claim_curvature_must_diverge=True,
            no_exclusion_of_topological_space_curved_global_geometry_or_different_readout=True,
            no_general_cognitive_three_dimensional_no_go=True,
            independent_of_parallel_round488=True,
            full_GR_goal_completed=False,phase_closure_triggered=False))
    assert json.loads(json.dumps(output))==output
    return output

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dry-run',action='store_true');parser.add_argument('--check',action='store_true')
    args=parser.parse_args();r=run()
    if args.check:assert r==json.loads(TARGET.read_text(encoding='utf-8'))
    elif not args.dry_run:
        with TARGET.open('x',encoding='utf-8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(r,ensure_ascii=False,indent=2))
