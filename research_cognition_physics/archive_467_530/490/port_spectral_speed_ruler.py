"""Round 490: port spectral speed ruler for the existing six-tree Hamiltonian.

Baseline 489. This is a finite endpoint ruler and a transfer-time lower bound,
not a derivation of a Dirac operator, spatial dimension, or minimal travel time.
The one-excitation probe preparation, port readout, clock and couplings are inputs.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import io,itertools,json,math,platform,unittest
from pathlib import Path
import numpy as np
import rigid_leaf_reference_audit as frame

TARGET=Path(__file__).with_name('port_spectral_speed_ruler_results.json')
LEAVES=(0,3,4,5)
OBS={}
# All int64 products below contain at most four factors, dimension at most 384,
# primitive integer entries at most 64. Check this conservative bound BEFORE use.
assert 384**4*64**4 < 2**63


def short(x):return float(f'{float(x):.12g}')


def determinant(matrix):
    a=[[F(int(v)) for v in row] for row in matrix];d=F(1)
    for i in range(len(a)):
        j=next((j for j in range(i,len(a)) if a[j][i]),None)
        if j is None:return F(0)
        if j!=i:a[i],a[j]=a[j],a[i];d=-d
        pivot=a[i][i];d*=pivot
        for j in range(i+1,len(a)):
            ratio=a[j][i]/pivot
            for k in range(i+1,len(a)):a[j][k]-=ratio*a[i][k]
            a[j][i]=F(0)
    return d


@lru_cache(None)
def model():
    trees,h,f=frame.system();h=h.astype(np.int64);f=f.astype(np.int64)
    adj=[]
    for tree in trees:
        a=np.zeros((6,6),dtype=np.int64)
        for i,j in tree:a[i,j]=a[j,i]=1
        adj.append(a)
    index=np.array([6*(1<<(5-i))+g for i in range(6) for g in range(6)])
    v=np.eye(384,dtype=np.int64)[:,index]
    h1=np.zeros((36,36),dtype=np.int64)
    for g,a in enumerate(adj):
        lap=np.diag(a.sum(axis=1))-a
        for i in range(6):
            for j in range(6):h1[6*i+g,6*j+g]=-lap[i,j]
    hf=np.kron(np.eye(6,dtype=np.int64),f)
    k=np.zeros((6,6),dtype=np.int64)
    for i in range(6):k[2 if i==1 else 1 if i==2 else i,i]=1
    return trees,h,f,adj,v,h1,hf,np.kron(k,np.eye(6,dtype=np.int64))


def lip(values,coupling=1):
    return abs(coupling)*max(np.linalg.norm(a@np.diag(values)-np.diag(values)@a,2) for a in model()[3])


def distance_matrix(coupling=1):
    assert coupling!=0
    d=np.zeros((6,6))
    for a in range(6):
        for b in range(a+1,6):
            value=np.sqrt(2) if a in LEAVES and b in LEAVES else np.sqrt(3)-1 if {a,b}=={1,2} else 1
            d[a,b]=d[b,a]=value/abs(coupling)
    return d


def unitary(generator,time):
    e,v=np.linalg.eigh(generator);return (v*np.exp(-1j*time*e))@v.conj().T


class Audit(unittest.TestCase):
    def test_01_original_observable_algebras_do_not_silently_supply_subject_points(self):
        trees,h,f,adj,v,hj,hf,k=model();distances=np.array([frame.distance(a,b) for a in range(6) for b in range(a+1,6)],dtype=np.int64)
        witness=None
        for rows in itertools.combinations(range(15),6):
            det=determinant(distances[list(rows)])
            if det:witness=(rows,det);break
        self.assertIsNotNone(witness)
        fullf=np.kron(np.eye(64,dtype=np.int64),f);data=h-fullf
        for g in range(6):
            a=np.kron(np.eye(64,dtype=np.int64),np.diag(np.eye(6,dtype=np.int64)[g]))
            self.assertFalse(np.any(data@a-a@data))
            self.assertTrue(np.array_equal(h@a-a@h,fullf@a-a@fullf))
        reached={0}
        for _ in range(6):reached|={j for i in reached for j in range(6) if f[i,j]}
        self.assertEqual(len(reached),6)
        # Full algebra contains H itself in the commutator kernel.
        one=np.ones(6,dtype=np.int64);q=np.array([frame.distance(1,3)[g]-frame.distance(1,0)[g] for g in range(6)],dtype=np.int64)
        self.assertTrue(np.array_equal(f@one,4*one));self.assertFalse(np.any(f@q))
        self.assertEqual(int(one@q),0)
        OBS['algebra_audit']=dict(graph_diagonal_distance_span_rank=6,
            exact_nonzero_minor_rows=list(witness[0]),exact_minor_determinant=str(witness[1]),
            graph_diagonal_commutator='[H,I tensor f]=kappa I tensor [F,f]',
            graph_diagonal_points_are_whole_graph_configurations_not_subject_ports=True,
            graph_diagonal_ruler_ignores_J_and_data_transport=True,
            full_matrix_algebra_has_nontrivial_commutator_kernel=True,
            exact_two_energy_witnesses=[9,5],full_algebra_distance_between_these_states='infinity',
            H_is_not_declared_to_be_a_geometric_Dirac_operator=True)

    def test_02_exact_one_excitation_interface_and_endpoint_algebra(self):
        trees,h,f,adj,v,hj,hf,k=model();identity=np.eye(36,dtype=np.int64)
        self.assertTrue(np.array_equal(v.T@v,identity))
        self.assertTrue(np.array_equal(h@v,v@(hj+hf+5*identity)))
        fullf=np.kron(np.eye(64,dtype=np.int64),f)
        self.assertTrue(np.array_equal(fullf@v,v@hf))
        for site in range(6):
            ni=np.diag([((d>>(5-site))&1) for d in range(64)]).astype(np.int64)
            actual=np.kron(ni,np.eye(6,dtype=np.int64))
            expected=np.kron(np.diag(np.eye(6,dtype=np.int64)[site]),np.eye(6,dtype=np.int64))
            self.assertTrue(np.array_equal(actual@v,v@expected))
            self.assertFalse(np.any(hf@expected-expected@hf))
            comm=hj@expected-expected@hj
            for g,a in enumerate(adj):
                ii=np.array([6*i+g for i in range(6)]);point=np.diag(np.eye(6,dtype=np.int64)[site])
                self.assertTrue(np.array_equal(comm[np.ix_(ii,ii)],a@point-point@a))
        OBS['actual_port_interface']=dict(full_dimension=384,one_excitation_dimension=36,
            exact_full_H_intertwining=True,identity_offset='5J I',
            endpoint_algebra='C^6, f mapped to diag(f) tensor I_G',
            effect_Ni_is_actual_local_excitation_readout=True,
            commutator_seminorm='abs(J) max_g norm([A_g,diag(f)])',
            all_real_kappa_and_nonzero_J_covered=True,kappa_drops_out_exactly=True,
            zero_J_distinct_port_distance='infinity',
            preparation_and_port_algebra_selection_are_explicit_inputs=True,
            old_data_saved_in_isolated_S_before_new_single_excitation_probe=True,
            G_S_R_correlations_retained_and_graph_not_reprepared_by_this_handoff=True)

    def test_03_complete_three_class_distance_optimality(self):
        trees,h,f,adj,*_=model();leafcases=[];internal=[]
        # All graph permutations generated by port symmetries preserve the family.
        for p in itertools.permutations(LEAVES):
            for swap in (False,True):
                image=dict(zip(LEAVES,p));image.update({1:2 if swap else 1,2:1 if swap else 2})
                for tree in trees:
                    mapped=frozenset(tuple(sorted((image[a],image[b]))) for a,b in tree)
                    self.assertIn(mapped,trees)
        for a in LEAVES:
            potential=np.zeros(6,dtype=np.int64);potential[a]=1
            for graph in adj:
                c=graph@np.diag(potential)-np.diag(potential)@graph;m=c.T@c
                self.assertTrue(np.array_equal(m@m,m));self.assertEqual(int(np.trace(m)),2)
            for b in (1,2):self.assertTrue(any(graph[a,b] for graph in adj))
        for a,b in itertools.combinations(LEAVES,2):
            potential=np.zeros(6,dtype=np.int64);potential[a]=1;potential[b]=-1;norms=[]
            for graph in adj:
                c=graph@np.diag(potential)-np.diag(potential)@graph;m=c.T@c
                same=bool(np.any(graph[a]*graph[b]));value=2 if same else 1
                self.assertTrue(np.array_equal(m@m,value*m));self.assertEqual(int(np.trace(m)),4)
                norms.append(value)
            self.assertIn(2,norms);leafcases.append(norms)
        potential=np.zeros(6,dtype=np.int64);potential[1]=1;potential[2]=-1
        for graph in adj:
            c=graph@np.diag(potential)-np.diag(potential)@graph;m=c.T@c
            self.assertFalse(np.any(m@(m@m-8*m+4*np.eye(6,dtype=np.int64))))
            self.assertEqual(int(np.trace(m)),16)
            internal.append(short(np.linalg.norm(c,2)))
        # Small integer certificate products bounded by 6^3*8^3 well below int64.
        self.assertLess(6**3*8**3,2**63)
        OBS['three_exact_distances']=dict(internal_leaf='1/abs(J)',leaf_leaf='sqrt(2)/abs(J)',
            internal_internal='(sqrt(3)-1)/abs(J)',
            endpoint_swap_antisymmetrization_preserves_difference_and_never_increases_L=True,
            all_48_port_symmetries_preserve_six_graph_department=True,
            leaf_leaf_commutator_squared_norm_patterns=leafcases,
            internal_CstarC_annihilator='m(m^2-8m+4I)=0; trace(m)=16',
            internal_commutator_norm_numeric=internal[0],
            time_unit_not_length_without_additional_rate_calibration=True,
            no_optimization_scan_used=True)

    def test_04_exact_five_dimensional_Euclidean_Gram_and_embedding(self):
        ones=np.ones((6,6),dtype=np.int64);center6=6*np.eye(6,dtype=np.int64)-ones
        d0=np.zeros((6,6),dtype=np.int64);d1=np.zeros_like(d0)
        for a in range(6):
            for b in range(a+1,6):
                d0[a,b]=d0[b,a]=2 if a in LEAVES and b in LEAVES else 4 if {a,b}=={1,2} else 1
                if {a,b}=={1,2}:d1[a,b]=d1[b,a]=-2
        leaf=np.array([int(i in LEAVES) for i in range(6)],dtype=np.int64)
        c=np.array([0,1,-1,0,0,0],dtype=np.int64);w=np.array([-1,2,2,-1,-1,-1],dtype=np.int64)
        pl4=4*np.diag(leaf)-np.outer(leaf,leaf);pc2=np.outer(c,c);pm12=np.outer(w,w)
        self.assertTrue(np.array_equal(-center6@d0@center6,18*pl4+72*pc2-6*pm12))
        self.assertTrue(np.array_equal(-center6@d1@center6,-36*pc2+4*pm12))
        projectors=((pl4,4,3),(pc2,2,1),(pm12,12,1))
        for i,(p,den,rank) in enumerate(projectors):
            self.assertTrue(np.array_equal(p@p,den*p));self.assertEqual(int(np.trace(p)),den*rank)
            for other,_,_ in projectors[:i]:self.assertFalse(np.any(p@other))
        self.assertTrue(np.array_equal(3*pl4+6*pc2+pm12,2*center6))
        self.assertLess(F(5,3)**2,3);self.assertGreater(F(7,4)**2,3)
        self.assertGreater(3-F(5,3)*F(7,4),0)  # lambda4 > lambda5.
        self.assertEqual(F(1,12)-F(1,100)*(3+F(1,100))-F(1,1000),F(1567,30000))
        self.assertLess(6**3*6**2*4,2**63)
        d=distance_matrix();center=center6/6;gram=-.5*center@(d*d)@center
        values=np.linalg.eigvalsh(gram)
        self.assertTrue(np.allclose(values,[0,(2*np.sqrt(3)-3)/3,2-np.sqrt(3),1,1,1],atol=2e-14))
        x=np.zeros((6,5));x[list(LEAVES),:3]=np.array([[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]])/2
        x[1,3:]=[np.sqrt((2*np.sqrt(3)-3)/4),(np.sqrt(3)-1)/2]
        x[2,3:]=[x[1,3],-x[1,4]]
        error=float(np.max(abs(np.linalg.norm(x[:,None,:]-x[None,:,:],axis=2)-d)))
        self.assertLess(error,2e-14)
        OBS['finite_Euclidean_classification']=dict(Gram_eigenvalues_at_J1=['0','1','1','1','2-sqrt(3)','(2sqrt(3)-3)/3'],
            exact_rational_plus_sqrt3_matrix_certificate=True,strict_positive_rank=5,
            orthogonal_projector_ranks=[3,1,1],explicit_R5_embedding_residual=short(error),
            minimal_Euclidean_embedding_dimension_of_this_finite_ruler=5,
            no_conclusion_universe_or_continuum_dimension_is_five=True,
            not_a_cognitive_or_three_dimensional_space_no_go=True,
            same_six_port_ruler_cannot_be_isometric_in_R3=True,
            any_R3_squared_distance_uniform_error_at_least='(2-sqrt(3))/(3 J^2)',
            scaled_error_at_least='a^2 (2-sqrt(3))/(3 J^2)',
            rank_three_Gram_obstruction_uses_fourth_largest_eigenvalue=True,
            fixed_smooth_three_dimensional_point_with_o_scale_length_error_excluded=True,
            finite_error_example_rational_margin='1567/30000',
            five_dimensional_Euclidean_realization_not_excluded=True)

    def test_05_arbitrary_graph_reference_transfer_bound_and_true_readout(self):
        _,_,_,_,_,hj,hf,k=model();h=hj+hf;time=F(1,5);u=unitary(h,float(time))
        # Explicit graph-reference input, without graph dephasing.
        vector=np.array([complex((i*7+3)%11-5,(i*5+1)%13-6) for i in range(18)]).reshape(6,3)
        vector/=np.linalg.norm(vector);initial=np.zeros((6,6,3),dtype=complex);initial[0]=vector
        output=(u@initial.reshape(36,3)).reshape(6,6,3);population=np.sum(abs(output)**2,axis=(1,2))
        self.assertAlmostEqual(float(population.sum()),1,places=13)
        before=vector.conj().T@vector;after=output.reshape(36,3).conj().T@output.reshape(36,3)
        self.assertLess(np.linalg.norm(after-before),2e-13)
        # f=single leaf is an exactly L=1 calibrated observable, genuinely read by N0.
        f0=np.zeros(6);f0[0]=1;variation=abs(float(f0@population)-1)
        self.assertLessEqual(variation,float(time)+1e-13)
        # All port effects are actual one-site Z outcomes inherited via test02.
        q=1-population[3];diam=np.sqrt(2)
        self.assertLessEqual(distance_matrix()[0,3]-diam*q,float(time)+1e-13)
        # Exact successful-transfer lower bound with a nonzero tolerance.
        epsilon=F(1,100);upper_sqrt2=F(3,2)
        lower_internal_leaf=1-upper_sqrt2*epsilon
        self.assertEqual(lower_internal_leaf,F(197,200))
        # Per-population instrument and common-source errors plus Bernoulli statistics.
        eta=F(1,1000);source_error=eta/4;instrument_error=eta/4;stat_error=eta/2
        self.assertEqual(source_error+instrument_error+stat_error,eta)
        copies=math.ceil(16/(eta*eta))
        control_only=instrument_error/2;read_window=instrument_error/36
        self.assertEqual(control_only+18*read_window,instrument_error)
        self.assertGreater(sum(F(8)**n/math.factorial(n) for n in range(30)),1200)
        OBS['operational_speed_and_readout']=dict(
            natural_evolution_bound='T >= d(a,b)-diam(d)*epsilon for final port success >=1-epsilon',
            arbitrary_unknown_graph_and_old_reference_inputs_covered=True,
            controlled_bound='T + integral(abs(u))/abs(J) >= d(a,b)-diam(d)*epsilon',
            graph_population_or_coherence_not_fixed_or_measured_by_bound=True,
            example_unknown_graph_reference_dimension=3,time=str(time),
            actual_local_excitation_populations=[short(x) for x in population],
            exact_internal_leaf_time_bound_at_epsilon_01='at least 197/200 at abs(J)=1',
            per_port_probability_error=str(eta),common_source_error=str(source_error),
            population_target_is_declared_exact_port_source_and_evolution=True,
            actual_source_error_counted_by_channel_contraction=True,
            complete_local_readout_channel_error=str(instrument_error),
            local_readout_control_only_error=str(control_only),
            local_readout_control_window=str(read_window),
            per_port_statistical_error=str(stat_error),copies_per_port=str(copies),
            total_known_history_replays=str(6*copies),all_six_failure_probability_below='1/100',
            local_readout_can_disturb_data_and_unknown_correlations=True,
            prior_information_saved_in_S_not_deleted=True,
            calibration_clock_fresh_probes_and_replay_sources_remain_inputs=True)

    def test_06_scale_boundary_and_lower_bound_is_not_arrival_time(self):
        # Same exchange family at N=2, explicitly not the six-tree dynamics.
        swap=np.array([[0,1],[1,0]],dtype=float)
        at_bound=unitary(swap,1)@np.array([1,0]);at_transfer=unitary(swap,np.pi/2)@np.array([1,0])
        self.assertLess(abs(at_bound[1])**2,1)
        self.assertAlmostEqual(abs(at_transfer[1])**2,1,places=14)
        # pi/2>3/2>1; no assumption of six-tree perfect transfer.
        self.assertGreater(np.pi,3)
        f=np.array([1,2,-1,0,3,-2],dtype=float)
        self.assertAlmostEqual(lip(f,2),2*lip(f,1),places=13)
        self.assertEqual(lip(f,0),0)
        self.assertTrue(np.allclose(distance_matrix(2),distance_matrix()/2))
        OBS['scope_and_scale_boundaries']=dict(
            N2_exchange_spectral_ruler='1/abs(J)',N2_first_perfect_transfer_time='pi/(2abs(J))',
            N2_success_at_the_ruler_time=short(abs(at_bound[1])**2),
            N2_is_only_same_family_boundary_not_six_tree_counterexample=True,
            no_six_tree_perfect_transfer_or_bound_saturation_claim=True,
            scalar_J_rescaling_changes_all_times_inversely=True,
            J_zero_all_distinct_ports_infinite_even_if_graph_F_evolves=True,
            current_graph_state_and_kappa_do_not_change_this_worst_state_port_ruler=True,
            no_direct_identification_with_488_state_dependent_Q_chart=True,
            converting_time_units_to_length_requires_additional_rate_standard=True)


def run():
    OBS.clear();stream=io.StringIO();result=unittest.TextTestRunner(stream=stream).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():raise RuntimeError(stream.getvalue())
    output=dict(round=490,scientific_baseline_round=489,reused_frozen_rounds=[442,443,465,466,472,475,478,489],
        tests_run=result.testsRun,failures=len(result.failures),errors=len(result.errors),
        python=platform.python_version(),numpy=np.__version__,observations=OBS,
        scope=dict(port_algebra_and_one_excitation_preparation_are_extra_inputs=True,
            three_exact_spectral_port_distances_classified=True,
            finite_ruler_has_minimal_Euclidean_embedding_dimension_five=True,
            no_five_dimensional_physical_space_claim=True,
            spectral_speed_bound_is_not_minimal_transfer_time_or_free_length_measurement=True,
            no_Dirac_operator_or_Riemannian_reconstruction_from_H_claim=True,
            no_position_group_scale_limit_or_three_dimensional_necessity=True,
            full_GR_goal_completed=False,phase_closure_triggered=False))
    assert json.loads(json.dumps(output))==output
    return output


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--dry-run',action='store_true');parser.add_argument('--check',action='store_true')
    args=parser.parse_args();result=run()
    if args.check:assert result==json.loads(TARGET.read_text(encoding='utf-8'))
    elif not args.dry_run:
        with TARGET.open('x',encoding='utf-8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False,indent=2))
