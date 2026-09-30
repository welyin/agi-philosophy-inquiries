"""Round 482: composite full-SWAP reference closure, scientific baseline 475.

Six original graph vertices now carry C4, without compressing private factors.
Exact evidence uses S6 group algebra x integer graph matrices, not floating ranks.
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
import rigid_leaf_reference_audit as old
import sequential_tree_distance_readout_audit as reader

TARGET=Path(__file__).with_name('composite_reference_closure_audit_results.json')
OBS={}
ORDER=10

# Gaussian arithmetic: pairs of genuine Python integers, never Python complex.
def gm(a,b):
    return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def ga(a,b):
    return a[0]+b[0],a[1]+b[1]
GI=(((1,0),(0,0)),((0,0),(1,0)))
GX=(((0,0),(1,0)),((1,0),(0,0)))
GY=(((0,0),(0,-1)),((0,1),(0,0)))
GZ=(((1,0),(0,0)),((0,0),(-1,0)))
GP=(GX,GY,GZ)
def madd(a,b):
    return tuple(tuple(ga(a[i][j],b[i][j]) for j in range(2)) for i in range(2))
def mm(a,b):
    return tuple(tuple(ga(gm(a[i][0],b[0][j]),gm(a[i][1],b[1][j])) for j in range(2)) for i in range(2))
def trace(a):
    return ga(a[0][0],a[1][1])
def local_numerators(axis):
    return (GI,GI if axis<0 else madd(GI,GP[axis]),GI,
            madd(GI,GX),madd(GI,GY),madd(GI,GZ))
def permutation_trace(p,axis):
    factors=local_numerators(axis)
    seen=set(); result=(1,0)
    for start in range(6):
        if start in seen: continue
        cycle=[];j=start
        while j not in seen:
            seen.add(j);cycle.append(j);j=p[j]
        product=GI
        for j in reversed(cycle): product=mm(product,factors[j])
        result=gm(result,trace(product))
    assert all(type(v) is int for v in result)
    return result

@lru_cache(None)
def algebra():
    permutations=list(itertools.permutations(range(6)))
    index={p:i for i,p in enumerate(permutations)}
    trees,_,fg=old.system()
    edges=sorted(set().union(*trees))
    masks=[np.array([int(e in g) for g in trees],dtype=np.int64) for e in edges]
    left=[];right=[]
    for a,b in edges:
        swap=list(range(6));swap[a],swap[b]=swap[b],swap[a]
        left.append(np.array([index[tuple(swap[p[k]] for k in range(6))] for p in permutations]))
        right.append(np.array([index[tuple(p[swap[k]] for k in range(6))] for p in permutations]))
    return permutations,index,edges,masks,left,right,fg

def ad_h(z):
    _,_,_,masks,left,right,fg=algebra()
    result=fg@z-z@fg
    for n,l,r in zip(masks,left,right):
        result[l]+=n[None,:,None]*z
        result[r]-=z*n[None,None,:]
    return result

@lru_cache(None)
def moments(leaf=3, axis_a=0, axis_b=0, layers=2):
    permutations,index,*_=algebra()
    weight=[]
    for p in permutations:
        value=permutation_trace(p,axis_a)
        if layers==2: value=gm(value,permutation_trace(p,axis_b))
        weight.append(value)
    # A regular S6 representation of H has absolute row sum 9.  Thus each
    # ad^k entry is <=2*18^k; even a 36-entry graph sum fits signed int64.
    assert 72*18**ORDER < 2**63
    z=np.zeros((720,6,6),dtype=np.int64)
    z[index[tuple(range(6))]]=np.diag(3*old.distance(0,leaf)-8)
    coefficients=[]; numerators=[]; peak=[]
    denominator=(64**layers)*6*3
    for k in range(ORDER+1):
        graph_sums=z.sum(axis=(1,2))
        total=(0,0)
        for j,value in enumerate(graph_sums):
            total=ga(total,(int(value)*weight[j][0],int(value)*weight[j][1]))
        rotated=gm(((1,0),(0,1),(-1,0),(0,-1))[k%4],total)
        assert type(rotated[0]) is int and type(rotated[1]) is int
        assert rotated[1]==0
        numerators.append(rotated[0])
        coefficients.append(F(rotated[0],denominator*math.factorial(k)))
        peak.append(int(np.max(np.abs(z))))
        if k<ORDER: z=ad_h(z)
    return coefficients,numerators,peak

def tail(t):
    x=18*t
    assert x<12
    return x**11/(2*math.factorial(11)*(1-x/12))

def finite_window():
    a=F(1,64);b=F(1,32)
    c=moments()[0]
    negative=sum(-v*b**(k-4) for k,v in enumerate(c) if k>=5 and v<0)
    positive=sum(v*b**(k-4) for k,v in enumerate(c) if k>=5 and v>0)
    # tail(t)/t^4 increases on this interval and down to t=0.
    remainder=tail(b)/b**4
    lower=c[4]-negative-remainder
    upper=c[4]+positive+remainder
    assert lower>F(7,100) and upper<F(2,25)
    return a,b,lower,upper,F(7,100)*a**4

def macro_swap(d):
    p=np.zeros((d*d,d*d),dtype=np.int64)
    for a in range(d):
        for b in range(d):p[b*d+a,a*d+b]=1
    return p

def gaussian_tensor(factors):
    real=np.array([[1]],dtype=np.int64);imag=np.array([[0]],dtype=np.int64)
    for factor in factors:
        ar=np.array([[v[0] for v in row] for row in factor],dtype=np.int64)
        ai=np.array([[v[1] for v in row] for row in factor],dtype=np.int64)
        real,imag=np.kron(real,ar)-np.kron(imag,ai),np.kron(real,ai)+np.kron(imag,ar)
    return real,imag

@lru_cache(None)
def full_d4_diagnostic():
    # Independent full 4^6 x 6 Hilbert-space action.  The 16 columns are only a
    # decomposition of the specified mixed preparation, not a private code.
    _,_,edges,masks,_,_,fg=algebra()
    labels=np.arange(4**6)
    digit=[labels//(4**(5-j))%4 for j in range(6)]
    swaps=[labels+(digit[a]-digit[b])*(4**(5-b)-4**(5-a)) for a,b in edges]
    def act(v):
        result=np.einsum('gh,dhi->dgi',fg,v)
        for permutation,n in zip(swaps,masks): result+=v[permutation]*n[None,:,None]
        return result
    # Independent Gaussian-integer pure columns and double-sided powers H^k.
    one=np.ones(4,dtype=np.int64);zero=np.zeros(4,dtype=np.int64)
    yy_r=np.array([1,0,0,-1],dtype=np.int64);yy_i=np.array([0,1,1,0],dtype=np.int64)
    zz=np.array([1,0,0,0],dtype=np.int64);eye4=np.eye(4,dtype=np.int64)
    column_r=[];column_i=[]
    for a,b in itertools.product(range(4),repeat=2):
        rr=np.array([1],dtype=np.int64);ii=np.array([0],dtype=np.int64)
        for ar,ai in ((eye4[a],zero),(one,zero),(eye4[b],zero),
                      (one,zero),(yy_r,yy_i),(zz,zero)):
            rr,ii=np.kron(rr,ar)-np.kron(ii,ai),np.kron(rr,ai)+np.kron(ii,ar)
        column_r.append(rr);column_i.append(ii)
    raw_r=np.repeat(np.stack(column_r,axis=1)[:,None,:],6,axis=1)
    raw_i=np.repeat(np.stack(column_i,axis=1)[:,None,:],6,axis=1)
    powers=[(raw_r,raw_i)]
    for k in range(4):powers.append((act(powers[-1][0]),act(powers[-1][1])))
    assert 4*24576*16*9**4 < 2**63
    o=(3*old.distance(0,3)-8)[None,:,None]
    independent_coeff=[]
    for k in range(5):
        value=F(0)
        for a in range(k+1):
            b=k-a;ar,ai=powers[a];br,bi=powers[b]
            real=int(np.sum(o*(ar*br+ai*bi)))
            imag=int(np.sum(o*(ar*bi-ai*br)))
            phase=((1,0),(0,1),(-1,0),(0,-1))[(a-b)%4]
            rotated=gm(phase,(real,imag))
            value+=F(rotated[0],3*6144*math.factorial(a)*math.factorial(b))
        independent_coeff.append(value)
    assert independent_coeff==moments()[0][:5]
    x=np.ones(4,dtype=complex);y=np.array([1,1j,1j,-1]);z=np.array([1,0,0,0],complex)
    columns=[]
    for a,b in itertools.product(range(4),repeat=2):
        value=np.array([1],complex)
        for factor in (np.eye(4)[a],x,np.eye(4)[b],x,y,z): value=np.kron(value,factor)
        columns.append(value)
    initial=np.stack(columns,axis=1)[:,None,:]*np.ones((1,6,1))/math.sqrt(6144)
    assert abs(np.linalg.norm(initial)**2-1)<1e-13
    observations=[];final_graph=None
    coeff=moments()[0]
    for t in (F(1,64),F(1,32)):
        term=initial.copy();out=initial.copy()
        for k in range(1,19):
            term=(-1j*float(t)/k)*act(term);out+=term
        graph=np.einsum('dgi,dhi->gh',out,out.conj())
        drift=float(np.diag(graph).real@old.distance(0,3)-8/3)
        polynomial=float(sum(coeff[k]*t**k for k in range(11)))
        assert abs(drift-polynomial)<float(tail(t))+5e-13
        xnorm=9*t
        unitary_tail=xnorm**19/(math.factorial(19)*(1-xnorm/20))
        observations.append(dict(time=str(t),distance_drift=drift,
            exact_polynomial=polynomial,absolute_difference=abs(drift-polynomial),
            exponent_tail=str(unitary_tail),state_norm_error=abs(np.trace(graph).real-1)))
        final_graph=graph
    return observations,final_graph,independent_coeff

class Audit(unittest.TestCase):
    def test_01_full_macro_exchange_and_exact_probe_embedding(self):
        s4=macro_swap(4);s2=macro_swap(2)
        self.assertTrue(np.array_equal(s4@s4,np.eye(16,dtype=np.int64)))
        # Fully arbitrary local density is covered analytically by tensor swap.
        c=np.array([[1,1j],[0,1],[1,-1j],[1j,0]],complex)
        tau=c@c.conj().T
        self.assertTrue(np.array_equal(s4@np.kron(tau,tau),np.kron(tau,tau)@s4))
        v=np.zeros((4,2),dtype=np.int64);v[0,0]=v[2,1]=1
        vv=np.kron(v,v)
        self.assertTrue(np.array_equal(s4@vv,vv@s2))
        # A retained seven-dimensional unknown symmetric B register is a
        # positive range control: every B transposition is exactly identity.
        dicke=np.zeros((64,7),dtype=np.int64)
        for x in range(64):dicke[x,x.bit_count()]=1
        for a,b in itertools.combinations(range(6),2):
            image=[]
            for x in range(64):
                da=(x>>(5-a))&1;db=(x>>(5-b))&1
                image.append(x+(da-db)*((1<<(5-b))-(1<<(5-a))))
            self.assertTrue(np.array_equal(dicke[image],dicke))
        permutations,index,edges,masks,left,right,fg=algebra()
        for l,r in zip(left,right):
            self.assertEqual(len(set(l.tolist())),720)
            self.assertEqual(len(set(r.tolist())),720)
        self.assertTrue(np.array_equal(fg.sum(axis=1),np.full(6,4)))
        self.assertTrue(np.array_equal(sum(masks),np.full(6,5)))
        OBS['model_contract']=dict(local_dimension=4,full_data_dimension=4096,
            graph_dimension=6,total_dimension=24576,H_norm_bound=9,
            generator='F + sum full_SWAP4_ab tensor n_ab',
            full_SWAP4_is_product_of_two_layer_swaps=True,
            not_sum_of_micro_exchange_generators=True,
            same_existing_G_not_merger_of_two_independent_graphs=True,
            no_unknown_private_factor_compression=True,
            collective_U4_and_leaf_S4_covariance=True,
            macro_product_Hamiltonian_is_explicit_model_input=True,
            symmetric_B_sector_dimension=7,
            symmetric_B_unknown_reference_retained_as_identity=True,
            symmetric_B_sector_is_extra_preparation_not_all_private_states=True,
            probe_V='|a> maps to |a>_A |0>_B',exact_probe_intertwiner=True)

    def test_02_permutation_trace_and_integer_safety(self):
        permutations,*_=algebra()
        for axis in (-1,0,1,2):
            real,imag=gaussian_tensor(local_numerators(axis))
            for p in permutations:
                image=[]
                for x in range(64):
                    image.append(sum(((x>>(5-j))&1)<<(5-p[j]) for j in range(6)))
                direct=(sum(int(real[j,image[j]]) for j in range(64)),
                        sum(int(imag[j,image[j]]) for j in range(64)))
                self.assertEqual(direct,permutation_trace(p,axis))
        coeff,nums,peak=moments()
        self.assertTrue(all(type(v) is int for v in nums))
        self.assertLess(72*18**ORDER,2**63)
        OBS['integer_certificate']=dict(group='S6',permutations=720,
            graph_block_shape=[6,6],local_trace_checks=4*720,
            cycle_order='reverse each cycle of old-to-new permutation',
            exact_Gaussian_arithmetic='pairs of Python integers',
            graph_recursion_int64_prior_bound=str(72*18**ORDER),
            signed_int64_limit=str(2**63),maximum_block_entries=peak,
            expectation_denominator_before_factorial=73728,
            exact_rotated_trace_numerators=nums,
            full_SWAP4_expectation_is_product_of_layer_cycle_traces=True)

    def test_03_true_first_nonzero_reference_drift(self):
        coeff,_,_=moments()
        expected=[F(0)]*11
        for k,value in {4:F(43,576),6:F(-7801,34560),8:F(1140511,3870720),10:F(-78899731,348364800)}.items():expected[k]=value
        self.assertEqual(coeff,expected)
        self.assertEqual(moments(layers=1)[0],[F(0)]*11)
        others={}
        for leaf in (4,5):
            c,_,_=moments(leaf=leaf)
            self.assertEqual(c[4],F(-43,1152))
            self.assertEqual(abs(c[5]),F(17,960))
            others[str(leaf)]={str(k):str(c[k]) for k in range(11) if c[k]}
        OBS['reference_counterexample']=dict(
            source='|+X>_A tensor |+X>_B at vertex 1',
            leaves={'0':'I4/4','3':'|+X,+X>','4':'|+Y,+Y>','5':'|0,0>'},
            internal_vertex_2='I4/4',graph='uniform coherent |s_G>',
            target='mean(D03)-8/3',coefficients=[str(v) for v in coeff],
            first_nonzero_order=4,fourth_derivative='43/24',
            other_two_leaf_distances=others,
            d2_control_zero_through_order_10=True,
            d2_all_time_protection_reused_from_475=True,
            no_generic_time_evenness_claim=True,
            witness_changes_whole_network_path_distance=True)

    def test_04_uniform_finite_time_theorem(self):
        a,b,lower,upper,gap=finite_window()
        OBS['finite_window']=dict(
            entire_punctured_interval=['0','1/32'],
            drift_lower='7*t^4/100',drift_upper='2*t^4/25',
            certificate_lower_ratio=str(lower),certificate_upper_ratio=str(upper),
            uniform_window=[str(a),str(b)],uniform_drift_gap=str(gap),
            centered_distance_norm='1/2',H_norm_bound=9,
            remainder='(18*t)^11/(2*11!*(1-18*t/12))',
            graph_reference_distribution_does_not_remain_fixed=True,
            whole_DG_and_untouched_R_information_preserved_by_unitarity=True,
            reference_constancy_not_claimed_for_arbitrary_initial_R_correlations=True)

    def test_05_independent_full_Hilbert_space_and_actual_instrument(self):
        observations,graph,independent_coeff=full_d4_diagnostic()
        for row in observations:
            t=F(row['time'])
            self.assertGreater(row['distance_drift'],float(F(7,100)*t**4)-1e-13)
            self.assertLess(row['distance_drift'],float(F(2,25)*t**4)+1e-13)
            self.assertLess(row['state_norm_error'],1e-12)
        # Original 472 uses actual exp(-i H2 v) Kraus maps; exact local V lifts
        # the entire instrument to the d4 generator, not just a mean formula.
        probe=F(1,65536)
        minus,plus,_=reader.instrument((1,0),probe)
        probabilities=[float(np.trace(np.einsum('ghij,ij->gh',op,graph)).real) for op in (minus,plus)]
        self.assertLess(abs(sum(probabilities)-1),1e-12)
        self.assertGreaterEqual(min(probabilities),0)
        trees=old.system()[0]
        true=float(np.diag(graph).real@np.array([int((0,1) in t) for t in trees]))
        estimate=(probabilities[1]-probabilities[0])/float(probe)
        bias=(18*probe)**2/(2*(1-6*probe)*probe)
        self.assertLess(abs(estimate-true),float(bias)+2e-9)
        pairs=[];products=[]
        for e1,e2 in (((1,0),(1,3)),((2,0),(2,3))):
            m1,p1,_=reader.instrument(e1,probe);m2,p2,_=reader.instrument(e2,probe)
            signed=np.einsum('ghij,ij->gh',p1-m1,graph)
            signed=np.einsum('ghij,ij->gh',p2-m2,signed)
            estimate_product=float(np.trace(signed).real)/float(probe**2)
            products.append(estimate_product)
            pairs.append(dict(ordered_edges=[list(e1),list(e2)],product_estimate=estimate_product))
        measured_distance=3-sum(products)
        true_distance=float(np.diag(graph).real@old.distance(0,3))
        relative=(18*probe)**2/(2*(1-6*probe)*probe)
        distance_bias=2*(2*relative+relative**2)
        self.assertLess(abs(measured_distance-true_distance),float(distance_bias)+2e-6)
        OBS['independent_numerics_and_CP_readout']=dict(
            full_24576_dimensional_H_action=True,mixed_pure_columns=16,
            independent_integer_coefficients_through_order4=[str(v) for v in independent_coeff],
            full_integer_density_denominator=6144,
            full_integer_double_sided_bound=str(4*24576*16*9**4),
            U_polynomial_degree=18,observations=observations,
            actual_472_instrument_lifted_by_exact_V=True,
            ordered_probe_edge=[1,0],measured_port=1,probe_wait=str(probe),
            probabilities=probabilities,true_n01=true,estimated_n01=estimate,
            diagnostic_bias_bound=str(bias),
            actual_two_stage_path_records=pairs,
            estimated_whole_graph_D03=measured_distance,true_whole_graph_D03=true_distance,
            two_path_distance_bias_bound=str(distance_bias),
            graph_coherence_not_dephased=True,
            old_full_C4_data_relocated_to_internal_storage=True,
            G_stays_active_and_final_CP_measurement_may_disturb_G=True)

    def test_06_complete_finite_distance_readout_budget(self):
        _,_,_,_,gap=finite_window()
        # On six graphs n03=0; D03=3 - (n10*n13+n20*n23).
        trees=old.system()[0]
        masks={e:np.array([int(e in t) for t in trees]) for e in sorted(set().union(*trees))}
        self.assertTrue(np.array_equal(old.distance(0,3),3-masks[(0,1)]*masks[(1,3)]-masks[(0,2)]*masks[(2,3)]))
        epsilon=gap/16
        # Two signed outcomes w1*w2 / v^2 estimate each product.  Sequential
        # Jordan ideal has the right traced product; it is not postselection.
        v=epsilon/F(1024)
        a=18*v
        rem=a*a/(2*(1-a/3))
        gamma=epsilon*v*v/32
        relative=rem/v
        physical_error=2*relative+relative**2
        implementation_error=2*gamma/(v*v)
        self.assertLess(physical_error+implementation_error,epsilon/2)
        source=gap/4
        copies=F(32,1)/(epsilon*epsilon*v**4)
        copies=(copies.numerator+copies.denominator-1)//copies.denominator
        self.assertGreaterEqual(copies*epsilon*epsilon*v**4/8,4)
        # Four means across two settings -> 8 exp(-8), so use exponent >=8.
        copies*=2
        self.assertGreaterEqual(copies*epsilon*epsilon*v**4/8,8)
        self.assertGreater(sum(F(8**j,math.factorial(j)) for j in range(14)),800)
        # Two settings: two product means each; source distance centered norm .5.
        remaining=gap-source-4*epsilon
        self.assertEqual(remaining,gap/2)
        OBS['resource_ledger']=dict(
            distance_identity='D03=3-n01*n13-n02*n23',
            settings=['prepared initial reference','after t=1/32'],
            each_setting_complete_source_and_control_diamond_budget=str(source),
            each_product_mean_total_error=str(epsilon),probe_wait=str(v),
            each_full_stage_instrument_diamond_budget=str(gamma),
            stages_per_trial=2,product_means_per_setting=2,
            independent_copies_per_product_mean=str(copies),total_copies=str(4*copies),
            joint_failure_probability_less_than='1/100',
            retained_observed_drift=str(remaining),
            initial_and_final_reference_compared_with_same_whole_G_readout=True,
            signed_records_not_free_joint_graph_measurement=True,
            source_waiting_preparation_swap_storage_and_instrument_errors_counted=True,
            old_record_feedback_excluded=True,
            same_actual_prepared_state_per_setting_used_for_both_product_means=True,
            unknown_single_seed_not_cloned=True,huge_sampling_not_executed=True)

def run():
    OBS.clear();stream=io.StringIO()
    result=unittest.TextTestRunner(stream=stream).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():raise RuntimeError(stream.getvalue())
    output=dict(round=482,scientific_baseline_round=475,tests_run=result.testsRun,
        failures=len(result.failures),errors=len(result.errors),python=platform.python_version(),
        numpy=np.__version__,observations=OBS,
        scope=dict(specified_macro_full_SWAP_reference_lift_falsified=True,
            not_a_refutation_of_all_composite_reference_architectures=True,
            no_universal_reference_statistics_protection_for_unknown_private_states=True,
            not_two_independently_positioned_sources=True,
            no_spatial_dimension_or_cognitive_axiom_no_go=True,
            full_GR_goal_completed=False,phase_closure_triggered=False))
    assert json.loads(json.dumps(output))==output
    return output

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dry-run',action='store_true');parser.add_argument('--check',action='store_true')
    args=parser.parse_args();result=run()
    if args.check:assert result==json.loads(TARGET.read_text(encoding='utf-8'))
    elif not args.dry_run:
        with TARGET.open('x',encoding='utf-8',newline='\n') as out:out.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False,indent=2))
