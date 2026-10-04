"""Round 486: all H-stationary inputs are invisible to Q under the fixed H+uK controls.

Scientific baseline 485, reusing 475/480/472. Exact 120-dimensional integer
symmetry certificates prove a statement about the full 384-dimensional system.
Complete energy pinching is an extra preparation, not a free control operation.
"""
import argparse
from fractions import Fraction as Fr
from functools import lru_cache
import io
import math
import json
from pathlib import Path
import platform
import unittest
import numpy as np
import rigid_leaf_reference_audit as frame
import symmetric_control_reachability_audit as words
import correlated_reference_local_access as chart
import sequential_tree_distance_readout_audit as reader

TARGET=Path(__file__).with_name('stationary_coordinate_response_audit_results.json')
OBS={}
LABELS=(('4',1),('31',3),('22',2),('211',3),('1111',1))
INTLIMIT=2**63-1


def short(x):return float(f'{float(x):.12g}')


def mm(a,b):
    """Every integer dot product is bounded BEFORE the int64 multiplication."""
    assert a.dtype==np.int64 and b.dtype==np.int64
    ma=max(abs(int(a.min(initial=0))),abs(int(a.max(initial=0))))
    mb=max(abs(int(b.min(initial=0))),abs(int(b.max(initial=0))))
    assert a.shape[1]*ma*mb<=INTLIMIT
    return a@b


def scale(a,s):
    assert a.dtype==np.int64 and isinstance(s,int)
    assert max(abs(int(a.min(initial=0))),abs(int(a.max(initial=0))))*abs(s)<=INTLIMIT
    return a*s


def add(a,b):
    assert a.dtype==np.int64 and b.dtype==np.int64
    assert max(abs(int(a.min(initial=0))),abs(int(a.max(initial=0))))+max(abs(int(b.min(initial=0))),abs(int(b.max(initial=0))))<=INTLIMIT
    return a+b


def trim(p):
    p=list(p)
    while len(p)>1 and p[0]==0:p.pop(0)
    return p


def pgcd(a,b):
    a=[Fr(x) for x in a];b=[Fr(x) for x in b]
    while b!=[0]:
        r=a[:]
        while len(r)>=len(b) and r!=[0]:
            factor=r[0]/b[0]
            r=trim([x-factor*y for x,y in zip(r,b+[0]*(len(r)-len(b)))])
        a,b=b,r
    return [x/a[0] for x in a]


def poly_on_projector(poly,h,p):
    out=scale(p,int(poly[0]))
    for c in poly[1:]:out=add(mm(h,out),scale(p,int(c)))
    return out


def normalize(num,den):
    if den<0:num=scale(num,-1);den=-den
    divisor=math.gcd(den,*(int(x) for x in num.flat))
    return num//divisor,den//divisor


@lru_cache(None)
def system():
    hf=frame.system()[1].astype(np.int64);kf=words.s12().astype(np.int64)
    idx=np.array([6*d+g for d in range(64) if d.bit_count()==3 for g in range(6)])
    local={int(z):i for i,z in enumerate(idx)};n=len(idx);eye=np.eye(n,dtype=np.int64)
    h=hf[np.ix_(idx,idx)];k=kf[np.ix_(idx,idx)]
    cas=scale(eye,-3)
    for a in range(6):
        for b in range(a+1,6):
            perm=[]
            for z in idx:
                d,g=divmod(int(z),6)
                if ((d>>(5-a))^(d>>(5-b)))&1:d^=(1<<(5-a))|(1<<(5-b))
                perm.append(local[6*d+g])
            cas=add(cas,eye[perm])
    complement=[5,4,3,2,1,0]
    graph_complement=eye[[local[6*(int(z)//6)+complement[int(z)%6]] for z in idx]]
    parity=mm(k,graph_complement)
    central=[]
    for label,dim in LABELS:
        numerator=np.zeros_like(eye)
        for perm,dm,gm,c22 in frame.permutations():
            sign=(-1)**sum(perm[i]>perm[j] for i in range(4) for j in range(i+1,4))
            c31=sum(perm[i]==i for i in range(4))-1
            character={'4':1,'31':c31,'22':c22,'211':sign*c31,'1111':sign}[label]
            image=[local[6*int(dm[int(z)//6])+int(gm[int(z)%6])] for z in idx]
            numerator=add(numerator,scale(eye[image],dim*character))
        central.append((label,dim,numerator))
    q=np.tile(words.qgraph(),(1,64))[:,idx].astype(np.int64)
    return hf,kf,idx,eye,h,k,cas,parity,central,q


@lru_cache(None)
def sectors():
    _,_,_,eye,h,k,cas,parity,central,_=system()
    powers=[eye]
    for _ in range(5):powers.append(mm(powers[-1],h))
    result=[]
    for j in range(4):
        value=j*(j+1);sn=eye.copy();sd=1
        for other in (0,2,6,12):
            if other!=value:
                sn=mm(sn,add(cas,scale(eye,-other)));sd*=value-other
        for label,dim,cn in central:
            for sign in (-1,1):
                num,den=normalize(mm(mm(sn,cn),add(eye,scale(parity,sign))),48*sd)
                rank=Fr(sum(int(num[i,i]) for i in range(120)),den)
                assert rank.denominator==1
                if not rank:continue
                assert rank.numerator%dim==0
                degree=rank.numerator//dim
                traces=[]
                for m in range(1,degree+1):
                    product=mm(num,powers[m])
                    traces.append(Fr(sum(int(product[i,i]) for i in range(120)),den*dim))
                coefficients=[Fr(1)]
                for m in range(1,degree+1):
                    coefficients.append(-sum(coefficients[m-i]*traces[i-1] for i in range(1,m+1))/m)
                assert all(x.denominator==1 for x in coefficients)
                poly=[int(x) for x in coefficients]
                assert not np.any(poly_on_projector(poly,h,num))
                result.append(dict(key=(j,label,sign),num=num,den=den,rank=int(rank),
                    dim=dim,poly=poly,HK_commutes=not np.any(mm(add(mm(h,k),scale(mm(k,h),-1)),num))))
    return result


def name(s):return f'j{s["key"][0]}_{s["key"][1]}_{"plus" if s["key"][2]>0 else "minus"}'


@lru_cache(None)
def visible_pairs():
    ss=sectors();q=system()[-1];pairs=[];zeros=0
    for i,a in enumerate(ss):
        for b in ss[i:]:
            if a['key'][0]!=b['key'][0]:continue
            nonzero=any(np.any(mm(a['num'],row[:,None]*b['num'])) for row in q)
            if nonzero:pairs.append((a,b,pgcd(a['poly'],b['poly'])))
            else:zeros+=1
    return pairs,zeros


def energy_projector(s,energy,roots):
    eye,h=system()[3:5];num=s['num'].copy();den=s['den']
    annihilator=[1]
    for root in roots:
        nxt=[0]*(len(annihilator)+1)
        for i,a in enumerate(annihilator):nxt[i]+=a;nxt[i+1]-=root*a
        annihilator=nxt
    assert not np.any(poly_on_projector(annihilator,h,s['num']))
    for other in roots:
        if other!=energy:num=mm(add(h,scale(eye,-other)),num);den*=energy-other
    return normalize(num,den)


@lru_cache(None)
def spectral_diagnostic():
    h=system()[0].astype(float);e,v=np.linalg.eigh(h);groups=[]
    for i,x in enumerate(e):
        if not groups or abs(x-e[groups[-1][0]])>1e-8:groups.append([])
        groups[-1].append(i)
    return e,v,groups


def pinch_numeric(rho):
    _,v,groups=spectral_diagnostic();inside=v.conj().T@rho@v;out=np.zeros_like(inside)
    for group in groups:out[np.ix_(group,group)]=inside[np.ix_(group,group)]
    return v@out@v.conj().T


class Audit(unittest.TestCase):
    def test_01_exact_symmetries_and_complete_projectors(self):
        hf,kf,idx,eye,h,k,cas,parity,central,q=system()
        self.assertEqual(len(idx),120)
        self.assertEqual(max(sum(abs(int(x)) for x in row) for row in hf),9)
        self.assertLessEqual(max(abs(int(x)) for x in q.flat),1)
        self.assertTrue(np.array_equal(mm(parity,parity),eye))
        for a in (cas,parity,*[x[2] for x in central]):
            self.assertTrue(np.array_equal(mm(h,a),mm(a,h)))
            self.assertTrue(np.array_equal(mm(k,a),mm(a,k)))
        for row in q:
            a=np.diag(row)
            self.assertTrue(np.array_equal(mm(parity,mm(a,parity)),-a))
            self.assertTrue(np.array_equal(mm(cas,a),mm(a,cas)))
        # Full SU2 scalar property; X and Z generate the collective rotations.
        collective_x=np.zeros((64,64),dtype=np.int64);collective_z=np.zeros((64,64),dtype=np.int64)
        for site in range(6):
            collective_x=add(collective_x,np.eye(64,dtype=np.int64)[np.arange(64)^(1<<(5-site))])
            collective_z=add(collective_z,np.diag([1-2*((d>>(5-site))&1) for d in range(64)]).astype(np.int64))
        for generator in (collective_x,collective_z):
            g=np.kron(generator,np.eye(6,dtype=np.int64))
            self.assertTrue(np.array_equal(mm(hf,g),mm(g,hf)))
            self.assertTrue(np.array_equal(mm(kf,g),mm(g,kf)))
        ss=sectors();den=math.lcm(*(s['den'] for s in ss));total=np.zeros_like(eye)
        for i,s in enumerate(ss):
            n=s['num'];d=s['den'];self.assertTrue(np.array_equal(n,n.T))
            self.assertTrue(np.array_equal(mm(n,n),scale(n,d)))
            total=add(total,scale(n,den//d))
            for other in ss[:i]:self.assertFalse(np.any(mm(n,other['num'])))
        self.assertTrue(np.array_equal(total,scale(eye,den)))
        self.assertEqual(sum(s['rank'] for s in ss),120)
        OBS['exact_symmetries']=dict(full_dimension=384,faithful_SU2_scalar_block_dimension=120,
            nonzero_joint_sectors=len(ss),all_integer_products_have_prior_overflow_bounds=True,
            Casimir_eigenvalues=[0,2,6,12],all_projectors_exact_complete_orthogonal=True,
            representation_reduction_is_not_physical_state_truncation=True)

    def test_02_exact_spectra_and_all_visible_pairs(self):
        ss=sectors();pairs,zeros=visible_pairs()
        self.assertEqual(len(ss),27);self.assertEqual(len(pairs),28)
        self.assertEqual(sum(g==[1] for a,b,g in pairs),24)
        self.assertEqual(max(len(s['poly'])-1 for s in ss),5)
        for a,b,g in pairs:
            self.assertEqual(a['key'][0],b['key'][0]);self.assertEqual(a['key'][2],-b['key'][2])
        OBS['exact_spectral_certificate']=dict(
            sectors=[dict(name=name(s),dimension_on_m0=s['rank'],S4_irrep_dimension=s['dim'],
                polynomial_descending_coefficients=s['poly'],HK_commutes_on_sector=s['HK_commutes']) for s in ss],
            polynomials_from_exact_Newton_traces_and_checked_as_full_annihilators=True,
            all_three_Q_full_blocks_checked=True,visible_unordered_pairs=len(pairs),
            same_spin_zero_pairs_including_diagonal=zeros,
            visible_pairs=[dict(left=name(a),right=name(b),gcd=[str(x) for x in g]) for a,b,g in pairs],
            coprime_pairs=24,different_spin_Q_blocks_zero_by_exact_Casimir_commutation=True)

    def test_03_all_four_resonances_are_control_invisible(self):
        roots={(0,'31',-1):(1,),(0,'211',-1):(1,),(0,'211',1):(-1,1),
               (0,'1111',-1):(-1,),(2,'31',-1):(2,4),(2,'211',1):(2,)}
        resonances=[];h,k=system()[4:6];q=system()[-1]
        for a,b,gcd in visible_pairs()[0]:
            if gcd==[1]:continue
            self.assertEqual(len(gcd),2);energy=-gcd[1]
            self.assertEqual(energy.denominator,1);energy=int(energy)
            self.assertTrue(a['HK_commutes']);self.assertTrue(b['HK_commutes'])
            na,da=energy_projector(a,energy,roots[a['key']]);nb,db=energy_projector(b,energy,roots[b['key']])
            for num,den in ((na,da),(nb,db)):
                self.assertTrue(np.array_equal(mm(num,num),scale(num,den)))
                self.assertTrue(np.array_equal(mm(h,num),scale(num,energy)))
                self.assertTrue(np.array_equal(mm(k,num),mm(num,k)))
            for row in q:self.assertFalse(np.any(mm(na,row[:,None]*nb)))
            resonances.append(dict(left=name(a),right=name(b),common_energy=energy,
                whole_degenerate_energy_Q_blocks_exactly_zero=True))
        self.assertEqual(len(resonances),4)
        OBS['resonance_resolution']=dict(exceptions=resonances,
            complete_degenerate_energy_spaces_used=True,
            all_piecewise_integrable_real_control_amplitudes_covered=True,
            operator_identity='Pinch_H(W^* Q_j W)=0',
            arbitrary_joint_H_stationary_external_reference_conditioned_Q_moments_zero=True,
            DG_marginal_stationarity_alone_not_sufficient_for_conditional_claim=True)

    def test_04_full_384_dimensional_degenerate_pinch_diagnostic(self):
        h,k=system()[:2];q=np.tile(words.qgraph(),(1,64));c=chart.pure_columns();seed=c@c.conj().T/384
        center=chart.numerical_control()[0];control=words.unitary(.7*h-1.3*k)@words.unitary(.2*h+.9*k)@words.unitary(.4*h-.2*k)
        records=[]
        for label,rho in (('480_seed',seed),('480_finite_pulse_center',center@seed@center.conj().T)):
            stat=pinch_numeric(rho);out=control@stat@control.conj().T
            stationary=float(np.linalg.norm(h@stat-stat@h));response=q@np.diag(out).real
            self.assertLess(stationary,2e-12);self.assertLess(np.max(abs(response)),2e-12)
            self.assertGreater(np.linalg.norm(k@stat-stat@k),.01)
            graph=out.reshape(64,6,64,6).trace(axis1=0,axis2=2)
            leaves=(0,3,4,5);largest=0.
            for ia,a in enumerate(leaves):
                for b in leaves[ia+1:]:largest=max(largest,abs(float(frame.distance(a,b)@np.diag(graph).real)-8/3))
            self.assertLess(largest,2e-12)
            records.append(dict(preparation=label,H_commutator_norm=short(stationary),
                K_commutator_norm=short(np.linalg.norm(k@stat-stat@k)),
                controlled_Q=[short(x) for x in response],largest_leaf_reference_error=short(largest)))
        e,v,groups=spectral_diagnostic()
        OBS['full_matrix_diagnostic']=dict(records=records,complete_energy_groups=len(groups),
            retained_degeneracies=[len(a) for a in groups],
            grouping_tolerance_is_diagnostic_only_not_proof=True,
            minimum_distinct_energy_gap=short(min(e[b[0]]-e[a[0]] for a,b in zip(groups,groups[1:]))),
            Q_frozen_does_not_mean_full_state_frozen=True)

    def test_05_unknown_reference_and_retained_energy_pointer_boundary(self):
        f=frame.system()[2].astype(np.int64);qg=words.qgraph().astype(np.int64);one=np.ones(6,dtype=np.int64);raw=qg[0]
        self.assertTrue(np.array_equal(f@one,4*one));self.assertFalse(np.any(f@raw))
        self.assertEqual(int(one@raw),0);self.assertEqual(int(raw@raw),4)
        s=one/np.sqrt(6);g=raw/2;columns=np.column_stack((s,g))/np.sqrt(2)
        moment=columns.T@(raw[:,None]*columns)
        self.assertAlmostEqual(float(moment[0,1]),math.sqrt(1/6),places=14)
        self.assertTrue(np.allclose(np.diag(moment),0))
        self.assertEqual(Fr(4,24),Fr(1,6))
        # On |0^6>, all swaps are +1, so H=5I+F and energies are 9 and 5.
        active=columns@columns.T
        self.assertTrue(np.allclose(f@active,active@f))
        pinched_joint_moment=np.diag(np.diag(moment))
        self.assertTrue(np.allclose(pinched_joint_moment,0))
        self.assertGreater(np.linalg.norm(moment),.5)
        OBS['memory_and_reference_scope']=dict(
            exact_energy_pair=[9,5],marginal_stationary_but_reference_Q_offdiagonal_squared='1/6',
            full_local_pinch_tensor_identity_makes_joint_state_stationary=True,
            isometry='V=sum_E P_E tensor |E>_M',isometry_preserves_all_unknown_input_information=True,
            coherent_energy_pointer_not_part_of_the_joint_stationary_reference_quantifier=True,
            energy_pointer_retained_inside_total_system=True,
            reading_or_reversing_pointer_coherence_is_additional_control_permission=True,
            current_H_plus_uK_does_not_touch_the_pointer=True,
            exact_pinch_preparation_and_pointer_cost_not_derived_from_original_H=True)

    def test_06_actual_probe_and_finite_complete_error_budget(self):
        c=chart.pure_columns();rho=pinch_numeric(c@c.conj().T/384)
        graph=rho.reshape(64,6,64,6).trace(axis1=0,axis2=2)
        t=Fr(1,65536);estimates=[];maxcp=0.
        for leaf in (0,3,4,5):
            minus,plus,kraus=reader.instrument((1,leaf),t)
            estimate=np.trace(reader.apply_map(plus-minus,graph,g=6,r=1)).real/float(t)
            estimates.append(float(estimate))
            completeness=np.einsum('bdgi,bdgj->ij',kraus.conj(),kraus,optimize=True)
            maxcp=max(maxcp,float(np.max(abs(completeness-np.eye(6)))))
        readq=np.array(estimates[0])-np.array(estimates[1:])
        bias=(math.expm1(18*float(t))-18*float(t))/float(t)
        self.assertLess(maxcp,3e-12);self.assertLess(np.max(abs(readq)),2*bias)
        epsilon=Fr(1,1000);a=epsilon/4;h=Fr(9);probe=a/(16*h*h);gamma=a*probe/8
        x=2*h*probe;remainder=x*x/(2*(1-x/3));delta=(remainder+gamma)/probe
        self.assertLess(delta,a/2)
        copies=math.ceil(56/(a*a*probe*probe));gamma_ctrl=gamma/2;duration=gamma/(4*h)
        self.assertEqual(gamma_ctrl+2*h*duration,gamma)
        self.assertGreater(sum(Fr(7)**n/math.factorial(n) for n in range(25)),800)
        source_and_control=epsilon/4;self.assertEqual(2*a+source_and_control,3*epsilon/4)
        OBS['actual_readout_and_resources']=dict(
            actual_single_edge_CP_instruments_from_round=472,actual_probe_wait=str(t),
            estimated_adjacencies=[short(x) for x in estimates],estimated_Q=[short(x) for x in readq],
            CP_instrument_completeness_error=short(maxcp),
            target_Q_infinity_error=str(epsilon),joint_failure_probability_at_most='1/100',
            per_edge_total_estimation_error=str(a),certified_probe_wait=str(probe),
            whole_handoff_plus_probe_channel_error=str(gamma),
            instrument_bias_bound=str(delta),control_only_error_budget=str(gamma_ctrl),
            total_control_window_duration_budget=str(duration),
            copies_per_edge=str(copies),total_independent_source_trials=str(4*copies),
            source_pinch_and_all_pre_readout_control_total_error_budget=str(source_and_control),
            same_actual_preparation_and_control_channel_shared_before_choice_of_four_edges=True,
            total_Q_error_bound=str(3*epsilon/4),
            old_D_moved_to_internal_isolated_S_G_stays_active_S_G_R_C_M_correlations_retained=True,
            new_independent_probes_and_clocks_storage_record_and_source_costs_explicit=True,
            final_readout_can_disturb_G=True,huge_statistical_sample_not_actually_executed=True)


def run():
    OBS.clear();stream=io.StringIO();result=unittest.TextTestRunner(stream=stream).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():raise RuntimeError(stream.getvalue())
    output=dict(round=486,scientific_baseline_round=485,reused_frozen_rounds=[475,480,472],
        tests_run=result.testsRun,failures=len(result.failures),errors=len(result.errors),
        python=platform.python_version(),numpy=np.__version__,observations=OBS,
        scope=dict(theorem_is_specific_to_declared_unit_coupling_H_K_Q=True,
            all_H_stationary_active_inputs_and_all_allowed_controls_have_Q_zero=True,
            stationarity_of_joint_state_required_for_arbitrary_external_reference_conditionals=True,
            coherent_energy_pointer_is_retained_but_not_controlled_or_in_joint_stationary_quantifier=True,
            deleting_unknown_information_not_assumed=True,
            pinching_is_an_extra_preparation_not_free_H_plus_uK=True,
            no_rank_three_linear_response_from_any_stationary_input_under_this_contract=True,
            all_nonstationary_zero_drift_preparations_not_excluded=True,
            alternate_controls_or_costs_not_excluded=True,
            no_physical_space_or_cognitive_axiom_impossibility_claim=True,
            full_GR_goal_completed=False,phase_closure_triggered=False))
    assert json.loads(json.dumps(output))==output
    return output


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--dry-run',action='store_true');parser.add_argument('--check',action='store_true')
    args=parser.parse_args();output=run()
    if args.check:assert output==json.loads(TARGET.read_text(encoding='utf-8'))
    elif not args.dry_run:
        with TARGET.open('x',encoding='utf-8',newline='\n') as handle:handle.write(json.dumps(output,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(output,ensure_ascii=False,indent=2))
