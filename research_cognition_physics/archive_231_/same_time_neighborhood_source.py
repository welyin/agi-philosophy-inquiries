"""Round 518: simultaneous radius-two source from four actual port records.

General bounds are record-space operator proofs in the note. Exact rational
certificates for the unchanged small model are checks, not extrapolations.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import io
import itertools
import json
import math
from pathlib import Path
import platform
import unittest

import numpy as np
import all_scale_monitored_graph_source as source517
import current_neighbor_detection as source500

HERE=Path(__file__).resolve().parent
TARGET=HERE/'same_time_neighborhood_source_results.json'
OBS={}
TIME=Q(1,2**20)
L=10368+1944


def edge(a,b):
    return tuple(sorted((a,b)))


def records(i,n,root):
    return [(u,v,w) for u in range(i) if u!=root
            for v,w in itertools.permutations([x for x in range(n) if x not in (root,u)],2)]


def good_masks(trees,rs,root):
    return np.array([[{edge(root,u),edge(u,v),edge(u,w)}<=g for g in trees]
                     for u,v,w in rs],dtype=bool)


def ball(tree,n,root,radius):
    adj=source517.old.adjacency(tree,n)
    seen={root}
    for _ in range(radius):
        seen |= {w for v in seen for w in adj[v]}
    return seen


def model(i):
    n,trees,f=source517.family(i)
    m=len(trees);dim=n*m
    h=np.zeros((dim,dim),dtype=np.int64)
    rows,cols,values=source517.sparse_h(trees,f,n)
    np.add.at(h,(rows,cols),values)
    h -= (n-1)*np.eye(dim,dtype=np.int64)
    return n,trees,f,h


def product(a,b):
    ar,ai=a;br,bi=b
    return ar@br-ai@bi, ar@bi+ai@br


def gram(k,mask=None):
    r,z=k
    if mask is not None:
        r=r[mask];z=z[mask]
    return r.T@r+z.T@z, r.T@z-z.T@r


def gershgorin(r,z,scale):
    low=[];high=[]
    for j in range(len(r)):
        radius=sum(abs(int(r[j,k]))+abs(int(z[j,k])) for k in range(len(r)) if k!=j)
        low.append(Q(int(r[j,j])-radius,scale))
        high.append(Q(int(r[j,j])+radius,scale))
    return min(low),max(high)


def exact_four_step():
    n,trees,f,h=model(2);m=len(trees);a=2
    rs=records(2,n,a);pi=good_masks(trees,rs,a)
    r,z,scale,tail,norm=source500.exact_taylor(h,denominator=TIME.denominator,degree=8)
    def block(v,u):
        ix=np.ix_(range(v*m,(v+1)*m),range(u*m,(u+1)*m))
        return r[ix],z[ix]
    er=np.zeros((m,m),dtype=object);ei=er.copy()
    br=er.copy();bi=er.copy();selected=[]
    for row,(u,v,w) in enumerate(rs):
        k=product(block(w,u),product(block(u,v),product(block(v,u),block(u,a))))
        gr,gi=gram(k);er+=gr;ei+=gi
        qr,qi=gram(k,~pi[row]);br+=qr;bi+=qi
        if not np.any(pi[row]):
            selected.append((row,k))
    # Complete non-stay record maps have norm <=3t; no number-of-records factor.
    a_bound=3*TIME
    amplitude_error=4*tail*(a_bound+tail)**3
    effect_error=2*(a_bound+tail)**4*amplitude_error+amplitude_error**2
    denominator=scale**8
    elo,ehi=gershgorin(er,ei,denominator)
    blo,bhi=gershgorin(br,bi,denominator)
    witness=None
    for row,k in selected:
        gr,gi=gram(k)
        for column in range(m):
            # The same norm bound also bounds a single recorded submap.
            lower=Q(int(gr[column,column]),denominator)-effect_error
            if lower>0:
                witness=(rs[row],column,lower)
                break
        if witness is not None:
            break
    return dict(n=n,trees=trees,root=a,records=rs,pi=pi,
        success_lower=elo-effect_error,success_upper=ehi+effect_error,
        bad_upper=bhi+effect_error,amplitude_tail=amplitude_error,
        effect_tail=effect_error,unitary_tail=tail,h_norm_bound=norm,witness=witness)


class Audit(unittest.TestCase):
    def test_01_leading_complete_record_effect_and_radius_two(self):
        rows=[]
        for i in (1,2,3,4):
            n,trees,f=source517.family(i);a=i
            rs=records(i,n,a);pi=good_masks(trees,rs,a)
            self.assertTrue(np.all(pi.sum(axis=0)==2))
            for row,(u,v,w) in enumerate(rs):
                for g in np.flatnonzero(pi[row]):
                    self.assertEqual(ball(trees[g],n,a,2),{a,u,v,w})
            rows.append(dict(internal=i,vertices=n,graphs=len(trees),
                accepted_record_labels=len(rs),nonzero_leading_records=int(np.count_nonzero(pi.any(axis=1))),
                leading_effect=2,same_time_ball_equals_recorded_set=True))
        OBS['exact_leading_source']=rows

    def test_02_original_H_and_size_free_local_commutators(self):
        n,trees,f,h=model(2)
        old_trees,old_h,_,_=source500.six_model()
        permutation=[1,2,0,3,4,5]
        lookup={g:j for j,g in enumerate(old_trees)}
        mapped=[lookup[frozenset(edge(permutation[a],permutation[b]) for a,b in g)] for g in trees]
        order=[permutation[v]*len(trees)+mapped[g] for v in range(n) for g in range(len(trees))]
        self.assertTrue(np.array_equal(h,old_h[np.ix_(order,order)]))
        rows=[]
        for i in (2,3,4):
            n,trees,f=source517.family(i);m=len(trees)
            adj=[source517.old.adjacency(g,n) for g in trees]
            a=i;rs=records(i,n,a);pi=good_masks(trees,rs,a)
            adjacency_row=0;good_row=0
            for g in range(m):
                for v in range(n):
                    adjacency_row=max(adjacency_row,sum(w*len(adj[g][v]^adj[gp][v]) for gp,w in f[g].items()))
                changes=np.zeros(len(rs),dtype=np.int64)
                for gp,w in f[g].items():
                    changes += w*(pi[:,gp]!=pi[:,g])
                good_row=max(good_row,int(changes.max()))
            self.assertLessEqual(adjacency_row,48)
            self.assertLessEqual(good_row,24)
            rows.append(dict(internal=i,commutator_F_A_absolute_row_max=adjacency_row,
                             recorded_full_star_commutator_row_max=good_row))
        OBS['original_model_and_local_bounds']=dict(original_six_H_equal=True,counts=rows)

    def test_03_uniform_operator_constants_and_finite_budget(self):
        cs=[48*(Q(9,2)-k)+18 for k in range(1,5)]
        self.assertEqual(27*sum(cs),L)
        self.assertLessEqual(L*TIME,Q(1,4))
        success=TIME**8;bad=L**2*TIME**10;eps=L**2*TIME**2
        self.assertEqual(success,Q(1,2**160))
        self.assertLess(eps,Q(1,5000))
        delta=success/10000
        robust=(bad+delta)/(success-delta)
        self.assertLess(robust,Q(1,4000))
        OBS['uniform_certificate']=dict(couplings=dict(J=1,kappa=1),time_per_leg=str(TIME),
            total_time=str(4*TIME),leg_error_coefficients=list(map(str,cs)),
            all_record_amplitude_error_coefficient=L,success_lower=str(success),
            conditional_bad_upper=str(eps),bad_below_one_per_5000=True,
            full_process_error_example=str(delta),conditional_bad_with_error_upper=str(robust),
            error_example_is_contract_not_constructed_hardware=True,
            complete_occupied_port_readings=4,uncompressed_occupancy_record_bits='4*N')

    def test_04_rational_actual_complete_source_certificate(self):
        got=exact_four_step()
        self.assertGreater(got['success_lower'],TIME**8)
        self.assertLess(got['success_upper'],3*TIME**8)
        self.assertLess(got['bad_upper'],L**2*TIME**10)
        self.assertIsNotNone(got['witness'])
        rec,col,lower=got['witness']
        OBS['finite_integer_certificate']=dict(active_dimension=36,source_graph_dimension=6,
            taylor_order=8,unitary_norm_bound=got['h_norm_bound'],
            certified_success_scaled_lower=str(got['success_lower']/TIME**8),
            certified_success_scaled_upper=str(got['success_upper']/TIME**8),
            certified_bad_scaled_upper=str(got['bad_upper']/TIME**10),
            unitary_remainder=str(got['unitary_tail']),
            complete_record_amplitude_remainder=str(got['amplitude_tail']),
            forbidden_star_record=list(rec),witness_initial_graph=col,
            forbidden_record_probability_lower=str(lower),
            forbidden_record_conditional_error_exactly_one=True,
            all_input_bounds_from_integer_Gershgorin=True)

    def test_05_full_actual_instrument_with_unknown_reference(self):
        n,trees,f,h=model(2);m=len(trees);a=2
        time=Q(1,8)
        u=source500.evolution(h,time)
        psi=np.zeros((m,2),complex)
        psi[0,0]=1/np.sqrt(2);psi[-1,1]=1j/np.sqrt(2)
        accepted=set(records(2,n,a))
        initial_r=psi.T@psi.conj()
        final_r=np.zeros((2,2),complex)
        success=bad=0.0;count=0
        def visit(last,history,state):
            nonlocal final_r,success,bad,count
            if len(history)==4:
                count+=1;final_r+=state.T@state.conj()
                b,c,d,e=history
                if b==d and (b,c,e) in accepted:
                    p=float(np.vdot(state,state).real);success+=p
                    mask=good_masks(trees,[(b,c,e)],a)[0]
                    bad+=float(np.vdot(state[~mask],state[~mask]).real)
                return
            for nxt in range(n):
                block=u[nxt*m:(nxt+1)*m,last*m:(last+1)*m]
                visit(nxt,history+[nxt],block@state)
        visit(a,[],psi)
        self.assertEqual(count,n**4)
        self.assertLess(float(np.max(np.abs(final_r-initial_r))),3e-12)
        self.assertGreater(success,0)
        self.assertGreater(bad,0)
        self.assertLess(bad,success)
        OBS['complete_history_diagnostic']=dict(histories=count,reference_dimension=2,
            time_per_leg=str(time),inside_uniform_small_time_certificate=False,
            success_probability=float(f'{success:.10g}'),
            conditional_same_time_bad_probability=float(f'{bad/success:.10g}'),
            full_unselected_reference_unchanged=True,failures_retained=True)

    def test_06_leaf_source_and_delivery_scope(self):
        n,trees,f=source517.family(2);a=0
        rs=records(2,n,a);pi=good_masks(trees,rs,a)
        counterexamples=0
        for row,(u,v,w) in enumerate(rs):
            for g in np.flatnonzero(pi[row]):
                self.assertFalse(ball(trees[g],n,a,2)<={a,u,v,w})
                counterexamples+=1
        self.assertGreater(counterexamples,0)
        delay=Q(1,2**22)
        delayed=(L*TIME+24*delay)**2
        self.assertLess(delayed,Q(1,5000))
        OBS['scope_boundaries']=dict(internal_root_counterexamples=counterexamples,
            leaf_root_required=True,delivery_delay_example=str(delay),
            delayed_bad_upper=str(delayed),
            delivery_facilities_assumed_not_generated=True,
            failure_does_not_reset_to_root=True,large_radius_buffer_or_predictor_claimed=False)


def run():
    OBS.clear();stream=io.StringIO()
    result=unittest.TextTestRunner(stream=stream,verbosity=0).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():raise AssertionError(stream.getvalue())
    return dict(date='2026-09-28',round=518,scientific_base_through_round=517,
        tests_run=result.testsRun,failures=len(result.failures),errors=len(result.errors),
        python=platform.python_version(),numpy=np.__version__,
        dependency_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in (
            'all_scale_monitored_graph_source.py','current_neighbor_detection.py',
            'research_note_496.md','research_note_499.md','research_note_500.md','research_note_502.md')},
        scope=dict(unchanged_440_H_and_declared_port_instrument=True,
            same_final_time_joint_radius_two_source=True,
            arbitrary_unknown_graph_and_passive_reference=True,
            size_independent_success_and_error_bounds=True,
            records_and_failure_instrument_retained=True,
            per_record_posterior_error_guarantee=False,
            automatic_failure_reset=False,arbitrary_large_buffer_generated=False,
            good_domain_predictor_generated=False,
            autonomous_readout_clock_or_delivery_generated=False,
            physical_three_dimensional_space_generated=False,full_GR_goal_completed=False),
        observations=OBS)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('round','tests_run','failures','errors')}))
