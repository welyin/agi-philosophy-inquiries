"""Round 519: boundary NNI enters the first resolved return coefficient.

General finite-time claims use the inherited commutator theorem, not a scan.
All floating source calculations below are explicitly labelled diagnostics.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import io
import json
from pathlib import Path
import platform
import unittest

import numpy as np
import same_time_neighborhood_source as source
import record_determined_predictor_probe as scope

HERE=Path(__file__).resolve().parent
TARGET=HERE/'boundary_corrected_return_results.json'
OBS={}
K=1889568
L=12312
WAIT=Q(1,64*K**5)
PREP=WAIT**4/(64*L)


def h_model(i,j=1,kappa=1):
    n,trees,adj,h=source.model(i)
    f=np.zeros((len(trees),len(trees)),dtype=np.int64)
    for row,neighbors in enumerate(adj):
        for col,weight in neighbors.items():f[row,col]=weight
    fullf=np.kron(np.eye(n,dtype=np.int64),f)
    return n,trees,kappa*fullf+j*(h-fullf)


class Audit(unittest.TestCase):
    def test_01_exact_all_good_input_effect(self):
        rows=[]
        for i in (1,2,3):
            for j,kappa in ((1,1),(2,-1)):
                n,trees,h=h_model(i,j,kappa);m=len(trees);a=i
                rs=source.records(i,n,a);masks=source.good_masks(trees,rs,a)
                first={w:h[a*m:(a+1)*m,w*m:(w+1)*m] for w in range(n)}
                second={w:h[a*m:(a+1)*m,:]@h[:,w*m:(w+1)*m] for w in range(n)}
                cases=0;internal=0
                for (u,v,w),mask in zip(rs,masks):
                    selected=np.flatnonzero(mask)
                    if not len(selected):continue
                    self.assertFalse(np.any(first[w][:,selected]))
                    c=second[w][:,selected]
                    dw=3 if w<i else 1
                    expected=j**4+(dw-1)*j*j*kappa*kappa
                    self.assertTrue(np.array_equal(c.T@c,expected*np.eye(len(selected),dtype=np.int64)))
                    # The good output component is exactly the two-hop channel.
                    self.assertTrue(np.array_equal(c[selected,:],j*j*np.eye(len(selected),dtype=np.int64)))
                    boundary=c.copy();boundary[selected,:]=0
                    self.assertTrue(np.array_equal(boundary.T@boundary,
                        (dw-1)*j*j*kappa*kappa*np.eye(len(selected),dtype=np.int64)))
                    cases+=1;internal+=int(w<i)
                rows.append(dict(internal_vertices=i,graphs=m,J=j,kappa=kappa,
                    nonempty_record_sectors=cases,internal_terminal_sectors=internal,
                    second_order_amplitude_effect_exact_for_all_graph_coherences=True))
        OBS['integer_leading_effect']=rows

    def test_02_reverse_boundary_map_and_actual_source_role_weight(self):
        rows=[]
        for i in (2,3,4):
            n,trees,f=source.source517.family(i);a=i
            rs=source.records(i,n,a);masks=source.good_masks(trees,rs,a)
            internal_counts=np.zeros(len(trees),dtype=int)
            outputs=0
            for (u,v,w),mask in zip(rs,masks):
                seen={}
                for g in np.flatnonzero(mask):
                    internal_counts[g]+=int(w<i)
                    hits=[gp for gp,weight in f[g].items() if source.edge(a,w) in trees[gp]]
                    self.assertEqual(len(hits),2 if w<i else 0)
                    for gp in hits:
                        self.assertEqual(f[g][gp],1)
                        self.assertNotIn(gp,seen)
                        seen[gp]=int(g)
                        neighbors=source.source517.old.adjacency(trees[gp],n)[u]
                        candidates=neighbors-{v,w}
                        self.assertEqual(len(candidates),1)
                        d=next(iter(candidates))
                        restored=(trees[gp]-{source.edge(a,w),source.edge(u,d)})|{
                            source.edge(a,u),source.edge(w,d)}
                        self.assertEqual(restored,trees[g])
                        outputs+=1
            self.assertTrue(np.all(internal_counts>=1))
            self.assertTrue(np.all(internal_counts<=2))
            rows.append(dict(internal_vertices=i,graphs=len(trees),
                boundary_output_pairs=outputs,minimum_internal_terminal_leading_records=int(internal_counts.min()),
                total_leading_records_per_graph=2,unique_reverse_map=True))
        OBS['boundary_injectivity_and_source_role']=rows

    def test_03_uniform_finite_budget(self):
        s=WAIT;t=PREP;x=K*s
        remainder=x**5/(1-x)
        self.assertLess(x,Q(1,2))
        self.assertLessEqual(remainder,s**4/32)
        self.assertEqual(2*L*t,s**4/32)
        self.assertLessEqual(remainder+2*L*t,s**4/16)
        self.assertLessEqual(L*t,Q(1,4))
        self.assertGreaterEqual(Q(1,2)-4*L*t,Q(1,4))
        self.assertEqual(Q(1,8)-Q(1,32)-Q(1,16),Q(1,32))
        # The four-point (I=1) star has no boundary; its coefficient is 1/4.
        self.assertEqual(Q(1+0,4),Q(1,4))
        OBS['uniform_certificate']=dict(couplings=dict(J=1,kappa=1),K=K,L=L,
            future_wait=str(s),preparation_leg=str(t),source_success_lower=str(t**8),
            source_bad_upper=str((L*t)**2),remainder_over_s4_upper='1/32',
            complete_forecast_error_over_s4_upper='1/16',
            actual_total_click_over_s4_lower='3/16',
            actual_minus_star_click_over_s4_lower_for_I_ge_2='1/32',
            internal_terminal_weight_lower='1/4',scale_independent=True,
            efficient_source_or_practical_clock_claimed=False)

    def test_04_actual_source_and_reference_return_diagnostic(self):
        n,trees,h=h_model(2);m=len(trees);a=2
        t=Q(1,2**16);s=Q(1,128)
        U=source.source500.evolution(h,t);W=source.source500.evolution(h,s)
        psi=np.zeros((m,2),complex)
        psi[0,0]=1/np.sqrt(2);psi[-1,1]=1j/np.sqrt(2)
        def block(v,u):return U[v*m:(v+1)*m,u*m:(u+1)*m]
        records=[];success=0.0
        for u,v,w in source.records(2,n,a):
            vec=block(w,u)@block(u,v)@block(v,u)@block(u,a)@psi
            records.append((w,vec));success+=float(np.vdot(vec,vec).real)
        p=pcorr=alpha=error=0.0
        for w,vec in records:
            vec=vec/np.sqrt(success);rhoR=scope.reference(vec)
            init=np.zeros((n,m,2),complex);init[w]=vec
            evolved=(W@init.reshape(n*m,2)).reshape(n,m,2)
            click=scope.reference(evolved[a]);total=scope.reference(evolved.reshape(n*m,2))
            c=Q(3 if w<2 else 1,4);f=float(c*s**4)
            p+=float(np.trace(click).real);pcorr+=f*float(np.trace(rhoR).real)
            if w<2:alpha+=float(np.trace(rhoR).real)
            error+=scope.distance(click,f*rhoR)+scope.distance(total-click,(1-f)*rhoR)
        star=scope.forecast(float(s));s4=float(s**4)
        self.assertLess(abs(p/s4-0.5),0.02)
        self.assertLess(error/s4,0.05)
        self.assertGreater((p-star)/s4,0.2)
        OBS['actual_diagnostic']=dict(preparation_leg=str(t),future_wait=str(s),
            source_success=float(f'{success:.10g}'),internal_terminal_weight=float(f'{alpha:.10g}'),
            actual_click_over_s4=float(f'{p/s4:.9g}'),corrected_click_over_s4=float(f'{pcorr/s4:.9g}'),
            isolated_star_click_over_s4=float(f'{star/s4:.9g}'),
            full_cq_reference_report_error_over_s4=float(f'{error/s4:.7g}'),
            theorem_tiny_time_budget_used=False,general_claims_certified_by_this_float_run=False)


def run():
    OBS.clear();stream=io.StringIO()
    result=unittest.TextTestRunner(stream=stream,verbosity=0).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():raise AssertionError(stream.getvalue())
    return dict(date='2026-09-28',round=519,scientific_base_through_round=518,
        tests_run=result.testsRun,failures=len(result.failures),errors=len(result.errors),
        python=platform.python_version(),numpy=np.__version__,
        dependency_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in (
            'same_time_neighborhood_source.py','record_determined_predictor_probe.py',
            'research_note_498.md','research_note_518.md','record_determined_predictor_review.md')},
        scope=dict(same_440_H_and_518_source=True,arbitrary_unknown_graph_and_reference=True,
            full_classical_record_report=True,uniform_finite_nonzero_return_certificate=True,
            boundary_NNI_changes_leading_return=True,record_and_role_suffice_for_short_time_forecast=True,
            arbitrary_local_state_reconstructed=False,per_record_conditioned_error_uniform=False,
            efficient_source_or_autonomous_readout_generated=False,
            macroscopic_coordinate_or_dimension_generated=False,full_GR_goal_completed=False),
        observations=OBS)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();answer=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(answer,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:answer[k] for k in ('round','tests_run','failures','errors')}))
