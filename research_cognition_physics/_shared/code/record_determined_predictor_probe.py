"""Unnumbered scope correction after 518; reuse 498/499, no new round.

Integer identities check the prepared local factor and star dynamics. Floating
point source checks are diagnostics; the note contains the general proof.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import io
import itertools
import json
from pathlib import Path
import platform
import unittest

import numpy as np
import same_time_neighborhood_source as source

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'record_determined_predictor_probe_results.json'
OBS = {}
K = 1889568


def distance(x, y):
    return float(np.abs(np.linalg.eigvalsh((x-y+(x-y).conj().T)/2)).sum()/2)


def density(v):
    v = v.reshape(-1)
    return np.outer(v, v.conj())


def reference(v):
    return v.T @ v.conj()


def star():
    # Order a,u,v,w, with u the centre. Full SWAP atoms, then one excitation.
    h = np.zeros((4,4), dtype=np.int64)
    for v in (0,2,3):
        swap = np.eye(4, dtype=np.int64)
        swap[1,1] = swap[v,v] = 0
        swap[1,v] = swap[v,1] = 1
        h += swap
    return h, 3*np.eye(4, dtype=np.int64)-h


def forecast(s):
    return float(abs(Q(1,4)-np.exp(1j*s)/3+np.exp(4j*s)/12)**2)


class Audit(unittest.TestCase):
    def test_01_exact_prepared_factor(self):
        rows = []
        for i in (1,2,3,4):
            n, trees, _ = source.source517.family(i)
            a = i
            rs = source.records(i,n,a)
            masks = source.good_masks(trees,rs,a)
            cases = 0
            for row,(u,v,w) in enumerate(rs):
                selected = np.flatnonzero(masks[row])
                A = {a,u,v,w}
                edges = {source.edge(a,u),source.edge(u,v),source.edge(u,w)}
                outside = []
                for g in selected:
                    inner = frozenset(e for e in trees[g] if set(e) <= A)
                    self.assertEqual(inner,edges)
                    self.assertEqual(source.ball(trees[g],n,a,2),A)
                    # Distinct good graphs differ solely in complement factors.
                    outside.append(frozenset(trees[g]-inner))
                    cases += 1
                self.assertEqual(len(set(outside)),len(outside))
            rows.append(dict(internal=i,graphs=len(trees),good_graph_record_pairs=cases,
                local_edge_bits_fixed=True,local_data_bits_fixed_by_single_excitation=True))
        OBS['exact_local_factor'] = rows

    def test_02_full_atom_star_and_spectral_formula(self):
        h,l = star()
        c = np.array([1,-3,1,1],dtype=np.int64)
        # All projectors scaled by 12, so the spectral check is exact integer.
        p0 = 3*np.ones((4,4),dtype=np.int64)
        p4 = np.outer(c,c)
        p1 = 12*np.eye(4,dtype=np.int64)-p0-p4
        projectors = [(0,p0),(1,p1),(4,p4)]
        for ev,p in projectors:
            self.assertTrue(np.array_equal(p@p,12*p))
            self.assertTrue(np.array_equal(l@p,ev*p))
        for (_,p),(_,q) in itertools.combinations(projectors,2):
            self.assertFalse(np.any(p@q))
        coefficients = [Q(int(p[0,3]),12) for _,p in projectors]
        self.assertEqual(coefficients,[Q(1,4),-Q(1,3),Q(1,12)])
        self.assertEqual(sum(coefficients),0)
        self.assertEqual(sum(ev*c for (ev,_),c in zip(projectors,coefficients)),0)
        self.assertEqual(sum(ev**2*c for (ev,_),c in zip(projectors,coefficients)),1)
        err = 0.0
        for s in (0.0,0.01,0.125,0.5,1.0):
            U = source.source500.evolution(h,s)
            err = max(err,abs(float(abs(U[0,3])**2)-forecast(s)))
        self.assertLess(err,2e-14)
        OBS['star_spectrum'] = dict(one_excitation_generator=h.tolist(),
            laplacian_eigenvalues=[0,1,1,4],leaf_to_other_leaf_projector_entries=list(map(str,coefficients)),
            leading_probability_coefficient='1/4',formula_numerically_cross_checked=True)

    def test_03_actual_source_reference_and_future_report(self):
        n,trees,_,h = source.model(2)
        m=len(trees);a=2;t=Q(1,64);s=Q(1,32)
        U=source.source500.evolution(h,t)
        W=source.source500.evolution(h,s)
        psi=np.zeros((m,2),complex)
        psi[0,0]=1/np.sqrt(2);psi[-1,1]=1j/np.sqrt(2)
        rs=source.records(2,n,a);masks=source.good_masks(trees,rs,a)
        states=[];projected=[]
        def block(v,u):
            return U[v*m:(v+1)*m,u*m:(u+1)*m]
        for (u,v,w),mask in zip(rs,masks):
            state=block(w,u)@block(u,v)@block(v,u)@block(u,a)@psi
            good=state*mask[:,None]
            states.append(state);projected.append(good)
        success=sum(float(np.vdot(v,v).real) for v in states)
        goodprob=sum(float(np.vdot(v,v).real) for v in projected)
        q=1-goodprob/success
        self.assertGreater(success,0);self.assertGreater(q,0);self.assertLess(q,1)
        dsource=dref=dfuture=pclick=0.0
        f=forecast(float(s))
        for (_,_,w),v,g in zip(rs,states,projected):
            v=v/np.sqrt(success);g=g/np.sqrt(goodprob)
            dsource+=distance(density(v),density(g))
            rhoR=reference(v);goodR=reference(g)
            dref+=distance(rhoR,goodR)
            init=np.zeros((n,m,2),complex);init[w]=v
            evolved=(W@init.reshape(n*m,2)).reshape(n,m,2)
            clicked=reference(evolved[a])
            pclick+=float(np.trace(clicked).real)
            total=reference(evolved.reshape(n*m,2))
            dfuture+=distance(clicked,f*rhoR)+distance(total-clicked,(1-f)*rhoR)
        self.assertLessEqual(dsource,np.sqrt(q)+1e-11)
        self.assertLessEqual(dref,dsource+1e-11)
        self.assertGreater(dref,1e-8) # Projection need not preserve the cq reference marginal.
        self.assertLessEqual(dfuture,1+1e-11)
        self.assertLessEqual(pclick,float(s*s)+1e-11)
        OBS['actual_source_diagnostic'] = dict(time_per_leg=str(t),forecast_wait=str(s),
            active_dimension=n*m,accepted_classical_labels=len(rs),reference_dimension=2,
            success_probability=float(f'{success:.10g}'),conditional_bad=float(f'{q:.10g}'),
            joint_source_to_normalized_good_distance=float(f'{dsource:.10g}'),
            cq_reference_marginal_change=float(f'{dref:.10g}'),
            future_classical_root_reference_report_distance=float(f'{dfuture:.10g}'),
            actual_future_root_click_probability=float(f'{pclick:.10g}'),
            within_518_uniform_small_time_window=False,numerics_prove_general_bound=False,
            saved_518_full_instrument_and_failure_analysis_reused=True)

    def test_04_information_scope_and_conservative_certificate(self):
        # A fixed graph / record does not determine unrestricted data input.
        zero=np.zeros(16,complex);zero[0]=1
        one=np.zeros(16,complex);one[1]=1
        self.assertAlmostEqual(distance(density(zero),density(one)),1)
        # A classical record marginal does not reproduce a coherent record+R.
        bell=np.array([1,0,0,1],complex)/np.sqrt(2)
        dephased=np.diag([0.5,0,0,0.5])
        self.assertAlmostEqual(distance(density(bell),dephased),0.5)
        # The leaf coupling to its complement has norm one on every graph.
        n,trees,_,h=source.model(2);m=len(trees);a=2
        inside=np.arange(a*m,(a+1)*m)
        outside=np.array([v for v in range(n*m) if v not in set(inside)])
        b=h[np.ix_(inside,outside)]
        self.assertTrue(np.array_equal(b@b.T,np.eye(m,dtype=np.int64)))
        # Concavity on z in [0,1] and 2*tanh(pi/2)>1 give E2>=min(1,K|s|).
        # With K>=1, this exceeds min(1,s*s), the actual vacuum error bound.
        self.assertGreaterEqual(K,1)
        eps=Q(12312,2**20)
        self.assertLess(eps,Q(1,80))
        OBS['scope_and_resolution'] = dict(arbitrary_data_same_record_can_have_distance=1,
            coherent_record_loss_example_distance='1/2',
            general_498_constant=K,zero_predictor_actual_error_bound='min(1,s^2)',
            exact_leaf_cross_block_BBt_identity=True,dominance_condition_exactly_checked='K >= 1',
            zero_click_predictor_is_not_excluded_by_inherited_tolerance=True,
            source_example_sqrt_bad_upper=str(eps),
            physical_signal_absent_claimed=False,macroscopic_space_or_dimension_selected=False)


def run():
    OBS.clear();stream=io.StringIO()
    result=unittest.TextTestRunner(stream=stream,verbosity=0).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise AssertionError(stream.getvalue())
    return dict(date='2026-09-28',scientific_baseline_round=518,
        numbered_round_created=False,numbered_scientific_test_increment=0,
        diagnostic_tests=result.testsRun,failures=len(result.failures),errors=len(result.errors),
        python=platform.python_version(),numpy=np.__version__,
        dependency_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in (
            'same_time_neighborhood_source.py','same_time_neighborhood_source_results.json',
            'research_note_498.md','research_note_499.md','research_note_518.md')},
        scope=dict(prepared_good_local_factor_is_record_determined=True,
            quantum_snapshot_not_required_for_this_forecast=True,
            joint_classical_record_and_reference_bound=True,
            corrected_generic_summary_requirement=True,
            full_coherent_record_channel_reproduced=False,
            arbitrary_unknown_local_state_classically_reconstructed=False,
            autonomous_forecaster_built=False,certified_nontrivial_spatial_resolution=False,
            physical_three_dimensional_space_generated=False,full_GR_goal_completed=False),
        observations=OBS)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();answer=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(answer,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:answer[k] for k in ('scientific_baseline_round','diagnostic_tests','failures','errors')}))
