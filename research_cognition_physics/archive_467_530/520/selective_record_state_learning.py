"""Round 520: retained port records learn the current finite relation state.

The all-size/all-input statement is analytic. Small matrices independently check
its word certificate and the full-reference report bound. Trajectories at a
larger diagnostic step do not certify the theorem's very small step or a rate.
"""
import argparse
from collections import Counter
from fractions import Fraction as Q
import hashlib
import io
import json
import math
from pathlib import Path
import platform
import unittest

import numpy as np
import all_scale_monitored_graph_source as source

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'selective_record_state_learning_results.json'
OBS = {}


def digest(name):
    return hashlib.sha256((HERE/name).read_bytes()).hexdigest()


def euler_walk(tree, n, root):
    adj = source.old.adjacency(tree, n)
    walk = [root]
    def visit(v, parent):
        for w in sorted(adj[v]):
            if w != parent:
                walk.append(w)
                visit(w, v)
                walk.append(v)
    visit(root, -1)
    return walk


def model(i, j=1, kappa=1):
    n, trees, flips = source.family(i)
    m = len(trees)
    h = np.zeros((n*m, n*m), dtype=np.int64)
    for g, tree in enumerate(trees):
        adj = source.old.adjacency(tree, n)
        for a in range(n):
            h[a*m+g, a*m+g] -= j*len(adj[a])
            for b in adj[a]:
                h[b*m+g, a*m+g] += j
            for gp, weight in flips[g].items():
                h[a*m+gp, a*m+g] += kappa*weight
    return n, trees, h


def taylor_unitary(h, time, order=32):
    term = np.eye(len(h), dtype=complex)
    result = term.copy()
    for k in range(1, order+1):
        term = (-1j*float(time)/k)*(h@term)
        result += term
    return result


def budget(i, j=Q(1), kappa=Q(1)):
    n = 2*i+2
    m = math.factorial(2*i)//2**i
    d, length = n*m, 2*(n-1)
    b = 4*abs(kappa)*(i-1)+6*abs(j)
    step = abs(j)/(16*d*length*b*b)
    x = b*b*step/abs(j)
    e = (1+x)**length-1
    effect = e*(2+e)
    return n, m, d, length, b, step, x, e, effect


def half_trace(x):
    return float(np.sum(np.abs(np.linalg.eigvalsh((x+x.conj().T)/2)))/2)


class Audit(unittest.TestCase):
    def test_01_all_root_tree_words_isolate_a_basis_state(self):
        rows = []
        for i in (1, 2, 3, 4):
            n, trees, _ = source.family(i)
            count = 0
            for g, tree in enumerate(trees):
                for a in range(n):
                    walk = euler_walk(tree, n, a)
                    counts = Counter(tuple(sorted((u,v))) for u,v in zip(walk,walk[1:]))
                    self.assertEqual(len(walk), 2*(n-1)+1)
                    self.assertEqual((walk[0],walk[-1]), (a,a))
                    self.assertEqual(set(counts), set(tree))
                    self.assertEqual(set(counts.values()), {2})
                    # For small families, explicitly form the product of all
                    # commuting edge effects on every graph basis vector.
                    if i <= 3:
                        surviving = [gp for gp,t in enumerate(trees)
                                     if all(e in t for e in counts)]
                        self.assertEqual(surviving, [g])
                    count += 1
            rows.append(dict(internal=i, vertices=n, graphs=len(trees), words=count,
                             length=2*(n-1), all_edges_twice=True,
                             explicit_all_graph_product=i<=3))
        OBS['integer_word_certificates'] = rows

    def test_02_one_step_budget_covers_all_words_and_projections(self):
        rows = []
        for i in range(1,13):
            for j,kappa in ((Q(1),Q(1)),(Q(2),Q(-1)),(Q(1),Q(0))):
                n,m,d,length,b,step,x,e,effect = budget(i,j,kappa)
                self.assertLessEqual(b*step, 1)
                self.assertEqual(length*x, Q(1,16*d))
                self.assertLessEqual(e, 2*length*x)
                self.assertLessEqual(e, Q(1,8*d))
                self.assertLessEqual(effect, Q(3,8*d))
                self.assertLess(effect, Q(1,d))
                # For every rank r>=2, q < r/(2d), without enumerating
                # uncountably many subspaces or testing floating ranks.
                self.assertLess(effect, Q(2,2*d))
                if j == kappa == 1 and i <= 4:
                    rows.append(dict(internal=i, dimension=d,
                                     step=str(step), word_length=length+1,
                                     effect_bound_upper=str(Q(3,8*d)),
                                     dark_projection_lower=str(Q(1,d))))
        OBS['exact_common_step_budgets'] = rows
        OBS['budget_cases_checked'] = 36

    def test_03_actual_kraus_words_and_rotated_basis(self):
        n,trees,h = model(2)
        m,d = len(trees),len(h)
        _,_,_,length,b,step,x,e,effect = budget(2)
        u = taylor_unitary(h,step)
        self.assertLess(np.linalg.norm(u.conj().T@u-np.eye(d)),1e-13)
        maximum = 0.0
        complete = np.zeros((d,d),complex)
        for a in range(n):
            first = u[a*m:(a+1)*m,:]
            for g,tree in enumerate(trees):
                walk=euler_walk(tree,n,a)
                prod=np.eye(m,dtype=complex)
                for s,t in zip(walk,walk[1:]):
                    block=u[t*m:(t+1)*m,s*m:(s+1)*m]/(-1j*float(step))
                    prod=block@prod
                v=prod@first
                effect_actual=v.conj().T@v
                row=first[g,:]
                projector=np.outer(row.conj(),row)
                complete+=projector
                error=float(np.linalg.norm(effect_actual-projector,2))
                maximum=max(maximum,error)
                self.assertLessEqual(error,float(effect)+1e-12)
        self.assertLess(np.linalg.norm(complete-np.eye(d)),1e-12)
        OBS['actual_word_effects']=dict(words=n*m, step=str(step),
            max_effect_error=maximum, analytic_effect_bound=float(effect),
            basis_completeness_error=float(np.linalg.norm(complete-np.eye(d))),
            taylor_order=32, b_step=float(b*step))

    def test_04_complete_history_and_unknown_reference_report(self):
        n,trees,h=model(2)
        m,d=len(trees),len(h)
        time=Q(1,8)
        u=taylor_unitary(h,time)
        rng=np.random.default_rng(520)
        psi=rng.normal(size=(d,2))+1j*rng.normal(size=(d,2))
        psi/=np.linalg.norm(psi)
        # A branch stores V_w : full initial DG -> current graph block.
        branches=[(a,u[a*m:(a+1)*m,:]) for a in range(n)]
        rows=[]
        for depth in range(1,5):
            mass_q=mass_p=impurity=report_error=0.0
            operator_sum=np.zeros((d,d),complex)
            max_dominance_error=0.0
            for a,v in branches:
                xx=v@v.conj().T/d
                q=float(np.trace(xx).real)
                ev=np.linalg.eigvalsh(xx)
                impurity+=q-float(ev[-1])
                sigma=xx/q
                z=v@psi
                p=float(np.vdot(z,z).real)
                rr=z.T@z.conj()
                actual=np.outer(z.reshape(-1),z.reshape(-1).conj())
                predicted=np.kron(sigma,rr)
                report_error+=half_trace(actual-predicted)
                operator_sum+=v.conj().T@v
                bad=float(max(0,-np.linalg.eigvalsh(d*xx-z@z.conj().T)[0]))
                max_dominance_error=max(max_dominance_error,bad)
                mass_q+=q;mass_p+=p
            self.assertAlmostEqual(mass_q,1,places=11)
            self.assertAlmostEqual(mass_p,1,places=11)
            self.assertLess(np.linalg.norm(operator_sum-np.eye(d)),1e-10)
            bound=2*math.sqrt(d*max(impurity,0))+d*impurity
            self.assertLessEqual(report_error,bound+1e-11)
            self.assertLess(max_dominance_error,1e-11)
            rows.append(dict(depth=depth,histories=len(branches),m_n=impurity,
                total_actual_probability=mass_p,total_estimator_probability=mass_q,
                full_reference_half_trace_error=report_error,
                analytic_bound_unclipped=bound,diagnostic_time=str(time)))
            if depth<4:
                branches=[(b,u[b*m:(b+1)*m,a*m:(a+1)*m]@v)
                          for a,v in branches for b in range(n)]
        OBS['finite_history_full_reference_checks']=rows

    def test_05_mixing_and_purity_do_not_mean_a_classical_fixed_graph(self):
        eye=np.eye(2,dtype=complex)
        x=np.array([[0,1],[1,0]],complex)
        y=np.array([[0,-1j],[1j,0]],complex)
        z=np.diag([1,-1]).astype(complex)
        rho=np.diag([0.8,0.2]).astype(complex)
        outputs=[a@rho@a.conj().T/4 for a in (eye,x,y,z)]
        self.assertLess(np.linalg.norm(sum(outputs)-eye/2),1e-14)
        self.assertTrue(all(abs(np.trace((4*r)@(4*r)).real-0.68)<1e-14
                            for r in outputs))
        n,trees,h=model(2)
        m=len(trees)
        # A pure graph superposition, under one actual port result, remains
        # pure but need not become any graph basis vector.
        u=taylor_unitary(h,Q(1,8))
        g=np.zeros(m,complex);g[0]=g[1]=1/math.sqrt(2)
        out=u[0:m,0:m]@g
        prob=float(np.vdot(out,out).real)
        out/=math.sqrt(prob)
        pop=np.abs(out)**2
        self.assertLess(float(np.max(pop)),0.9)
        OBS['scope_counterexamples']=dict(
            primitive_pauli_channel_posterior_purity=0.68,
            primitive_channel_alone_does_not_imply_purification=True,
            actual_pure_graph_no_jump_probability=prob,
            largest_graph_basis_probability=float(np.max(pop)),
            pure_posterior_does_not_mean_classical_graph=True)

    def test_06_nonvacuous_reference_bound_for_actual_rare_words(self):
        n,trees,h=model(2)
        m,d=len(trees),len(h)
        _,_,_,length,_,step,_,word_error,_=budget(2)
        # Singular-value perturbation: every non-leading singular value of
        # B is <= word_error; its largest is >= 1-word_error.
        eta_upper=Q(m-1)*word_error**2/(1-word_error)**2
        u=taylor_unitary(h,step)
        rng=np.random.default_rng(52006)
        psi=rng.normal(size=(d,2))+1j*rng.normal(size=(d,2))
        psi/=np.linalg.norm(psi)
        errors=[];bounds=[];log_probs=[]
        for a in range(n):
            for tree in trees:
                walk=euler_walk(tree,n,a)
                v=u[a*m:(a+1)*m,:].copy()
                for s,t in zip(walk,walk[1:]):
                    v=(u[t*m:(t+1)*m,s*m:(s+1)*m]/(-1j*float(step)))@v
                xx=v@v.conj().T/d
                q=float(np.trace(xx).real)
                sigma=xx/q
                z=v@psi
                p=float(np.vdot(z,z).real)
                normalized=z/math.sqrt(p)
                true=np.outer(normalized.reshape(-1),normalized.reshape(-1).conj())
                rr=normalized.T@normalized.conj()
                error=half_trace(true-np.kron(sigma,rr))
                bound=2*math.sqrt(d*q*float(eta_upper)/p)+float(eta_upper)
                self.assertLessEqual(error,bound+1e-12)
                self.assertLess(bound,1)
                errors.append(error);bounds.append(bound)
                log_probs.append(2*length*math.log10(float(step))+math.log10(p))
        OBS['rare_word_reference_diagnostic']=dict(
            words=len(errors),step=str(step),max_conditioned_error=max(errors),
            max_conditioned_bound=max(bounds),
            analytic_purity_defect_upper=float(eta_upper),
            log10_actual_probability_min=min(log_probs),
            log10_actual_probability_max=max(log_probs),
            all_probabilities_retained=True,
            not_a_typical_purification_rate=True)


def run():
    OBS.clear()
    suite=unittest.defaultTestLoader.loadTestsFromTestCase(Audit)
    stream=io.StringIO()
    result=unittest.TextTestRunner(stream=stream,verbosity=2).run(suite)
    if not result.wasSuccessful():
        raise AssertionError(stream.getvalue())
    return dict(date='2026-09-28',round=520,scientific_base_through_round=519,
        tests_run=result.testsRun,failures=len(result.failures),errors=len(result.errors),
        python=platform.python_version(),numpy=np.__version__,
        dependency_sha256={name:digest(name) for name in (
            'all_scale_monitored_graph_source.py','research_note_517.md',
            'research_note_518.md','research_note_519.md',
            'uniform_tree_continuum_review.md')},
        scope=dict(same_retained_DG_and_original_H=True,
            fixed_finite_size_common_step_no_dark_subspace=True,
            all_unknown_initial_states_and_passive_references=True,
            retained_complete_classical_history_estimator=True,
            uniform_mean_full_reference_report_bound=True,
            graph_resets_or_new_target_dependent_measurements=False,
            exact_autonomous_reader_or_computer_generated=False,
            uniform_in_size_or_explicit_mixing_rate=False,
            permanently_fixed_classical_graph=False,
            macroscopic_positions_or_three_dimensions_derived=False,
            full_GR_goal_completed=False),observations=OBS.copy())


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dry-run',action='store_true')
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    answer=json.loads(json.dumps(run()))
    if args.check:
        assert answer==json.loads(TARGET.read_text(encoding='utf8'))
    elif not args.dry_run:
        with TARGET.open('x',encoding='utf8',newline='\n') as handle:
            json.dump(answer,handle,ensure_ascii=False,indent=2)
            handle.write('\n')
    print(json.dumps(answer,ensure_ascii=False,indent=2))
