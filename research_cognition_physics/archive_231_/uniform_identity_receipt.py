"""Round 503: the same identity marking Hamiltonian for every root.

All sites mark, including the querying site. Self labels are explicit failures.
The role register of round 501 is absent, while identity and control inputs remain.
"""
import argparse
from fractions import Fraction as F
import hashlib
import io
import itertools
import json
from pathlib import Path
import platform
import unittest
import numpy as np
import local_neighbor_receipt as old501
import local_return_recycling as old502
import current_neighbor_detection as old500
import branching_tree_distance_audit as graphs

HERE = Path(__file__).resolve().parent
TARGET = HERE/'uniform_identity_receipt_results.json'
OBS = {}


def model(trees, flip, n):
    """No root argument: port, graph, query or any of N reply colours."""
    ng, nc = len(trees), n+1
    base = np.kron(np.eye(n, dtype=np.int64), flip)
    for a, tree in enumerate(trees):
        for u, v in tree:
            x, y = u*ng+a, v*ng+a
            base[x, y] += 1; base[y, x] += 1
            base[x, x] -= 1; base[y, y] -= 1
    h = np.kron(base, np.eye(nc, dtype=np.int64))
    for v in range(n):
        for a in range(ng):
            q = (v*ng+a)*nc
            h[q, q+v+1] += 1; h[q+v+1, q] += 1
    return h


def root_interface(trees, n, root):
    ng, nc = len(trees), n+1
    source = (root*ng+np.arange(ng))*nc
    success = np.zeros(n*ng*nc, bool)
    bad = success.copy(); self_reply = success.copy()
    for a, tree in enumerate(trees):
        self_reply[(root*ng+a)*nc+root+1] = True
        for label in range(n):
            if label == root:
                continue
            row = (root*ng+a)*nc+label+1
            success[row] = True
            bad[row] = tuple(sorted((root, label))) not in tree
    return source, success, bad, self_reply


def cycle_blocks(u, source, n, ng, root, cap):
    nc = n+1
    away = np.ones(len(u), bool)
    away[root*ng*nc:(root+1)*ng*nc] = False
    survivor = np.eye(len(u), dtype=complex)[:, source]
    accepted, reset = [], []
    for k in range(1, cap+1):
        evolved = u@survivor
        for colour in range(nc):
            block = evolved[(root*ng+np.arange(ng))*nc+colour]
            is_accept = k == 1 and colour not in (0, root+1)
            (accepted if is_accept else reset).append((k, colour, block))
        survivor = evolved*away[:, None]
    return accepted, reset, survivor


class Checks(unittest.TestCase):
    def test_01_full_tensor_without_roles(self):
        n, d = 3, 5
        trees = sorted({graphs.prufer_tree(n, (v,)) for v in range(n)}, key=lambda g:sorted(g))
        ng = len(trees)
        h = model(trees, np.zeros((ng, ng), dtype=np.int64), n)
        words = list(itertools.product(range(d), repeat=n))
        index = {w:k for k,w in enumerate(words)}
        full = np.zeros((len(words)*ng, len(words)*ng), dtype=np.int64)
        for col, word in enumerate(words):
            for a, tree in enumerate(trees):
                for u,v in tree:
                    target=list(word); target[u],target[v]=target[v],target[u]
                    full[index[tuple(target)]*ng+a,col*ng+a]+=1
                for v in range(n):
                    if word[v] in (1, v+2):
                        target=list(word); target[v]=v+2 if word[v]==1 else 1
                        full[index[tuple(target)]*ng+a,col*ng+a]+=1
        selected=[]
        for v in range(n):
            for a in range(ng):
                for colour in range(n+1):
                    word=[0]*n; word[v]=colour+1
                    selected.append(index[tuple(word)]*ng+a)
        outside=sorted(set(range(len(full)))-set(selected))
        self.assertFalse(np.any(full[np.ix_(outside,selected)]))
        self.assertTrue(np.array_equal(full[np.ix_(selected,selected)]-2*np.eye(len(h),dtype=np.int64),h))
        OBS['full_tensor']=dict(N=n,full_dimension=len(full),active_dimension=len(h),
            all_N_colours_retained=True,role_registers=0,exact_restriction=True)

    def test_02_all_roots_same_generator(self):
        trees,_,_,flip=old500.six_model()
        n,ng,nc=6,len(trees),7
        h=model(trees,flip,n)
        fg=np.kron(np.kron(np.eye(n,dtype=np.int64),flip),np.eye(nc,dtype=np.int64))
        comm=fg@h-h@fg
        self.assertLessEqual(np.max(np.sum(abs(comm),axis=1)),48)
        cases=[]
        for root in range(n):
            source,success,bad,self_reply=root_interface(trees,n,root)
            power=np.eye(len(h),dtype=np.int64)[:,source]
            degree=np.array([sum(root in edge for edge in tree) for tree in trees])
            expected=np.zeros_like(power)
            for a,tree in enumerate(trees):
                for label in range(n):
                    if tuple(sorted((root,label))) in tree:
                        expected[(root*ng+a)*nc+label+1,a]=1
            for order in range(4):
                if order < 3:
                    self.assertFalse(np.any(power[success]))
                else:
                    self.assertTrue(np.array_equal(power*success[:,None],expected))
                    self.assertTrue(np.array_equal(expected.T@expected,np.diag(degree)))
                    self.assertFalse(np.any(power[bad]))
                if order == 1:
                    self.assertTrue(np.array_equal(power[self_reply],np.eye(ng,dtype=np.int64)))
                power=h@power
            cases.append(dict(root=root,degrees=sorted(set(map(int,degree))),
                first_self_response_exact=True,third_nonself_response_exact=True))
        OBS['root_independence']=dict(active_dimension=len(h),one_H_for_all_roots=True,
            root_dependent_H_terms=0,commutator_row_bound=int(np.max(np.sum(abs(comm),axis=1))),cases=cases)

    def test_03_all_root_rational_certificates(self):
        trees,_,_,flip=old500.six_model()
        h=model(trees,flip,6)
        certs=[]
        for root in range(6):
            source,success,bad,_=root_interface(trees,6,root)
            data=old501.source_taylor(h,source)
            low,high=old501.effect_bound(data,success)
            _,wrong=old501.effect_bound(data,bad)
            self.assertGreater(low,F(1,10**13))
            self.assertLess(wrong/low,F(1,10000))
            certs.append(dict(root=root,success_lower=str(low),bad_mass_upper=str(wrong),
                joint_error_upper=str(wrong/low),readable_success_lower=float(low),
                readable_joint_error_upper=float(wrong/low),norm_bound=data[-1],tail=str(data[-2])))
        OBS['rational_all_root_certificates']=dict(time='1/64',Taylor_degree=12,
            uniform_finite_success_lower='1/10000000000000',uniform_finite_error_upper='1/10000',
            all_graph_inputs_and_passive_references=True,certificates=certs)

    def test_04_complete_recycling_at_two_roots(self):
        trees,_,_,flip=old500.six_model()
        h=model(trees,flip,6); u=old500.evolution(h,F(1,64))
        psi=np.array([[1,1j],[2j,-1],[1+1j,2],[-1j,1],[2,-2j],[1,3j]],complex)
        psi/=np.linalg.norm(psi)
        initial=np.outer(psi.ravel(),psi.ravel().conj())
        diagnostics=[]
        for root in [0,1]:
            source,_,_,_=root_interface(trees,6,root)
            fresh,reset,timeout=cycle_blocks(u,source,6,len(trees),root,8)
            gram=timeout.conj().T@timeout
            for _,_,k in fresh+reset:gram+=k.conj().T@k
            self.assertLess(np.linalg.norm(gram-np.eye(6)),2e-12)
            active=initial.copy();reference=np.zeros((2,2),complex)
            accept=bad=timed=0.0
            for cycle in range(2):
                rho4=active.reshape(6,2,6,2)
                for _,colour,k in fresh:
                    kr=np.kron(k,np.eye(2));out=kr@active@kr.conj().T
                    accept+=float(np.trace(out).real)
                    mask=np.repeat([tuple(sorted((root,colour-1))) not in t for t in trees],2)
                    bad+=float(np.trace(out[np.ix_(mask,mask)]).real)
                    reference+=old502.reference_after(k,rho4)
                rt=old502.reference_after(timeout,rho4);timed+=float(np.trace(rt).real);reference+=rt
                following=np.zeros_like(active)
                for _,_,k in reset:
                    kr=np.kron(k,np.eye(2));following+=kr@active@kr.conj().T
                active=following
            remaining=float(np.trace(active).real)
            reference+=np.einsum('grgs->rs',active.reshape(6,2,6,2))
            self.assertAlmostEqual(accept+timed+remaining,1,places=12)
            self.assertLess(np.linalg.norm(reference-psi.T@psi.conj()),2e-12)
            self.assertLess(bad/accept,1/10000)
            diagnostics.append(dict(root=root,accepted=accept,timeout=timed,exhausted=remaining,
                conditional_current_bad=bad/accept,self_return_in_reset=True,
                complete_residual=float(np.linalg.norm(gram-np.eye(6))),
                reference_residual=float(np.linalg.norm(reference-psi.T@psi.conj()))))
        OBS['complete_recycling']=dict(same_matrix_used=True,cycles=2,cap=8,diagnostics=diagnostics)

    def test_05_self_reply_must_be_excluded(self):
        trees=[frozenset({(0,1)})]
        h=model(trees,np.zeros((1,1),dtype=np.int64),2)
        source,success,_,own=root_interface(trees,2,0)
        data=old501.source_taylor(h,source)
        own_low,_=old501.effect_bound(data,own)
        _,nonself_high=old501.effect_bound(data,success)
        self.assertGreater(own_low/(own_low+nonself_high),F(99,100))
        OBS['self_reply_counterexample']=dict(N=2,time='1/64',
            self_mass_lower=str(own_low),nonself_mass_upper=str(nonself_high),
            self_fraction_among_all_replies_strict_lower='99/100',
            physical_graph_has_no_self_edges=True,actual_instrument_rejects_self=True)

    def test_06_relabel_covariance_and_resource_change(self):
        n=3
        trees=sorted({graphs.prufer_tree(n,(v,)) for v in range(n)},key=lambda g:sorted(g))
        h=model(trees,np.zeros((3,3),dtype=np.int64),n)
        for perm in itertools.permutations(range(n)):
            mapping=[]
            for v in range(n):
                for tree in trees:
                    mapped=frozenset(tuple(sorted((perm[a],perm[b]))) for a,b in tree)
                    gi=trees.index(mapped)
                    for colour in range(n+1):
                        c=0 if colour==0 else perm[colour-1]+1
                        mapping.append((perm[v]*3+gi)*4+c)
            self.assertTrue(np.array_equal(h[np.ix_(mapping,mapping)],h))
        OBS['scope_resources']=dict(exact_joint_relabellings=6,
            N6_active_dimension=252,N6_previous_active_dimension=216,
            N6_listed_register_qubits=72,N6_previous_listed_register_qubits=78,
            root_role_qubits_removed=6,all_data_identity_graph_registers_retained=True,
            root_readout_reachable_outcomes=8,root_readout_bits=3,
            different_roots_share_dynamics=True,simultaneous_multiple_packets_certified=False,
            local_self_comparison_uses_prepared_identity_dictionary=True)


def run():
    OBS.clear()
    stream=io.StringIO()
    r=unittest.TextTestRunner(stream=stream,verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Checks))
    if not r.wasSuccessful():raise AssertionError(stream.getvalue())
    deps=['local_neighbor_receipt_results.json','local_return_recycling_results.json']
    return dict(round=503,scientific_baseline_round=502,tests_run=r.testsRun,
        failures=len(r.failures),errors=len(r.errors),python=platform.python_version(),numpy=np.__version__,
        dependency_results_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in deps},
        scope=dict(global_root_role_configuration_required=False,all_sites_always_mark=True,
            same_fixed_H_for_every_query_root=True,one_packet_source_still_required=True,
            trusted_identity_and_marking_interaction_still_inputs=True,
            controls_and_clocks_generated=False,multi_packet_concurrency_proved=False,
            current_random_label_joint_bound_not_per_label=True,
            actual_position_displacements_generated=False,full_GR_goal_completed=False,
            phase_closure_triggered=False),observations=OBS.copy())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--write-results',action='store_true');args=p.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf-8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False,indent=2))
