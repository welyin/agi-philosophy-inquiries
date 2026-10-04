"""Unnumbered mature-theorem bridge after 519, not a new geometry experiment.

Check the exact finite sampling law before applying a CRT limit theorem.
No simulation or finite fitted exponent is used as proof of a dimension.
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

HERE=Path(__file__).resolve().parent
TARGET=HERE/'uniform_tree_continuum_probe_results.json'
OBS={}


def catalan(i):
    return math.comb(2*i,i)//(i+1)


def rooted_children(tree,n,a):
    adj=source.old.adjacency(tree,n)
    root=next(iter(adj[a]));children={}
    def visit(v,parent):
        children[v]=sorted(adj[v]-{parent})
        for w in children[v]:visit(w,v)
    visit(root,a)
    return root,children


def ordered_shape(root,children,mask):
    def visit(v):
        cs=children[v]
        if not cs:return '.'
        assert len(cs)==2
        if mask & (1<<v):cs=cs[::-1]
        return '('+visit(cs[0])+visit(cs[1])+')'
    return visit(root)


def unordered_shape(root,children):
    def visit(v):
        cs=children[v]
        if not cs:return '.'
        return '('+''.join(sorted(visit(w) for w in cs))+')'
    return visit(root)


def trace_distance(x,y):
    z=(x-y+(x-y).conj().T)/2
    return float(np.abs(np.linalg.eigvalsh(z)).sum()/2)


class Audit(unittest.TestCase):
    def test_01_exact_plane_sampling_law(self):
        rows=[]
        for i in (1,2,3,4):
            n,trees,_=source.family(i);a=i
            counts=Counter();unordered=Counter()
            for tree in trees:
                root,children=rooted_children(tree,n,a)
                self.assertEqual(len(children),2*i+1)
                self.assertEqual(sum(len(v)==2 for v in children.values()),i)
                unordered[unordered_shape(root,children)]+=1
                for mask in range(1<<i):
                    counts[ordered_shape(root,children,mask)]+=1
            self.assertEqual(len(counts),catalan(i))
            self.assertEqual(set(counts.values()),{math.factorial(i)*math.factorial(i+1)})
            self.assertEqual(len(trees)*(1<<i),math.factorial(2*i))
            if i==3:
                self.assertEqual(sorted(unordered.values()),[18,72])
            rows.append(dict(internal_vertices=i,original_graphs=len(trees),
                oriented_shapes=len(counts),multiplicity_per_plane_shape=next(iter(counts.values())),
                rooted_unordered_shape_counts=sorted(unordered.values()),
                exact_uniform_Catalan_law=True))
        OBS['finite_sampling_law']=rows

    def test_02_exact_depth_recurrence_and_scale_diagnostic(self):
        # S_i sums all vertex depths over every plane full binary tree.
        sums=[0]
        for i in range(1,33):
            val=sum(sums[k]*catalan(i-1-k)+catalan(k)*sums[i-1-k]
                +2*i*catalan(k)*catalan(i-1-k) for k in range(i))
            sums.append(val)
            self.assertEqual(val,2*4**i-2*(2*i+1)*catalan(i))
        self.assertEqual(Q(1,2)*0+Q(1,2)*2,1)
        self.assertEqual(Q(1,2)*(0-1)**2+Q(1,2)*(2-1)**2,1)
        rows=[]
        for i in (4,16,64,256,1024):
            mean=Q(2*4**i,(2*i+1)*catalan(i))-2
            rows.append(dict(internal_vertices=i,total_vertices_after_root_leaf_removal=2*i+1,
                mean_depth_over_sqrt_n=float(f'{float(mean)/math.sqrt(2*i+1):.10g}')))
        OBS['normalization_diagnostic']=dict(offspring_support=[0,2],mean=1,variance=1,
            allowed_total_sizes='odd only',exact_depth_recurrence_cases=32,values=rows,
            asymptotic_mean_depth_over_sqrt_n='sqrt(pi/2)',
            fitted_exponent_used_as_proof=False)

    def test_03_finite_source_error_pushforward_with_reference(self):
        n,trees,_=source.family(3);m=len(trees);a=3
        classes={}
        for g,tree in enumerate(trees):
            root,children=rooted_children(tree,n,a)
            classes.setdefault(unordered_shape(root,children),[]).append(g)
        # A correlated source perturbation; the reference marginal stays I/2.
        eps=Q(1,32)
        psi=np.zeros((m,2),complex);psi[0,0]=1/math.sqrt(2);psi[-1,1]=1j/math.sqrt(2)
        ideal=np.eye(2*m)/(2*m)
        rho=(1-float(eps))*ideal+float(eps)*np.outer(psi.reshape(-1),psi.reshape(-1).conj())
        inp=trace_distance(rho,ideal)
        out=0.0
        for indices in classes.values():
            block=np.zeros((2,2),complex)
            for g in indices:block+=rho[2*g:2*g+2,2*g:2*g+2]
            ref=len(indices)*np.eye(2)/(2*m)
            out+=trace_distance(block,ref)
        self.assertAlmostEqual(inp,float(eps)*(1-1/(2*m)),places=12)
        self.assertLessEqual(out,inp+1e-12)
        OBS['source_error_pushforward_diagnostic']=dict(graphs=m,reference_dimension=2,
            input_half_trace_distance=float(f'{inp:.10g}'),
            shape_reference_report_distance=float(f'{out:.10g}'),
            preserved_passive_reference=True,physical_joint_graph_reader_built=False)


def run():
    OBS.clear();stream=io.StringIO()
    result=unittest.TextTestRunner(stream=stream,verbosity=0).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():raise AssertionError(stream.getvalue())
    return dict(date='2026-09-28',scientific_baseline_round=519,
        numbered_round_created=False,numbered_scientific_test_increment=0,
        diagnostic_tests=result.testsRun,failures=len(result.failures),errors=len(result.errors),
        python=platform.python_version(),numpy=np.__version__,
        dependency_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in (
            'all_scale_monitored_graph_source.py','research_note_467.md','research_note_469.md',
            'research_note_491.md','research_note_517.md','research_note_519.md')},
        source_urls=dict(periodic_binary_limit='https://arxiv.org/pdf/math/0511515',
            metric_normalization='https://www.stats.ox.ac.uk/~goldschm/PIMSminicoursev2.pdf',
            dimension='https://arxiv.org/pdf/math/0501079'),
        scope=dict(exact_sampling_law_matched=True,mature_CRT_limit_applied=True,
            odd_size_subsequence_explicit=True,finite_source_error_transferred=True,
            new_numbered_dimension_theorem_claimed=False,
            measured_dynamic_propagation_metric_identified=False,
            graph_metric_readout_autonomously_generated=False,
            single_growing_universe_constructed=False,
            macroscopic_three_dimensional_space_refuted_in_general=False,
            full_GR_goal_completed=False),observations=OBS)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();answer=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(answer,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:answer[k] for k in ('scientific_baseline_round','diagnostic_tests','failures','errors')}))
