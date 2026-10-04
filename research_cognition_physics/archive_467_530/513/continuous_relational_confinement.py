"""Round 513: finite continuous protection of relative memory sectors.

The size-uniform theorem reuses the old local NNI commutator bound 48.
X labels are conserved inputs, not dynamically created memberships.
"""
import argparse
from fractions import Fraction as Q
from functools import lru_cache
import hashlib
import io
import json
from pathlib import Path
import platform
import unittest

import numpy as np
import coherent_graph_mean_obstruction as old
from distributed_role_reader import exponential, opnorm
import historical_membership_interface as history

HERE = Path(__file__).resolve().parent
TARGET = HERE/'continuous_relational_confinement_results.json'
OBS = {}


@lru_cache(None)
def blocks():
    a = history.model()
    hf = old.model(2)[3]
    out = []
    for sign in a['signs']:
        c = np.repeat(sign, a['G'])
        same = c[:,None] == c[None,:]
        hd = np.where(same,a['h'],0)
        ho = a['h']-hd
        out.append((c,hd,ho))
    return a,hf,out


def budget(t,g):
    assert g != 0
    return (Q(3)+Q(93,2)*abs(t))/abs(g)


@lru_cache(None)
def propagators(t,g):
    a,_,data = blocks()
    original,comparison = [],[]
    for c,hd,_ in data:
        original.append(exponential(float(t)*(a['h']+float(g)*np.diag(c))))
        comparison.append(exponential(float(t)*(hd+float(g)*np.diag(c))))
    return np.array(original),np.array(comparison)


def readout(t,g,anchor):
    a,_,data = blocks()
    us,_ = propagators(t,g)
    columns = slice(anchor*a['G'],(anchor+1)*a['G'])
    branches = []
    for s,u in enumerate(us):
        for v in range(a['n']):
            k = u[v*a['G']:(v+1)*a['G'],columns]/np.sqrt(a['modes'])
            branches.append((s,v,k))
    return branches


class Audit(unittest.TestCase):
    def test_01_exact_decomposition_and_intertwiner_residual(self):
        a,_,data = blocks()
        g=8
        largest = 0
        for c,hd,ho in data:
            self.assertTrue(np.all(c*c==1))
            self.assertTrue(np.array_equal(c[:,None]*hd,hd*c[None,:]))
            self.assertTrue(np.array_equal(c[:,None]*ho,-ho*c[None,:]))
            hh=a['h']+g*np.diag(c)
            h0=hd+g*np.diag(c)
            scaled_w=2*g*np.eye(a['d'],dtype=np.int64)-c[:,None]*ho
            expected=c[:,None]*(ho@ho-(hd@ho-ho@hd))
            residual=hh@scaled_w-scaled_w@h0-expected
            largest=max(largest,int(np.max(np.abs(residual))))
        self.assertEqual(largest,0)
        OBS['exact_intertwiner'] = dict(memory_sectors=a['modes'],block_dimension=a['d'],
            integer_residual_max=largest,comparison_W_is_mathematical_not_applied=True,
            H_form_unchanged_from_510=True)

    def test_02_local_bounds_without_global_graph_norm(self):
        a,hf,data = blocks()
        maxima = dict(off_row_sum=0,graph_commutator_row_sum=0,diagonal_commutator_row_sum=0)
        for c,hd,ho in data:
            off = int(np.max(np.sum(abs(ho),axis=1)))
            fcomm=hf@ho-ho@hf
            dcomm=hd@ho-ho@hd
            maxima['off_row_sum']=max(maxima['off_row_sum'],off)
            maxima['graph_commutator_row_sum']=max(maxima['graph_commutator_row_sum'],int(np.max(np.sum(abs(fcomm),axis=1))))
            maxima['diagonal_commutator_row_sum']=max(maxima['diagonal_commutator_row_sum'],int(np.max(np.sum(abs(dcomm),axis=1))))
            self.assertLessEqual(off,3)
            self.assertLessEqual(np.max(np.sum(abs(fcomm),axis=1)),48)
            self.assertLessEqual(np.max(np.sum(abs(dcomm),axis=1)),84)
            self.assertLessEqual(opnorm(hd-hf),6+1e-12)
        OBS['local_bound_crosscheck']=dict(observed=maxima,
            theorem_off_norm=3,theorem_F_commutator=48,theorem_hd_commutator=84,
            global_F_norm_assumed_uniform=False,
            arbitrary_size_step_is_analytic_reuse_of_round496=True)

    def test_03_full_propagator_error_and_sector_leakage(self):
        a,_,data = blocks()
        rows=[]
        for t,g in ((Q(1,1024),Q(65536)),(Q(1),Q(1024))):
            us,vs=propagators(t,g)
            error=max(opnorm(u-v) for u,v in zip(us,vs))
            leakage=0.
            for anchor in (0,2):
                for s,(u,v) in enumerate(zip(us,vs)):
                    mask=data[s][0]==a['signs'][s,anchor]
                    self.assertLess(np.linalg.norm(v[np.ix_(~mask,mask)]),1e-12)
                    if np.any(~mask):
                        leakage=max(leakage,opnorm(u[np.ix_(~mask,mask)]))
            bound=budget(t,g)
            self.assertLess(error,float(bound))
            self.assertLess(leakage,float(bound))
            rows.append(dict(time=str(t),coupling=str(g),operator_bound=str(bound),
                all_input_operator_error=old.short(error),maximum_sector_leak_amplitude=old.short(leakage),
                leakage_probability_bound=str(min(Q(1),bound**2))))
        OBS['finite_protection']=rows

    def test_04_nontrivial_propagation_within_relative_class(self):
        a,_,data=blocks()
        t,g=Q(1,1024),Q(65536)
        anchor=2
        columns=slice(anchor*a['G'],(anchor+1)*a['G'])
        first_gram=np.zeros((a['G'],a['G']),dtype=np.int64)
        for c,hd,ho in data:
            # Source is a leaf; all first-step destinations are internal.
            b=hd[:2*a['G'],columns]
            first_gram+=b.T@b
        self.assertTrue(np.array_equal(first_gram,32*np.eye(a['G'],dtype=np.int64)))
        lower_amplitude=Q(2,3)*t-48*t*t-budget(t,g)
        self.assertGreater(lower_amplitude,0)
        effect=np.zeros((a['G'],a['G']),complex)
        for s,v,k in readout(t,g,anchor):
            if v<2 and a['signs'][s,v]==a['signs'][s,anchor]:
                effect+=k.conj().T@k
        eigen=np.linalg.eigvalsh(effect)
        self.assertGreater(eigen[0],float(lower_amplitude**2))
        OBS['nontrivial_same_class_reception']=dict(anchor=anchor,
            exact_uniform_memory_first_order_gram='(1/2) I_G',
            time=str(t),coupling=str(g),strict_amplitude_lower=str(lower_amplitude),
            strict_probability_lower=str(lower_amplitude**2),
            actual_effect_eigenvalue_min=old.short(eigen[0]),
            actual_effect_eigenvalue_max=old.short(eigen[-1]),
            unknown_graph_reference_covered=True,
            memory_source_is_independent_Z_vacuum=True)

    def test_05_conserved_source_and_nonspatial_partition_counterexample(self):
        a,_,data=blocks()
        anchor=2
        size=[int(np.sum(sign==sign[anchor])) for sign in a['signs']]
        counts={str(k):size.count(k) for k in range(1,7)}
        self.assertEqual(list(counts.values()),[2,10,20,20,10,2])
        # Four leaves same as anchor, both internal vertices opposite: not connected.
        sign=np.array([-1,-1,1,1,1,1])
        allowed={v for v in range(6) if sign[v]==sign[anchor]}
        self.assertEqual(allowed,{2,3,4,5})
        self.assertTrue(all(not any(u in allowed and v in allowed for u,v in tree) for tree in a['trees']))
        # Two degenerate classes: all sites, or just the anchor, are permitted too.
        self.assertTrue(any(k==1 for k in size))
        self.assertTrue(any(k==6 for k in size))
        OBS['source_and_geometry_boundary']=dict(class_size_histogram=counts,
            total_sign_strings=64,expected_anchor_class_size='7/2',
            disconnected_four_leaf_partition_exact=True,
            all_X_labels_commute_with_both_generators=True,
            partition_dynamically_created_claimed=False,
            same_statistics_after_X_dephasing_for_X_block_menu=True,
            original_Z_history_labels_frozen_claimed=False)

    def test_06_complete_actual_record_instrument(self):
        a,_,_=blocks()
        t,g,anchor=Q(1,1024),Q(65536),2
        effects=[np.zeros((a['G'],a['G']),complex) for _ in range(2)]
        sign_effects=[np.zeros((a['G'],a['G']),complex) for _ in range(a['modes'])]
        total=np.zeros_like(effects[0])
        branch_count=0
        for s,v,k in readout(t,g,anchor):
            value=k.conj().T@k
            total+=value
            sign_effects[s]+=value
            same=int(a['signs'][s,v]==a['signs'][s,anchor])
            effects[same]+=value
            branch_count+=1
        error=opnorm(total-np.eye(a['G']))
        self.assertLess(error,2e-12)
        opposite_max=float(np.linalg.eigvalsh(effects[0])[-1])
        self.assertLess(opposite_max,float(budget(t,g)**2))
        source_error=max(opnorm(e-np.eye(a['G'])/a['modes']) for e in sign_effects)
        self.assertLess(source_error,1e-13)
        OBS['actual_readout']=dict(final_full_X_strings_and_port_retained=True,
            Kraus_branches=branch_count,input_graph_columns=a['G'],
            completeness_error=old.short(error),
            maximum_opposite_class_probability=old.short(opposite_max),
            X_record_alone_independent_of_unknown_graph_effect_error=old.short(source_error),
            one_final_readout_with_no_repeated_projection=True,
            reader_clock_and_record_comparison_remain_inputs=True)


def run():
    OBS.clear()
    stream=io.StringIO()
    result=unittest.TextTestRunner(stream=stream,verbosity=0).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise AssertionError(stream.getvalue())
    return dict(date='2026-09-28',round=513,scientific_baseline_round=512,
        tests_run=result.testsRun,failures=len(result.failures),errors=len(result.errors),
        python=platform.python_version(),numpy=np.__version__,
        dependency_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in (
            'historical_membership_interface.py','coherent_graph_mean_obstruction.py',
            'research_note_496.md','research_note_510.md','research_note_512.md')},
        scope=dict(same_continuous_510_generator=True,
            finite_time_size_uniform_relative_sector_bound=True,
            all_unknown_memory_graph_and_reference_for_protection=True,
            nontrivial_actual_same_class_reception=True,
            conserved_label_source_boundary_proved=True,
            original_Z_history_membership_permanently_frozen=False,
            membership_partition_dynamically_created=False,
            connected_macroscopic_spatial_regions_generated=False,
            infinite_time_protection_at_fixed_g=False,
            autonomous_readout_and_controller_generated=False,
            dimension_three_generated=False,full_GR_goal_completed=False,
            phase_closure_triggered=False),observations=OBS)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf-8',newline='\n') as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('round','tests_run','failures','errors')},ensure_ascii=False))
