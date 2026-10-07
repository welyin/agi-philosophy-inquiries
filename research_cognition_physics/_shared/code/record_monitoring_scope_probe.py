"""Unnumbered follow-up to 514: check an existing monitoring corollary.

This is not a new spatial construction or a new discovery of the Zeno effect.
All readouts are complete classical instruments, with the initial pinch paid.
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
import historical_membership_interface as old
from distributed_role_reader import exponential, opnorm

HERE=Path(__file__).resolve().parent
TARGET=HERE/'record_monitoring_scope_probe_results.json'
OBS={}


def trace_rank_two(k,l):
    """Trace norm of normalized Choi branch difference; no diamond assertion."""
    d=k.shape[1]
    _,r=np.linalg.qr(np.column_stack((k.ravel(),l.ravel()))/np.sqrt(d))
    return float(np.abs(np.linalg.eigvalsh(r@np.diag([1.,-1.])@r.conj().T)).sum())


def symmetry_orbits_and_actions():
    a=old.model();total=a['modes']*a['d'];parent=list(range(total))
    tree_index={tuple(sorted(t)):i for i,t in enumerate(a['trees'])}
    def find(v):
        while parent[v]!=v:
            parent[v]=parent[parent[v]];v=parent[v]
        return v
    actions=[]
    for left,right in ((0,1),(2,3),(3,4),(4,5)):
        perm=list(range(a['n']));perm[left],perm[right]=perm[right],perm[left]
        gp=[tree_index[tuple(sorted(tuple(sorted((perm[u],perm[v]))) for u,v in t))]
            for t in a['trees']]
        mp=[sum(((bits>>v)&1)<<perm[v] for v in range(a['n'])) for bits in range(a['modes'])]
        action=np.array([mp[m]*a['d']+perm[v]*a['G']+gp[g]
            for m in range(a['modes']) for v in range(a['n']) for g in range(a['G'])])
        for j,k in enumerate(action):parent[find(j)]=find(int(k))
        actions.append(action)
    groups={}
    for j in range(total):groups.setdefault(find(j),[]).append(j)
    return list(groups.values()),actions


class Audit(unittest.TestCase):
    def test_01_complete_history_bound_without_enumerating_bad_records(self):
        a=old.model();time=Q(1,2);rows=[]
        for n in (8,16,32):
            delta=time/n
            blocks=old.history(delta)
            u=exponential(float(delta)*a['h'])
            k=blocks[0]
            step_return=opnorm(k-u)
            step_bad=np.eye(a['d'])-k.conj().T@k
            self.assertLess(step_return,float(delta**2/2)+1e-13)
            self.assertLess(np.linalg.eigvalsh(step_bad)[-1],float(delta**2)+1e-13)
            good=np.linalg.matrix_power(k,n)
            target=exponential(float(time)*a['h'])
            bad=np.eye(a['d'])-good.conj().T@good
            badmax=float(np.linalg.eigvalsh(bad)[-1])
            diamond_upper=2*opnorm(good-target)+badmax
            analytic=float(2*time*delta)
            self.assertLess(diamond_upper,analytic+2e-12)
            self.assertLess(badmax,float(time*delta)+1e-12)
            choi=trace_rank_two(good,target)+float(np.trace(bad).real/a['d'])
            rows.append(dict(steps=n,delta=str(delta),duration=str(time),
                step_return_error=old.old.short(step_return),
                all_changed_history_probability_max=old.old.short(badmax),
                complete_history_operator_diamond_upper=old.old.short(diamond_upper),
                complete_classical_record_choi_trace=old.old.short(choi),
                analytic_diamond_upper=str(2*time*delta),
                monitoring_classical_record_bits=n*a['n'],
                initial_memory_record_bits=a['n'],
                total_with_initial_record_bits=(n+1)*a['n']))
        OBS['complete_histories']=rows

    def test_02_actual_formation_then_member_read_and_monitoring(self):
        a=old.model();formation=old.history(Q(1,2))
        time=Q(1,2);n=16;delta=time/n
        k=np.linalg.matrix_power(old.history(delta)[0],n//2)
        u=exponential(float(time/2)*a['h'])
        good_effect=np.zeros((a['d'],a['d']),complex)
        ideal_effect=np.zeros_like(good_effect)
        choi=0.;multipart_effect=np.zeros_like(good_effect)
        for bits,source in enumerate(formation):
            if bits.bit_count()>=2:
                multipart_effect+=source.conj().T@source
            q=np.diag(a['q'][bits])
            for r in (0,1):
                projector=q if r else np.eye(a['d'])-q
                actual=k@projector@k@source
                target=u@projector@u@source
                good_effect+=actual.conj().T@actual
                ideal_effect+=target.conj().T@target
                choi+=trace_rank_two(actual,target)
        bad=np.eye(a['d'])-good_effect
        choi+=float(np.trace(bad).real/a['d'])
        self.assertLess(opnorm(ideal_effect-np.eye(a['d'])),1e-12)
        self.assertGreater(np.linalg.eigvalsh(bad)[0],-1e-12)
        self.assertLess(choi,float(2*time*delta))
        self.assertLess(np.linalg.eigvalsh(bad)[-1],float(time*delta)+1e-12)
        OBS['formed_membership_then_read']=dict(formation_time='1/2',steps=n,
            duration=str(time),endpoint_read_time=str(time/2),
            full_ideal_instrument_error=old.old.short(opnorm(ideal_effect-np.eye(a['d']))),
            source_two_or_more_labels_probability_min=old.old.short(np.linalg.eigvalsh(multipart_effect)[0]),
            complete_record_choi_trace=old.old.short(choi),
            analytic_diamond_upper=str(2*time*delta),
            changed_history_probability_max=old.old.short(np.linalg.eigvalsh(bad)[-1]),
            every_source_label_and_endpoint_outcome_retained=True,
            source_graph_not_reset=True,endpoint_projective_instrument_is_declared_input=True,
            finite_time_endpoint_reader_implementation_not_tested=True)

    def test_03_quantum_environment_scope_counterexample(self):
        delta=Q(1,64)
        # Fixed occupied site: coherent premeasurement has amplitude O(delta).
        s=np.sin(float(delta))
        coherent_distance=2*abs(s)
        classical_distance=2*s*s
        self.assertLess(classical_distance,float(2*delta**2))
        self.assertGreater(2*(delta-delta**3/6),2*delta**2)
        self.assertGreater(coherent_distance,float(2*delta**2)*50)
        OBS['readout_environment_boundary']=dict(delta=str(delta),
            joint_coherent_premeasurement_trace=old.old.short(coherent_distance),
            full_classical_instrument_trace=old.old.short(classical_distance),
            classical_quadratic_upper=str(2*delta**2),
            strict_sine_lower=str(delta-delta**3/6),
            coherent_environment_quadratic_claim_refuted=True,
            unknown_memory_coherence_not_implicitly_dephased=True)


    def test_04_fine_source_cannot_flow_into_invariant_vector_code(self):
        a=old.model();groups,actions=symmetry_orbits_and_actions()
        # Exact sparse integer H, without forming a dense 2304-square matrix.
        assert np.array_equal(a['h'],np.rint(a['h']))
        rows,cols=np.nonzero(a['h']);entries={}
        for m in range(a['modes']):
            for row,col in zip(rows,cols):
                entries[(m*a['d']+int(row),m*a['d']+int(col))]=int(a['h'][row,col])
            for dg in range(a['d']):
                v=dg//a['G']
                entries[(m*a['d']+dg,(m^(1<<v))*a['d']+dg)]=1
        for action in actions:
            transformed={(int(action[row]),int(action[col])):value
                for (row,col),value in entries.items()}
            self.assertEqual(entries,transformed)
        cases=[]
        for anchor,record in ((0,1),(0,2),(2,0),(2,3)):
            ka=2 if anchor<2 else 4;kb=2 if record<2 else 4
            orbit_size=ka*(kb-int((anchor<2)==(record<2)))
            ratio=max(Q(sum(j//a['d']==1<<record and (j%a['d'])//a['G']==anchor
                for j in group),len(group)) for group in groups)
            self.assertEqual(ratio,Q(1,orbit_size))
            cases.append(dict(anchor=anchor,record=record,pointer_orbit_size=orbit_size,
                exact_maximum_invariant_code_weight=str(ratio),
                half_trace_distance_lower_to_any_code_state=str(1-ratio)))
        OBS['fine_source_to_symmetric_code']=dict(active_dimension=a['modes']*a['d'],
            joint_role_permutation_orbits=len(groups),integer_H_nonzero_entries=len(entries),
            all_four_permutation_generator_commutators_exactly_zero=True,
            cases=cases,code_is_invariant_vector_space_not_covariant_density_space=True,
            no_general_collective_endpoint_impossibility_claimed=True)


def run():
    OBS.clear();stream=io.StringIO()
    result=unittest.TextTestRunner(stream=stream,verbosity=0).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():raise AssertionError(stream.getvalue())
    return dict(date='2026-09-28',scientific_baseline_round=514,
        numbered_round_created=False,numbered_scientific_test_increment=0,
        diagnostic_tests=result.testsRun,failures=len(result.failures),errors=len(result.errors),
        python=platform.python_version(),numpy=np.__version__,
        dependency_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in (
            'historical_membership_interface.py','distributed_role_reader.py',
            'research_note_494.md','research_note_506.md','research_note_510.md','research_note_514.md')},
        scope=dict(existing_measurement_tool_corollary=True,
            full_classical_record_and_unknown_reference_bound=True,
            same_real_H_during_monitoring_waits=True,
            initial_memory_measurement_counted=True,
            fine_source_to_symmetric_vector_code_boundary_proved=True,
            global_quantum_environment_same_bound=False,
            natural_unobserved_record_stability_proved=False,
            autonomous_instruments_or_timers_generated=False,
            shared_spatial_endpoint_generated=False,dimension_three_generated=False,
            full_GR_goal_completed=False),observations=OBS)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('scientific_baseline_round',
        'numbered_round_created','diagnostic_tests','failures','errors')},ensure_ascii=False))
