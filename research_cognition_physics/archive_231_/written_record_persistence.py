"""Round 514: a generated remote Z record with a continuous retention window.

All branches and source cost remain explicit. No switch removes the real g term.
"""
import argparse
from fractions import Fraction as Q
from functools import lru_cache
import hashlib
import io
import json
import math
from pathlib import Path
import platform
import unittest

import numpy as np
import historical_membership_interface as old
import coherent_graph_mean_obstruction as graph
from distributed_role_reader import exponential, opnorm

HERE=Path(__file__).resolve().parent
TARGET=HERE/'written_record_persistence_results.json'
OBS={}
FORMATION=Q(1,2048)
RETENTION=Q(1,8)


def tail(x,start):
    return x**start/Q(math.factorial(start))/(1-x/Q(start+1))


def uniform_source_remainder(t,g=Q(1)):
    return 24*abs(g)*t**4+tail((6+abs(g))*t,4)


def initial(anchor,dtype=complex):
    a=old.model()
    out=np.zeros((a['modes'],a['d'],a['G']),dtype=dtype)
    out[0,anchor*a['G']:(anchor+1)*a['G']]=np.eye(a['G'],dtype=dtype)
    return out


@lru_cache(None)
def source(anchor):
    out=initial(anchor)
    term=out.copy()
    for k in range(1,19):
        term=(-1j*float(FORMATION)/k)*old.apply_historical(term)
        out+=term
    return out


def selected(anchor,state):
    a=old.model()
    return {b:state[1<<b,anchor*a['G']:(anchor+1)*a['G']]
            for b in range(a['n']) if b!=anchor}


class Audit(unittest.TestCase):
    def test_01_exact_third_order_returned_remote_record(self):
        a=old.model();rows=[]
        for anchor in (0,2):
            power=initial(anchor,np.int64)
            selected_powers=[]
            for k in range(1,4):
                power=old.apply_historical(power)
                selected_powers.append(selected(anchor,power))
            self.assertTrue(all(not np.any(v) for row in selected_powers[:2] for v in row.values()))
            gram=np.zeros((a['G'],a['G']),dtype=np.int64)
            for b,matrix in selected_powers[2].items():
                incidence=np.diag([tuple(sorted((anchor,b))) in tree for tree in a['trees']]).astype(np.int64)
                self.assertTrue(np.array_equal(matrix,incidence))
                gram+=matrix.T@matrix
            degree=3 if anchor<2 else 1
            self.assertTrue(np.array_equal(gram,degree*np.eye(a['G'],dtype=np.int64)))
            rows.append(dict(anchor=anchor,first_nonzero_order=3,integer_gram_factor=degree))
        OBS['exact_formation']=rows

    def test_02_size_uniform_rational_source_certificate(self):
        t=FORMATION
        remainder=uniform_source_remainder(t)
        self.assertLess(remainder,t**3/12)
        self.assertLess(remainder,128*t**4)
        local_tail=tail(11*t,19)
        self.assertLess(local_tail,t**3/Q(10**30))
        OBS['source_certificate']=dict(duration=str(t),coupling=1,
            F_interaction_local_norm=7,old_F_adjacency_commutator_bound=48,
            third_order_time_dependence_remainder=str(24*t**4),
            higher_order_tail=str(tail(7*t,4)),total_remainder=str(remainder),
            amplitude_lower=str(t**3/12),probability_lower=str(t**6/144),
            finite_six_graph_Taylor_tail=str(local_tail),
            source_bound_uniform_in_network_size=True,
            no_F_global_norm_in_uniform_proof=True)

    def test_03_full_source_instrument_unknown_graph_columns(self):
        a=old.model();rows=[]
        for anchor in (0,2):
            state=source(anchor)
            flat=state.reshape(-1,a['G'])
            full=flat.conj().T@flat
            self.assertLess(opnorm(full-np.eye(a['G'])),1e-13)
            scale=float(FORMATION**3)
            matrices=selected(anchor,state)
            scaled_effect=sum((k/scale).conj().T@(k/scale) for k in matrices.values())
            spectrum=np.linalg.eigvalsh(scaled_effect)
            self.assertGreater(spectrum[0],1/144)
            # Independent X-sector implementation, also including every failure.
            spectral=old.history(FORMATION)[:,:,anchor*a['G']:(anchor+1)*a['G']]
            difference=opnorm((state-spectral).reshape(-1,a['G']))
            self.assertLess(difference,2e-13)
            self.assertGreater(np.linalg.eigvalsh(np.eye(a['G'])-scale**2*scaled_effect)[0],0)
            rows.append(dict(anchor=anchor,complete_instrument_error=graph.short(opnorm(full-np.eye(a['G']))),
                independent_X_sector_difference=graph.short(difference),
                success_effect_scaled_by_T6_min=graph.short(spectrum[0]),
                success_effect_scaled_by_T6_max=graph.short(spectrum[-1]),
                source_probability_min=graph.short(scale**2*spectrum[0]),
                source_probability_max=graph.short(scale**2*spectrum[-1]),
                labels_summed_not_each_label_guaranteed=True))
        OBS['full_source_instrument']=rows

    def test_04_retention_full_operator_bound(self):
        a=old.model();u=RETENTION;bound=Q(3,2)*u*u
        rows=[]
        for record in (0,2):
            outside=np.repeat(np.arange(a['n'])!=record,a['G'])
            largest=0.
            for sign in a['signs']:
                real=a['h']+np.diag(np.repeat(sign,a['G']))
                omitted=sign.copy();omitted[record]=0
                comparison=a['h']+np.diag(np.repeat(omitted,a['G']))
                error=opnorm((exponential(float(u)*real)-exponential(float(u)*comparison))[:,outside])
                largest=max(largest,error)
            self.assertLess(largest,float(bound))
            rows.append(dict(record_site=record,all_memory_X_blocks=64,
                observed_full_operator_error=graph.short(largest),operator_upper=str(bound)))
        OBS['retention']=dict(duration=str(u),operator_bound=str(bound),
            erasure_upper=str(bound**2),ordinary_joint_channel_error_upper=str(2*bound),
            cases=rows,initial_packet_outside_record_required=True,
            comparison_missing_term_is_not_physical_switch=True)

    def test_05_generated_source_then_same_H_retention(self):
        a=old.model();future=old.history(RETENTION)
        beta=Q(3,2)*RETENTION**2;rows=[]
        for anchor in (0,2):
            succ=np.zeros((a['G'],a['G']),complex)
            bad=np.zeros_like(succ)
            continue_source=future[:,:,anchor*a['G']:(anchor+1)*a['G']]
            for b,k in selected(anchor,source(anchor)).items():
                k=k/float(FORMATION**3)
                succ+=k.conj().T@k
                continued=np.einsum('zpg,gh->zph',continue_source,k)
                # Since X_b commutes with H, starting e_b means XOR b on outputs.
                erased=np.array([bool(z & (1<<b)) for z in range(a['modes'])])
                flat=continued[erased].reshape(-1,a['G'])
                bad+=flat.conj().T@flat
            eig,vec=np.linalg.eigh(succ)
            invroot=(vec*(1/np.sqrt(eig)))@vec.conj().T
            ratios=np.linalg.eigvalsh(invroot@bad@invroot)
            margin=float(np.linalg.eigvalsh(float(beta**2)*succ-bad)[0])
            self.assertGreater(margin,0)
            self.assertLess(ratios[-1],float(beta**2))
            rows.append(dict(anchor=anchor,conditional_erasure_min=graph.short(ratios[0]),
                conditional_erasure_max=graph.short(ratios[-1]),
                unknown_graph_operator_inequality_margin=graph.short(margin),
                source_H_equals_retention_H=True,source_measurement_backaction_retained=True))
        OBS['composed_source_and_retention']=rows

    def test_06_occupied_counterexample_and_conditioning_cost(self):
        a=old.model();u=Q(1,128);record=2
        future=old.history(u)[:,:,record*a['G']:(record+1)*a['G']]
        erased=np.array([bool(z&(1<<record)) for z in range(a['modes'])])
        f=future[erased].reshape(-1,a['G'])
        erasure_min=float(np.linalg.eigvalsh(f.conj().T@f)[0])
        false_bound=(Q(3,2)*u*u)**2
        self.assertGreater(erasure_min,float(false_bound)*100)
        pmin=FORMATION**6/144
        eta=pmin/1000
        beta2=(Q(3,2)*RETENTION**2)**2
        conditional_upper=beta2+(1+beta2)*(eta/2)/(pmin-eta/2)
        self.assertGreater(conditional_upper,beta2)
        self.assertLess(conditional_upper,Q(2,1000))
        OBS['scope_and_precision']=dict(occupied_counterexample_time=str(u),
            occupied_record_erasure_min=graph.short(erasure_min),
            outside_support_bound_wrongly_applied=str(false_bound),
            ideal_source_probability_lower=str(pmin),
            illustrative_full_instrument_diamond_error=str(eta),
            conditional_error_upper_with_that_input=str(conditional_upper),
            finite_precision_implementation_not_constructed=True,
            heralding_is_not_free_retry=True)


def run():
    OBS.clear();stream=io.StringIO()
    result=unittest.TextTestRunner(stream=stream,verbosity=0).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise AssertionError(stream.getvalue())
    return dict(date='2026-09-28',round=514,scientific_baseline_round=513,
        tests_run=result.testsRun,failures=len(result.failures),errors=len(result.errors),
        python=platform.python_version(),numpy=np.__version__,
        dependency_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in (
            'historical_membership_interface.py','coherent_graph_mean_obstruction.py',
            'research_note_496.md','research_note_510.md','research_note_513.md')},
        scope=dict(same_H_for_formation_and_retention=True,
            true_remote_Z_record_source=True,size_uniform_source_and_retention_bounds=True,
            arbitrary_unknown_graph_reference_covered=True,
            all_failure_branches_and_conditioning_error_accounted=True,
            retention_clock_starts_at_source_instrument_end=True,
            unchanged_real_memory_coupling=True,
            deterministic_or_high_probability_source=False,
            source_success_for_each_specified_remote_site=False,
            permanent_record_or_entire_member_table_protected=False,
            unmeasured_natural_history_proved=False,
            autonomous_readout_controller_generated=False,
            source_record_collection_deadline_proved=False,
            dimension_three_generated=False,full_GR_goal_completed=False,
            phase_closure_triggered=False),observations=OBS)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf-8',newline='\n') as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('round','tests_run','failures','errors')},ensure_ascii=False))
