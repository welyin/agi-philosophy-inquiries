"""689: physical Weyl contact symbol and actual equal-time CAR projection.

The formal source identity and polynomial continuation are analytical results.
This code checks original matrices, noncommuting multi-time group insertions,
singular backgrounds, and the exact projected expectation inherited from688.
"""
import argparse
from fractions import Fraction
import hashlib
import importlib.util
import json
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_physical_car_projector_results.json'
spec=importlib.util.spec_from_file_location('source_entry689',HERE/'round689_drafts/physical_car_source_probe.py')
entry=importlib.util.module_from_spec(spec);spec.loader.exec_module(entry)


def all_time_insertions():
    rows=[]
    for nt in (2,3,4):
        single=entry.single
        links=[single.gauge.rep(*single.gauge.group(689120+20*nt+x,.27)) for x in range(nt)]
        f=entry.original_fields(links);Q=f['Q'];I=np.eye(len(Q));H=np.zeros_like(Q)
        for x in range(nt):
            h=np.kron(single.gauge.rep(*single.gauge.group(689230+20*nt+x,.43)),np.eye(2))
            H[32*x:32*x+32,32*x:32*x+32]=h
        B=(H-I)/2;D=(H+I)/2
        # The first block uses the actual original spinor projectionsK,A.
        polynomial=np.block([[2*f['K'],B],[2*f['A'],D]])
        actual=entry.detvalue(polynomial)
        finite_lattice=entry.detvalue(I-H@Q)
        closed_word=-np.linalg.matrix_power(H@Q,nt)[:32,:32]
        fock=entry.detvalue(np.eye(32)+closed_word)
        relative=float(max(abs(actual-finite_lattice),abs(actual-fock))/max(1.,abs(fock)))
        assert relative<2e-11
        # Gauge covariance at actual separate time sites, original group only.
        gauge=np.zeros_like(Q)
        for x in range(nt):
            g=np.kron(single.gauge.rep(*single.gauge.group(689390+20*nt+x,.35)),np.eye(2))
            gauge[32*x:32*x+32,32*x:32*x+32]=g
        transformed=entry.detvalue(I-(gauge@H@gauge.conj().T)@(gauge@Q@gauge.conj().T))
        covariance=float(abs(transformed-finite_lattice)/max(1.,abs(finite_lattice)))
        assert covariance<2e-11
        rows.append(dict(N=nt,original_all_time_polynomial_error=relative,
            independent_local_gauge_error=covariance,
            actual_word_determinant_real=float(fock.real),actual_word_determinant_imag=float(fock.imag)))
    return rows


def projected_expectation():
    prior=json.loads((HERE/'joint_gauss_support_marginal_results.json').read_text('utf8'))
    modes=prior['exact_projected_support']['rows'];n0=modes[0]['Gauss_multiplicity']
    rows=[]
    for saved in prior['CAR_marginal']['rows']:
        N=saved['N']
        denominator=sum(r['Gauss_multiplicity']*Fraction(r['mu'])**N for r in modes)
        # First Haar integral over the actual physical insertion removes CAR
        # character at each original background: int_h chiF(g h) = dim F^G.
        numerator=n0*sum(r['auxiliary_singlets']*Fraction(r['mu'])**N for r in modes)
        p=numerator/denominator
        assert p==Fraction(saved['CAR_singlet_probability'])
        assert 1-p==Fraction(saved['trace_distance_to_CAR_only_Gauss_state'])
        rows.append(dict(N=N,original_physical_projector_expectation=str(p),
            original_physical_expectation_float=float(p),
            gap_to_CAR_only_reference=str(1-p),
            distinct_time_repeated_projector_expectation=str(p) if N>=2 else None))
    return dict(rows=rows,full_group_Haar_used_analytically=True,
        zero_weight_backgrounds_kept_polynomial=True,
        ordinary_same_time_Grassmann_product_not_identified_as_CAR_product=True,
        original_full_HF_not_replaced_by_CAR_only_reference=True)


def run():
    probe=entry.run()
    assert probe==json.loads(entry.TARGET.read_text('utf8'))
    deps=('research_note_653.md','research_note_654.md','research_note_646.md','research_note_661.md',
        'research_note_670.md','research_note_679.md','research_note_688.md',
        'joint_gauss_support_marginal.py','joint_gauss_support_marginal_results.json',
        'round689_drafts/physical_car_source_probe.py','round689_drafts/physical_car_source_probe_results.json')
    return dict(date='2026-10-02',round=689,tests_run=2,failures=0,errors=0,
        physical_source_and_group_insertions=dict(entry=probe,multi_time=all_time_insertions()),
        original_projector=projected_expectation(),
        dependency_hashes={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in deps},
        scope=dict(inherited653_time_transfer_not_new=True,
            original_pure_time_physical_Y_projector_dictionary_proved=True,
            contact_prescription_fixed_by_original_matrices=True,
            scalar_measurement_instrument_or_new_cognitive_design=False,
            general_dynamic_Gauss_Q0_RP_decided=False,original_full_HF_identity_proved=False,
            common_spacetime_limit_or_quantum_GR_completed=False,
            old_space_and_full_goal_unchanged=True))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true');args=parser.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(round=689,all_checks_passed=True,
        max_multi_time_error=max(r['original_all_time_polynomial_error'] for r in result['physical_source_and_group_insertions']['multi_time']),
        original_projector_at_N2=result['original_projector']['rows'][1]['original_physical_expectation_float'])))
