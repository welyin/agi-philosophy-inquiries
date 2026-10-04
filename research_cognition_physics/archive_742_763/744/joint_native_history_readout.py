"""744: exact single-read symmetry and native parity-history time jets.

The complete-state symmetry is analytic. Parity-history numerics are for the
original full onsite H5/CAR Hamiltonian on a declared one-vertex graph, with a
new normal even Gaussian ready wavefunction. No multi-vertex or continuum
claim, exact time propagator, interval sign certificate or ideal n instrument.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'round744_drafts'))
import native_instrument_selection_entry as entry
import history_jet_probe as fourth
import history_sixth_jet_probe as sixth
TARGET=HERE/'joint_native_history_readout_results.json'

def history_check():
    first=sixth.more_jets(0);second=sixth.more_jets((1<<30)|(1<<31))
    saved4=json.loads((HERE/'round744_drafts/history_jet_probe_results.json').read_text('utf8'))
    saved6=json.loads((HERE/'round744_drafts/history_sixth_jet_probe_results.json').read_text('utf8'))
    rows=[]
    for n in (48,64,80):
        a=sixth.measure(first,n);b=sixth.measure(second,n)
        raw=next(r for r in saved6['rows'] if r['order']==n)
        assert a==raw['empty'] and b==raw['pair']
        lowa=fourth.measure(first[:5],n);lowb=fourth.measure(second[:5],n)
        if n in (48,64):
            prev=next(r for r in saved4['rows'] if r['order']==n)
            assert lowa==prev['empty'] and lowb==prev['pair']
        delta=b['sixth']-a['sixth']
        assert -.56<delta<-.54
        assert max(abs(a['identity']),abs(b['identity']))<1e-8
        assert max(abs(lowb['derivatives'][str(j)][0]-lowa['derivatives'][str(j)][0]) for j in range(5))<1e-9
        rows.append(dict(order=n,difference_sixth=delta,parity_expectation_t6_coefficient=delta/720,
            even_parity_probability_t6_coefficient=delta/1440,
            lower_derivative_difference=[lowb['derivatives'][str(j)][0]-lowa['derivatives'][str(j)][0] for j in range(5)],
            identity_sixth_residual=max(abs(a['identity']),abs(b['identity']))))
    changes=[abs(rows[i+1]['difference_sixth']-rows[i]['difference_sixth']) for i in range(2)]
    assert max(changes)<2e-9
    return dict(graph_vertices=1,target_coordinates=5,physical_CAR_modes=32,
        quantum_boson_H_retained=True,original_potential_and_all_mass_couplings_retained=True,
        alpha=fourth.ALPHA,volume_weight=fourth.W,normal_ready_state='exp(-|x|^2/2), normalized in actual H5 measure',
        rows=rows,quadrature_changes=changes,
        no_claim_of_analytic_vanishing_of_all_lower_derivatives=True,
        no_interval_sign_certificate=True,no_multi_vertex_continuum_or_full_propagation_claim=True,
        actual_two_read_effect='Sum(r*s E_r E_s)=sin(s_ball)^2/4; original square-root instruments, zero separation at endpoint.',
        preparation_is_not_frozen_round743_compact_packet=True)

def run():
    a=entry.run();assert a==json.loads(entry.TARGET.read_text('utf8'))
    b=history_check()
    deps=('research_note_577.md','research_note_598.md','research_note_623.md','research_note_624.md',
          'research_note_717.md','research_note_718.md','research_note_743.md',
          'round744_drafts/native_instrument_selection_entry.py',
          'round744_drafts/native_instrument_selection_entry_results.json',
          'round744_drafts/history_jet_probe.py','round744_drafts/history_sixth_jet_probe.py')
    return dict(round=744,tests_run=2,failures=0,errors=0,
        checks=['original_unknown_input_instrument_selection','native_full_onsite_history_sixth_jet'],
        results=dict(selection=a,parity_history=b),
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope='Exact full-H joint symmetry constrains the native single-read unknown-code POVM. True two-read histories permit population sensitivity. Full native onsite quantum Hamiltonian jets give reproducible sixth-order numerical evidence for that signal on a declared normal Gaussian preparation. It is not a certified nonzero connected-graph instrument, exact n measurement, autonomous detector, continuum SM or gravity construction.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');a=p.parse_args();r=run()
    if a.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps({'round':744,'tests_run':2,'history':r['results']['parity_history']['rows'][-1]},ensure_ascii=False))
