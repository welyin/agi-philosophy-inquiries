"""592: original read histories, finite-band obstruction, and source control.

Full-model estimates are analytic. The added numerical history diagnostic
uses actual Gauss packets and original Kraus products with zero waiting time,
not a replacement Hamiltonian or a claimed full spectral simulation.
"""
import argparse
import hashlib
import importlib.util
import itertools
import json
from pathlib import Path
import numpy as np
import joint_geometry_work_noise as prior

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_record_band_control_results.json'
spec=importlib.util.spec_from_file_location('band_entry592',HERE/'round592_drafts/record_band_entry.py')
entry=importlib.util.module_from_spec(spec);spec.loader.exec_module(entry)


def trace_norm_from_gram(gram):
    values,vectors=np.linalg.eigh(gram)
    assert values.min()>-2e-13
    root=(vectors*np.sqrt(np.maximum(values,0)))@vectors.conj().T
    signed=root@np.diag([-1.,.5,.5])@root
    return float(np.sum(abs(np.linalg.eigvalsh(signed))))


def actual_history_check():
    h,s,p,k=prior.packet_data(128)
    n=16;alpha=.4;weight=alpha/n**2
    phases=np.exp(1j*s[...,None]*np.array([0,n,-n]))
    L,_,_=entry.instrument(s)
    def gram(density):
        return np.einsum('xy,xyi,xyj->ij',density,phases.conj(),phases)
    initial=.5*weight*trace_norm_from_gram(gram(p))
    rows=[]
    for count in (1,3,5):
        distance=0.;mass=0.;largest_probability_gap=0.
        for history in itertools.product((0,1),repeat=count):
            amplitude=np.ones_like(s)
            for r in history:amplitude*=L[...,r]
            G=gram(p*amplitude**2)
            distance+=.5*weight*trace_norm_from_gram(G)
            probability=float(G[0,0].real)
            mixture_probability=(1-weight)*probability+.5*weight*float(G[1,1].real+G[2,2].real)
            largest_probability_gap=max(largest_probability_gap,abs(mixture_probability-probability))
            mass+=probability
        assert abs(mass-1)<2e-14 and largest_probability_gap<2e-15
        assert distance<=initial+2e-14
        rows.append(dict(reads=count,total_history_probability=mass,
                         full_classical_quantum_trace_distance=distance,
                         maximum_same_history_probability_error=largest_probability_gap))
    return dict(original_phase_mixture=True,initial_trace_distance=initial,rows=rows,
                original_packet_and_primitive_instruments=True,waiting_times_are_zero=True,
                same_reports_do_not_imply_same_conditional_quantum_states=True,
                arbitrary_waiting_times_are_covered_analytically_not_simulated=True,
                actual_spectral_projection_not_computed=True)


def run():
    old=entry.run();assert old==json.loads(entry.TARGET.read_text('utf8'))
    evidence=dict(old['evidence']);evidence['actual_history_check']=actual_history_check()
    names=('research_note_555.md','research_note_556.md','research_note_590.md','research_note_591.md',
           'joint_geometry_work_noise.py','joint_full_spatial_metric.py',
           'round592_drafts/record_band_entry.py','round592_drafts/record_band_entry_results.json',
           'round592_drafts/record_band_entry.md','round592_drafts/record_band_entry_checks.json')
    return dict(round=592,tests_run=4,failures=0,errors=0,evidence=evidence,
                dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names},
                scope='Original fixed finite graph, fixed positive geometry, and primitive reads. No nonzero finite-dimensional exact instrument-invariant sector. A propagated H^2 input budget controls finite-history CQ states and form-source means of an explicitly defined approximation. Spectral reset is a mathematical comparison, not a free physical operation. No moving geometry, original spectrum numerics, full source noise, Einstein dynamics, or fixed-hbar continuum.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args();r=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==r
    print(json.dumps(r,ensure_ascii=False,indent=2))
