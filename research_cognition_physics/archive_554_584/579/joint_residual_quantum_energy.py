"""579: the unchanged full positive potential survives the xi=1/5 spectral shift.
Finite-volume lower bounds apply to the whole physical state space by analytic forms.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_residual_quantum_energy_results.json'
spec=importlib.util.spec_from_file_location('residual_entry',HERE/'round579_drafts/residual_energy_entry_probe.py')
entry=importlib.util.module_from_spec(spec);spec.loader.exec_module(entry)


def finite_node_certificate_check():
    radius=entry.threshold_radius()
    wstar=entry.HBAR/(np.sqrt(8*entry.A)*radius*np.exp(entry.BETA*radius/2))
    volume=2.3
    minimum_nodes=int(np.ceil(2*volume/wstar))
    rows=[]
    for count in (minimum_nodes,2*minimum_nodes,4*minimum_nodes):
        cap=2*volume/count
        rc=entry.matching_radius(cap)
        logupper=radius+2/entry.BETA*np.log(wstar/cap)
        sharp=count/2*entry.lower_bound(cap)
        explicit=entry.HBAR**2*count**2/(32*volume*logupper**2)
        assert entry.RHO0<rc<=logupper*(1+1e-14)
        assert explicit<=sharp*(1+1e-14)
        # Three large cells are outside the small-w regime; retain only positivity there.
        w=np.full(count,volume/(4*(count-3)))
        w[:3]=volume/4
        small=w<=cap
        assert np.count_nonzero(small)>=count/2
        eligible=w<=wstar
        assert np.count_nonzero(~eligible)==3
        sum_bound=sum(entry.lower_bound(float(v)) for v in w[eligible])
        assert sum_bound>=sharp
        residual_budget=explicit-1.
        rows.append(dict(nodes=count,volume=float(w.sum()),root=rc,log_root_bound=float(logupper),
            implicit_total_bound=sharp,explicit_total_bound=float(explicit),
            very_unequal_sum_bound=sum_bound,large_cells_ignored=int(np.count_nonzero(~eligible)),
            necessary_extra_scalar_subtraction_for_energy_at_most_one=float(residual_budget)))
    return dict(radius_threshold=radius,cell_volume_threshold=float(wstar),minimum_node_count=minimum_nodes,
                rows=rows,large_cell_energy_replaced_only_by_proved_nonnegativity=True)


def run():
    original=json.loads((HERE/'round579_drafts/residual_energy_entry_results.json').read_text('utf8'))
    data=entry.run();assert data==original and data['checks_passed']==4
    evidence=dict(hardy_form_check=data['hardy_identity'],
        full_original_potential_check=dict(constants=data['constants'],
            sampled_min_ratio=data['potential_bound_min_sample_ratio'],
            global_bound_is_analytic_not_sampling=True),
        node_bound_check=data['node_rows'],
        nonuniform_volume_check=data['graph_rows'],
        finite_node_certificate_check=finite_node_certificate_check())
    deps=('joint_curved_quantum_source.py','research_note_574.md','research_note_578.md',
          'research_round_578_checks.json','round579_drafts/residual_energy_entry_probe.py',
          'round579_drafts/residual_energy_entry_results.json')
    return dict(round=579,tests_run=len(evidence),failures=0,errors=0,checks=list(evidence),
        evidence=evidence,dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope='original H5 full positive potential after only xi=1/5 target-curvature spectral shift; analytic all-state residual energy lower bound of order N^2/log(N)^2 at fixed hbar and bounded volume; conditional obstruction to unchanged direct regular Einstein matching, not actual ground-state asymptotics or a no-go for renormalized theory')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        payload=json.dumps(result,ensure_ascii=False,indent=2)+'\n'
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(payload)
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps({k:result[k] for k in ('round','tests_run','failures','errors')}))
