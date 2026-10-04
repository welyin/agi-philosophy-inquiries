"""726 entry: original point, canonical normalization and unexchanged limits."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent;ARCHIVE=HERE.parent
sys.path.insert(0,str(ARCHIVE))
import joint_material_boundary_transport as wall
import joint_quantum_reference_forms as old
TARGET=HERE/'common_wall_limit_entry_results.json'


def run():
    N=32;data=wall.old.fields(N);c,lam,_=wall.base_data();idx=(0,N//4,N//8)
    phi=data['q']['phi'][idx];eps=2*np.pi/N;psi=float(data['psi'][idx])
    volume=eps**3*psi**6;momentum=eps**3*data['q']['p'][idx]
    h=float(np.linalg.norm(phi[:4]));df=np.r_[-lam*phi[:4]/h,1.]
    inverse=old.original.inverse(phi);Gamma=float(df@inverse@df)
    velocity=float(df@inverse@momentum/volume)
    spatial=float(psi**-4*np.sum((data['ds'][idx]-lam*data['dh'][idx])**2))
    expected=c['b']**2*psi**-4
    assert abs(velocity)<1e-14 and abs(spatial-expected)<1e-14
    # Actual652 prescribed forward differences, without identifying them with the continuum.
    field=data['q']['phi'][...,4]-lam*np.linalg.norm(data['q']['phi'][...,:4],axis=-1)
    differences=np.array([((np.roll(field,-1,axis=axis)-field)/eps)[idx] for axis in range(3)])
    Wdiscrete=float(psi**-4*(differences@differences))
    critical=volume*np.sqrt(Wdiscrete)/Gamma
    phase_rows=[]
    for ratio in (-1.5,-.5,.5,1.5):
        kick=ratio*critical
        direct=Wdiscrete-(df@inverse@(momentum+kick*df)/volume)**2
        predicted=Wdiscrete-(velocity+kick*Gamma/volume)**2
        assert abs(direct-predicted)<1e-12
        assert (direct>0)==(abs(ratio)<1)
        phase_rows.append(dict(kick_over_symbolic_boundary=ratio,normal_symbol=float(direct),
                               independent_identity_error=float(abs(direct-predicted))))
    # The1/(4 w²) coefficient is actual725. This is symbol-level scale audit,
    # not the finite-hbar normal measurement's spectral probability.
    coefficient=Gamma**2/(4*volume**2)
    dependencies=('research_note_574.md','research_note_575.md','research_note_579.md',
                  'research_note_651.md','research_note_725.md','joint_material_boundary_transport.py',
                  'joint_gravity_material_coordinates.py')
    return dict(entry_round=726,new_formal_round=False,tests_run=2,failures=0,errors=0,
                result=dict(N=N,epsilon=float(eps),psi=psi,node_volume=float(volume),
                    original_lambda=float(lam),original_wall_normal_velocity=velocity,
                    continuum_normal_squared=spatial,prescribed_graph_normal_squared=Wdiscrete,
                    finite_graph_continuum_difference=float(Wdiscrete-spatial),
                    original_Gamma=Gamma,critical_canonical_kick=float(critical),
                    Q_backaction_coefficient_per_hbar_over_sigma_squared=float(coefficient),
                    phase_rows=phase_rows,
                    measurement_decomposition_label_is_not_pointer_record=True,
                    no_finite_hbar_signature_probability_claim=True),
                dependencies={n:hashlib.sha256((ARCHIVE/n).read_bytes()).hexdigest() for n in dependencies},
                scope='Deduplication and scale audit only. Reuse574 exact Gauss packets and651 original source. Gaussian position backaction uses canonical node volume; fixed-graph hbar/sigma convergence does not give joint refinement. Classical phase symbols do not prove the actual Q spectral or conditional-record concentration.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');a=p.parse_args();r=run()
    if a.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
