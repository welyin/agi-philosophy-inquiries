"""729 entry: restore574's actual weak and circle links in the727 star."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent;ARCHIVE=HERE.parent
sys.path.insert(0,str(ARCHIVE))
import joint_matter_ground_source as old
TARGET=HERE/'nonflat_reference_entry_results.json'


def original_star():
    N=32;q=old.matter.original.shared_source(N)
    field=old.old.wall.old.fields(N)
    assert np.max(abs(q['phi']-field['q']['phi']))<1e-14
    W=old.matter.original.lattice.su2(q['a'],q['eps'])
    z=np.exp(1j*q['eps']*q['a0'])
    idx=(0,N//4,N//8)
    points=[idx,(1,N//4,N//8),(0,N//4+1,N//8),(0,N//4,N//8+1)]
    pairs=[old.matter.mass_matrices(q['phi'][p]) for p in points]
    h=old.block([p[0] for p in pairs]);delta=old.block([p[1] for p in pairs]);hop=np.zeros_like(h)
    directions=[];higgs=[q['f']['X'][idx]]
    for axis,alpha in enumerate(old.chiral.kinetic_matrices()):
        weak=W[idx+(axis,)];phase=z[idx+(axis,)]
        R=old.matter.representation(np.eye(3),weak,phase)
        J=-1j*alpha[:32,:32]@R/(2*q['eps'])
        sl=slice(32*(axis+1),32*(axis+2))
        hop[:32,sl]=J;hop[sl,:32]=J.conj().T
        higgs.append(phase**3*weak@q['f']['X'][points[axis+1]])
        directions.append(dict(axis=axis,weak_identity_distance=float(np.linalg.norm(weak-np.eye(2))),
                               circle_phase=[float(phase.real),float(phase.imag)]))
    return old.bdg(h+hop,delta),old.bdg(-hop,np.zeros_like(hop)),directions,np.array(higgs).T


def run():
    B,G,links,higgs=original_star();flat,Gflat=old.star()
    d=old.data(B,G);df=old.data(flat,Gflat)
    delta=float(np.linalg.norm(B-flat,2))
    ranks=np.linalg.svd(higgs,compute_uv=False)
    assert min(ranks)>.001 and delta>.013
    assert d['gap']>.001
    # Retain original fiber source; uniform gamma scales these restored edges.
    step=2e-5
    plus=B+(np.exp(-step)-1)*(-G);minus=B+(np.exp(step)-1)*(-G)
    fd=(old.data(plus,G)['E']-old.data(minus,G)['E'])/(2*step)
    source_error=abs(fd-d['mean']);assert source_error<1e-6
    deps=('research_note_574.md','research_note_604.md','research_note_727.md',
          'research_note_728.md','joint_curved_quantum_source.py',
          'joint_gauss_continuum_sampling.py','joint_matter_ground_source.py')
    return dict(entry_round=729,new_formal_round=False,tests_run=2,failures=0,errors=0,
        results=dict(original_N=32,nodes=4,physical_modes=128,links=links,
                     transported_Higgs_singular_values=ranks.tolist(),
                     restored_link_perturbation_norm=delta,
                     flat_star_gap=df['gap'],original_link_star_gap=d['gap'],
                     original_link_ground_energy=d['E'],flat_ground_energy=df['E'],
                     original_link_geometry_source_mean=d['mean'],
                     original_link_geometry_source_variance=d['noise'],
                     same_source_derivative_error=source_error,
                     entire_periodic_graph_not_diagonalized=True,
                     original_full_Gauss_momenta_or_geometry_not_declared_solved=True),
        dependencies={n:hashlib.sha256((ARCHIVE/n).read_bytes()).hexdigest() for n in deps},
        scope='Same727 conditional star with574 original weak/circle midpoint links restored. Finite spectra and source comparison only. The unit-link rational bound cannot simply absorb this perturbation. No full periodic graph gap, curved Dirac spin connection, original Einstein backreaction or new formal round claimed.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');a=p.parse_args();r=run()
    if a.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(r['results'],ensure_ascii=False,indent=2))
