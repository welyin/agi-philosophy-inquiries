"""Unnumbered760 tool-applicability diagnostic on the original periodic lepton sector."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path.insert(0,str(ROOT))
import joint_curved_periodic_reference as old
import joint_reference_constraint_strata as backgrounds
TARGET=HERE/'transport_source_entry_results.json'

def run():
    N=4
    h,d,K,_,psi_old=old.matrices(N,old.LEPTON)
    q,psi_new,*_=backgrounds.completed(N,1.)
    scalar=old.background(N)[0]
    assert np.max(abs(q['phi']-scalar['phi']))<1e-13
    # Replace only the old fixed geometry weights by the actual753 geometry.
    # The complete lepton sector is color neutral, so753's new color links act
    # trivially here. Quarks are not discarded from the full analytic model.
    w=np.repeat((psi_old/psi_new).ravel(),8)
    hop=w[:,None]*K*w[None,:];mass=h-K
    B=old.old.bdg(mass+hop,d);M=old.old.bdg(mass,d)
    C=B@M-M@B
    frob2=float(np.vdot(C,C).real)
    row_bound=float(np.max(np.sum(abs(B),axis=1)))
    bound=frob2/(4*len(B)*row_bound**2)
    assert np.max(abs(B-B.conj().T))<1e-13 and np.max(abs(M-M.conj().T))<1e-13
    assert frob2>1 and bound>1e-5
    deps=('research_note_727.md','research_note_729.md','research_note_753.md',
          'research_note_759.md','joint_curved_periodic_reference.py',
          'joint_reference_constraint_strata.py')
    return dict(kind='unnumbered_transport_applicability_diagnostic',latest_scientific_round=759,
        N=N,physical_lepton_modes=8*N**3,Nambu_dimension=len(B),
        geometry_min=float(psi_new.min()),geometry_max=float(psi_new.max()),
        source_mass_commutator_Frobenius_squared=frob2,
        Hermitian_operator_norm_row_sum_upper=row_bound,
        normalized_Hilbert_Schmidt_off_energy_block_lower=bound,
        lower_bound_scope='Analytic pinching inequality applied to the represented Nambu matrices; floating coefficients and products are diagnostics, not interval certification. This number is not the original physical Fock trace variance.',
        full_original_periodic_lepton_sector=True,old_original_mass_weak_circle_and_Weyl_links_retained=True,
        new_color_trivial_on_leptons=True,whole_quark_model_not_simulated=True,
        phase_space_transport_and_physical_time_limit_proven=False,
        no_new_numbered_tests=True,
        dependency_hashes={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in deps},
        all_entry_checks_passed=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args();r=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps({k:r[k] for k in ('latest_scientific_round','source_mass_commutator_Frobenius_squared','normalized_Hilbert_Schmidt_off_energy_block_lower','all_entry_checks_passed')}))

