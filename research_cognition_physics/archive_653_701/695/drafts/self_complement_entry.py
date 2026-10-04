"""695 entry: unused-spin complement and phase at the original fixed background.
No Monte Carlo integration or inference of positivity from sampledE fields.
"""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path.insert(0,str(ROOT))
import joint_gauss_boundary_functional as base
import joint_rational_physical_limit as soft
import joint_critical_source_certificate as cert
TARGET=HERE/'self_complement_entry_results.json'


def run():
    sp=base.internal.spin;k4=sp.G5@sp.GAMMA[1]
    assert np.array_equal(k4@k4,-np.eye(4))
    assert np.array_equal(k4.conj().T,-k4)
    pairing=k4.T@base.internal.B@k4
    sign=1 if np.array_equal(pairing,base.internal.B) else -1
    assert np.array_equal(pairing,sign*base.internal.B)
    assert np.array_equal(k4@sp.G5+sp.G5@k4,np.zeros((4,4)))
    for gamma in (sp.GAMMA[0],sp.GAMMA[3]):
        assert np.array_equal(k4@(sp.G5@gamma)+(sp.G5@gamma)@k4,np.zeros((4,4)))
    k=np.kron(np.kron(np.eye(4),k4),np.eye(16))
    full_links,erandom,phis=soft.fixture()
    # Original group path of694, same anchor; reuse existing entry module.
    import importlib.util
    spec=importlib.util.spec_from_file_location('entry694',ROOT/'round694_drafts/critical_source_entry.py')
    previous=importlib.util.module_from_spec(spec);spec.loader.exec_module(previous)
    anchor_links=previous.links_at(4*np.arctan(float(cert.LEFT))/np.pi)
    e0=np.tile(np.eye(10)[0],(4,1));ealt=e0.copy();ealt[[1,2]]*=-1
    rows=[]
    for label,links in [('original_full_group_fixture',full_links),('certified_flux_anchor',anchor_links)]:
        u,v,d,h,gap=base.kernel(links);vk=k@u
        complement_error=float(np.max(abs(vk@vk.conj().T-v@v.conj().T)))
        assert np.max(abs(k@h@k.conj().T+h))==0
        # The anchor has gap~4e-6: comparing independent numerical eigenspaces
        # amplifies roundoff by1/gap. The theorem uses exact Clifford identities.
        complement_tolerance=20*len(h)*np.finfo(float).eps/gap
        assert complement_error<complement_tolerance
        detframe=np.linalg.det(np.column_stack((u,vk)));values=[]
        for name,e in [('common_axis',e0),('original_fixture_E',erandom),('alternating_axis',ealt)]:
            mat=base.fixed_matrices(e,phis);a=base.pf(u.T@mat[2]@u);av=base.pf(vk.T@mat[2]@vk)
            assert abs(a)>1e-24
            relative=float(abs(a-detframe*av.conjugate())/abs(a))
            fixed_phase=float(abs(av-a)/abs(a))
            assert max(relative,fixed_phase)<3e-10
            n=base.regular(u,v,d,mat,lam=0);weight=n['weight']
            values.append(dict(auxiliary_configuration=name,weight=base.old.cpair(weight),
                auxiliary_phase_identity_error=relative,paired_auxiliary_identity_error=fixed_phase))
        reference=complex(*values[0]['weight']);assert abs(reference)>1e-24
        for r in values:r['weight_relative_to_common_axis']=base.old.cpair(complex(*r['weight'])/reference)
        rows.append(dict(background=label,Wilson_gap=gap,complement_projector_error=complement_error,
            complement_diagnostic_gap_scaled_tolerance=complement_tolerance,
            all_original_channels_retained=True,auxiliary_rows=values))
    deps=('research_note_674.md','research_note_676.md','research_note_692.md','research_note_694.md',
          'joint_critical_source_certificate_results.json')
    return dict(date='2026-10-02',entry_round=695,latest_formal_round=694,not_formal_round=True,
        exact_spin_K_square=-1,exact_spin_pairing_sign=sign,
        K_maps_original_H_negative_to_positive_for_all_retained_links=True,
        auxiliary_phase_independent_of_E_modulo_sign_in_balanced_rank=True,rows=rows,
        global_auxiliary_nonnegativity_not_proved=True,full_S9_or_Haar_or_Hb_integral_not_evaluated=True,
        original_physical_RP_not_proved=True,
        dependency_hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in deps})


if __name__=='__main__':
    out=run()
    if TARGET.exists():assert out==json.loads(TARGET.read_text('utf8'))
    else:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(dict(entry_round=695,pairing_sign=out['exact_spin_pairing_sign'],
        ratios=[(r['background'],[v['weight_relative_to_common_axis'] for v in r['auxiliary_rows']]) for r in out['rows']])) )
