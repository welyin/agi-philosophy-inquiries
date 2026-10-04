"""693 entry: screen exact finite Peter-Weyl closure using frozen676 support.

No new full averaged sign claim. Exact open zero support is inherited, not
re-proved by numerical null-value thresholds.
"""
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT))
import joint_nonflat_seed_support as support
import joint_auxiliary_orbit_quadrature as orbit
base=support.base
TARGET=HERE/'finite_character_entry_results.json'


def run():
    cert=support.load('exact_flux_seed.py').run()
    cert['mod4_multiplicities']={str(k):v for k,v in cert['mod4_multiplicities'].items()}
    assert cert==json.loads((ROOT/'round676_drafts/exact_flux_seed_results.json').read_text('utf8'))
    free=np.tile(np.eye(16,dtype=complex),(2,4,1,1))
    e=np.tile(np.eye(10)[0],(4,1))
    free_source=orbit.entry.evaluate(free,e,base.mass.car.PHI)
    normalized=free_source['source']/free_source['weight']
    expected=12*((np.sqrt(5)-1)/2)**3
    assert abs(normalized-expected)<2e-12 and free_source['weight'].real>0
    flux=free.copy()
    for i,(t,x) in enumerate(base.prior.SITES):
        for mu,z in enumerate(((-1.)**t if x==1 else 1.,1j if x==1 else 1.)):
            flux[mu,i]=base.prior.rep(np.eye(3),np.eye(2),z)
    s4=1j*base.internal.spin.GAMMA[1]@base.internal.spin.GAMMA[2]
    ev,vec=np.linalg.eigh(s4)
    frames=[np.kron(np.kron(np.eye(4),vec[:,ev*sg>0]),np.eye(16)) for sg in (-1,1)]
    rows=[]
    for scale in (0.,.002):
        links=flux.copy()
        if scale:
            for mu in range(2):
                for i in range(4):links[mu,i]=base.prior.rep(*base.prior.group(69320+4*mu+i,scale))@links[mu,i]
        u,_,_,h,gap=base.kernel(links)
        counts=[int(sum(np.linalg.eigvalsh(f.conj().T@h@f)<0)) for f in frames]
        assert counts==[72,56]
        mat=base.fixed_matrices(e,base.mass.car.PHI)
        null=int(sum(np.linalg.svd(u.T@mat[2]@u,compute_uv=False)<1e-10))
        assert null>=16
        rows.append(dict(perturbation_scale=scale,negative_spin_counts=counts,
            Wilson_gap=gap,numerical_auxiliary_nullity_diagnostic=null))
    deps=('research_note_676.md','round676_drafts/exact_flux_seed.py',
          'round676_drafts/exact_flux_seed_results.json','research_note_657.md',
          'research_note_691.md','research_note_692.md','joint_auxiliary_orbit_quadrature_results.json')
    return dict(date='2026-10-02',entry_round=693,latest_formal_round=692,new_formal_round=False,
        exact676_certificate_reproduced=True,free_original_weight=base.old.cpair(free_source['weight']),
        free_original_B_source=base.old.cpair(free_source['source']),
        free_normalized_B_norm=float(normalized.real),old_formula_error=float(abs(normalized-expected)),
        inherited_zero_support_diagnostics=rows,
        finite_exact_global_Peter_Weyl_representation_screened_out=True,
        controlled_infinite_or_approximate_representation_not_excluded=True,
        full_Haar_integration_or_RP_not_decided=True,
        dependency_hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in deps})


if __name__=='__main__':
    result=run()
    if TARGET.exists():assert result==json.loads(TARGET.read_text('utf8'))
    else:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(dict(entry_round=693,exact676_certificate_reproduced=True,
        free_B_norm=result['free_normalized_B_norm'],all_checks_passed=True)))
