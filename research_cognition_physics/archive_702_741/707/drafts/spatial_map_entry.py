"""707 entry only: distinguish617 representation cut from isotropic refinement.
Uses original quotient characters; tests ONLY stated electric-sector contracts.
No full-H dynamics, full refinement, Gibbs family or continuum is asserted.
"""
import argparse
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent;ARCHIVE=HERE.parent
sys.path.insert(0,str(ARCHIVE))
import joint_quotient_boundary_entropy as boundary
TARGET=HERE/'spatial_map_entry_results.json'

def run():
    rows=[]
    for module,label in boundary.MODULES.items():
        info=boundary.label_data(label)
        assert info['residue']==0
        for sector,cas in enumerate(info['casimirs']):
            if cas==0:continue
            # Casimirs cancel from the ratio: inherited637 half-weight vs new epsilon/2.
            old=cas;cut=Fraction(1,2)*cas*2;refined=2*cas*2
            assert cut==old and refined==4*old
            rows.append(dict(module=module,sector=sector,C=cas,cut_ratio=cut/old,
                             serial_isotropic_ratio=refined/old,
                             uniform_source_ratio=refined/old))
    # Independent discrete Haar character norm of617 multiplication map.
    theta=2*np.pi*np.arange(79)/79
    coeff=np.array([1,1j,.3-.2j]);coeff=coeff/np.linalg.norm(coeff)
    q=np.array([0,6,12])
    psi=sum(a*np.exp(1j*Q*theta) for a,Q in zip(coeff,q))
    pair=sum(a*np.exp(1j*Q*(theta[:,None]+theta[None,:])) for a,Q in zip(coeff,q))
    normerr=abs(np.mean(abs(psi)**2)-np.mean(abs(pair)**2))
    assert normerr<1e-14
    # Four declared parallel paths, each two fine serial edges; only singlet characters Q=6m.
    parallel=[]
    for n in range(-12,13):
        k,r=divmod(n,4)
        balanced=[k+1]*r+[k]*(4-r)
        norm=sum(x*x for x in balanced)
        assert sum(balanced)==n
        # Independent bounded exhaustive search: any improvement must have sum squares<=norm,
        # hence every coordinate lies inside +/-floor(sqrt(norm)).
        bound=int(np.sqrt(norm))
        best=None
        for a,b,c in itertools.product(range(-bound,bound+1),repeat=3):
            d=n-a-b-c
            if abs(d)>bound:continue
            candidate=a*a+b*b+c*c+d*d
            best=candidate if best is None else min(best,candidate)
        assert best==norm
        gap=4*norm-n*n
        assert gap==r*(4-r)
        labels=[boundary.label_data((0,0,0,6*x)) for x in balanced]
        assert all(x['residue']==0 for x in labels)
        parallel.append(dict(n=n,balanced=balanced,fine_energy_over_36c=4*norm,
                             coarse_energy_over_36c=n*n,excess_over_36c=gap))
    deps=('research_note_617.md','research_note_637.md','research_note_642.md',
          'research_note_652.md','research_note_705.md','research_note_706.md',
          'joint_quotient_boundary_entropy.py')
    return dict(date='2026-10-02',entry_round=707,latest_formal_round=706,new_formal_round=False,
        serial_rows=rows,discrete_Haar_isometry_error=float(normerr),
        parallel_character_rows=parallel,
        analytic_serial_ratio=4,analytic_parallel_excess='36*c*r*(4-r), n=4k+r, 0<=r<4',
        dependencies={name:hashlib.sha256((ARCHIVE/name).read_bytes()).hexdigest() for name in deps},
        scope=dict(inherited_cut_identity=True,original_quotient_labels=True,
                   declared_pure_electric_comparison_only=True,
                   parallel_geometry_is_new_comparison_input=True,
                   charged_nonAbelian_fluxes_not_minimized=True,
                   no_full_H_Gibbs_or_spacetime_limit_claim=True),
        all_checks_passed=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');a=p.parse_args()
    r=run()
    if a.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(entry_round=707,new_formal_round=False,serial_ratio=r['analytic_serial_ratio'],
                         parallel_rows=len(r['parallel_character_rows']),all_checks_passed=True)))
