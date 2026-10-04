"""694 entry: actual original-source behavior across the676 flux path.

Numerical diagnostic, not yet an exact root/nonzero-jump certificate or Haar/RP test.
"""
import hashlib
import json
import sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT))
import joint_gauss_boundary_functional as base
import joint_auxiliary_orbit_quadrature as orbit
import joint_rational_physical_limit as soft
TARGET=HERE/'critical_source_entry_results.json'
CHARGES=(1,-4,2,-3,6,0)
MULT=(6,3,3,2,1,1)


def charge_h(s,q):
    sp=base.internal.spin;dw=np.zeros((16,16),complex)
    for mu in (0,1):
        shift=np.zeros_like(dw)
        for i,(t,x) in enumerate(base.prior.SITES):
            target=((1-t,x) if mu==1 else (t,1-x));j=base.prior.SITES.index(target)
            angle=(np.pi*s/2 if mu else np.pi*s*t) if x else 0.
            phase=np.exp(1j*q*angle)*(-1 if mu==1 and t else 1)
            shift[4*i:4*i+4,4*j:4*j+4]=phase*np.eye(4)
        gamma=np.kron(np.eye(4),sp.GAMMA[3 if mu else 0])
        dw+=np.eye(16)-(shift+shift.conj().T)/2+gamma@(shift-shift.conj().T)/2
    return np.kron(np.eye(4),sp.G5)@(dw-np.eye(16))


def sector(s,q,sign=-1):
    spin=np.array([[1,0],[-sign,0],[0,1],[0,-sign]])/np.sqrt(2)
    frame=np.kron(np.eye(4),spin)
    return frame.conj().T@charge_h(s,q)@frame


def spectrum_counts(s):
    return [sum(m*int(sum(np.linalg.eigvalsh(sector(s,q,sign))<0))
                for q,m in zip(CHARGES,MULT)) for sign in (-1,1)]


def first_crossings():
    rows=[]
    for q in CHARGES[:-1]:
        last_s=0.;last=int(sum(np.linalg.eigvalsh(sector(0,q))<0))
        for s in np.linspace(0,1,257)[1:]:
            count=int(sum(np.linalg.eigvalsh(sector(s,q))<0))
            if count!=last:
                lo,hi=last_s,float(s);left=last
                for _ in range(44):
                    mid=(lo+hi)/2
                    if int(sum(np.linalg.eigvalsh(sector(mid,q))<0))==left:lo=mid
                    else:hi=mid
                rows.append(dict(charge=q,root_midpoint=(lo+hi)/2,
                    numerical_bracket=[lo,hi],negative_counts=[left,count]))
                break
            last_s=float(s);last=count
    assert rows
    return sorted(rows,key=lambda r:r['root_midpoint'])


def links_at(s):
    links=np.tile(np.eye(16,dtype=complex),(2,4,1,1))
    for i,(t,x) in enumerate(base.prior.SITES):
        for mu in (0,1):
            angle=(np.pi*s/2 if mu else np.pi*s*t) if x else 0.
            links[mu,i]=base.prior.rep(np.eye(3),np.eye(2),np.exp(1j*angle))
    return links


def actual_sources(s,e,phis):
    u,v,d,h,gap=base.kernel(links_at(s));mat=orbit.entry.matrices(e,phis)
    full=base.regular(u,v,d,mat,lam=0)
    _,ry,_=orbit.entry.reflection.reflection_matrices(mat,np.eye(4)[[2,3,0,1]])
    source=sum(c*orbit.entry.integral_coefficient(full,z) for c,z,_ in orbit.entry.source_rows(ry))
    # Independent charge decomposition of the same original256-dimensional H.
    predicted=np.sort(np.concatenate([np.tile(np.linalg.eigvalsh(charge_h(s,q)),m)
                                       for q,m in zip(CHARGES,MULT)]))
    spectral_error=float(np.max(abs(np.linalg.eigvalsh(h)-predicted)))
    assert spectral_error<2e-12
    return dict(path_parameter=s,Wilson_gap=gap,spin_negative_counts=spectrum_counts(s),
        scalar_weight=base.old.cpair(full['weight']),B_source=base.old.cpair(source),
        charge_to_full_spectral_error=spectral_error)


def run():
    crossings=first_crossings();star=crossings[0]['root_midpoint']
    assert spectrum_counts(0)==[64,64] and spectrum_counts(1)==[72,56]
    _,random_e,phis=soft.fixture()
    # Reuse originalE sphere configurations; no E interpreted as physical record.
    choices=[('original_fixture',random_e),('common_axis',np.tile(np.eye(10)[0],(4,1)))]
    rows=[]
    for label,e in choices:
        for offset in (-1e-3,-1e-5,-1e-7,1e-7,1e-5,1e-3):
            row=actual_sources(star+offset,e,phis)
            row.update(auxiliary_choice=label,offset_from_numerical_root=offset)
            rows.append(row)
    deps=('research_note_673.md','research_note_676.md','research_note_679.md','research_note_693.md',
          'joint_gauss_boundary_functional.py','round692_drafts/odd_gauss_source_probe.py')
    return dict(date='2026-10-02',entry_round=694,latest_formal_round=693,not_formal_round=True,
        path='spatial z(t,x=1)=exp(i*pi*s*t); temporal z(t,x=1)=exp(i*pi*s/2); others identity',
        first_charge_sector_crossings=crossings,rows=rows,
        endpoints_reuse676=True,all_original_16_channels_retained=True,
        exact_root_and_nonzero_limit_certificate=False,
        fixed_background_diagnostic_not_complete_Haar_or_RP=True,
        source_uses_no_inverse_or_weight_division=True,
        dependency_hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in deps})


if __name__=='__main__':
    out=run()
    if TARGET.exists():assert out==json.loads(TARGET.read_text('utf8'))
    else:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(dict(first_crossings=out['first_charge_sector_crossings'],
        rows=[{k:r[k] for k in ('auxiliary_choice','offset_from_numerical_root',
              'spin_negative_counts','Wilson_gap','scalar_weight','B_source')} for r in out['rows']])))
