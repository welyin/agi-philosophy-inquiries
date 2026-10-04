"""673 entry: boundary holonomy changes the original signed sewing weight.

Reuses670 covariance and643 Gauss endpoint rule. Not Haar integration or a
new formal round. In particular a conditioned complex weight is not a no-go.
"""
import hashlib
import json
import sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;BASE=HERE.parent;sys.path.insert(0,str(BASE))
import joint_nonflat_mass_measure as prior
old=prior.old


def run():
    angles=(0.,.12)
    spatial=[prior.rep(np.eye(3),np.eye(2),np.exp(1j*t)) for t in angles]
    links=np.tile(np.eye(16,dtype=complex),(2,4,1,1))
    links[0,:2]=spatial[0];links[0,2:]=spatial[1]
    e=np.tile(np.eye(10)[0],(4,1));phis=np.tile(prior.mass.car.PHI[0],(4,1))
    base=prior.spatial_data(links,e,phis)
    omegas=[prior.rep(*prior.group(67331+i,.055)) for i in range(2)]
    twisted=links.copy()
    twisted[1,2]=omegas[0];twisted[1,3]=omegas[1]
    q=prior.spatial_data(twisted,e,phis)
    # Independent local gauge frames, applied to both ends, including the seam.
    groups=[prior.group(67341+i,.045) for i in range(4)]
    rs=[prior.rep(*g) for g in groups]
    changed=twisted.copy()
    sites=prior.SITES
    for mu in range(2):
        for i,(t,x) in enumerate(sites):
            dest=((t+1)%2,x) if mu==1 else (t,(x+1)%2)
            j=sites.index(dest)
            changed[mu,i]=rs[i]@twisted[mu,i]@rs[j].conj().T
    et=np.array([prior.vector_rotation(r)@ei for r,ei in zip(rs,e)])
    pt=[]
    for (_,weak,z),ph in zip(groups,phis):
        h=z**3*weak@(ph[:2]+1j*ph[2:4])
        pt.append(np.r_[h.real,h.imag,ph[4]])
    transformed=prior.spatial_data(changed,et,np.array(pt))
    identity_error=float(abs(transformed['weight']/q['weight']-1))
    assert identity_error<3e-11
    shift=q['weight']/base['weight']
    assert abs(shift-1)>.01
    traces=[old.cpair(np.trace(u@u)) for u in spatial]
    assert abs(complex(*traces[0])-complex(*traces[1]))>.1
    closure_errors=[]
    for x in range(2):
        a=sites.index((0,x));b=sites.index((1,x))
        closed=changed[1,a]@changed[1,b]
        closure_errors.append(old.err(closed-rs[a]@omegas[x]@rs[a].conj().T))
    assert max(closure_errors)<3e-13
    return dict(date='2026-10-02',entry_for_round=673,formal_round_complete=False,
        original_signed_weight=old.cpair(base['weight']),twisted_weight=old.cpair(q['weight']),
        twisted_to_identity_ratio=old.cpair(shift),minimum_Wilson_gap=min(base['gap'],q['gap'],transformed['gap']),
        full_gauge_covariance_error=identity_error,temporal_holonomy_conjugation_errors=closure_errors,
        spatial_Wilson_loop_traces=traces,
        orbit_distinction_is_not_a_Gauss_projected_negative_witness=True,
        normalized_physical_two_point_change=old.norm(q['covariance']-base['covariance']),
        no_Haar_projection_no_bosonic_twisted_endpoint_integral=True,
        all_entry_checks_passed=True,
        dependency_hashes={p:hashlib.sha256((BASE/p).read_bytes()).hexdigest() for p in
            ('joint_nonflat_mass_measure.py','joint_gauge_history_kernel.py','research_note_643.md','research_note_672.md')})


if __name__=='__main__':
    result=run();target=HERE/'boundary_holonomy_probe_results.json'
    if '--write-results' in sys.argv:
        with target.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(target.read_text('utf8'))
    print(json.dumps(result,ensure_ascii=False))
