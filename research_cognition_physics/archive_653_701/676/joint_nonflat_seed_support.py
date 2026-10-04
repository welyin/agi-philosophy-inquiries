"""676: fixed boundary seed, exact nonflat ranks and auxiliary null support."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import numpy as np
import joint_gauss_boundary_functional as base
HERE=Path(__file__).resolve().parent;TARGET=HERE/'joint_nonflat_seed_support_results.json'
old=base.old

def load(name):
    path=HERE/'round676_drafts'/name
    spec=importlib.util.spec_from_file_location(name[:-3],path)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    return m

def transfer_and_certificate():
    probe=load('nonflat_boundary_probe.py').run()
    cert=load('exact_flux_seed.py').run()
    cert['mod4_multiplicities']={str(k):v for k,v in cert['mod4_multiplicities'].items()}
    # Clarify the preliminary entry label: the open neighborhood is in ALL
    # retained two-direction links of the original full internal group,
    # not arbitrary extra spacetime links or all four-dimensional backgrounds.
    cert.pop('positive_full_gauge_neighborhood_by_continuity')
    cert['positive_full_retained_two_direction_gauge_neighborhood']=True
    cert['extra_spacetime_directions_not_in_this_witness']=True
    return dict(transfer=probe,certificate=cert)

def support_check():
    sp=base.internal.spin;s4=1j*sp.GAMMA[1]@sp.GAMMA[2]
    s=np.kron(np.kron(np.eye(4),s4),np.eye(16))
    local_identity=old.err(s4.T@base.internal.B+base.internal.B@s4)
    assert local_identity==0
    eig,b=np.linalg.eigh(s4)
    frames=[np.kron(np.kron(np.eye(4),b[:,eig*sg>0]),np.eye(16)) for sg in (-1,1)]
    rng=np.random.default_rng(67641);e=rng.normal(size=(4,10));e/=np.linalg.norm(e,axis=1)[:,None]
    mat=base.fixed_matrices(e,base.mass.car.PHI);jm,jp,m,mb,pair=mat
    rows=[]
    for seed in (67350,67351,67352,None):
        if seed is None:
            links=np.tile(np.eye(16,dtype=complex),(2,4,1,1))
            for i,(t,x) in enumerate(base.prior.SITES):
                for mu,z in enumerate(((-1.)**t if x==1 else 1.,1j if x==1 else 1.)):
                    links[mu,i]=base.prior.rep(np.eye(3),np.eye(2),z)
        else:
            links=np.array([[base.prior.rep(*base.prior.group(seed+8*mu+i,.9)) for i in range(4)] for mu in range(2)])
        u,v,d,h,gap=base.kernel(links)
        counts=[int(sum(np.linalg.eigvalsh(f.conj().T@h@f)<0)) for f in frames]
        imbalance=abs(counts[0]-counts[1])
        a=u.T@m@u
        au=int(sum(np.linalg.svd(a,compute_uv=False)<1e-10))
        full=base.regular(u,v,d,mat)
        nu=int(sum(np.linalg.svd(full['N'],compute_uv=False)<1e-10))
        su=u.conj().T@s@u
        error=old.err(su.T@a+a@su);assert error<3e-13
        assert au>=imbalance and nu>=2*imbalance
        rows.append(dict(seed='exact_flux' if seed is None else seed,Wilson_gap=gap,
            negative_spin_counts=counts,predicted_auxiliary_nullity_lower=imbalance,
            predicted_full_N_nullity_lower=2*imbalance,
            measured_auxiliary_nullity=au,measured_full_N_nullity=nu,
            skew_spin_identity_error=error,raw_weight_not_used_as_sign=old.cpair(full['weight']),
            all_pure_physical_source_coefficients_zero_by_auxiliary_block=bool(imbalance)))
    return dict(local_Clifford_pair_identity_exact_zero=local_identity,rows=rows,
        exact_rank_certificate_available_for_flux_only=True,
        no_candidate_positive_or_negative_Gauss_integral_claim=True)

def run():
    deps=('research_note_659.md','research_note_673.md','research_note_675.md',
          'joint_gauss_boundary_functional.py','round676_drafts/nonflat_boundary_probe.py',
          'round676_drafts/exact_flux_seed.py','round676_drafts/seed_sector_probe_results.json')
    return dict(date='2026-10-02',round=676,tests_run=2,failures=0,errors=0,
        positive_fifth_transfer_and_fixed_seed_obstruction=transfer_and_certificate(),
        auxiliary_support_and_physical_sources=support_check(),
        dependency_hashes={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in deps},
        scope='Original two-direction finite box: fifth transfer remains positive for all retained unitary gauge links, but659 fixed-Jminus image cannot cover an exactly certified nonflat72/56 sector with64/64 seed. Same unused-spin symmetry forces auxiliary and pure-physical-source weight zero on that sector. Explains prior near-zero cases; no refutation of physical marginal, domain-wall theory, RP, original HF, continuum or quantum GR.',
        all_checks_passed=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(round=676,tests_run=2,all_checks_passed=True)))
