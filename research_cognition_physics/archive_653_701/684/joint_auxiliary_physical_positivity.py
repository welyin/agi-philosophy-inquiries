"""684: scoped copy-repair obstruction and an exactly compensated spectator.
The original interacting physical Gauss/S9 RP problem remains unresolved.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import numpy as np
import joint_local_source_lift as body
HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_auxiliary_physical_positivity_results.json'
base=body.base;err=body.err

def copy_class():
    path=HERE/'round684_drafts/compensation_involution_probe.py'
    spec=importlib.util.spec_from_file_location('entry684',path)
    entry=importlib.util.module_from_spec(spec);spec.loader.exec_module(entry)
    inherited=entry.run()
    assert inherited==json.loads(entry.TARGET.read_text('utf8'))
    J=np.array([[0.,1.],[-1.,0.]])
    a=.2;beta=a/(4+a*a)**2;rng=np.random.default_rng(68411);rows=[]
    for copies in (2,4,6):
        o,_=np.linalg.qr(rng.normal(size=(copies,copies)))
        A=o@np.kron(np.eye(copies//2),J)@o.T
        assert err(A.T@A-np.eye(copies))<2e-14
        assert err(A@A+np.eye(copies))<2e-14
        # Restriction to original W=0 mode and one fixed internal component.
        gram=beta*np.kron(A,J);ev=np.linalg.eigvalsh(gram)
        assert err(gram-gram.T)<1e-14
        assert err(ev[:copies]+beta)<2e-14 and err(ev[copies:]-beta)<2e-14
        rows.append(dict(copies=copies,real_copy_rotation=o.tolist(),
            antisymmetric_orthogonal_copy_matrix=A.tolist(),
            negative_modes=int(np.sum(ev<0)),positive_modes=int(np.sum(ev>0)),
            minimum_eigenvalue=float(ev[0]),analytic_negative_eigenvalue=-beta))
    return dict(executed_entry=inherited,real_orthogonal_factorized_copy_class=rows,
        arbitrary_complex_or_background_dependent_reflections_not_excluded=True)

def unit_spectator():
    links,e,phis=body.prior.fixture();mat=base.fixed_matrices(e,phis)
    _,_,_,h,_=base.kernel(links);n=len(h)
    g5=mat[1]@mat[1].T-mat[0]@mat[0].T
    k,_,_=body.blocks(g5@h,g5,.2,1);G=k.conj().T@k;m=len(G)
    # Two complex bosons with G, and two independent complex Grassmann pairs G.
    # The per-pair Berezin orientation is fixed to integral=det G.
    anti=np.block([[np.zeros_like(G),-G.T],[G,np.zeros_like(G)]])
    phase,logpf=body.logpf(anti);sign=(-1)**(m*(m+1)//2)
    G2=np.kron(np.eye(2),G);sg2,ld2=np.linalg.slogdet(G2)
    factor=(phase/sign)**2/sg2*np.exp(2*logpf-ld2)
    assert abs(factor-1)<3e-10
    eps=body.prior.regulate(h,g5,.23,1)[0]
    original=body.prior.soft(eps,mat,.37)
    rng=np.random.default_rng(68412)
    z=rng.normal(size=(n,4))+1j*rng.normal(size=(n,4));z/=np.linalg.norm(z,axis=0)
    rows=[]
    for sources in (0,2,4):
        old=base.pf(body.augmented(original['N'],original['Phi'],z[:,:sources]))
        kept=old*factor
        residual=float(abs(kept-old)/max(abs(old),1e-250))
        assert residual<3e-10
        rows.append(dict(sources=sources,original=base.old.cpair(old),
            fully_compensated=base.old.cpair(kept),relative_identity_error=residual))
    # Separate reflection test in the original FREE physical branch, where661/658
    # already prove a normalized physical state. The spectator still has a bad full algebra.
    free=np.tile(np.eye(16,dtype=complex),(2,4,1,1));_,_,_,hf,_=base.kernel(free)
    kf,_,_=body.blocks(g5@hf,g5,.2,1);Cf=np.linalg.inv(kf.conj().T@kf)
    permutation=np.eye(4)[[2,3,0,1]]
    ell=np.kron(permutation,np.kron(base.internal.spin.GAMMA[3],np.eye(16)))
    R=np.kron(np.array([[0.,1.],[1.,0.]]),g5@ell)
    pos=np.r_[np.arange(128,256),np.arange(128,256)+256]
    B=(R@Cf)[np.ix_(pos,pos)]
    Q=np.block([[np.zeros_like(B),B],[-B,np.zeros_like(B)]])
    u=np.zeros(128);u[[0,64]]=1/np.sqrt(2)
    v=np.zeros(512);v[:128]=u/np.sqrt(2);v[384:]=-u/np.sqrt(2)
    norm=np.vdot(v,Q@v);beta=.2/(4+.2**2)**2
    assert abs(norm+beta)<1e-14
    return dict(original_nonflat_full_16_channels=True,extra_boson_copies=2,
        extra_Grassmann_pairs=2,Gram_dimension=m,
        independently_computed_spectator_factor=base.old.cpair(factor),source_rows=rows,
        full_S9_free_physical_state_positivity_inherited_from_658_661=True,
        normalized_extended_free_auxiliary_norm=base.old.cpair(norm),
        physical_all_source_identity_pointwise_in_background=True,
        no_new_physical_degrees_of_freedom_claimed=True,
        no_dynamic_Gauss_physical_negative_witness_claimed=True)

def run():
    deps=('research_note_658.md','research_note_661.md','research_note_675.md',
          'research_note_678.md','research_note_679.md','research_note_680.md',
          'research_note_682.md','research_note_683.md','joint_local_source_lift.py',
          'round684_drafts/compensation_involution_probe.py',
          'round684_drafts/compensation_involution_probe_results.json',
          'round684_drafts/compensation_involution_entry.md','round684_drafts/entry_checks.json')
    return dict(date='2026-10-02',round=684,tests_run=2,failures=0,errors=0,
        copy_reflection_class=copy_class(),complete_compensated_extension=unit_spectator(),
        dependency_hashes={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in deps},
        scope=dict(real_factorized_finite_copy_repairs_excluded=True,
            all_reflections_or_all_extensions_excluded=False,
            all_physical_sources_preserved_by_declared_unit_spectator=True,
            full_auxiliary_RP_not_necessary_for_original_physical_RP=True,
            necessity_counterexample_uses_original_proven_free_branch=True,
            original_dynamic_Gauss_RP_decided=False,
            original_HF_continuum_or_quantum_GR_completed=False,
            old_space_and_full_goal_unchanged=True))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(round=684,tests_run=2,all_checks_passed=True,
        spectator_factor=result['complete_compensated_extension']['independently_computed_spectator_factor'])))
