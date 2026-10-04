"""682: actual678 auxiliary deformations and the complete boundary-source test.
Berezin inverse blocks are never interpreted as positive classical variances.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import numpy as np
import joint_local_source_lift as body
HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_bulk_fluctuation_matching_results.json'
base=body.base
err=body.err

def entry(name,file):
    spec=importlib.util.spec_from_file_location(name,HERE/'round682_drafts'/file)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module

probe=entry('body_probe682','body_fluctuation_probe.py')
control=entry('conditional_probe682','conditional_bulk_probe.py')

def mixed_deformation():
    links,e,phis=body.prior.fixture();mat=base.fixed_matrices(e,phis)
    n=len(mat[2]);g5=mat[1]@mat[1].T-mat[0]@mat[0].T
    _,_,_,h,_=base.kernel(links)
    lifted=body.lift(g5@h,mat,lam=.37,layers=1)
    k=lifted['K'];O=lifted['out'];m=len(k)
    V=np.vstack((np.eye(n),np.eye(n)));W=np.vstack((np.eye(n),np.zeros((n,n))))
    delta=.09;Gamma=mat[3]
    sigma=delta*(V@Gamma@W.T+W@Gamma@V.T)
    R=k@sigma@k.T
    one_leg=float(np.linalg.norm(O@sigma,2))
    visible=float(np.linalg.norm(O@sigma@O.T,2))
    assert one_leg>.1 and visible<2e-14
    full=lifted['N'].copy();full[2*n+m:,2*n+m:]+=R
    Y=lifted['Y'];Q=full[2*n:,2*n:];Qi=np.linalg.inv(Q)
    mix=full[:2*n,2*n:];Yq=Y[:,2*n:]
    eff=full[:2*n,:2*n]+mix@Qi@mix.T
    phi=Y[:,:2*n]+Yq@Qi@mix.T
    contact=Yq@Qi@Yq.T;ref=lifted['reference']
    errors=dict(action=err(eff-ref['N']),source=err(phi-ref['Phi']),
                contact=err(contact))
    assert max(errors.values())<4e-12
    rng=np.random.default_rng(68261)
    sources=rng.normal(size=(n,2))+1j*rng.normal(size=(n,2))
    sources/=np.linalg.norm(sources,axis=0)
    source_rows=[]
    for count in (0,2):
        z=sources[:,:count]
        observed=probe.corrected_large_coefficient(body.augmented(full,Y,z),k)
        target=base.pf(body.augmented(ref['N'],ref['Phi'],z))
        rel=float(abs(observed-target)/max(abs(observed),abs(target),1e-250))
        assert rel<3e-9
        source_rows.append(dict(source_count=count,full_compensated=base.old.cpair(observed),
            original=base.old.cpair(target),relative_error=rel))
    # Actual independent local transformations of color, weak and hypercharge.
    groups=[base.prior.group(68240+i,.37) for i in range(4)]
    lt,et,pt,rs=base.transform(links,e,phis,groups)
    mt=base.fixed_matrices(et,pt)
    _,_,_,ht,_=base.kernel(lt)
    kt,_,ot=body.blocks(g5@ht,g5,.23,1)
    r4=np.zeros((n,n),complex)
    for i,r in enumerate(rs):r4[64*i:64*(i+1),64*i:64*(i+1)]=np.kron(np.eye(4),r)
    rm=np.kron(np.eye(2),r4)
    sigmat=delta*(V@mt[3]@W.T+W@mt[3]@V.T)
    Rt=kt@sigmat@kt.T
    gauge=dict(K=err(kt-rm@k@rm.conj().T),
        O=err(ot-r4@O@rm.conj().T),
        Sigma=err(sigmat-rm@sigma@rm.T),
        eta_pair=err(Rt-rm@R@rm.T),
        invariant_visible_zero=err(ot@sigmat@ot.T))
    assert max(gauge.values())<2e-12
    return dict(mass=.37,one_leg_inverse_block=one_leg,
        two_leg_visible_inverse_block=visible,Schur_errors=errors,
        source_rows=source_rows,original_Gauss_covariance_errors=gauge,
        nonzero_one_leg_block_does_not_break_conditional_physical_identity=True)

def run():
    inherited_probe=probe.run()
    assert inherited_probe==json.loads(probe.TARGET.read_text('utf8'))
    bulk_control=control.run()
    assert bulk_control==json.loads(control.TARGET.read_text('utf8'))
    deps=('research_note_643.md','research_note_673.md','research_note_675.md',
          'research_note_677.md','research_note_678.md','research_note_679.md',
          'research_note_680.md','research_note_681.md','joint_local_source_lift.py',
          'round682_drafts/body_fluctuation_probe.py',
          'round682_drafts/body_fluctuation_probe_results.json',
          'round682_drafts/conditional_bulk_probe.py',
          'round682_drafts/conditional_bulk_probe_results.json',
          'round682_drafts/conditional_bulk_entry.md')
    return dict(date='2026-10-02',round=682,tests_run=2,failures=0,errors=0,
        original_auxiliary_body=dict(direct_probe=inherited_probe,mixed=mixed_deformation()),
        complete_positive_bulk_and_boundary_control=bulk_control,
        dependency_hashes={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in deps},
        scope=dict(original_678_full_16_channel_finite_body=True,
            conditional_matching_criterion_O_Sigma_O_transpose_zero=True,
            matching_for_all_physical_sources_and_finite_mass=True,
            pointwise_invisible_family_preserves_original_complete_average=True,
            visible_deformation_not_original_model_failure=True,
            no_Grassmann_to_positive_noise_inference=True,
            positivity_and_original_HF_identity_not_proved=True,
            whole_goal_and_old_spatial_results_unchanged=True))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(round=682,tests_run=2,all_checks_passed=True,
        mixed=result['original_auxiliary_body']['mixed']['Schur_errors'])))

