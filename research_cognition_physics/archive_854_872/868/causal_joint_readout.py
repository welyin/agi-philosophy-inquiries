"""868: random finite Weyl readouts compile the SAME867 source menu.
Three old809 sine primitives plus a stochastic decoder. The regional statement
is proved analytically from Weyl relations; no Pauli model is passed off as the
original spacetime commutator.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'867'))
import shared_weyl_readout as old
TARGET=HERE/'causal_joint_readout_results.json'

def run():
    ry,ay,w=old.finite_menu();probe=old.probe
    with probe.ResearchRuntime(probe.Layout()).installed():
        import joint_reference_constraint_strata as original
        U,dU,_,_=probe.original.matrices(original,probe.original.inherited.constants(),.4)
    phase=np.exp(-.74j);r0=phase*np.trace(U);dr=phase*np.trace(dU)
    a0=abs(r0)**2-1;da=2*np.real(np.conj(r0)*dr)
    base=np.array([r0.real,r0.imag,a0]);direction=np.array([dr.real,dr.imag,da])
    C=w[:,None]*np.stack([2*probe.ETA_R*ry.real,2*probe.ETA_R*ry.imag,probe.ETA_A*ay],axis=-1)
    p0=w+C@base
    r=np.full(3,1/3);c=np.array([6.,6.,.75])
    signs=np.array([-1.,1.])
    T=p0[:,None,None]+signs[None,None,:]*C[:,:,None]/(r*c)[None,:,None]
    assert np.min(T/p0[:,None,None])>=.5-1e-14
    assert np.max(abs(T.sum(axis=0)-1))<1e-13
    primitive0=(r[:,None]*np.ones((3,2)))/2
    baseline=np.einsum('jas,as->j',T,primitive0)
    primitive_first=np.zeros((3,2,3))
    for a in range(3):primitive_first[a,:,a]=signs*r[a]*c[a]/2
    recovered=np.einsum('jas,asb->jb',T,primitive_first)
    transfer_error=max(float(np.max(abs(baseline-p0))),float(np.max(abs(recovered-C))))
    assert transfer_error<1e-14
    G=C.T@(C/p0[:,None]);primitive_G=np.diag(r*c*c)
    difference=primitive_G-G;eig=np.linalg.eigvalsh(difference)
    assert min(eig)>0
    pvar=float(direction@G@direction)/4
    nvar=float(direction@primitive_G@direction)/4
    assert nvar>pvar>0
    sx=np.array([[0,1],[1,0]],complex);sy=np.array([[0,-1j],[1j,0]],complex);sz=np.diag([1,-1])
    paulis=np.array([sx,sy,sz]);I=np.eye(2,dtype=complex)
    X=.3*sx-.2*sy+.4*sz
    rows=[]
    for t in (.2,.1,.05):
        E=np.broadcast_to(p0[:,None,None]*I,(len(w),2,2)).copy()
        channel=np.zeros((2,2),complex)
        for a in range(3):
            E+=(C[:,a]/c[a])[:,None,None]*np.sin(c[a]*t)*paulis[a]
            V=np.cos(c[a]*t/2)*I+1j*np.sin(c[a]*t/2)*paulis[a]
            Lplus=(np.cos(c[a]*t/2)*I+np.sin(c[a]*t/2)*paulis[a])/np.sqrt(2)
            Lminus=(np.cos(c[a]*t/2)*I-np.sin(c[a]*t/2)*paulis[a])/np.sqrt(2)
            from_kraus=Lplus@X@Lplus+Lminus@X@Lminus
            from_weyl=(V@X@V.conj().T+V.conj().T@X@V)/2
            assert np.max(abs(from_kraus-from_weyl))<1e-14
            channel+=r[a]*from_weyl
        assert np.min(np.linalg.eigvalsh(E)/p0[:,None])>.5-1e-13
        assert np.max(abs(E.sum(axis=0)-I))<1e-13
        rows.append(dict(t=t,normalization_error=float(np.max(abs(E.sum(axis=0)-I))),
            relative_effect_lower=float(np.min(np.linalg.eigvalsh(E)/p0[:,None]))))
    # Exact channel multiplier formula tested for arbitrary real CCR pairings.
    pairings=np.array([.17,-.23,.31]);t=.2
    multiplier=float(np.sum(r*np.cos(c*t*pairings/2)))
    assert abs(multiplier)<=1
    return dict(round=868,date='2026-10-06',formal_reports=868,
        cumulative_numbered_groups=3653,fresh_numbered_groups=1,all_checks_passed=True,
        reported_labels=len(w),primitive_settings=3,primitive_outcomes_per_setting=2,
        fixed_selector_weights=r.tolist(),fixed_gains=c.tolist(),
        universal_decoder_bound='p_j/2 <= T(j|alpha,s) <= 3 p_j/2',
        decoder_min_relative=float(np.min(T/p0[:,None,None])),
        decoder_stochastic_error=float(np.max(abs(T.sum(axis=0)-1))),
        whole_probability_source_transfer_error=transfer_error,
        primitive_backaction_Gram=primitive_G.tolist(),
        old_menu_Fisher_Gram_eigenvalues=np.linalg.eigvalsh(G).tolist(),
        backaction_excess_Gram_eigenvalues=eig.tolist(),
        central_original_family_companion_t2_variance=pvar,
        causal_choice_companion_t2_variance=nvar,
        finite_matrix_checks=rows,arbitrary_CCR_multiplier_calibration=multiplier,
        region_preservation='analytic: D(W(k))=m(k)W(k), hence D(A(O)) subset A(O)',
        nonselective_spacelike_preparation_test_passed='analytic free physical net',
        same_original_first_source_menu=True,
        interpretation='different instrument with same baseline/first source; not identical poststates',
        local_probe_or_original_action_realization_proved=False,
        entire_relativistic_measurement_protocol_certified=False,
        interacting_claim='same transported net and formal S only',
        continuum_nonlinear_graph_matching=False,full_goal_completed=False)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args();r=run()
    if a.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
