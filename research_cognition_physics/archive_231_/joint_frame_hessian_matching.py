"""583: same-action scalar/conformal Hessians and field-dependent freezing.

This is a restricted tangent calculation, not a gauge-fixed gravity determinant.
The Euclidean conformal direction is negative; no convergent Gaussian is asserted.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_curved_quantum_source as original

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_frame_hessian_matching_results.json'


def vacuum():
    L,u,_=original.lattice.scalar.parameters()
    phi=np.array([0.,np.sqrt(u[0]),0.,0.,np.sqrt(u[1])])
    f=float(original.F(phi));K=original.metric(phi)
    hessian=np.zeros((5,5))
    hessian[1,1]=2*L[0,0]*u[0]/f**2
    hessian[4,4]=2*L[1,1]*u[1]/f**2
    hessian[1,4]=hessian[4,1]=2*L[0,1]*np.sqrt(u[0]*u[1])/f**2
    return phi,f,K,hessian,-phi/(6*f)


def vacuum_hessian_check():
    phi,f,K,H,c=vacuum();rows=[];eye=np.eye(5)
    assert original.node_potential(phi)<1e-28 and f>0
    for step in (2e-3,1e-3,5e-4):
        numeric=np.empty((5,5))
        for a in range(5):
            for b in range(5):
                ea=step*eye[a];eb=step*eye[b]
                numeric[a,b]=(original.node_potential(phi+ea+eb)
                    -original.node_potential(phi+ea-eb)-original.node_potential(phi-ea+eb)
                    +original.node_potential(phi-ea-eb))/(4*step**2)
        rows.append(dict(step=step,hessian_error=float(np.max(abs(numeric-H)))))
    assert rows[-1]['hessian_error']<rows[0]['hessian_error']/12
    assert np.linalg.matrix_rank(H,tol=1e-10)==2
    assert np.max(abs(K-6*np.outer(c,c)-np.eye(5)/f))<1e-14
    return dict(phi=phi.tolist(),F=f,M=original.M,K=K.tolist(),potential_hessian=H.tolist(),
        c=c.tolist(),Jordan_background_scale=1/f,finite_difference=rows,
        potential_value=float(original.node_potential(phi)),Hessian_rank=2)


def blocks(p2):
    _,f,K,H,c=vacuum()
    DE=p2*K+H;CJ=p2*np.eye(5)/f+H
    HE=np.zeros((6,6));HE[0,0]=-6*p2;HE[1:,1:]=DE
    HJ=np.zeros((6,6));HJ[0,0]=-6*p2
    HJ[0,1:]=HJ[1:,0]=-6*p2*c;HJ[1:,1:]=CJ
    J=np.eye(6);J[0,1:]=c
    return HE,HJ,J,DE,CJ


def joint_blocks_check():
    rows=[]
    for p2 in (.001,.01,.1,1.,10.):
        HE,HJ,J,DE,CJ=blocks(p2)
        congruence=float(np.max(abs(HJ-J.T@HE@J)))
        schur=HJ[1:,1:]-np.outer(HJ[1:,0],HJ[0,1:])/HJ[0,0]
        schur_error=float(np.max(abs(schur-DE)))
        determinant_error=float(abs(np.linalg.det(HJ)/np.linalg.det(HE)-1))
        inverse_error=float(np.max(abs(np.linalg.inv(HJ)[1:,1:]-np.linalg.inv(DE))))
        omitted_mixing=float(np.linalg.norm(np.linalg.inv(CJ)-np.linalg.inv(DE)))
        assert max(congruence,schur_error,determinant_error)<1e-12
        assert inverse_error<2e-10 and omitted_mixing>1e-4
        inertia=[int(np.sum(np.linalg.eigvalsh(HJ)<0)),int(np.sum(np.linalg.eigvalsh(HJ)>0))]
        assert inertia==[1,5]
        rows.append(dict(p_squared=p2,congruence_error=congruence,Schur_error=schur_error,
            determinant_relative_error=determinant_error,scalar_inverse_block_error=inverse_error,
            scalar_inverse_error_if_mixing_omitted=omitted_mixing,negative_positive_inertia=inertia))
    return dict(rows=rows,zero_momentum_excluded=True,no_gravity_path_integral_claim=True)


def normalized_scalar_check():
    _,f,K,H,c=vacuum()
    eig,O=np.linalg.eigh(K);inverse_root=(O/np.sqrt(eig))@O.T
    ME=inverse_root@H@inverse_root;MJ=f*H
    massE=np.linalg.eigvalsh(ME);massJ=np.linalg.eigvalsh(MJ)
    rows=[]
    for p2 in (1e-5,.001,.01,.1,1.,100.,1e5):
        _,_,_,DE,CJ=blocks(p2)
        raw=float(np.exp(np.linalg.slogdet(DE)[1]-np.linalg.slogdet(CJ)[1]))
        lemma=float(1+6*p2*c@np.linalg.solve(CJ,c))
        normalized=float(np.exp(np.linalg.slogdet(p2*np.eye(5)+ME)[1]
            -np.linalg.slogdet(p2*np.eye(5)+MJ)[1]))
        assert abs(raw-lemma)<1e-11 and abs(normalized-f/original.M*raw)<1e-11
        assert 1<raw<original.M/f+1e-12
        assert f/original.M-1e-12<normalized<1
        rows.append(dict(p_squared=p2,raw_scalar_determinant_ratio=raw,
            lemma_error=abs(raw-lemma),kinetically_normalized_ratio=normalized))
    assert abs(rows[-1]['kinetically_normalized_ratio']-1)<1e-6
    assert np.max(abs(massE-massJ))>1e-3
    return dict(Einstein_partial_mass_squared=massE.tolist(),Jordan_partial_mass_squared=massJ.tolist(),
        raw_high_momentum_limit=original.M/f,normalized_high_momentum_limit=1,
        normalized_low_momentum_limit=f/original.M,rows=rows,
        partial_scalar_loop_trace_M_squared_E=float(np.trace(ME@ME)),
        partial_scalar_loop_trace_M_squared_J=float(np.trace(MJ@MJ)),
        spectra_are_not_full_physical_poles=True)


def profiles(amplitude,direction,mode=2,nodes=2048):
    phi0,_,_,_,_=vacuum()
    x=np.arange(nodes)*2*np.pi/nodes
    wave=np.cos(mode*x);wavep=-mode*np.sin(mode*x);wavepp=-mode**2*wave
    phi=phi0+amplitude*wave[:,None]*direction
    phip=amplitude*wavep[:,None]*direction
    phipp=amplitude*wavepp[:,None]*direction
    f=original.F(phi);fp=-np.sum(phi*phip,axis=-1)/3
    fpp=-np.sum(phip*phip+phi*phipp,axis=-1)/3
    assert np.min(f)>0
    return wave,wavep,wavepp,phi,phip,f,fp,fpp


def actions(amplitude,direction,metric_amplitude=0.,mode=2,nodes=2048):
    # a is the coefficient of eta*cos(nx) in sigma_J, about gJ0=delta/F0.
    _,f0,_,_,_=vacuum()
    wave,wp,wpp,phi,phip,f,fp,fpp=profiles(amplitude,direction,mode,nodes)
    sigmaJ=-.5*np.log(f0)+amplitude*metric_amplitude*wave
    sJp=amplitude*metric_amplitude*wp;sJpp=amplitude*metric_amplitude*wpp
    sigmaE=sigmaJ+.5*np.log(f)
    sEp=sJp+.5*fp/f;sEpp=sJpp+.5*(fpp/f-(fp/f)**2)
    e2E=np.exp(2*sigmaE);e2J=np.exp(2*sigmaJ)
    scalar_gradient=np.einsum('...a,...ab,...b->...',phip,original.metric(phi),phip)
    U=original.node_potential(phi)
    L,u,_=original.lattice.scalar.parameters()
    d=np.stack((np.sum(phi[:,:4]**2,axis=-1)-u[0],phi[:,4]**2-u[1]),axis=-1)
    V=np.einsum('...a,ab,...b->...',d,L,d)/4
    integral=lambda a:float(2*np.pi*np.mean(a))
    einstein=integral(3*e2E*(sEpp+sEp*sEp)+.5*e2E*scalar_gradient+e2E**2*U)
    einstein_ibp=integral(-3*e2E*sEp*sEp+.5*e2E*scalar_gradient+e2E**2*U)
    jordan=integral(3*f*e2J*(sJpp+sJp*sJp)+.5*e2J*np.sum(phip*phip,axis=-1)+e2J**2*V)
    fixedE=integral(.5*scalar_gradient+U)
    fixedJ=integral(np.sum(phip*phip,axis=-1)/(2*f0)+V/f0**2)
    # Pull back fixed gE=delta: sigmaJ=-log(F(phi))/2, not constant sigmaJ.
    sigp=-.5*fp/f;sigpp=-.5*(fpp/f-(fp/f)**2)
    pulled_fixedE=integral(3*(sigpp+sigp*sigp)
        +np.sum(phip*phip,axis=-1)/(2*f)+V/f**2)
    return dict(E=einstein,E_integrated_by_parts=einstein_ibp,J=jordan,
        fixed_E=fixedE,fixed_J=fixedJ,pulled_fixed_E=pulled_fixedE,min_F=float(np.min(f)))


def frozen_action_check():
    _,_,_,_,c=vacuum();v=np.array([.13,.8,-.11,.07,.6]);mode=2
    _,_,_,DE,CJ=blocks(mode**2)
    expectedE=float(np.pi*v@DE@v);expectedJ=float(np.pi*v@CJ@v)
    gap=float(-6*np.pi*mode**2*(c@v)**2)
    rows=[]
    for eta in (.04,.02,.01,.005):
        plus=actions(eta,v);minus=actions(-eta,v)
        numericE=(plus['fixed_E']+minus['fixed_E'])/eta**2
        numericJ=(plus['fixed_J']+minus['fixed_J'])/eta**2
        for r in (plus,minus):
            assert abs(r['E']-r['E_integrated_by_parts'])<1e-13
            assert abs(r['E']-r['J'])<1e-13
            assert abs(r['J']-r['fixed_J'])<1e-13
            assert abs(r['pulled_fixed_E']-r['fixed_E'])<1e-13
        rows.append(dict(amplitude=eta,second_E=numericE,second_J=numericJ,
            error_E=abs(numericE-expectedE),error_J=abs(numericJ-expectedJ),
            gap_error=abs(numericJ-numericE-gap),nonlinear_frame_action_error=abs(plus['E']-plus['J']),
            pullback_freeze_error=abs(plus['pulled_fixed_E']-plus['fixed_E'])))
    assert rows[-1]['error_E']<rows[0]['error_E']/50
    assert rows[-1]['error_J']<rows[0]['error_J']/50
    assert abs(expectedJ-expectedE-gap)<1e-13 and abs(gap)>.05
    return dict(direction=v.tolist(),mode=mode,expected_second_E=expectedE,
        expected_second_J=expectedJ,expected_difference=gap,rows=rows,
        periodic_integral_transverse_coordinate_volume=1)


def joint_action_check():
    rows=[];v=np.array([.13,.8,-.11,.07,.6]);mode=2
    HE,HJ,J,_,_=blocks(mode**2)
    for a in (.17,-.21):
        vec=np.r_[a,v];expected=float(np.pi*vec@HJ@vec)
        assert abs(expected-np.pi*(J@vec)@HE@(J@vec))<1e-13
        values=[]
        for eta in (.04,.02,.01,.005):
            plus=actions(eta,v,a);minus=actions(-eta,v,a)
            curvature=(plus['J']+minus['J'])/eta**2
            frame_error=max(abs(r['E']-r['J']) for r in (plus,minus))
            ibp_error=max(abs(r['E']-r['E_integrated_by_parts']) for r in (plus,minus))
            assert max(frame_error,ibp_error)<1e-13
            values.append(dict(amplitude=eta,actual_second_variation=curvature,
                error=abs(curvature-expected),frame_action_error=frame_error))
        assert values[-1]['error']<values[0]['error']/50
        rows.append(dict(Jordan_conformal_amplitude=a,predicted_second_variation=expected,rows=values))
    return dict(rows=rows,nonlinear_original_Jordan_action_checked=True)


def run():
    evidence=dict(vacuum_hessian=vacuum_hessian_check(),joint_blocks=joint_blocks_check(),
        normalized_scalar=normalized_scalar_check(),frozen_action=frozen_action_check(),
        joint_action=joint_action_check())
    deps=('joint_curved_quantum_source.py','joint_gauss_einstein_initial_data.py',
        'research_note_553.md','research_note_580.md','research_note_582.md','research_round_582_checks.json')
    return dict(round=583,tests_run=len(evidence),failures=0,errors=0,checks=list(evidence),evidence=evidence,
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope='same four-dimensional leading scalar-tensor action and original H5 matter at its common flat vacuum; restricted scalar/conformal Hessian, normalized frozen scalar blocks and nonlinear Weyl substitution; not a full gauge-fixed gravitational determinant, convergent Euclidean integral, physical pole spectrum or graph-continuum quantum matching')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run();payload=json.dumps(result,ensure_ascii=False,indent=2)+'\n'
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(payload)
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps({k:result[k] for k in ('round','tests_run','failures','errors')}))
