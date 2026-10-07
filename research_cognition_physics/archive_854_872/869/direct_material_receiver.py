"""869: original loop source menu -> declared local receiver coupling.
Finite differential/probe matrices calibrate ordering and response identities.
They are not a numerical curved-spacetime field theory or continuum CP proof.
"""
from pathlib import Path
import argparse,itertools,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'867'))
import shared_weyl_readout as prior
TARGET=HERE/'direct_material_receiver_results.json'
I=np.eye(2,dtype=complex)
X=np.array([[0,1],[1,0]],complex)
Y=np.array([[0,-1j],[1j,0]],complex)
Z=np.diag([1.,-1.]).astype(complex)
PAULIS=np.array([X,Y,Z])
GAINS=np.array([6.,6.,.75])

def kronall(*items):
    out=np.array([[1.]],complex)
    for a in items:out=np.kron(out,a)
    return out

def source_calibration():
    # Same oblique finite differential projector used by785.
    n=7;eye=np.eye(n)
    d=(np.roll(eye,1,axis=1)-np.roll(eye,-1,axis=1))/2
    d=np.diag(np.linspace(.6,1.4,n))@d+np.diag(np.linspace(-.2,.3,n))
    k=np.vstack((d,eye));f=np.hstack((eye,-d))
    l=np.hstack((np.zeros((n,n)),eye))+(.17*eye+.09*d)@f
    pi=np.eye(2*n)-k@l
    rng=np.random.default_rng(869)
    j=f.T@rng.normal(size=(n,3))
    chi=-2*(1.2+.15*np.cos(np.arange(n)[:,None]+np.arange(3)[None,:]))
    q=rng.normal(size=(n,3))
    recovered=[];ward=[];naive=[];wrong=[];duality=[]
    for a in range(3):
        b=j[:,a]/np.tile(chi[:,a],2)
        B=np.vstack((np.diag(b[:n]),np.diag(b[n:])))
        recovered.append(pi.T@B@chi[:,a])
        full=pi.T@B@q[:,a]
        ward.append(float(np.max(abs(k.T@full))))
        naive.append(float(np.max(abs(k.T@(B@q[:,a])))))
        # Moving the susceptibility outside the differential adjoint is wrong.
        bad=np.tile(chi[:,a],2)*(pi.T@b)
        wrong.append(float(np.linalg.norm(bad-j[:,a])))
        u=rng.normal(size=2*n)
        duality.append(float(abs(u@full-(pi@u)@(B@q[:,a]))))
    recovered=np.stack(recovered,axis=1)
    error=float(np.max(abs(recovered-j)))
    assert error<2e-14 and max(ward)<2e-14 and max(duality)<2e-14
    assert min(naive)>.05 and min(wrong)>.01
    return dict(source_recovery_error=error,gauge_Ward_error=max(ward),
        source_field_duality_error=max(duality),
        unprojected_local_coupling_Ward_defects=naive,
        misplaced_susceptibility_source_errors=wrong,
        nonzero_susceptibility_margin=float(np.min(abs(chi))),
        scope='finite jet-ordering calibration, not original continuum PDE')

def original_menu():
    ry,ay,w=prior.finite_menu();probe=prior.probe
    with probe.ResearchRuntime(probe.Layout()).installed():
        import joint_reference_constraint_strata as old
        U,dU,_,_=probe.original.matrices(old,probe.original.inherited.constants(),.4)
    phase=np.exp(-.74j);r=phase*np.trace(U);dr=phase*np.trace(dU)
    base=np.array([r.real,r.imag,abs(r)**2-1])
    tangent=np.array([dr.real,dr.imag,2*np.real(np.conj(r)*dr)])
    C=w[:,None]*np.stack((2*probe.ETA_R*ry.real,2*probe.ETA_R*ry.imag,probe.ETA_A*ay),axis=1)
    p=w+C@base
    values=np.stack((ry.real/probe.ETA_R,ry.imag/probe.ETA_R,ay/probe.ETA_A,1+ay/probe.ETA_A))
    return C,p,values,tangent

def unitary(h,t):
    ev,u=np.linalg.eigh(h)
    return (u*np.exp(-1j*t*ev))@u.conj().T

def run():
    projection=source_calibration()
    C,p,values,tangent=original_menu()
    bits=np.array(list(itertools.product((1.,-1.),repeat=3)))
    decoder=p[:,None]+(C/GAINS)@bits.T
    assert np.min(decoder/p[:,None])>.5
    stochastic=float(np.max(abs(decoder.sum(axis=0)-1)))
    # Three actual prepared flavor records; baseline is uniform on eight words.
    first_word=(bits@(GAINS[:,None]*PAULIS.reshape(3,4)))/8
    first_effect=(decoder@first_word).reshape(len(p),2,2)
    target=np.einsum('ja,amn->jmn',C,PAULIS)
    all_source_error=float(np.max(abs(first_effect-target)))
    probability_derivative=decoder@(bits@(GAINS*tangent)/8)
    original_probability_error=float(np.max(abs(probability_derivative-C@tangent)))
    assert stochastic<2e-13 and all_source_error<2e-15 and original_probability_error<1e-18
    fisher=float(np.sum(probability_derivative**2/p))
    assert fisher>0
    # A two-pulse diagnostic of receiver dynamics with the same first sources.
    # Extra time-ordering dependence at finite strength must not be called868.
    kick=np.array([.09*Z,.06*X,.07*Y])
    parts=(.5*PAULIS+kick,.5*PAULIS-kick)
    hs=[]
    for part in parts:
        h=np.zeros((16,16),complex)
        for a in range(3):
            local=[I,I,I];local[a]=-Y/2
            h+=GAINS[a]*kronall(part[a],*local)
        hs.append(h)
    eta=np.ones((8,1),complex)/np.sqrt(8)
    isometry=np.kron(I,eta)
    def output(t):
        U=unitary(hs[1],t)@unitary(hs[0],t)
        out=(U@isometry).reshape(2,8,2).transpose(1,0,2)
        word_effects=out.conj().transpose(0,2,1)@out
        effects=np.einsum('jz,zmn->jmn',decoder,word_effects)
        ideal=p[:,None,None]*I
        for a in range(3):
            ideal=ideal+(C[:,a]/GAINS[a])[:,None,None]*prior.herm_fun(GAINS[a]*t*PAULIS[a],np.sin)
        return effects,float(np.max(abs(U.conj().T@U-np.eye(16)))),float(np.max(abs(effects-ideal)))
    rows=[]
    for t in (.08,.04,.02,.01):
        e,ue,diff=output(t);em,_,_=output(-t)
        derivative=(e-em)/(2*t)
        err=float(np.max(abs(np.einsum('kj,jmn->kmn',values,derivative-target))))
        norm=float(np.max(abs(e.sum(axis=0)-I)))
        lower=float(np.min(np.linalg.eigvalsh(e)/p[:,None]))
        assert norm<3e-13 and ue<3e-14 and lower>.5
        assert diff>1e-12
        rows.append(dict(strength=t,all_menu_first_response_error=err,
            normalization_error=norm,relative_effect_lower=lower,
            unitarity_error=ue,difference_from_exact868_effects=diff))
    assert all(3.7<rows[i]['all_menu_first_response_error']/rows[i+1]['all_menu_first_response_error']<4.3 for i in range(3))
    return dict(round=869,date='2026-10-06',formal_reports=869,
        cumulative_numbered_groups=3654,fresh_numbered_groups=1,all_checks_passed=True,
        differential_source_calibration=projection,
        same_original_final_labels=len(p),native_binary_receivers=3,native_joint_records=8,
        finite_decoder_stochastic_error=stochastic,
        complete_first_source_transfer_error=all_source_error,
        original_physical_family_probability_derivative_error=original_probability_error,
        original_central_anchor_probability_Fisher=fisher,
        receiver_dynamics_calibration=rows,
        continuum_result='analytic: lowest read-coupling and sqrt(hbar) response is c_alpha Phi(j_alpha)',
        nonlinear_completion='original785 formal relational slice; every fixed field order finite jet',
        added_species=0,added_interaction_and_calibrated_preparation=True,
        background_and_free_hessian_preserved=True,
        exact868_sine_instrument_or_poststate_realized=False,
        all_order_first_coupling_response_matched=False,
        autonomy_or_finite_coupling_continuum_realization=False,
        full_goal_completed=False)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');a=ap.parse_args()
    r=run()
    if a.write:
        assert not TARGET.exists()
        TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
