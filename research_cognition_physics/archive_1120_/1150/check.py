"""1150 finite identities and diagnostic history; not a continuum field simulation."""
from pathlib import Path
import argparse, hashlib, importlib.util, json, math
import numpy as np
BASE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('density1150',BASE.parent/'_shared/after1148_local_density_check.py')
q=importlib.util.module_from_spec(spec);spec.loader.exec_module(q)

def material_group():
    v=q.radial(128)
    beta=np.diag([1,1,-1,-1])
    sig=[np.array([[0,1],[1,0]],complex),np.array([[0,-1j],[1j,0]],complex),np.diag([1,-1])]
    rows=[]
    for r,theta,phi in [(.13,.31,.7),(.4,1.3,1.7),(.81,2.1,-.2),(.97,2.8,.8)]:
        n=[math.sin(theta)*math.cos(phi),math.sin(theta)*math.sin(phi),math.cos(theta)]
        sn=sum(a*b for a,b in zip(n,sig))
        F=float(q.j0([v['q']*r])[0]);G=float(v['b']*q.j1([v['q']*r])[0])
        U=math.sqrt(v['N_squared'])*np.vstack([F*np.eye(2),1j*G*sn])
        V=math.sqrt(v['N_squared'])*np.vstack([-1j*G*sn,F*np.eye(2)])
        both=np.column_stack([U,V])
        rho=v['N_squared']*(F*F-G*G)
        err=float(np.linalg.norm(both.conj().T@beta@both-rho*beta,2))
        assert err<2e-14
        rows.append(dict(radius=r,scalar_block_residual=err,pair_block_norm=float(np.linalg.norm(U.conj().T@beta@V))))
    # Both spinor columns obey the same MIT boundary.
    F=float(q.j0([v['q']])[0]);G=float(v['b']*q.j1([v['q']])[0])
    sn=sig[2]
    alpha=np.block([[np.zeros((2,2)),sn],[sn,np.zeros((2,2))]])
    B=np.eye(4)+1j*beta@alpha
    spinors=np.column_stack([np.vstack([F*np.eye(2),1j*G*sn]),np.vstack([-1j*G*sn,F*np.eye(2)])])
    bc=float(np.linalg.norm(B@spinors))
    assert bc<1e-13
    return dict(kind='actual_positive_negative_ground_scalar_blocks',passed=True,
                rows=rows,MIT_boundary_residual=bc,pair_vertices_vanish=True,
                strict_subball_locality_claimed=False)

def budget_group():
    lam=1/600000
    amin=3/20
    rho=math.sqrt(11)
    new=2*math.sqrt(2)*lam*amin*rho
    old=4/39999
    remainder=math.pi/(2e5*(1-1e-8))
    one=new+old+remainder
    two=2*(new+old)+remainder
    assert one<1.19e-4 and two<2.21e-4
    assert 27*math.pi**2/25<11
    delta_bound=1.5+.5*amin**2*11
    Qbound=8*delta_bound+3
    lower=64*amin**2*11+Qbound
    assert Qbound<16 and lower<32
    theta=math.pi/4
    U=math.cos(theta)*np.eye(2)-1j*math.sin(theta)*np.array([[0,1],[1,0]])
    P=np.diag([0,1])
    noncomm=float(np.linalg.norm(U@P-P@U,2))
    assert abs(noncomm-1/math.sqrt(2))<1e-14
    return dict(kind='fixed_original_coupling_and_duration_budget',passed=True,
                lambda_upper=lam,a_upper=amin,density_L2_squared_upper=11,
                new_cloud_upper=new,old_cloud_upper=old,original_time_phase_upper=remainder,
                unperturbed_comparison_upper=one,one_insertion_history_upper=two,
                matching_Q_norm_upper=Qbound,energy_lower_shift_coefficient_upper=lower,
                noncommuting_rotation_norm=noncomm,original_T0_reselected=False)

def coherent(alpha,n):
    vecs=[]
    for z in alpha:
        v=np.empty(n,complex);v[0]=np.exp(-abs(z)**2/2)
        for k in range(1,n):v[k]=v[k-1]*z/math.sqrt(k)
        vecs.append(v)
    return np.kron(*vecs)

def history_group():
    # Independent finite two-oscillator identity diagnostic, NOT physical 1147 parameters.
    n=10
    freq=np.array([1.2,1.9])
    ga=np.array([.12,.08]);gb=np.array([.07,-.09])
    tau=.7;T=1.6
    theta=math.pi/4
    RA=math.cos(theta)*np.eye(2)-1j*math.sin(theta)*np.array([[0,1],[1,0]])
    R=np.kron(RA,np.eye(2))
    ann=np.diag(np.sqrt(np.arange(1,n)),1)
    a=[np.kron(ann,np.eye(n)),np.kron(np.eye(n),ann)]
    hfield=sum(w*(x.conj().T@x) for w,x in zip(freq,a))
    num=hfield.shape[0]
    energies=np.array([.021,.448,.443,.870])
    labels=[(0,0),(0,1),(1,0),(1,1)]
    matrices=[];couplings=[];decomposed=[]
    for i,(aa,bb) in enumerate(labels):
        g=aa*ga+bb*gb
        V=sum(t*(x+x.conj().T) for t,x in zip(g,a))
        h=energies[i]*np.eye(num)+hfield+V
        matrices.append(h);couplings.append(g);decomposed.append(V)
    diag=[np.linalg.eigh(h) for h in matrices]
    def evolve(state,t):
        out=np.empty_like(state)
        for i,(e,U) in enumerate(diag):
            out[i]=(U@(np.exp(-1j*t*e)[:,None]*(U.conj().T@state[i].T))).T
        return out
    rng=np.random.default_rng(1150)
    coefficients=rng.normal(size=(4,2))+1j*rng.normal(size=(4,2))
    coefficients/=np.linalg.norm(coefficients)
    initial=np.zeros((4,2,num),complex);initial[:,:,0]=coefficients
    before=evolve(initial,tau)
    after=np.einsum('ij,jrn->irn',R,before)
    final=evolve(after,T-tau)
    def evolve_coherent(alpha,g,t):
        z=g/freq
        a2=(alpha+z)*np.exp(-1j*freq*t)
        phase=float(np.sum(z*z*freq*t+np.imag(z*np.conj(alpha))-np.imag(z*np.conj(a2))))
        return a2-z,phase
    direct=np.zeros_like(final)
    for start in range(4):
        alpha1,p1=evolve_coherent(np.zeros(2,complex),couplings[start],tau)
        for finish in range(4):
            if abs(R[finish,start])==0:continue
            alpha2,p2=evolve_coherent(alpha1,couplings[finish],T-tau)
            ph=np.exp(1j*(p1+p2)-1j*(energies[start]*tau+energies[finish]*(T-tau)))
            direct[finish]+=R[finish,start]*ph*coefficients[start,:,None]*coherent(alpha2,n)
    err=float(np.linalg.norm(final-direct))
    assert err<1e-10
    def expect(state,ops):
        return float(sum(np.vdot(state[i],(ops[i]@state[i].T).T).real for i in range(4)))
    work=expect(after,matrices)-expect(before,matrices)
    free=[energies[i]*np.eye(num) for i in range(4)]
    work_free=expect(after,free)-expect(before,free)
    work_int=expect(after,decomposed)-expect(before,decomposed)
    work_field=expect(after,[hfield]*4)-expect(before,[hfield]*4)
    residual=abs(work-work_free-work_int)
    conserved=abs(expect(final,matrices)-expect(after,matrices))
    assert residual<2e-13 and abs(work_field)<2e-13 and conserved<2e-12
    return dict(kind='noncommuting_history_and_complete_kick_work_identity_diagnostic',
                passed=True,oscillator_cutoff=n,reference_dimension=2,unknown_input_seed=1150,
                total_time=T,pulse_time=tau,physical_continuum_parameters=False,
                field_reset_at_pulse=False,coherent_vs_spectral_full_state_residual=err,
                pulse_work=work,material_including_matching_work=work_free,
                interaction_storage_work=work_int,free_field_jump=work_field,
                complete_work_residual=residual,post_pulse_energy_conservation_residual=conserved)

def compute():
    manifest=json.loads((BASE/'source_hashes.json').read_text(encoding='utf-8'))
    for n,h in manifest.items():
        assert hashlib.sha256((BASE/n).read_bytes()).hexdigest()==h,n
    groups=[material_group(),budget_group(),history_group()]
    return dict(round=1150,all_passed=True,numerical_check_groups=3,groups=groups,
                frozen_source_count=len(manifest),scientific_calibration_increment=0,
                new_adopted_axioms=0,full_continuum_simulation=False,goal_complete=False)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--write-results',action='store_true')
    args=ap.parse_args();v=compute();s=json.dumps(v,ensure_ascii=False,indent=2)+'\n'
    f=BASE/'results.json'
    if args.write_results:
        if f.exists():assert f.read_text(encoding='utf-8')==s
        else:f.write_text(s,encoding='utf-8')
    else:assert json.loads(f.read_text(encoding='utf-8'))==v
    print(json.dumps(dict(status='PASS',groups=3,
                         history_error_upper=v['groups'][1]['one_insertion_history_upper'],
                         diagnostic_full_state_residual=v['groups'][2]['coherent_vs_spectral_full_state_residual'])))
if __name__=='__main__':main()
