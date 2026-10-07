"""927: same finite autonomous process, instantaneous vs future encoding contracts.
The reference shift is reused; no continuum, control synthesis, or parameter optimisation.
"""
from pathlib import Path
import argparse,json,hashlib
import numpy as np
HERE=Path(__file__).resolve().parent
I=np.eye(2,dtype=complex)
X=np.array([[0,1],[1,0]],complex);Y=np.array([[0,-1j],[1j,0]],complex);Z=np.diag([1.,-1.]).astype(complex)
S=(X,Y,Z)
def op(a):return float(np.linalg.norm(a,2))
def com(a,b):return a@b-b@a
def kron(*xs):
    a=np.array([[1]],complex)
    for x in xs:a=np.kron(a,x)
    return a
def expH(h,t):
    e,u=np.linalg.eigh(h);return (u*np.exp(-1j*t*e))@u.conj().T
def expectation(psi,a):return float(np.vdot(psi,a@psi).real)
def basis(d,j):
    a=np.zeros(d,complex);a[j]=1;return a

def build(gvalues):
    N=5;L=4;r=2;a=.7;kappa=4.;b=1.;e0=2.;d=2*2*r*N
    def ix(s,z,q,n):return ((s*2+z)*r+q)*N+n
    def cx(s,z,q,n):return ((s*2+z)*r+q)*L+n
    E=np.zeros((d,2*2*r*L),complex)
    h0=np.zeros(d);end=np.zeros(d);ref=np.zeros(d);qdiag=np.zeros(d);pb=np.zeros(d)
    J=np.zeros((d,d),complex);JG=np.zeros((d,d),complex)
    for s in range(2):
        for z in range(2):
            for q in range(r):
                for n in range(N):
                    i=ix(s,z,q,n);end[i]=e0+kappa*q-s*a*q;ref[i]=b+a*n
                    h0[i]=end[i]+ref[i];qdiag[i]=q;pb[i]=s
                for n in range(L):E[ix(s,z,q,n+s*q),cx(s,z,q,n)]=1
    for z in range(2):
        for q in range(r):
            for n in range(N-q):
                i=ix(1,z,q,n+q);j=ix(0,z,q,n)
                J[i,j]=1;JG[i,j]=gvalues[q]
    alpha=[kron(Z,s,np.eye(r),np.eye(N)) for s in S]
    h0=np.diag(h0);he=np.diag(end);hr=np.diag(ref);Q=np.diag(qdiag);PB=np.diag(pb)
    hint=JG+JG.conj().T;H=h0+hint
    ca=[kron(Z,s,np.eye(r),np.eye(L)) for s in S]
    resource=kron(np.eye(4),np.diag(np.arange(r)*kappa),np.eye(L))+kron(np.eye(4*r),np.diag(b+e0+a*np.arange(L)))
    mix=kron(X,I,np.diag(gvalues),np.eye(L));hc=resource+mix
    return dict(N=N,L=L,r=r,a=a,kappa=kappa,E=E,alpha=alpha,h0=h0,he=he,hr=hr,Q=Q,PB=PB,H=H,hc=hc,ca=ca,mix=mix,resource=resource,ix=ix,cx=cx)

def run():
    T=np.pi/2;t=np.pi/8
    results=[]
    for name,gs in [('uniform',[1.,1.]),('distinct_same_arrival',[1.,5.])]:
        m=build(gs);E=m['E'];H=m['H'];dc=E.shape[1];d=len(H);L=m['L'];r=m['r']
        P=E@E.conj().T
        assert op(E.conj().T@E-np.eye(dc))<1e-14
        code_error=op(H@E-E@m['hc'])
        assert code_error<1e-13
        assert op(com(H,P))<1e-13
        timeU=expH(H,T)
        half=dc//2;VA=E[:,:half];VB=E[:,half:]
        rhs=-1j*VB@expH(m['resource'][:half,:half],T)
        transfer_error=op(timeU@VA-rhs)
        assert transfer_error<5e-13
        parity=kron(X,I,np.eye(r),np.eye(L))
        parity_error=op(com(m['hc'],parity))
        rotation_error=max(op(com(m['hc'],kron(I,s/2,np.eye(r),np.eye(L)))) for s in S)
        current=iPB=1j*com(H,m['PB'])
        ref_flow_error=op(1j*com(H,m['hr'])-m['a']*m['Q']@current)
        endpoint_flow_error=op(1j*com(H,m['he'])+m['a']*m['Q']@current)
        assert max(ref_flow_error,endpoint_flow_error)<1e-13
        conservation_error=op(com(m['h0'],H))
        assert conservation_error<1e-13
        instant=0.;future=0.;derivative=0.;encoding_energies=[]
        momentum=np.array([.13,-.17,.19])
        hp=H+sum(momentum[i]*m['alpha'][i] for i in range(3))
        hc_p=m['hc']+sum(momentum[i]*m['ca'][i] for i in range(3))
        momentum_code_error=op(hp@E-E@hc_p)
        assert momentum_code_error<1e-13
        for u in (X,Z):
            uc=kron(np.eye(4),u,np.eye(L));up=E@uc@E.conj().T+np.eye(d)-P
            assert op(up.conj().T@up-np.eye(d))<1e-14
            for ac,ap in zip(m['ca'],m['alpha']):
                instant=max(instant,op(com(ap,up)))
                derivative=max(derivative,op(com(1j*com(m['hc'],ac),uc)))
                v=expH(hc_p,.23);future=max(future,op(com(v.conj().T@ac@v,uc)))
            encoding_energies.append(op(com(up,m['h0'])))
        psi=(basis(d,m['ix'](0,0,0,0))+basis(d,m['ix'](0,0,1,0)))/np.sqrt(2)
        out=timeU@psi
        flow={k:expectation(out,m[k])-expectation(psi,m[k]) for k in ('hr','he','H')}
        assert abs(flow['hr']-.35)<1e-12 and abs(flow['he']+.35)<1e-12 and abs(flow['H'])<1e-12
        shortU=expH(H,t);probs=[]
        effect=(np.eye(d)+m['alpha'][2])/2
        for q in range(r):
            psi_q=basis(d,m['ix'](0,0,q,0));evolved=shortU@psi_q
            probs.append(expectation(evolved,effect))
        difference=abs(probs[0]-probs[1])
        assert instant<1e-13 and parity_error<1e-13 and rotation_error<1e-13
        if name=='uniform':assert derivative<1e-13 and future<1e-12 and difference<1e-12
        else:assert abs(derivative-8)<1e-12 and abs(difference-1/np.sqrt(2))<1e-12 and future>.1
        results.append(dict(name=name,mixing_rates=gs,physical_dimension=d,code_dimension=dc,
          positive_total_H_minimum=float(np.min(np.linalg.eigvalsh(H))),code_intertwining_error=code_error,
          momentum_code_intertwining_error=momentum_code_error,all_unknown_input_transfer_error=transfer_error,
          bare_energy_conservation_error=conservation_error,resource_current_identity_error=ref_flow_error,
          endpoint_current_identity_error=endpoint_flow_error,role_exchange_error_on_code=parity_error,
          rotation_error_on_code=rotation_error,instant_encoding_commutator_norm=instant,
          first_future_derivative_encoding_defect=derivative,finite_future_encoding_commutator_norm=future,
          encoding_bare_energy_commutators=encoding_energies,mean_energy_changes=flow,
          direction_probabilities_at_pi_over_8=probs,finite_probability_difference=difference))
    # Calibration of the analytic Schur result for the strengthened future contract.
    ca=[kron(Z,s,I) for s in S]
    herm=[I,X,Y,Z];cols=[]
    for g in herm:
        h=kron(X,I,g);col=[]
        for u in (X,Z):
            for a in ca:col.extend(com(1j*com(h,a),kron(I,I,u)).reshape(-1))
        col=np.array(col);cols.append(np.r_[col.real,col.imag])
    rank=int(np.count_nonzero(np.linalg.svd(np.column_stack(cols),compute_uv=False)>1e-10))
    assert 4-rank==1
    return dict(round=927,date='2026-10-06',all_scientific_checks_passed=True,arrival_time=T,
      comparison_time=t,joint_candidate_family=results,strong_future_contract_mixing_dimension=4-rank,
      one_process_carries_records_direction_and_resource=True,unknown_record_kept_in_joint_reference_code=True,
      reference_preparation_common_to_unknown_inputs=True,resource_shift_not_claimed_new=True,
      uniform_coupling_not_assumed_as_cognitive_requirement=True,
      strong_future_invariance_is_extra_input=True,active_encodings_not_claimed_costless=True,
      autonomous_reencoding_controller_constructed=False,spacetime_dimension_generated=False,
      gauge_group_or_GR_generated=False,all_physics_inputs_selected=False,whole_stage_completed=False,full_goal_completed=False,
      source_hashes={str(Path('927')/Path(__file__).name):hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    out=run();p=HERE/'joint_encoding_resource_selection_results.json'
    if args.write:
        assert not p.exists();p.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in out.items() if k!='source_hashes'},ensure_ascii=False,indent=2))
