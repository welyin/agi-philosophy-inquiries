"""606: field-dependent chiral fibres and the inherited electric quadratic form.
Uses a declared Euclidean Wilson/GW projector only as an interface diagnostic.
It is not a completed chiral Hamiltonian or a proof of operational FUCP rights.
"""
import argparse,hashlib,itertools,json
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_chiral_fibre_source_results.json'
S=[np.array([[0,1],[1,0]],complex),np.array([[0,-1j],[1j,0]]),np.diag([1,-1]).astype(complex)]
Z=np.zeros((2,2),complex);I2=np.eye(2)
GAMMA=[np.block([[Z,-1j*s],[1j*s,Z]]) for s in S]+[np.block([[Z,I2],[I2,Z]])]
G5=GAMMA[0]@GAMMA[1]@GAMMA[2]@GAMMA[3]
SITES=list(itertools.product(range(2),repeat=4))
def wilson(theta,charge):
    n=len(SITES);D=np.zeros((4*n,4*n),complex);dD=np.zeros_like(D)
    for mu in range(4):
        T=np.zeros((n,n),complex)
        for row,x in enumerate(SITES):
            y=list(x);y[mu]=(y[mu]+1)%2
            T[row,SITES.index(tuple(y))]=np.exp(1j*charge*theta) if mu==0 else 1
        dagger=T.conj().T
        D+=np.kron(np.eye(n)-.5*(T+dagger),np.eye(4))+.5*np.kron(T-dagger,GAMMA[mu])
        if mu==0:
            dT=1j*charge*T;ddagger=dT.conj().T
            dD+=-.5*np.kron(dT+ddagger,np.eye(4))+.5*np.kron(dT-ddagger,GAMMA[mu])
    g5=np.kron(np.eye(n),G5)
    return g5@(D-np.eye(4*n)),g5@dD,g5
def projector(theta,charge,positive=True):
    H,dH,g5=wilson(theta,charge)
    E,V=np.linalg.eigh(H);gap=float(min(abs(E)))
    assert gap>.8 and np.max(abs(H-H.conj().T))<1e-13
    signs=(E>0).astype(float)
    if not positive:signs=1-signs
    P=(V*signs)@V.conj().T
    delta=E[:,None]-E[None,:]
    numerator=signs[:,None]-signs[None,:]
    factor=np.divide(numerator,delta,out=np.zeros_like(delta),where=abs(delta)>1e-10)
    dP=V@(factor*(V.conj().T@dH@V))@V.conj().T
    Phi=P@dP@dP@P
    return P,dP,Phi,gap,H,dH,g5
def projector_check():
    rows=[]
    for q in (0,1,-2,-3,4,-6):
        P,dP,Phi,gap,H,dH,g5=projector(.001,q)
        eps=2e-6
        fd=(projector(.001+eps,q)[0]-projector(.001-eps,q)[0])/(2*eps)
        error=float(np.max(abs(fd-dP)))
        E,V=np.linalg.eigh(H);sign=(V*np.sign(E))@V.conj().T
        D=np.eye(64)+g5@sign
        GW=float(np.max(abs(g5@D+D@g5-D@g5@D)))
        assert max(error,GW)<3e-8
        assert np.max(abs(P@P-P))<1e-13
        assert np.max(abs(P@dP@P))<1e-13
        assert np.linalg.eigvalsh(Phi).min()>-1e-12
        bound=float(np.linalg.norm(dH,2)/(2*gap))
        assert np.linalg.norm(dP,2)<=bound+1e-12
        rows.append(dict(integer_charge=q,rank=int(round(np.trace(P).real)),gap=gap,
            finite_difference_error=error,GW_residual=GW,
            geometric_trace=float(np.trace(Phi).real),
            derivative_norm=float(np.linalg.norm(dP,2)),derivative_bound=bound))
    return dict(dimension=64,sites=16,input_Euclidean_dimensions=4,m0_input=1,
        flat_plaquettes_exactly_identity=True,rows=rows,
        temporal_Hamiltonian_identification_not_assumed=True)
def combined_source_check():
    modules={'Q':(6,1,True),'u':(3,4,False),'d':(3,-2,False),
             'L':(2,-3,True),'e':(1,-6,False),'nu':(1,0,False)}
    one=float(np.trace(projector(0,1)[2]).real)
    rows=[];total=0.;charge_index=0;cubic=0
    for name,(multiplicity,q,left) in modules.items():
        value=float(np.trace(projector(0,q,left)[2]).real)
        assert abs(value-q*q*one)<1e-10
        total+=multiplicity*value;charge_index+=multiplicity*q*q
        cubic+=(1 if left else -1)*multiplicity*q**3
        rows.append(dict(module=name,multiplicity=multiplicity,charge=q,
                         geometric_trace=value))
    assert charge_index==120 and cubic==0 and abs(total-120*one)<1e-9
    theta=.12;kappa=np.exp(-2*theta);eps=1e-6
    source=-2*kappa*total
    fd=(np.exp(-2*(theta+eps))*total-np.exp(-2*(theta-eps))*total)/(2*eps)
    assert abs(fd-source)<1e-6
    return dict(rows=rows,single_unit_charge_trace=one,total_geometric_trace=total,
        quadratic_charge_index=charge_index,signed_cubic_anomaly=0,
        declared_electric_scale=kappa,filled_fibre_energy=kappa*total,
        scale_source=source,source_derivative_residual=abs(fd-source),
        source_is_not_cancelled_by_anomaly=True,
        filled_state_is_diagnostic_not_selected_vacuum=True)
def annihilators(n):
    out=[]
    for j in range(n):
        a=np.zeros((2**n,2**n),complex)
        for state in range(2**n):
            if state>>j&1:
                a[state^(1<<j),state]=(-1)**((state&((1<<j)-1)).bit_count())
        out.append(a)
    return out
def exterior_form_check():
    P,dP,_,_,_,_,_=projector(0,1)
    val,V=np.linalg.eigh(P);inside=V[:,val>.5];outside=V[:,val<.5]
    sv=np.linalg.svd(outside.conj().T@dP@inside,compute_uv=False)[:2]
    B=np.diag(sv);zero=np.zeros((2,2));K=np.block([[zero,-B],[B,zero]])
    c=annihilators(4);KF=sum(K[i,j]*c[i].conj().T@c[j] for i in range(4) for j in range(4))
    embed=np.eye(16)[:,:4]
    geometric=(KF@embed).conj().T@(KF@embed)
    logical=annihilators(2)
    expected=sum((B.conj().T@B)[i,j]*logical[i].conj().T@logical[j] for i in range(2) for j in range(2))
    error=float(np.max(abs(geometric-expected)));assert error<1e-13
    # Local frame change: connection and horizontal contribution are distinct.
    A=np.array([[.2j,.3+.1j],[-.3+.1j,-.4j]])
    V0=np.eye(4,dtype=complex)[:,:2];dV=np.vstack([A,B])
    f=np.array([.4+.2j,-.3+.7j]);df=np.array([-.1+.4j,.2-.1j])
    direct=float(np.linalg.norm(V0@df+dV@f)**2)
    split=float(np.linalg.norm(df+A@f)**2+np.linalg.norm(B@f)**2)
    assert abs(direct-split)<1e-13
    return dict(actual_projector_offdiagonal_singular_values=sv.tolist(),
        four_ambient_mode_Fock_check=True,Fock_geometric_error=error,
        logical_geometric_eigenvalues=np.linalg.eigvalsh(expected).tolist(),
        local_frame_form_residual=abs(direct-split),
        exterior_identity_for_general_rank_proved_in_note=True)
def run():
    data=dict(projector=projector_check(),common_source=combined_source_check(),exterior=exterior_form_check())
    deps=('research_note_598.md','research_note_603.md','research_note_605.md',
          'cognitive_foundation_bridge_605.md','cognitive_foundation_bridge_605_results.json',
          'round606_drafts/foundation_inheritance_entry.md')
    return dict(round=606,tests_run=3,failures=0,errors=0,**data,
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope=dict(electric_connection_and_positive_term_derived_conditionally=True,
            actual_GW_projector_diagnostic=True,original_charge_multiplicities_retained=True,
            no_chiral_Hamiltonian_or_measure_completion=True,
            no_FUCP_operational_rights_inferred_from_complex_fibres=True,
            no_lattice_continuum_or_GR_claim=True))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');a=p.parse_args();r=run()
    if a.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==r
    print(json.dumps(dict(round=606,tests=3,all_passed=True,source=r['common_source'])))
