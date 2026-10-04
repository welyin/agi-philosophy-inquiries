"""751 entry: source information retained inside a fixed original record.
Finite original128 local-symbol audit; no continuum/GR implementation claim.
"""
import sys, json
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent; ROOT=HERE.parent
sys.path.insert(0,str(ROOT)); sys.path.insert(0,str(ROOT/'round750_drafts'))
import marked_source_entry as old

def destroy(n,j):
    out=np.zeros((2**n,2**n),complex)
    for state in range(2**n):
        if (state>>j)&1:
            out[state^(1<<j),state]=(-1)**((state&((1<<j)-1)).bit_count())
    return out

def run():
    bg=old.sources.background(0.)
    B,_=old.sources.matrix_and_source(*bg,*[np.zeros_like(x) for x in bg])
    labels,Js=old.full_sources(B)
    # Original event mode e; all modes are sterile/gauge-singlet.
    e=np.zeros(64);e[30]=e[62]=1/np.sqrt(2)
    d=np.zeros(64);d[30]=1/np.sqrt(2);d[62]=-1/np.sqrt(2)
    f=np.zeros(64);f[31]=1.
    V=np.column_stack((e,d,f))
    W=np.block([[V,np.zeros_like(V)],[np.zeros_like(V),V.conj()]])
    cs=[destroy(3,j) for j in range(3)]
    eta=cs+[c.conj().T for c in cs]
    N=cs[0].conj().T@cs[0]
    code=np.eye(8)[:,[0,6]] # |vac>, d^dagger f^dagger |vac>, even and n_e=0
    S=np.array([[[0,1],[1,0]],[[0,-1j],[1j,0]],[[1,0],[0,-1]]],complex)
    small=[]; opfull=[]; vecs=[]; source_errors=[]; parity_errors=[]
    full_parity=np.diag([(-1)**i.bit_count() for i in range(8)])
    for J in Js:
        j=W.conj().T@J@W
        O=sum(.5*j[a,b]*eta[a].conj().T@eta[b] for a in range(6) for b in range(6))
        A=code.conj().T@O@code;small.append(A);opfull.append(O)
        vecs.append([float(np.trace(A@s).real/2) for s in S])
        parity_errors.append(np.max(abs(O@full_parity-full_parity@O)))
    vecs=np.array(vecs)
    A=small[0]; ev,U=np.linalg.eigh(A)
    inputstates=[code@U[:,0],code@U[:,-1]]
    scalar=labels.index('scalar_4')
    rho_delta=np.outer(inputstates[1],inputstates[1].conj())-np.outer(inputstates[0],inputstates[0].conj())
    # In this convention P_ab=<eta_b^dag eta_a>, hence <O>=Tr(J P)/2.
    def covariance(psi):
        return np.array([[np.vdot(psi,eta[b].conj().T@eta[a]@psi) for b in range(6)] for a in range(6)])
    Ps=[covariance(v) for v in inputstates];delta=W@(Ps[1]-Ps[0])@W.conj().T
    differences=[]
    for J,O in zip(Js,opfull):
        by_fock=float(np.trace(rho_delta@O).real);by_full=float(np.trace(J@delta).real/2)
        source_errors.append(abs(by_fock-by_full));differences.append(by_full)
    record=[]
    for psi in inputstates:
        rho=np.outer(psi,psi.conj());branch=N@rho@N;post=(np.eye(8)-N)@rho@(np.eye(8)-N)+branch
        record.append(dict(probability_one=float(np.trace(branch).real),poststate_error=float(np.max(abs(post-rho)))))
    # Gauge generators annihilate the three sterile modes in original species dictionary.
    gauge_error=0.
    for _,_,H in old.sources.physical_generators():
        h=np.block([[H,np.zeros_like(H)],[np.zeros_like(H),H]])
        gauge_error=max(gauge_error,float(np.max(abs(h@V))))
    source_gap=float(ev[-1]-ev[0]);bracket=A@small[scalar]-small[scalar]@A
    bracket_norm=float(np.linalg.norm(bracket,2))
    errors=dict(mode_orthogonality=float(np.max(abs(V.conj().T@V-np.eye(3)))),
        code_occupation=float(np.max(abs(N@code))),code_even_parity=float(np.max(abs(full_parity@code-code))),
        gauge_singlet=gauge_error,all24_full_covariance_source_agreement=max(source_errors),
        even_source=max(map(float,parity_errors)))
    assert max(errors.values())<1e-12 and source_gap>.01 and bracket_norm>.001
    assert all(x['probability_one']<1e-14 and x['poststate_error']<1e-14 for x in record)
    return dict(round=751,original_Nambu_dimension=128,source_labels=labels,
        encoded_source_Pauli_vectors=vecs.tolist(),encoded_energy_eigenvalues=ev.tolist(),
        energy_gap=source_gap,record_only_energy_minimax_error=source_gap/2,
        singlet_source_minimax_error=float(np.linalg.norm(vecs[scalar])),
        energy_singlet_source_commutator_norm=bracket_norm,
        encoded_energy_singlet_Pauli_rank=int(np.linalg.matrix_rank(vecs[[0,scalar]],tol=1e-12)),
        identical_record_inputs=record,all24_source_differences=differences,errors=errors,
        arbitrary_input_audit_not_original_fixed_reference=True,
        encoded_sector_not_claimed_invariant_under_full_H=True,
        original_gauge_singlet_and_even_parity=True,
        not_full_quantum_Einstein_constraint_states=True,
        source_calibration_not_continuum_geometry_solution=True)

if __name__=='__main__':
    r=run()
    with (HERE/'record_source_sufficiency_results.json').open('x',encoding='utf8') as f:
        json.dump(r,f,ensure_ascii=False,indent=2)
    print(json.dumps({k:r[k] for k in ('energy_gap','record_only_energy_minimax_error','singlet_source_minimax_error','energy_singlet_source_commutator_norm','errors')},ensure_ascii=False))

