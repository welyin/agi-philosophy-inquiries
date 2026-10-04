"""751: original record erases a source coherence required by internal balance."""
import sys,json
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import record_source_sufficiency as base

def run():
    old=base.old
    bg=old.sources.background(0.)
    B,_=old.sources.matrix_and_source(*bg,*[np.zeros_like(x) for x in bg])
    labels,Js=old.full_sources(B)
    e=np.zeros(64);e[30]=e[62]=1/np.sqrt(2)
    f=np.zeros(64);f[31]=1.
    V=np.column_stack((e,f))
    W=np.block([[V,np.zeros_like(V)],[np.zeros_like(V),V.conj()]])
    cs=[base.destroy(2,j) for j in range(2)]
    eta=cs+[c.conj().T for c in cs]
    n=cs[0].conj().T@cs[0];ls=[np.eye(4)-n,n];code=np.eye(4)[:,[0,3]]
    operators=[]
    for J in Js:
        j=W.conj().T@J@W
        O=sum(.5*j[a,b]*eta[a].conj().T@eta[b] for a in range(4) for b in range(4))
        operators.append(O)
    H=operators[0];D=sum(l@H@l for l in ls)-H
    gap_ev,U=np.linalg.eigh(code.conj().T@D@code)
    states=[code@U[:,0],code@U[:,-1]]
    rhos=[np.outer(v,v.conj()) for v in states]
    cq=[[l@rho@l for l in ls] for rho in rhos]
    cqerror=max(float(np.max(abs(cq[0][r]-cq[1][r]))) for r in range(2))
    differences=[];full_errors=[];costs=[]
    def covariance(rho):
        return np.array([[np.trace(rho@eta[b].conj().T@eta[a]) for b in range(4)] for a in range(4)])
    for rho,out in zip(rhos,cq):
        delta=W@(covariance(sum(out))-covariance(rho))@W.conj().T
        cc=[]
        for J,O in zip(Js,operators):
            a=float(np.trace((sum(out)-rho)@O).real);b=float(np.trace(J@delta).real/2)
            full_errors.append(abs(a-b));cc.append(a)
        costs.append(cc)
    before=[float(np.trace(rho@H).real) for rho in rhos]
    after=[float(np.trace(sum(out)@H).real) for out in cq]
    # A minimal exact dilation puts original r in an orthogonal physical register.
    Wiso=np.zeros((8,4),complex)
    for r,l in enumerate(ls):
        Wiso+=np.kron(l,np.eye(2)[:,r:r+1])
    def trace_env(rho):return np.einsum('aibi->ab',rho.reshape(4,2,4,2))
    def trace_sys(rho):return np.einsum('aiaj->ij',rho.reshape(4,2,4,2))
    outstates=[Wiso@v for v in states]
    global_rhos=[np.outer(v,v.conj()) for v in outstates]
    sm=[trace_env(r) for r in global_rhos];am=[trace_sys(r) for r in global_rhos]
    marginal_error=max(float(np.max(abs(sm[0]-sm[1]))),float(np.max(abs(am[0]-am[1]))))
    arbitrary_HA=np.array([[.2,.3+.1j],[.3-.1j,1.4]])
    additive=np.kron(H,np.eye(2))+np.kron(np.eye(4),arbitrary_HA)
    totals=[float(np.trace(r@additive).real) for r in global_rhos]
    overlap=float(abs(np.vdot(outstates[0],outstates[1])))
    # Same marginals do not mean the joint state is the same; global correlation keeps the phase.
    joint_trace_distance=float(np.linalg.norm(global_rhos[0]-global_rhos[1],ord='nuc')/2)
    phi=bg[1];F=old.native.old.matter.original.F(phi);Y=old.native.old.matter.Y['s']
    analytic_kappa=abs(Y*phi[4])/np.sqrt(2*F)
    kappa_error=abs(analytic_kappa-abs((code.conj().T@H@code)[0,1]))
    assert kappa_error<1e-13
    assert cqerror<1e-12 and marginal_error<1e-12 and abs(totals[1]-totals[0])<1e-12
    assert max(full_errors)<1e-12 and abs(before[1]-before[0])>.01
    assert overlap<1e-12 and abs(joint_trace_distance-1)<1e-12
    return dict(round=751,original_sources=labels,record_probabilities=[float(np.trace(r).real) for r in cq[0]],
        identical_full_cq_output_error=cqerror,
        energy_before=before,energy_after=after,energy_injections=[x[0] for x in costs],
        all24_source_injections=costs,all24_full_covariance_agreement=max(full_errors),
        record_only_compensation_minimax_error=float((gap_ev[-1]-gap_ev[0])/2),
        original_energy_offdiagonal_abs=float(abs((code.conj().T@H@code)[0,1])),
        analytic_majorana_kappa=float(analytic_kappa),analytic_kappa_error=float(kappa_error),
        minimal_dilation_equal_marginals_error=marginal_error,
        minimal_dilation_joint_trace_distance=joint_trace_distance,
        minimal_dilation_additive_final_energies=totals,
        no_claim_selected_minimal_dilation_energy_conserving=True,
        general_exact_dilation_obstruction_analytic=True,
        scope='Independent fixed apparatus; exact original Luders operation; conserved additive energy and zero endpoint interaction; finite expectation domain. No rejection of approximate records or interacting endpoints.',
        diagnostic_not_full_GR_state_or_continuum_solution=True)

if __name__=='__main__':
    r=run()
    with (HERE/'erased_source_balance_results.json').open('x',encoding='utf8') as f:
        json.dump(r,f,ensure_ascii=False,indent=2)
    print(json.dumps({k:r[k] for k in ('record_probabilities','identical_full_cq_output_error','energy_before','energy_after','energy_injections','record_only_compensation_minimax_error','minimal_dilation_joint_trace_distance')},ensure_ascii=False))
