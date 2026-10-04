"""733: finite record embedding does not remove geometric reference response.

Complete original local matrices calibrate a continuum principal-symbol
obstruction. No bare-cutoff quantity is called a renormalized stress prediction.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
import joint_dynamic_continuum_reference as reference
import joint_relative_source_development as prior

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'round732_drafts'))
import source_feedback_entry as family
TARGET=HERE/'joint_reference_polarization_boundary_results.json'


def maxabs(x):
    return float(np.max(abs(x)))


def matrix(g,k):
    p=reference.point();ev,V=np.linalg.eigh(g);frame=(V/np.sqrt(ev))@V.T
    gen=reference.gauge_h(p['a'],p['a0'])
    kinetic=sum(reference.GAMMA[j]*sum(frame[j,i]*k[i] for i in range(3)) for j in range(3))
    gauge=sum(reference.GAMMA[j][:32,:32]@sum(frame[j,i]*gen[i] for i in range(3)) for j in range(3))
    return (kinetic+reference.old.bdg(gauge,np.zeros_like(gauge))
            +reference.old.bdg(*reference.old.matter.mass_matrices(p['phi'])))


def geometries():
    base=reference.point()['g'];gamma=.12
    # Original initial spatial metric is conformally Euclidean at this point.
    assert maxabs(base-np.trace(base)/3*np.eye(3))<1e-13
    changed=np.diag(np.exp([2*gamma,-2*gamma,0.]))@base
    n=np.array([1.,1.,0.])/np.sqrt(2)
    direction=np.exp([-gamma,gamma,0.])*n;direction/=np.linalg.norm(direction)
    symbol0=(np.eye(64)-sum(n[i]*reference.GAMMA[i] for i in range(3)))/2
    symbol1=(np.eye(64)-sum(direction[i]*reference.GAMMA[i] for i in range(3)))/2
    cosine=np.cosh(gamma)/np.sqrt(np.cosh(2*gamma))
    exact=np.sqrt((1-cosine)/2)
    return base,changed,n,symbol0,symbol1,float(exact)


def principal_reference_check():
    base,changed,n,s0,s1,exact=geometries()
    eigen=np.linalg.eigvalsh(s1-s0)
    assert maxabs(abs(eigen)-exact)<3e-15 and exact>.01
    rows=[]
    for scale in (10.,100.,1000.):
        P0=reference.projector(matrix(base,scale*n));P1=reference.projector(matrix(changed,scale*n))
        norm=float(np.linalg.norm(P1-P0,2))
        error=float(np.linalg.norm((P1-P0)-(s1-s0),2))
        rows.append(dict(momentum_scale=scale,projector_difference=norm,principal_symbol_error=error,
                         min_singular_difference=float(np.linalg.svd(P1-P0,compute_uv=False).min())))
    assert rows[-1]['principal_symbol_error']<rows[0]['principal_symbol_error']/80
    assert abs(rows[-1]['projector_difference']-exact)<2e-3
    # A conformal scale has the same leading normalized direction; no false
    # claim that every metric change gives a nonzero order-zero symbol.
    conformal_direction=.73*n;conformal_direction/=np.linalg.norm(conformal_direction)
    conformal_error=float(np.linalg.norm(conformal_direction-n))
    assert conformal_error<1e-14
    return dict(gamma=.12,exact_asymptotic_operator_norm=exact,
                all64_symbol_singular_values_error=maxabs(abs(eigen)-exact),rows=rows,
                conformal_leading_direction_error=conformal_error,
                local_symbol_calibration_not_a_continuum_covariance_computation=True)


def finite_patch_boundary_check():
    base,changed,n,_,_,exact=geometries();Bs=[];Cs=[]
    for scale in (800.,1000.,1200.):
        Bs.append(prior.assemble([matrix(base,scale*n),matrix(base,-scale*n)]))
        Cs.append(prior.assemble([matrix(changed,scale*n),matrix(changed,-scale*n)]))
    B0=reference.old.block(Bs);B1=reference.old.block(Cs)
    P0=reference.projector(B0);Pstar=reference.projector(B1)
    Cpair=np.block([[np.zeros((64,64)),np.eye(64)],[np.eye(64),np.zeros((64,64))]])
    conjugation=reference.old.block([Cpair]*3)
    f=np.zeros(384,complex)
    for i,phase in enumerate((1.,1j,np.exp(.4j))):f[i*128+30]=phase/np.sqrt(3)
    cf=conjugation@f.conj();Q=np.outer(f,f.conj())+np.outer(cf,cf.conj());R=np.eye(384)-2*Q
    spanning=np.column_stack((f,cf,P0@f,P0@cf));U,s,_=np.linalg.svd(spanning,full_matrices=False)
    basis=U[:,s>1e-12];K=basis@basis.conj().T;Z=np.eye(384)-K
    patched=K@P0@K+Z@Pstar@Z
    finite=patched-Pstar;remaining=patched-P0
    rank=int(np.count_nonzero(abs(np.linalg.eigvalsh(finite))>1e-10));assert rank<=8
    residual_singular=np.linalg.svd(remaining,compute_uv=False)
    persistent=int(np.count_nonzero(residual_singular>exact/2))
    source_difference_error=maxabs((R@patched@R-patched)/2-(R@P0@R-P0)/2)
    reality=maxabs(conjugation@patched.conj()@conjugation+patched-np.eye(384))
    assert source_difference_error<2e-13 and reality<3e-13
    assert persistent>=384-rank and residual_singular.max()>exact*.9
    eig=np.linalg.eigvalsh(patched);assert eig.min()>-3e-13 and eig.max()<1+3e-13
    return dict(Nambu_dimension=384,record_block_dimension=len(basis.T),patch_rank=rank,
        unchanged_record_difference_error=source_difference_error,reality_error=reality,
        residual_operator_norm=float(residual_singular.max()),
        singular_values_above_half_asymptotic_bound=persistent,
        rank_bound_required_count=384-rank,
        discarded_reference_change_cannot_be_finite_record_patch=True)


def second_order_reference_response_check():
    P0,_,_,_,C=prior.initial_modes();e=np.zeros(128,complex);e[30]=e[62]=1/np.sqrt(2)
    ce=C@e.conj();Q=np.outer(e,e.conj())+np.outer(ce,ce.conj());R=np.eye(128)-2*Q
    U,s,_=np.linalg.svd(np.column_stack((e,ce,P0@e,P0@ce)),full_matrices=False)
    basis=U[:,s>1e-12];K=basis@basis.conj().T;Z=np.eye(128)-K
    delta0=(R@P0@R-P0)/2
    P,dPs,B,Ds=family.flow(0.,0.,derivatives=True)
    assert maxabs(P-P0)<1e-12
    G0=Ds[0];force0=prior.source(P0,G0)
    # d_gamma G=-G at the endpoint. Reference patch retains its state response.
    response=prior.source(P0,-G0)+prior.source(Z@dPs[0]@Z,G0)
    rows=[]
    for eps in (.04,.02,.01):
        P,_,B,Ds=family.flow(eps,0.,derivatives=False);G=Ds[0]
        patched=K@P0@K+Z@P@Z
        fref=prior.source(patched,G);fplus=prior.source(patched+delta0,G)
        frecord=prior.source(delta0,G)
        omitted=eps*(fref-force0)
        identity=(fplus-force0)-frecord-(fref-force0)
        assert abs(identity)<1e-12 and abs(omitted)>1e-5
        rows.append(dict(epsilon=eps,reference_source=fref,record_source=frecord,
            omitted_weighted_reference_term=omitted,term_over_epsilon_squared=omitted/eps**2,
            exact_branch_decomposition_error=abs(identity),
            reference_patch_source=prior.source(patched-P,G)))
    error=abs(rows[-1]['term_over_epsilon_squared']-response)
    assert error<abs(rows[0]['term_over_epsilon_squared']-response)/3
    h=2e-5;values=[]
    for gamma in (h,-h):
        P,_,_,Ds=family.flow(gamma,0.,derivatives=False)
        values.append(prior.source(K@P0@K+Z@P@Z,Ds[0]))
    fd=(values[0]-values[1])/(2*h)
    assert abs(fd-response)<2e-8
    return dict(base_reference_geometric_source=force0,
        exact_reference_response_coefficient=response,finite_difference=fd,
        derivative_error=abs(fd-response),rows=rows,
        finite_matrix_sources_not_renormalized_continuum_predictions=True)


def run():
    checks=('principal_reference_check','finite_patch_boundary_check','second_order_reference_response_check')
    results={name:globals()[name]() for name in checks}
    deps=('research_note_730.md','research_note_731.md','research_note_732.md',
          'joint_relative_source_development.py','joint_relative_source_development_results.json',
          'round732_drafts/source_feedback_entry.py','round732_drafts/source_feedback_entry_results.json',
          'round733_drafts/reference_embedding_entry.py','round733_drafts/reference_embedding_entry_results.json')
    return dict(round=733,tests_run=3,failures=0,errors=0,checks=list(checks),results=results,
        dependency_hashes={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in deps},
        scope='Original record block admits Hadamard embedding, but a nonconformal metric variation changes the Cauchy covariance principal symbol in the fixed operational identification. Finite smooth record patches cannot absorb that reference change. Nonlinear branch comparison contains background-reference polarization, already order epsilon squared in a perturbative response. This closes only the finite-record-only absolute feedback shortcut, not the full unified programme, relative models, or renormalized quantum gravity.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');a=p.parse_args();r=run()
    if a.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps({k:r[k] for k in ('round','tests_run','failures','errors')}))
