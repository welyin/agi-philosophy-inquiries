"""706: local boundary-preserving regularizers in original real histories.
Full Gauss/region theorem is analytic; numerics are an explicitly frozen-link,
aligned two-node radial/neutral-Majorana diagnostic, not a full gauge spectrum.
Original256 dimensional process is kept during every waiting interval.
"""
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import numpy as np
import joint_preparation_history_limit as jets
import joint_gibbs_preserving_transfer as thermal
HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_local_history_source_limit_results.json'
KEYS=jets.KEYS;HBAR=.7;BETA=1.2;TIMES=(.12,.19)


def data():
    original=thermal.original;n=2;w=.8;eps=w**(1/3)
    hs=np.linspace(.15,1.55,n+2)[1:-1];ss=np.linspace(-1.1,1.1,n+2)[1:-1]
    h,s=np.meshgrid(hs,ss,indexing='ij');h=h.ravel();s=s.ravel();r=len(h)
    F=original.M-(h*h+s*s)/6;measure=np.sqrt(original.M)*h**3/F**3
    def derivative(step):return (np.diag(np.ones(n-1),1)-np.diag(np.ones(n-1),-1))/(2*step)
    ri=1/np.sqrt(measure)
    dh=np.kron(derivative(hs[1]-hs[0]),np.eye(n))*ri[None,:]
    ds=np.kron(np.eye(n),derivative(ss[1]-ss[0]))*ri[None,:]
    ghh=F*(1-h*h/(6*original.M));gss=F*(1-s*s/(6*original.M));ghs=-F*h*s/(6*original.M)
    T=HBAR**2/(2*w)*(dh.T@((measure*ghh)[:,None]*dh)+ds.T@((measure*gss)[:,None]*ds)
       +dh.T@((measure*ghs)[:,None]*ds)+ds.T@((measure*ghs)[:,None]*dh))
    phi=np.zeros((r,5));phi[:,1]=h;phi[:,4]=s
    W=np.diag(w*original.node_potential(phi));e,v=np.linalg.eigh(T+W)
    localT=np.kron(T,np.eye(4));localW=np.kron(W,np.eye(4))
    B=np.zeros_like(localT,dtype=complex)
    for j,p in enumerate(phi):
        _,d=thermal.inherited.old.mass_matrices(p)
        B[4*j:4*j+4,4*j:4*j+4]=thermal.car.fock(np.zeros((2,2)),d[np.ix_([30,31],[30,31])])
    nloc=4*r;I=np.eye(nloc)
    parity=np.tile(np.array([1.,-1.,-1.,1.]),r);P=np.diag(parity)
    annih=[]
    for mode in range(2):
        a=np.zeros((4,4))
        for occupied in range(4):
            target,sign=thermal.car.word(occupied,[(mode,False)])
            if sign:a[target,occupied]=sign
        annih.append(np.kron(np.eye(r),a))
    ca=[np.kron(a,I) for a in annih]
    cb=[np.kron(P,a) for a in annih]
    allc=ca+cb
    car_error=max(float(np.linalg.norm(a@b.conj().T+b.conj().T@a
                   -(np.eye(nloc*nloc) if i==j else 0),2))
                  for i,a in enumerate(allc) for j,b in enumerate(allc))
    assert car_error<1e-13
    hop=.23*sum(a.conj().T@b+b.conj().T@a for a,b in zip(ca,cb))
    distances=original.distance_squared(phi[:,None,:],phi[None,:,:])
    grad_values=eps/2*distances
    grad_diag=np.broadcast_to(grad_values[:,None,:,None],(r,4,r,4)).ravel()
    parts=[np.kron(localT,I)+np.kron(I,localT),
           np.kron(localW,I)+np.kron(I,localW),
           np.diag(grad_diag),np.kron(B,I)+np.kron(I,B),hop]
    q=np.kron(v,np.eye(4));Q=np.kron(q,q)
    parts=[Q.conj().T@a@Q for a in parts];powers=(-6,6,2,0,-2)
    H=sum(parts);G=sum(k*a for k,a in zip(powers,parts));C=sum(k*k*a for k,a in zip(powers,parts))
    effects=[q.conj().T@np.diag(np.repeat(np.sqrt(.5+sgn*np.sin(s)/4),4))@q for sgn in (1,-1)]
    local_energy=np.repeat(e-e[0]+1,4);full_energy=(local_energy[:,None]+local_energy[None,:]).ravel()
    local_effects=[[np.kron(L,I) for L in effects],[np.kron(I,L) for L in effects]]
    assert np.linalg.norm(H@np.kron(P,P)-np.kron(P,P)@H)<1e-12
    assert np.linalg.norm(parts[-1]@np.kron(P,I)-np.kron(P,I)@parts[-1])>1e-3
    return dict(parts=parts,powers=powers,H=H,G=G,C=C,L=local_effects,nloc=nloc,
        e=e,local_energy=local_energy,full_energy=full_energy,parity=parity,
        CAR_error=car_error,minimum_F=float(F.min()),eps=eps,r=r)


def keep_data(d,bands):
    rank=4*bands
    while rank<d['nloc'] and abs(d['local_energy'][rank]-d['local_energy'][rank-1])<1e-10:rank+=1
    keep=np.arange(d['nloc'])<rank
    ground={1:0,-1:1}
    tails=[(j,ground[int(d['parity'][j])]) for j in np.where(~keep)[0]]
    return keep,tails


def local_map(x,d,bands,side):
    keep,tails=keep_data(d,bands);l=d['nloc'];t=x.reshape(l,l,l,l)
    if side==0:
        out=t*keep[:,None,None,None]*keep[None,None,:,None]
        for j,g in tails:out[g,:,g,:]+=t[j,:,j,:]
    else:
        out=t*keep[None,:,None,None]*keep[None,None,None,:]
        for j,g in tails:out[:,g,:,g]+=t[:,j,:,j]
    return out.reshape(x.shape)


def make_jets(d,plus=1.,minus=-.4):
    prep,logs=jets.geometry.statejets(jets.expjet(d['H'],d['G'],d['C'],BETA))
    U=[(jets.expjet(d['H'],plus*d['G'],plus*plus*d['C'],1j*t/HBAR),
        jets.expjet(d['H'],minus*d['G'],minus*minus*d['C'],1j*t/HBAR)) for t in TIMES]
    return prep,U


def histories_jets(d,bands,prepared=None,plus=1.,minus=-.4,apply_after_reads=True):
    prep,U=prepared if prepared is not None else make_jets(d,plus,minus)
    allrows=[];states=[]
    for hist in itertools.product((0,1),repeat=2):
        x={key:prep[key[0]].copy() if key[1]==0 else np.zeros_like(prep[0]) for key in KEYS}
        for side,(r,(up,um)) in enumerate(zip(hist,U)):
            y={}
            for j,k in KEYS:
                total=np.zeros_like(prep[0])
                for a in range(k+1):
                    for b in range(k-a+1):
                        c=k-a-b
                        factor=math.factorial(k)/(math.factorial(a)*math.factorial(b)*math.factorial(c))
                        total+=factor*up[a]@x[(j,b)]@um[c].conj().T
                L=d['L'][side][r];total=L@total@L.conj().T
                y[(j,k)]=local_map(total,d,bands,side) if apply_after_reads else total
            x=y
        allrows.append([np.trace(x[key]) for key in KEYS]);states.append(x[(0,0)])
    return np.array(allrows),states


def histories_value(d,bands,xi,epsilon,plus=1.,minus=-.4):
    def h(u):return sum(np.exp(k*u)*a for k,a in zip(d['powers'],d['parts']))
    rho=thermal.exp_h(h(xi),BETA);rho/=np.trace(rho)
    U=[(thermal.exp_h(h(plus*epsilon),1j*t/HBAR),
        thermal.exp_h(h(minus*epsilon),1j*t/HBAR)) for t in TIMES]
    values=[]
    for hist in itertools.product((0,1),repeat=2):
        x=rho.copy()
        for side,(r,(up,um)) in enumerate(zip(hist,U)):
            L=d['L'][side][r];x=L@up@x@um.conj().T@L.conj().T
            x=local_map(x,d,bands,side)
        values.append(np.trace(x))
    return np.array(values)


def local_model_check(d):
    prep,_=make_jets(d);rho=prep[0];rows=[]
    l=d['nloc'];P=np.diag(d['parity']);I=np.eye(l)
    for bands in (1,2,3,4):
        keep,tails=keep_data(d,bands)
        K=[np.diag(keep.astype(float))]
        for j,g in tails:
            k=np.zeros((l,l));k[g,j]=1;K.append(k)
        tp=float(np.linalg.norm(sum(k.conj().T@k for k in K)-I))
        even=max(float(np.linalg.norm(k@P-P@k)) for k in K)
        ref_before=np.einsum('abad->bd',rho.reshape(l,l,l,l))
        after=local_map(rho,d,bands,0)
        ref_after=np.einsum('abad->bd',after.reshape(l,l,l,l))
        remote=float(np.linalg.norm(ref_before-ref_after))
        expected=np.zeros_like(rho)
        for k in K:
            full=np.kron(k,I);expected+=full@rho@full.conj().T
        formula=float(np.linalg.norm(expected-after))
        m2_before=float(np.dot(d['full_energy']**2,np.diag(rho).real))
        m2_after=float(np.dot(d['full_energy']**2,np.diag(after).real))
        assert max(tp,even,remote,formula)<2e-12 and m2_after<=m2_before+1e-9
        rows.append(dict(radial_bands=bands,local_main_rank=int(keep.sum()),
            trace_preservation_error=tp,even_Kraus_error=even,remote_state_error=remote,
            independent_Kraus_formula_error=formula,
            comparison_second_moment_before=m2_before,comparison_second_moment_after=m2_after))
    step=2e-5
    hp=sum(np.exp(k*step)*a for k,a in zip(d['powers'],d['parts']))
    hm=sum(np.exp(-k*step)*a for k,a in zip(d['powers'],d['parts']))
    contact_error=float(np.linalg.norm((hp+hm-2*d['H'])/(step*step)-d['C']))
    assert contact_error<2e-4
    return dict(rows=rows,CAR_error=d['CAR_error'],contact_independent_difference_error=contact_error,
        neutral_hopping_norm=float(np.linalg.norm(d['parts'][-1],2)),
        original_geodesic_edge_norm=float(np.linalg.norm(d['parts'][2],2)),
        contact_includes_gradient_and_hopping=True)


def history_check(d):
    prepared=make_jets(d);full,states=histories_jets(d,4,prepared)
    rows=[]
    for bands in (1,2,3,4):
        z,out=histories_jets(d,bands,prepared)
        assert abs(z[:,0].sum()-1)<2e-11
        cq=sum(thermal.norm1(a-b) for a,b in zip(out,states))
        rows.append(dict(radial_bands=bands,source_jet_max_branch_errors=np.max(abs(z-full),axis=0).tolist(),
            cq_state_error=cq,mixed_total=[float(z[:,4].sum().real),float(z[:,4].sum().imag)]))
    assert max(rows[-1]['source_jet_max_branch_errors'])<1e-11
    assert rows[1]['source_jet_max_branch_errors'][-1]>1e-3
    target,_=histories_jets(d,2,prepared);checks=[]
    for step in (.0008,.0004):
        pp=histories_value(d,2,step,step);pm=histories_value(d,2,step,-step)
        mp=histories_value(d,2,-step,step);mm=histories_value(d,2,-step,-step)
        fd=(pp-pm-mp+mm)/(4*step**2)
        checks.append(dict(step=step,mixed_error=float(np.max(abs(fd-target[:,4])))))
    assert checks[1]['mixed_error']<checks[0]['mixed_error']/3
    same,_=histories_jets(d,2,plus=1.,minus=1.)
    norm=float(np.max(abs(same.sum(axis=0)-np.array([1,0,0,0,0,0]))))
    assert norm<3e-9
    # Erasing source derivatives of the regulator's post-read state is NOT allowed.
    bare,_=histories_jets(d,2,prepared,apply_after_reads=False)
    changed=float(np.max(abs(bare[:,5]-target[:,5])));assert changed>1e-3
    return dict(jet_order=[list(k) for k in KEYS],rows=rows,
        full_mixed_total=[float(full[:,4].sum().real),float(full[:,4].sum().imag)],
        independent_mixed_checks=checks,equal_history_joint_normalization_error=norm,
        ignoring_local_regulation_changes_second_source_by=changed,
        process_uses_original_full_diagnostic_H_between_every_record=True,
        no_claim_of_uniform_cutoff_rate_or_monotone_intermediate_errors=True)


def run():
    d=data()
    deps=('research_note_623.md','research_note_625.md','research_note_637.md','research_note_704.md',
        'research_note_705.md','round706_drafts/local_domain_entry.md',
        'round706_drafts/parity_completion_705.md','joint_gibbs_preserving_transfer.py',
        'joint_preparation_history_limit.py','joint_full_spatial_metric.py','research_note_667.md')
    return dict(date='2026-10-02',round=706,tests_run=2,failures=0,errors=0,
        calibration=dict(total_dimension=len(d['H']),node_dimension=d['nloc'],beta=BETA,hbar=HBAR,
            times=list(TIMES),radial_nodes_per_axis=2,site_volume=.8,lattice_spacing=d['eps'],
            minimum_F=d['minimum_F'],declared_neutral_hopping=.23,
            frozen_aligned_Higgs_and_gauge_link_conditional_diagnostic=True,
            not_full_Gauss_or_spatial_continuum_simulation=True),
        local_model=local_model_check(d),history=history_check(d),
        dependency_hashes={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in deps},
        scope=dict(original_complete_fixed_graph_theorem=True,
            graded_boundary_retaining_local_regularization=True,
            original_real_time_H_and_actual_record_instruments=True,
            scalar_joint_preparation_and_dynamic_sources_through_total_order_two=True,
            original_Gibbs_preparation_not_effective_regulator_Gibbs=True,
            full_process_remains_infinite_dimensional=True,
            no_uniform_spatial_limit_or_quantum_gravity_claim=True))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(round=706,checks=2,history=result['history'],all_checks_passed=True)))

