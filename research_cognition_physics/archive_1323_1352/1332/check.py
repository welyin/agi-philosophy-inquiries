"""1332: scalar reference transport channels and finite A/W predictions.

All parent reference directions retained. Frozen 1330 CAR/vertex implementation
is a declared dependency. Exact Z conservation alone permits dynamics reduction.
"""
from __future__ import annotations
import argparse
import functools
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
DEPENDENCY=HERE.parent/'1330'/'check.py'
spec=importlib.util.spec_from_file_location('frozen1330',DEPENDENCY)
old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)


def norm(x):return float(np.linalg.norm(x,2)) if x.size else 0.


@functools.lru_cache(None)
def refblock(r,s):
    inp=[(i,j) for i in range(r+1) for j in range(s+1)]
    out=[(i,j) for i in range(r) for j in range(s+2)]
    oi={x:k for k,x in enumerate(out)}
    t=np.zeros((len(out),len(inp)),complex)
    for col,(i,j) in enumerate(inp):
        if r-i:t[oi[(i,j)],col]+=math.sqrt((r-i)*(s-j+1))
        if i:t[oi[(i-1,j+1)],col]+=math.sqrt(i*(j+1))
    u,sv,vh=np.linalg.svd(t,full_matrices=False)
    keep=sv>1e-10
    w=u[:,keep]@vh[keep]
    a=t/math.sqrt(r*(s+1))
    return dict(inp=inp,out=out,T=t,A=a,W=w,sv=sv,rank=int(keep.sum()))


def reference_J2(r,s):
    basis=[(i,j) for i in range(r+1) for j in range(s+1)]
    lookup={x:k for k,x in enumerate(basis)};d=len(basis)
    jp=np.zeros((d,d),complex);z=np.zeros(d)
    for col,(i,j) in enumerate(basis):
        z[col]=(r+s)/2-i-j
        if i:jp[lookup[(i-1,j)],col]=math.sqrt((r-i+1)*i)
        if j:jp[lookup[(i,j-1)],col]=math.sqrt((s-j+1)*j)
    return np.diag(z*z)+(jp@jp.conj().T+jp.conj().T@jp)/2


def channel_checks():
    rows=[]
    for r in range(1,7):
        for s in range(7):
            b=refblock(r,s);jmax=(r+s)/2;j2=reference_J2(r,s)
            gram=b['T'].conj().T@b['T'];expected=(r*(s+1)-jmax*(jmax+1))*np.eye(len(gram))+j2
            casimir=norm(gram-expected)
            identity=norm(b['W'].conj().T@b['W']@b['W'].conj().T@b['W']-b['W'].conj().T@b['W'])
            top=max(norm(b[k][:,0]-np.eye(len(b['out']))[:,0]) for k in ['A','W'])
            input_js=np.arange(abs(r-s)/2,jmax+.1)
            predicted=[]
            for jj in input_js:
                lam=r*(s+1)-jmax*(jmax+1)+jj*(jj+1)
                predicted += [lam]*int(2*jj+1)
            spectral=float(np.max(np.abs(np.linalg.eigvalsh(gram)-np.sort(predicted))))
            missing=int(r<=s);actual_missing=(r+1)*(s+1)-b['rank']
            expected_missing=(s-r+1) if missing else 0
            marker_error=0.
            if r>1 and s>0:
                col=b['inp'].index((1,0));row=b['out'].index((0,1))
                aa=1/math.sqrt(r*(s+1));ww=1/(math.sqrt(r*(s+1))+math.sqrt((r-1)*s))
                marker_error=max(abs(b['A'][row,col]-aa),abs(b['W'][row,col]-ww))
            passed=max(casimir,identity,top,spectral,marker_error)<2e-9 and actual_missing==expected_missing
            rows.append(dict(r=r,s=s,input_dimension=(r+1)*(s+1),output_dimension=r*(s+2),
                common_rank=b['rank'],input_kernel_dimension=actual_missing,casimir_error=casimir,
                polar_projector_error=identity,top_gain_error=top,spectrum_error=spectral,
                one_excitation_marker_error=float(marker_error),passed=bool(passed)))
    return rows


def scalar_parent(q,Ns):
    a=old.physical(q,Ns);d=len(a['basis']);look=a['lookup']
    edges={menu:{s:np.zeros((d,d),complex) for s in ['D','L','U']} for menu in ['A','W']}
    z=np.zeros(d,int);jp=np.zeros((d,d),complex)
    low=np.zeros((d,d),complex)
    for col,(s,h0,h1) in enumerate(a['basis']):
        r,b=int(a['bs'][col,0]),int(a['bs'][col,1])
        z[col]=h0+h1+sum(((s>>(4*x+3))&1)+((s>>(4*x+1))&1) for x in range(2))
        for x,hh in enumerate([h0,h1]):
            bb=int(a['bs'][col,x])
            if hh:
                hs=[h0,h1];hs[x]-=1;jp[look[(s,*hs)],col]+=math.sqrt((bb-hh+1)*hh)
            moved=old.hop(s,4*x+2,4*x+3)
            if moved:
                ss,sign=moved;jp[look[(ss,h0,h1)],col]+=sign
        for sec,mm,w in [('D',[0],0),('L',[2,3],1),('U',[1],2)]:
            for m in mm:
                moved=old.hop(s,m,4+m)
                if not moved:continue
                ss,sign=moved
                if not w:
                    for menu in edges:edges[menu][sec][look[(ss,h0,h1)],col]+=sign
                    continue
                first=refblock(r,b);icol=first['inp'].index((h0,h1))
                for menu in edges:
                    mat=first[menu];out=first['out']
                    if w==2:
                        second=refblock(r-1,b+1);mat=second[menu]@mat;out=second['out']
                    for row,(k0,k1) in enumerate(out):
                        value=mat[row,icol]
                        if abs(value)>1e-14:edges[menu][sec][look[(ss,k0,k1)],col]+=sign*value
    # This is a model diagnostic, not a freely supplied low-spin preparation or readout.
    for s in old.masks(q):
        indices=[look[(s,h0,h1)] for h0 in range(int(a['bs'][look[(s,0,0)],0])+1)
                 for h1 in range(int(a['bs'][look[(s,0,0)],1])+1)]
        r,b=map(int,a['bs'][indices[0]])
        j2=reference_J2(r,b);ev,vec=np.linalg.eigh(j2);top=(r+b)/2*((r+b)/2+1)
        keep=ev<top-1e-8
        low[np.ix_(indices,indices)]=vec[:,keep]@vec[:,keep].conj().T
    a.update(scalar_edges=edges,Z=z,Jplus=jp,Jx=(jp+jp.conj().T)/2,
             Jy=(jp-jp.conj().T)/(2j),Jz=sum(Ns)/2-z,reference_lowspin=low)
    return a


def assemble(a,menu,md,mu,jcoef,p):
    d=len(a['basis']);v=np.zeros((d,d),complex);edge=np.zeros_like(v)
    for x in range(2):
        t=md[x]*a['D'][x]+mu[x]*a['U'][x];v+=t+t.conj().T
    for sec in jcoef:
        t=jcoef[sec]*a['scalar_edges'][menu][sec];edge+=t+t.conj().T
    h0=np.diag(sum(p[x]*(a['occ'][:,x,0]+a['occ'][:,x,1]-a['occ'][:,x,2]) for x in range(2)))
    return v+edge+h0,v,edge,h0


def calibrate(md,mu,jcoef,p):
    a=scalar_parent(1,[3,3]);l=old.logical_model(1,[3,3]);out={}
    for menu in ['A','W']:
        h,_,_,_=assemble(a,menu,md,mu,jcoef,p);comp=a['E'].conj().T@h@a['E'];li=l['lookup']
        recovered_md=[];recovered_mu=[]
        for x in range(2):
            o=4*x
            recovered_md.append(comp[li[1<<o],li[1<<(o+2)]].conjugate())
            recovered_mu.append(-comp[li[1<<(o+1)],li[1<<(o+3)]].conjugate()/math.sqrt(2/3))
        recovered_j=comp[li[1<<2],li[1<<6]]
        err=max(*(abs(recovered_md[x]-md[x]) for x in range(2)),
                *(abs(recovered_mu[x]-mu[x]) for x in range(2)),abs(recovered_j-jcoef['L']))
        out[menu]=dict(md=[float(x.real) for x in recovered_md],mu=[float(x.real) for x in recovered_mu],
                       JL=float(recovered_j.real),recovery_error=float(err))
    return dict(q=1,Ns=[3,3],menus=out,all_balance_top_gain_is_candidate_rule=True,
                finite_calibration_does_not_prove_all_balance_rule=True,after_calibration_refits=0)


def polynomial_amplitude(h,psi,marker,t,order,bound):
    term=psi.copy();total=psi.copy()
    for k in range(1,order+1):
        term=(-1j*t/k)*(h@term);total+=term
    # Hermitian integral Taylor remainder uses unitarity, no exponential prefactor.
    tail=(abs(t)*bound)**(order+1)/math.factorial(order+1)
    amp=total[marker];radius=tail+2e-12
    lo=max(0.,abs(amp)-radius)**2;hi=(abs(amp)+radius)**2
    return dict(order=order,analytic_tail=float(tail),numerical_allowance=2e-12,
                machine_interval_proof=False,
                amplitude_real=float(amp.real),amplitude_imag=float(amp.imag),
                probability_interval=[float(lo),float(hi)])


def evolve(eigen,initial,t):
    w,v=eigen;return (v*np.exp(-1j*t*w))@(v.conj().T@initial)


def run_case(Ns,md,mu,jcoef,p):
    q=2;a=scalar_parent(q,Ns);d=len(a['basis']);l=old.logical_model(q,Ns)
    support=np.flatnonzero(a['Z']<=q);out=np.flatnonzero(a['Z']>q);look=a['lookup']
    initial_mask=(1<<3)|(1<<6);marker_mask=(1<<0)|(1<<2)
    initial_index=look[(initial_mask,0,0)];marker_index=look[(marker_mask,0,1)]
    ip=int(np.flatnonzero(support==initial_index)[0]);mp=int(np.flatnonzero(support==marker_index)[0])
    cols=[look[(s,0,0)] for s in old.masks(q)
          if all(old.counts(s)[x,0]==old.counts(s)[x,1]==0 and old.counts(s)[x,2]==1 for x in range(2))]
    initial=np.eye(d,dtype=complex)[:,cols];probe_col=cols.index(initial_index)
    inputp=initial[support];psip=inputp[:,probe_col]
    htarget,_,_=old.assemble(l,md,mu,jcoef,p)
    cs,ls=old.residual_parts(l,md,mu);e=a['E'];prediction=e@sum(cs)+sum(a['EX'][x]@ls[x] for x in range(2))
    eta=old.bound(q,Ns,md,mu)
    norm_bound=2*q*sum(abs(x) for x in md+mu)+2*(abs(jcoef['D'])+2*abs(jcoef['L'])+abs(jcoef['U']))+q*max(abs(x) for x in p)
    times=[.65,1.1,1.65];menus={};operators={}
    # Per-site detailed pointer counts are reported separately, without free compression.
    local_counts=[len({((s>>(4*x))&15,int(a['bs'][i,x]-hs[x]),int(hs[x]))
                       for i,(s,*hs) in enumerate(a['basis'])}) for x in range(2)]
    for menu in ['A','W']:
        h,v,edge,h0=assemble(a,menu,md,mu,jcoef,p);hp=h[np.ix_(support,support)]
        ep=edge[np.ix_(support,support)];eig=np.linalg.eigh(hp);operators[menu]=hp
        su2=norm(h@a['Jplus']-a['Jplus']@h)
        zres=float(np.linalg.norm(h*(a['Z'][None,:]-a['Z'][:,None])))
        inv=norm(h[np.ix_(out,support)])
        residual=norm(h@e-e@htarget-prediction)
        source_rows=[]
        for sec in ['D','L','U']:
            b=a['scalar_edges'][menu][sec];logical=l['edges'][sec].conj().T
            real=b+b.conj().T;imag=1j*(b-b.conj().T)
            source_rows.append(dict(sector=sec,real_source_intertwiner=norm(real@e-e@(logical+logical.conj().T)),
                                    imag_source_intertwiner=norm(imag@e-e@(1j*(logical-logical.conj().T)))))
        h2=(h@h[:,initial_index])[marker_index]
        r,s=Ns[0],Ns[1]-1
        k=1/math.sqrt(r*(s+1)) if menu=='A' else 1/(math.sqrt(r*(s+1))+math.sqrt((r-1)*s))
        predicted_h2=md[0]*jcoef['L']/math.sqrt(r)*k
        low_seed=float(np.real(v[:,initial_index].conj()@a['reference_lowspin']@v[:,initial_index]))
        holdouts=[]
        for t in times:
            sp=evolve(eig,inputp,t);full=np.zeros_like(initial);full[support]=sp;psi=full[:,probe_col]
            probs=np.abs(psi)**2;prob=float(probs[marker_index])
            cert=polynomial_amplitude(hp,psip,mp,t,60,norm_bound)
            low=float(np.real(psi.conj()@a['reference_lowspin']@psi))
            energy=norm(full.conj().T@h@full-initial.conj().T@h@initial)
            work=-full.conj().T@(h-h0)@full
            final=full.conj().T@h0@full
            workres=norm(final-initial.conj().T@h@initial-work)
            gram=norm(full.conj().T@full-np.eye(len(cols)))
            # Complete Fock occupation readout keeps Jz and changes transverse system moments.
            record_sources={name:float(norm(full.conj().T@a[name]@full)) for name in ['Jx','Jy']}
            refs=[float(np.real(psi.conj()@(a['bs'][:,x]*psi))) for x in range(2)]
            balances=[norm(full.conj().T@((a['bs'][:,x]+a['occ'][:,x,2]+2*a['occ'][:,x,1])[:,None]*full)
                           -Ns[x]*np.eye(len(cols))) for x in range(2)]
            source_expect={}
            for sec in ['D','L','U']:
                b=a['scalar_edges'][menu][sec]
                source_expect[sec]=dict(real=float(np.real(psi.conj()@(b+b.conj().T)@psi)),
                    phase=float(np.real(psi.conj()@(1j*jcoef[sec]*(b-b.conj().T))@psi)))
            holdouts.append(dict(time=t,marker_probability=prob,finite_Taylor_certificate=cert,
                reference_lowspin_probability=low,reference_numbers=refs,edge_sources=source_expect,
                complete_record_probabilities=probs.tolist(),record_probability_sum=float(probs.sum()),
                record_gram_error=gram,record_system_orientation_change_operator_norm=record_sources,
                switch_off_work_operator_norm=norm(work),source_energy_conservation_error=energy,
                complete_switch_energy_error=workres,local_balance_error=max(balances),
                local_blank_pointer_dimensions=local_counts,
                local_record_capacity_bits=sum(math.ceil(math.log2(x)) for x in local_counts),
                all_outcomes_retained=True,reference_or_record_reset=False,
                passed=bool(cert['probability_interval'][0]-1e-14<=prob<=cert['probability_interval'][1]+1e-14
                            and max(gram,energy,workres,*balances,abs(probs.sum()-1))<3e-9)))
        maxsource=max(max(x['real_source_intertwiner'],x['imag_source_intertwiner']) for x in source_rows)
        passed=max(su2,zres,inv,residual,abs(h2-predicted_h2),maxsource)<3e-9 and low_seed>0 and norm(hp)<=norm_bound
        menus[menu]=dict(global_SU2_commutator_error=su2,Z_commutator_error=zres,
             exact_invariant_support_error=inv,common_residual_error=residual,source_checks=source_rows,
             marker_H2=float(h2.real),predicted_marker_H2=float(predicted_h2),marker_H2_error=float(abs(h2-predicted_h2)),
             actual_point_generated_lowspin_seed=low_seed,holdouts=holdouts,
             passed=bool(passed and all(x['passed'] for x in holdouts)))
    comparisons=[]
    for ar,wr in zip(menus['A']['holdouts'],menus['W']['holdouts']):
        ai=ar['finite_Taylor_certificate']['probability_interval'];wi=wr['finite_Taylor_certificate']['probability_interval']
        gap=ai[0]-wi[1]
        assert gap>0,(Ns,ar['time'],ai,wi)
        eps=gap/8;systematic_each=gap/8;failure_per_mean=.005
        shots=math.ceil(math.log(2/failure_per_mean)/(2*eps*eps))
        comparisons.append(dict(time=ar['time'],probability_gap=ar['marker_probability']-wr['marker_probability'],
            numerically_checked_probability_gap=gap,absolute_precision_each=eps,
            systematic_instrument_budget_each=systematic_each,combined_tolerance=2*(eps+systematic_each),
            combined_failure_probability_bound=2*failure_per_mean,sufficient_repetitions_per_menu=shots,
            total_repreparations=2*shots,initial_reference_quanta_per_trial=sum(Ns)-2,
            finite_effective_counting_accuracy_is_registered_input=True,
            separation_margin=gap-2*(eps+systematic_each),passed=True))
    # Separate exact analytic small-time theorem; no eigensolver/Taylor evaluation needed.
    small_t=1e-5;tail2=(norm_bound*small_t)**3/6
    small_intervals={menu:[max(0.,menus[menu]['predicted_marker_H2']*small_t**2/2-tail2)**2,
                          (menus[menu]['predicted_marker_H2']*small_t**2/2+tail2)**2]
                     for menu in ['A','W']}
    small_gap=small_intervals['A'][0]-small_intervals['W'][1]
    assert small_gap>0
    return dict(q=q,Ns=Ns,parent_dimension=d,invariant_support_dimension=len(support),
        actual_unknown_input_dimension=len(cols),reference_lowspin_prepared=False,
        marker=dict(initial=[initial_mask,0,0],record=[marker_mask,0,1]),norm_bound=norm_bound,
        inherited_effective_process_generator_bound=eta,menus=menus,comparisons=comparisons,
        analytic_small_time_certificate=dict(time=small_t,norm_bound=norm_bound,
            exact_symbolic_second_coefficients='MD0*JL/sqrt(N0) times k_A or k_W from proof (12)',
            second_order_remainder_bound=tail2,probability_intervals=small_intervals,
            probability_gap_lower_bound=small_gap,uses_eigensolver_or_measured_fit=False),
        record_basis=[list(x) for x in a['basis']],passed=bool(all(x['passed'] for x in menus.values())))


def compute():
    md=[.23,.17];mu=[.11,.13];jcoef={'D':.19,'L':.4,'U':.21};p=[.04,-.025]
    channels=channel_checks();cal=calibrate(md,mu,jcoef,p)
    cases=[run_case(Ns,md,mu,jcoef,p) for Ns in [[5,5],[6,6]]]
    paths=[HERE/'proof.md',HERE/'check.py',DEPENDENCY,HERE.parent/'1330'/'proof.md',
           HERE.parent/'_admission'/'round1332_admission.md']
    provenance={str(p.relative_to(HERE.parent)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    passed=all(x['passed'] for x in channels+cases) and all(x['recovery_error']<1e-10 for x in cal['menus'].values())
    return dict(round=1332,status='PASS' if passed else 'FAIL',provenance=provenance,
                parameters=dict(MD=md,MU=mu,J=jcoef,p=p),channel_checks=channels,calibration=cal,cases=cases,
                new_axioms=0,empirical_experiments=0,all_possible_SU2_parents_classified=False,
                unique_winner_selected=False,stage_complete=False,long_goal_complete=False)


def compare(a,b,path='root'):
    if isinstance(a,dict):
        assert set(a)==set(b),(path,'keys')
        for k in a:compare(a[k],b[k],path+'.'+k)
    elif isinstance(a,list):
        assert len(a)==len(b),(path,'length')
        for i,(x,y) in enumerate(zip(a,b)):compare(x,y,path+f'[{i}]')
    elif isinstance(a,float):assert math.isclose(a,b,rel_tol=1e-7,abs_tol=2e-9),(path,a,b)
    elif isinstance(a,int) and 'repetition' in path:assert abs(a-b)<=max(1,int(abs(a)*1e-7)),(path,a,b)
    else:assert a==b,(path,a,b)


def main():
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group();g.add_argument('--check',action='store_true');g.add_argument('--write-new',action='store_true')
    args=p.parse_args();r=compute();assert r['status']=='PASS',r
    dst=HERE/'results.json'
    if args.write_new:
        payload=json.dumps(r,ensure_ascii=False,indent=2)+'\n'
        with dst.open('x',encoding='utf-8') as f:f.write(payload)
    else:compare(r,json.loads(dst.read_text(encoding='utf-8')))
    print(json.dumps(dict(status='PASS',round=1332,reference_channel_cases=42,full_parent_candidates=4,
                         finite_record_comparisons=6,results_written=args.write_new)))


if __name__=='__main__':main()
