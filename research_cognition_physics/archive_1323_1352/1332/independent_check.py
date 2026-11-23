"""Independent 1332 reference-channel and actual two-fermion marker audit.

Uses only our frozen independent 1331 tuple-CAR and sparse parent builder.
Reference polar transport is obtained by rectangular SVD, not author code.
All reference directions and zero singular channels are retained in parent H.
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
DEP=HERE.parent/'1331'/'independent_check.py'
EXPECTED='6508090952adda9b2d5da29aa6a81dd93dbde6251e7d50c4c3972805a7b642d8'
assert hashlib.sha256(DEP.read_bytes()).hexdigest()==EXPECTED
spec=importlib.util.spec_from_file_location('our_frozen1331',DEP)
old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
Op=old.Op


def spin(r,s):
    states=list(itertools_product(r,s));ix={a:i for i,a in enumerate(states)}
    jp=np.zeros((len(states),len(states)),complex)
    for col,(hx,hy) in enumerate(states):
        if hx:jp[ix[(hx-1,hy)],col]=math.sqrt((r-hx+1)*hx)
        if hy:jp[ix[(hx,hy-1)],col]=math.sqrt((s-hy+1)*hy)
    jz=np.diag([(r+s)/2-hx-hy for hx,hy in states])
    casimir=jz@jz+(jp@jp.conj().T+jp.conj().T@jp)/2
    return jp,jz,casimir


def itertools_product(r,s):
    return ((x,y) for x in range(r+1) for y in range(s+1))


@functools.lru_cache(None)
def transport(r,s):
    assert r>=1 and s>=0
    ins=list(itertools_product(r,s));outs=list(itertools_product(r-1,s+1));oi={v:i for i,v in enumerate(outs)}
    t=np.zeros((len(outs),len(ins)),complex)
    for col,(hx,hy) in enumerate(ins):
        if r-hx:t[oi[(hx,hy)],col]+=math.sqrt((r-hx)*(s-hy+1))
        if hx:t[oi[(hx-1,hy+1)],col]+=math.sqrt(hx*(hy+1))
    u,singular,vh=np.linalg.svd(t,full_matrices=False)
    keep=singular>1e-11
    w=u[:,keep]@vh[keep,:]
    a=t/math.sqrt(r*(s+1))
    return dict(T=t,A=a,W=w,ins=ins,outs=outs,singular=singular,rank=int(keep.sum()))


@functools.lru_cache(None)
def kernel(kind,r,s,weight):
    if weight==0:return np.eye((r+1)*(s+1)),list(itertools_product(r,s)),list(itertools_product(r,s))
    first=transport(r,s)
    if weight==1:return first[kind],first['ins'],first['outs']
    assert weight==2 and r>=2
    second=transport(r-1,s+1)
    return second[kind]@first[kind],first['ins'],second['outs']


def build(N,kind,md,mu,edges,p):
    zero=np.zeros((1,1),complex)
    a=old.model(1,2,N,[np.array([[v]]) for v in md],[np.array([[v]]) for v in mu],
                [zero]*4,p)
    triples=[]
    sources={role:[] for role in range(4)}
    for col,(occupation,hx,hy) in enumerate(a['states']):
        r,s=map(int,a['bs'][col])
        for role,weight in enumerate((0,1,1,2)):
            moved=old.move(occupation,role,4+role)
            if moved is None:continue
            changed,sign=moved
            k,ins,outs=kernel(kind,r,s,weight);source_col=ins.index((hx,hy))
            for row,(hxx,hyy) in enumerate(outs):
                coefficient=k[row,source_col]
                if abs(coefficient)>2e-14:
                    rr=a['ix'][(changed,hxx,hyy)]
                    sources[role].append((rr,col,sign*coefficient))
                    triples.append((rr,col,edges[role]*sign*coefficient))
    forward=Op(len(a['states']),triples);edge=forward.add(forward.adj())
    a.update(H=a['H0'].add(a['V'],edge),edge=edge,edge_forward=forward,
             edge_sources={k:Op(len(a['states']),v) for k,v in sources.items()})
    return a


def norm(a):return old.norm(a)


def compute():
    checks=[]
    def eq(name,x,y,tol=3e-9):
        r=float(np.linalg.norm(np.asarray(x)-np.asarray(y)))
        checks.append(dict(name=name,residual=r,tolerance=tol))
        if r>tol:raise AssertionError((name,r,tol))
    def le(name,x,y,tol=3e-9):
        r=max(0.,float(x-y));checks.append(dict(name=name,residual=r,tolerance=tol,value=float(x),bound=float(y)))
        if r>tol:raise AssertionError((name,x,y))

    channels=[]
    for r in range(1,6):
        for s in range(5):
            b=transport(r,s);t,a,w=b['T'],b['A'],b['W'];jp,jz,c=spin(r,s);op,oz,oc=spin(r-1,s+1)
            maximum=(r+s)/2;tag=f'ref.r{r}.s{s}'
            eq(tag+'.Casimir',t.conj().T@t,r*(s+1)*np.eye(t.shape[1])-maximum*(maximum+1)*np.eye(t.shape[1])+c)
            eq(tag+'.A_scalar',op@a,a@jp)
            eq(tag+'.W_scalar',op@w,w@jp)
            eq(tag+'.W_reverse_scalar',jp@w.conj().T,w.conj().T@op)
            eq(tag+'.W_partial_isometry',w@w.conj().T@w,w)
            input_eigen=[];coeff=[]
            for j2 in range(abs(r-s),r+s+1,2):
                j=j2/2;lam=(j+(r-s)/2)*(j-(r-s)/2+1)
                input_eigen.extend([lam]*(j2+1))
                coeff.append(dict(j=j,lambda_j=lam,A=math.sqrt(max(0,lam)/(r*(s+1))),W=int(lam>0)))
            eq(tag+'.channel_eigenvalues',np.sort(np.linalg.eigvalsh(t.conj().T@t)),np.sort(input_eigen))
            expected_missing=s-r+1 if r<=s else 0
            eq(tag+'.missing_input_dimension',t.shape[1]-b['rank'],expected_missing)
            eq(tag+'.highest_calibration_A',a[0,0],1)
            eq(tag+'.highest_calibration_W',w[0,0],1)
            le(tag+'.A_norm',norm(a),1.)
            le(tag+'.W_norm',norm(w),1.)
            if r>1 and s>0:
                row=b['outs'].index((0,1));col=b['ins'].index((1,0))
                eq(tag+'.actual_one_blue_A',a[row,col],1/math.sqrt(r*(s+1)))
                eq(tag+'.actual_one_blue_W',w[row,col],1/(math.sqrt(r*(s+1))+math.sqrt((r-1)*s)))
            if r>=2:
                for kind in ('A','W'):
                    two,_,_=kernel(kind,r,s,2);second=transport(r-1,s+1)[kind]
                    eq(tag+'.'+kind+'_two_real_steps',two,second@b[kind])
                    eq(tag+'.'+kind+'_true_reverse',two.conj().T,b[kind].conj().T@second.conj().T)
            channels.append(dict(r=r,s=s,input_dimension=t.shape[1],output_dimension=t.shape[0],rank=b['rank'],coefficients=coeff))

    md=[.28+.03j,.23-.02j];mu=[.19-.01j,.17+.025j]
    edges=[.21,.55,.55,.27];p=[.04,-.02]
    cases=[]
    for N in [(5,5),(6,5)]:
        models={k:build(N,k,md,mu,edges,p) for k in ('A','W')}
        a=models['A'];dim=len(a['states']);P=a['P'];selected=a['selected']
        initial_index=a['ix'][((2,5),0,0)]
        marker_index=a['ix'][((0,1),0,1)]
        unknown_cols=[a['ix'][(tuple(sorted((i,4+j))),0,0)] for i in (1,2) for j in (1,2)]
        F=np.eye(dim,dtype=complex)[:,unknown_cols]
        actual_marker_input=unknown_cols.index(initial_index)
        E=np.zeros((dim,28),complex)
        for col,s in enumerate(sorted({x[0] for x in a['states']})):E[a['ix'][(s,0,0)],col]=1
        eq(f'N{N}.edge_same_whole_polarized_code',models['A']['edge'].act(E),models['W']['edge'].act(E))
        rolecal=[]
        for role in range(4):
            sa=models['A']['edge_sources'][role];sw=models['W']['edge_sources'][role]
            eq(f'N{N}.same_calibrated_edge_source{role}',sa.act(E),sw.act(E))
            eq(f'N{N}.same_calibrated_reverse_source{role}',sa.adj().act(E),sw.adj().act(E))
            rolecal.append(float(norm(sa.act(E))))
        outputs={};details={}
        for kind,b in models.items():
            tag=f'N{N}.{kind}';hP=b['H'].act(P);H=P.conj().T@hP
            eq(tag+'.invariant_Z',b['H'].v*(b['z'][b['H'].r]-b['z'][b['H'].c]),0)
            eq(tag+'.actual_support',hP-P@H,0)
            eq(tag+'.full_SU2_columns',b['H'].act(b['Jx'].act(P))-b['Jx'].act(hP),0)
            eq(tag+'.reference_balance',b['bs']+b['counts'][:,:,1]+b['counts'][:,:,2]+2*b['counts'][:,:,3],np.asarray(N)[None,:]+np.zeros((dim,2)))
            # Nontrivial full-parent probe, including channels outside actual P.
            rng=np.random.default_rng(1332);probe=rng.normal(size=(dim,3))+1j*rng.normal(size=(dim,3));probe/=np.linalg.norm(probe,axis=0)
            eq(tag+'.full_SU2_parent',b['H'].act(b['Jx'].act(probe))-b['Jx'].act(b['H'].act(probe)),0)
            # All CAR paths to the declared record, without fitting the result.
            first=b['H'].act(np.eye(dim,dtype=complex)[:,initial_index]);back=b['H'].act(np.eye(dim,dtype=complex)[:,marker_index]).conj()
            paths=back*first
            nonzero=np.flatnonzero(np.abs(paths)>1e-12)
            assert len(nonzero)==1
            intermediate=b['states'][nonzero[0]]
            assert intermediate==((0,5),1,0)
            expected_reference=1/math.sqrt(N[0]*N[1]) if kind=='A' else 1/(math.sqrt(N[0]*N[1])+math.sqrt((N[0]-1)*(N[1]-1)))
            h2=complex(paths.sum());expected=md[0].conjugate()*edges[1]/math.sqrt(N[0])*expected_reference
            eq(tag+'.all_paths_H2_marker',h2,expected)
            eq(tag+'.no_first_order_marker',first[marker_index],0)
            eig,vec=np.linalg.eigh(H);hn=float(np.max(np.abs(eig)));Fp=P.conj().T@F
            # Each reference transport is contractive, each CAR branch has norm
            # at most one. This rational ceiling is independent of eigensolvers.
            analytic_H_bound=7.1
            raw_bound=4*sum(abs(z.real)+abs(z.imag) for z in md+mu)+2*sum(abs(z) for z in edges)+2*max(map(abs,p))
            le(tag+'.analytic_parent_norm_ceiling',hn,raw_bound)
            le(tag+'.rational_parent_ceiling',raw_bound,analytic_H_bound)
            times=[];out=[]
            for t in (1.,1.7,2.4):
                evolved=P@((vec*np.exp(-1j*eig*t))@(vec.conj().T@Fp))
                psi=evolved[:,actual_marker_input]
                probability=float(abs(psi[marker_index])**2)
                # Fixed degree-120 polynomial with a separate analytic truncation bound.
                v=Fp[:,actual_marker_input].copy();polynomial=v.copy()
                for k in range(1,121):v=(-1j*t/k)*(H@v);polynomial+=v
                approx=complex((P@polynomial)[marker_index]);z=analytic_H_bound*t
                tail=math.exp(z+(121*math.log(z) if z else -math.inf)-math.lgamma(122))
                amplitude_error=abs(approx-psi[marker_index]);le(tag+f'.t{t}.finite_remainder',amplitude_error,tail,tol=3e-12)
                radius=tail+3e-12 # numerical guard is not formal interval arithmetic
                interval=[max(0.,abs(approx)-radius)**2,min(1.,abs(approx)+radius)**2]
                assert interval[0]<=probability+1e-13 and probability<=interval[1]+1e-13
                eq(tag+f'.t{t}.unknown_partner_record_gram',evolved.conj().T@evolved,np.eye(4))
                work=-evolved.conj().T@b['V'].add(b['edge']).act(evolved)
                after=evolved.conj().T@(b['hdiag'][:,None]*evolved)
                eq(tag+f'.t{t}.quench_complete_record_energy',after-F.conj().T@b['H'].act(F),work)
                for x in range(2):
                    balance=b['bs'][:,x]+b['counts'][:,x,1]+b['counts'][:,x,2]+2*b['counts'][:,x,3]
                    eq(tag+f'.t{t}.actual_reference_ledger{x}',evolved.conj().T@(balance[:,None]*evolved),N[x]*np.eye(4))
                source_values={}
                for role,source in b['edge_sources'].items():
                    for suffix,op in [('real',source.add(source.adj())),('imag',source.add(source.adj().scaled(-1)).scaled(1j))]:
                        value=evolved.conj().T@op.act(evolved)
                        source_values[f'{role}_{suffix}']=[float(np.trace(value).real/4),float(norm(value))]
                jx=evolved.conj().T@b['Jx'].act(evolved)
                jy=evolved.conj().T@b['Jy'].act(evolved)
                eq(tag+f'.t{t}.closed_parent_Jx',jx,F.conj().T@b['Jx'].act(F))
                eq(tag+f'.t{t}.closed_parent_Jy',jy,F.conj().T@b['Jy'].act(F))
                # Raw local occupation record destroys system transverse coherence,
                # assigned to measurement source; all branches and pointer retained.
                probabilities=np.abs(psi)**2
                eq(tag+f'.t{t}.all_record_probability',probabilities.sum(),1.)
                times.append(dict(time=t,marker_probability=probability,marker_amplitude=[psi[marker_index].real,psi[marker_index].imag],
                    certificate_interval=interval,taylor_tail_bound=tail,numerical_guard=3e-12,
                    all_records_probability_sum=float(probabilities.sum()),marker_branch_is_not_postselected=True,
                    mean_switch_work=float(np.trace(work).real/4),record_source_Jx_operator_norm=norm(jx),record_source_Jy_operator_norm=norm(jy),
                    complete_edge_sources=source_values))
                out.append((probability,interval))
            details[kind]=dict(physical_dimension=dim,actual_support_dimension=P.shape[1],unknown_input_dimension=4,
                full_parent_H_norm_on_support=hn,analytic_parent_H_bound=analytic_H_bound,
                second_order_amplitude=[h2.real,h2.imag],single_intermediate_state=[list(intermediate[0]),intermediate[1],intermediate[2]],
                expected_reference_element=expected_reference,held_out=times)
            outputs[kind]=out
        # A separate mathematical finite witness uses only exact second-order
        # coefficients and the analytic third-order tail; no numerical expm.
        small_t=1e-5;small_tail=math.exp(7.1*small_t)*(7.1*small_t)**3/6
        small_intervals={}
        for kind in ('A','W'):
            second=complex(*details[kind]['second_order_amplitude'])
            leading=abs(second)*small_t**2/2
            small_intervals[kind]=[max(0.,leading-small_tail)**2,(leading+small_tail)**2]
        small_gap=small_intervals['A'][0]-small_intervals['W'][1]
        assert small_gap>0
        comparisons=[]
        for index,t in enumerate((1.,1.7,2.4)):
            a_gap=outputs['A'][index][1][0]-outputs['W'][index][1][1]
            w_gap=outputs['W'][index][1][0]-outputs['A'][index][1][1]
            gap=max(a_gap,w_gap)
            assert gap>0,(N,t,outputs)
            accuracy=1e-6;failure=.01
            systematic=1e-7
            repetitions=math.ceil(math.log(4/failure)/(2*accuracy*accuracy))
            local_outcomes=[len(set((tuple(i-4*x for i in s if 4*x<=i<4*x+4),hs[x]) for s,*hs in a['states'])) for x in range(2)]
            bits=sum(math.ceil(math.log2(v)) for v in local_outcomes)
            margin=gap-2*(accuracy+systematic)
            assert margin>0
            comparisons.append(dict(time=t,higher_probability_candidate='A' if a_gap>w_gap else 'W',
                certified_probability_gap=gap,single_mean_precision=accuracy,
                registered_single_mean_instrument_bias=systematic,joint_failure_probability=failure,
                source_or_statistics_free=False,separation_margin=margin,repetitions_per_candidate=repetitions,
                reference_quanta_per_trial=sum(N)-2,total_reference_quanta=2*repetitions*(sum(N)-2),
                local_record_dimensions=local_outcomes,blank_bits_per_trial=bits,total_blank_bits=2*repetitions*bits))
        cases.append(dict(N=list(N),calibration_refits=0,initial_record=[['L1_x','L0_y'],[N[0]-1,N[1]-1]],
            marker='D_x and L0_x, no matter at y, reference b1_x=0,b1_y=1',models=details,comparisons=comparisons,
            analytic_small_time_certificate=dict(time=small_t,H_bound=7.1,amplitude_remainder=small_tail,
                exact_second_order_radicals_verified=True,probability_intervals=small_intervals,
                probability_gap=small_gap,floating_values_display_exact_formula_certificate=True)))
    return dict(round=1332,status='PASS',new_axioms=0,stage_complete=False,
        dependency_sha256={str(DEP.relative_to(HERE.parent)).replace('\\','/'):EXPECTED},
        classification='Reference-only SU2-scalar transport times fixed matter bilinear; not all total scalar parents.',
        channels=channels,cases=cases,check_count=len(checks),checks=checks,max_residual=max(x['residual'] for x in checks),
        same_highest_gain_is_rule_contract_not_inferred_from_one_shot=True,
        U_uses_two_updated_steps=True,zero_channels_kept_in_unitary_parent=True,
        author_imports=False,preferred_transport_selected=False,
        after_this_round_stop_repairing_transport=True)


def compare(a,b,path='root'):
    if isinstance(a,dict):
        assert a.keys()==b.keys(),path
        for k in a:compare(a[k],b[k],path+'.'+k)
    elif isinstance(a,list):
        assert len(a)==len(b),path
        for i,(x,y) in enumerate(zip(a,b)):compare(x,y,path+f'[{i}]')
    elif isinstance(a,float):assert math.isfinite(b) and abs(a-b)<=1e-8*(1+abs(a)),(path,a,b)
    else:assert a==b,(path,a,b)


def main():
    p=argparse.ArgumentParser();p.add_argument('--write-new',action='store_true');p.add_argument('--check',action='store_true')
    args=p.parse_args();r=compute();dest=HERE/'independent_results.json'
    if args.write_new:
        with dest.open('x',encoding='utf-8') as f:json.dump(r,f,ensure_ascii=False,indent=2);f.write('\n')
    else:compare(r,json.loads(dest.read_text(encoding='utf-8')))
    print(json.dumps({'round':1332,'status':r['status'],'checks':r['check_count'],'max_residual':r['max_residual'],
        'cases':[{'N':x['N'],'comparisons':x['comparisons']} for x in r['cases']]},ensure_ascii=False))


if __name__=='__main__':main()
