"""1330: one two-site CAR/reference parent with conversion, compensated motion and records.

NumPy only. --write-new creates results.json exclusively; default/--check is read-only.
Every reference direction and unilateral-shift boundary is kept in the physical parent.
"""
from __future__ import annotations
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
WEIGHTS={'D':0,'L':1,'U':2}
MODES={'D':[0],'L':[2,3],'U':[1]}


def norm(a): return float(np.linalg.norm(a,2))
def frob(a): return float(np.linalg.norm(a))
def masks(q):
    return [sum(1<<i for i in c) for c in itertools.combinations(range(8),q)]


def hop(s,i,j):
    if not (s>>j)&1:return None
    sign=(-1)**((s&((1<<j)-1)).bit_count());t=s^(1<<j)
    if (t>>i)&1:return None
    sign*=(-1)**((t&((1<<i)-1)).bit_count())
    return t|(1<<i),sign


def counts(s):
    return np.array([[(s>>(4*x))&1,(s>>(4*x+1))&1,
                       ((s>>(4*x+2))&3).bit_count()] for x in range(2)])


def add_hop(matrix,lookup,s,col,i,j,coefficient=1):
    moved=hop(s,i,j)
    if moved:
        ss,sg=moved;matrix[lookup[ss],col]+=sg*coefficient


def logical_model(q,Ns):
    basis=masks(q);li={s:i for i,s in enumerate(basis)};d=len(basis)
    ds=[np.zeros((d,d),complex) for _ in range(2)]
    us=[np.zeros((d,d),complex) for _ in range(2)]
    ld=[np.zeros((d,d),complex) for _ in range(2)]
    lu=[np.zeros((d,d),complex) for _ in range(2)]
    edges={k:np.zeros((d,d),complex) for k in WEIGHTS}
    for col,s in enumerate(basis):
        for x in range(2):
            off=4*x
            add_hop(ds[x],li,s,col,off+2,off)
            add_hop(us[x],li,s,col,off+3,off+1,-1)
            add_hop(ld[x],li,s,col,off,off+3)
            add_hop(lu[x],li,s,col,off+2,off+1)
        for sector,mm in MODES.items():
            for m in mm:add_hop(edges[sector],li,s,col,4+m,m)
    occ=np.array([counts(s) for s in basis])
    bs=np.asarray(Ns)[None,:]-occ[:,:,2]-2*occ[:,:,1]
    return dict(basis=basis,lookup=li,D=ds,U=us,LD=ld,LU=lu,edges=edges,occ=occ,bs=bs)


def physical(q,Ns):
    lb=masks(q);basis=[]
    for s in lb:
        occ=counts(s);bs=np.asarray(Ns)-occ[:,2]-2*occ[:,1]
        basis += [(s,h0,h1) for h0 in range(bs[0]+1) for h1 in range(bs[1]+1)]
    lookup={v:i for i,v in enumerate(basis)};d=len(basis)
    ds=[np.zeros((d,d),complex) for _ in range(2)]
    us=[np.zeros((d,d),complex) for _ in range(2)]
    edges={k:np.zeros((d,d),complex) for k in WEIGHTS}
    occ=np.array([counts(s) for s,_,_ in basis])
    bs=np.asarray(Ns)[None,:]-occ[:,:,2]-2*occ[:,:,1]
    for col,(s,h0,h1) in enumerate(basis):
        hs=[h0,h1]
        for x in range(2):
            off=4*x;b=int(bs[col,x]);h=hs[x]
            for alpha in range(2):
                occupied=b-h if alpha==0 else h
                moved=hop(s,off+2+alpha,off)
                if moved and occupied:
                    ss,sg=moved;out=hs.copy();out[x]-=(alpha==1)
                    ds[x][lookup[(ss,*out)],col]+=sg*math.sqrt(occupied/b)
                beta=1-alpha;occupied=b-h if beta==0 else h
                moved=hop(s,off+2+alpha,off+1)
                if moved:
                    ss,sg=moved;out=hs.copy();out[x]+=(beta==1)
                    us[x][lookup[(ss,*out)],col]+=sg*(1 if alpha==0 else -1)*math.sqrt((occupied+1)/(b+2))
        for sector,mm in MODES.items():
            weight=WEIGHTS[sector]
            for m in mm:
                moved=hop(s,4+m,m)
                if moved and bs[col,1]-h1>=weight:
                    ss,sg=moved
                    edges[sector][lookup[(ss,h0,h1)],col]+=sg
    e=np.zeros((d,len(lb)),complex)
    ex=[np.zeros_like(e) for _ in range(2)]
    for col,s in enumerate(lb):
        e[lookup[(s,0,0)],col]=1
        ex[0][lookup[(s,1,0)],col]=1;ex[1][lookup[(s,0,1)],col]=1
    return dict(basis=basis,lookup=lookup,D=ds,U=us,edges=edges,E=e,EX=ex,
                occ=occ,bs=bs,leak=np.array([h0+h1>0 for _,h0,h1 in basis]))


def assemble(a,md,mu,jcoef,p):
    d=len(a['D'][0]);v=np.zeros((d,d),complex)
    for x in range(2):
        t=md[x]*a['D'][x]+mu[x]*a['U'][x];v+=t+t.conj().T
    for s in WEIGHTS:
        t=jcoef[s]*a['edges'][s];v+=t+t.conj().T
    occ=a['occ'];diag=sum(p[x]*(occ[:,x,0]+occ[:,x,1]-occ[:,x,2]) for x in range(2))
    return v+np.diag(diag),v,np.diag(diag)


def residual_parts(l,md,mu):
    cs=[];ls=[]
    for x in range(2):
        b=l['bs'][:,x];gap=np.diag(1-np.sqrt((b+1)/(b+2)))
        y=-mu[x]*l['U'][x]
        cs.append(y@gap+gap@y.conj().T)
        ls.append(md[x].conjugate()*l['LD'][x]@np.diag((b+1.)**-.5)
                  +mu[x]*l['LU'][x]@np.diag((b+2.)**-.5))
    return cs,ls


def bound(q,Ns,md,mu):
    return sum((abs(md[x])+abs(mu[x]))/math.sqrt(Ns[x]-2*q+1)
               +2*abs(mu[x])/(Ns[x]-2*q+1) for x in range(2))


def evo(eigen,columns,t):
    w,v=eigen
    return (v*np.exp(-1j*t*w))@(v.conj().T@columns)


def calibration(md0,mu0,kappa,gamma,p):
    q,Ns=1,[3,3];a=physical(q,Ns);l=logical_model(q,Ns)
    h,_,_=assemble(a,md0,mu0,{s:gamma*kappa[s] for s in WEIGHTS},p)
    compressed=a['E'].conj().T@h@a['E'];li=l['lookup']
    md=[];mu=[]
    for x in range(2):
        off=4*x
        md.append(compressed[li[1<<off],li[1<<(off+2)]].conjugate())
        mu.append(-compressed[li[1<<(off+1)],li[1<<(off+3)]].conjugate()/math.sqrt(2/3))
    scale=compressed[li[1<<6],li[1<<2]]/kappa['L']
    err=max(max(abs(md[x]-md0[x]) for x in range(2)),
            max(abs(mu[x]-mu0[x]) for x in range(2)),abs(scale-gamma))
    return md,mu,float(scale.real),dict(q=q,Ns=Ns,recovery_error=float(err),
        gamma=float(scale.real),edge_ratios_phases_are_registered_inputs=True,refits_after_calibration=0)


def history_check(eigen,leigen,e,f,embed,h,v,h0,lh,lv,lh0,physical_tags,logical_tags,occ,locc,bs,q,total_eta):
    t1,t2=.8,1.1
    a=evo(eigen,f,t1);z=evo(leigen,embed,t1)
    blocks=[rows[:,None]*a for rows in physical_tags]
    lblocks=[rows[:,None]*z for rows in logical_tags]
    first_work=-a.conj().T@v@a
    restart_work=sum(k.conj().T@v@k for k in blocks)
    # Propagate the real record branches jointly; no reference or record reset.
    joined=np.column_stack(blocks);lj=np.column_stack(lblocks)
    continued=evo(eigen,joined,t2);lc=evo(leigen,lj,t2)
    m=f.shape[1];finals=np.hsplit(continued,len(blocks));lfin=np.hsplit(lc,len(blocks))
    stacked=np.vstack(finals);target=np.vstack([e@x for x in lfin])
    second_work=-sum(x.conj().T@v@x for x in finals)
    final_energy=sum(x.conj().T@h0@x for x in finals)
    initial_energy=f.conj().T@h@f
    work_res=norm(final_energy-initial_energy-first_work-restart_work-second_work)
    gram=norm(stacked.conj().T@stacked-np.eye(m))
    reference_res=[]
    for x in range(2):
        total=bs[:,x]+occ[:,x,2]+2*occ[:,x,1]
        reference_res.append(norm(sum(y.conj().T@(total[:,None]*y) for y in finals)-f.conj().T@(total[:,None]*f)))
    receive=occ[:,1,:].sum(axis=1)
    rng=np.random.default_rng(1330+q);probe=rng.normal(size=m)+1j*rng.normal(size=m);probe/=np.linalg.norm(probe)
    at_first=float(np.real((a@probe).conj()@(receive*(a@probe))))
    at_second=float(sum(np.real((k@probe).conj()@(receive*(k@probe))) for k in finals))
    error=norm(stacked-target);estimate=min(2.,(t1+t2)*total_eta)
    return dict(first_time=t1,continuation_time=t2,complete_history_error=error,error_bound=estimate,
                receiver_occupation_first=at_first,receiver_occupation_after_record_and_continuation=at_second,
                reference_reset=False,record_reset=False,restart_work_operator_norm=norm(restart_work),
                max_resource_residual=max(work_res,gram,*reference_res),
                passed=bool(error<=estimate+1e-10 and max(work_res,gram,*reference_res)<1e-9))


def run_case(q,Ns,md,mu,jcoef,p):
    a=physical(q,Ns);l=logical_model(q,Ns)
    h,v,h0=assemble(a,md,mu,jcoef,p);lh,lv,lh0=assemble(l,md,mu,jcoef,p)
    e=a['E'];cs,ls=residual_parts(l,md,mu)
    r=h@e-e@lh;pred=e@sum(cs)+sum(a['EX'][x]@ls[x] for x in range(2))
    exact_res=norm(r-pred)
    gram_res=norm(r.conj().T@r-sum(cs).conj().T@sum(cs)-sum(x.conj().T@x for x in ls))
    total_eta=bound(q,Ns,md,mu)
    occ,locc,bs=a['occ'],l['occ'],a['bs']
    qs=occ.sum(axis=2);lqs=locc.sum(axis=2)
    edge_checks=[];currents={};lcurrents={};vedge=np.zeros_like(h)
    for s in WEIGHTS:
        b=jcoef[s]*a['edges'][s];bl=jcoef[s]*l['edges'][s]
        current=1j*(b.conj().T-b);lc=1j*(bl.conj().T-bl)
        currents[s]=current;lcurrents[s]=lc;vedge+=b+b.conj().T
        real_source=a['edges'][s]+a['edges'][s].conj().T
        real_target=l['edges'][s]+l['edges'][s].conj().T
        imag_source=1j*(a['edges'][s]-a['edges'][s].conj().T)
        imag_target=1j*(l['edges'][s]-l['edges'][s].conj().T)
        checks=[norm(a['edges'][s]@e-e@l['edges'][s]),
                norm(a['edges'][s].conj().T@e-e@l['edges'][s].conj().T),
                norm(real_source@e-e@real_target),norm(imag_source@e-e@imag_target),
                norm(current@e-e@lc)]
        edge_checks.append(dict(sector=s,weight=WEIGHTS[s],max_source_intertwiner=max(checks),passed=max(checks)<1e-11))
    itotal=sum(currents.values());litotal=sum(lcurrents.values())
    continuity=frob(1j*h*(qs[:,1][None,:]-qs[:,1][:,None])-itotal)
    ref_edge=frob(1j*vedge*(bs[:,1][None,:]-bs[:,1][:,None])+
                  sum(WEIGHTS[s]*currents[s] for s in WEIGHTS))
    # Retain independent local E_D,E_L,E_U,omega before using the work dictionary.
    independent=np.zeros(len(h));omega=[.4,.3];base=.8
    for x in range(2):
        independent+=(base+p[x])*occ[:,x,0]+(base+2*omega[x]+p[x])*occ[:,x,1]
        independent+=(base+omega[x]-p[x])*occ[:,x,2]+omega[x]*bs[:,x]
    energy_dictionary=float(np.max(np.abs(independent-(base*q+sum(omega[x]*Ns[x] for x in range(2))+np.diag(h0)))))
    number_res=max(norm(bs[:,x,None]*e-e*l['bs'][:,x][None,:]) for x in range(2))
    local_sources=[]
    for x in range(2):
        for sector in ['D','U']:
            smd=[0j,0j];smu=[0j,0j]
            if sector=='D':smd[x]=1j
            else:smu[x]=1.
            source=(1j*a['D'][x]-1j*a['D'][x].conj().T if sector=='D'
                    else a['U'][x]+a['U'][x].conj().T)
            target=(1j*l['D'][x]-1j*l['D'][x].conj().T if sector=='D'
                    else l['U'][x]+l['U'][x].conj().T)
            sc,sl=residual_parts(l,smd,smu)
            sr=source@e-e@target
            err=norm(sr-e@sum(sc)-sum(a['EX'][y]@sl[y] for y in range(2)))
            local_sources.append(dict(site=x,sector=sector,decomposition_residual=err,passed=err<1e-11))
    distribution=[1,0] if q==1 else [1,1]
    cols=[i for i,c in enumerate(locc) if all(c[x,0]==0 and c[x,1]==0 and c[x,2]==distribution[x] for x in range(2))]
    f=e[:,cols];embed=np.eye(len(lh))[:,cols]
    tags=sorted(set(tuple(x) for x in occ.reshape(len(occ),6).tolist()))
    physical_tags=[np.all(occ.reshape(len(occ),6)==tag,axis=1) for tag in tags]
    logical_tags=[np.all(locc.reshape(len(locc),6)==tag,axis=1) for tag in tags]
    record_identity=max(norm(rows[:,None]*e-e*lr[None,:]) for rows,lr in zip(physical_tags,logical_tags))
    eigen=np.linalg.eigh(h);leigen=np.linalg.eigh(lh)
    rng=np.random.default_rng(1330+q);probe=rng.normal(size=len(cols))+1j*rng.normal(size=len(cols));probe/=np.linalg.norm(probe)
    init=f@probe
    nodes,weights=np.polynomial.legendre.leggauss(12)
    holdouts=[]
    for time in [.7,1.8,3.4]:
        actual_all=evo(eigen,e,time);logical_all=evo(leigen,np.eye(len(lh)),time)
        actual=evo(eigen,f,time);logical=evo(leigen,embed,time)
        target=e@logical
        error=norm(actual-target);full_error=norm(actual_all-e@logical_all)
        estimate=min(2.,time*total_eta)
        recorded=np.vstack([rows[:,None]*actual for rows in physical_tags])
        rtarget=np.vstack([e@(rows[:,None]*logical) for rows in logical_tags])
        record_error=norm(recorded-rtarget)
        gram=norm(recorded.conj().T@recorded-np.eye(len(cols)))
        partner=norm(recorded.conj().T@recorded/len(cols)-np.eye(len(cols))/len(cols))
        energy=norm(actual.conj().T@h@actual-f.conj().T@h@f)
        off_work=-actual.conj().T@v@actual
        record_energy=sum((rows[:,None]*actual).conj().T@h0@(rows[:,None]*actual) for rows in physical_tags)
        work=norm(record_energy-f.conj().T@h@f-off_work)
        exact_balance=max(norm(actual.conj().T@((bs[:,x]+occ[:,x,2]+2*occ[:,x,1])[:,None]*actual)-
                               f.conj().T@((bs[:,x]+occ[:,x,2]+2*occ[:,x,1])[:,None]*f)) for x in range(2))
        integral=np.zeros((len(cols),len(cols)),complex)
        for node,weight in zip(nodes,weights):
            state=evo(eigen,f,(node+1)*time/2)
            integral+=weight*time/2*(state.conj().T@itotal@state)
        flow=norm(integral-(actual.conj().T@(qs[:,1,None]*actual)-f.conj().T@(qs[:,1,None]*f)))
        psi=actual@probe;lp=logical@probe
        reception=float(np.real(psi.conj()@(qs[:,1]*psi)))
        predicted_reception=float(np.real(lp.conj()@(lqs[:,1]*lp)))
        reference=[float(np.real(psi.conj()@(bs[:,x]*psi)-init.conj()@(bs[:,x]*init))) for x in range(2)]
        probabilities=[float(np.linalg.norm(psi[rows])**2) for rows in physical_tags]
        edge_values=[dict(sector=s,actual=float(np.real(psi.conj()@currents[s]@psi)),
                          target=float(np.real(lp.conj()@lcurrents[s]@lp))) for s in WEIGHTS]
        checks=[gram,partner,energy,work,exact_balance,flow,abs(record_error-error),abs(sum(probabilities)-1)]
        holdouts.append(dict(time=time,full_code_error=full_error,natural_input_error=error,error_bound=estimate,
             receiver_occupation=reception,target_receiver_occupation=predicted_reception,
             reference_transfers=reference,edge_currents=edge_values,
             transverse_reference_probability=float(np.linalg.norm(psi[a['leak']])**2),
             complete_record_error=record_error,
             records=[dict(site_counts_D_U_L=list(tag),probability=prob) for tag,prob in zip(tags,probabilities)],
             current_integral_identity_residual=flow,max_resource_residual=max(checks),
             passed=bool(max(checks)<1e-9 and full_error<=estimate+1e-10)))
    history=history_check(eigen,leigen,e,f,embed,h,v,h0,lh,lv,lh0,physical_tags,logical_tags,occ,locc,bs,q,total_eta)
    passed=(max(exact_res,gram_res,continuity,ref_edge,energy_dictionary,number_res,record_identity)<1e-10 and
            norm(r)<=total_eta+1e-11 and all(x['passed'] for x in edge_checks+local_sources+holdouts) and history['passed'])
    return dict(q=q,Ns=Ns,parent_dimension=len(h),logical_dimension=len(lh),actual_input_dimension=len(cols),
        initial_local_occupations=distribution,observed_generator_residual=norm(r),generator_bound=total_eta,
        exact_residual_decomposition=exact_res,orthogonal_gram_residual=gram_res,
        occupation_continuity_residual=continuity,edge_reference_current_residual=ref_edge,
        energy_dictionary_residual=energy_dictionary,reference_source_identity_residual=number_res,
        record_intertwiner_residual=record_identity,edge_sources=edge_checks,local_sources=local_sources,
        holdouts=holdouts,recorded_continuation=history,passed=bool(passed))


def compute():
    md0=[.16+.015j,.16+.015j];mu0=[.09-.01j,.09-.01j];p=[.025,-.015]
    kappa={'D':.9*np.exp(.2j),'L':1.1*np.exp(-.15j),'U':.7*np.exp(.35j)}
    md,mu,gamma,cal=calibration(md0,mu0,kappa,.15,p)
    jc={s:gamma*kappa[s] for s in WEIGHTS}
    cases=[run_case(q,Ns,md,mu,jc,p) for q,Ns in [(1,[3,3]),(1,[5,4]),(2,[5,5])]]
    paths=[HERE/'proof.md',HERE/'check.py',HERE.parent/'_admission'/'round1330_admission.md',
           HERE.parent/'1329'/'proof.md']
    hashes={str(p.relative_to(HERE.parent)).replace('\\','/'):
            hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    passed=cal['recovery_error']<1e-12 and all(x['passed'] for x in cases)
    return dict(round=1330,status='PASS' if passed else 'FAIL',new_axioms=0,provenance=hashes,
        calibration=cal,cases=cases,full_reference_boundaries_retained=True,
        local_K_is_Gauss=False,local_SU2_certified=False,empirical_experiments=0,
        actual_inputs='known local left occupations, arbitrary internal state and partner',
        arbitrary_spatial_code_preparation_certified=False,stage_complete=False,long_goal_complete=False)


def compare(a,b,path='root'):
    if isinstance(a,dict):
        assert set(a)==set(b),(path,'keys')
        for k in a:compare(a[k],b[k],path+'.'+k)
    elif isinstance(a,list):
        assert len(a)==len(b),(path,'length')
        for i,(x,y) in enumerate(zip(a,b)):compare(x,y,path+f'[{i}]')
    elif isinstance(a,float):
        assert math.isclose(a,b,rel_tol=1e-8,abs_tol=1e-9),(path,a,b)
    else:assert a==b,(path,a,b)


def main():
    parser=argparse.ArgumentParser();group=parser.add_mutually_exclusive_group()
    group.add_argument('--check',action='store_true');group.add_argument('--write-new',action='store_true')
    args=parser.parse_args();r=compute();assert r['status']=='PASS',r
    destination=HERE/'results.json'
    if args.write_new:
        payload=json.dumps(r,ensure_ascii=False,indent=2)+'\n'
        with destination.open('x',encoding='utf-8') as f:
            f.write(payload)
    else:compare(r,json.loads(destination.read_text(encoding='utf-8')))
    print(json.dumps(dict(status='PASS',round=1330,complete_two_site_parents=3,time_holdouts=9,
                         recorded_continuations=3,results_written=args.write_new),ensure_ascii=False))


if __name__=='__main__':main()
