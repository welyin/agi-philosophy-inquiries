"""Independent 1330 full two-cell CAR/reference parent audit (no author imports).

Actual mode order in each cell: D,L0,L1,U. A common nonzero edge matrix for all
roles is a declared sample input, not a derived universal coupling law.
Default/--check compares read-only; --write-new creates the output exclusively.
"""
from __future__ import annotations
import argparse
import itertools
import json
import math
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent


def move(occupation,destination,origin):
    if origin not in occupation:return None
    position=occupation.index(origin)
    remainder=occupation[:position]+occupation[position+1:]
    if destination in remainder:return None
    insertion=sum(i<destination for i in remainder)
    return remainder[:insertion]+(destination,)+remainder[insertion:],(-1)**(position+insertion)


def counts(occupation,n):
    output=[]
    for x in range(2):
        local=[i-x*4*n for i in occupation if x*4*n<=i<(x+1)*4*n]
        output.append((sum(i<n for i in local),sum(n<=i<3*n for i in local),sum(i>=3*n for i in local)))
    return output


def build(n,q,N,MD,MU,J,p):
    logical=list(itertools.combinations(range(8*n),q))
    states=[]
    for occupation in logical:
        numbers=counts(occupation,n)
        balance=[N[x]-numbers[x][1]-2*numbers[x][2] for x in range(2)]
        states += [(occupation,h0,h1) for h0 in range(balance[0]+1) for h1 in range(balance[1]+1)]
    index={s:i for i,s in enumerate(states)};li={s:i for i,s in enumerate(logical)}
    size,lsize=len(states),len(logical)
    on=np.zeros((size,size),complex);edge=np.zeros_like(on)
    lon=np.zeros((lsize,lsize),complex);ledge=np.zeros_like(lon)
    physical_counts=np.array([counts(s[0],n) for s in states])
    logical_counts=np.array([counts(s,n) for s in logical])
    references=np.array([[N[x]-c[x,1]-2*c[x,2] for x in range(2)] for c in physical_counts])
    for column,(occupation,h0,h1) in enumerate(states):
        blue=(h0,h1)
        for x in range(2):
            b=references[column,x];red=b-blue[x];offset=x*4*n
            for i in range(n):
                for j in range(n):
                    transitions=[(offset+n+i,offset+j,MD[x][i,j]*math.sqrt(red/b),0),
                                 (offset+2*n+i,offset+j,MD[x][i,j]*math.sqrt(blue[x]/b),-1),
                                 (offset+n+i,offset+3*n+j,MU[x][i,j]*math.sqrt((blue[x]+1)/(b+2)),1),
                                 (offset+2*n+i,offset+3*n+j,-MU[x][i,j]*math.sqrt((red+1)/(b+2)),0)]
                    for destination,origin,coefficient,change in transitions:
                        if coefficient==0:continue
                        result=move(occupation,destination,origin)
                        if result:
                            out,sign=result;updated=list(blue);updated[x]+=change
                            on[index[(out,*updated)],column]+=sign*coefficient
        # Actual y=1 -> x=0 edge, with reverse direction supplied by the adjoint.
        for role,weight in enumerate((0,1,1,2)):
            if references[column,0]-h0<weight:continue
            for i in range(n):
                for j in range(n):
                    result=move(occupation,role*n+i,4*n+role*n+j)
                    if result:
                        out,sign=result
                        edge[index[(out,h0,h1)],column]+=sign*J[i,j]
    for column,occupation in enumerate(logical):
        for x in range(2):
            offset=x*4*n
            for i in range(n):
                for j in range(n):
                    for destination,origin,coefficient in [(offset+n+i,offset+j,MD[x][i,j]),(offset+2*n+i,offset+3*n+j,-MU[x][i,j])]:
                        result=move(occupation,destination,origin)
                        if result:
                            out,sign=result;lon[li[out],column]+=sign*coefficient
        for role in range(4):
            for i in range(n):
                for j in range(n):
                    result=move(occupation,role*n+i,4*n+role*n+j)
                    if result:
                        out,sign=result;ledge[li[out],column]+=sign*J[i,j]
    E=np.zeros((size,lsize),complex)
    for column,occupation in enumerate(logical):E[index[(occupation,0,0)],column]=1
    H0=np.diag(sum(p[x]*(physical_counts[:,x,0]+physical_counts[:,x,2]-physical_counts[:,x,1]) for x in range(2)))
    h0=np.diag(sum(p[x]*(logical_counts[:,x,0]+logical_counts[:,x,2]-logical_counts[:,x,1]) for x in range(2)))
    V=on+on.conj().T+edge+edge.conj().T
    return dict(states=states,logical=logical,index=index,E=E,on=on,edge=edge,H=H0+V,V=V,H0=H0,
                h=h0+lon+lon.conj().T+ledge+ledge.conj().T,logical_edge=ledge,logical_on=lon,
                counts=physical_counts,logical_counts=logical_counts,references=references)


def opnorm(a):return float(np.linalg.norm(a,ord=2))
def nuclear(a):return float(np.linalg.svd(a,compute_uv=False).sum())


def evolve_columns(eigen,columns,time):
    values,vectors=eigen
    return (vectors*np.exp(-1j*time*values))@(vectors.conj().T@columns)


def compute():
    checks=[];cases=[]
    def eq(name,a,b,tolerance=8e-10):
        residual=float(np.linalg.norm(np.asarray(a)-np.asarray(b)))
        checks.append(dict(name=name,residual=residual,tolerance=tolerance))
        if residual>tolerance:raise AssertionError((name,residual,tolerance))
    def bound(name,value,limit,tolerance=8e-10):
        residual=max(0.0,float(value-limit))
        checks.append(dict(name=name,residual=residual,tolerance=tolerance,value=float(value),bound=float(limit)))
        if residual>tolerance:raise AssertionError((name,value,limit))

    original_D=[np.array([[.17+.02j,.03-.01j],[.02j,.13]]),np.array([[.11-.01j,.02j],[.015,.19+.02j]])]
    original_U=[np.array([[.09,.025j],[-.01+.02j,.14]]),np.array([[.16,.01-.02j],[.03j,.10-.01j]])]
    original_J=np.array([[.21+.04j,.03-.015j],[.02j,.18-.03j]])
    zero2=np.zeros((2,2),complex)
    calibration=build(2,1,(3,3),original_D,original_U,original_J,(.07,-.03))
    compression=calibration['E'].conj().T@calibration['H']@calibration['E']
    MD=[];MU=[]
    for x in range(2):
        o=x*8
        MD.append(compression[o+2:o+4,o:o+2])
        MU.append(-compression[o+4:o+6,o+6:o+8]/math.sqrt(2/3))
        eq('calibration.MD'+str(x),MD[-1],original_D[x])
        eq('calibration.MU'+str(x),MU[-1],original_U[x])
    J=compression[2:4,10:12]
    eq('calibration.common_edge_from_left',J,original_J)
    for n,q,N in [(1,1,(3,4)),(1,2,(5,5)),(2,1,(3,3))]:
        A=[x[:n,:n] for x in MD];B=[x[:n,:n] for x in MU];G=J[:n,:n]
        a=build(n,q,N,A,B,G,(.07,-.03));E,H,h=a['E'],a['H'],a['h']
        dim,ldim=E.shape;tag=f'n{n}.q{q}.N{N[0]}_{N[1]}'
        eta=sum((nuclear(A[x])+nuclear(B[x]))/math.sqrt(N[x]-2*q+1)+2*nuclear(B[x])/(N[x]-2*q+1) for x in range(2))
        eq(tag+'.isometry',E.conj().T@E,np.eye(ldim))
        eq(tag+'.exact_edge_forward',a['edge']@E,E@a['logical_edge'])
        eq(tag+'.exact_edge_backward',a['edge'].conj().T@E,E@a['logical_edge'].conj().T)
        residual=H@E-E@h
        on_residual=(a['on']+a['on'].conj().T)@E-E@(a['logical_on']+a['logical_on'].conj().T)
        eq(tag+'.shared_residual_is_onsite_only',residual,on_residual)
        bound(tag+'.full_generator_bound',opnorm(residual),eta)
        for x in range(2):
            K=a['references'][:,x]+a['counts'][:,x,1]+2*a['counts'][:,x,2]
            eq(tag+f'.full_local_K{x}',K,N[x]*np.ones(dim))
            for category in range(3):
                eq(tag+f'.count_source{x}_{category}',a['counts'][:,x,category,None]*E,E*a['logical_counts'][:,x,category][None,:])
            qx=np.diag(a['counts'][:,x,:].sum(axis=1))
            current=1j*(H@qx-qx@H)
            phase_source=1j*(a['edge']-a['edge'].conj().T)
            eq(tag+f'.edge_continuity{x}',current,(-1 if x==0 else 1)*phase_source)
        # The phase source and every real/imaginary edge matrix direction are independent sources.
        zero=np.zeros((n,n),complex)
        for i in range(n):
            for j in range(n):
                for part,z in [('real',1.),('imag',1j)]:
                    dG=np.zeros((n,n),complex);dG[i,j]=z
                    source=build(n,q,N,[zero,zero],[zero,zero],dG,(0.,0.))
                    eq(tag+f'.edge_source{i}{j}.{part}',source['V']@E,E@source['h'])
        # Full-boundary witness: a direction-1-only destination cannot supply the a0 shift.
        blocked=0
        for col,(occupation,h0,h1) in enumerate(a['states']):
            for role,weight in enumerate((0,1,1,2)):
                if weight and a['references'][col,0]-h0<weight:
                    for i in range(n):
                        for j in range(n):
                            out=move(occupation,role*n+i,4*n+role*n+j)
                            if out:
                                dest,sign=out
                                if (dest,h0,h1) not in a['index']:
                                    blocked+=1
        assert blocked>0
        # Known local left particle numbers; all internal directions/flavors remain unknown.
        local_q=(1,0) if q==1 else (1,1)
        selected=[j for j,s in enumerate(a['logical']) if all(c[0]==0 and c[2]==0 and c[1]==local_q[x] for x,c in enumerate(counts(s,n)))]
        F=E[:,selected];lf=np.eye(ldim)[:,selected]
        assert len(selected)==math.prod(math.comb(2*n,local_q[x]) for x in range(2))
        for x in range(2):eq(tag+f'.actual_prepared_reference{x}',a['references'][:,x,None]*F,(N[x]-local_q[x])*F)
        occupation_tags=sorted(set(tuple(row.flatten()) for row in a['counts']))
        pmasks=[np.all(a['counts'].reshape(dim,6)==c,axis=1) for c in occupation_tags]
        lmasks=[np.all(a['logical_counts'].reshape(ldim,6)==c,axis=1) for c in occupation_tags]
        for i,(pm,lm) in enumerate(zip(pmasks,lmasks)):eq(tag+f'.record_dictionary{i}',pm[:,None]*E,E*lm[None,:])
        heigen,leigen=np.linalg.eigh(H),np.linalg.eigh(h)
        m=len(selected);time_rows=[]
        for time in (.43,1.2,2.3):
            actual=evolve_columns(heigen,F,time);ideal=E@evolve_columns(leigen,lf,time)
            error=opnorm(actual-ideal)
            bound(tag+f'.t{time}.complete_unknown_input',error,min(2,time*eta))
            branch=np.concatenate([mask[:,None]*actual for mask in pmasks],axis=0)
            ideal_branch=np.concatenate([mask[:,None]*ideal for mask in pmasks],axis=0)
            eq(tag+f'.t{time}.record_norm',branch.conj().T@branch,np.eye(m))
            eq(tag+f'.t{time}.record_same_error',opnorm(branch-ideal_branch),error)
            off_work=-actual.conj().T@a['V']@actual
            postenergy=sum((mask[:,None]*actual).conj().T@a['H0']@(mask[:,None]*actual) for mask in pmasks)
            eq(tag+f'.t{time}.record_and_quench',postenergy-F.conj().T@H@F,off_work)
            refchanges=[]
            for x in range(2):
                dref=actual.conj().T@(a['references'][:,x,None]*actual)-F.conj().T@(a['references'][:,x,None]*F)
                weight=a['counts'][:,x,1]+2*a['counts'][:,x,2]
                dweight=actual.conj().T@(weight[:,None]*actual)-F.conj().T@(weight[:,None]*F)
                eq(tag+f'.t{time}.local_reference_balance{x}',dref+dweight,np.zeros((m,m)))
                refchanges.append(float(np.trace(dref).real/m))
            receiver=float(np.trace(actual.conj().T@(a['counts'][:,1,:].sum(axis=1)[:,None]*actual)).real/m)
            # Fresh inputs were not supplied after the first segment.
            middle=evolve_columns(heigen,F,time/2)
            continued=evolve_columns(heigen,middle,time/2)
            eq(tag+f'.t{time}.actual_successor_without_reset',continued,actual)
            time_rows.append(dict(time=time,actual_input_error=error,receiver_mean_count=receiver,
                                  reference_changes=refchanges,mean_quench_work=float(np.trace(off_work).real/m)))
        # A true intermediate record followed by further simultaneous onsite/edge dynamics.
        t1,t2=.43,.71
        first=evolve_columns(heigen,F,t1);lfirst=evolve_columns(leigen,lf,t1)
        second=[];desired=[]
        for pm,lm in zip(pmasks,lmasks):
            second.append(evolve_columns(heigen,pm[:,None]*first,t2))
            desired.append(E@evolve_columns(leigen,lm[:,None]*lfirst,t2))
        second=np.concatenate(second,axis=0);desired=np.concatenate(desired,axis=0)
        eq(tag+'.record_successor_norm',second.conj().T@second,np.eye(m))
        bound(tag+'.record_successor_error',opnorm(second-desired),min(2,(t1+t2)*eta))
        # All V terms change at least one recorded occupancy, so reopening after the record costs zero in expectation.
        reopened=sum((pm[:,None]*first).conj().T@a['V']@(pm[:,None]*first) for pm in pmasks)
        eq(tag+'.reopening_control_work',reopened,np.zeros((m,m)))
        cases.append(dict(n=n,q=q,N=list(N),physical_dimension=dim,logical_dimension=ldim,
                          actual_input_dimension=m,boundary_blocked_transitions=blocked,generator_error=opnorm(residual),
                          generator_bound=eta,onsite_and_edge_nonzero=True,times=time_rows,
                          recorded_successor_error=opnorm(second-desired)))
    return dict(round=1330,status='PASS',new_axioms=0,stage_complete=False,
                construction='Complete tuple-CAR basis with both complete reference directions, actual a0 unit shifts and their adjoints; no author imports.',
                calibration=dict(q0=1,N0=[3,3],later_refits=0,edge_matrix_from_left_hop=True),
                check_count=len(checks),max_residual=max(x['residual'] for x in checks),cases=cases,checks=checks)


def compare(a,b,path='root'):
    if isinstance(b,dict):
        assert a.keys()==b.keys(),path
        for key in b:compare(a[key],b[key],path+'.'+key)
    elif isinstance(b,list):
        assert len(a)==len(b),path
        for i,(x,y) in enumerate(zip(a,b)):compare(x,y,path+f'[{i}]')
    elif isinstance(b,float):assert math.isfinite(a) and abs(a-b)<=1e-9*(1+abs(b)),path
    else:assert a==b,path


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-new',action='store_true');parser.add_argument('--check',action='store_true')
    args=parser.parse_args();result=compute();destination=HERE/'independent_results.json'
    if args.write_new:
        with destination.open('x',encoding='utf-8') as handle:json.dump(result,handle,ensure_ascii=False,indent=2);handle.write('\n')
    else:compare(json.loads(destination.read_text(encoding='utf-8')),result)
    print(json.dumps({key:result[key] for key in ('round','status','check_count','max_residual','cases')},ensure_ascii=False))


if __name__=='__main__':main()
