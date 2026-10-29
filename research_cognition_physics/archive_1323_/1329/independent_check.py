"""Independent 1329 full finite CAR parent audit; no author implementation imports.

Actual mode order D,L0,L1,U. The U target has its fixed epsilon sign (-M_U).
--write-new creates output exclusively; default/--check compares read-only.
"""
from __future__ import annotations
import argparse
import itertools
import json
import math
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent


def hop(occupied,destination,origin):
    if origin not in occupied:
        return None
    location=occupied.index(origin)
    rest=occupied[:location]+occupied[location+1:]
    if destination in rest:
        return None
    insertion=sum(x<destination for x in rest)
    return rest[:insertion]+(destination,)+rest[insertion:], (-1)**(location+insertion)


def counts(occupied,n):
    nd=sum(x<n for x in occupied)
    nl=sum(n<=x<3*n for x in occupied)
    nu=sum(x>=3*n for x in occupied)
    return nd,nl,nu


def build(n,q,N,MD,MU,p=0.0):
    fermions=list(itertools.combinations(range(4*n),q))
    states=[]
    for occupied in fermions:
        nd,nl,nu=counts(occupied,n)
        remaining=N-nl-2*nu
        if remaining>=0:
            states += [(occupied,blue) for blue in range(remaining+1)]
    index={x:i for i,x in enumerate(states)}
    logical_index={x:i for i,x in enumerate(fermions)}
    T=np.zeros((len(states),len(states)),complex)
    logical_T=np.zeros((len(fermions),len(fermions)),complex)
    nds,nls,nus,nb=[],[],[],[]
    for column,(occupied,blue) in enumerate(states):
        nd,nl,nu=counts(occupied,n)
        balance=N-nl-2*nu
        red=balance-blue
        nds.append(nd);nls.append(nl);nus.append(nu);nb.append(balance)
        for i in range(n):
            for j in range(n):
                if red:
                    moved=hop(occupied,n+i,j)
                    if moved:
                        out,sign=moved
                        T[index[(out,blue)],column]+=sign*MD[i,j]*math.sqrt(red/balance)
                if blue:
                    moved=hop(occupied,2*n+i,j)
                    if moved:
                        out,sign=moved
                        T[index[(out,blue-1)],column]+=sign*MD[i,j]*math.sqrt(blue/balance)
                moved=hop(occupied,n+i,3*n+j)
                if moved:
                    out,sign=moved
                    T[index[(out,blue+1)],column]+=sign*MU[i,j]*math.sqrt((blue+1)/(balance+2))
                moved=hop(occupied,2*n+i,3*n+j)
                if moved:
                    out,sign=moved
                    T[index[(out,blue)],column]-=sign*MU[i,j]*math.sqrt((red+1)/(balance+2))
    for column,occupied in enumerate(fermions):
        for i in range(n):
            for j in range(n):
                for destination,origin,coefficient in [(n+i,j,MD[i,j]),(2*n+i,3*n+j,-MU[i,j])]:
                    moved=hop(occupied,destination,origin)
                    if moved:
                        out,sign=moved
                        logical_T[logical_index[out],column]+=sign*coefficient
    E=np.zeros((len(states),len(fermions)),complex)
    for i,occupied in enumerate(fermions):
        if (occupied,0) in index:
            E[index[(occupied,0)],i]=1
    PD,PL,PU,NB=[np.diag(x) for x in (nds,nls,nus,nb)]
    ld,ll,lu=[np.diag([counts(s,n)[i] for s in fermions]) for i in range(3)]
    V=T+T.conj().T
    lv=logical_T+logical_T.conj().T
    return dict(states=states,fermions=fermions,index=index,E=E,T=T,V=V,
                H=V+p*(PD+PU-PL),logical_H=lv+p*(ld+lu-ll),logical_V=lv,
                ND=PD,NL=PL,NU=PU,NB=NB,ld=ld,ll=ll,lu=lu)


def operator_norm(a):
    return float(np.linalg.norm(a,ord=2))


def nuclear_norm(a):
    return float(np.linalg.svd(a,compute_uv=False).sum())


def propagator(H,time):
    values,vectors=np.linalg.eigh(H)
    return (vectors*np.exp(-1j*time*values))@vectors.conj().T


def compute():
    checks=[]
    cases=[]

    def eq(name,a,b,tolerance=5e-10):
        residual=float(np.linalg.norm(np.asarray(a)-np.asarray(b)))
        checks.append(dict(name=name,residual=residual,tolerance=tolerance))
        if residual>tolerance:raise AssertionError((name,residual,tolerance))

    def bound(name,value,ceiling,tolerance=5e-10):
        residual=max(0.0,float(value-ceiling))
        checks.append(dict(name=name,residual=residual,tolerance=tolerance,value=float(value),bound=float(ceiling)))
        if residual>tolerance:raise AssertionError((name,value,ceiling))

    original_D=np.array([[.31+.04j,.07-.03j],[-.02+.05j,.22-.01j]])
    original_U=np.array([[.13-.02j,.05+.06j],[.08-.04j,.26+.03j]])
    n,N0=2,3
    calibration=build(n,1,N0,original_D,original_U)
    compressed=calibration['E'].conj().T@calibration['H']@calibration['E']
    MD=compressed[n:2*n,:n]
    MU=-compressed[2*n:3*n,3*n:4*n]/math.sqrt((N0-1)/N0)
    eq('one_calibration.D',MD,original_D)
    eq('one_calibration.U',MU,original_U)
    assert operator_norm(MD@MU-MU@MD)>0.005
    configuration=[(1,1,3),(1,2,5),(2,1,3),(2,2,5),(2,2,9),(2,2,17),(2,3,7),(2,4,9)]
    for n,q,N in configuration:
        A=MD[:n,:n];B=MU[:n,:n]
        p=.11
        a=build(n,q,N,A,B,p)
        E,H,h=a['E'],a['H'],a['logical_H']
        dim=len(a['states']);ldim=len(a['fermions'])
        minimum=N-2*q
        assert minimum>=1
        label=f'n{n}.q{q}.N{N}'
        sD,sU=nuclear_norm(A),nuclear_norm(B)
        declared=(sD+sU)/math.sqrt(minimum+1)+2*sU/(minimum+1)
        delta=H@E-E@h
        leakage=delta-E@(E.conj().T@delta)
        code_defect=E.conj().T@delta
        eq(label+'.isometry',E.conj().T@E,np.eye(ldim))
        bound(label+'.generator_bound',operator_norm(delta),declared)
        bound(label+'.leakage_bound',operator_norm(leakage),(sD+sU)/math.sqrt(minimum+1))
        bound(label+'.code_defect_bound',operator_norm(code_defect),2*sU/(minimum+1))
        eq(label+'.full_fermion_count',a['ND']+a['NL']+a['NU'],q*np.eye(dim))
        eq(label+'.full_K',a['NB']+a['NL']+2*a['NU'],N*np.eye(dim))
        for key,logical in [('ND',a['ld']),('NL',a['ll']),('NU',a['lu'])]:
            eq(label+'.source_'+key,a[key]@E,E@logical)
        eq(label+'.reference_source',a['NB']@E,E@(N*np.eye(ldim)-a['ll']-2*a['lu']))
        eq(label+'.p_source',(a['ND']+a['NU']-a['NL'])@E,E@(a['ld']+a['lu']-a['ll']))
        # Every real and imaginary independent portal parameter is tested, not only a scalar rate.
        zero=np.zeros((n,n),complex)
        selected_source=None
        for department in ('D','U'):
            for i in range(n):
                for j in range(n):
                    for component,z in [('real',1.0),('imag',1j)]:
                        differential=np.zeros((n,n),complex);differential[i,j]=z
                        dA=differential if department=='D' else zero
                        dB=differential if department=='U' else zero
                        source=build(n,q,N,dA,dB)
                        eta=(nuclear_norm(dA)+nuclear_norm(dB))/math.sqrt(minimum+1)+2*nuclear_norm(dB)/(minimum+1)
                        error=operator_norm(source['V']@E-E@source['logical_V'])
                        bound(label+f'.Msource.{department}{i}{j}.{component}',error,eta)
                        if department=='U' and i==0 and j==n-1 and component=='imag':
                            selected_source=(source['V'],source['logical_V'],error)
        # No state outside E is projected away in the parent evolution.
        eq(label+'.hermitian_parent',H,H.conj().T)
        allowed_left=[j for j,occ in enumerate(a['fermions']) if all(n<=m<3*n for m in occ)]
        assert len(allowed_left)==math.comb(2*n,q)
        physical_left=E[:,allowed_left]
        logical_left=np.eye(ldim)[:,allowed_left]
        eq(label+'.actual_all_left_reference',a['NB']@physical_left,(N-q)*physical_left)
        eq(label+'.actual_all_left_unknown',physical_left.conj().T@physical_left,np.eye(len(allowed_left)))
        E0,omega=2.7,.63
        H0=(E0+p)*a['ND']+(E0+omega-p)*a['NL']+(E0+2*omega+p)*a['NU']+omega*a['NB']
        Hfull=H0+a['V']
        eq(label+'.all_independent_energy_sources',Hfull,(E0*q+omega*N)*np.eye(dim)+H)
        # Occupancy triples are an actual complete record, including every failure and leakage.
        triples=sorted(set(zip(np.diag(a['ND']).astype(int),np.diag(a['NL']).astype(int),np.diag(a['NU']).astype(int))))
        physical_masks=[];logical_masks=[]
        for triple in triples:
            pm=(np.diag(a['ND'])==triple[0])&(np.diag(a['NL'])==triple[1])&(np.diag(a['NU'])==triple[2])
            lm=(np.diag(a['ld'])==triple[0])&(np.diag(a['ll'])==triple[1])&(np.diag(a['lu'])==triple[2])
            physical_masks.append(pm);logical_masks.append(lm)
            eq(label+'.record_dictionary_'+str(triple),pm[:,None]*E,E*lm[None,:])
        time_rows=[]
        for time in (.29,1.3,3.1):
            U,u=propagator(H,time),propagator(h,time)
            full_error=operator_norm(U@E-E@u)
            bound(label+f'.t{time}.full_process',full_error,min(2,time*declared))
            bound(label+f'.t{time}.measured_residual_duhamel',full_error,min(2,time*operator_norm(delta)))
            actual=U@physical_left
            ideal=E@u@logical_left
            allowed_error=operator_norm(actual-ideal)
            bound(label+f'.t{time}.qualified_all_left',allowed_error,min(2,time*declared))
            source,target_source,source_error=selected_source
            source_observable=(U@E).conj().T@source@(U@E)-u.conj().T@target_source@u
            bound(label+f'.t{time}.evolved_source',operator_norm(source_observable),2*operator_norm(source)*full_error+source_error)
            # Full maps cover arbitrary entangled partners; the normalized identity sample is maximally entangled.
            record_actual=np.concatenate([m[:,None]*actual for m in physical_masks],axis=0)
            record_ideal=np.concatenate([m[:,None]*ideal for m in physical_masks],axis=0)
            eq(label+f'.t{time}.full_record',record_actual.conj().T@record_actual,np.eye(len(allowed_left)))
            eq(label+f'.t{time}.record_same_error',operator_norm(record_actual-record_ideal),allowed_error)
            out=actual/math.sqrt(len(allowed_left))
            initial=physical_left/math.sqrt(len(allowed_left))
            work_off=-float(np.vdot(out,a['V']@out).real)
            energy_record=sum(np.vdot(m[:,None]*out,H0@(m[:,None]*out)).real for m in physical_masks)
            eq(label+f'.t{time}.quench_and_records',energy_record-np.vdot(out,Hfull@out).real,work_off)
            eq(label+f'.t{time}.energy_conservation',np.vdot(out,Hfull@out),np.vdot(initial,Hfull@initial))
            change_reference=np.vdot(out,a['NB']@out)-np.vdot(initial,a['NB']@initial)
            change_left=np.vdot(out,a['NL']@out)-np.vdot(initial,a['NL']@initial)
            change_up=np.vdot(out,a['NU']@out)-np.vdot(initial,a['NU']@initial)
            eq(label+f'.t{time}.real_reference_flow',change_reference+change_left+2*change_up,0)
            leaked=actual-E@(E.conj().T@actual)
            leakage_probability=float(np.linalg.norm(leaked)**2/len(allowed_left))
            time_rows.append(dict(time=time,whole_code_error=full_error,qualified_left_error=allowed_error,
                                  actual_left_leakage=leakage_probability,reference_transfer=float(change_reference.real),
                                  work_off=work_off,bound=min(2,time*declared)))
        cases.append(dict(n=n,q=q,N=N,bmin=minimum,physical_dimension=dim,logical_dimension=ldim,
                          actual_left_input_dimension=len(allowed_left),generator_error=operator_norm(delta),
                          declared_bound=declared,code_defect=operator_norm(code_defect),leakage_norm=operator_norm(leakage),times=time_rows))
    return dict(round=1329,status='PASS',new_axioms=0,stage_complete=False,
                construction='Complete tuple-CAR parent, direct zero-mode reference embedding, fixed negative U target convention; no author imports.',
                calibration=dict(N0=N0,q0=1,subsequent_refits=0),check_count=len(checks),
                max_residual=max(c['residual'] for c in checks),cases=cases,checks=checks)


def compare(a,b,path='root'):
    if isinstance(b,dict):
        assert a.keys()==b.keys(),path
        for key in b:compare(a[key],b[key],path+'.'+key)
    elif isinstance(b,list):
        assert len(a)==len(b),path
        for i,(x,y) in enumerate(zip(a,b)):compare(x,y,path+f'[{i}]')
    elif isinstance(b,float):
        assert math.isfinite(a) and abs(a-b)<=1e-9*(1+abs(b)),path
    else:assert a==b,path


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-new',action='store_true')
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args();result=compute()
    output=HERE/'independent_results.json'
    if args.write_new:
        with output.open('x',encoding='utf-8') as handle:
            json.dump(result,handle,ensure_ascii=False,indent=2);handle.write('\n')
    else:compare(json.loads(output.read_text(encoding='utf-8')),result)
    print(json.dumps({key:result[key] for key in ('round','status','check_count','max_residual')},ensure_ascii=False))


if __name__=='__main__':main()
