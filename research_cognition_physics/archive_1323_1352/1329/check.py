"""1329: two opposite reference ports and a complete many-fermion parent.

NumPy only. Default/--check compares with saved results read-only.
--write-new creates results.json exclusively. No favorable reference rows are discarded.
"""
from __future__ import annotations
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent


def norm(a):
    return float(np.linalg.norm(a, 2))


def nuclear(a):
    return float(np.linalg.svd(a, compute_uv=False).sum())


def masks(n, q):
    return [sum(1 << i for i in c) for c in itertools.combinations(range(n), q)]


def hop(s, i, j):
    if not (s >> j) & 1:
        return None
    sign = (-1) ** ((s & ((1 << j)-1)).bit_count())
    t = s ^ (1 << j)
    if (t >> i) & 1:
        return None
    sign *= (-1) ** ((t & ((1 << i)-1)).bit_count())
    return t | (1 << i), sign


def counts(s, n):
    bits = (1 << n)-1
    nd = (s & bits).bit_count()
    nu = ((s >> n) & bits).bit_count()
    nl = (s >> (2*n)).bit_count()
    return nd, nu, nl


def physical(n, q, N, md, mu, p):
    logical = masks(4*n, q)
    basis = []
    for s in logical:
        _, nu, nl = counts(s, n)
        b = N-nl-2*nu
        basis += [(s, h) for h in range(b+1)]
    lookup = {x:i for i,x in enumerate(basis)}
    size = len(basis)
    td = np.zeros((size,size), complex)
    tu = np.zeros_like(td)
    occ = np.array([counts(s,n) for s,_ in basis])
    nb = N-occ[:,2]-2*occ[:,1]
    for col,(s,h) in enumerate(basis):
        b = int(nb[col])
        for alpha in range(2):
            occupied = b-h if alpha==0 else h
            beta = 1-alpha
            created = (b-h if beta==0 else h)+1
            for i in range(n):
                for j in range(n):
                    a = hop(s,(2+alpha)*n+i,j)
                    if a and occupied:
                        ss,sg = a
                        td[lookup[(ss,h-(alpha==1))],col] += sg*md[i,j]*math.sqrt(occupied/b)
                    a = hop(s,(2+alpha)*n+i,n+j)
                    if a:
                        ss,sg = a
                        tu[lookup[(ss,h+(beta==1))],col] += (
                            sg*(1 if alpha==0 else -1)*mu[i,j]*math.sqrt(created/(b+2)))
    v = td+tu+td.conj().T+tu.conj().T
    h0 = np.diag(p*(occ[:,0]+occ[:,1]-occ[:,2]))
    e0 = np.zeros((size,len(logical)),complex)
    e1 = np.zeros_like(e0)
    for col,s in enumerate(logical):
        e0[lookup[(s,0)],col]=1
        e1[lookup[(s,1)],col]=1
    return dict(basis=basis,lookup=lookup,logical=logical,E0=e0,E1=e1,
                H=v+h0,V=v,H0=h0,occupations=occ,nb=nb,
                leak_rows=np.array([h!=0 for _,h in basis]))


def target(n,q,N,md,mu,p):
    basis=masks(4*n,q);lookup={s:i for i,s in enumerate(basis)}
    size=len(basis)
    x=np.zeros((size,size),complex);y=np.zeros_like(x)
    ld=np.zeros_like(x);lu=np.zeros_like(x)
    for col,s in enumerate(basis):
        for i in range(n):
            for j in range(n):
                for matrix,a,b,coefficient in [
                    (x,2*n+i,j,md[i,j]),(y,3*n+i,n+j,mu[i,j]),
                    (ld,j,3*n+i,md[i,j].conjugate()),(lu,2*n+i,n+j,mu[i,j])]:
                    moved=hop(s,a,b)
                    if moved:
                        ss,sg=moved;matrix[lookup[ss],col]+=sg*coefficient
    occ=np.array([counts(s,n) for s in basis])
    b=N-occ[:,2]-2*occ[:,1]
    f=np.diag(np.sqrt((b+1)/(b+2)))
    gap=np.eye(size)-f
    compression=y@gap+gap@y.conj().T
    leak=ld@np.diag((b+1.)**-.5)+lu@np.diag((b+2.)**-.5)
    h=x+x.conj().T-y-y.conj().T+np.diag(p*(occ[:,0]+occ[:,1]-occ[:,2]))
    return dict(basis=basis,H=h,C=compression,L=leak,occupations=occ,b=b)


def evolve_columns(eigen,columns,t):
    w,v=eigen
    return (v*np.exp(-1j*t*w))@(v.conj().T@columns)


def evolve(h,t):
    return evolve_columns(np.linalg.eigh(h),np.eye(len(h)),t)


def eta(nuclear_d,nuclear_u,bmin):
    return (nuclear_d+nuclear_u)/math.sqrt(bmin+1)+2*nuclear_u/(bmin+1)


def once_calibration(md0,mu0,p):
    n,N0=2,3
    a=physical(n,1,N0,md0,mu0,p)
    compressed=a['E0'].conj().T@a['H']@a['E0']
    li={s:i for i,s in enumerate(a['logical'])}
    md=np.zeros_like(md0);mu=np.zeros_like(mu0)
    # Read the inverse amplitudes of actual natural-left preparations.
    for i in range(n):
        for j in range(n):
            md[i,j]=compressed[li[1<<j],li[1<<(2*n+i)]].conjugate()
            mu[i,j]=-compressed[li[1<<(n+j)],li[1<<(3*n+i)]].conjugate()/math.sqrt((N0-1)/N0)
    return md,mu,dict(q=1,N=N0,recovery_error=max(norm(md-md0),norm(mu-mu0)),
                      physical_input='natural left; full outputs include transverse leakage',
                      later_parameter_refits=0)


def run_case(n,q,N,md,mu,p):
    a=physical(n,q,N,md,mu,p)
    b=target(n,q,N,md,mu,p)
    e0,e1,h,hh=a['E0'],a['E1'],a['H'],b['H']
    residual=h@e0-e0@hh
    predicted=e0@b['C']+e1@b['L']
    observed_eta=norm(residual)
    bound=eta(nuclear(md),nuclear(mu),N-2*q)
    exact_res=max(norm(residual-predicted),
                  norm(residual.conj().T@residual-b['C'].conj().T@b['C']-b['L'].conj().T@b['L']),
                  norm(e0.conj().T@e1),norm(e0.conj().T@e0-np.eye(len(hh))))
    # Original independent sources are retained before the working energy relation.
    occ=a['occupations'];nb=a['nb'];locc=b['occupations'];lb=b['b']
    ed,el,eu,omega=.8+p,.8+.35-p,.8+.7+p,.35
    h0_independent=ed*occ[:,0]+eu*occ[:,1]+el*occ[:,2]+omega*nb
    energy_dictionary_res=float(np.max(np.abs(h0_independent-(.8*q+omega*N+np.diag(a['H0'])))))
    source_number_res=max(norm(occ[:,k,None]*e0-e0*locc[:,k][None,:]) for k in range(3))
    source_number_res=max(source_number_res,norm(nb[:,None]*e0-e0*lb[None,:]))
    coupling_sources=[]
    for sector in ['D','U']:
        dm=np.zeros((n,n),complex);um=np.zeros_like(dm)
        if sector=='D':dm[0,1]=1j
        else:um[1,0]=1
        pa=physical(n,q,N,dm,um,0.)
        ta=target(n,q,N,dm,um,0.)
        r=pa['V']@e0-e0@ta['H']
        decomp=norm(r-e0@ta['C']-e1@ta['L'])
        bound_source=eta(nuclear(dm),nuclear(um),N-2*q)
        coupling_sources.append(dict(sector=sector,observed_residual=norm(r),bound=bound_source,
                                     decomposition_error=decomp,passed=norm(r)<=bound_source+1e-11 and decomp<1e-11))
    # Every unknown all-left q-fermion state has one fixed, actual reference preparation.
    left_columns=[i for i,s in enumerate(a['logical']) if counts(s,n)[:2]==(0,0)]
    f0=e0[:,left_columns]
    embed=np.eye(len(hh))[:,left_columns]
    eigen=np.linalg.eigh(h);leigen=np.linalg.eigh(hh)
    tags=sorted(set(tuple(x) for x in occ.tolist()))
    physical_tag=[np.all(occ==tag,axis=1) for tag in tags]
    logical_tag=[np.all(locc==tag,axis=1) for tag in tags]
    # Exact record intertwiner, unlike the 1328 CG-direction recording question.
    record_intertwiner=max(norm(rows[:,None]*e0-e0*cols[None,:])
                           for rows,cols in zip(physical_tag,logical_tag))
    rng=np.random.default_rng(1329+q)
    probe=rng.normal(size=len(left_columns))+1j*rng.normal(size=len(left_columns))
    probe/=np.linalg.norm(probe)
    initial=f0@probe
    holdouts=[]
    for time in [.5,1.3,3.1]:
        actual_all=evolve_columns(eigen,e0,time)
        logical_all=evolve_columns(leigen,np.eye(len(hh)),time)
        actual=evolve_columns(eigen,f0,time)
        logical=evolve_columns(leigen,embed,time)
        desired=e0@logical
        full_error=norm(actual_all-e0@logical_all)
        natural_error=norm(actual-desired)
        amp_bound=min(2.,abs(time)*bound)
        epsilon=min(1.,abs(time)*bound)
        recorded=np.vstack([rows[:,None]*actual for rows in physical_tag])
        target_recorded=np.vstack([e0@(cols[:,None]*logical) for cols in logical_tag])
        recorded_error=norm(recorded-target_recorded)
        m=len(left_columns)
        record_gram=norm(recorded.conj().T@recorded-np.eye(m))
        partner_res=norm(recorded.conj().T@recorded/m-np.eye(m)/m)
        # Source and backreaction identities on the complete unknown input family.
        weight=occ[:,2]+2*occ[:,1]
        initial_energy=f0.conj().T@h@f0
        energy_res=norm(actual.conj().T@h@actual-initial_energy)
        flow_res=norm(
            actual.conj().T@((nb+weight)[:,None]*actual)-f0.conj().T@((nb+weight)[:,None]*f0))
        source_flow_error=norm(actual.conj().T@(weight[:,None]*actual)-
                               logical.conj().T@((locc[:,2]+2*locc[:,1])[:,None]*logical))
        off_work=-actual.conj().T@a['V']@actual
        recorded_energy=sum((rows[:,None]*actual).conj().T@a['H0']@(rows[:,None]*actual)
                            for rows in physical_tag)
        record_energy_res=norm(recorded_energy-actual.conj().T@a['H0']@actual)
        work_res=norm(recorded_energy-initial_energy-off_work)
        psi=actual@probe
        probabilities=[float(np.linalg.norm(psi[rows])**2) for rows in physical_tag]
        reference_delta=float(np.real(psi.conj()@(nb*psi)-initial.conj()@(nb*initial)))
        reference_leak_probability=float(np.linalg.norm(psi[a['leak_rows']])**2)
        number_predictions=[]
        for k,name in enumerate(['D','U','L']):
            physical_value=float(np.real(psi.conj()@(occ[:,k]*psi)))
            logical_psi=logical@probe
            target_value=float(np.real(logical_psi.conj()@(locc[:,k]*logical_psi)))
            number_predictions.append(dict(name=name,actual=physical_value,target=target_value))
        checks=[energy_res,flow_res,record_gram,partner_res,record_energy_res,work_res,
                abs(recorded_error-natural_error),abs(sum(probabilities)-1)]
        passed=(max(checks)<1e-9 and full_error<=amp_bound+1e-10 and
                source_flow_error<=2*q*epsilon+1e-10)
        holdouts.append(dict(time=time,full_conditional_code_error=full_error,
                             natural_left_error=natural_error,amplitude_bound=amp_bound,
                             reference_change_error=source_flow_error,reference_change_bound=2*q*epsilon,
                             complete_record_error=recorded_error,record_results=[
                                 dict(counts_D_U_L=list(tag),probability=prob) for tag,prob in zip(tags,probabilities)],
                             number_predictions=number_predictions,reference_transfer=reference_delta,
                             transverse_reference_probability=reference_leak_probability,
                             max_record_resource_residual=max(checks),passed=passed))
    passed=(exact_res<1e-10 and observed_eta<=bound+1e-11 and energy_dictionary_res<1e-10 and
            source_number_res<1e-10 and record_intertwiner<1e-11 and
            all(x['passed'] for x in coupling_sources+holdouts))
    return dict(n=n,q=q,N=N,parent_dimension=len(h),logical_dimension=len(hh),
                natural_left_dimension=len(left_columns),bmin=N-2*q,p=p,
                observed_generator_residual=observed_eta,generator_bound=bound,
                exact_orthogonal_decomposition_residual=exact_res,
                source_number_identity_residual=source_number_res,energy_dictionary_residual=energy_dictionary_res,
                record_intertwiner_residual=record_intertwiner,coupling_sources=coupling_sources,
                holdouts=holdouts,passed=passed)


def compute():
    md0=np.array([[.18,.03+.014j],[-.006j,.11]])
    mu0=np.array([[.08,.016j],[.024,.13]])
    p=.07
    md,mu,calibration=once_calibration(md0,mu0,p)
    cases=[run_case(2,q,N,md,mu,p) for q,N in [(1,3),(2,5),(2,9),(3,7),(3,11),(4,9)]]
    paths=[HERE/'proof.md',HERE/'check.py',HERE.parent/'_admission'/'round1329_admission.md',
           HERE.parent/'1327'/'proof.md',HERE.parent/'1328'/'proof.md']
    hashes={str(x.relative_to(HERE.parent)).replace('\\','/'):
            hashlib.sha256(x.read_bytes()).hexdigest() for x in paths}
    passed=calibration['recovery_error']<1e-12 and all(x['passed'] for x in cases)
    return dict(round=1329,status='PASS' if passed else 'FAIL',new_axioms=0,
                provenance=hashes,calibration=calibration,cases=cases,
                exact_complete_parent=True,empirical_experiments=0,
                actual_preparation='all unknown left q-fermion states and fixed reference N-q,0',
                arbitrary_logical_code_preparation_certified=False,
                uniform_all_time_small_error_claim=False,stage_complete=False,long_goal_complete=False)


def compare(a,b,path='root'):
    if isinstance(a,dict):
        assert set(a)==set(b),(path,'keys')
        for k in a:compare(a[k],b[k],path+'.'+k)
    elif isinstance(a,list):
        assert len(a)==len(b),(path,'length')
        for i,(x,y) in enumerate(zip(a,b)):compare(x,y,path+f'[{i}]')
    elif isinstance(a,float):
        assert math.isclose(a,b,rel_tol=1e-8,abs_tol=1e-9),(path,a,b)
    else:
        assert a==b,(path,a,b)


def main():
    parser=argparse.ArgumentParser()
    choice=parser.add_mutually_exclusive_group()
    choice.add_argument('--check',action='store_true')
    choice.add_argument('--write-new',action='store_true')
    args=parser.parse_args()
    r=compute();assert r['status']=='PASS',r
    result_path=HERE/'results.json'
    if args.write_new:
        with result_path.open('x',encoding='utf-8') as f:
            json.dump(r,f,ensure_ascii=False,indent=2);f.write('\n')
    else:compare(r,json.loads(result_path.read_text(encoding='utf-8')))
    print(json.dumps(dict(status='PASS',round=1329,complete_parents=len(r['cases']),
                         time_holdouts=sum(len(x['holdouts']) for x in r['cases']),
                         calibrated_once=True,results_written=args.write_new),ensure_ascii=False))


if __name__=='__main__':
    main()
