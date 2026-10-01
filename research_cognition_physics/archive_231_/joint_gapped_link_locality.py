"""610: common gap controls actual Wilson link derivatives and Fock coefficients.

Four Euclidean grid dimensions are an input. No unbounded-gauge LR claim.
Trace-norm estimates, not one-particle norm times particle count, control dGamma.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np
import joint_chiral_fibre_source as prior
HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_gapped_link_locality_results.json'
M=12.;JROW=11.;GAP=.8;RANGE=1.
MU=float(np.log1p(GAP/(4*JROW)))
CONTOUR=2*(M+GAP)/(2*np.pi)
CP=CONTOUR*4/GAP

def tn(A,kind='hermitian'):
    if kind=='anti':
        return float(np.sum(abs(np.linalg.eigvalsh(1j*A))))
    if kind=='hermitian':
        return float(np.sum(abs(np.linalg.eigvalsh((A+A.conj().T)/2))))
    return float(np.sum(np.linalg.svd(A,compute_uv=False)))
def op(A):
    return float(np.linalg.norm(A,2))
def kernel(L,phases=None,derivative=None,charge=1):
    sites=list(itertools.product(range(L),repeat=4));index={x:i for i,x in enumerate(sites)}
    phases={} if phases is None else phases
    n=len(sites);D=np.zeros((4*n,4*n),complex);dD=np.zeros_like(D)
    for mu in range(4):
        T=np.zeros((n,n),complex);dT=np.zeros_like(T)
        for i,x in enumerate(sites):
            y=list(x);y[mu]=(y[mu]+1)%L
            phase=np.exp(1j*charge*phases.get((x,mu),0.))
            T[i,index[tuple(y)]]=phase
            if derivative==(x,mu):dT[i,index[tuple(y)]]=1j*charge*phase
        D+=np.kron(np.eye(n)-.5*(T+T.conj().T),np.eye(4))+.5*np.kron(T-T.conj().T,prior.GAMMA[mu])
        dD+=-.5*np.kron(dT+dT.conj().T,np.eye(4))+.5*np.kron(dT-dT.conj().T,prior.GAMMA[mu])
    gamma5=np.kron(np.eye(n),prior.G5)
    return gamma5@(D-np.eye(4*n)),gamma5@dD,sites
def projection(H,D):
    E,V=np.linalg.eigh(H);select=(E>0).astype(float)
    p=(V*select)@V.conj().T
    denominator=E[:,None]-E[None,:];numerator=select[:,None]-select[None,:]
    divided=np.divide(numerator,denominator,out=np.zeros_like(denominator),where=abs(denominator)>1e-10)
    b=V@(divided*(V.conj().T@D@V))@V.conj().T
    c=p@b-b@p;s=b@b
    return dict(p=p,b=b,c=c,s=s,gap=float(min(abs(E))),maximum=float(max(abs(E))))
def distance(x,y,L):
    return sum(min(abs(a-b),L-abs(a-b)) for a,b in zip(x,y))
def radii(sites,S,L):
    return np.array([min(distance(x,y,L) for y in S) for x in sites])
def constants(D):
    d1=tn(D)
    cb=CONTOUR*(4/GAP)**2*d1
    cc=2*CP*cb;cs=cb*cb
    db=64*CONTOUR*M*d1/GAP**3
    dc=2*db+32*CONTOUR*M*cb/GAP**2
    ds=2*cb*db
    return dict(D_trace_norm=d1,projection_bound=CP,derivative_bound=cb,
                connection_bound=cc,positive_bound=cs,
                remote_derivative_constant=db,remote_connection_constant=dc,
                remote_positive_constant=ds)

def actual_and_derivative_check():
    cache={};rows=[]
    ell=((0,0,0,0),0);phases={ell:.02}
    # Full-direction phase reconstructs the inherited 606 kernel exactly.
    all_phases={(x,0):.013 for x in itertools.product(range(2),repeat=4)}
    H0,_,_=kernel(2,all_phases,ell)
    Hprev,_,_=prior.wilson(.013,1)
    inherited=op(H0-Hprev);assert inherited<1e-12
    for L in (2,3):
        H,D,sites=kernel(L,phases,ell);data=projection(H,D);const=constants(D)
        eps=2e-6
        plus=kernel(L,{ell:.02+eps},ell)[0]
        minus=kernel(L,{ell:.02-eps},ell)[0]
        fd=(projection(plus,D)['p']-projection(minus,D)['p'])/(2*eps)
        error=op(fd-data['b'])
        assert error<2e-8 and data['gap']>GAP and data['maximum']<M
        rowmax=0.
        for x in range(len(sites)):
            rowmax=max(rowmax,sum(op(H[4*x:4*x+4,4*y:4*y+4]) for y in range(len(sites))))
        assert rowmax<=JROW+1e-12
        assert op(data['p']@data['p']-data['p'])<1e-12
        assert op(data['b']+data['c']@data['p']-data['p']@data['c'])<1e-12
        assert np.linalg.eigvalsh(data['s'])[0]>-1e-12
        S=(ell[0],(1,0,0,0));rad=radii(sites,S,L)
        cache[L]=(H,D,sites,data,const,rad)
        rows.append(dict(L=L,sites=len(sites),dimension=len(H),gap=data['gap'],
            row_block_norm_bound=rowmax,derivative_difference_error=error,
            derivative_operator_norm=op(data['b']),derivative_trace_norm=tn(data['b']),
            positive_trace=tn(data['s']),constants=const))
    # The local perturbation is gauge covariant, including its derivative.
    L=2;sites=cache[L][2];angle={x:.07*np.sin(i) for i,x in enumerate(sites)}
    transformed={}
    for x in sites:
        for mu in range(4):
            y=list(x);y[mu]=(y[mu]+1)%L;y=tuple(y)
            transformed[(x,mu)]=phases.get((x,mu),0)+angle[x]-angle[y]
    Hg,Dg,_=kernel(L,transformed,ell);original=cache[L]
    U=np.kron(np.diag([np.exp(1j*angle[x]) for x in sites]),np.eye(4))
    gauge=max(op(Hg-U@original[0]@U.conj().T),op(Dg-U@original[1]@U.conj().T))
    assert gauge<1e-12
    return dict(inherited_kernel_error=inherited,gauge_covariance_error=gauge,rows=rows),cache

def weighted_and_truncation_check(cache):
    rows=[]
    for L,(H,D,sites,data,k,rad) in cache.items():
        weight=np.repeat(np.exp(MU*rad),4)
        wp=weight[:,None]*data['p']/weight[None,:]
        assert op(wp)<=CP+1e-10
        measured={}
        for name,key,kind in [('b','derivative_bound','hermitian'),
                              ('c','connection_bound','anti'),('s','positive_bound','hermitian')]:
            weighted=weight[:,None]*data[name]*weight[None,:]
            value=tn(weighted,kind);assert value<=k[key]+1e-9
            measured[name]=value
        tails=[]
        for r in range(int(max(rad))+1):
            mask=np.repeat(rad<=r,4)
            values={}
            for name,key,kind in [('b','derivative_bound','hermitian'),
                                  ('c','connection_bound','anti'),('s','positive_bound','hermitian')]:
                truncated=mask[:,None]*data[name]*mask[None,:]
                value=tn(data[name]-truncated,kind)
                bound=2*k[key]*np.exp(-MU*r)
                assert value<=bound+1e-9
                values[name]=value
            tails.append(dict(radius=r,tail_trace_norms=values))
        assert tails[-1]['tail_trace_norms']['b']<1e-12
        rows.append(dict(L=L,weighted_projection_norm=op(wp),
            weighted_trace_norms=measured,tails=tails))
    return dict(rows=rows,mu=MU,volume_independence_proved_not_fitted=True,
        dGamma_error_bounded_by_trace_norm_not_particle_number=True,
        fermion_index_support_alone_not_full_gauge_locality=True)

def remote_configuration_check(cache):
    ell=((0,0,0,0),0);rows=[]
    for L,(H,D,sites,old,k,rad) in cache.items():
        far=max(sites,key=lambda x:distance(x,ell[0],L))
        remote=(far,1)
        phases={ell:.02,remote:.03}
        H1,D1,_=kernel(L,phases,ell);new=projection(H1,D1)
        assert op(D1-D)<1e-12 and new['gap']>GAP
        delta=H1-H
        support=np.array([i for i in range(len(sites)) if np.linalg.norm(delta[4*i:4*i+4,:])>1e-12])
        r=float(min(rad[support]));assert r>0
        assert op(delta)<=2*M
        changed={}
        for name,key,kind in [('b','remote_derivative_constant','hermitian'),
                              ('c','remote_connection_constant','anti'),
                              ('s','remote_positive_constant','hermitian')]:
            value=tn(new[name]-old[name],kind)
            assert value<=k[key]*np.exp(-MU*r)+1e-9
            changed[name]=value
        assert changed['c']>1e-8
        rows.append(dict(L=L,remote_link=[list(far),1],distance_to_changed_kernel=r,
            second_gap=new['gap'],changes_trace_norm=changed,
            same_local_derivative_error=op(D1-D)))
    return dict(rows=rows,both_configurations_must_satisfy_same_gap=True,
        coefficient_depends_on_remote_gauge_field_but_has_conditional_tail_bound=True,
        no_claim_of_cylinder_extension_preserving_gap=True)

def run():
    actual,cache=actual_and_derivative_check()
    result=dict(actual=actual,locality=weighted_and_truncation_check(cache),
                gauge_dependence=remote_configuration_check(cache))
    deps=('research_note_603.md','research_note_605.md','research_note_606.md',
          'research_note_608.md','research_note_609.md','joint_chiral_fibre_source.py',
          'joint_fock_covariant_completion.py')
    return dict(round=610,tests_run=3,failures=0,errors=0,**result,
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope=dict(uniform_gap_finite_range_local_derivative_contract=True,
            shared_projector_connection_positive_Fock_bounds=True,
            same_gap_remote_configuration_stability=True,
            actual_Euclidean_Wilson_link_not_global_holonomy=True,
            full_original_thermal_support_not_replaced_by_hard_domain=True,
            no_unbounded_gauge_LR_no_SM_mass_dictionary_no_GR=True))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps(dict(round=610,tests=3,all_passed=True,mu=MU,
        derivative_rows=result['actual']['rows'],remote=result['gauge_dependence']['rows'])))
