"""1325 independent representation, occupied-parent and reference-budget checks.

No author code is imported. The full small parent and a closed occupation
subspace of the large parent are separately constructed from ladder rules.
"""
from __future__ import annotations
import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
RESULT=HERE/'independent_results.json'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def op(x):
    return float(np.linalg.norm(x,2))


class Checks:
    def __init__(self): self.rows=[]
    def near(self,name,x,y=0.,tol=5e-10):
        r=float(np.linalg.norm(np.asarray(x)-np.asarray(y)))
        if not np.isfinite(r) or r>tol: raise AssertionError((name,r,tol))
        self.rows.append(dict(name=name,residual=r,tolerance=tol))
    def le(self,name,x,y,tol=5e-10):
        if not (np.isfinite(x) and np.isfinite(y) and x<=y+tol): raise AssertionError((name,x,y))
        self.rows.append(dict(name=name,value=float(x),upper=float(y),tolerance=tol))
    def yes(self,name,yes):
        if not yes: raise AssertionError(name)
        self.rows.append(dict(name=name,passed=True))


def cg(a,b):
    """Twice-spin notation; exact irreducible multiplicities."""
    return Counter({j:1 for j in range(abs(a-b),a+b+1,2)})


def direction(K):
    W=np.zeros((2*K,K+1),complex)
    for k in range(K+1):
        if k: W[k-1,k]=math.sqrt(k/K)
        if k<K: W[K+k,k]=math.sqrt((K-k)/K)
    J=np.zeros((2*K,max(0,K-1)),complex)
    for k in range(K-1):
        J[k,k]=math.sqrt((K-1-k)/K)
        J[K+k+1,k]=-math.sqrt((k+1)/K)
    return W,J


def representations(c):
    data=[]
    for K in [2,3,4,5,8]:
        target={K,K-2}
        allowed_D=[];allowed_U=[]
        for j in range(9):
            if set(cg(K,j)).issubset(target): allowed_D.append(j)
            if set(cg(K-2,j)).issubset(target): allowed_U.append(j)
        c.yes(f'K{K}.D_nontrivial_excluded',allowed_D==[0])
        c.yes(f'K{K}.U_role',allowed_U==([0,2] if K<4 else [0]))
        W,J=direction(K)
        c.near(f'K{K}.W_isometry',W.conj().T@W,np.eye(K+1))
        c.near(f'K{K}.J_isometry',J.conj().T@J,np.eye(K-1))
        c.near(f'K{K}.no_cross_port',W.conj().T@J)
        c.near(f'K{K}.whole_left',W@W.conj().T+J@J.conj().T,np.eye(2*K))
        # An equivariant full cover must match each distinct spin multiplicity.
        permitted=[(d,u,l) for d in range(1,5) for u in range(1,5) for l in range(1,5)
                   if Counter({K:d,K-2:u})==Counter({K:l,K-2:l})]
        c.yes(f'K{K}.single_K_total_multiplicity',permitted==[(n,n,n) for n in range(1,5)])
        data.append(dict(K=K,twice_spins_checked=list(range(9)),D_allowed=allowed_D,U_allowed=allowed_U))
    low=cg(3,0)+cg(1,0)+cg(1,2)
    c.yes('K3.full_cover_nontrivial_U_counterexample',low==Counter({3:2,1:2}))
    c.yes('K4.same_U_fails',0 in cg(2,2))
    W,J=direction(1)
    c.near('K1.no_U_boundary',W@W.conj().T,np.eye(2))
    c.yes('K1.J_empty',J.shape==(2,0))
    return dict(finite_checks=data,
                all_spins_supported_by_analytic_highest_weight_and_lowest_spin_argument=True,
                trivial_copies_absorbed_into_total_multiplicities=True,
                K3_counterexample=dict(WD_twice_spins=[0],WU_twice_spins=[0,2],nD=1,nU=1,nL=2),
                full_coverage_is_additional_contract=True)


def weights(nu,L):
    logs=np.array([k*math.log(nu)-math.lgamma(k+1) for k in range(L+1)])
    p=np.exp(logs-logs.max());return p/p.sum()


def moments(p):
    L=len(p)-1
    s=float(np.sqrt(p[:-1]*p[1:]).sum())
    z=float(np.dot(np.sqrt(p[:-1]*p[1:]),np.sqrt(np.arange(1,L+1)/np.arange(2,L+2))))
    cc=float(np.dot(p,1/(np.arange(L+1)+1)))
    kap=float(np.dot(np.sqrt(p[:-2]*p[2:]),np.sqrt(np.arange(1,L)/np.arange(2,L+1)))) if L>=2 else 0.
    return dict(s=s,z=z,c=cc,q=1-cc-z*z,kappa=kap)


def hermitian_evolve(H,t):
    vals,vec=np.linalg.eigh(H)
    return (vec*np.exp(-1j*t*vals))@vec.conj().T


def matrix_on_states(states,MD,MU):
    """Build each occupied-space conversion directly; one CAR bilinear sign is +."""
    n=len(MD);index={state:i for i,state in enumerate(states)}
    A=np.zeros((len(states),len(states)),complex)
    for col,(loc,j,a,b) in enumerate(states):
        if loc==0 and a+b:
            for alpha,count in enumerate((a,b)):
                if count:
                    for i in range(n):
                        A[index[alpha+1,i,a-int(alpha==0),b-int(alpha==1)],col]=MD[i,j]*math.sqrt(count/(a+b))
        elif loc==3:
            for alpha,beta,sign in [(0,1,1),(1,0,-1)]:
                count=(a,b)[beta]
                for i in range(n):
                    A[index[alpha+1,i,a+int(beta==0),b+int(beta==1)],col]=sign*MU[i,j]*math.sqrt((count+1)/(a+b+2))
    return A+A.conj().T


def full_parent(L,MD,MU):
    n=len(MD);cap=L+2
    states=[(loc,j,a,b) for loc in range(4) for j in range(n)
            for total in range(cap-[0,1,1,2][loc]+1)
            for a in range(total+1) for b in [total-a]]
    return states,matrix_on_states(states,MD,MU)


def embed(states,n,p,phi,u):
    E=np.zeros((len(states),4*n),complex)
    for row,(loc,j,a,b) in enumerate(states):
        k=a+b
        if k<len(p):
            E[row,loc*n+j]=math.sqrt(p[k]*math.comb(k,a))*np.exp(1j*phi*k)*u[0]**a*u[1]**b
    return E


def h_eff(MD,MU,m,phi,u):
    n=len(MD);h=np.zeros((4*n,4*n),complex)
    down=m['s']*np.exp(1j*phi)*u
    up=m['z']*np.exp(-1j*phi)*np.array([u[1].conjugate(),-u[0].conjugate()])
    for alpha in range(2):
        rows=slice((alpha+1)*n,(alpha+2)*n)
        h[rows,:n]=down[alpha]*MD
        h[rows,3*n:]=up[alpha]*MU
    return h+h.conj().T


def expected_gram(MD,MU,p,m,u):
    n=len(MD);Ad=MD@MD.conj().T;Au=MU@MU.conj().T
    P=np.outer(u,u.conj());G=np.zeros((4*n,4*n),complex)
    G[:n,:n]=(1-p[0]-m['s']**2)*(MD.conj().T@MD)
    G[n:3*n,n:3*n]=np.kron(P,(1-m['s']**2)*Ad)+np.kron(np.eye(2)-P,m['c']*Ad+m['q']*Au)
    G[3*n:,3*n:]=(1-m['z']**2)*(MU.conj().T@MU)
    return G


def pure_distance(a,b):
    return math.sqrt(max(0.,1-abs(np.vdot(a,b))**2))


def observables(states):
    nb=np.array([a+b for loc,j,a,b in states],float)
    left=np.array([loc in (1,2) for loc,j,a,b in states],float)
    up=np.array([loc==3 for loc,j,a,b in states],float)
    K=nb+left+2*up
    return nb,left,up,K


def ledger(c,name,states,E,A,VA,VE):
    nb,pl,pu,K=observables(states)
    dn=A.conj().T@(nb[:,None]*A)-E.conj().T@(nb[:,None]*E)
    dl=A.conj().T@(pl[:,None]*A)-E.conj().T@(pl[:,None]*E)
    du=A.conj().T@(pu[:,None]*A)-E.conj().T@(pu[:,None]*E)
    c.near(name+'.complete_reference_identity',dn+dl+2*du)
    H0=1*(1-pl-pu)+1.4*pl+1.8*pu+.4*nb
    c.near(name+'.resonant_free_source',H0,1+.4*K)
    c.near(name+'.free_energy',A.conj().T@(H0[:,None]*A),E.conj().T@(H0[:,None]*E),tol=3e-9)
    c.near(name+'.interaction_energy',A.conj().T@VA,E.conj().T@VE)
    c.near(name+'.total_energy',A.conj().T@(H0[:,None]*A+VA),E.conj().T@(H0[:,None]*E+VE),tol=3e-9)
    return dn,dl,du


def small_checks(c):
    MD=np.array([[.7+.1j,.2j],[.13-.17j,.4]])
    MU=np.array([[.2+.3j,-.21j],[.1,.52-.12j]])
    directions=[np.array([1,0],complex),np.array([math.sqrt(.37)*np.exp(.27j),math.sqrt(.63)*np.exp(-.51j)])]
    rng=np.random.default_rng(1325);result=[]
    for L,nu in [(1,.6),(2,1.3),(3,2.2)]:
        states,V=full_parent(L,MD,MU);n=len(MD);d=len(states)
        nb,pl,pu,K=observables(states)
        c.yes(f'small.L{L}.full_dimension',d==n*(2*(L+2)**2+2*(L+2)+1))
        c.near(f'small.L{L}.K_conserved',(K[:,None]-K[None,:])*V)
        c.near(f'small.L{L}.parent_norm',op(V),max(op(MD),op(MU)))
        c.yes(f'small.L{L}.only_vacuum_global_dark',int(np.sum(abs(np.linalg.eigvalsh(V))<1e-9))==n)
        p=weights(nu,L);m=moments(p)
        c.near(f'small.L{L}.kappa_finite_identity',m['kappa'],1-m['c']-L/(L+1)*p[-1])
        c.yes(f'small.L{L}.positive_up_variance',m['q']>=-1e-13)
        for iu,u in enumerate(directions):
            name=f'small.L{L}.u{iu}';phi=.4
            E=embed(states,n,p,phi,u);h=h_eff(MD,MU,m,phi,u);VE=V@E
            G=expected_gram(MD,MU,p,m,u);D=VE-E@h;eta=math.sqrt(float(np.linalg.eigvalsh(G)[-1]))
            c.near(name+'.input_isometry',E.conj().T@E,np.eye(4*n))
            c.near(name+'.mean_generator',E.conj().T@VE,h)
            c.near(name+'.whole_error_Gram',D.conj().T@D,G)
            c.near(name+'.reference_orthogonal_error',E.conj().T@D)
            c.near(name+'.eta',op(D),eta)
            c.near(name+'.D_U_cross_Gram',(D.conj().T@D)[:n,3*n:])
            t=.47;A=hermitian_evolve(V,t)@E;I=E@hermitian_evolve(h,t)
            c.near(name+'.evolved_isometry',A.conj().T@A,np.eye(4*n))
            c.le(name+'.joint_Duhamel',op(A-I),t*eta)
            c.le(name+'.full_partner',pure_distance(A/math.sqrt(4*n),I/math.sqrt(4*n)),min(1,t*eta))
            z=rng.normal(size=(4*n,3))+1j*rng.normal(size=(4*n,3));z/=np.linalg.norm(z)
            c.le(name+'.arbitrary_partner',pure_distance(A@z,I@z),min(1,t*eta))
            ledger(c,name,states,E,A,V@A,VE)
            leakage=A-E@(E.conj().T@A)
            c.le(name+'.actual_reference_Taylor',abs(op(leakage)/t-eta),t*op(V)**2/2)
            result.append(dict(L=L,nu=nu,direction=iu,full_dimension=d,eta=eta,
                               actual_reference_leak_worst=float(op(leakage)**2)))
    return result


def closed_parent(L,MD,MU):
    """Closed axial-input occupation space, retaining both common-left arrows."""
    n=len(MD);groups=[[(0,0,0)]]
    for K in range(1,L+3):
        groups.append([(0,K,0),(1,K-1,0)])
        shared=[(0,K-1,1),(2,K-1,0)]
        if K>=2: shared += [(1,K-2,1),(3,K-2,0)]
        groups.append(shared)
    nodes=sorted({s for group in groups for s in group})
    states=[(loc,j,a,b) for loc,a,b in nodes for j in range(n)]
    index={s:i for i,s in enumerate(states)}
    blocks=[]
    for group in groups:
        local=[(loc,j,a,b) for loc,a,b in group for j in range(n)]
        H=matrix_on_states(local,MD,MU)
        blocks.append((np.array([index[s] for s in local]),H))
    # Every actual nonzero ladder transition stays inside this common space.
    tests=[]
    for loc,j,a,b in states:
        if loc==0 and a+b:
            for alpha,count in enumerate((a,b)):
                if count: tests.extend((alpha+1,i,a-int(alpha==0),b-int(alpha==1)) in index for i in range(n))
        if loc==3:
            for alpha,beta in [(0,1),(1,0)]:
                tests.extend((alpha+1,i,a+int(beta==0),b+int(beta==1)) in index for i in range(n))
        if loc in (1,2):
            alpha=loc-1
            tests.extend((0,i,a+int(alpha==0),b+int(alpha==1)) in index for i in range(n))
            beta=1-alpha
            if (a,b)[beta]:
                tests.extend((3,i,a-int(beta==0),b-int(beta==1)) in index for i in range(n))
    return states,blocks,all(tests)


def apply_blocks(blocks,X,time=None):
    Y=np.zeros_like(X)
    for idx,H in blocks:
        Y[idx]=(H if time is None else hermitian_evolve(H,time))@X[idx]
    return Y


def reduced_cross(states,n,a,b):
    refs=sorted({(x,y) for loc,j,x,y in states});idx={r:i for i,r in enumerate(refs)}
    A=np.zeros((4*n,len(refs)),complex);B=np.zeros_like(A)
    for row,(loc,j,x,y) in enumerate(states):
        A[loc*n+j,idx[x,y]]=a[row];B[loc*n+j,idx[x,y]]=b[row]
    return A@B.conj().T


def coherent_witness(c):
    result=[];gamma=.1;T=math.pi/(2*gamma);phi=.37
    MD=MU=np.array([[gamma]])
    for nu,L in [(16.,64),(64.,160)]:
        name=f'coherent.nu{nu:g}.L{L}';p=weights(nu,L);m=moments(p)
        states,blocks,closed=closed_parent(L,MD,MU);E=embed(states,1,p,phi,np.array([1,0]))
        c.yes(name+'.occupation_closure',closed)
        c.yes(name+'.kept_dimension',len(states)==6*L+11)
        h=h_eff(MD,MU,m,phi,np.array([1,0]));VE=apply_blocks(blocks,E)
        D=VE-E@h;eta=gamma*math.sqrt(1-m['z']**2)
        c.near(name+'.all_input_Gram',D.conj().T@D,expected_gram(MD,MU,p,m,np.array([1,0])))
        c.near(name+'.max_error',op(D),eta)
        A=apply_blocks(blocks,E,T);I=E@hermitian_evolve(h,T)
        c.near(name+'.complete_isometry',A.conj().T@A,np.eye(4))
        c.le(name+'.full_joint_Duhamel',op(A-I),T*eta)
        c.le(name+'.full_partner',pure_distance(A/2,I/2),min(1,T*eta))
        cross=reduced_cross(states,1,A[:,0],A[:,3])
        target=np.zeros((4,4),complex);target[1,2]=-np.exp(2j*phi)*m['kappa']
        c.near(name+'.actual_cross_port_coherence',cross,target)
        psi=(A[:,0]+A[:,3])/math.sqrt(2)
        rho=reduced_cross(states,1,psi,psi)
        X=np.zeros((4,4),complex);X[1,2]=-np.exp(2j*phi);X[2,1]=-np.exp(-2j*phi)
        readout=float(np.trace(rho@X).real)
        c.near(name+'.actual_X_readout',readout,m['kappa'])
        c.near(name+'.kappa_recursion',m['kappa'],1-m['c']-L/(L+1)*p[-1])
        c.near(name+'.kappa_other_finite_sum',m['kappa'],nu*np.dot(p[:-2],1/np.arange(2,L+1)))
        dn,dl,du=ledger(c,name,states,E,A,apply_blocks(blocks,A),VE)
        c.near(name+'.D_probability',dl[0,0],1-p[0])
        c.near(name+'.U_probability',dl[3,3],1.)
        c.near(name+'.D_reference_loss',dn[0,0],-(1-p[0]))
        c.near(name+'.U_reference_gain',dn[3,3],1.)
        # Same parent, a fixed-number preparation: trace over reference kills the cross-port readout.
        fixed=np.zeros(L+1);fixed[4]=1.
        EF=embed(states,1,fixed,phi,np.array([1,0]));AF=apply_blocks(blocks,EF,T)
        c.near(name+'.fixed_number_cross_dephased',reduced_cross(states,1,AF[:,0],AF[:,3]))
        phases=2*math.pi*np.arange(L+1)/(L+1)
        c.near(name+'.phase_twirl_cross_zero',np.exp(2j*phases).mean()*m['kappa'])
        result.append(dict(nu=nu,L=L,occupation_dimension=len(states),gamma=gamma,time=T,
                           s=m['s'],z=m['z'],c=m['c'],kappa=m['kappa'],
                           measured_X=readout,eta=eta,Duhamel_bound=min(1,T*eta),
                           actual_joint_difference=op(A-I),
                           reference_change_D=float(dn[0,0].real),reference_change_U=float(dn[3,3].real)))
    return result


def budget_witness(c):
    nu,L=16.,64;n=3;p=weights(nu,L);m=moments(p);phi=.19
    # All singular values are odd multiples of .01: same exact handoff clock.
    MD=np.diag([.13,.03,.01]).astype(complex)
    base=np.diag([.11,.05,.01]).astype(complex)
    perm=np.array([[0,0,1],[0,1,0],[1,0,0]],complex)
    beta=.0337;tau=.1;Thandoff=math.pi/.02
    outputs=[]
    for tag,MU in [('aligned',base),('reversed',perm@base)]:
        name='budget.'+tag
        states,blocks,closed=closed_parent(L,MD,MU)
        c.yes(name+'.closure',closed)
        c.yes(name+'.complete_kept_dimension',len(states)==n*(6*L+11))
        E=embed(states,n,p,phi,np.array([1,0]));h=h_eff(MD,MU,m,phi,np.array([1,0]))
        VE=apply_blocks(blocks,E);D=VE-E@h;G=expected_gram(MD,MU,p,m,np.array([1,0]))
        eta=op(D);gmax=max(op(MD),op(MU));shared=m['c']*(MD@MD.conj().T)+m['q']*(MU@MU.conj().T)
        other=math.sqrt(max((1-m['s']**2)*op(MD)**2,(1-m['z']**2)*op(MU)**2))
        c.near(name+'.common_Gram',D.conj().T@D,G)
        c.near(name+'.initial_reference_projection',(VE-E@(E.conj().T@VE)),D)
        c.near(name+'.mean_error_orthogonal',E.conj().T@D)
        c.near(name+'.budget_shared_maximum',eta,math.sqrt(float(np.linalg.eigvalsh(shared)[-1])))
        c.yes(name+'.other_blocks_below_budget',other<beta)
        c.near(name+'.same_U_singular_values',np.linalg.svd(MU,compute_uv=False),np.linalg.svd(base,compute_uv=False))
        selected=2*n
        c.near(name+'.same_fixed_input_attains_worst',np.linalg.norm(D[:,selected]),eta)
        A=apply_blocks(blocks,E,tau);leak=A-E@(E.conj().T@A)
        rate_fixed=float(np.linalg.norm(leak[:,selected])/tau);rate_worst=op(leak)/tau
        remainder=tau*gmax*gmax/2
        c.le(name+'.finite_Taylor_worst',abs(rate_worst-eta),remainder)
        c.le(name+'.finite_Taylor_fixed',abs(rate_fixed-eta),remainder)
        if tag=='aligned':
            c.yes(name+'.reject_asymptotic_budget',eta>beta)
            c.yes(name+'.finite_same_input_violation',rate_fixed>beta)
            c.yes(name+'.separation_even_with_remainder',eta-remainder>beta)
        else:
            c.yes(name+'.passes_asymptotic_budget',eta<beta)
            c.yes(name+'.finite_all_input_pass',rate_worst<beta)
            c.yes(name+'.separation_even_with_remainder',eta+remainder<beta)
        ledger(c,name,states,E,A,apply_blocks(blocks,A),VE)
        # Exact all-unknown fixed-K task on a complete occupied parent, not just chi.
        K=4;full,H=full_parent(K-2,MD,MU)
        sector=np.array([i for i,(loc,j,a,b) in enumerate(full) if a+b+[0,1,1,2][loc]==K])
        HH=H[np.ix_(sector,sector)];UU=hermitian_evolve(HH,Thandoff)
        locs=[full[i][0] for i in sector]
        right=np.array([i for i,x in enumerate(locs) if x in (0,3)])
        left=np.array([i for i,x in enumerate(locs) if x in (1,2)])
        output=UU[:,right]
        c.near(name+'.actual_same_clock_right_failure',output[right])
        c.near(name+'.actual_same_clock_full_unknown_isometry',output[left].conj().T@output[left],np.eye(len(right)))
        c.near(name+'.actual_same_clock_full_left_cover',output[left]@output[left].conj().T,np.eye(len(left)))
        outputs.append(dict(orientation=tag,eta=eta,Gamma=eta*eta,other_blocks_bound=other,
                            actual_probability_reference_leaves=float(np.linalg.norm(leak[:,selected])**2),
                            actual_rate_fixed=rate_fixed,actual_rate_worst=rate_worst,
                            remainder=remainder,fixed_input_system_index=selected,
                            relative_invariant=float(np.trace((MD@MD.conj().T)@(MU@MU.conj().T)).real),
                            exact_K4_handoff_time=Thandoff))
    c.yes('budget.relative_invariant_changes',outputs[0]['relative_invariant']!=outputs[1]['relative_invariant'])
    return dict(nu=nu,L=L,beta=beta,time=tau,probability_threshold=(beta*tau)**2,
                D_singular_values=[.13,.03,.01],U_singular_values=[.11,.05,.01],
                shared_exact_handoff_time=Thandoff,examples=outputs,
                budget_is_additional_not_natural_constant=True)


def compute():
    c=Checks();rep=representations(c);small=small_checks(c);coh=coherent_witness(c);budget=budget_witness(c)
    paths=[HERE/'independent_check.py',HERE/'proof.md',STAGE/'_admission'/'round1325_admission.md',
           STAGE/'1323'/'proof.md',STAGE/'1324'/'proof.md']
    return dict(schema='round1325_independent_v1',status='PASS_INDEPENDENT_CHECKS_ONLY',
                scientific_acceptance=False,author_code_imported=False,
                checks=len(c.rows),max_absolute_residual=max(x['residual'] for x in c.rows if 'residual' in x),
                representation=rep,small_full_parent=small,coherent_joint_witness=coh,
                operational_reference_budget=budget,evidence=c.rows,
                scope=dict(all_unknown_fixed_K_task_distinct_from_product_reference_task=True,
                           coherent_reference_preparation_and_projection_permission_are_inputs=True,
                           large_closed_subspace_contains_all_declared_input_columns=True,
                           inactive_parent_complement_not_deleted=True,
                           complete_local_field_theory_not_certified=True,new_axioms=0),
                source_bindings={p.relative_to(STAGE).as_posix():digest(p) if p.exists() else None for p in paths})


def compare(a,b,path='result'):
    if isinstance(a,dict):
        if not isinstance(b,dict) or a.keys()!=b.keys(): raise AssertionError(path+': keys')
        for k in a: compare(a[k],b[k],path+'.'+k)
    elif isinstance(a,list):
        if not isinstance(b,list) or len(a)!=len(b): raise AssertionError(path+': list')
        for i,(x,y) in enumerate(zip(a,b)): compare(x,y,path+f'[{i}]')
    elif isinstance(a,float):
        if isinstance(b,bool) or not isinstance(b,(float,int)) or not np.isfinite(b) or not math.isclose(a,b,rel_tol=5e-10,abs_tol=5e-12):
            raise AssertionError((path,a,b))
    elif type(a) is not type(b) or a!=b: raise AssertionError((path,a,b))


def main():
    parser=argparse.ArgumentParser();modes=parser.add_mutually_exclusive_group()
    modes.add_argument('--check',action='store_true');modes.add_argument('--write-new',action='store_true')
    modes.add_argument('--preview',action='store_true');args=parser.parse_args()
    if args.write_new and RESULT.exists(): raise FileExistsError(RESULT)
    if not args.write_new and not args.preview and not RESULT.is_file(): raise FileNotFoundError(RESULT)
    r=compute()
    if args.write_new:
        if any(x is None for x in r['source_bindings'].values()): raise AssertionError('Unfinished source binding')
        with RESULT.open('x',encoding='utf-8',newline='\n') as f:
            json.dump(r,f,ensure_ascii=False,indent=2,allow_nan=False);f.write('\n')
    elif not args.preview: compare(r,json.loads(RESULT.read_text(encoding='utf-8')))
    print(json.dumps(dict(status='PASS',mode='preview' if args.preview else 'write-new' if args.write_new else 'check',
                          checks=r['checks'],max_residual=r['max_absolute_residual'],
                          reference_budget_rates=[x['actual_rate_fixed'] for x in r['operational_reference_budget']['examples']]),indent=2))


if __name__=='__main__': main()
