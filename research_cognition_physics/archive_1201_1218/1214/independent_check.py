"""Independent 1214 checks. No author implementation is imported.
Default mode only compares and emits compact PASS JSON; --write never overwrites.
Finite examples audit the analytic statements; they do not prove all-n classification.
"""
import argparse
import hashlib
import itertools
import json
import math
from fractions import Fraction as F
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
STAGE=HERE.parent

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def entropy(rho):
    p=np.linalg.eigvalsh((rho+rho.conj().T)/2)
    assert p.min()>-1e-11
    p=p[p>1e-14]
    return float(-np.dot(p,np.log(p)))
def trace_distance(a,b):
    return float(np.sum(np.abs(np.linalg.eigvalsh((a-b+(a-b).conj().T)/2)))/2)
def partial(rho,da,db,keep):
    a=rho.reshape(da,db,da,db)
    return np.einsum('abcb->ac',a) if keep==0 else np.einsum('abad->bd',a)
def relative(rho,sigma):
    vals,vec=np.linalg.eigh((sigma+sigma.conj().T)/2)
    assert vals.min()>0
    log_sigma=(vec*np.log(vals))@vec.conj().T
    return float(-entropy(rho)-np.trace(rho@log_sigma).real)

def first_hit():
    # Two 3-level blocks; first local hit of {0,2}; unit level spacing.
    strings=list(itertools.product(range(3),repeat=2));n=len(strings)
    idx={s:i for i,s in enumerate(strings)}
    perm=[];flag=[]
    for s in strings:
        t=list(s);hit=next((j for j,x in enumerate(s) if x in (0,2)),None)
        flag.append(int(hit is not None and s[hit]==0))
        if hit is not None:t[hit]=2-t[hit]
        perm.append(idx[tuple(t)])
    perm=np.array(perm);energy=np.array([sum(s) for s in strings])
    assert np.array_equal(perm[perm],np.arange(n))
    L=5;nc=L+2;nk=2;da=nc*nk
    labels=list(itertools.product(range(n),range(nc),range(nk)))
    loc={t:i for i,t in enumerate(labels)};full=np.arange(len(labels))
    for s in range(n):
        t=int(perm[s])
        if energy[s]>=energy[t]:continue
        assert energy[t]-energy[s]==2
        for k in range(1,L+2):
            for rec in range(nk):
                a=loc[(s,k,rec)];b=loc[(t,k-1,rec)]
                assert full[a]==a and full[b]==b
                full[a]=b;full[b]=a
    assert np.array_equal(full[full],np.arange(len(labels)))
    total=np.array([energy[s]+2*k for s,k,j in labels])
    assert np.array_equal(total[full],total)
    # Apply QND output flag after the lift. This full permutation need not be involutive.
    copy=np.array([loc[(s,k,j^flag[s])] for s,k,j in labels])
    combined=copy[full]
    assert len(set(map(int,combined)))==len(labels)
    assert np.array_equal(total[combined],total)
    # Columns retain integer amplitudes before a common 1/sqrt(L).
    vi=np.zeros((n,da,n),dtype=np.int64)
    for s in range(n):
        for k in range(1,L+1):
            out=labels[int(combined[loc[(s,k,0)]])]
            vi[out[0],out[1]*nk+out[2],s]+=1
    flat=vi.reshape(n*da,n)
    assert np.array_equal(flat.T@flat,L*np.eye(n,dtype=np.int64))
    image_i=np.einsum('ias,jas->ij',vi,vi)
    assert np.array_equal(image_i,L*np.eye(n,dtype=np.int64))
    aux=np.array([vi[int(perm[s]),:,s] for s in range(n)])
    assert all(np.count_nonzero(vi[:,:,s])==L for s in range(n))
    gram=aux@aux.T
    assert np.array_equal(np.diag(gram),L*np.ones(n,dtype=np.int64))
    assert np.linalg.eigvalsh(gram).min()>-1e-12
    rng=np.random.default_rng(1214001)
    sr=rng.normal(size=(n,3))+1j*rng.normal(size=(n,3));sr/=np.linalg.norm(sr)
    out=(flat@sr/math.sqrt(L)).reshape(n,da,3)
    reduced=np.einsum('iar,jar->ij',out,out.conj())
    rho=sr@sr.conj().T
    predicted=(gram/L*rho)[perm][:,perm]
    reduced_error=float(np.linalg.norm(reduced-predicted))
    joint_sr=np.einsum('iar,jas->irjs',out,out.conj())
    input_sr=np.einsum('ir,js->irjs',sr,sr.conj())
    expected_sr=(input_sr*(gram/L)[:,None,:,None])[perm,:,:,:][:,:,perm,:]
    gram_error=float(np.linalg.norm(joint_sr-expected_sr))
    assert max(gram_error,reduced_error)<1e-12
    base=[F(1,2),F(1,3),F(1,6)]
    pop=[base[a]*base[b] for a,b in strings]
    pout=[pop[int(perm[j])] for j in range(n)]
    drift=sum((abs(a-b) for a,b in zip(pop,pout)),F(0))/2
    px=sum((pop[s] for s in range(n) if flag[s]),F(0))
    py=sum((pout[s] for s in range(n) if flag[s]),F(0))
    assert drift==px-py==F(4,9)
    ro=np.diag(list(map(float,pop)))
    output=flat@ro@flat.T/L
    rho_aux=partial(output,n,da,1)
    assert np.linalg.matrix_rank(rho_aux,tol=1e-12)<=3
    return {'block_length':1,'block_count':2,'material_dimension':n,
        'ladder_L':L,'auxiliary_dimension_including_record':da,
        'full_energy_permutation_conserves_labels':True,
        'TP_and_unital_integer_Gram_exact':True,
        'unknown_SR_Gram_channel_residual':gram_error,
        'nonthermal_diagonal_drift_exact':str(drift),
        'record_event_before_after':[str(px),str(py)],
        'complete_auxiliary_entropy_diagnostic':entropy(rho_aux)},(perm,gram/L,energy)

def gibbs_bound(data):
    perm,gram,energy=data;rows=[]
    for beta in [-1.2,0.,.7,2.]:
        p=np.exp(-beta*energy);p/=p.sum();q=p[perm]
        ds=float(-np.dot(q,np.log(q))+np.dot(p,np.log(p)))
        de=float(np.dot(q-p,energy));da=-de
        drel=float(np.dot(q,np.log(q/p)))
        dtr=float(np.abs(q-p).sum()/2)
        assert abs(drel-(beta*de-ds))<1e-12
        assert drel<=abs(beta)*abs(da)+1e-12
        assert dtr<=math.sqrt(abs(beta)*abs(da)/2)+1e-12
        if beta==0:assert dtr==0
        rows.append({'beta':beta,'material_trace_distance':dtr,
                     'auxiliary_energy_change':da,'relative_entropy':drel,
                     'uniform_upper':min(1.,math.sqrt(abs(beta)*abs(da)/2))})
    return rows

def clausius_and_subsystem():
    hs=np.diag([0.,1.]);hb=hs.copy();beta=math.log(2)
    tau=np.diag([2/3,1/3]);rho=np.array([[.4,.18j],[-.18j,.6]])
    omega=np.kron(rho,tau)
    swap=np.array([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]],complex)
    phase=np.diag([1,1,-1,-1]);kraus=[math.sqrt(.6)*swap,math.sqrt(.4)*phase]
    assert np.linalg.norm(sum(k@k.conj().T for k in kraus)-np.eye(4))<1e-14
    out=sum(k@omega@k.conj().T for k in kraus)
    a=partial(out,2,2,0);b=partial(out,2,2,1)
    sigma=entropy(out)-entropy(omega)
    mutual=entropy(a)+entropy(b)-entropy(out)
    rel=relative(b,tau)
    de=float(np.trace(hb@(b-tau)).real)
    lhs=entropy(a)-entropy(rho)+beta*de
    assert abs(lhs-(sigma+mutual+rel))<1e-12
    assert min(sigma,mutual,rel,lhs)>-1e-12
    # Global SWAP unital, thermal partial channel rho->tau not unital.
    reduced_identity=partial(swap@np.kron(np.eye(2)/2,tau)@swap,2,2,0)
    defect=trace_distance(reduced_identity,np.eye(2)/2)
    assert abs(defect-1/6)<1e-14
    return {'full_AB_unital':True,'target_only_normalized_defect':defect,
            'Clausius_LHS':lhs,'total_entropy_increase':sigma,
            'final_mutual_information':mutual,'final_bath_relative_entropy':rel,
            'identity_residual':abs(lhs-sigma-mutual-rel)}

def sharp_faithfulness():
    D=6;rank=2;lam=F(1,24);rows=[]
    initial=[lam]*D;initial[0]+=1-D*lam
    assert sum(initial)==1 and min(initial)==lam
    for amount in [F(0),F(1,5),F(1)]:
        replacement=[F(1,rank) if i<rank else F(0) for i in range(D)]
        out=[(1-amount)*p+amount*r for p,r in zip(initial,replacement)]
        image_uniform=[(1-amount)/D+amount*r for r in replacement]
        epsilon=sum((abs(p-F(1,D)) for p in image_uniform),F(0))/2
        failure=sum(out[rank:],F(0))
        lower=lam*D*(1-F(rank,D)-epsilon)
        assert failure==lower
        rows.append({'reset_mixture':str(amount),'lambda':str(lam),
                     'normalized_unital_defect':str(epsilon),
                     'failure_equals_lower_bound':str(failure)})
    return {'full_dimension':D,'target_rank':rank,'exact_saturation':rows}

def native_bridge():
    cold=list(map(F,json.loads((STAGE/'1202/results.json').read_text(encoding='utf-8'))['one_marginal_rational']))
    hot=list(map(F,json.loads((STAGE/'1204/results.json').read_text(encoding='utf-8'))['new_reference']['one_marginal_rational']))
    lam=min(cold)*min(hot)
    h=np.array([[0,-.2,0,-math.sqrt(3)/80],[-.2,1,0,0],[0,0,1,0],[-math.sqrt(3)/80,0,0,-.025]])
    q=np.zeros((4,4));q[1,2]=q[2,1]=1
    en,ev=np.linalg.eigh(h);ident=np.eye(4)
    h0=np.kron(h,ident)+np.kron(ident,h)
    full=h0+.1*np.kron(q,q)
    e,u=np.linalg.eigh(full)
    projectors=[np.kron(ev[:,j:j+1]@ev[:,j:j+1].T,ident) for j in range(4)]
    ground=projectors[0]
    def evolve(H,time):
        a,v=np.linalg.eigh(H);return (v*np.exp(-1j*time*a))@v.conj().T
    T=F(606);mass=F(10**24);width=F(1,10**10)
    bt=width*(1+T*T/(8*mass*mass*width**4))
    rows=[]
    for speed in [1,2]:
        duration=600/speed
        pre=evolve(h0,1/speed);middle=(u*np.exp(-1j*duration*e))@u.conj().T
        post=evolve(h0,606-601/speed)
        operators=[post@P@middle@pre for P in projectors]
        tp=float(np.linalg.norm(sum(k.conj().T@k for k in operators)-np.eye(16)))
        uni=float(np.linalg.norm(sum(k@k.conj().T for k in operators)-np.eye(16)))
        assert max(tp,uni)<2e-13
        eps=(150+F(3,2)*speed*F(22,7))*T*bt+F(1,10000*speed)
        bound=lam*(12-16*eps)
        assert bound>F(2693,10000)
        outputs=[]
        for orientation,left,right in [('CH',cold,hot),('HC',hot,cold)]:
            ra=ev@np.diag(list(map(float,left)))@ev.T
            rb=ev@np.diag(list(map(float,right)))@ev.T
            omega=np.kron(ra,rb)
            out=sum(k@omega@k.conj().T for k in operators)
            failure=1-float(np.trace(ground@out).real)
            assert failure+1e-12>=float(12*lam)
            outputs.append({'orientation':orientation,'ideal_failure_diagnostic':failure,
                            'ideal_total_entropy_gain_diagnostic':entropy(out)-entropy(omega)})
        rows.append({'speed':speed,'QQ_duration':duration,'epsilon_exact':str(eps),
                     'epsilon_diagnostic':float(eps),'TP_residual':tp,'unital_residual':uni,
                     'actual_failure_lower_exact':str(bound),
                     'actual_failure_lower_diagnostic':float(bound),'ideal_outputs':outputs})
    return {'lambda_exact':str(lam),'lambda_diagnostic':float(lam),
            'cold_min_index':cold.index(min(cold)),'hot_min_index':hot.index(min(hot)),
            'ideal_failure_lower_exact':str(12*lam),'rows':rows,
            'scope':'1203 arbitrary ABR, fixed independent C/K; not the correlated 1205 preparation section.'}

def deletion():
    d=4;ground=np.diag([1.,0,0,0])
    reset_defect=trace_distance(ground,np.eye(d)/d)
    assert reset_defect==.75
    # If the pure inventory is included as material, the complete channel is SWAP,
    # but its actual initial material state has zero minimum eigenvalue.
    thermal=np.diag([.4,.3,.2,.1]);complete=np.kron(thermal,ground)
    assert np.linalg.eigvalsh(complete).min()==0
    swap=np.zeros((d*d,d*d))
    for a in range(d):
        for b in range(d):swap[b*d+a,a*d+b]=1
    assert np.array_equal(swap@swap.T,np.eye(d*d))
    out=swap@complete@swap.T
    assert np.linalg.norm(partial(out,d,d,0)-ground)<1e-14
    return {'1213_reset_normalized_unital_defect_exact':'3/4',
            'pure_inventory_included_material_min_eigenvalue':0,
            'pure_inventory_as_aux_target_reset_defect_exact':'3/4',
            'whole_material_SWAP_unital':True,'exact_reset_target_probability':1}

def derive():
    witness,data=first_hit()
    groups={'first_hit_lift_with_complete_pointer':witness,
            'positive_zero_negative_beta':gibbs_bound(data),
            'full_material_Clausius_and_target_distinction':clausius_and_subsystem(),
            'faithfulness_and_normalized_unital_defect':sharp_faithfulness(),
            'native_independent_contact_bridge':native_bridge(),
            'reset_and_pure_inventory_deletions':deletion()}
    sources=['1202/results.json','1203/proof.md','1204/results.json','1209/proof.md','1213/proof.md']
    return {'schema':'1214-independent-v1','status':'PASS','all_checks_passed':True,
            'check_groups':len(groups),'groups':groups,
            'source_sha256':{p:sha(STAGE/p) for p in sources},
            'independent_code_sha256':sha(Path(__file__)),
            'not_certified':['No all-n classification is proved by finite examples.',
             'Unital role is not derived from FUCP/CC or imposed on all work processes.',
             'Negative beta remains in the mathematical Gibbs class.',
             'The window references remain non-Gibbs; their initial relative-entropy defects remain.',
             'All auxiliary budgets must be assessed jointly; channel composition alone does not add entropy-return budgets.',
             'No new controller, preparation, physical area, temperature instrument, or universal third law is certified.']}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    target=HERE/'independent_results.json'
    if args.write and target.exists():raise SystemExit('Refusing --write: results already exist.')
    result=derive()
    if args.write:
        with target.open('x',encoding='utf-8') as out:json.dump(result,out,indent=2,sort_keys=True,ensure_ascii=False);out.write('\n')
    saved=json.loads(target.read_text(encoding='utf-8'))
    if saved!=result:raise SystemExit('Read-only comparison failed; no file overwritten.')
    print(json.dumps({'status':'PASS','all_checks_passed':True,'results_match':True,'check_groups':len(result['groups'])}))

if __name__=='__main__':main()
