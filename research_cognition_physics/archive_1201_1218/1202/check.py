"""1202: finite energy-window preparation, kept records and native QQ predictions."""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal as D,localcontext
from bisect import bisect_right
import hashlib,json,math,random
import numpy as np
HERE=Path(__file__).resolve().parent; OLD=HERE.parent/'1201'
N=100; DEN=10**18
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def poly(x):return 32000*x**3-31200*x*x*DEN-2095*x*DEN**2-17*DEN**3
def bracket(a,b):
    a,b=a*(DEN//10**6),b*(DEN//10**6);s=poly(a)
    assert s*poly(b)<0
    while b-a>1:
        c=(a+b)//2
        if s*poly(c)>0:a,s=c,poly(c)
        else:b=c
    return a,b
def setup():
    boxes=[bracket(-54072,-54071),bracket(-9461,-9460),(DEN,DEN),bracket(1038532,1038533)]
    fac=[math.factorial(i) for i in range(N+1)]
    rows=[];total=0;one=[0]*4;two=[[0]*4 for _ in range(4)];margin=10**30
    for a in range(N+1):
        for b in range(N-a+1):
            for c in range(N-a-b+1):
                n=(a,b,c,N-a-b-c)
                lo=sum(n[j]*boxes[j][0] for j in range(4));hi=sum(n[j]*boxes[j][1] for j in range(4))
                inside=2*lo>=39*DEN and 2*hi<=41*DEN
                outside=2*hi<39*DEN or 2*lo>41*DEN
                assert inside or outside
                margin=min(margin,*[abs(2*x-y*DEN) for x in (lo,hi) for y in (39,41)])
                if not inside:continue
                w=fac[N]//math.prod(fac[x] for x in n)
                rows.append((n,w,total));total+=w
                for i in range(4):
                    one[i]+=w*n[i]
                    for j in range(4):two[i][j]+=w*n[i]*(n[j]-(i==j))
    r=[F(x,N*total) for x in one];pair=[[F(x,N*(N-1)*total) for x in row] for row in two]
    assert len(rows)==1728 and sum(r)==1
    assert all(sum(pair[i])==r[i] for i in range(4))
    return boxes,rows,total,r,pair,F(margin,2*DEN)
def mult(n):return math.factorial(sum(n))//math.prod(math.factorial(v) for v in n)
def unrank_permutation(n,k):
    n=list(n);word=[]
    while sum(n):
        for i in range(4):
            if not n[i]:continue
            n[i]-=1;w=mult(n)
            if k<w:word.append(i);break
            k-=w;n[i]+=1
        else:raise AssertionError('rank outside block')
    assert k==0
    return word
def rank_permutation(word):
    n=[word.count(i) for i in range(4)];k=0
    for v in word:
        for i in range(v):
            if n[i]:
                n[i]-=1;k+=mult(n);n[i]+=1
        n[v]-=1
    return k
def mm(a,b):return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]
def tr(a):return list(map(list,zip(*a)))
def mv(a,v):return [sum(x*y for x,y in zip(row,v)) for row in a]
def outer(a,b):return [[x*y for y in b] for x in a]
def tv(a,b):
    if isinstance(a[0],list):return sum(abs(x-y) for ar,br in zip(a,b) for x,y in zip(ar,br))/2
    return sum(abs(x-y) for x,y in zip(a,b))/2
def floats(v):
    return [floats(x) if isinstance(x,list) else float(x) for x in v]
def dec(x):return D(x.numerator)/D(x.denominator)
def ent(p):
    with localcontext() as c:
        c.prec=70
        return float(-sum(dec(x)*dec(x).ln() for x in p if x))
def entropy_matrix(a):
    v=np.linalg.eigvalsh((a+a.conj().T)/2);v=v[v>1e-14]
    return float(-np.sum(v*np.log(v)))
def joint_entropy_diagnostic(pair):
    h=np.array([[0,-.2,0,-np.sqrt(3)/80],[-.2,1,0,0],[0,0,1,0],[-np.sqrt(3)/80,0,0,-.025]])
    q=np.zeros((4,4));q[1,2]=q[2,1]=1
    e,R=np.linalg.eigh(h);H=np.kron(h,np.eye(4))+np.kron(np.eye(4),h)+.1*np.kron(q,q)
    ev,V=np.linalg.eigh(H);basis=np.kron(R,R);channels={}
    for t in [600,300]:
        U=basis.T@((V*np.exp(-1j*t*ev))@V.T)@basis
        a=[];b=[]
        for i in range(4):
            psi=U[:,8+i].reshape(4,4)
            a.append(psi@psi.conj().T);b.append(psi.T@psi.conj())
        channels[t]=(a,b)
    AA=np.zeros((16,16),complex);BB=AA.copy();pairf=np.array(floats(pair))
    for i in range(4):
        for j in range(4):
            AA+=pairf[i,j]*np.kron(channels[600][0][i],channels[300][0][j])
            BB+=pairf[i,j]*np.kron(channels[600][1][i],channels[300][1][j])
    sa,sb=entropy_matrix(AA),entropy_matrix(BB);sin=ent([x for row in pair for x in row])
    Eb=np.kron(np.diag(e),np.eye(4))+np.kron(np.eye(4),np.diag(e))
    initial_E=float(np.sum(pairf*(e[:,None]+e[None,:])));final_E=float(np.trace(BB@Eb).real)
    beta=1.2;logZ=np.log(np.exp(-beta*e).sum())
    din=-sin+beta*initial_E+2*logZ;dout=-sb+beta*final_E+2*logZ
    mutual=sa+sb-sin;lhs=sa+beta*(final_E-initial_E);rhs=mutual+dout-din
    assert abs(lhs-rhs)<1e-12
    return {'S_AA':sa,'S_BB':sb,'S_AB_global_unchanged':sin,'I_AA_BB':mutual,
        'D_initial_pair':din,'D_final_pair':dout,'delta_E_two_sources':final_E-initial_E,
        'clausius_lhs':lhs,'identity_rhs':rhs,'identity_error':abs(lhs-rhs),
        'pointer_not_discarded_in_identity':True,'rigorous_numeric_entropy_certificate':False}
def main():
    boxes,rows,Z,r,pair,margin=setup()
    starts=[row[2] for row in rows];lookup={row[0]:row for row in rows}
    rng=random.Random(1202)
    samples=sorted(set([0,Z-1,*starts,*[s+w-1 for _,w,s in rows],*[rng.randrange(Z) for _ in range(32)]]))
    for idx in samples:
        n,w,start=rows[bisect_right(starts,idx)-1];word=unrank_permutation(n,idx-start)
        assert tuple(word.count(i) for i in range(4))==n
        assert lookup[n][2]+rank_permutation(word)==idx
    # Exhaust every short multiset permutation for the same bijection rule.
    tiny=0
    for a in range(6):
        for b in range(6-a):
            for c in range(6-a-b):
                n=(a,b,c,5-a-b-c);w=mult(n);seen=set()
                for idx in range(w):
                    word=unrank_permutation(n,idx);assert rank_permutation(word)==idx
                    seen.add(tuple(word));tiny+=1
                assert len(seen)==w
    assert tiny==4**5
    k=Z.bit_length();L=10;q=F(Z,2**k);failure=(1-q)**L
    assert k==173 and failure<F('0.000000211')
    # First valid finite register selection is exactly uniform on each rank.
    per_rank=sum((1-q)**j/F(2**k) for j in range(L))
    assert per_rank==(1-failure)/Z
    old=json.loads((OLD/'results.json').read_text(encoding='utf-8'))
    M=[[F(x) for x in row] for row in old['contacts']['600']['K']]
    T=[[F(x) for x in row] for row in old['contacts']['300']['K']]
    joint=mm(mm(M,pair),tr(T));product=outer(mv(M,r),mv(T,r));diff=tv(joint,product)
    event=[(i,j) for i in range(4) for j in range(4) if joint[i][j]>product[i][j]]
    eventgap=sum(joint[i][j]-product[i][j] for i,j in event)
    # Rounded K columns may differ slightly from normalization: safe difference bound.
    matrix_error=F('1e-10')
    joint_entry_error=3*matrix_error
    separation_error=16*joint_entry_error
    assert diff-separation_error>F('.00069502')
    energy={}
    for t in ['600','300']:
        energy[t]={}
        for name in ['delta_EA_coefficients','delta_EB_coefficients','switch_work_coefficients']:
            coeff=list(map(F,old['contacts'][t][name]))
            energy[t][name]={'expectation':float(sum(a*b for a,b in zip(coeff,r))),'numerical_error_upper':'1e-10'}
    with localcontext() as ctx:
        ctx.prec=70
        e=[D(a+b)/D(2*DEN) for a,b in boxes];beta=D('1.2')
        weights=[(-beta*x).exp() for x in e];tau=[F(x/sum(weights)) for x in weights]
        # Exact roots are inside 1e-18 brackets; .5e-18*1.2 exponent changes
        # and correctly rounded Decimal exp are dominated by this common radius.
        tauerr=F('1e-12');assert min(tau)>F('.1')
        local_trace_upper=tv(r,tau)+2*tauerr
        qp=outer(tau,tau);qpair_lower=[[(tau[i]-tauerr)*(tau[j]-tauerr) for j in range(4)] for i in range(4)]
        chi1=sum(r[i]**2/(tau[i]-tauerr) for i in range(4))-1
        chi2=sum(pair[i][j]**2/qpair_lower[i][j] for i in range(4) for j in range(4))-1
        D1=sum(dec(r[i])*(dec(r[i])/dec(tau[i])).ln() for i in range(4))
        D2=sum(dec(pair[i][j])*(dec(pair[i][j])/dec(qp[i][j])).ln() for i in range(4) for j in range(4))
        pool_initial_energy=N*sum(e[i]-e[0] for i in range(4))
        pool_final_energy=N*sum(dec(r[i])*(e[i]-e[0]) for i in range(4))
        exact_energy=sum(dec(r[i])*e[i] for i in range(4))
    assert local_trace_upper<F('.01001') and chi1<F('.000405') and chi2<F('.000912'), (float(local_trace_upper),float(chi1),float(chi2))
    # Independent readout probability errors: local preparation trace distance
    # controls any one-material channel, pair trace distance controls any two.
    pair_product_tv=tv(pair,outer(r,r))
    pair_gibbs_upper=pair_product_tv+2*local_trace_upper
    result={'schema':'energy_window_reference_1202_v1','status':'PASS_CONDITIONAL_PREPARATION_AND_USE',
        'preparation':{'N':N,'energy_window':['39/2','41/2'],'root_denominator':DEN,'root_boxes':boxes,
            'window_boundary_margin_lower':float(margin),'selected_types':len(rows),'dimension':str(Z),
            'seed_bits_each':k,'seed_registers':L,'seed_bits_total':k*L,'seed_success_probability':float(q),
            'all_rejected_probability':float(failure),'declared_all_rejected_upper':'0.000000211',
            'uniform_rank_probability_identity_exact':True,'rank_unrank_checked_indices':len(samples),
            'exhaustive_short_words':tiny,
            'materials_per_pool':{'targets_ground':100,'known_level_inventory':400},
            'postselection_over_random_100_materials_not_used':True,
            'conditional_on_success':True,'all_seed_selection_inventory_records_retained':True},
        'one_marginal_rational':[str(x) for x in r],
        'two_marginal_rational':[[str(x) for x in row] for row in pair],
        'one_marginal':floats(r),'two_marginal':floats(pair),
        'energies':{'initial_inventory_excitation':float(pool_initial_energy),
            'successful_target_excitation':float(pool_final_energy),
            'successful_target_mean_bare_energy':float(exact_energy),
            'inventory_loss_equals_target_gain_exact':True,
            'preparation_switch_net_work':0,
            'sequential_peak_gate_norm':'pi/tau','parallel_peak_gate_norm_upper':'100*pi/tau',
            'controller_implementation_and_cost_completed':False},
        'effective_temperature_comparison':{'beta':'6/5','not_used_to_choose_preparation':True,
            'gibbs_population':floats(tau),'local_trace_distance_upper':float(local_trace_upper),
            'D_one_diagnostic':float(D1),'D_two_diagnostic':float(D2),
            'D_one_certified_chi_square_upper':float(chi1),'D_two_certified_chi_square_upper':float(chi2),
            'two_source_trace_distance_to_product_gibbs_upper':float(pair_gibbs_upper),
            'exact_single_temperature_claimed':False},
        'native_contacts':{'g':'1/10','times':[600,300],'one_outputs':{'600':floats(mv(M,r)),'300':floats(mv(T,r))},
            'joint_outputs':floats(joint),'joint_entry_numerical_error_upper':str(joint_entry_error),
            'product_same_marginals_outputs':floats(product),'correlation_TV_diagnostic':float(diff),
            'correlation_TV_lower':float(diff-separation_error),
            'correlation_event':[list(x) for x in event],'correlation_event_gap':float(eventgap),
            'pair_initial_TV_from_product':float(pair_product_tv),'energy_exchange':energy},
        'entropy':{'pool_marginal_entropy':math.log(Z),'local_entropy':ent(r),
            'pair_entropy':ent([x for row in pair for x in row]),
            'exact_identity':'DeltaS_A+beta*DeltaE_B=I_AB_out+D(B_out||tau^2)-D(B_in||tau^2), for A initially independent of B',
            'certified_two_source_clausius_lower':-float(chi2),
            'initial_A_B_correlation_requires_subtracting_I_initial':True},
        'floating_full_quantum_entropy_diagnostic':joint_entropy_diagnostic(pair),
        'remaining_inputs':['uniform distinguishable seed labels or coherent seed with kept orthogonal record',
            'pure known energy inventory and ground targets','finite reversible rank program',
            'controlled same-material SWAP permission and bounded timed gate schedule',
            'original QQ switch, timing and pointer contracts','actual spatial region and clock identity'],
        'not_claimed':['unconditional success','prepared reference independent of inventory/records',
            'autonomous thermalization or exact Gibbs state','free cooling from maximally mixed stock',
            'full internal clock/controller construction','area entropy, Unruh or GR'],
        'source_sha256':{str((OLD/'results.json').relative_to(HERE.parents[2])):sha(OLD/'results.json'),
            str((OLD/'acceptance.json').relative_to(HERE.parents[2])):sha(OLD/'acceptance.json')}}
    (HERE/'results.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':result['status'],'failed_probability':float(failure),
        'local_trace_upper':float(local_trace_upper),'D1_upper':float(chi1),'D2_upper':float(chi2),
        'joint_corr_lower':float(diff-separation_error),'resources':result['preparation']['materials_per_pool']},ensure_ascii=False))
if __name__=='__main__':main()
