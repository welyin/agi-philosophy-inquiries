"""Independent 1209 finite resource checks; no imports of author code.

Default is read-only and compares the saved independent result. --write saves it.
Exact arithmetic certifies permutations, probabilities and large-instance bounds.
NumPy checks reduced states, entropy and finite isometries, not universal claims.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
import argparse
import hashlib
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent


def entropy(a):
    w = np.linalg.eigvalsh((a + a.conj().T) / 2)
    assert w.min() > -2e-12
    return float(-sum(x * math.log(x) for x in w if x > 1e-14))


def trace_distance(a, b):
    return float(np.abs(np.linalg.eigvalsh((a-b + (a-b).conj().T)/2)).sum()/2)


def reduced(a, ds, da):
    t = a.reshape(ds, da, ds, da)
    return np.einsum('iaja->ij', t), np.einsum('iaib->ab', t)


def tidy(x):
    return float(f'{x:.12g}')


def dec(x, places=24, upper=False):
    scale = 10**places
    n = -((-x.numerator*scale)//x.denominator) if upper else x.numerator*scale//x.denominator
    return ('-' if n < 0 else '') + str(abs(n)//scale) + '.' + str(abs(n)%scale).zfill(places)


def log_upper(x, terms=28):
    """Strict rational upper bound for log(x), x>=1, via binary reduction."""
    assert x >= 1
    m = 0
    while x >= 2:
        x /= 2
        m += 1

    def bound(y):
        z = (y-1)/(y+1)
        s = 2*sum((z**(2*k+1)/F(2*k+1) for k in range(terms)), F(0))
        tail = 2*z**(2*terms+1)/((2*terms+1)*(1-z*z))
        return s+tail
    return m*bound(F(2))+bound(x)


def first_hit_case(blocks, length, p, record=False, beta=None, autonomous=False):
    """One-symbol x=0,y=1, inactive=2; E=(0,1,2), energy gap one."""
    words = list(product(range(3), repeat=blocks))
    index = {s:i for i,s in enumerate(words)}
    ds, db, dp = len(words), length+2, (2 if record else 1)
    da, dim = db*dp, len(words)*db*dp
    energy_s = [sum(s) for s in words]
    probabilities = [math.prod(p[j] for j in s) for s in words]
    target, hit, low = [], [], []
    for s in words:
        places = [j for j,v in enumerate(s) if v in (0,1)]
        if not places:
            target.append(index[s]); hit.append(False); low.append(False)
        else:
            j = places[0]
            t = list(s); t[j] = 1-t[j]
            target.append(index[tuple(t)]); hit.append(True); low.append(s[j] == 0)
    assert all(target[target[s]] == s for s in range(ds))

    def ix(s,k,b=0):
        return (s*db+k)*dp+b

    perm = list(range(dim))
    # Low system energy with battery k has the same total as high with k-1.
    for s in range(ds):
        if not hit[s] or not low[s]:
            continue
        t = target[s]
        assert energy_s[t] == energy_s[s]+1
        for k in range(1, length+2):
            for b in range(dp):
                a,c = ix(s,k,b),ix(t,k-1,b)
                perm[a],perm[c] = c,a
    if record:
        # Full controlled-X extension; output target first-x event is copied.
        copy = list(range(dim))
        for s in range(ds):
            if hit[s] and low[s]:
                for k in range(db):
                    copy[ix(s,k,0)],copy[ix(s,k,1)] = ix(s,k,1),ix(s,k,0)
        perm = [copy[j] for j in perm]
    assert sorted(perm) == list(range(dim))
    qenergy = [energy_s[s]+k for s in range(ds) for k in range(db) for _ in range(dp)]
    assert all(qenergy[i] == qenergy[perm[i]] for i in range(dim))

    battery = np.zeros(db); battery[1:length+1] = 1/math.sqrt(length)
    aux = np.zeros(da); aux[::dp] = battery
    J = np.zeros((dim,ds)); J0 = np.zeros_like(J)
    for s in range(ds):
        out_pointer = int(record and hit[target[s]] and low[target[s]])
        for k in range(1,length+1):
            J[perm[ix(s,k)],s] = 1/math.sqrt(length)
            J0[ix(target[s],k,out_pointer),s] = 1/math.sqrt(length)
    assert np.linalg.norm(J.T@J-np.eye(ds)) < 2e-14
    iso = np.linalg.norm(J-J0,2)
    assert abs(iso-math.sqrt(2/length)) < 2e-14
    pf = np.array([float(x) for x in probabilities])
    rho = (J*pf)@J.T
    rs, ra = reduced(rho,ds,da)
    targetp = [F(0) for _ in range(ds)]
    for s in range(ds):
        targetp[target[s]] += probabilities[s]
    assert np.linalg.norm(rs-np.diag([float(v) for v in targetp])) < 2e-14

    # Independent exact marginal accumulation from basis paths, no eigensolver.
    exact_a = [[F(0) for _ in range(da)] for _ in range(da)]
    for s in range(ds):
        paths = [divmod(perm[ix(s,k)],da) for k in range(1,length+1)]
        assert {out_s for out_s,_ in paths} == {target[s]}
        for _,a in paths:
            for _,b in paths:
                exact_a[a][b] += probabilities[s]/length
    exact_as_float = np.array([[float(x) for x in row] for row in exact_a])
    assert np.linalg.norm(ra-exact_as_float) < 2e-14
    d_s = sum((abs(targetp[i]-probabilities[i]) for i in range(ds)),F(0))/2
    q = p[0]+p[1]
    formula = (1-(1-q)**blocks)*abs(p[0]-p[1])/q
    assert d_s == formula
    first_x_in = sum((probabilities[s] for s in range(ds) if hit[s] and low[s]), F(0))
    first_x_out = sum((targetp[s] for s in range(ds) if hit[s] and low[s]), F(0))
    assert abs(first_x_out-first_x_in) == d_s
    if record:
        pointer_one = sum((exact_a[k*dp+1][k*dp+1] for k in range(db)), F(0))
        assert pointer_one == first_x_out
    delta_e = sum((F(k)*(exact_a[k*dp+b][k*dp+b]-F(1,length)*int(1<=k<=length and b==0))
                   for k in range(db) for b in range(dp)),F(0))
    delta_s_energy = sum((F(energy_s[s])*(targetp[s]-probabilities[s]) for s in range(ds)),F(0))
    assert delta_e == -delta_s_energy and abs(delta_e) <= 1
    s_aux = entropy(ra)
    rank = int(np.count_nonzero(np.linalg.eigvalsh(ra)>1e-12))
    assert rank <= 3
    s_in = float(-sum(float(x)*math.log(float(x)) for x in probabilities))
    mutual = entropy(rs)+s_aux-entropy(rho)
    assert abs(entropy(rs)-s_in)<2e-12 and abs(mutual-s_aux)<2e-12
    overlap = float(aux@ra@aux)
    if record:
        # Here the low-energy x is the high-population word; record outputs 1
        # only from initially y. Probabilities are tested without a sign shortcut.
        py = (1-(1-q)**blocks)*p[1]/q
        nohit = (1-q)**blocks
        px = (1-(1-q)**blocks)*p[0]/q
        exact_overlap = nohit+px*F(length-1,length)**2
    else:
        exact_overlap = (1-q)**blocks + (1-(1-q)**blocks)*F(length-1,length)**2
    assert abs(overlap-float(exact_overlap))<2e-14
    fail = 1-exact_overlap
    td_aux = trace_distance(ra,np.outer(aux,aux))
    assert td_aux <= math.sqrt(float(fail))+2e-14
    if fail <= F(2,3):
        ff = float(fail)
        bound = -ff*math.log(ff)-(1-ff)*math.log1p(-ff)+ff*math.log(2) if ff else 0
        assert s_aux <= bound+2e-13
    else:
        bound = math.log(3)
    rel = sum(float(targetp[i])*math.log(float(targetp[i]/probabilities[i])) for i in range(ds))
    identity_error = None
    if beta is not None:
        rhs = -beta*float(delta_e)+s_aux-mutual
        identity_error = abs(rel-rhs)
        assert identity_error < 2e-12
        assert rel <= abs(beta)*abs(float(delta_e))+abs(s_aux)+2e-12

    autonomous_result = None
    if autonomous:
        # Hermitian logarithm of the complete permutation, including record.
        G = np.zeros((dim,dim),complex)
        visited = set()
        for start in range(dim):
            if start in visited: continue
            cyc=[]; j=start
            while j not in visited:
                visited.add(j);cyc.append(j);j=perm[j]
            assert j == start
            ell=len(cyc)
            for m in range(ell):
                phase = 0 if m==0 else 2*math.pi*(1-m/ell)
                v=np.zeros(dim,complex)
                for k,j in enumerate(cyc): v[j]=np.exp(-2j*math.pi*m*k/ell)/math.sqrt(ell)
                G += phase*np.outer(v,v.conj())
        Qtotal=np.diag(qenergy)
        assert np.linalg.norm(G@Qtotal-Qtotal@G)<2e-12
        energies,vec=np.linalg.eigh(Qtotal+G)
        actual=(vec*np.exp(-1j*energies))@vec.conj().T
        W=np.zeros((dim,dim));W[perm,np.arange(dim)]=1
        expected=np.exp(-1j*np.array(qenergy))[:,None]*W
        error=np.linalg.norm(actual-expected,2)
        assert error<1e-12 and np.linalg.eigvalsh(G).min()>-1e-12
        rin=np.kron(np.diag(pf),np.outer(aux,aux))
        rout=actual@rin@actual.conj().T
        dg=float(np.trace(G@(rout-rin)).real)
        assert abs(dg)<2e-12
        autonomous_result={
            'full_log_instead_of_involution_log':True,
            'norm_G':tidy(np.linalg.norm(G,2)),
            'free_phase_identity_error':tidy(error),
            'interaction_energy_change':tidy(dg),
            'endpoint_interaction_is_not_switched_off':True}
    return {
        'blocks':blocks,'L':length,'pointer_included':record,'joint_dimension':dim,
        'exact_material_distance':str(d_s),'exact_auxiliary_energy_change':str(delta_e),
        'first_x_read_attains_material_distance_exactly':True,
        'auxiliary_entropy':tidy(s_aux),'auxiliary_rank':rank,
        'auxiliary_trace_distance':tidy(td_aux),'auxiliary_initial_overlap':tidy(overlap),
        'rank3_entropy_upper':tidy(bound),'mutual_information':tidy(mutual),
        'full_isometry_error':tidy(iso),'gibbs_identity_error':None if identity_error is None else tidy(identity_error),
        'autonomous':autonomous_result}


def generic_gibbs_identity():
    records=[]
    rng=np.random.default_rng(1209003)
    es=np.arange(3,dtype=float);ea=np.arange(4,dtype=float)
    et=(es[:,None]+ea[None,:]).ravel()
    W=np.eye(12,dtype=complex)
    for e in sorted(set(et)):
        ids=np.where(et==e)[0]
        z=rng.normal(size=(len(ids),len(ids)))+1j*rng.normal(size=(len(ids),len(ids)))
        q,r=np.linalg.qr(z);W[np.ix_(ids,ids)]=q
    assert np.linalg.norm(W.conj().T@np.diag(et)@W-np.diag(et))<2e-14
    b=np.array([.4,.3,.2,.1]);b=np.diag(b)
    for beta in (-.7,0.,.8):
        p=np.exp(-beta*es);p/=sum(p)
        rin=np.kron(np.diag(p),b);rout=W@rin@W.conj().T
        rs,ra=reduced(rout,3,4)
        de=float(np.trace(np.diag(ea)@(ra-b)).real)
        ds=entropy(ra)-entropy(b)
        mutual=entropy(rs)+entropy(ra)-entropy(rout)
        rel=-entropy(rs)-float(np.dot(np.diag(rs).real,np.log(p)))
        error=abs(rel-(-beta*de+ds-mutual))
        assert error<2e-12
        bound=abs(beta)*abs(de)+abs(ds)
        td=trace_distance(rs,np.diag(p))
        assert rel<=bound+2e-12 and 2*td*td<=rel+2e-12
        records.append({'beta':beta,'delta_E_aux':tidy(de),'delta_S_aux':tidy(ds),
                        'relative_entropy':tidy(rel),'mutual_information':tidy(mutual),
                        'identity_error':tidy(error),'pinsker_slack':tidy(rel-2*td*td)})
    return records


def one_sided_entropy_deletion():
    q=(math.sqrt(97)-1)/12
    tau=np.array([1,q,q*q]);tau/=sum(tau)
    r=np.array([.4,.4,.2]);e=np.arange(3)
    de=float(np.dot(e,r-tau))
    ds=entropy(np.diag(r))-entropy(np.diag(tau))
    td=float(abs(r-tau).sum()/2)
    assert abs(de)<1e-14 and ds<-.01 and td>.07
    return {'same_mean_energy_residual':tidy(abs(de)),
            'swap_auxiliary_entropy_change':tidy(ds),'swap_auxiliary_state_distance':tidy(td),
            'one_sided_entropy_contract_allows_zero_energy_zero_upper_entropy_budget':True,
            'absolute_entropy_and_state_return_contract_rejects_this_swap':True}


def native_bounds():
    old=json.loads((HERE.parent/'1208/results.json').read_text(encoding='utf-8'))
    prep=json.loads((HERE.parent/'1202/results.json').read_text(encoding='utf-8'))
    source=old['native_parity_witness']
    L=10**6
    qfail=F(2,L)-F(1,L*L)
    # Layer 1: preserve the accepted 1208 witness; no large matrix construction.
    eta_hi=F(source['energy_interval'][1])
    s_hi=F(2,L)*(1+log_upper(F(L)))
    eps=F('.001415')
    assert eps*eps>qfail and s_hi<F('.00003')
    lower=F(source['trace_distance_strict_lower_exact'])
    assert lower>F('.47')
    layer1={
        'L':L,'auxiliary_dimension':L+2,'pointer_included':False,
        'blocks':source['blocks'],'material_count':source['total_materials'],
        'target_distance_strict_lower':dec(lower),
        'absolute_energy_change_strict_upper':dec(eta_hi,upper=True),
        'absolute_entropy_change_strict_upper':dec(s_hi,upper=True),
        'auxiliary_state_distance_strict_upper':str(eps),
        'initial_auxiliary_energy_upper':dec(eta_hi*F(L+1,2),upper=True),
        'initial_auxiliary_energy_variance_upper':dec(eta_hi**2*F(L*L-1,12),upper=True)}

    # Layer 2: the original 1207 counts, both multiplied by eight.
    p=[F(v) for v in prep['one_marginal_rational']]
    x=[0,120,0,80];y=[112,0,88,0]
    assert sum(x)==sum(y)==200
    a=math.prod(p[i]**x[i] for i in range(4));b=math.prod(p[i]**y[i] for i in range(4))
    assert a>b>0
    pair=a+b;asym=(a-b)/pair
    M=(pair.denominator+pair.numerator-1)//pair.numerator
    assert (M-1)*pair<1<=M*pair
    # exp(1)>163/60 gives 1-exp(-1)>103/163.
    assert sum((F(1,math.factorial(k)) for k in range(7)),F(0))>F(163,60)
    hit_lower=F(103,163)
    dist_lower=hit_lower*asym
    e_raw=F(2,L)+b/pair
    e_upper=F(dec(e_raw,places=16,upper=True))
    entropy_upper=e_upper*(1+log_upper(2/e_upper))
    eps2=F('.022')
    # Energy vector is four times the accepted doubled-z vector.
    gap_lo=4*F(source['energy_interval'][0]);gap_hi=4*eta_hi
    assert dist_lower>F('.63') and 0<gap_lo<gap_hi<F('.00339')
    assert e_upper<eps2**2 and entropy_upper<F('.0045')
    layer2={
        'L':L,'auxiliary_dimension':2*(L+2),'pointer_included':True,
        'word_x_counts':x,'word_y_counts':y,'block_length':200,
        'blocks':str(M),'material_count':str(200*M),
        'pair_probability_approx':format(float(pair),'.12e'),
        'low_probability_conditional_upper':dec(b/pair,upper=True),
        'auxiliary_overlap_loss_upper':dec(e_upper,upper=True),
        'auxiliary_rank_upper':3,
        'target_distance_strict_lower':dec(dist_lower),
        'absolute_energy_change_strict_upper':dec(gap_hi,upper=True),
        'absolute_entropy_change_strict_upper':dec(entropy_upper,upper=True),
        'auxiliary_state_distance_strict_upper':str(eps2),
        'initial_auxiliary_energy_upper':dec(gap_hi*F(L+1,2),upper=True),
        'initial_auxiliary_energy_variance_upper':dec(gap_hi**2*F(L*L-1,12),upper=True),
        'blank_pointer_initial_entropy':0,
        'unknown_input_target_is_T_then_record_not_T_tensor_blank':True}
    return {'battery_only_1208':layer1,'complete_record_first_hit':layer2,
            'uniform_initial_resource_bound_claimed':False,
            'native_QQ_or_geometric_implementation_claimed':False,
            'huge_matrix_enumeration_performed':False}


def compute():
    source_paths=['1202/results.json','1207/proof.md','1208/results.json','1208/acceptance.json']
    hashes={p:hashlib.sha256((HERE.parent/p).read_bytes()).hexdigest() for p in source_paths}
    checks=[]
    nonthermal=(F(1,2),F(1,3),F(1,6))
    for n in (1,2):
        for L in (3,5,9):
            checks.append(first_hit_case(n,L,nonthermal,record=False))
            checks.append(first_hit_case(n,L,nonthermal,record=True,autonomous=(n==1 and L==3)))
    for p,beta in (((F(4,7),F(2,7),F(1,7)),math.log(2)),
                   ((F(1,7),F(2,7),F(4,7)),-math.log(2)),
                   ((F(1,3),)*3,0.)):
        checks.append(first_hit_case(2,5,p,record=True,beta=beta))
    return {'schema':'round1209_independent_three_budget_v1',
            'status':'INDEPENDENT_CALCULATIONS_PASS_NOT_ACCEPTANCE',
            'source_sha256':hashes,
            'independent_code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'finite_exact_and_numpy_cases':checks,
            'generic_energy_conserving_gibbs_identity':generic_gibbs_identity(),
            'one_sided_entropy_deletion':one_sided_entropy_deletion(),
            'native_finite_bounds':native_bounds(),
            'all_checks_passed':True,
            'scope':{'analytic_classification_proved_by_this_code':False,
                     'full_auxiliary_state_and_absolute_entropy_budgets_checked':True,
                     'all_auxiliary_records_retained_in_record_examples':True,
                     'existing_QQ_control_permission_inferred':False}}


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true')
    args=parser.parse_args();result=compute();dest=HERE/'independent_results.json'
    if args.write:
        dest.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:
        assert json.loads(dest.read_text(encoding='utf-8'))==result, 'Saved independent result differs.'
    print(json.dumps({'status':'PASS','finite_cases':len(result['finite_exact_and_numpy_cases']),
                      'native_bounds':result['native_finite_bounds'],
                      'mode':'write' if args.write else 'read_only'},ensure_ascii=False,indent=2))
