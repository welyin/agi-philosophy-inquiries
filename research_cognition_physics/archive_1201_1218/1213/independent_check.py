"""Independent 1213 audit. Default: derive, read-only compare, print a compact PASS JSON.
No author module is imported. --write is the only result-writing mode.
Finite clock tests are supported-input regressions, not a continuous-spectrum proof.
"""
import argparse
import hashlib
import json
import math
from fractions import Fraction as F
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def poly(x):
    return 32000*x**3 - 31200*x**2 - 2095*x - 17

def interval_add(x, y):
    return (x[0]+y[0], x[1]+y[1])

def interval_mul(x, y):
    vals = [a*b for a in x for b in y]
    return min(vals), max(vals)

def interval_scale(q, x):
    return interval_mul((q, q), x)

def decimal_interval(x):
    return {'lower_rational':str(x[0]), 'upper_rational':str(x[1]),
            'midpoint_diagnostic':float((x[0]+x[1])/2)}

def inner(f, g):
    # Piecewise constant, real two-branch energy wavefunctions.
    out = F(0)
    for s,a,b,u in f:
        for t,c,d,v in g:
            if s == t:
                out += max(F(0), min(b,d)-max(a,c))*u*v
    return out

def shift(f, delta):
    return [(s,a+delta,b+delta,u) for s,a,b,u in f]

def adjoint_shift(f, delta):
    return [(s,max(a,delta)-delta,b-delta,u)
            for s,a,b,u in f if b > delta]

def restrict(f, delta):
    return [(s,max(a,delta),b,u) for s,a,b,u in f if b > delta]

def moment(f, power):
    # f uses disjoint intervals on each branch.
    return sum((u*u*(b**(power+1)-a**(power+1))/F(power+1)
                for s,a,b,u in f), F(0))

def equal_l2(f,g):
    return inner(f,f)+inner(g,g)-2*inner(f,g) == 0

def packet_tests():
    f = [(1,F(0),F(1,3),F(2)), (1,F(1),F(5,3),F(-1)),
         (-1,F(1,4),F(3,4),F(3))]
    rows=[]
    for delta in map(F,['0','1/8','1/2','1','2']):
        tf=shift(f,delta)
        assert inner(tf,tf)==inner(f,f)
        assert equal_l2(adjoint_shift(tf,delta), f)
        assert equal_l2(shift(adjoint_shift(f,delta),delta),restrict(f,delta))
        assert moment(tf,1)==moment(f,1)+delta*moment(f,0)
        assert moment(tf,2)==moment(f,2)+2*delta*moment(f,1)+delta**2*moment(f,0)
        rows.append({'shift':str(delta),'norm_squared':str(inner(tf,tf)),
                     'shifted_first_moment':str(moment(tf,1)),
                     'shifted_second_moment':str(moment(tf,2)),
                     'range_projection_norm_squared':str(inner(restrict(f,delta),restrict(f,delta)))})
    return rows

def entropy(rho):
    vals=np.linalg.eigvalsh((rho+rho.conj().T)/2)
    assert vals.min()>-1e-12
    vals=vals[vals>1e-14]
    return float(-np.dot(vals,np.log(vals)))

def finite_clock_test(gaps, last, support):
    # Explicit finite discrete analogue; no claim of a global truncated T isometry.
    d=4; nc=last+1; nk=4
    labels=[(s,c,k) for s in range(d) for c in range(nc) for k in range(nk)]
    index={label:i for i,label in enumerate(labels)}
    perm=np.arange(len(labels)); pairs=[]
    for s in range(1,d):
        for c in range(nc):
            if c+gaps[s] < nc:
                i=index[(s,c,0)];j=index[(0,c+gaps[s],s)]
                assert perm[i]==i and perm[j]==j
                perm[i]=j;perm[j]=i;pairs.append((i,j))
    q=np.array([gaps[s]+c for s,c,k in labels],dtype=np.int64)
    assert np.array_equal(perm[perm],np.arange(len(labels)))
    assert np.array_equal(q[perm],q)
    # Each maximal-clock excited state is an explicit boundary failure.
    boundary=index[(3,last,0)]
    assert perm[boundary]==boundary
    assert all(c+max(gaps)<=last for c in support)
    phi=np.zeros(nc,dtype=complex)
    for j,c in enumerate(support):phi[c]=complex(j+1,1-j)
    phi/=np.linalg.norm(phi)
    nr=5
    sr=np.array([[complex((s+1)*(r+2)%7-3,(2*s+r)%5-2)
                  for r in range(nr)] for s in range(d)],dtype=complex)
    sr/=np.linalg.norm(sr)
    inp=np.zeros((d,nc,nk,nr),complex)
    inp[:,:,0,:]=sr[:,None,:]*phi[None,:,None]
    arr=inp.reshape(len(labels),nr)
    tau=.73
    phase=np.exp(-1j*tau*q)
    actual=phase[:,None]*arr[perm]
    expected=np.zeros_like(inp)
    shifted=[]
    for s in range(d):
        psi=np.zeros(nc,complex)
        for c in support:psi[c+gaps[s]]=phi[c]
        shifted.append(psi)
        expected[0,:,s,:]=np.exp(-1j*tau*np.arange(nc))[:,None]*psi[:,None]*sr[s][None,:]
    error=float(np.linalg.norm(actual-expected.reshape(len(labels),nr)))
    assert error<1e-13
    out=actual.reshape(d,nc,nk,nr)
    ground=float(np.sum(np.abs(out[0])**2))
    assert abs(ground-1)<1e-13
    rho_s=sr@sr.conj().T
    ck=out[0].reshape(nc*nk,nr)
    rho_ck=ck@ck.conj().T
    entropy_error=abs(entropy(rho_s)-entropy(rho_ck))
    assert entropy_error<1e-12
    # Reduced K keeps the overlapping shifted-packet coherences.
    rho_k=np.einsum('ckr,clr->kl',out[0],out[0].conj())
    rho_k_expected=np.empty((nk,nk),complex)
    for i in range(d):
        for j in range(d):
            rho_k_expected[i,j]=rho_s[i,j]*np.vdot(shifted[j],shifted[i])
    kerror=float(np.linalg.norm(rho_k-rho_k_expected))
    assert kerror<1e-13
    # Intermediate times include the resident G and Q; verify both energy accounts.
    p0=float(rho_s[0,0].real)
    mean_g_expected=math.pi*(1-p0)/(2*tau)
    mean_q_initial=float(np.sum(q[:,None]*abs(arr)**2))
    ledger=[]
    for ratio in [0,.25,.5,1,2]:
        t=tau*ratio;z=np.exp(-1j*math.pi*ratio)
        a=(1+z)/2;b=(1-z)/2
        state=np.exp(-1j*t*q)[:,None]*(a*arr+b*arr[perm])
        gmean=float((math.pi/(2*tau)*np.vdot(state,state-state[perm])).real)
        qmean=float(np.sum(q[:,None]*abs(state)**2))
        assert abs(gmean-mean_g_expected)<1e-12
        assert abs(qmean-mean_q_initial)<1e-12
        ground_now=float(np.sum(abs(state.reshape(d,nc,nk,nr)[0])**2))
        if ratio==2:assert abs(ground_now-p0)<1e-12
        ledger.append({'t_over_tau':ratio,'Q_mean':qmean,'G_mean':gmean,
                       'ground_probability':ground_now})
    # Initial p times a diagonal clock, mixed K: exact Fraction deletion check.
    p=[F(1,10),F(2,10),F(3,10),F(4,10)]
    wc={c:F(j+1,sum(range(1,len(support)+1))) for j,c in enumerate(support)}
    probability=F(0)
    for s in range(d):
        for c,weight in wc.items():
            for k in range(nk):
                sout=labels[int(perm[index[(s,c,k)]])][0]
                if sout==0:probability+=p[s]*weight/F(4)
    threshold=sum((sum((v for c,v in wc.items() if c<gaps[i]),F(0))
                   for i in range(1,d)),F(0))
    formula=F(1,4)+p[0]*threshold/F(4)
    assert probability==formula<1
    # Erasing the K label would violate an input Gram matrix in this analogue too.
    overlaps=[abs(np.vdot(shifted[i],shifted[j]))
              for i in range(d) for j in range(i+1,d)]
    assert max(overlaps)>0
    return {'target_gaps':gaps,'clock_dimension':nc,'record_dimension':nk,
            'supported_clock_indices':support,'pair_count':len(pairs),
            'exact_involution_and_energy_labels':True,
            'arbitrary_entangled_SR_endpoint_norm_error':error,
            'K_coherence_norm_error':kerror,'CK_entropy_error':entropy_error,
            'mixed_K_ground_probability_rational':str(probability),
            'max_overlap_without_record':float(max(overlaps)),
            'unpaired_upper_boundary_is_not_reset':True,'resident_interaction_ledger':ledger}

def native_ledger():
    eigen=json.loads((STAGE/'1201/results.json').read_text(encoding='utf-8'))['eigenvalues']
    cold=json.loads((STAGE/'1202/results.json').read_text(encoding='utf-8'))['one_marginal_rational']
    p=list(map(F,cold));assert sum(p)==1 and min(p)>0
    margin=F(1,10**60)
    roots=[]
    for i,s in enumerate(eigen):
        x=F(s)
        if i==2:roots.append((x,x));continue
        lo,hi=x-margin,x+margin
        assert poly(lo)*poly(hi)<0
        roots.append((lo,hi))
    assert all(roots[i][1]<roots[i+1][0] for i in range(3))
    gaps=[(F(0),F(0))]+[(lo-roots[0][1],hi-roots[0][0]) for lo,hi in roots[1:]]
    assert gaps[-1][1]<F(1094,1000)
    mu=F(17,24);nu=F(355,576) # M=3,p_c=2,s=1/2.
    # Recompute moments via Gaussian binomial moments, not the displayed formulas.
    raw=[]
    for power in [2,4]:
        total=F(0)
        for k in range(0,power+1,2):
            double_factorial=math.prod(range(1,k,2))
            total+=math.comb(power,k)*F(2)**(power-k)*F(1,2)**k*double_factorial
        raw.append(total)
    assert raw[0]/6==mu and raw[1]/36==nu
    dm=(F(0),F(0));dm2=dm
    for prob,gap in zip(p,gaps):
        dm=interval_add(dm,interval_scale(prob,gap))
        dm2=interval_add(dm2,interval_scale(prob,interval_mul(gap,gap)))
    mout=interval_add((nu,nu),interval_add(interval_scale(2*mu,dm),dm2))
    upper=nu+2*mu*gaps[-1][1]+gaps[-1][1]**2
    assert mout[1]<=upper
    return {'native_root_intervals_certified_by_exact_cubic_sign':True,
            'root_radius_rational':str(margin),'gap_width':decimal_interval(gaps[-1]),
            'cold_ground_probability_rational':str(p[0]),
            'gaussian_parameters':{'M':3,'p_c':2,'s':'1/2'},
            'gaussian_initial_H_mean_rational':str(mu),
            'gaussian_initial_H2_mean_rational':str(nu),
            'cold_energy_to_clock':decimal_interval(dm),
            'cold_final_clock_H2':decimal_interval(mout),
            'cold_G_mean_coefficient_of_pi_over_2tau':str(1-p[0]),
            'cold_mixed_K_ground_probability_strict_upper_rational':str(F(1,4)+3*p[0]/4),
            'cold_CK_entropy_diagnostic':float(-sum(float(x)*math.log(float(x)) for x in p))}

def derive():
    packet=packet_tests()
    finite=[finite_clock_test([0,1,3,4],8,[0,1,2]),
            finite_clock_test([0,2,3,5],11,[0,1,2,3])]
    native=native_ledger()
    amplitude=F(-1,4);density=2*amplitude;derivative=2*(amplitude-1)
    assert density> -1 and derivative< -1
    return {'schema':'1213-independent-audit-v1','status':'PASS',
            'all_checks_passed':True,'check_groups':6,
            'source_sha256':{p:sha(STAGE/p) for p in ['1201/results.json','1202/results.json']},
            'independent_code_sha256':sha(Path(__file__)),
            'half_line_packet_exact_fraction_tests':packet,
            'finite_supported_input_regressions':finite,
            'native_exact_interval_energy_ledger':native,
            'gaussian_threshold_exponents':{'amplitude':str(amplitude),
                'density':str(density),'squared_derivative':str(derivative),
                'normalizable_but_not_DX':True},
            'scope':['Finite clock regressions do not establish a global truncated energy-shift isometry.',
                     'Continuous-spectrum strong commutation and domains are reviewed analytically.',
                     'No original QQ, local geometry, physical trigger, released cold product, or all-resource bound is certified.',
                     'Initial pure clock and blank record are resources, not generated preparations.']}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true')
    args=parser.parse_args();result=derive();target=HERE/'independent_results.json'
    output=json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True)+'\n'
    if args.write:
        if target.exists():raise SystemExit('--write requires an absent results file; existing results were not modified.')
        target.write_text(output,encoding='utf-8')
    saved=json.loads(target.read_text(encoding='utf-8'))
    if saved!=result:raise SystemExit('Independent result differs; no files modified.')
    print(json.dumps({'status':'PASS','all_checks_passed':True,'results_match':True,'check_groups':6}))

if __name__=='__main__':main()
