"""1277 author group: exact waiting bounds and inherited finite channel.
Default: read-only recomputation. --write: exclusive first result creation.
No imported historical/independent checker and no new physical permissions.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import math
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in HERE.parents if (p / 'research_cognition_physics').is_dir())
RESULT = HERE / 'results.json'

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def rat(x):
    x = F(x)
    return f'{x.numerator}/{x.denominator}'

def rounded(x):
    return round(float(x), 13)

def peval(p, x):
    return sum(a * x ** i for i, a in enumerate(p))

def qadd(a, b):
    return (a[0] + b[0], a[1] + b[1])

def qmul(a, b):
    return (a[0]*b[0]+2*a[1]*b[1], a[0]*b[1]+a[1]*b[0])

def qeval(poly, x):
    out = (F(0), F(0))
    for a in reversed(poly):
        out = qadd(qmul(out, x), (a, F(0)))
    return out

def trace_distance(a, b):
    d = a-b
    return float(np.abs(np.linalg.eigvalsh((d+d.conj().T)/2)).sum()/2)

def entropy(a):
    v = np.linalg.eigvalsh((a+a.conj().T)/2)
    return float(-sum(x*math.log(x) for x in v if x > 1e-14))

def kraus(u):
    p0, p1, p2 = u*u, 2*u*(1-u), (1-u)**2
    k0, k1 = np.zeros((7,2), complex), np.zeros((7,2), complex)
    k0[:2] = math.sqrt(p0)*np.eye(2)
    k1[2:4] = math.sqrt(p1)*np.eye(2)
    out = [k0, k1]
    for j in range(5):
        for a in range(2):
            k = np.zeros((7,2), complex)
            k[2+j,a] = math.sqrt(p2/5)
            out.append(k)
    return out

def apply(ks, rho, anc=1):
    return sum((np.kron(k,np.eye(anc)) @ rho @ np.kron(k,np.eye(anc)).conj().T for k in ks),
               np.zeros((7*anc,7*anc),complex))

def formula(u):
    return abs(u*u-.5)/2 + abs(1+6*u-8*u*u)/10 + 3*abs(1-4*u+2*u*u)/20

def build():
    baseline = json.loads((HERE/'source_baseline.json').read_text(encoding='utf-8-sig'))
    for name, digest in baseline['files_sha256'].items():
        assert sha(ROOT/name) == digest, ('source changed', name)
    bindings = {name:sha(HERE/name) for name in ('admission.md','source_baseline.json','proof.md','check.py')}
    # Exact weighted absolute-value algebra, without numerical minimization.
    pol = [(F(-1,2),F(0),F(1)),(F(1),F(6),F(-8)),(F(1),F(-4),F(2))]
    weights = [F(1,2),F(1,10),F(3,20)]
    signs = [(-1,1,1),(-1,1,-1),(1,1,-1),(1,-1,-1)]
    expected = [(F(1,2),F(0),F(-1)),(F(1,5),F(6,5),F(-8,5)),
                (F(-3,10),F(6,5),F(-3,5)),(F(-1,2),F(0),F(1))]
    for sg, ex, sample in zip(signs,expected,[F(1,10),F(1,2),F(4,5),F(19,20)]):
        got = tuple(sum(s*w*p[i] for s,w,p in zip(sg,weights,pol)) for i in range(3))
        assert got == ex
        assert all(s*peval(p,sample)>0 for s,p in zip(sg,pol))
    rt2 = (F(1414213562373095,10**15),F(1414213562373096,10**15))
    rt17 = (F(4123105625617660,10**15),F(4123105625617661,10**15))
    assert rt2[0]**2<2<rt2[1]**2 and rt17[0]**2<17<rt17[1]**2
    qlo,qhi=1/rt2[1],1/rt2[0]
    rlo,rhi=1-qhi,1-qlo
    slo,shi=(3+rt17[0])/8,(3+rt17[1])/8
    assert 0<rlo<rhi<qlo<qhi<slo<shi<1
    q=(F(0),F(1,2)); r=(F(1),F(-1,2))
    dstar=(F(-3,5),F(3,5))
    assert qeval(expected[0],r)==qeval(expected[1],r)==(F(-1),F(1))
    assert qeval(expected[1],q)==qeval(expected[2],q)==dstar
    assert F(3,5)*(rt2[0]-1)>F(497,2000)
    assert F(3,5)*(rt2[0]-1)-F(1,500)>F(493,2000)
    # Rank-three target mass and robust uniform-channel bound are exact.
    assert 3*F(1,8)==F(3,8)
    robust=F(3,8)-2*F(1,1000)
    assert robust==F(373,1000)

    target=np.diag([.25]*2+[.1]*5).astype(complex)
    sigma=np.eye(2)/2
    sx=np.array([[0,1],[1,0]],complex)
    sy=np.array([[0,-1j],[1j,0]],complex)
    sz=np.diag([1,-1]).astype(complex)
    paulis=[np.eye(2),sx,sy,sz]
    rng=np.random.default_rng(1277)
    states=[sigma,np.diag([1.,0.]),np.array([[.5,.5],[.5,.5]])]
    for _ in range(12):
        a=rng.normal(size=(2,2))+1j*rng.normal(size=(2,2))
        rho=a@a.conj().T; states.append(rho/np.trace(rho))
    us=sorted(set([j/20 for j in range(21)]+[1/math.sqrt(2)]))
    residuals=[]; min_gap=1.
    for u in us:
        ks=kraus(u)
        residuals.append(np.linalg.norm(sum(k.conj().T@k for k in ks)-np.eye(2),ord=2))
        got=apply(ks,sigma)
        residuals.append(abs(trace_distance(got,target)-formula(u)))
        for rho in states:
            out=apply(ks,rho)
            min_gap=min(min_gap,trace_distance(out,target)-formula(u))
            twirl=sum(v@rho@v.conj().T for v in paulis)/4
            residuals.append(np.linalg.norm(twirl-sigma))
            for v in paulis:
                ext=np.eye(7,dtype=complex);ext[:2,:2]=v;ext[2:4,2:4]=v
                residuals.append(np.linalg.norm(apply(ks,v@rho@v.conj().T)-ext@out@ext.conj().T))
    assert min_gap>-2e-12
    opt=apply(kraus(1/math.sqrt(2)),sigma)
    residuals.append(abs(np.trace(opt[:2,:2])-.5))
    dopt=trace_distance(opt,target)
    residuals.append(abs(dopt-3*(math.sqrt(2)-1)/5))
    residuals.append(abs(np.trace(opt[2:4,2:4]-target[2:4,2:4]).real-dopt))
    endmix=.5*apply(kraus(1),sigma)+.5*apply(kraus(0),sigma)
    residuals.append(np.linalg.norm(endmix-target))
    joint=np.zeros((14,14),complex)
    joint[:7,:7]=.5*apply(kraus(1),sigma)
    joint[7:,7:]=.5*apply(kraus(0),sigma)
    mutual=math.log(2)+entropy(endmix)-entropy(joint)
    residuals.append(abs(mutual-math.log(2)))

    bell=np.array([1,0,0,1],complex)/math.sqrt(2)
    bellrho=np.outer(bell,bell.conj())
    bellout=np.zeros(14,complex);bellout[:4]=bell
    effect=np.zeros((14,14),complex);effect[:4,:4]=np.eye(4)-bellrho
    desired=np.kron(target,np.eye(2)/2)
    reference_prob=np.trace(effect@desired).real
    residuals.append(abs(reference_prob-3/8))
    mixtures=[[(1.,0.)],[(1.,1.)],[(.5,0.),(.5,1.)],[(.2,.13),(.3,.5),(.5,.91)]]
    bell_probs=[]; bell_distances=[]
    for mixture in mixtures:
        actual=sum(w*apply(kraus(u),bellrho,2) for w,u in mixture)
        prob=np.trace(effect@actual).real
        bell_probs.append(rounded(prob));bell_distances.append(rounded(trace_distance(actual,desired)))
        residuals.append(abs(prob))
        assert trace_distance(actual,desired)>=3/8-2e-12
    max_half=0.
    for p in np.linspace(0,1,41):
        psi=np.array([math.sqrt(p),0,0,math.sqrt(1-p)],complex)
        rho=np.outer(psi,psi.conj()); rr=np.diag([p,1-p])
        diff=rho-np.kron(np.eye(2)/2,rr)
        eig=np.linalg.eigvalsh(diff)
        root=math.sqrt(1+12*p*(1-p))
        analytic=sorted([-p/2,-(1-p)/2,(1-root)/4,(1+root)/4])
        residuals.append(np.max(np.abs(eig-analytic)))
        dist=trace_distance(rho,np.kron(np.eye(2)/2,rr))
        residuals.append(abs(dist-(1+root)/4))
        actual=.5*apply(kraus(0),rho,2)+.5*apply(kraus(1),rho,2)
        wanted=np.kron(target,rr)
        got=trace_distance(actual,wanted)
        residuals.append(abs(got-dist/2));max_half=max(max_half,got)
    assert max(residuals)<2e-12 and abs(max_half-3/8)<2e-12
    return {
      'round':1277,'verdict':'PASS','author_groups':1,
      'source_bindings':bindings,'historical_sources_checked':len(baseline['files_sha256']),
      'exact':{
        'piece_coefficients_low_to_high':[[rat(x) for x in p] for p in expected],
        'sqrt2_enclosure':[rat(x) for x in rt2],
        'sqrt17_enclosure':[rat(x) for x in rt17],
        'ordered_breakpoints_certified':True,
        'optimum':'3*(sqrt(2)-1)/5','unique_u':'1/sqrt(2)',
        'optimum_strict_lower_bound':'497/2000','optimal_area_distribution':['1/2','1/2'],
        'universal_random_waiting_half_diamond':'3/8','bell_effect_rank':3,
        'robust_two_errors_each_1e_minus3':'373/1000'},
      'diagnostics':{
        'u_samples':len(us),'input_states_per_u':len(states),'schmidt_samples':41,
        'fixed_time_optimum':rounded(dopt),'all_input_distance_gap_min':rounded(min_gap),
        'endpoint_known_input_distance':rounded(trace_distance(endmix,target)),
        'retained_time_label_mutual_information':rounded(mutual),
        'bell_effect_actual_probabilities':bell_probs,'bell_distances':bell_distances,
        'endpoint_channel_half_diamond_max':rounded(max_half),
        'maximum_matrix_residual':float(f'{max(residuals):.6e}')},
      'scope':{
        'same_original_1237_parent_and_1274_reduced_channel':True,
        'all_times_and_inputs_supported_by_analytic_proof_not_sampling':True,
        'random_scheduler_physically_certified':False,'new_measurement_effect_authorized':False,
        'full_parent_gibbs_preparation_certified':False,
        'universal_replacer_is_separate_stronger_task':True,
        'all_TR_references_required_to_be_replacers':False,
        'geometry_or_Jacobson_interface_completed':False}}

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    result=build()
    if args.write:
        with RESULT.open('x',encoding='utf-8',newline='\n') as out:
            json.dump(result,out,ensure_ascii=False,indent=2,sort_keys=True)
            out.write('\n')
        status='FIRST_WRITE_PASS'
    else:
        saved=json.loads(RESULT.read_text(encoding='utf-8'))
        assert saved==result,'saved result differs; read-only mode does not write'
        status='READ_ONLY_MATCH'
    print(json.dumps({'status':status,'verdict':result['verdict'],
                      'd_star':result['diagnostics']['fixed_time_optimum'],
                      'universal_half_diamond':'3/8',
                      'maximum_residual':result['diagnostics']['maximum_matrix_residual']},ensure_ascii=False))

if __name__=='__main__':
    main()
