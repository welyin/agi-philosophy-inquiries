"""624: original record probabilities and source-conditioned response.

Exact generator diagnostics use the original full Hamiltonian commutator:
all multiplicative potentials, CAR mass and spectator terms cancel for f(s).
No full graph propagator or new measuring apparatus is simulated.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np
import joint_geometry_work_noise as original

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_record_source_response_results.json'
PHASE=3.


def conditional_rows(nodes):
    h,s,p,kss=original.packet_data(nodes)
    velocity=original.HBAR*PHASE*kss/original.W
    rows=[]
    for sign in (1,-1):
        effect=.5+sign*np.sin(s)/4
        derivative=sign*np.cos(s)/4
        prob=float(np.sum(p*effect))
        dp=float(np.sum(p*derivative*velocity))
        mean=float(np.sum(p*effect*s)/prob)
        numerator_derivative=float(np.sum(p*(effect+s*derivative)*velocity))
        actual=(numerator_derivative-mean*dp)/prob
        post_state_only=float(np.sum(p*effect*velocity)/prob)
        correction=float(np.sum(p*(s-mean)*derivative*velocity)/prob)
        assert abs(actual-post_state_only-correction)<2e-14
        rows.append(dict(sign=sign,probability=prob,probability_derivative=dp,
                         mean=mean,conditional_mean_derivative=actual,
                         post_state_only=post_state_only,selection_correction=correction))
    unconditional=float(np.sum(p*velocity))
    weighted=sum(r['probability']*r['conditional_mean_derivative'] for r in rows)
    weight_term=sum(r['probability_derivative']*r['mean'] for r in rows)
    assert abs(weighted+weight_term-unconditional)<2e-14
    assert abs(weight_term)>1e-4
    assert all(abs(r['selection_correction']) > 100*np.spacing(1+abs(r['conditional_mean_derivative'])) for r in rows)
    fisher=sum(r['probability_derivative']**2/r['probability'] for r in rows)
    return dict(nodes=nodes,rows=rows,unconditional_mean_derivative=unconditional,
                weighted_conditional_derivative=weighted,probability_weight_term=weight_term,
                classical_record_fisher=fisher)


def conditional_check():
    rows=[conditional_rows(n) for n in (64,96,128)]
    gap=max(abs(rows[-1][key]-rows[-2][key]) for key in
            ('unconditional_mean_derivative','probability_weight_term','classical_record_fisher'))
    assert gap<2e-12
    return dict(quadrature=rows,convergence_error=float(gap),
                original_compact_Gauss_packet_with_singlet_phase=True,
                original_full_H_instantaneous_generator_not_time_propagation=True)


def commutator_check():
    h,s=.537,.419
    metric,drift=original.metric_radial(h,s)
    alpha=original.HBAR**2/(2*original.W)
    def psi(x,y):return original.amplitude(x,y)*np.exp(1j*PHASE*y)
    def effect(y):return .5+np.sin(y)/4
    def lap(f,step):
        center=f(h,s)
        dh=(f(h+step,s)-f(h-step,s))/(2*step)
        ds=(f(h,s+step)-f(h,s-step))/(2*step)
        hh=(f(h+step,s)-2*center+f(h-step,s))/step**2
        ss=(f(h,s+step)-2*center+f(h,s-step))/step**2
        hs=(f(h+step,s+step)-f(h+step,s-step)-f(h-step,s+step)+f(h-step,s-step))/(4*step**2)
        return metric[0,0]*hh+2*metric[0,1]*hs+metric[1,1]*ss+drift[0]*dh+drift[1]*ds
    expected=original.HBAR*PHASE*metric[1,1]*np.cos(s)/(4*original.W)
    rows=[]
    for step in (.0005,.00025,.000125):
        comm=-alpha*(lap(lambda x,y:effect(y)*psi(x,y),step)-effect(s)*lap(psi,step))
        actual=float(np.real(1j*comm/(original.HBAR*psi(h,s))))
        rows.append(dict(step=step,current=actual,error=abs(actual-expected)))
    assert rows[-1]['error']<2e-6 and rows[0]['error']>10*rows[-1]['error']
    return dict(expected_current=float(expected),rows=rows,
                target_metric_and_measure_drift_retained=True,
                all_other_full_H_terms_commute_with_original_effect=True)


def coarse_check():
    h,s,p,kss=original.packet_data(128)
    velocity=original.HBAR*PHASE*kss/original.W
    rows=[]
    for count in (1,2,3,4):
        data=[]
        for signs in itertools.product((1,-1),repeat=count):
            effects=[.5+sign*np.sin(s)/4 for sign in signs]
            E=np.prod(effects,axis=0)
            logarithmic=sum(sign*np.cos(s)/(4*e) for sign,e in zip(signs,effects))
            prob=float(np.sum(p*E))
            dp=float(np.sum(p*E*logarithmic*velocity))
            data.append((int(np.prod(signs)),prob,dp))
            assert prob>=4.**(-count)-1e-14
        full=sum(dp*dp/prob for _,prob,dp in data)
        coarse=0.;lost=0.
        for parity in (1,-1):
            group=[r for r in data if r[0]==parity]
            prob=sum(r[1] for r in group);dp=sum(r[2] for r in group)
            coarse+=dp*dp/prob
            lost+=sum(r[1]*(r[2]/r[1]-dp/prob)**2 for r in group)
        assert abs(full-coarse-lost)<2e-13
        assert abs(sum(r[1] for r in data)-1)<2e-14
        assert abs(sum(r[2] for r in data))<2e-14
        if count>1:assert lost>1e-3
        rows.append(dict(reads=count,full_record_fisher=full,parity_fisher=coarse,
                         lost_score_variance=lost,min_history_probability=min(r[1] for r in data)))
    return dict(rows=rows,all_reads_after_same_infinitesimal_original_H_flow=True,
                zero_wait_between_reads=True,coarse_reporting_not_a_new_direct_instrument=True,
                score_variance_identity_is_analytic=True)


def run():
    names=('joint_geometry_work_noise.py','research_note_586.md','research_note_592.md',
           'research_note_593.md','research_note_598.md','research_note_623.md')
    return dict(round=624,tests_run=3,failures=0,errors=0,
                conditional_response=conditional_check(),
                original_generator=commutator_check(),coarse_record_response=coarse_check(),
                dependency_hashes={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in names},
                scope='Same original fixed graph Hamiltonian and primitive singlet instrument. Full finite-history source expansion is analytic under common-domain and A2 conditions; numerics test original instantaneous generator on a normal Gauss packet. No full graph evolution simulation, new hardware, continuum or gravity completion.')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true');args=parser.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False,indent=2))
