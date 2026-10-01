"""584: original finite-graph quantum measure, records and background response.

The finite graph, positive spatial weights and Laplace--Beltrami ordering are
inherited inputs. A change of coordinate measure is not a new gravity measure.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_curved_quantum_source as original
from round584_drafts import half_density_entry_probe as entry

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_quantum_measure_records_results.json'
HBAR=.7
WEIGHT=.8
HCENTER,HRADIUS=.52,.10
SCENTER,SRADIUS=.40,.15


def compact_packet(nodes):
    z,w=np.polynomial.legendre.leggauss(nodes)
    x,y=np.meshgrid(z,z,indexing='ij')
    h=HCENTER+HRADIUS*x;s=SCENTER+SRADIUS*y
    chi=np.exp(-1/(1-x*x)-1/(1-y*y))
    weights=w[:,None]*w[None,:]*HRADIUS*SRADIUS*h**3
    Z=float(np.sum(weights*chi**2));prob=weights*chi**2/Z
    assert Z>0 and np.all(prob>=0)
    f=original.M-(h*h+s*s)/6
    mu=np.sqrt(original.M)*f**-3
    psiK=chi/np.sqrt(mu)
    norm_K=float(np.sum(weights*mu*psiK**2)/Z)
    q=2.5-(h*h+s*s)/(3*original.M)
    grad_chi=np.stack((-2*x/(HRADIUS*(1-x*x)**2),-2*y/(SRADIUS*(1-y*y)**2)),axis=-1)
    coords=np.stack((h,s),axis=-1)
    G=f[...,None,None]*(np.eye(2)-coords[..., :,None]*coords[...,None,:]/(6*original.M))
    grad_psi=grad_chi-coords/(2*f[...,None])
    kinetic_b=float(np.sum(prob*np.einsum('...a,...ab,...b->...',grad_chi,G,grad_chi)))
    kinetic_K=float(np.sum(prob*np.einsum('...a,...ab,...b->...',grad_psi,G,grad_psi)))
    qmean=float(np.sum(prob*q))
    acceleration=HBAR**2*f*f*s*np.cos(s)/(3*original.M**2*WEIGHT**2)
    phi=np.zeros(h.shape+(5,));phi[...,1]=h;phi[...,4]=s
    Umean=float(np.sum(prob*original.node_potential(phi)))
    return dict(nodes=nodes,normalizer_without_common_S3_area=Z,norm_curved=norm_K,
        kinetic_bare=kinetic_b,kinetic_curved=kinetic_K,q_mean=qmean,
        quadratic_form_identity_error=abs(kinetic_K-kinetic_b-qmean),
        mean_acceleration=float(np.sum(prob*acceleration)),
        probability_t_squared_coefficient=float(np.sum(prob*acceleration))/8,
        initial_probability=float(np.sum(prob*(.5+np.sin(s)/4))),potential_mean=Umean)


def compact_gauss_source_check():
    rows=[compact_packet(n) for n in (48,80,128)]
    fmin=original.M-((HCENTER+HRADIUS)**2+(SCENTER+SRADIUS)**2)/6
    lower=HBAR**2*fmin**2*(SCENTER-SRADIUS)*np.cos(SCENTER+SRADIUS)/(3*original.M**2*WEIGHT**2)
    assert lower>0 and all(r['mean_acceleration']>lower for r in rows)
    assert max(abs(r['norm_curved']-1) for r in rows)<1e-14
    assert abs(rows[-1]['mean_acceleration']-rows[-2]['mean_acceleration'])<1e-11
    assert rows[-1]['quadratic_form_identity_error']<1e-9
    phi=np.array([.26,.34,.18,.27,.4]);X=phi[:2]+1j*phi[2:4]
    rotation=original.lattice.old.group_exp(np.array([.4,-.2,.3]),2)
    Y=np.exp(.19j)**3*(rotation@X)
    transformed=np.r_[Y.real,Y.imag,phi[4]]
    gauge_error=max(abs(float(original.F(phi)-original.F(transformed))),
        abs(float(entry.qterm(phi)-entry.qterm(transformed))),
        abs(float(np.dot(phi[:4],phi[:4])-np.dot(transformed[:4],transformed[:4]))))
    assert gauge_error<1e-14
    return dict(rows=rows,support=dict(h=[HCENTER-HRADIUS,HCENTER+HRADIUS],
        s=[SCENTER-SRADIUS,SCENTER+SRADIUS],F_lower=fmin),
        analytic_acceleration_lower=float(lower),probability_coefficient_lower=float(lower/8),
        gauge_invariant_source_error=gauge_error,
        original_577_effect='E=1/2+sin(s)/4; probability difference coefficient is <A>/8',
        finite_time_result_is_analytic_not_a_propagation_simulation=True)


def geometric_response_check():
    packet=compact_packet(128);eps=.8;psi0=1.13
    tK,tb,q,V=(packet[n] for n in ('kinetic_curved','kinetic_bare','q_mean','potential_mean'))
    def energy(psi):
        weight=eps**3*psi**6
        return np.array([HBAR**2*tK/(2*weight)+weight*V,
            HBAR**2*(tb+q)/(2*weight)+weight*V,
            HBAR**2*tb/(2*weight)+weight*V])
    weight=eps**3*psi0**6;initial=energy(psi0)
    predicted_difference=HBAR**2*q/(2*weight)
    predicted_response=-3*HBAR**2*q/(eps**3*psi0**7)
    predicted_correct_response=6*weight/psi0*(-HBAR**2*tK/(2*weight**2)+V)
    assert abs(initial[0]-initial[1])<1e-9
    assert abs(initial[0]-initial[2]-predicted_difference)<1e-9
    rows=[]
    for step in (.004,.002,.001):
        plus,minus=energy(psi0+step),energy(psi0-step)
        derivative=(plus-minus)/(2*step)
        delta=float(derivative[0]-derivative[2])
        rows.append(dict(step=step,curved_response=float(derivative[0]),
            correct_flat_response=float(derivative[1]),bare_flat_response=float(derivative[2]),
            response_dictionary_error=float(abs(derivative[0]-derivative[1])),
            missing_q_response=delta,missing_q_response_error=abs(delta-predicted_response),
            correct_response_error=float(abs(derivative[0]-predicted_correct_response))))
    assert rows[-1]['missing_q_response_error']<rows[0]['missing_q_response_error']/12
    assert rows[-1]['correct_response_error']<rows[0]['correct_response_error']/12
    assert max(r['response_dictionary_error'] for r in rows)<1e-8
    return dict(epsilon=eps,spatial_conformal_factor=psi0,node_volume=weight,
        energy_curved=float(initial[0]),energy_correct_flat=float(initial[1]),energy_bare_flat=float(initial[2]),
        exact_missing_q_energy=predicted_difference,exact_missing_q_background_response=predicted_response,
        exact_correct_background_response=predicted_correct_response,rows=rows,
        fixed_corresponding_state_partial_derivative=True,
        node_kinetic_and_original_onsite_potential_diagnostic_not_complete_graph_stress=True)


def run():
    initial=entry.run()
    assert initial==json.loads(entry.TARGET.read_text('utf8'))
    evidence=dict(full_five_coordinate_operator=initial['operators'],
        local_record_double_commutator=initial['record_difference'],
        compact_Gauss_source=compact_gauss_source_check(),geometric_response=geometric_response_check())
    deps=('joint_curved_quantum_source.py','research_note_574.md','research_note_577.md',
        'research_note_578.md','research_note_583.md','research_round_583_checks.json',
        'round584_drafts/half_density_entry_probe.py','round584_drafts/half_density_entry_results.json')
    return dict(round=584,tests_run=len(evidence),failures=0,errors=0,checks=list(evidence),evidence=evidence,
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope='same fixed finite graph and original Laplace-Beltrami matter; unitary coordinate-measure dictionary with inherited closed domain, Gauss and source/effect mapping; bounded omitted-q competitor has a same-source short-time record difference and changed background-weight response; not a full Jordan-gravity measure, continuum quantum stress, dynamical geometry or dimensional selection')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run();payload=json.dumps(result,ensure_ascii=False,indent=2)+'\n'
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(payload)
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps({k:result[k] for k in ('round','tests_run','failures','errors')}))
