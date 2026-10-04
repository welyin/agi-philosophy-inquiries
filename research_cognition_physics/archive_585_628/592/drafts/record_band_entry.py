"""592 entry: original-read commutator and common energy-moment budgets.

Checks differential identities in the original H5 model. No finite surrogate
Hamiltonian, original eigenvalues, or spectral-cut simulation is used.
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
import joint_geometry_work_noise as prior
import joint_full_spatial_metric as model
TARGET=HERE/'record_band_entry_results.json'
HBAR,W,M=prior.HBAR,prior.W,model.original.M


def instrument(s):
    e=np.stack((.5+np.sin(s)/4,.5-np.sin(s)/4),axis=-1)
    de=np.stack((np.cos(s)/4,-np.cos(s)/4),axis=-1)
    dde=np.stack((-np.sin(s)/4,np.sin(s)/4),axis=-1)
    L=np.sqrt(e)
    return L,de/(2*L),dde/(2*L)-de**2/(4*L**3)


def constants():
    A1,A2=1/16,9/32
    b0=np.sqrt(6*M)/3
    Bdelta=2*M*M*A2+2*b0*b0*A1
    c1=4*HBAR**2*M*A1/W
    c0=HBAR**4/(2*W**2)*Bdelta
    return dict(A1=A1,A2=A2,laplacian_s_bound=b0,laplacian_L_square_bound=Bdelta,
                commutator_energy_coefficient=c1,commutator_norm_coefficient=c0,
                single_read_energy_gain_bound=HBAR**2*M/(32*W))


def commutator_check():
    c=constants();s=np.linspace(-np.sqrt(6*M),np.sqrt(6*M),1201)
    L,dL,ddL=instrument(s)
    assert np.max(np.sum(dL*dL,axis=-1))<=c['A1']*(1+1e-14)
    assert np.max(np.sum(ddL*ddL,axis=-1))<c['A2']
    q=np.array([.14,.53,-.17,.1,.42]);g=model.original.inverse(q)
    drift=-model.original.F(q)*q/(3*M)
    neighbor=np.array([.23,.41,-.09,.14,.31])
    phase=np.array([.03,.13,-.02,0,.09])
    def wave(x):return np.exp(-.21*np.dot(x,x)+1j*np.dot(phase,x))
    def potential(x):return W*float(model.original.node_potential(x))+float(model.original.distance_squared(x,neighbor))/2
    L,dL,ddL=instrument(q[-1]);grad=(-.42*q+1j*phase)*wave(q)
    target=-HBAR**2/(2*W)*(2*dL*np.dot(g[-1],grad)+(g[-1,-1]*ddL+drift[-1]*dL)*wave(q))
    def H(fn,step):
        e=np.eye(5)*step;f=fn(q);lap=0j
        for a in range(5):
            lap+=g[a,a]*(fn(q+e[a])-2*f+fn(q-e[a]))/step**2
            lap+=drift[a]*(fn(q+e[a])-fn(q-e[a]))/(2*step)
            for b in range(a):
                cross=(fn(q+e[a]+e[b])-fn(q+e[a]-e[b])-fn(q-e[a]+e[b])+fn(q-e[a]-e[b]))/(4*step**2)
                lap+=2*g[a,b]*cross
        return -HBAR**2/(2*W)*lap+potential(q)*f
    errors=[]
    for step in (.02,.01,.005):
        direct=np.array([H(lambda x,r=r:instrument(x[-1])[0][r]*wave(x),step)-L[r]*H(wave,step) for r in range(2)])
        errors.append(float(np.max(abs(direct-target))))
    assert errors[-1]<errors[0]/12
    grad_energy=HBAR**2/(2*W)*float(np.vdot(grad,g@grad).real)
    lhs=float(np.sum(abs(target)**2))
    rhs=c['commutator_energy_coefficient']*grad_energy+c['commutator_norm_coefficient']*abs(wave(q))**2
    assert lhs<rhs
    return dict(constants=c,original_five_coordinate_commutator_errors=errors,
                pointwise_commutator_square=lhs,pointwise_certified_upper=rhs,
                original_potential_included_and_cancels=True,full_other_graph_terms_commute_analytically=True)


def source_tail_check():
    h,s,p,kss=prior.packet_data(128)
    a=HBAR**2/(2*W)*float(np.sum(p*kss));alpha=.4
    rows=[]
    for n in (10,30,100,300):
        weight=alpha/n**2
        overlap=np.sum(p*np.exp(1j*n*s))
        trace_distance_upper=weight*np.sqrt(max(0,1-abs(overlap)**2))
        rows.append(dict(n=n,phase_mixture_weight=weight,trace_distance_upper=float(trace_distance_upper),
                         unchanged_energy_increment=alpha*a,unchanged_conformal_source_increment=-6*alpha*a))
    assert rows[-1]['trace_distance_upper']<rows[0]['trace_distance_upper']/100
    return dict(original_compact_Gauss_phase_family=True,rows=rows,
                trace_bound_uses_convexity_and_pure_state_overlap=True,
                arbitrary_future_fixed_channel_cannot_increase_trace_distance=True,
                finite_mean_energy_does_not_control_source_continuity=True)


def budget_check():
    c=constants();energy=1000.;moment=2e6;steps=4;rows=[]
    # These are assumed input budgets, not ground-energy estimates.
    for j in range(steps+1):
        rows.append(dict(reads=j,mean_energy_bound=energy,second_moment_bound=moment))
        moment=2*moment+2*(c['commutator_energy_coefficient']*energy+c['commutator_norm_coefficient'])
        energy+=c['single_read_energy_gain_bound']
    worstE,worstB=rows[-1]['mean_energy_bound'],rows[-1]['second_moment_bound']
    approximation=[]
    # Initial compression plus four record steps. Cut values are purely
    # analysis parameters and are not an actual model energy spectrum.
    for R in (1e6,1e7,1e8):
        distance=sum(np.sqrt(r['second_moment_bound'])/R+r['second_moment_bound']/(2*R**2) for r in rows)
        T=R**(2/3)  # dimensionless energy units used in this diagnostic
        source_bound=12*T*distance+12*(worstB/T+2*np.sqrt(worstE*worstB/T))
        approximation.append(dict(spectral_regulator=R,source_regulator=T,
                                  full_history_trace_distance_bound=distance,
                                  mean_conformal_source_difference_bound=source_bound))
    assert approximation[-1]['full_history_trace_distance_bound']<approximation[0]['full_history_trace_distance_bound']/90
    assert approximation[-1]['mean_conformal_source_difference_bound']<approximation[0]['mean_conformal_source_difference_bound']
    return dict(input_budgets_are_assumptions=True,no_nonempty_budget_claim=True,
                finite_history_moment_rows=rows,conservative_approximation_bounds=approximation,
                source_bounds_not_optimized=True,actual_spectral_truncation_not_computed=True)


def run():
    evidence={f.__name__:f() for f in (commutator_check,source_tail_check,budget_check)}
    names=('research_note_555.md','research_note_590.md','research_note_591.md',
           'joint_geometry_work_noise.py','joint_full_spatial_metric.py')
    return dict(status='592 entry, not yet a completed numbered round',checks_passed=True,evidence=evidence,
                dependency_hashes={n:hashlib.sha256((HERE.parent/n).read_bytes()).hexdigest() for n in names},
                scope='Original static finite graph and primitive binary read. Analytic finite-band obstruction and moment propagation; numerical local diagnostics only. No Einstein dynamics, autonomous reset, moving geometry, or fixed-hbar continuum.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args();r=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==r
    print(json.dumps(r,ensure_ascii=False,indent=2))
