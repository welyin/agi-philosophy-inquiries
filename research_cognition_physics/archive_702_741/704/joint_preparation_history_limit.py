"""704: own preparation Gibbs state and dynamical-source mixed jets.
Original full finite-graph proof is analytic. Numerics retain702's declared
radial/neutral-Majorana diagnostic and625's actual CP record compression.
"""
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import numpy as np
import joint_gibbs_preserving_transfer as thermal
import joint_geometry_thermal_limit as geometry

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_preparation_history_limit_results.json'
KEYS=((0,0),(1,0),(0,1),(2,0),(1,1),(0,2))
BETA=1.2
HBAR=.7
TIMES=(.12,.19)


def expjet(h,g,c,t):
    scale=max(0,int(np.ceil(np.log2(max(1.,8*abs(t)*np.linalg.norm(h,np.inf))))))
    x=[-t*m/(2**scale) for m in (h,g,c)]
    term=[np.eye(len(h),dtype=complex),np.zeros_like(h,dtype=complex),np.zeros_like(h,dtype=complex)]
    total=[m.copy() for m in term]
    for k in range(1,25):
        term=[m/k for m in geometry.muljet(term,x)]
        total=[a+b for a,b in zip(total,term)]
    assert max(np.linalg.norm(m,np.inf) for m in term)<1e-17
    for _ in range(scale):total=geometry.muljet(total,total)
    return total


def trace_norm(a):return float(np.linalg.svd(a,compute_uv=False).sum())


def data():
    K,W,B,L=thermal.interacting_radial_diagnostic();u=.03
    K=np.exp(-6*u)*K;W=np.exp(6*u)*W
    H=K+W+B;e,v=np.linalg.eigh(H)
    conv=lambda m:v.conj().T@m@v
    K,W,B=map(conv,(K,W,B));L=[conv(x) for x in L]
    G=-6*K+6*W;C=36*(K+W);A=e-e[0]+1
    ranks=[]
    for desired in (24,48,96,160,224,256):
        rank=desired
        while rank<len(e) and abs(e[rank]-e[rank-1])<1e-9:rank+=1
        ranks.append(rank)
    return dict(K=K,W=W,B=B,E=e,G=G,C=C,L=L,A=A,ranks=sorted(set(ranks)),dimension=len(e))


def preparation(d,rank):
    h=np.diag(d['E'][:rank]);g=d['G'][:rank,:rank];c=d['C'][:rank,:rank]
    jets=expjet(h,g,c,BETA)
    rho,logs=geometry.statejets(jets)
    assert abs(np.trace(rho[0])-1)<2e-12
    assert max(abs(np.trace(rho[j])) for j in (1,2))<5e-10
    return rho,logs


def record_data(d,rank):
    L=[x[:rank,:rank] for x in d['L']]
    defects=[(x.conj().T@x)[:rank,:rank]-y.conj().T@y for x,y in zip(d['L'],L)]
    cp=min(float(np.linalg.eigvalsh((x+x.conj().T)/2)[0]) for x in defects)
    tp=float(np.linalg.norm(sum(x.conj().T@x+y for x,y in zip(L,defects))-np.eye(rank),2))
    assert cp>-5e-13 and tp<1e-12
    return L,defects,cp,tp


def histories_jets(d,rank,rho,plus=1.,minus=-.4):
    h=np.diag(d['E'][:rank]);g=d['G'][:rank,:rank];c=d['C'][:rank,:rank]
    L,defects,cp,tp=record_data(d,rank)
    unitary=[(expjet(h,plus*g,plus**2*c,1j*t/HBAR),
              expjet(h,minus*g,minus**2*c,1j*t/HBAR)) for t in TIMES]
    for p,m in unitary:
        assert np.linalg.norm(p[0].conj().T@p[0]-np.eye(rank))<2e-11
    result=[]
    for hist in itertools.product((0,1),repeat=len(TIMES)):
        x={key:rho[key[0]].copy() if key[1]==0 else np.zeros_like(rho[0]) for key in KEYS}
        for r,(up,um) in zip(hist,unitary):
            y={}
            for k,q in KEYS:
                total=np.zeros_like(rho[0])
                for left in range(q+1):
                    for middle in range(q-left+1):
                        right=q-left-middle
                        factor=math.factorial(q)/(math.factorial(left)*math.factorial(middle)*math.factorial(right))
                        total+=factor*up[left]@x[(k,middle)]@um[right].conj().T
                after=L[r]@total@L[r].conj().T
                after[0,0]+=np.trace(defects[r]@total)
                y[(k,q)]=after
            x=y
        result.append([complex(np.trace(x[key])) for key in KEYS])
    out=np.array(result)
    assert abs(out[:,0].sum()-1)<2e-11
    assert max(abs(out[:,j].sum()) for j in (1,3))<2e-9
    return out,cp,tp


def histories_value(d,rank,xi,epsilon,plus=1.,minus=-.4,reset=True):
    def h(u):return np.exp(-6*u)*d['K'][:rank,:rank]+np.exp(6*u)*d['W'][:rank,:rank]+d['B'][:rank,:rank]
    rho=thermal.exp_h(h(xi),BETA);rho/=np.trace(rho)
    L,defects,_,_=record_data(d,rank)
    def unitary(u,t):
        e,v=np.linalg.eigh(h(u))
        return (v*np.exp(-1j*t*e/HBAR))@v.conj().T
    steps=[(unitary(plus*epsilon,t),unitary(minus*epsilon,t)) for t in TIMES]
    out=[]
    for hist in itertools.product((0,1),repeat=len(TIMES)):
        x=rho.copy()
        for r,(up,um) in zip(hist,steps):
            x=up@x@um.conj().T;tail=np.trace(defects[r]@x)
            x=L[r]@x@L[r].conj().T
            if reset:x[0,0]+=tail
        out.append(np.trace(x))
    return np.array(out)


def preparation_weight_check(d):
    full,logfull=preparation(d,d['dimension']);rows=[]
    for rank in d['ranks']:
        rho,logs=preparation(d,rank);errors=[];sandwiched=[]
        for j in range(3):
            embedded=np.zeros_like(full[j]);embedded[:rank,:rank]=rho[j]
            errors.append(trace_norm(d['A'][:,None]*(embedded-full[j])))
            sandwiched.append(trace_norm(d['A'][:,None]*(embedded-full[j])*d['A'][None,:]))
        rows.append(dict(rank=rank,weighted_A_preparation_jet_errors=errors,
            sandwiched_A_preparation_jet_errors=sandwiched,
            log_partition_jet_errors=[abs(logs[j]-logfull[j]) for j in range(3)]))
    assert max(rows[-1]['weighted_A_preparation_jet_errors'])<1e-11
    assert max(rows[-1]['sandwiched_A_preparation_jet_errors'])<1e-10
    # Relative-operator Neumann identity on actual original conformal coefficients.
    radius=.01;probe=.006+.005j;rank=d['ranks'][2]
    h0=np.diag(d['E'][:rank]);ar=d['A'][:rank]
    vz=(np.exp(-6*probe)-1)*d['K'][:rank,:rank]+(np.exp(6*probe)-1)*d['W'][:rank,:rank]
    relative=vz/ar[None,:];q=float(np.linalg.norm(relative,2));assert q<1
    direct=ar[:,None]*np.linalg.inv(np.diag(ar)+vz)
    factored=np.linalg.inv(np.eye(rank)+relative)
    residual=float(np.linalg.norm(direct-factored,2));assert residual<3e-13
    return dict(rows=rows,complex_probe=[probe.real,probe.imag],probe_rank=rank,
        actual_relative_perturbation_norm=q,Neumann_identity_error=residual,
        fixed_baseline_projection_and_weight=True,
        general_small_neighborhood_proof_analytic_not_inferred_from_one_probe=True)


def mixed_history_check(d):
    fullrho,_=preparation(d,d['dimension'])
    full,_,_=histories_jets(d,d['dimension'],fullrho)
    rows=[]
    for rank in d['ranks']:
        rho,_=preparation(d,rank);actual,cp,tp=histories_jets(d,rank,rho)
        errors=np.max(abs(actual-full),axis=0)
        rows.append(dict(rank=rank,source_jet_max_branch_errors=errors.tolist(),
            preparation_dynamics_mixed_sum=[float(actual[:,4].sum().real),float(actual[:,4].sum().imag)],
            minimum_CP_defect=cp,trace_preservation_error=tp))
    assert max(rows[-1]['source_jet_max_branch_errors'])<2e-10
    mixed=full[:,4];assert np.max(abs(mixed))>1e-4
    rank=d['ranks'][2];rho,_=preparation(d,rank);jets,_,_=histories_jets(d,rank,rho)
    checks=[]
    for h in (.0008,.0004):
        pp=histories_value(d,rank,h,h);pm=histories_value(d,rank,h,-h)
        mp=histories_value(d,rank,-h,h);mm=histories_value(d,rank,-h,-h)
        fd=(pp-pm-mp+mm)/(4*h*h)
        checks.append(dict(step=h,mixed_derivative_error=float(np.max(abs(fd-jets[:,4])))))
    assert checks[-1]['mixed_derivative_error']<checks[0]['mixed_derivative_error']/3
    # Probability normalization on equal histories, including the compressed CP defect.
    same,_,_=histories_jets(d,rank,rho,plus=1.,minus=1.)
    normalization=float(np.max(abs(same.sum(axis=0)-np.array([1,0,0,0,0,0]))))
    assert normalization<3e-9
    noreset=histories_value(d,rank,0.,0.,plus=1.,minus=1.,reset=False)
    loss=float(1-noreset.sum().real);assert loss>1e-8
    return dict(jet_order=[list(k) for k in KEYS],rows=rows,
        full_mixed_coefficients=[[float(z.real),float(z.imag)] for z in mixed],
        full_mixed_total=[float(mixed.sum().real),float(mixed.sum().imag)],
        independent_four_point_mixed_checks=checks,
        same_history_joint_jet_normalization_error=normalization,
        omitted_original_CP_defect_probability_loss=loss,
        frozen_initial_state_would_drop_entire_mixed_derivative=True,
        source_jets_are_derivatives_not_Taylor_coefficients=True,
        exact_jet_recursion_primary_finite_difference_only_cross_check=True)


def run():
    d=data()
    deps=('research_note_603.md','research_note_623.md','research_note_624.md','research_note_625.md',
        'research_note_702.md','research_note_703.md','joint_source_preserving_compression.py',
        'joint_gibbs_preserving_transfer.py','joint_geometry_thermal_limit.py',
        'round704_drafts/own_gibbs_compression_entry.md')
    return dict(date='2026-10-02',round=704,tests_run=2,failures=0,errors=0,
        calibration=dict(dimension=d['dimension'],beta=BETA,hbar=HBAR,times=list(TIMES),
            inherited702_radial_neutral_sector=True,full_graph_CAR_not_numerically_diagonalized=True),
        preparation_weight=preparation_weight_check(d),mixed_history=mixed_history_check(d),
        dependency_hashes={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in deps},
        scope=dict(original_fixed_graph_full_Gauss_theorem=True,
            one625_Galerkin_family_with_own_Gibbs_and_actual_CP_records=True,
            joint_preparation_and_dynamical_jets_through_total_order_two=True,
            not_a_parameter_derivative_theorem_for703_log_transfer=True,
            no_uniform_spatial_limit_or_free_physical_reset=True,
            quantum_gravity_chiral_unification_and_observational_match_still_open=True))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(round=704,checks=2,mixed=result['mixed_history']['full_mixed_total'],
        independent=result['mixed_history']['independent_four_point_mixed_checks'],all_checks_passed=True)))
