"""1045: actual 1002 post-write state and source continuation.

The source/lapse contract is explicitly adopted; this is not a
dynamical-gravity or natural-coupling-selection result. Default: read-only
numerical replay. --write is exclusively for the initial result file.
"""
from pathlib import Path
import argparse
import json
from fractions import Fraction
import hashlib
import numpy as np

HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parents[1]
OLD = RESEARCH / 'archive_990_1008/1002/thermal_record_reuse_results.json'
CERTIFICATE = RESEARCH / 'archive_956_989/980/finite_thermal_records_results.json'
TARGET = HERE / 'postwrite_source_bridge_results.json'
HISTORY = ('archive_956_989/research_note_956.md',
           'archive_956_989/research_note_957.md',
           'archive_956_989/research_note_976.md',
           'archive_956_989/research_note_977.md',
           'archive_956_989/research_note_980.md',
           'archive_956_989/980/finite_thermal_records_results.json',
           'archive_990_1008/research_note_1002.md',
           'archive_990_1008/1002/thermal_record_reuse.py',
           'archive_990_1008/1002/thermal_record_reuse_results.json',
           'archive_990_1008/1002/organization_lifecycle_adoption_v1.md',
           'archive_1009_1043/research_note_1043.md')


def trace_distance(a, b):
    return float(np.linalg.eigvalsh((a-b + (a-b).conj().T)/2).__abs__().sum()/2)


def propagator(H, t):
    e, v = np.linalg.eigh(H)
    return (v*np.exp(-1j*t*e)) @ v.conj().T


def marginals(rho):
    a = rho.reshape(2, 4, 5, 2, 4, 5)
    return [np.einsum('ijkljk->il', a),
            np.einsum('ijkimk->jm', a),
            np.einsum('ijkijn->kn', a)]


def joint_output(H, effect, rho, lambdas, wait, refdim=1):
    """Population outcome plus complete path/reference output, path starts +."""
    d = len(H)
    us = [np.kron(propagator(H, lam*wait), np.eye(refdim)) for lam in lambdas]
    outputs = []
    # The population effect is diagonal in every basis used here.
    for weights in [1-np.diag(effect).real, np.diag(effect).real]:
        out = np.zeros((2*refdim, 2*refdim), complex)
        for l in range(2):
            for m in range(2):
                a = (us[l] @ rho @ us[m].conj().T).reshape(d, refdim, d, refdim)
                block = np.einsum('i,iris->rs', weights, a)/2
                out[l*refdim:(l+1)*refdim, m*refdim:(m+1)*refdim] = block
        outputs.append(out)
    return outputs


def calculate():
    old = json.loads(OLD.read_text(encoding='utf-8-sig'))
    pars = old['parameters']
    e = np.array(pars['energies']); p = np.array(pars['thermal_populations'])
    aux = np.array(pars['auxiliary_energies']); g = pars['g']
    hx = np.array(pars['source_energies']); rev = np.arange(3, -1, -1)
    def ix(x, i, b): return (x*4+i)*5+b
    energies = (hx[:,None,None]+e[None,:,None]+aux[None,None,:]).ravel()
    H0 = np.diag(energies); V = np.zeros((40, 40))
    for i in range(4):
        a, b = ix(1,i,0), ix(1,rev[i],i+1)
        V[a,b] = V[b,a] = 1
    H = H0 + g*V
    E = np.diag(np.tile(np.repeat([1.,1.,0.,0.],5),2))
    rest = 20.
    M = rest*np.eye(40)+H
    # Exact physical invariant sector: four x=0 points and four x=1 pairs.
    columns = [ix(0,i,0) for i in range(4)]
    for i in range(4): columns.extend([ix(1,i,0), ix(1,rev[i],i+1)])
    S = np.eye(40)[:,columns]
    h = S.T@H@S; m = S.T@M@S; effect = S.T@E@S
    blocks = [[i] for i in range(4)] + [[4+2*i,5+2*i] for i in range(4)]
    def pinch(rho, refdim=1):
        answer = np.zeros_like(rho)
        for indices in blocks:
            index = [q*refdim+r for q in indices for r in range(refdim)]
            answer[np.ix_(index,index)] = rho[np.ix_(index,index)]
        return answer
    residuals = {
        'H0_V_commutator': float(np.linalg.norm(H0@V-V@H0,2)),
        'source_intertwining': float(np.linalg.norm(M@S-S@m,2)),
        'effect_intertwining': float(np.linalg.norm(E@S-S@effect,2)),
    }
    eig = np.linalg.eigvalsh(h)
    assert np.min(np.diff(eig)) > 1e-3
    # Rational disjointness certificate inherited from 980's root isolation.
    cert = json.loads(CERTIFICATE.read_text(encoding='utf-8-sig'))
    cert = cert['eigensystem_certificate']['coarse_analytic_certificate']
    intervals = [[Fraction(v) for v in pair] for pair in cert['root_brackets']]
    intervals.insert(2, [Fraction(1), Fraction(1)])
    spectral_intervals = [[a+2,b+2] for a,b in intervals]
    for a,b in intervals:
        for sign in [-1,1]:
            shift=Fraction(11,5)+sign*Fraction(1,100)
            spectral_intervals.append([a+shift,b+shift])
    spectral_intervals.sort()
    spectral_gap = min(spectral_intervals[k+1][0]-spectral_intervals[k][1]
                       for k in range(11))
    assert spectral_gap>0
    thermal_q_upper=sum(Fraction(str(pair[1])) for pair in cert['thermal_intervals'][2:])
    assert thermal_q_upper<Fraction(188903,1000000)
    rho0 = np.zeros((40,40), complex)
    for i in range(4): rho0[ix(1,i,0),ix(1,i,0)] = p[i]
    # Stay inside the *original* 1002 certified population-writing window.
    times = [np.pi/(2*g)-10, np.pi/(2*g)+10]
    histories = [propagator(H,t)@rho0@propagator(H,t).conj().T for t in times]
    marginal_gaps = [trace_distance(a,b) for a,b in zip(marginals(histories[0]),marginals(histories[1]))]
    en, basis = np.linalg.eigh(H)
    source_weights = [np.diag(basis.conj().T@rho@basis).real for rho in histories]
    source_gap = float(np.max(np.abs(source_weights[0]-source_weights[1])))
    assert max(marginal_gaps+[source_gap]) < 1e-11
    lambdas = [.99,.97]
    wait = 25*np.pi/g
    result = [joint_output(M,E,rho,lambdas,wait) for rho in histories]
    neutral = [joint_output(M,E,rho,[.98,.98],wait) for rho in histories]
    cq_gap = sum(trace_distance(a,b) for a,b in zip(*result))
    baseline_gap = sum(trace_distance(a,b) for a,b in zip(*neutral))
    contrast = 2*float(p[0]+p[1])-1
    witness_gap=contrast*np.sin(.2)
    assert abs(cq_gap-witness_gap) < 2e-9
    assert baseline_gap < 2e-9
    record_gap = max(abs(np.trace(result[0][b])-np.trace(result[1][b])) for b in range(2))
    probe_gap = trace_distance(sum(result[0]),sum(result[1]))
    assert max(record_gap,probe_gap) < 2e-9
    # Check the complete extended output on genuinely unknown X/reference input.
    rng = np.random.default_rng(10441002)
    expanded_errors = []
    for t, ls, s in [(0.0,[.99,.97],1.), (73.,[.98,.96],17.), (199.,[.99,.97],wait)]:
        z = rng.normal(size=(4,4))+1j*rng.normal(size=(4,4))
        xr = z@z.conj().T; xr /= np.trace(xr)
        rho = np.zeros((80,80),complex)
        for i in range(4):
            embedding = np.zeros((40,2))
            embedding[ix(0,i,0),0] = embedding[ix(1,i,0),1] = 1
            a = np.kron(embedding,np.eye(2))
            rho += p[i]*(a@xr@a.conj().T)
        u = np.kron(propagator(H,t),np.eye(2)); rho = u@rho@u.conj().T
        sr = np.kron(S,np.eye(2)); compressed = pinch(sr.conj().T@rho@sr,2)
        full = joint_output(M,E,rho,ls,s,2)
        reduced = joint_output(m,effect,compressed,ls,s,2)
        expanded_errors.append(sum(trace_distance(a,b) for a,b in zip(full,reduced)))
    assert max(expanded_errors)<2e-9
    # The same dictionary is sufficient on the full invariant-sector domain,
    # not merely on the special thermal writing histories.
    full_sector_errors = []
    for _ in range(3):
        z = rng.normal(size=(24,24))+1j*rng.normal(size=(24,24))
        rho = z@z.conj().T; rho /= np.trace(rho)
        ls=rng.uniform(.96,.995,2); s=float(rng.uniform(1,200))
        full=joint_output(m,effect,rho,ls,s,2)
        reduced=joint_output(m,effect,pinch(rho,2),ls,s,2)
        full_sector_errors.append(sum(trace_distance(a,b) for a,b in zip(full,reduced)))
    assert max(full_sector_errors)<2e-12
    # Independent finite rank diagnostic of the 20-dimensional visible space.
    # Analytic exponential independence, not this sample rank, proves equality.
    rows = []
    for _ in range(72):
        # Actual bounded positive-lapse menu, not an unrelated u/v square.
        s=float(rng.uniform(200,8000))
        u,v = s*rng.uniform(.96,.995,2)
        A = propagator(h,-u)@effect@propagator(h,v)
        rows.append(A.ravel())
        rows.append(propagator(h,u-v).ravel())
    sv = np.linalg.svd(np.array(rows),compute_uv=False)
    rank = int(np.sum(sv>1e-9))
    assert rank==20
    delta_rational = Fraction(13,14000)+Fraction(1,10000)
    delta = float(delta_rational)
    # q < 188903/1e6 is inherited from the earlier certified thermal intervals.
    # sin(1/5) > 1/5-(1/5)^3/6.  No numerical sine is needed for this bound.
    sine_lower = Fraction(1,5)-Fraction(1,5)**3/6
    ideal_gap_lower = (1-2*Fraction(188903,1000000))*sine_lower
    reset_gap_lower = ideal_gap_lower-2*delta_rational
    assert reset_gap_lower>Fraction(12,100)
    return {
        'round':1045,
        'date':'2026-10-08',
        'status':'author_complete_pending_independent_review',
        'new_science_groups':1,
        'new_cognitive_axioms':0,
        'stage_goal_completed_here':False,
        'all_scientific_checks_passed':True,
        'parameters':{'energies':e.tolist(),'populations':p.tolist(),'g':g,'rest':rest,
                      'postwrite_times':times,'lapses':lambdas,'feedback_wait':wait},
        'positive_full_H_min':float(np.linalg.eigvalsh(H)[0]),
        'positive_feedback_min':float(min(lambdas)*np.linalg.eigvalsh(M)[0]),
        'invariant_dimension':12,
        'visible_complex_dimension':rank,
        'visible_algebra':'C^4 plus four full M2 blocks',
        'restricted_spectrum_min_separation':float(np.min(np.diff(eig))),
        'restricted_spectrum_rational_gap_lower':str(spectral_gap),
        'inherited_identity_residuals':residuals,
        'actual_history_marginal_half_trace_gaps':marginal_gaps,
        'complete_source_spectral_probability_gap':source_gap,
        'current_population_probabilities':[float(np.trace(E@r).real) for r in histories],
        'feedback_population_path_distributions':[[np.diag(q).real.tolist() for q in pair] for pair in result],
        'joint_output_half_trace_gap':cq_gap,
        'thermal_contrast_C':contrast,
        'analytic_joint_gap_C_sin_point_two':witness_gap,
        'original_write_window_equal_prior_errors':[
            float((1-p[0]-p[1])+(contrast/2)*np.cos(g*t)**2) for t in times],
        'same_average_lapse_no_differential_gap':baseline_gap,
        'after_feedback_record_marginal_gap':float(record_gap.real),
        'after_feedback_full_path_gap':probe_gap,
        'unknown_input_reference_cp_dictionary_errors':expanded_errors,
        'arbitrary_sector_reference_cp_dictionary_errors':full_sector_errors,
        'rank_20th_singular_value':float(sv[19]),
        'rank_21st_singular_value':float(sv[20]),
        'old_reset_error_per_history':delta,
        'certified_ideal_joint_gap_rational_lower':str(ideal_gap_lower),
        'certified_thermal_and_reset_gap_rational_lower':str(reset_gap_lower),
        'certified_thermal_and_reset_gap_lower':float(reset_gap_lower),
        'ideal_common_summary_minimax_error_lower':float(ideal_gap_lower/2),
        'reset_scope':'Per-history reset error gives output separation and approximately equal summaries; it does not imply exactly identical summaries or a reset-adjusted common-summary minimax bound.',
        'historical_sha256':{p:hashlib.sha256((RESEARCH/p).read_bytes()).hexdigest() for p in HISTORY},
        'limits': ['static prescribed lapse, not dynamical Einstein feedback',
                   '1002 engineered writer and auxiliary preparation remain adopted',
                   'different actual writing histories; explicit retained clock/history can distinguish them',
                   'fixed final population effect and all path POVMs; no claim for arbitrary later material operations',
                   'full-sector necessity and original reachable-history necessity are different statements',
                   'no natural-coupling selection or whole-stage completion claim'],
    }


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    result=calculate()
    encoded=json.dumps(result,ensure_ascii=False,indent=2)+'\n'
    if args.write:
        with TARGET.open('x',encoding='utf-8') as f: f.write(encoded)
    else:
        saved=json.loads(TARGET.read_text(encoding='utf-8'))
        assert result==saved, 'Read-only replay differs from saved result.'
    print(json.dumps({'round':1045,'all_scientific_checks_passed':True,
                      'mode':'exclusive_write' if args.write else 'read_only_compare',
                      'visible_dimension':result['visible_complex_dimension'],
                      'joint_gap':result['joint_output_half_trace_gap'],
                      'certified_reset_gap_lower':result['certified_thermal_and_reset_gap_lower']},
                     ensure_ascii=False))
