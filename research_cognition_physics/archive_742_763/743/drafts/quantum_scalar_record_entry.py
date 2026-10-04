"""743 entry: original scalar fluctuations generate a nonGaussian jet.

The short-time coefficient uses the full original mass matrices and reference.
Finite conditional branches independently check that coefficient only. They
are a covariance cubature, not a physical stochastic controller or a solved
bosonic/gravitational process.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/'round742_drafts'))
import gaussian_record_entry as old
import joint_local_source_normalization as local
import joint_retarded_reference_response as response
TARGET=HERE/'quantum_scalar_record_entry_results.json'

def derivative_covariance(dP):
    n=len(dP)//2;I=np.eye(n);T=np.block([[I,I],[-1j*I,1j*I]])
    value=-1j*T@dP@T.conj().T
    assert np.max(abs(value.imag))<2e-12
    return value.real

def evolve(P,B,t):
    ev,V=np.linalg.eigh(B);U=(V*np.exp(-1j*t*ev))@V.conj().T
    return U@P@U.conj().T

def run():
    P,_,_=response.flow(0.);G=old.covariance(P);n=len(P)//2
    f=np.zeros(n);f[30]=f[62]=1/np.sqrt(2)
    u=np.r_[f,np.zeros(n)];v=np.r_[np.zeros(n),f]
    outside=np.eye(2*n)-np.outer(u,u)-np.outer(v,v)
    a=float(u@G@v);kappa=np.sqrt(1-a*a)
    W=np.array([u,v,u@G@outside/kappa,v@G@outside/kappa])
    B=response.matrices(0.)[0]
    operators=np.array([response.prior.assemble([A,A]) for A in local.basis()])
    # All five coordinate variances positive: compatible with a normal packet.
    variances=np.array([1e-4]*4+[.04])
    jets=np.array([W@derivative_covariance(-1j*(A@P-P@A))@W.T for A in operators])
    coefficients=np.array([-old.pf4(D) for D in jets])
    coefficient=float(variances@coefficients)
    assert abs(coefficient)>1e-8
    Q=np.outer(f,f);N=np.block([[Q,np.zeros_like(Q)],[np.zeros_like(Q),-Q]])
    commutators=[float(np.linalg.norm(A@N-N@A)) for A in operators]
    overlaps=[float(np.trace(A.conj().T@N).real) for A in operators]
    assert max(abs(np.array(overlaps)))<1e-13 and commutators[4]>1e-3
    shifts=[s*np.sqrt(5*variances[j])*operators[j] for j in range(5) for s in (-1,1)]
    def defect(t):
        gammas=[W@old.covariance(evolve(P,B+shift,t))@W.T for shift in shifts]
        mean=sum(gammas)/10
        return float(-sum(old.pf4(g) for g in gammas)/10+old.pf4(mean))
    rows=[]
    for t in (.08,.04,.02):
        positive=defect(t);negative=defect(-t)
        estimate=(positive+negative)/(2*t*t)
        rows.append(dict(time=t,positive_time_defect=positive,centered_quadratic_coefficient=estimate,
                         coefficient_error=abs(estimate-coefficient)))
    ratios=[rows[i]['coefficient_error']/rows[i+1]['coefficient_error'] for i in range(2)]
    assert all(3.6<r<4.4 for r in ratios),ratios
    names=('research_note_598.md','research_note_623.md','research_note_735.md','research_note_742.md',
           'joint_local_source_normalization.py','joint_retarded_reference_response.py',
           'round742_drafts/gaussian_record_entry.py')
    return dict(entry_round=743,latest_completed_round=742,formal_test_count_unchanged=3422,
        original_Nambu_dimension=len(P),original_record_p=(1+a)/2,
        mass_coordinate_variances=variances.tolist(),per_unit_variance_Wick_coefficients=coefficients.tolist(),
        original_scalar_fluctuation_Wick_coefficient=coefficient,
        original_mass_number_commutator_norms=commutators,
        original_mass_number_Hilbert_Schmidt_overlaps=overlaps,
        number_operator_relative_distance_from_mass_span=1.,rows=rows,coefficient_error_ratios=ratios,
        conditional_covariance_cubature_only=True,
        no_claim_of_original_instrument_realization_or_full_Gauss_state=True,
        no_claim_of_boson_quantum_continuum_or_joint_gravity_solution=True,
        dependencies={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names},all_checks_passed=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');a=p.parse_args();r=run()
    if a.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
