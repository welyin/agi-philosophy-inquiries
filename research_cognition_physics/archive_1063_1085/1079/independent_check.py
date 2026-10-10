"""1079 independent coarse-instrument audit. Default stdout JSON only.

Recompute the actual source numerically through frozen 1077. Retain Q-A-R
and classical branch/result outputs in the SWAP realization. The reduced raw
Q instrument is the original final Z measure-and-prepare, not a Luders model.
Exact source statements reuse the saved, outward rational 1078 boxes.
"""
from fractions import Fraction as Fraction
from pathlib import Path
import hashlib
import importlib.util
import json
import numpy as np

BASE=Path(__file__).resolve().parent
GROUPS=[]
I=np.eye(2,dtype=complex)
PAULI=np.array([[[0,1],[1,0]],[[0,-1j],[1j,0]],[[1,0],[0,-1]]],complex)

def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def checked(name,condition):
    if not bool(condition):
        raise AssertionError(name)
    GROUPS.append(name)

def norm1(a):
    return float(np.sum(np.abs(np.linalg.eigvalsh((a+a.conj().T)/2))))

def trace_sub(a,dims,keep):
    dims=list(dims); current=list(range(len(dims)))
    x=a.reshape(dims+dims)
    for original in sorted(set(current)-set(keep),reverse=True):
        axis=current.index(original)
        x=np.trace(x,axis1=axis,axis2=axis+len(current))
        current.pop(axis)
    d=int(np.prod([dims[k] for k in keep]))
    return x.reshape(d,d)

def prepare_qar(rho_qr):
    return np.einsum("irjs,ab->iarjbs",rho_qr.reshape(2,2,2,2),I/2).reshape(8,8)

def swap_qa():
    s=np.zeros((8,8),complex)
    for q in range(2):
        for a in range(2):
            for r in range(2):
                s[(a*2+q)*2+r,(q*2+a)*2+r]=1
    return s

SWAP=swap_qa()

def raw_measure_prepare_kraus(effect):
    """Original final node-Z read, reduced to the measured Q output."""
    result=[]
    for outcome,es in enumerate((effect,I-effect)):
        vals,vecs=np.linalg.eigh(es)
        if min(vals)<-1e-12:
            raise AssertionError("invalid raw effect")
        out=np.eye(2,dtype=complex)[:,outcome]
        result.append([np.sqrt(max(0.,v))*np.outer(out,vecs[:,k].conj())
                       for k,v in enumerate(vals)])
    return result

def apply_raw(rho_qar,kraus,outcome):
    result=np.zeros((8,8),complex)
    for k in kraus[outcome]:
        a=np.kron(k,np.eye(4))
        result+=a@rho_qar@a.conj().T
    return result

def blocks(rho_qr,kraus,coin=.5,raw_error=0.,reset_error=0.):
    """Four retained (coin,reported-result) QAR blocks.

Noise is a mixture of full declared channels, so its budget covers storage,
not just the working-qubit marginal. Returned blocks include branch weights.
"""
    initial=prepare_qar(rho_qr)
    reset=(1-reset_error)*(SWAP@initial@SWAP.conj().T)+reset_error*initial
    result=[]
    for branch in range(2):
        state=initial if branch==0 else reset
        weight=coin if branch==0 else 1-coin
        for reported in range(2):
            actual=reported^branch
            value=(1-raw_error)*apply_raw(state,kraus,actual)
            value+=raw_error*apply_raw(state,kraus,1-actual)
            result.append(weight*value)
    return result

def references():
    bell=np.array([1,0,0,1],complex)/np.sqrt(2)
    psi=np.array([np.sqrt(.7),0,0,np.exp(.31j)*np.sqrt(.3)],complex)
    return [np.outer(bell,bell.conj()),
            (6*np.outer(psi,psi.conj())+np.eye(4)/4)/7]

def full_process_checks(effect):
    ks=raw_measure_prepare_kraus(effect)
    c=float(np.trace(effect).real/2)
    tf=effect-c*I
    coarse=I/2+tf/2
    fine=[effect/2,(I-effect)/2,(1-c)*I/2,c*I/2]
    calibration=[I/2]+[(I+p)/2 for p in PAULI]
    worst_probability=worst_ref=worst_storage=worst_tp=0.
    smallest=1.
    retained_fine=[]
    for rho_q in calibration:
        rho=np.kron(rho_q,np.array([[.4,.1j],[-.1j,.6]]))
        out=blocks(rho,ks)
        probs=np.array([np.trace(a).real for a in out])
        expected=np.array([np.trace(f@rho_q).real for f in fine])
        worst_probability=max(worst_probability,float(np.max(np.abs(probs-expected))))
        worst_probability=max(worst_probability,abs(probs[0]+probs[2]-np.trace(coarse@rho_q).real))
        retained_fine.append(probs.tolist())
    for rho in references():
        out=blocks(rho,ks)
        worst_tp=max(worst_tp,abs(sum(np.trace(a).real for a in out)-1))
        smallest=min(smallest,min(np.linalg.eigvalsh(a)[0] for a in out))
        r_before=trace_sub(rho,(2,2),[1])
        r_after=sum(trace_sub(a,(2,2,2),[2]) for a in out)
        worst_ref=max(worst_ref,float(np.max(abs(r_after-r_before))))
        # Each reset result is independent of the old Q-R data; it retains it
        # on A-R rather than claiming it stays on the measured Q-R output.
        for a in out[2:]:
            probability=float(np.trace(a).real)
            if probability>1e-14:
                stored=trace_sub(a,(2,2,2),[1,2])/probability
                worst_storage=max(worst_storage,float(np.max(abs(stored-rho))))
    # The two-input Kraus matrices give complete CP behavior with any reference.
    complete=sum(k.conj().T@k for branch in ks for k in branch)
    checked("complete QAR instrument, retained fine probabilities, and reference trace",
            max(worst_probability,worst_ref,worst_storage,worst_tp,
                float(np.max(abs(complete-I))))<2e-13 and smallest>-2e-13)
    checked("fine coin record still retains raw trace bias",
            abs(np.trace(fine[0]).real-c)<1e-14
            and abs(np.trace(fine[2]).real-(1-c))<1e-14
            and abs(c-.5)>1e-3)
    # Direct raw measure-and-prepare is not substituted by sqrt(E) rho sqrt(E).
    vals,vecs=np.linalg.eigh(effect)
    root=(vecs*np.sqrt(np.maximum(vals,0)))@vecs.conj().T
    j_actual=sum(np.outer(k.ravel(order="F"),k.ravel(order="F").conj()) for k in ks[0])
    vr=root.ravel(order="F")
    gap=float(np.linalg.norm(j_actual-np.outer(vr,vr.conj())))
    checked("raw measured-Q output is distinguished from a Luders surrogate",gap>1e-3)
    # Error bound checks on entangled inputs and full retained output blocks.
    ei,delta,ed=.031,.017,.023
    distances=[]
    for rho in references():
        exact=blocks(rho,ks)
        noisy=blocks(rho,ks,coin=.5+delta,raw_error=ei,reset_error=ed)
        distances.append(sum(norm1(a-b) for a,b in zip(exact,noisy))/2)
    checked("finite-error full-output witnesses obey epsilon_I + delta + epsilon_D",
            max(distances)<=ei+abs(delta)+ed+1e-13)
    return {
      "coarse_effect_eigenvalues":np.linalg.eigvalsh(coarse).tolist(),
      "max_probability_residual":worst_probability,
      "max_reference_marginal_residual":worst_ref,
      "max_reset_branch_stored_QR_residual":worst_storage,
      "max_trace_preservation_residual":worst_tp,
      "smallest_branch_eigenvalue_diagnostic":smallest,
      "fine_effect_traces":[float(np.trace(f).real) for f in fine],
      "four_calibration_fine_probabilities":retained_fine,
      "raw_measure_prepare_vs_luders_choi_frobenius":gap,
      "finite_error":{"epsilon_I":ei,"coin_total_variation":delta,
        "epsilon_complete_reset_with_storage":ed,
        "conservative_half_diamond_upper":ei+abs(delta)+ed,
        "full_QAR_reference_trace_distance_witnesses":distances,
        "samples_do_not_certify_diamond_supremum":True},
    }

def main():
    old=load("frozen1077_independent_physics",BASE.parent/"1077/independent_check.py")
    source=old.run_independent()
    center=np.asarray(source["effect_center"],float)
    j=np.asarray(source["effect_Jacobian"],float)
    c,b=center[0],center[1:]
    effect=c*I+np.einsum("i,ijk->jk",b,PAULI)
    new_center=np.r_[.5,b/2]
    new_j=np.vstack((np.zeros(3),j[1:]/2))
    checked("independent original propagation and source continuity data reused",
            source["scope"]["all_graph_coherences_retained"]
            and not source["scope"]["coordinate_to_Pauli_feedback"]
            and source["graph_coherence_frobenius"]>.38)
    process=full_process_checks(effect)
    g=2*b@j[1:]
    _,_,vh=np.linalg.svd(g[None,:])
    tangent=vh[1:].T
    singular=np.linalg.svd(j[1:]@tangent/2,compute_uv=False)
    checked("centered Bloch response and angular tangent have diagnostic ranks three and two",
            np.min(np.linalg.svd(new_j[1:],compute_uv=False))>1e-4
            and np.min(singular)>1e-4 and abs(g[1]/4)>.0025
            and np.linalg.norm(g@tangent)<1e-15)
    saved=json.loads((BASE.parent/"1078/results.json").read_text(encoding="utf-8"))
    scalar_box=tuple(Fraction(a)/4 for a in saved["uniform_spectral_derivative_boxes"][1][1])
    bz_box=tuple(Fraction(a)/2 for a in saved["uniform_effect_boxes"][3])
    checked("frozen rational boxes prove nonzero spectral gradient and Bloch amplitude",
            scalar_box[1]<-Fraction(1,400) and bz_box[0]>Fraction(2,25))
    lower=Fraction(saved["reused_bloch_embedding_lower"])/2
    checked("centered map inherits positive strict source embedding budget",
            lower==Fraction(59,240000)>0)
    # Small instrument error controls statistics, not exact level-set topology.
    # An unfair coin retains raw c: c_delta=1/2-delta+2*delta*c,
    # b_delta=(1/2+delta)*b. The old nonzero minor survives multiplied by
    # 2*delta*(1/2+delta)^2, for arbitrarily small nonzero delta.
    coin_bias=Fraction(1,1000)
    bias_factor=2*coin_bias*(Fraction(1,2)+coin_bias)**2
    bias_minor=tuple(bias_factor*Fraction(x) for x in saved["uniform_minor_enclosure"])
    checked("small coin bias need not preserve the exact two-dimensional isospectral level",
            bias_factor>0 and bias_minor[1]<0)
    hashes={"independent_check":sha(Path(__file__))}
    for rel in ("1077/independent_check.py","1077/proof.md","1078/check.py",
                "1078/results.json","1078/proof.md"):
        hashes[rel]=sha(BASE.parent/rel)
    # Binding only; no author algorithms are used for new instrument checks.
    for name in ("check.py","results.json","proof.md","source_scope.md"):
        if (BASE/name).exists():
            hashes["author_"+name]=sha(BASE/name)
    result={"round":1079,"status":"passed","assertion_groups":len(GROUPS),"groups":GROUPS,
      "original_center_effect":center.tolist(),"centered_effect":new_center.tolist(),
      "centered_effect_jacobian":new_j.tolist(),
      "raw_norm_squared_gradient":g.tolist(),
      "centered_norm_squared_gradient":(g/4).tolist(),
      "tangent_source_columns_diagnostic":tangent.tolist(),
      "angular_tangent_singular_values_diagnostic":singular.tolist(),
      "strict_centered_spectral_u_derivative_box":[str(v) for v in scalar_box],
      "strict_centered_bz_box":[str(v) for v in bz_box],
      "strict_centered_bloch_embedding_lower":str(lower),
      "exactness_boundary":{"example_coin_bias":str(coin_bias),
        "raw_spectral_minor_multiplier":str(bias_factor),
        "perturbed_strict_spectral_minor_box":[str(x) for x in bias_minor],
        "statistical_error_control_implies_exact_level_set_dimension_stability":False},
      "process":process,"hashes":hashes,
      "scope":{"actual_coarse_binary_effect_realized_with_added_fixed_operations":True,
        "raw_source_propagation_changed":False,"coordinate_estimation_or_feedback":False,
        "maximally_mixed_blank_and_fair_coin_are_declared_resources":True,
        "coin_original_labels_and_storage_are_retained":True,
        "whole_fine_instrument_is_isospectral":False,
        "raw_afterstate_assumed_luders":False,
        "old_QR_stays_on_working_QR_in_reset_branch":False,
        "coarse_effect_unique_afterstate_claimed":False,
        "strict_gradient_comes_from_reused_outward_rational_boxes":True,
        "numerical_tangent_rank_used_as_proof":False,
        "actual_rotations_reverse_or_spatial_dimension_derived":False,
        "general_arbitrary_Q_source_correlations_allowed_for_effect_formula":False,
        "complete_environment_state_of_raw_384dim_process_recomputed_here":False}}
    print(json.dumps(result,ensure_ascii=False,indent=2))

if __name__=="__main__":
    main()
