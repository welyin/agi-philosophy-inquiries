"""Count-zero anchored direction bridge: independent finite audit.

Default stdout JSON only. No source propagation, recurrence search, or new
control model. Eight classical branches retain the one unknown Q, blank M,
and reference R. Frozen 1080 source effects are algebraic diagnostic inputs;
the second source point is not inside the certified tiny source ball.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import importlib.util
import json
import numpy as np

BASE=Path(__file__).resolve().parent
ARCHIVE=BASE.parent.parent
GROUPS=[]

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def load(p):
    spec=importlib.util.spec_from_file_location('frozen1079_process_helpers',p)
    m=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m

def check(name,condition):
    if not bool(condition):
        raise AssertionError(name)
    GROUPS.append(name)

def partial_effect(effect,rho_qr):
    # Unnormalized reference state for this effect, independent of output-Q map.
    return np.einsum('ji,irjs->rs',effect,rho_qr.reshape(2,2,2,2))

def full_eight_branch(helper,effects,rho_qr):
    initial=helper.prepare_qar(rho_qr)
    reset=helper.SWAP@initial@helper.SWAP.conj().T
    kraus=[helper.raw_measure_prepare_kraus(e) for e in effects]
    result=[]
    for r in range(2):
        state=initial if r==0 else reset
        for j in range(2):
            for y in range(2):
                z=y^j^r
                output=helper.apply_raw(state,kraus[j],z)/4
                if r==0:
                    e=(effects[j] if z==0 else helper.I-effects[j])/4
                else:
                    c=float(np.trace(effects[j]).real/2)
                    e=((c if z==0 else 1-c)/4)*helper.I
                result.append((r,j,y,z,e,output))
    return result

def run():
    GROUPS.clear()
    helpers_path=ARCHIVE/'1079/independent_check.py'
    saved_source=ARCHIVE/'1080/independent_results.json'
    source_certificate=ARCHIVE/'1077/results.json'
    source_proof=ARCHIVE/'1077/proof.md'
    old_budget=BASE.parent/'history_aware_endpoint_after1075/results.json'
    h=load(helpers_path)
    saved=json.loads(saved_source.read_text(encoding='utf-8'))
    coefficients=[np.array(s['effect_coefficients'],float) for s in saved['sources']]
    # x=p2, o=p1. Both are actual old fixed-reader results, used only to check
    # this finite instrument algebra, not to claim membership of the tiny ball.
    co,cx=coefficients
    effect_o=co[0]*h.I+np.einsum('i,ijk->jk',co[1:],h.PAULI)
    effect_x=cx[0]*h.I+np.einsum('i,ijk->jk',cx[1:],h.PAULI)
    effects=[effect_x,effect_o]
    eigs=[np.linalg.eigvalsh(e) for e in effects]
    check('frozen actual source effects are valid and distinct',
          saved['status']=='passed' and all(min(s)>0 and max(s)<1 for s in eigs)
          and np.linalg.norm(effect_x-effect_o)>1e-4)
    delta=(h.I+effect_x-effect_o)/2
    coarse=h.I/2+np.einsum('i,ijk->jk',(cx[1:]-co[1:])/4,h.PAULI)
    check('two-stage effect algebra is positive and has fixed trace',
          np.max(abs(coarse-(h.I/2+(delta-np.trace(delta)*h.I/2)/2)))<2e-15
          and np.linalg.eigvalsh(delta)[0]>0 and np.linalg.eigvalsh(delta)[-1]<1
          and abs(np.trace(coarse)-1)<2e-15)
    calibration=[h.I/2]+[(h.I+p)/2 for p in h.PAULI]
    states=[np.kron(r,np.array([[.4,.1j],[-.1j,.6]])) for r in calibration]
    states+=h.references()
    prob_res=ref_res=storage_res=tp_res=coarse_ref_res=post_q_res=0.
    smallest=1.
    complete_res=0.
    fine_probabilities=[]
    for e in effects:
        ks=h.raw_measure_prepare_kraus(e)
        complete=sum(k.conj().T@k for branch in ks for k in branch)
        complete_res=max(complete_res,float(np.max(abs(complete-h.I))))
    for rho in states:
        out=full_eight_branch(h,effects,rho)
        rho_q=h.trace_sub(rho,(2,2),[0])
        probabilities=[]
        for r,j,y,z,e,a in out:
            p=float(np.trace(a).real)
            probabilities.append(p)
            prob_res=max(prob_res,abs(p-float(np.trace(e@rho_q).real)))
            smallest=min(smallest,float(np.linalg.eigvalsh(a)[0]))
            if p>1e-14:
                outq=h.trace_sub(a,(2,2,2),[0])/p
                expected=np.zeros((2,2),complex);expected[z,z]=1
                post_q_res=max(post_q_res,float(np.max(abs(outq-expected))))
                if r==1:
                    stored=h.trace_sub(a,(2,2,2),[1,2])/p
                    storage_res=max(storage_res,float(np.max(abs(stored-rho))))
        total=sum(a for r,j,y,z,e,a in out)
        tp_res=max(tp_res,abs(float(np.trace(total).real)-1))
        ref_res=max(ref_res,float(np.max(abs(h.trace_sub(total,(2,2,2),[2])
                                              -h.trace_sub(rho,(2,2),[1])))))
        positive=sum(a for r,j,y,z,e,a in out if y==0)
        coarse_ref_res=max(coarse_ref_res,float(np.max(abs(
            h.trace_sub(positive,(2,2,2),[2])-partial_effect(coarse,rho)))))
        fine_probabilities.append(probabilities)
    check('eight retained QMR branches reproduce fine and coarse probabilities',
          max(prob_res,coarse_ref_res,complete_res)<3e-13 and smallest>-3e-13)
    check('full process preserves reference trace and stores unknown QR on reset branch',
          max(ref_res,storage_res,tp_res,post_q_res)<3e-13)
    # The anchor coincidence is neutral only for the selected coarse task.
    same=full_eight_branch(h,[effect_o,effect_o],h.references()[0])
    same_effect=sum(e for r,j,y,z,e,a in same if y==0)
    fine_trace=[float(np.trace(e).real) for r,j,y,z,e,a in same]
    check('coincident source and anchor are coarse neutral but fine task is retained',
          np.max(abs(same_effect-h.I/2))<2e-15
          and max(fine_trace)-min(fine_trace)>1e-3
          and len(same)==8)

    # Exact arithmetic reuses the frozen source certificate and old compiler.
    cert=json.loads(source_certificate.read_text(encoding='utf-8'))
    budget=json.loads(old_budget.read_text(encoding='utf-8'))
    r0=F(cert['reused_480_control_radius']);B=F(30000);L=F(2052)
    r=r0/2;epsilon=r/(16*B);delta_target=r/(4*B);a=r0/(32*B)
    source_radius=r0/4
    source_contraction=B*L*source_radius
    source_map_radius=source_contraction*source_radius+B*a
    compiler_contraction=F(1,4)+76*B*epsilon
    check('exact source ball and one-request compiler budgets cover every rotation',
          F(cert['bloch_inverse_infinity_certified_upper'])<B
          and a==delta_target/4==epsilon and 2*a<delta_target
          and source_contraction==F(1,16)
          and source_map_radius<source_radius
          and F(budget['r'])==r and F(budget['epsilon'])==epsilon
          and F(budget['delta'])==delta_target
          and compiler_contraction==F(budget['contraction_lambda'])<1)
    # Orthogonal R preserves |w|=a. Every pair in the closed ball differs by
    # at most 2a in Euclidean and hence infinity norm. This analytic inequality,
    # not finite sampling below, gives the all-R/all-shell request budget.
    rx=((1,0,0),(0,0,-1),(0,1,0))
    rz=((0,-1,0),(1,0,0),(0,0,1))
    def mv(m,v):return tuple(sum(F(m[i][j])*v[j] for j in range(3)) for i in range(3))
    ez=(F(0),F(0),F(1));ez_real=tuple(a*x for x in ez)
    zx=mv(rz,mv(rx,ez_real));xz=mv(rx,mv(rz,ez_real))
    difference=tuple(x-y for x,y in zip(zx,xz))
    antipode_effect_gap=a/2
    order_effect_gap_squared=sum(x*x for x in difference)/16
    plus_x_probability_difference=difference[0]/4
    check('exact antipode contrast and quarter-turn order witness use the new normalization',
          zx==(a,F(0),F(0)) and xz==(F(0),-a,F(0))
          and antipode_effect_gap>0 and order_effect_gap_squared==a*a/8
          and plus_x_probability_difference==a/4>0 and a<2)
    # Reflection acts on the rich source history, not on an arbitrary unknown Q.
    # Synthetic identity Pauli actions only verify the already mature covariance.
    ux=np.cos(np.pi/4)*h.I-1j*np.sin(np.pi/4)*h.PAULI[0]
    uz=np.cos(np.pi/4)*h.I-1j*np.sin(np.pi/4)*h.PAULI[2]
    covariance_res=max(float(np.max(abs(u@h.PAULI[k]@u.conj().T
        -np.einsum('i,ijk->jk',np.asarray(m,dtype=float)[:,k],h.PAULI))))
        for u,m in ((ux,rx),(uz,rz)) for k in range(3))
    check('mature Pauli covariance control has the declared rotation convention',covariance_res<2e-15)
    dependencies={
        '1079/independent_check.py':helpers_path,
        '1080/independent_results.json':saved_source,
        '1077/results.json':source_certificate,
        '1077/proof.md':source_proof,
        '_admission/history_aware_endpoint_after1075/results.json':old_budget,
    }
    return {
        'status':'passed','scientific_round_increment':0,'assertion_groups':len(GROUPS),
        'groups':GROUPS,
        'source_diagnostic':{
            'anchor_parameters':saved['sources'][0]['source_parameters'],
            'target_parameters':saved['sources'][1]['source_parameters'],
            'anchor_effect_coefficients':co.tolist(),'target_effect_coefficients':cx.tolist(),
            'raw_effect_eigenvalues':[v.tolist() for v in eigs],
            'difference_effect_eigenvalues':np.linalg.eigvalsh(delta).tolist(),
            'coarse_effect_coefficients':np.r_[.5,(cx[1:]-co[1:])/4].tolist(),
            'coarse_effect_eigenvalues':np.linalg.eigvalsh(coarse).tolist(),
            'second_source_is_small_ball_point':False,
            'old_384_dimensional_propagation_recomputed':False,
        },
        'retained_instrument':{
            'branch_order':'reset r, source j (0 target, 1 anchor), reported y',
            'actual_raw_result':'z = y XOR j XOR r','branch_weight':'1/4',
            'quantum_order':'Q,M,R; M is initially maximally mixed',
            'fine_branch_count':8,'input_states_checked':len(states),
            'max_fine_probability_residual':prob_res,
            'max_coarse_conditioned_reference_residual':coarse_ref_res,
            'max_full_reference_residual':ref_res,
            'max_reset_branch_stored_QR_residual':storage_res,
            'max_full_trace_preservation_residual':tp_res,
            'max_raw_instrument_TP_residual':complete_res,
            'max_original_Z_output_residual':post_q_res,
            'minimum_branch_eigenvalue_diagnostic':smallest,
            'fine_probabilities':fine_probabilities,
            'coincident_anchor_fine_effect_traces':fine_trace,
        },
        'exact_budget':{
            'r0':str(r0),'B':str(B),'source_parameter_radius':str(source_radius),
            'source_contraction':str(source_contraction),
            'source_map_radius_bound':str(source_map_radius),
            'request_delta':str(delta_target),'shell_radius_a':str(a),
            'all_pair_request_distance_upper':str(2*a),
            'history_recurrence_epsilon':str(epsilon),
            'compiler_contraction':str(compiler_contraction),
            'antipode_effect_operator_gap':str(antipode_effect_gap),
            'quarter_turn_order_effect_gap_squared':str(order_effect_gap_squared),
            'quarter_turn_plus_X_probability_difference':str(plus_x_probability_difference),
            'quarter_turn_endpoints':[[str(z) for z in v] for v in (zx,xz)],
            'all_rotations_budget_basis':'||R w-w||_infinity <= ||R w-w||_2 <= 2a < delta',
            'complete_source_ball_basis':'Frozen1077 inverse and Hessian bounds; contraction on C_(r0/4).',
        },
        'pauli_covariance_residual':covariance_res,
        'fixed_dependency_sha256':{k:sha(p) for k,p in dependencies.items()},
        'independent_code_sha256':sha(BASE/'independent_check.py'),
        'scope':{
            'new_bridge':'Combine same actual source, physical anchor comparison and known-history requests into a complete relative task shell.',
            'new_cognitive_principle':False,'new_formal_round':False,
            'unknown_probe_cloned':False,
            'fresh_QR_independent_of_source_and_anchor_resources_required':True,
            'same_calibrated_probe_interface_for_both_sources_required':True,
            'anchor_preparation_or_refresh_resource_is_input':True,
            'persistent_undisturbed_catalytic_anchor_proved':False,
            'full_fine_process_isospectral_or_covariant_claimed':False,
            'arbitrary_read_then_update_histories_covered':False,
            'limited_known_unitary_history_endpoint_quotient':True,
            'same_b_implies_original_Q_or_all_future_identity':False,
            'actual_spatial_dimension_derived':False,
            'tiny_shell_contrasts_resolved_by_floating_differences':False,
            'uniform_finite_resource_or_runtime_bound_derived':False,
            'finite_matrix_checks_replace_universal_CP_or_continuity_proof':False,
        },
    }

if __name__=='__main__':
    print(json.dumps(run(),ensure_ascii=False,indent=2))
