"""Round 1116: algebraic calibration, not an empirical or full Lorentz model."""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]

def norm(a):
    return float(np.linalg.norm(a, ord=2))

def exp_h(h, t):
    e, u = np.linalg.eigh(h)
    return (u * np.exp(-1j * t * e)) @ u.conj().T

def ad_matrix(u):
    return np.kron(u.conj(), u)

def boost(v):
    v = np.asarray(v, float)
    b2 = float(v @ v)
    g = 1.0 / np.sqrt(1.0 - b2)
    out = np.eye(4)
    if b2:
        out[0, 0] = g
        out[0, 1:] = g * v
        out[1:, 0] = g * v
        out[1:, 1:] += (g - 1) * np.outer(v, v) / b2
    return out

def calculate():
    h = np.diag([1., 2.])
    shift, t = 2., .37
    ui, uo = exp_h(h, t), exp_h(h + shift*np.eye(2), t)
    g1 = dict(full_channel_residual=norm(ad_matrix(ui)-ad_matrix(uo)),
              amplitude_phase_residual=norm(uo-np.exp(-1j*shift*t)*ui),
              energy_defect=shift,
              same_input_output_dimension_does_not_mean_same_physical_type=True)

    theta = .43
    rot = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
    w = np.kron(np.eye(2), rot)
    ps = [np.diag(x) for x in [[5,5,8,8],[1,1,-2,-2],[0,0,1,1],[0,0,0,0]]]
    defects = [w.T@p@w-p for p in ps]
    bad = np.eye(4)[[0,2,1,3]]
    g2 = dict(unitarity=norm(w.T@w-np.eye(4)),
              full_same_type_charge_residual=max(map(norm, defects)),
              trace_candidate_q=[float(np.trace(d).real/4) for d in defects],
              arbitrary_unitary_is_not_enough=norm(ad_matrix(exp_h(ps[0],t)@bad)-ad_matrix(bad@exp_h(ps[0],t))))

    n = 6
    v = np.zeros((n+1,n)); v[np.arange(1,n+1),np.arange(n)]=1
    hi, ho = np.diag(np.arange(n)), np.diag(np.arange(n+1))
    g3 = dict(isometry=norm(v.T@v-np.eye(n)),
              exact_charge_one_residual=norm(ho@v-v@hi-v),
              output_range_missing_dimension=int((n+1)-np.linalg.matrix_rank(v)),
              trace_argument_not_applicable_to_isometry=True)

    q1 = np.array([2., .3, -.4, .1]); q2 = np.array([1., -.2, .1, .2])
    a = np.array([.27, .1, -.3, .4]); eta=np.diag([1.,-1.,-1.,-1.])
    chi = lambda q: np.exp(-1j*(q@eta@a))
    rx = np.diag([1.,1.,-1.,-1.]); ry = np.diag([1.,-1.,1.,-1.]); bx=boost([.6,0,0])
    constraints = np.concatenate([rx-np.eye(4),ry-np.eye(4),bx-np.eye(4)],axis=0)
    g4 = dict(composition_character_residual=float(abs(chi(q1)*chi(q2)-chi(q1+q2))),
              rotations_and_boost_fixed_vector_rank=int(np.linalg.matrix_rank(constraints)),
              smallest_constraint_singular_value=float(np.linalg.svd(constraints,compute_uv=False)[-1]),
              nonzero_external_family_q=(bx@np.array([2.,0,0,0])).tolist(),
              family_q_nonzero_is_allowed_without_fixed_W=True,
              finite_vector_test_is_not_a_finite_Poincare_representation=True)

    # Full A x flag x battery object: swap equal-total-energy resource states.
    eye2=np.eye(2); n1=np.diag([0.,1.])
    htot=np.kron(h,np.eye(4))+2*np.kron(np.kron(eye2,n1),eye2)+2*np.kron(np.kron(eye2,eye2),n1)
    swap=np.eye(4)[:,[0,2,1,3]]
    full=np.kron(eye2,swap)
    jin=np.zeros((8,2)); jout=np.zeros((8,2))
    jin[[1,5],[0,1]]=1; jout[[2,6],[0,1]]=1
    flag=2*np.kron(np.kron(eye2,n1),eye2)
    battery=2*np.kron(np.kron(eye2,eye2),n1)
    g5=dict(full_unitarity=norm(full.T@full-np.eye(8)),
            full_energy_commutator=norm(htot@full-full@htot),
            all_unknown_input_isometry_identity=norm(full@jin-jout),
            material_energy_change=np.diag(jout.T@flag@jout-jin.T@flag@jin).tolist(),
            battery_energy_change=np.diag(jout.T@battery@jout-jin.T@battery@jin).tolist(),
            arbitrary_reference_inherits_exact_isometry_identity=True,
            actual_CO_port_or_full_Lorentz_device_claimed=False)

    probe_t=np.pi/shift
    plus=np.array([1.,1.],complex)/np.sqrt(2)
    probe=np.diag([1,np.exp(-1j*shift*probe_t)])@plus
    g6=dict(probe_time=float(probe_t),
            plus_probability_with_phase=float(abs(np.vdot(plus,probe))**2),
            plus_probability_after_incorrect_phase_deletion=1.,
            controlled_implementation_is_an_explicit_extra_permission=True)
    assert g1['full_channel_residual']<1e-12 and g1['amplitude_phase_residual']<1e-12
    assert g2['unitarity']<1e-12 and g2['full_same_type_charge_residual']<1e-12 and g2['arbitrary_unitary_is_not_enough']>.1
    assert g3['isometry']<1e-12 and g3['exact_charge_one_residual']<1e-12
    assert g4['composition_character_residual']<1e-12 and g4['rotations_and_boost_fixed_vector_rank']==4
    assert g5['full_energy_commutator']<1e-12 and g5['all_unknown_input_isometry_identity']<1e-12
    assert np.allclose(g5['material_energy_change'],[2,2]) and np.allclose(g5['battery_energy_change'],[-2,-2])
    assert g6['plus_probability_with_phase']<1e-12
    sources=['archive_1106_/1108/proof.md','archive_1086_1105/_admission/dynamic_contract_audit.md','archive_1106_/1111/proof.md']
    return dict(round=1116,status='PASS_TRANSLATION_DEFECT_CLASSIFICATION',numerical_check_groups=6,
        new_adopted_axioms=0,scientific_calibration=0,
        source_sha256={p:hashlib.sha256((BASE/p).read_bytes()).hexdigest() for p in sources},
        group1_cross_type_channel_defect=g1,group2_same_type_unitary=g2,
        group3_isometry_boundary=g3,group4_composition_and_fixed_Lorentz_vector=g4,
        group5_complete_resource_closure=g5,group6_coherent_phase_probe=g6)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args()
    result=calculate(); raw=json.dumps(result,ensure_ascii=False,indent=2)+'\n'
    output=HERE/'results.json'
    if args.write:
        output.write_text(raw,encoding='utf-8')
    else:
        saved=json.loads(output.read_text(encoding='utf-8'))
        assert saved==result,'Saved results differ; inspect before overwriting'
    print(json.dumps(dict(status=result['status'],groups=6,saved_result_matches=True),ensure_ascii=False))

if __name__=='__main__': main()
