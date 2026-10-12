"""Non-numbered qualification checks after 1146.

Three groups. These are exact finite algebra / rational diagnostics, not a
finite-mode approximation theorem for the Dirac or mediator dynamics.
Default: recompute and compare the saved result without writes.
--write-results: explicitly save after1146_joint_results.json.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ARCHIVE = HERE.parent
RESULT = HERE / "after1146_joint_results.json"
CAR_PATH = ARCHIVE / "1128" / "check.py"
_spec = importlib.util.spec_from_file_location("round1128_reused_car", CAR_PATH)
car1128 = importlib.util.module_from_spec(_spec)
assert _spec.loader is not None
_spec.loader.exec_module(car1128)

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Z = np.diag([1, -1]).astype(complex)
P1 = np.diag([0, 1]).astype(complex)


def opnorm(a):
    return float(np.linalg.norm(a, 2))


def fraction(x):
    return {"exact": str(x), "numerator": x.numerator,
            "denominator": x.denominator, "float": float(x)}


def qs_add(a, b):
    return a[0]+b[0], a[1]+b[1]


def qs_mul(a, b):
    # Exact arithmetic in Q[sqrt(3)].
    return a[0]*b[0]+3*a[1]*b[1], a[0]*b[1]+a[1]*b[0]


def group_rational():
    poly = [F(0), F(0), F(30), F(-60), F(30)]
    area = sum((c/F(n+1) for n, c in enumerate(poly)), F(0))
    weighted = sum((c*(F(1,n+1)-F(1,n+2)) for n,c in enumerate(poly)), F(0))
    assert area == 1 and weighted == F(1,2)
    # Critical points of f' are u=(3 +- sqrt(3))/6. Endpoints give f'=0.
    u = F(1,2), -F(1,6)
    derivative = F(0), F(0)
    power = F(1), F(0)
    dpoly = [F(0), F(60), F(-180), F(120)]
    for c in dpoly:
        derivative = qs_add(derivative, (c*power[0], c*power[1]))
        power = qs_mul(power, u)
    assert derivative == (F(0), F(10,3))
    squared = qs_mul(derivative, derivative)
    assert squared == (F(100,3), F(0))
    assert squared[0] < 36

    mass, delta, epsilon = 10**10, F(1,1000), F(1,10**16)
    p0 = F(3*mass,4)
    projection_tail = (20+delta)/(2*mass)
    kernel = 2*epsilon*(p0+20)
    vacuum_pair = F(3,5)-kernel-2*projection_tail
    # sqrt(2) < 3/2 bounds b+2sqrt(2)b by 4b.
    assert F(3,2)**2 > 2
    prepared_pair = F(3,5)-kernel-4*projection_tail
    certified = F(5999984,10**7)
    assert vacuum_pair > certified and prepared_pair > certified
    pair_energy_upper = F(25000000050)
    actual_energy_upper = F(5*mass,2)+40+delta
    assert actual_energy_upper < pair_energy_upper
    response = 2*certified**2/pair_energy_upper
    assert response > F(2879,10**14)
    return {
        "pass": True,
        "scope": "exact rational identities and rigorous finite material thresholds",
        "slow_switch_area": fraction(area),
        "slow_switch_weighted_area": fraction(weighted),
        "slow_switch_Lipschitz_squared": fraction(squared[0]),
        "normal_projection_amplitude_upper": fraction(projection_tail),
        "spatial_kernel_error_upper": fraction(kernel),
        "vacuum_pair_lower": fraction(vacuum_pair),
        "prepared_code_pair_lower": fraction(prepared_pair),
        "published_common_pair_lower": fraction(certified),
        "pair_energy_upper": fraction(pair_energy_upper),
        "zero_frequency_response_lower": fraction(response),
        "no_UV_limit_taken": True,
    }


def group_car():
    # Closed spin sector of a genuine Dirac alpha_z/beta representation.
    # Four Cauchy modes = two spin components x two flavors; all 16 CAR states.
    alpha, beta = X, Z
    mass, momentum, energy = 4.0, 3.0, 5.0
    h = momentum*alpha+mass*beta
    r, v = 4/5, 3/5
    u = np.array([3.,1.], complex)/np.sqrt(10)
    w = (beta-r*I2)@u/v
    assert np.linalg.norm(h@u-energy*u) < 1e-14
    assert np.linalg.norm(h@w+energy*w) < 1e-14
    assert abs(np.vdot(u,w)) < 1e-14
    assert abs(np.linalg.norm(w)-1) < 1e-14
    e0, e1 = np.eye(2, dtype=complex)
    f0, f1 = np.kron(u,e0), np.kron(u,e1)
    hm0, hm1 = np.kron(w,e0), np.kron(w,e1)
    ann = car1128.annihilators(4)
    create = lambda f: car1128.creation(f, ann)
    ident = np.eye(16, dtype=complex)
    car_residual = max(opnorm(ann[i]@ann[j].conj().T +
                             ann[j].conj().T@ann[i] -
                             (ident if i == j else 0*ident))
                       for i in range(4) for j in range(4))
    assert car_residual < 1e-14

    empty = np.zeros(16, complex)
    empty[0] = 1
    sea = create(hm0)@create(hm1)@empty
    assert abs(np.linalg.norm(sea)-1) < 1e-14
    psi1 = create(f1)@sea
    assert abs(np.linalg.norm(psi1)-1) < 1e-14
    h4 = np.kron(h,I2)
    scalar = np.kron(beta,X)
    vector = np.kron(I2,X)
    qscalar = car1128.quadratic(scalar, ann)
    qvector = car1128.quadratic(vector, ann)
    normal_H = car1128.quadratic(h4, ann)+2*energy*ident
    assert np.linalg.norm(normal_H@sea) < 2e-14
    assert np.linalg.eigvalsh(normal_H)[0] > -2e-14

    probe = create(hm1)@create(f0).conj().T
    n0 = create(f0)@create(f0).conj().T
    n1 = create(f1)@create(f1).conj().T
    both = n0@n1
    pair_read = abs(np.vdot(psi1,probe@qscalar@psi1))
    pair_leakage = float(np.linalg.norm(both@qscalar@psi1))
    assert abs(pair_read-v) < 2e-14
    assert pair_leakage >= v-2e-14
    assert abs(np.vdot(sea,qscalar@sea)) < 1e-14
    scalar_variance = float(np.linalg.norm(qscalar@sea)**2)
    assert abs(scalar_variance-2*v*v) < 2e-14

    pp = np.kron((I2+h/energy)/2,I2)
    pm = np.eye(4)-pp
    scalar_off = opnorm(pp@scalar@pm)
    vector_off = opnorm(pp@vector@pm)
    assert abs(scalar_off-v) < 2e-14
    assert vector_off < 2e-14
    assert np.linalg.norm(qvector@sea) < 2e-14
    vector_leakage = float(np.linalg.norm(both@qvector@psi1))
    assert vector_leakage < 2e-14

    es, ws = np.linalg.eigh(normal_H)
    inv = np.zeros_like(es)
    inv[es>1e-9] = 1/es[es>1e-9]
    Hinv = (ws*inv)@ws.conj().T
    fluctuation = qscalar@sea
    response = float((2*np.vdot(fluctuation,Hinv@fluctuation)).real)
    assert abs(response-2*v*v/energy) < 2e-14
    return {
        "pass": True,
        "scope": "closed center spin sector; all four-mode CAR occupation sectors",
        "reuse": "1128/check.py annihilators, creation, quadratic",
        "dimension": 16, "Cauchy_modes": 4,
        "mass": mass, "momentum_z": momentum, "center_energy": energy,
        "CAR_residual": car_residual,
        "normal_sea_energy_residual": float(np.linalg.norm(normal_H@sea)),
        "prepared_code_probe_amplitude": float(pair_read),
        "two_code_occupation_leakage": pair_leakage,
        "scalar_vacuum_variance": scalar_variance,
        "scalar_positive_negative_block_norm": scalar_off,
        "vector_positive_negative_block_norm": vector_off,
        "vector_vacuum_source_norm": float(np.linalg.norm(qvector@sea)),
        "vector_two_code_occupation_leakage": vector_leakage,
        "finite_center_static_response": response,
        "not_a_spatial_Dirac_dynamics_approximation": True,
    }


def three(a,b,c):
    return np.kron(np.kron(a,b),c)


def group_cross():
    # Exactly defined algebraic diagnostic. Mediator Hilbert space here is
    # vacuum plus two labels, not a certified truncation of a physical field.
    delta_a, delta_b = 5.0, 7.0
    ha = np.diag([0.,0.,delta_a,delta_a])
    hb = np.diag([0.,0.,delta_b,delta_b])
    frequencies = [F(2),F(3)]
    rho_a = [F(1),F(1,2)]
    rho_b = [F(3,4),-F(2,3)]
    cross_coefficient = sum((ra*rb/(w*w) for ra,rb,w
                             in zip(rho_a,rho_b,frequencies)),F(0))
    assert cross_coefficient == F(65,432)
    hc = np.diag([0.,float(frequencies[0]),float(frequencies[1])])
    phi_a = np.zeros((3,3),complex)
    phi_b = np.zeros((3,3),complex)
    for j,(ra,rb,w) in enumerate(zip(rho_a,rho_b,frequencies),1):
        phi_a[0,j]=phi_a[j,0]=float(ra)/np.sqrt(2*float(w))
        phi_b[0,j]=phi_b[j,0]=float(rb)/np.sqrt(2*float(w))
    i4, i3 = np.eye(4), np.eye(3)
    h0 = three(ha,i4,i3)+three(i4,hb,i3)+three(i4,i4,hc)
    diagonal = np.diag(h0).real
    pin = np.zeros((48,4),complex)
    for a in range(2):
        for b in range(2):
            pin[(a*4+b)*3,a*2+b]=1
    assert opnorm(pin.conj().T@pin-np.eye(4)) < 1e-15
    rinv = np.zeros_like(diagonal)
    rinv[diagonal>0]=1/diagonal[diagonal>0]
    resolvent = np.diag(rinv)
    assert opnorm(h0@pin) < 1e-14

    cases = [
        (F(1,5), np.array([[1,1],[0,1]],complex)),
        (F(9,20),np.array([[1,0.5j],[1/3,-1]],complex)),
    ]
    target = -float(cross_coefficient)*np.kron(P1,X)
    rows = []
    outputs = []
    for kk, mat in cases:
        kappa = float(kk)
        mb = mat.conj().T
        ja = np.block([[P1,kappa*mat],[kappa*mat.conj().T,P1]])
        jb = np.block([[X,kappa*mb],[kappa*mb.conj().T,X]])
        assert opnorm(ja-ja.conj().T)+opnorm(jb-jb.conj().T)<1e-14
        va, vb = three(ja,i4,phi_a),three(i4,jb,phi_b)
        cross = -pin.conj().T@(va@resolvent@vb+vb@resolvent@va)@pin
        self_a = -pin.conj().T@va@resolvent@va@pin
        self_b = -pin.conj().T@vb@resolvent@vb@pin
        theory_a = np.zeros((2,2),complex)
        theory_b = np.zeros((2,2),complex)
        for ra,rb,w in zip(rho_a,rho_b,frequencies):
            wf = float(w)
            theory_a -= float(ra*ra)/(2*wf)*(
                P1@P1/wf+kappa*kappa*(mat@mat.conj().T)/(delta_a+wf))
            theory_b -= float(rb*rb)/(2*wf)*(
                X@X/wf+kappa*kappa*(mb@mb.conj().T)/(delta_b+wf))
        a_formula = opnorm(self_a-np.kron(theory_a,I2))
        b_formula = opnorm(self_b-np.kron(I2,theory_b))
        cross_error = opnorm(cross-target)
        commutator = opnorm(theory_a@P1-P1@theory_a)
        assert cross_error < 3e-14
        assert a_formula < 3e-14 and b_formula < 3e-14
        assert commutator > 1e-4
        full_second = -pin.conj().T@(va+vb)@resolvent@(va+vb)@pin
        decomposition = opnorm(full_second-cross-self_a-self_b)
        assert decomposition < 3e-14
        rows.append({
            "kappa": str(kk),
            "M_real": mat.real.tolist(), "M_imag": mat.imag.tolist(),
            "cross_identity_residual": cross_error,
            "self_A_formula_residual": a_formula,
            "self_B_formula_residual": b_formula,
            "self_A_commutator_with_a": commutator,
            "self_A_eigenvalues": np.linalg.eigvalsh(theory_a).tolist(),
            "self_B_eigenvalues": np.linalg.eigvalsh(theory_b).tolist(),
            "full_second_order_decomposition_residual": decomposition,
        })
        outputs.append((cross,self_a,self_b))
    cross_change = opnorm(outputs[1][0]-outputs[0][0])
    self_a_change = opnorm(outputs[1][1]-outputs[0][1])
    self_b_change = opnorm(outputs[1][2]-outputs[0][2])
    assert cross_change < 3e-14
    assert self_a_change > 1e-3 and self_b_change > 1e-3
    return {
        "pass": True,
        "scope": "exact second-order algebra only; not field cutoff or full-time record",
        "dimension": 48,
        "material_gaps": [delta_a,delta_b],
        "mediator_frequencies": [str(w) for w in frequencies],
        "rho_A": [str(x) for x in rho_a], "rho_B": [str(x) for x in rho_b],
        "cross_coefficient": fraction(cross_coefficient),
        "cases": rows,
        "cross_change_between_cases": cross_change,
        "self_A_change_between_cases": self_a_change,
        "self_B_change_between_cases": self_b_change,
        "local_self_terms_retained": True,
        "free_packet_not_claimed_to_have_isolated_low_energy_band": True,
    }


def dependency_hashes():
    files = [
        CAR_PATH,
        ARCHIVE/"1146"/"packet_bound.md",
        HERE/"after1146_material_mediator_screen.md",
        HERE/"after1146_mediator_domain_screen.md",
    ]
    return {p.relative_to(ARCHIVE).as_posix():
            hashlib.sha256(p.read_bytes()).hexdigest() for p in files}


def run():
    groups = [group_rational(),group_car(),group_cross()]
    return {
        "kind": "non_numbered_qualification_after_1146",
        "formal_round_increment": 0,
        "scientific_check_count_increment": 0,
        "check_groups": 3,
        "all_passed": all(g["pass"] for g in groups),
        "dependency_hashes": dependency_hashes(),
        "groups": groups,
    }


def compare(old,new,path="root"):
    if isinstance(new,dict):
        assert isinstance(old,dict) and old.keys()==new.keys(), path
        for key in new:
            compare(old[key],new[key],path+"."+key)
    elif isinstance(new,list):
        assert isinstance(old,list) and len(old)==len(new), path
        for i,(a,b) in enumerate(zip(old,new)):
            compare(a,b,path+f"[{i}]")
    elif isinstance(new,float):
        assert isinstance(old,(int,float)) and np.isclose(
            old,new,rtol=1e-10,atol=1e-12), (path,old,new)
    else:
        assert old==new, (path,old,new)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    result=run()
    if args.write_results:
        RESULT.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
        status="explicitly_saved"
    else:
        if RESULT.exists():
            compare(json.loads(RESULT.read_text(encoding="utf-8")),result)
            status="saved_results_match"
        else:
            status="recomputed_without_saved_results"
    print(json.dumps({
        "check_groups":result["check_groups"],"all_passed":result["all_passed"],
        "status":status,"formal_round_increment":0,
        "scientific_check_count_increment":0,
        "pair_amplitude":result["groups"][1]["prepared_code_probe_amplitude"],
        "cross_case_change":result["groups"][2]["cross_change_between_cases"],
        "self_A_change":result["groups"][2]["self_A_change_between_cases"],
    },ensure_ascii=False,indent=2))
