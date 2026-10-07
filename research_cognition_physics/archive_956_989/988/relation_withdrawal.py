"""988: distinguish stopping formation from removing future influence.

Finite conditional-operation contract; not an autonomous feedback Hamiltonian.
The 4D material is the exact invariant sector inherited from round 968.
"""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
ROOT = HERE.parents[2]
OUT = HERE / "relation_withdrawal_results.json"


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def kron(*xs):
    out = np.ones((1, 1))
    for x in xs:
        out = np.kron(out, x)
    return out


def norm(x):
    return float(np.linalg.norm(x, 2))


def comm(a, b):
    return a @ b - b @ a


def expectation(psi, op):
    return float(np.vdot(psi, op @ psi).real)


def propagator(H, time):
    values, V = np.linalg.eigh(H)
    gram = norm(V.conj().T @ V - np.eye(len(H)))
    residual = norm(H @ V - V * values)
    # Conservative short-time arithmetic allowance, standard double model.
    # The defining h and Q are rational; F contains sqrt(3). Dimension <= 256.
    budget = gram + abs(time) * (residual + 1e-12) * norm(V) + 1e-10
    assert budget < 1e-8
    U = (V * np.exp(-1j * time * values)) @ V.conj().T
    return U, dict(gram_defect=gram, eigen_residual=residual,
                   operator_error_budget=budget)


def run():
    old = load("material968", STAGE / "968/internal_relay.py")
    _, h, q, f, _, material_checks = old.material()
    I = np.eye(4)
    basis = np.eye(4)
    T, dplus, dminus = basis[:, 3], basis[:, 1], basis[:, 2]
    P = np.diag([1., 1., 1., 0.])
    kappa, eta, formation_time, future_time = .2, .05, 20., 5.
    Hfree = kron(h, I, I) + kron(I, h, I) + kron(I, I, h)
    contact = kappa * (kron(q, q, I) + kron(I, q, q))
    FB = kron(I, f, I)
    H0 = Hfree + contact
    H1 = H0 + eta * FB
    PB = kron(I, P, I)
    EC = kron(I, I, np.outer(dminus, dminus))
    XA = I.copy()
    XA[:, [1, 3]] = XA[:, [3, 1]]
    X = kron(XA, I, I)
    # Exact rational counterexample: H0 = K/5. Noncommutation of the
    # Heisenberg receiving effect with a source-only unitary proves that
    # H0 is not a universally severed A->C channel. No floating time is used.
    K = np.rint(5*H0).astype(np.int64)
    assert np.max(abs(5*H0-K)) < 1e-14
    nested = EC.astype(np.int64)
    exact_rows = []
    source_int = X.astype(np.int64)
    for order in range(1,4):
        nested = K@nested-nested@K
        commutator = nested@source_int-source_int@nested
        exact_rows.append(dict(order=order,nonzero_entries=int(np.count_nonzero(commutator)),
            entry_17_38=int(commutator[17,38]),denominator=5**order))
    assert exact_rows[0]["nonzero_entries"] == exact_rows[1]["nonzero_entries"] == 0
    assert exact_rows[2]["entry_17_38"] == 1
    initial = np.kron(np.kron(T, T), dplus)
    Uform, cform = propagator(H1, formation_time)
    Ufuture, cfuture = propagator(H0, future_time)
    formed = Uform @ initial
    before_p = expectation(formed, PB)
    no_park = []
    for action in (np.eye(64), X):
        start = action @ formed
        end = Ufuture @ start
        no_park.append(dict(receiver_probability=expectation(end, EC),
            active_relay_before=expectation(start, PB),
            active_relay_after=expectation(end, PB),
            intervention_work=expectation(start, H0)-expectation(formed, H0),
            energy_drift=expectation(end, H0)-expectation(start, H0)))
    contrast = abs(no_park[0]["receiver_probability"]-no_park[1]["receiver_probability"])
    # A probability difference involves two histories of two propagations.
    contrast_budget = 4*(cform["operator_error_budget"]+cfuture["operator_error_budget"])+1e-9

    # Full 4-site carrier A,B,C,D. D is a declared isolated, same-material blank.
    Hfree4 = sum(kron(*[h if j == i else I for j in range(4)]) for i in range(4))
    contact4 = kappa*(kron(q,q,I,I)+kron(I,q,q,I))
    Hafter = Hfree4 + contact4
    S = np.zeros((256,256))
    embedding = np.zeros((256,64))
    V = np.zeros((256,64))
    for a in range(4):
        for b in range(4):
            for c in range(4):
                col = 16*a+4*b+c
                embedding[64*a+16*b+4*c+3,col] = 1
                V[64*a+48+4*c+b,col] = 1
                for d in range(4):
                    S[64*a+16*d+4*c+b,64*a+16*b+4*c+d] = 1
    checks = dict(swap_isometry=norm(S.T@S-np.eye(256)),
        swap_recovers=norm(S@S-np.eye(256)),
        exact_state_relocation=norm(S@embedding-V),
        arbitrary_input_isometry=norm(V.T@V-np.eye(64)),
        dark_sector_intertwining=norm(Hafter@V-V@Hfree),
        unchanged_local_energy_operator=norm(S.T@Hfree4@S-Hfree4),
        contact_work_identity=norm(V.T@Hafter@V-H0+contact),
        stopped_update_conserves_active_population=norm(comm(H0,PB)),
        formation_changes_active_population=norm(comm(H1,PB)))
    assert max(v for key,v in checks.items() if key != "formation_changes_active_population") < 1e-13
    assert checks["formation_changes_active_population"] > .02
    Upark, cpark = propagator(Hafter, future_time)
    E4 = kron(I,I,np.outer(dminus,dminus),I)
    X4 = kron(XA,I,I,I)
    parked = V@formed
    park_rows = []
    for action in (np.eye(256),X4):
        state=action@parked
        end=Upark@state
        park_rows.append(dict(receiver_probability=expectation(end,E4),
            intervention_work=expectation(state,Hafter)-expectation(parked,Hafter),
            energy_drift=expectation(end,Hafter)-expectation(state,Hafter)))
    park_contrast=abs(park_rows[0]["receiver_probability"]-park_rows[1]["receiver_probability"])
    assert park_contrast < 1e-12
    # Check recovery on the whole 64D input space, including arbitrary references.
    Ufree, cfree = propagator(Hfree, future_time)
    transport_defect=norm(Upark@V-V@Ufree)
    assert transport_defect < 1e-12
    work = dict(stop_update=expectation(formed,H0)-expectation(formed,H1),
        expected_stop_update=-eta*expectation(formed,FB),
        park=expectation(parked,Hafter)-expectation(formed,H0),
        expected_park=-expectation(formed,contact),
        parking_work_operator_norm=norm(contact),
        analytic_parking_work_abs_bound=2*kappa,
        analytic_update_switch_work_abs_bound=3*eta/4)
    assert abs(work["park"]-work["expected_park"])<1e-13
    assert abs(work["stop_update"]-work["expected_stop_update"])<1e-13
    assert max(abs(row["energy_drift"]) for row in no_park+park_rows)<1e-12
    # A complete process must fund these operation exchanges; no battery here.
    assert contrast > contrast_budget
    files=[STAGE/"968/internal_relay.py",STAGE/"968/internal_relay_results.json",
        STAGE/"research_note_932.md",STAGE/"research_note_962.md",
        STAGE/"research_note_968.md",STAGE/"research_note_972.md",
        STAGE/"research_note_987.md",STAGE/"987/drafts/mechanism_map_v0_3.md",
        HERE/"drafts/STATUS.md",HERE/"drafts/withdrawal_decision.md"]
    return dict(round=988,all_scientific_checks_passed=True,
        parameters=dict(U=1,v=.1,kappa=kappa,eta=eta,
            formation_time=formation_time,future_time=future_time,
            initial_material_labels=["T","T","D+"],receiver_effect="D- population",
            control_actions=["identity on A","swap D+ and T on A"]),
        reuse_material_checks=material_checks,operator_checks=checks,
        exact_heisenberg_source_commutators=exact_rows,
        formed_active_population=before_p,stop_update_only=no_park,
        new_intervention_contrast=contrast,numerical_contrast_budget=contrast_budget,
        positive_contrast_lower=contrast-contrast_budget,
        park_then_stop=park_rows,park_numerical_contrast=park_contrast,
        park_exact_future_contrast=0,all_input_free_transport_defect=transport_defect,
        certificates=dict(formation=cform,future=cfuture,parked=cpark,free=cfree),
        work_ledger=work,
        resources=dict(additional_same_material_blank_dimension=4,
            historical_input_dimension=64,exact_relocation=True,
            blank_reusable_without_recovery=False,
            includes_autonomous_control_and_battery=False),
        decision=dict(reject_stop_formation_as_universal_withdrawal=True,
            retain_internal_relocation_as_conditional_withdrawal=True,
            direct_contact_control_is_an_alternative=True,
            no_claim_unique_mechanism=True,full_goal_completed=False),
        source_hashes={str(p.relative_to(ROOT)):sha(p) for p in files})


def compare(a,b,path=""):
    if isinstance(a,dict):
        assert a.keys()==b.keys(),path
        for k in a:compare(a[k],b[k],path+"/"+k)
    elif isinstance(a,list):
        assert len(a)==len(b),path
        for i,(x,y) in enumerate(zip(a,b)):compare(x,y,path+f"/{i}")
    elif isinstance(a,float):
        assert math.isclose(a,b,rel_tol=1e-8,abs_tol=1e-11),(path,a,b)
    else:assert a==b,(path,a,b)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--write",action="store_true")
    args=parser.parse_args()
    result=run()
    if args.write:
        with OUT.open("x",encoding="utf-8") as f:
            json.dump(result,f,ensure_ascii=False,indent=2);f.write("\n")
    else:compare(result,json.loads(OUT.read_text("utf-8")))
    print(json.dumps({k:result[k] for k in ("round","all_scientific_checks_passed",
        "formed_active_population","new_intervention_contrast","positive_contrast_lower",
        "park_numerical_contrast","work_ledger")},ensure_ascii=False,indent=2))
