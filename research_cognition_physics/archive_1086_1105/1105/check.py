"""Round 1105: exact finite recoil ledger; no whole-model certification."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SOURCES = [
    "research_cognition_physics/archive_1086_/1095/proof.md",
    "research_cognition_physics/archive_1086_/1103/proof.md",
    "research_cognition_physics/archive_1086_/1104/proof.md",
    "research_cognition_physics/archive_935_955/research_note_941.md",
    "research_cognition_physics/archive_370_428/research_note_421.md",
    "research_cognition_physics/archive_1086_/_admission/after_1102_atomic_macro_scope.md",
]


def cs(k):
    a = F(2) ** k
    return (a + 1/a)/2, (a - 1/a)/2


def exact_string(x):
    return str(x)


def compute():
    groups = {}
    m, parent = F(2), F(5)
    hin, pin = [parent*cs(k)[0] for k in (0, 1)], [parent*cs(k)[1] for k in (0, 1)]
    hout = [m*(cs(a)[0]+cs(b)[0]) for a in (1, 2) for b in (-1, 0)]
    pout = [m*(cs(a)[1]+cs(b)[1]) for a in (1, 2) for b in (-1, 0)]
    assert hin == [F(5), F(25,4)] and pin == [0, F(15,4)]
    assert hout == [5, F(9,2), F(27,4), F(25,4)]
    assert pout == [0, F(3,2), F(9,4), F(15,4)]
    assert [hout[i] for i in (0,3)] == hin
    assert [pout[i] for i in (0,3)] == pin
    # A specified passive boost of the whole ledger is a further exact check.
    # This does not certify an active reference preparation.
    gamma, beta = F(5,4), F(3,5)
    boosted_h = [gamma*(h-beta*p) for h,p in zip(hout,pout)]
    boosted_p = [gamma*(p-beta*h) for h,p in zip(hout,pout)]
    assert [boosted_h[i] for i in (0,3)] == [gamma*(h-beta*p) for h,p in zip(hin,pin)]
    assert [boosted_p[i] for i in (0,3)] == [gamma*(p-beta*h) for h,p in zip(hin,pin)]
    groups["exact_vertex_intertwining"] = {
        "H_in": list(map(exact_string,hin)), "P_in": list(map(exact_string,pin)),
        "H_out": list(map(exact_string,hout)), "P_out": list(map(exact_string,pout)),
        "exact_energy_and_momentum_intertwining": True,
        "passive_ledger_boost_checked": True,
        "actual_all_reference_preparation_claimed": False,
    }

    w = np.zeros((4,2),complex)
    w[0,0] = w[3,1] = 1
    joint_w = np.kron(w,np.eye(4))  # motion, Q, passive R
    rng = np.random.default_rng(1105)
    x = rng.normal(size=(8,8)) + 1j*rng.normal(size=(8,8))
    rho = x@x.conj().T
    rho /= np.trace(rho)
    out = joint_w@rho@joint_w.conj().T
    recovered = joint_w.conj().T@out@joint_w
    before_qr = np.einsum("aiaj->ij",rho.reshape(2,4,2,4))
    after_qr = np.einsum("aiaj->ij",out.reshape(4,4,4,4))
    identity_error = float(np.max(np.abs(joint_w.conj().T@joint_w-np.eye(8))))
    recovery_error = float(np.max(np.abs(recovered-rho)))
    payload_error = float(np.max(np.abs(before_qr-after_qr)))
    plus = np.ones(2)/np.sqrt(2)
    bell = w@plus
    full = np.outer(bell,bell.conj()).reshape(2,2,2,2)
    reduced = np.einsum("abcb->ac",full)
    clean_target = np.outer(plus,plus)
    coherence_loss = float(np.linalg.norm(reduced-clean_target,ord="nuc")/2)
    assert max(identity_error,recovery_error,payload_error) < 1e-12
    assert np.isclose(coherence_loss,.5,atol=1e-12)
    groups["complete_unknown_input_and_recoil"] = {
        "joint_input_dimension":8, "joint_output_dimension":16,
        "isometry_residual":identity_error, "full_recovery_residual":recovery_error,
        "unknown_Q_R_marginal_residual":payload_error,
        "motion_coherence_loss_if_recoil_discarded":coherence_loss,
        "all_passive_references_by_operator_identity": True,
        "whole_unknown_motion_worldline_chain_certified": False,
    }

    q, mu = F(5,2), F(1)
    ledgers, endpoints = [], []
    delta, tau = F(1,4), F(2)
    for n in range(1,13):
        m0 = mu*q**n
        sides = [(mu*q**(n-j), j-2, 2**(n-j)) for j in range(1,n+1)]
        objects = sides+[(mu,n,1)]
        e = sum(mass*cs(k)[0] for mass,k,_ in objects)
        p = sum(mass*cs(k)[1] for mass,k,_ in objects)
        cores = sum(count for _,_,count in objects)
        stored = m0-cores*mu
        retained = sum(mass-count*mu for mass,_,count in objects)
        kinetic = sum(mass*(cs(k)[0]-1) for mass,k,_ in objects)
        assert e == m0 and p == 0 and cores == 2**n
        assert stored == retained+kinetic and retained >= 0 and kinetic > 0
        cn,sn = cs(n)
        speed = sn/cn
        assert speed == F(4**n-1,4**n+1) and 0 < speed < 1
        ledgers.append({
            "N":n, "cores":cores, "initial_rest_mass":str(m0),
            "initial_stored_energy_c_equals_1":str(stored),
            "unreleased_side_energy":str(retained), "kinetic_energy":str(kinetic),
            "forward_speed_over_c":str(speed),
            "exact_total_energy":True, "exact_total_momentum_zero":True,
        })
        events = []
        t,xp = F(0),F(0)
        for k in range(n):
            ck,sk = cs(k)
            dt,dx = delta*ck,delta*sk
            assert dt*dt-dx*dx == delta*delta
            t,xp = t+dt,xp+dx
            events.append((t,xp))
        assert t == delta*(F(2)**n+1-F(2)**(1-n))/2
        assert xp == delta*(F(2)**n-3+F(2)**(1-n))/2
        dt,dx = tau*cn,tau*sn
        assert dt*dt-dx*dx == tau*tau
        T = t+dt
        side_positions = []
        for j,(mass,k,count) in enumerate(sides,1):
            tj,xj = events[j-1]
            ck,sk = cs(k)
            elapsed = T-tj
            pos = xj+sk/ck*elapsed
            assert elapsed > 0 and elapsed**2-(pos-xj)**2 > 0
            side_positions.append(str(pos))
        margin = dt-dx
        assert margin == tau/F(2)**n > 0
        et,ex = margin/8,margin/8
        robust_margin = (dt-et)-(dx+ex)
        assert robust_margin == 3*margin/4 > 0
        endpoints.append({
            "N":n, "preparation_t":str(t), "preparation_x":str(xp),
            "co_moving_output_dt":str(dt), "co_moving_output_dx":str(dx),
            "timelike_margin":str(margin),
            "conditional_each_readout_error":str(et),
            "remaining_margin":str(robust_margin),
            "retained_side_positions_at_output_time":side_positions,
        })
    groups["finite_tree_resource_ledger"] = {
        "chi":"log(2)", "mu":1, "c":1, "rows":ledgers,
        "uniform_N_budget_asserted":False,
        "live_clock_and_control_energy_certified":False,
    }
    groups["actual_finite_worldline_endpoints"] = {
        "proper_preparation_delay":str(delta), "proper_co_moving_time":str(tau),
        "rows":endpoints,
        "Lorentz_position_and_clock_laws_input":True,
        "fixed_receiver_takeover_certified":False,
        "uniform_all_N_error_budget_asserted":False,
    }
    return {
        "round":1105, "status":"PASS_PARTIAL_PREPARATION_CONTRACT",
        "groups":groups,
        "source_sha256": {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCES},
        "both_branches_share_only_K":True,
        "new_universal_axioms":0, "new_scientific_experiments":0,
        "scientific_count_increment":0, "scientific_count_total":3860,
        "full_joint_model_certified":False, "full_consensus_conjecture_refuted":False,
        "Lorentz_structure_derived_from_cognition":False,
        "goal_complete":False, "empirical_data":False,
    }


def compare(a,b):
    if isinstance(a,dict):
        assert a.keys() == b.keys()
        for key in a:
            compare(a[key],b[key])
    elif isinstance(a,list):
        assert len(a) == len(b)
        for x,y in zip(a,b):
            compare(x,y)
    elif isinstance(a,float):
        assert np.isclose(a,b,atol=1e-12,rtol=1e-12), (a,b)
    else:
        assert a == b, (a,b)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write",action="store_true")
    args = parser.parse_args()
    result = compute()
    if args.write:
        with (HERE/"results.json").open("x",encoding="utf-8") as stream:
            json.dump(result,stream,ensure_ascii=False,indent=2)
            stream.write("\n")
    else:
        compare(json.loads((HERE/"results.json").read_text(encoding="utf-8")),result)
    print(json.dumps({"status":result["status"],"groups":len(result["groups"]),
                      "saved_results_match":not args.write}))
