"""Round 1320. Complete finite CAR model and conditional local-current dictionary.
Default/--check reads only. --write-new creates results.json exclusively.
"""
from __future__ import annotations
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
CHECKS = []


def check(name, value, tol=2e-9):
    err = float(np.max(np.abs(value))) if np.size(value) else 0.
    CHECKS.append({"name": name, "error": err, "tolerance": tol})
    if not np.isfinite(err) or err > tol:
        raise AssertionError((name, err, tol))


def truth(name, value):
    check(name, 0 if bool(value) else 1, 0)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def exponential(h, t):
    e, u = np.linalg.eigh(h)
    return (u*np.exp(-1j*t*e))@u.conj().T


def wedge(u):
    pairs = np.array(list(itertools.combinations(range(len(u)), 2)))
    i, j = pairs[:, 0, None], pairs[:, 1, None]
    k, l = pairs[None, :, 0], pairs[None, :, 1]
    return u[i, k]*u[j, l]-u[i, l]*u[j, k]


def calculate():
    CHECKS.clear()
    rng = np.random.default_rng(1320)
    pairs = list(itertools.combinations(range(28), 2))
    dim = 1+len(pairs)+3
    positions = {(1 << p)+(1 << q): k+1 for k, (p, q) in enumerate(pairs)}

    def create(mask, mode):
        if mask & (1 << mode): return mask, 0
        return mask | (1 << mode), (-1)**((mask & ((1 << mode)-1)).bit_count())

    def pair(p, q):
        mask, phase = create(0, q)
        mask, phase2 = create(mask, p)
        out = np.zeros(dim, complex)
        if phase*phase2: out[positions[mask]] = phase*phase2
        return out

    eps = np.array([[0, 1], [-1, 0]])
    qcols = []
    for i, j in [(0, 1), (0, 2), (1, 2)]:
        out = np.zeros(dim, complex)
        for a in range(2):
            for alpha in range(2):
                out += eps[a, 1-a]*eps[alpha, 1-alpha]*pair(4*i+2*a+alpha, 4*j+2*(1-a)+1-alpha)/2
        qcols.append(out)
    c0 = np.array([[0, 0, 1], [0, -1, 0], [1, 0, 0]], complex)
    wa = np.column_stack(qcols)@c0.T
    ccols = []
    for k in range(3):
        out = np.zeros(dim, complex)
        for a in range(2):
            for alpha in range(2):
                out += eps[a, 1-a]*eps[alpha, 1-alpha]*pair(12+4*k+2*a+alpha, 24+2*(1-a)+1-alpha)/2
        ccols.append(out)
    wc = np.column_stack(ccols)
    ws = np.eye(dim, dtype=complex)[:, -3:]
    pa, pc, ps = wa@wa.conj().T, wc@wc.conj().T, ws@ws.conj().T
    w9 = np.column_stack((wa, ws, wc))
    ta, tc = ws@wa.conj().T, ws@wc.conj().T
    nq = np.array([0]+[sum(i < 12 for i in p) for p in pairs]+[0]*3)
    nqb = np.array([0]+[sum(12 <= i < 24 for i in p) for p in pairs]+[0]*3)
    nlb = np.array([0]+[sum(i >= 24 for i in p) for p in pairs]+[0]*3)
    ns = np.array([0]*379+[1]*3)
    e_q, e_qb, e_lb, es = .7, .6, .8, 1.4
    ga, gc = .3+.4j, .4-.3j
    pars = np.array([e_q, e_qb, e_lb, es, ga.real, ga.imag, gc.real, gc.imag])

    def father(x=pars):
        eq, eqb, elb, esm, ar, ai, cr, ci = x
        aa, cc = ar+1j*ai, cr+1j*ci
        return np.diag(eq*nq+eqb*nqb+elb*nlb+esm*ns)+aa*ta+aa.conjugate()*ta.conj().T+cc*tc+cc.conjugate()*tc.conj().T

    h = father()
    h0 = np.diag(e_q*nq+e_qb*nqb+e_lb*nlb+es*ns)
    charge = np.diag(nq+nqb+nlb+2*ns)
    ff = np.zeros((dim, dim)); ff[1:379, 1:379] = np.eye(378)
    dark = ff-pa-pc
    check("CAR_QQ_orthonormal", wa.conj().T@wa-np.eye(3))
    check("CAR_antiQL_orthonormal", wc.conj().T@wc-np.eye(3))
    check("all_nine_codes_orthonormal", w9.conj().T@w9-np.eye(9))
    truth("complete_382_dimension", dim == 382 and len(pairs) == 378)
    check("complete_dark_dimension_372", np.trace(dark)-372)
    check("dark_not_discarded", (h-h0)@dark)
    check("K_conserved", h@charge-charge@h)
    check("Hamiltonian_Hermitian", h-h.conj().T)
    check("vacuum_inert", h[:, 0])
    # QQ combined color-weak Yukawa coefficient symmetric, not flavor antisymmetric.
    color_eps = np.zeros((3, 3)); color_eps[0, 1] = 1; color_eps[1, 0] = -1
    y_combined = np.kron(color_eps, eps)
    check("QQ_combined_Yukawa_symmetric", y_combined-y_combined.T)
    check("QQ_full_CAR_coefficient_antisymmetric", np.kron(y_combined, eps)+np.kron(y_combined, eps).T)
    truth("same_flavor_QQ_nonzero", np.linalg.norm(wa[:, 2]) > .99)

    # Complete action of a common internal/weak/spin change on retained states.
    raw = rng.normal(size=(3, 3))+1j*rng.normal(size=(3, 3))
    u, rr = np.linalg.qr(raw); u = u@np.diag(np.diag(rr)/np.abs(np.diag(rr)))
    weak = np.array([[np.cos(.21), np.sin(.21)], [-np.sin(.21), np.cos(.21)]])
    spin = np.diag(np.exp(1j*np.array([.17, -.17])))
    f = np.zeros((28, 28), complex)
    f[:12, :12] = np.kron(u, np.kron(weak, spin))
    f[12:24, 12:24] = np.kron(u.conj(), np.kron(weak.conj(), spin.conj()))
    f[24:, 24:] = np.linalg.det(u)*np.kron(weak.conj(), spin.conj())
    rep = np.zeros((dim, dim), complex); rep[0, 0] = 1
    rep[1:379, 1:379] = wedge(f); rep[-3:, -3:] = np.linalg.det(u)*u.conj()
    b_rep = np.linalg.det(u)*u.conj()
    check("QQ_role_intertwiner", rep@wa-wa@b_rep)
    check("antiQL_role_intertwiner", rep@wc-wc@b_rep)
    check("same_parent_internal_weak_spin", rep@h-h@rep)
    fermion_r = np.array([1j]*12+[-1j]*16)
    rdiag = np.array([1]+[fermion_r[p]*fermion_r[q] for p, q in pairs]+[-1]*3)
    check("old_r_same_parent", rdiag[:, None]*h-h*rdiag[None, :])

    # Character conclusions and exact deletion examples.
    omega = np.exp(2j*np.pi/3)
    ub = np.diag([1, omega, omega**2])
    twist = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]], complex)
    check("Z3_full_rank_nontrivial_self_twist", twist@(omega*ub)-ub@twist)
    check("Z3_twist_full_rank", twist.conj().T@twist-np.eye(3))
    truth("Z3_character_not_equal", abs(omega-1) > 1)
    theta = .43
    ub1 = np.diag([np.exp(-1j*theta), 1, 1]); dchar = np.exp(1j*theta)
    partial = np.zeros((3, 3), complex); partial[0, 1] = 1
    check("connected_rank_one_nontrivial_intertwiner", partial@ub1-dchar*ub1@partial)
    truth("partial_intertwiner_not_full", np.linalg.matrix_rank(partial) == 1)
    capacities = []
    for w in range(1, 6):
        lower = 1-1/w
        capacities.append({"w": w, "C_dimension": 3*w, "roundtrip_lower_bound": lower})
        truth(f"capacity_w_{w}", (lower < .5) == (w == 1))
    subcode = np.zeros((6, 3)); subcode[:3, :] = np.eye(3)
    check("w2_selected_subcode_only", subcode.T@subcode-np.eye(3))
    truth("w2_rest_not_delivered", np.linalg.matrix_rank(subcode@subcode.T) == 3)

    def h3(a, c, delta=0):
        return np.array([[0, a.conjugate(), 0], [a, delta, c], [0, c.conjugate(), 0]], complex)

    e0 = 2*e_q
    check("actual_three_channel_parent", w9.conj().T@h@w9-e0*np.eye(9)-np.kron(h3(ga, gc), np.eye(3)))
    strength = np.sqrt(abs(ga)**2+abs(gc)**2)
    tstar = np.pi/strength
    ufull = exponential(h, tstar)
    check("balanced_S_returns_empty", ws.conj().T@ufull@wa)
    forward, reverse = wc.conj().T@ufull@wa, wa.conj().T@ufull@wc
    target = -2*gc.conjugate()*ga/strength**2*np.exp(-1j*e0*tstar)
    check("balanced_actual_full_transfer", forward-target*np.eye(3))
    check("balanced_actual_reverse", reverse-(-2*ga.conjugate()*gc/strength**2*np.exp(-1j*e0*tstar))*np.eye(3))
    check("balanced_complete_roundtrip", (reverse@forward).conj().T@(reverse@forward)-np.eye(3))
    darkchain = np.array([gc, 0, -ga])/strength
    check("actual_chain_dark_vector", h3(ga, gc)@darkchain)
    halfway = exponential(h3(ga, gc), tstar/2)[:, 0]
    check("S_not_exclusive_archive", abs(halfway[1])**2-.5)
    truth("endpoints_do_not_exit_midway", abs(halfway[0])**2+abs(halfway[2])**2 > .49)
    check("all_finite_dark_evolution", ufull@dark-np.diag(np.exp(-1j*tstar*np.diag(h0)))@dark)

    ratio_cases = []
    for label, a, c in [("balanced", ga, gc), ("admissible_imbalance", .5+0j, .4j), ("outside_contract", .5+0j, .1j)]:
        uu, vv = abs(a)**2, abs(c)**2
        gg = np.sqrt(uu+vv); tp = np.pi/gg
        pp = 4*uu*vv/(uu+vv)**2
        eps_tr, eps_rt = 1-pp, 1-pp**2
        evo = exponential(h3(a, c), tp)
        check(label+"_peak_probability", abs(evo[2, 0])**2-pp)
        check(label+"_S_empty", evo[1, 0])
        recovered = evo[0, 2]*evo[2, 0]
        check(label+"_two_pass_probability", abs(recovered)**2-pp**2)
        check(label+"_rt_tr_dictionary", eps_rt-(2*eps_tr-eps_tr**2))
        delta_ratio = (uu-vv)/(uu+vv)
        check(label+"_imbalance_exact", delta_ratio**2-(1-np.sqrt(max(0, 1-eps_rt))))
        phi = np.eye(3).ravel()/np.sqrt(3)
        ideal = np.zeros((12, 12), complex); ideal[:9, :9] = np.outer(phi, phi)
        output = pp**2*ideal; output[9:, 9:] = eps_rt*np.eye(3)/3
        check(label+"_partner_complete_erasure", np.sum(np.abs(np.linalg.eigvalsh(output-ideal)))/2-eps_rt)
        ratio_cases.append({"label": label, "u_over_v": uu/vv, "peak_probability": pp, "transfer_erasure": eps_tr, "roundtrip_erasure": eps_rt, "within_half_contract": eps_rt < .5})
    truth("unbalanced_common_nonempty", ratio_cases[1]["within_half_contract"])
    truth("strong_imbalance_failure", not ratio_cases[2]["within_half_contract"])
    # Deleting the simultaneous-constant-control rule invalidates a universal balance claim.
    sequential_a, sequential_c = .5+0j, .1j
    first = exponential(h3(sequential_a, 0j), np.pi/(2*abs(sequential_a)))
    second = exponential(h3(0j, sequential_c), np.pi/(2*abs(sequential_c)))
    sequential = second@first
    check("independent_portal_schedule_exact_transfer", abs(sequential[2, 0])**2-1)
    check("independent_schedule_S_empty", sequential[1, 0])
    for k, (delta, t) in enumerate([(0, .2), (.3, 1.1), (-.4, 2), (.7, 5)]):
        a, c = .5+0j, .4j
        pmax = 4*abs(a*c)**2/(abs(a)**2+abs(c)**2)**2
        pp = abs(exponential(h3(a, c, delta), t)[2, 0])**2
        truth(f"all_detuning_peak_upper_{k}", pp <= pmax+1e-12)

    # Whole sources, before restricting to the resonant parameter surface.
    sources = [np.diag(nq), np.diag(nqb), np.diag(nlb), np.diag(ns), ta+ta.conj().T, 1j*(ta-ta.conj().T), tc+tc.conj().T, 1j*(tc-tc.conj().T)]
    step = 2e-6
    for j in range(8):
        plus, minus = pars.copy(), pars.copy(); plus[j] += step; minus[j] -= step
        check(f"full_independent_parent_source_{j}", (father(plus)-father(minus))/(2*step)-sources[j], 3e-8)
    psi = wa@np.array([1, 1j, -1])/np.sqrt(3)
    out = ufull@psi
    expect = lambda op, state: float(np.vdot(state, op@state).real)
    check("full_energy_conserved", expect(h, out)-expect(h, psi))
    check("resonant_free_interaction_commute", h0@(h-h0)-(h-h0)@h0)
    terminal = pa@np.outer(out, out.conj())@pa+pc@np.outer(out, out.conj())@pc+ps@np.outer(out, out.conj())@ps
    check("terminal_mean_readout_energy", np.trace(h@terminal).real-expect(h, out))
    check("terminal_Nq", expect(np.diag(nq), out))
    check("terminal_NantiQ", expect(np.diag(nqb), out)-1)
    check("terminal_NantiL", expect(np.diag(nlb), out)-1)
    check("terminal_NS", expect(np.diag(ns), out))
    ephys = .31
    endpoint = np.column_stack((wa, wc))
    vrow = np.column_stack((ga*np.eye(3), gc*np.eye(3)))
    finite_sigma = -vrow.conj().T@vrow/(es-ephys)
    active3 = np.kron(h3(ga, gc), np.eye(3))+e0*np.eye(9)
    ends = [0, 1, 2, 6, 7, 8]; meds = [3, 4, 5]
    schur = active3[np.ix_(ends, ends)]-ephys*np.eye(6)-active3[np.ix_(ends, meds)]@np.linalg.inv(active3[np.ix_(meds, meds)]-ephys*np.eye(3))@active3[np.ix_(meds, ends)]
    check("same_parent_energy_Schur", schur-(e0-ephys)*np.eye(6)-finite_sigma)

    # Fixed spatial-overlap dictionary, not an identification of energy g and dimensionless lambda.
    kappa = np.array([2+0j, 3*np.exp(.2j)])
    lam = np.array([ga, gc])/kappa
    y_d, y_l = lam[0], lam[1].conjugate()
    mass = 1.4; m2 = mass**2
    def kernel(l=lam, mm=m2, momentum2=0):
        return np.outer(l.conj(), l)/(mm-momentum2)
    k0 = kernel()
    check("fixed_overlap_dictionary", kappa*lam-[ga, gc])
    check("field_Gram_Hermitian", k0-k0.conj().T)
    truth("field_Gram_positive", np.linalg.eigvalsh(k0)[0] > -1e-12)
    check("field_Gram_rank_one", np.linalg.det(k0))
    check("field_cross_modulus_identity", abs(k0[0, 1])**2-k0[0, 0]*k0[1, 1])
    check("QQQL_conjugation_dictionary", k0[1, 0]-y_l*y_d/m2)
    check("task_field_ratio", abs(kappa[0])**2*k0[0, 0]/(abs(kappa[1])**2*k0[1, 1])-abs(ga/gc)**2)
    check("exact_task_balanced_Gram_weight", abs(kappa[0])**2*k0[0, 0]-abs(kappa[1])**2*k0[1, 1])
    for j, momentum2 in enumerate([-.2*m2, .1*m2, .25*m2]):
        err = np.linalg.norm(kernel(momentum2=momentum2)-k0, 2)
        bound = abs(momentum2)/(m2-abs(momentum2))*np.linalg.norm(k0, 2)
        truth(f"mass_squared_finite_window_{j}", err <= bound+1e-12)
    dl = np.array([.17-.09j, -.04+.11j]); dm2 = .23
    differential = (np.outer(dl.conj(), lam)+np.outer(lam.conj(), dl))/m2-np.outer(lam.conj(), lam)*dm2/m2**2
    numeric = (kernel(lam+step*dl, m2+step*dm2)-kernel(lam-step*dl, m2-step*dm2))/(2*step)
    check("field_complete_source", numeric-differential, 3e-8)
    dk = np.array([.07+.03j, -.08+.01j])
    dg = kappa*dl+dk*lam
    numeric_g = ((kappa+step*dk)*(lam+step*dl)-(kappa-step*dk)*(lam-step*dl))/(2*step)
    check("variable_overlap_source_chain", numeric_g-dg, 3e-8)
    check("rest_scalar_energy_mass_source", (np.sqrt(m2+step*dm2)-np.sqrt(m2-step*dm2))/(2*step)-dm2/(2*mass), 3e-8)
    # Two independent mediators and direct contacts explicitly delete the rank-one equality.
    second = np.array([.13+.02j, -.03+.19j])
    k_two = k0+np.outer(second.conj(), second)/(2.1**2)
    truth("two_mediators_strict_Gram_inequality", np.linalg.det(k_two).real > 1e-7)
    contact = np.diag([-.2, 0])
    truth("free_contact_can_break_positivity", np.linalg.eigvalsh(k0+contact)[0] < 0)
    # This is a bound on a conditional projection, not a claimed full-QFT implementation.
    eta, total_time = 1e-4, 2*tstar
    projection_error_budget = total_time*eta
    truth("conditional_projection_budget_finite", projection_error_budget < .01)
    return {
        "round": 1320, "status": "author_checks_passed_not_acceptance", "new_axioms": 0,
        "bindings": {"proof_sha256": digest(HERE/"proof.md"), "check_sha256": digest(Path(__file__)), "admission_sha256": digest(HERE.parent/"_admission"/"round1320_admission.md")},
        "scope": "Finite CAR parent and conditional fixed-mode local-current bridge only. No complete QFT instrument, full S1 action, actual unique gauge group, family selection or observed parameters certified.",
        "finite_parent": {"dimension": dim, "fermion_modes": 28, "all_two_fermion_states": 378, "scalar_states": 3, "vacuum_states": 1, "dark_two_fermion_dimension": 372},
        "capacity_bounds": capacities, "ratio_cases": ratio_cases,
        "balanced_point": {"eq": e_q, "eqbar": e_qb, "ellbar": e_lb, "ES": es, "M": mass, "gA": [ga.real, ga.imag], "gC": [gc.real, gc.imag], "peak_time": float(tstar), "kappa": [[z.real, z.imag] for z in kappa], "lambda": [[z.real, z.imag] for z in lam]},
        "conditional_projection_example": {"eta_assumed_not_proven": eta, "two_pass_bound": float(projection_error_budget)},
        "family_scope": "Only one declared bright flavor was constructed. Candidate family counts 3/4/5 were not numerically certified as full joint theories in this round.",
        "checks": CHECKS, "check_count": len(CHECKS), "maximum_error": max(c["error"] for c in CHECKS),
    }


def equivalent(a, b, path="root"):
    if isinstance(a, dict):
        if not isinstance(b, dict) or a.keys() != b.keys(): raise AssertionError(path+" keys")
        for key in a: equivalent(a[key], b[key], path+"."+key)
    elif isinstance(a, (list, tuple)):
        if not isinstance(b, list) or len(a) != len(b): raise AssertionError(path+" length")
        for i, (x, y) in enumerate(zip(a, b)): equivalent(x, y, f"{path}[{i}]")
    elif isinstance(a, (int, float)) and not isinstance(a, bool):
        if not np.isclose(a, b, rtol=1e-11, atol=1e-12): raise AssertionError((path, a, b))
    elif a != b: raise AssertionError((path, a, b))


def main():
    parser = argparse.ArgumentParser()
    flags = parser.add_mutually_exclusive_group()
    flags.add_argument("--check", action="store_true")
    flags.add_argument("--write-new", action="store_true")
    args = parser.parse_args()
    result = calculate()
    destination = HERE/"results.json"
    if args.write_new:
        with destination.open("x", encoding="utf-8", newline="\n") as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2, allow_nan=False); stream.write("\n")
        mode = "exclusive_new_result"
    else:
        equivalent(result, json.loads(destination.read_text(encoding="utf-8")))
        mode = "read_only_recalculation"
    print(json.dumps({"status": "PASS", "mode": mode, "checks": result["check_count"], "maximum_error": result["maximum_error"]}, ensure_ascii=False))


if __name__ == "__main__":
    main()
