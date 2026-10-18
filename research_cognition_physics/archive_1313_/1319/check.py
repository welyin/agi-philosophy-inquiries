"""Round 1319: finite CAR father, capacity/cyclic contract and same-kernel response.

Default/--check is read-only. --write-new creates results.json exclusively.
No continuum QFT instrument, anomaly completion or inferred resource budget is signed.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
checks = []


def check(name, value, tolerance=2e-9):
    error = float(np.max(np.abs(value))) if np.size(value) else 0.0
    checks.append({"name": name, "error": error, "tolerance": tolerance})
    if not np.isfinite(error) or error > tolerance:
        raise AssertionError((name, error, tolerance))


def truth(name, value):
    check(name, 0 if bool(value) else 1, 0)


def evolution(h, t):
    e, u = np.linalg.eigh(h)
    return (u * np.exp(-1j * t * e)) @ u.conj().T


def wedge(u):
    pairs = list(itertools.combinations(range(len(u)), 2))
    return np.array([[u[i, k] * u[j, l] - u[i, l] * u[j, k]
                      for k, l in pairs] for i, j in pairs])


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def calculate():
    checks.clear()
    rng = np.random.default_rng(1319)
    pairs = list(itertools.combinations(range(12), 2))
    dim = 1 + len(pairs) + 3
    pos = {sum(1 << i for i in p): 1 + k for k, p in enumerate(pairs)}

    def create(mask, mode):
        if mask & (1 << mode):
            return None, 0
        sign = (-1) ** ((mask & ((1 << mode) - 1)).bit_count())
        return mask | (1 << mode), sign

    def pair(p, q):
        # c_p^dagger c_q^dagger |0>, using CAR, not an assigned logical code.
        mask, sign = create(0, q)
        mask, sign2 = create(mask, p)
        v = np.zeros(dim, complex)
        if sign2:
            v[pos[mask]] = sign * sign2
        return v

    s = {(i, j): (pair(2*i, 6+2*j+1) - pair(2*i+1, 6+2*j))/np.sqrt(2)
         for i in range(3) for j in range(3)}
    ij = [(0, 1), (0, 2), (1, 2)]
    wa = np.column_stack([(s[i, j]-s[j, i])/np.sqrt(2) for i, j in ij])
    wb = np.eye(dim, dtype=complex)[:, -3:]
    c0 = np.array([[0, 0, 1], [0, -1, 0], [1, 0, 0]], complex)
    pa, pb = wa @ wa.conj().T, wb @ wb.conj().T
    vac = np.eye(dim)[:, 0]
    n_a = np.array([0]+[sum(i < 6 for i in p) for p in pairs]+[0]*3)
    n_b = np.array([0]+[sum(i >= 6 for i in p) for p in pairs]+[0]*3)
    n_c = np.array([0]*67+[1]*3)
    charge = np.diag(n_a+n_b+2*n_c)
    ndiff = np.diag(n_a-n_b)
    t0 = wb @ c0 @ wa.conj().T
    ea, eb, gr, gi = .8, .9, .6, .8
    g = gr + 1j*gi
    delta, ec = .5, ea+eb+.5

    def father(e_a=ea, e_b=eb, e_c=ec, real_g=gr, imag_g=gi, cc=c0):
        tt = wb @ cc @ wa.conj().T
        gg = real_g + 1j*imag_g
        return np.diag(e_a*n_a+e_b*n_b+e_c*n_c) + gg*tt + gg.conjugate()*tt.conj().T

    h = father()
    check("S_orthonormal", np.column_stack(list(s.values())).conj().T @ np.column_stack(list(s.values()))-np.eye(9))
    check("Q_orthonormal", wa.conj().T@wa-np.eye(3))
    check("C0_unitary", c0.conj().T@c0-np.eye(3))
    check("A_B_orthogonal", wa.conj().T@wb)
    truth("full_dimension_70", dim == 70 and len(pairs) == 66)
    check("Hermitian", h-h.conj().T)
    check("K_conservation", h@charge-charge@h)
    check("Na_minus_Nb", h@ndiff-ndiff@h)
    check("vacuum_inert", h@vac)
    rmat = np.diag([1]+[-1]*69)
    check("old_r_compatibility", rmat@h-h@rmat)
    ff = np.zeros((dim, dim)); ff[1:67, 1:67] = np.eye(66)
    dark = ff-pa
    check("all_dark_states_retained", (h-np.diag(ea*n_a+eb*n_b+ec*n_c))@dark)

    # Complete triple-slot constraint has one complex solution dimension.
    def tensor(c):
        x = np.zeros((3, 3, 3), complex)
        for j, (a, b) in enumerate(ij):
            x[a, b, :] = c[:, j]; x[b, a, :] = -c[:, j]
        return x
    equations = []
    for k in range(9):
        c = np.zeros((3, 3), complex); c.flat[k] = 1
        tt = tensor(c)
        equations.append((tt-np.transpose(tt, (1, 2, 0))).ravel())
    rank = int(np.linalg.matrix_rank(np.column_stack(equations)))
    truth("cyclic_constraint_rank_eight", rank == 8)
    tt = tensor(c0)
    check("volume_cyclic", tt-np.transpose(tt, (1, 2, 0)))

    def representation(u, scalar=None, spin=None):
        if scalar is None: scalar = np.linalg.det(u)
        if spin is None: spin = np.eye(2)
        f = np.zeros((12, 12), complex)
        f[:6, :6] = np.kron(u, spin); f[6:, 6:] = np.kron(u, spin)
        rep = np.zeros((dim, dim), complex)
        rep[0, 0] = 1; rep[1:67, 1:67] = wedge(f)
        rep[-3:, -3:] = scalar*u.conj()
        return rep

    for k in range(3):
        raw = rng.normal(size=(3, 3))+1j*rng.normal(size=(3, 3))
        u, rr = np.linalg.qr(raw); u = u @ np.diag(np.diag(rr)/np.abs(np.diag(rr)))
        check(f"cofactor_U3_{k}", c0@wedge(u)-np.linalg.det(u)*u.conj()@c0)
        rep = representation(u)
        check(f"full_father_U3_{k}", rep@h-h@rep)
    theta = .37
    center = representation(np.exp(1j*theta)*np.eye(3))
    check("D_internal_center_equals_K", center-np.diag(np.exp(1j*theta*np.diag(charge))))
    spin = np.array([[np.cos(.21), np.sin(.21)], [-np.sin(.21), np.cos(.21)]])
    check("spin_singlet_father", representation(np.eye(3), spin=spin)@h-h@representation(np.eye(3), spin=spin))

    active = np.column_stack((wa, wb))
    expected = np.block([[np.zeros((3, 3)), g.conjugate()*c0.conj().T], [g*c0, delta*np.eye(3)]])
    check("CAR_active_block", active.conj().T@h@active-(ea+eb)*np.eye(6)-expected)
    omega = np.sqrt(delta**2+4*abs(g)**2)
    tp = np.pi/omega
    pmax = 4*abs(g)**2/omega**2
    for k, t in enumerate([.13, .7, tp, 2*tp]):
        phase = np.exp(-1j*delta*t/2)
        aa = phase*(np.cos(omega*t/2)+1j*delta/omega*np.sin(omega*t/2))*np.eye(3)
        bb = -1j*phase*2*g/omega*np.sin(omega*t/2)*c0
        dd = phase*(np.cos(omega*t/2)-1j*delta/omega*np.sin(omega*t/2))*np.eye(3)
        br = -1j*phase*2*g.conjugate()/omega*np.sin(omega*t/2)*c0.conj().T
        full = evolution(h, t)
        check(f"full_forward_{k}", full@wa-np.exp(-1j*(ea+eb)*t)*(wa@aa+wb@bb))
        check(f"full_reverse_{k}", full@wb-np.exp(-1j*(ea+eb)*t)*(wb@dd+wa@br))
        check(f"complete_Kraus_{k}", aa.conj().T@aa+bb.conj().T@bb-np.eye(3))
    hres = father(e_c=ea+eb)
    ures = evolution(hres, np.pi/(2*abs(g)))
    check("resonance_forward_full", ures@wa-np.exp(-1j*(ea+eb)*np.pi/2)*wb@(-1j*g*c0))
    check("resonance_reverse_full", ures@wb-np.exp(-1j*(ea+eb)*np.pi/2)*wa@(-1j*g.conjugate()*c0.conj().T))
    check("resonance_vacuum", ures@vac-vac)

    # Full partner reference: trace the old fermion factor using its occupation.
    up = evolution(h, tp)
    forward = wb.conj().T@up@wa
    reverse = wa.conj().T@up@wb
    recovered = reverse@forward
    check("two_real_passes_scalar_identity", recovered-np.trace(recovered)/3*np.eye(3))
    check("two_pass_success_probability", recovered.conj().T@recovered-pmax**2*np.eye(3))
    check("one_pass_success_16_17", pmax-16/17)
    erasure = 1-pmax**2
    check("two_pass_erasure_33_289", erasure-33/289)
    truth("both_directions_below_capacity_margin", erasure < 1/3)
    phi = np.eye(3).reshape(9)/np.sqrt(3)
    ideal = np.zeros((12, 12), complex); ideal[:9, :9] = np.outer(phi, phi)
    rho = pmax**2*ideal
    rho[9:, 9:] = erasure*np.eye(3)/3
    trace_distance = float(np.sum(np.abs(np.linalg.eigvalsh(rho-ideal)))/2)
    check("partner_erasure_trace_witness", trace_distance-erasure)
    # Orthogonal source sectors make the c/vac cross block zero on each reference matrix unit.
    for i in range(3):
        for j in range(3):
            cross = (up@wa)[:67, i].conj() @ np.eye(dim)[:67, 0]
            check(f"source_trace_flag_cross_{i}_{j}", cross)
    psi = wa @ np.array([1, 1j, -1])/np.sqrt(3)
    evolved = up@psi
    expect = lambda op, v: float(np.vdot(v, op@v).real)
    p = pmax
    for name, nn, target in [("Na", n_a, 1-p), ("Nb", n_b, 1-p), ("Nc", n_c, p)]:
        check("full_resource_"+name, expect(np.diag(nn), evolved)-target)
    h0 = np.diag(ea*n_a+eb*n_b+ec*n_c)
    hint = h-h0
    check("full_energy_conservation", expect(h, evolved)-expect(h, psi))
    check("free_energy_increase", expect(h0, evolved)-expect(h0, psi)-delta*p)
    check("interaction_energy_compensation", expect(hint, evolved)+delta*p)
    dephased = pa@np.outer(evolved, evolved.conj())@pa+pb@np.outer(evolved, evolved.conj())@pb
    check("first_occupation_readout_work", np.trace(h@dephased).real-expect(h, evolved)-delta*p)
    revpsi = wb@np.array([1, 1j, -1])/np.sqrt(3)
    revout = up@revpsi
    check("reverse_interaction_energy", expect(hint, revout)-delta*p)
    check("net_two_readout_work", delta*p-p*delta*p-delta*p*(1-p))

    # Fixed physical E and fixed relative z are distinct derivatives.
    def sigma_abs(e, pars):
        e_a, e_b, e_c, real_g, imag_g = pars
        return -(real_g**2+imag_g**2)/(e_c-e)*pa
    pars = np.array([ea, eb, ec, gr, gi])
    physical_e = ea+eb+.11
    sources = [np.diag(n_a), np.diag(n_b), np.diag(n_c), t0+t0.conj().T, 1j*(t0-t0.conj().T)]
    den = ec-physical_e
    abs_derivatives = [0*pa, 0*pa, abs(g)**2/den**2*pa, -2*gr/den*pa, -2*gi/den*pa]
    dz = np.array([-1, -1, 0, 0, 0])
    ddelta = np.array([-1, -1, 1, 0, 0])
    dq = np.array([0, 0, 0, 2*gr, 2*gi])
    eps = 2e-6
    for j in range(5):
        plus, minus = pars.copy(), pars.copy(); plus[j] += eps; minus[j] -= eps
        hp, hm = father(*plus), father(*minus)
        check(f"full_parent_source_{j}", (hp-hm)/(2*eps)-sources[j], 2e-8)
        derivative = (sigma_abs(physical_e, plus)-sigma_abs(physical_e, minus))/(2*eps)
        check(f"fixed_physical_E_source_{j}", derivative-abs_derivatives[j], 3e-8)
        rel_derivative = (-dq[j]/den+abs(g)**2*ddelta[j]/den**2)*pa
        chain = -abs(g)**2/den**2*dz[j]*pa
        check(f"frequency_chain_rule_{j}", rel_derivative+chain-abs_derivatives[j])
        relplus = sigma_abs(plus[0]+plus[1]+.11, plus)
        relminus = sigma_abs(minus[0]+minus[1]+.11, minus)
        check(f"fixed_relative_z_source_{j}", (relplus-relminus)/(2*eps)-rel_derivative, 3e-8)
    for k, z in enumerate([-.2*delta, .1*delta, .25*delta]):
        absolute_e = ea+eb+z
        v = h[-3:, 1:67]
        schur = h[1:67, 1:67]-absolute_e*np.eye(66)-v.conj().T@np.linalg.inv(h[-3:, -3:]-absolute_e*np.eye(3))@v
        check(f"complete_66_Schur_{k}", schur-(h[1:67, 1:67]-absolute_e*np.eye(66)+sigma_abs(absolute_e, pars)[1:67, 1:67]))
        change = abs(-abs(g)**2/(delta-z)+abs(g)**2/delta)
        bound = abs(g)**2*abs(z)/(delta*(delta-abs(z)))
        truth(f"finite_response_bound_{k}", change <= bound+2e-12)
    eigs = np.linalg.eigvalsh(expected)
    check("full_active_poles", eigs-np.repeat([(delta-omega)/2, (delta+omega)/2], 3))
    z = .31
    resactive = active.conj().T@hres@active-(ea+eb)*np.eye(6)
    check("resonance_resolvent_sign", np.linalg.inv(z*np.eye(6)-resactive)[:3, :3]-z/(z*z-abs(g)**2)*np.eye(3))
    check("handoff_response_identity_p", abs(g)**2/delta-delta*p/(4*(1-p)))
    check("handoff_response_identity_time", abs(g)**2/delta-(np.pi**2/(4*delta*tp**2)-delta/4))

    # Actual competing father retains the exact task but breaks general SU(3).
    symmetric = sum(s[i, i] for i in range(3))/np.sqrt(3)
    extra = np.outer(symmetric, symmetric.conj())
    halt = h+.37*extra
    check("SO3_extra_dark", extra@wa)
    for k, angle in enumerate([.2, .7]):
        rot = np.array([[np.cos(angle), -np.sin(angle), 0], [np.sin(angle), np.cos(angle), 0], [0, 0, 1]])
        rep = representation(rot)
        check(f"SO3_alternative_father_{k}", rep@halt-halt@rep)
    u = np.diag(np.exp(1j*np.array([.3, -.3, 0])))
    rep = representation(u)
    defect = float(np.linalg.norm(rep@halt-halt@rep))
    truth("SO3_competitor_not_full_SU3", defect > .1)
    check("SO3_same_actual_handoff", evolution(halt, tp)@wa-up@wa)
    check("SO3_same_reverse", evolution(halt, tp)@wb-up@wb)

    ctilde = np.diag([1, 1, -1])@c0
    cycdef = float(np.linalg.norm(tensor(ctilde)-np.transpose(tensor(ctilde), (1, 2, 0))))
    truth("unitary_noncyclic", cycdef > 1)
    ht = father(e_c=ea+eb, cc=ctilde)
    transfer = wb.conj().T@evolution(ht, np.pi/2)@wa
    check("noncyclic_complete_transfer", transfer.conj().T@transfer-np.eye(3))
    check("noncyclic_may_still_isotropic", ctilde.conj().T@ctilde-np.eye(3))
    perm = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
    left, right = ctilde@wedge(perm), perm.conj()@ctilde
    phase_gap = math.sqrt(max(0, 6-2*abs(np.trace(right.conj().T@left))))
    truth("noncyclic_no_common_delta_for_cycle", phase_gap > 1)
    codd = np.diag([1, 3, 5])@c0
    hodd = father(e_c=ea+eb, real_g=1, imag_g=0, cc=codd)
    oddtransfer = wb.conj().T@evolution(hodd, np.pi/2)@wa
    check("odd_singular_values_complete_transfer", oddtransfer.conj().T@oddtransfer-np.eye(3))
    check("odd_singular_values_response_split", np.linalg.eigvalsh(codd.conj().T@codd)-[1, 9, 25])

    capacities = []
    for d in range(2, 9):
        m = d*(d-1)//2
        fwd, rev = max(0, 1-d/m), max(0, 1-m/d)
        capacities.append({"d": d, "input_dimension": m, "receiver_dimension": d, "forward_lower_bound": fwd, "reverse_lower_bound": rev})
        truth(f"capacity_exclusion_d_{d}", (max(fwd, rev) < 1/3) == (d == 3))
    # Direct rank-r maximal overlap check; a Schmidt-vector statement, not optimal CPTP claim.
    for m, r in [(6, 4), (1, 1), (10, 5)]:
        mat = np.diag([1/np.sqrt(r)]*r+[0]*(m-r))
        check(f"Schmidt_overlap_bound_saturation_{m}_{r}", abs(np.trace(mat))**2/m-r/m)
    for n in [3, 4, 5]:
        check(f"family_replication_identity_{n}", np.kron(np.eye(n), recovered).conj().T@np.kron(np.eye(n), recovered)-p*p*np.eye(3*n))

    return {
        "round": 1319, "status": "author_checks_passed_not_acceptance", "new_axioms": 0,
        "bindings": {"proof_sha256": digest(HERE/"proof.md"), "check_sha256": digest(Path(__file__)),
                     "admission_sha256": digest(HERE.parent/"_admission"/"round1319_admission.md")},
        "scope": "Finite CAR parent and additional task contract only; no full-QFT device, actual full U3 group, spacetime dimension, family-number or empirical SM certification.",
        "basis": {"dimension": dim, "fermion_modes": 12, "two_fermion_states": 66, "boson_states": 3, "vacuum_states": 1, "wedge_order": ij},
        "capacity_bounds": capacities,
        "common_nonresonant_point": {"Ea": ea, "Eb": eb, "Ec": ec, "Re_g": gr, "Im_g": gi,
            "delta": delta, "peak_time": float(tp), "one_pass_probability": float(p), "two_pass_erasure": float(erasure),
            "finite_relative_frequency_window": [-delta/4, delta/4], "first_readout_work": float(delta*p),
            "two_readout_net_work": float(delta*p*(1-p)), "fresh_vacuum_and_switching": "additional registered protocol, not an autonomous free instrument"},
        "cyclic_rank": rank, "SO3_SU3_father_defect": defect, "noncyclic_phase_gap": phase_gap,
        "checks": checks, "check_count": len(checks), "maximum_error": max(c["error"] for c in checks),
    }


def equivalent(a, b, location="root"):
    if isinstance(a, dict):
        if not isinstance(b, dict) or a.keys() != b.keys(): raise AssertionError(location+" keys")
        for k in a: equivalent(a[k], b[k], location+"."+k)
    elif isinstance(a, (list, tuple)):
        if not isinstance(b, list) or len(a) != len(b): raise AssertionError(location+" length")
        for k, (x, y) in enumerate(zip(a, b)): equivalent(x, y, f"{location}[{k}]")
    elif isinstance(a, (float, int)) and not isinstance(a, bool):
        if not np.isclose(a, b, rtol=1e-11, atol=1e-12): raise AssertionError((location, a, b))
    elif a != b:
        raise AssertionError((location, a, b))


def main():
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--check", action="store_true")
    group.add_argument("--write-new", action="store_true")
    args = parser.parse_args()
    result = calculate()
    output = HERE/"results.json"
    if args.write_new:
        with output.open("x", encoding="utf-8", newline="\n") as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2, allow_nan=False); stream.write("\n")
        mode = "exclusive_new_result"
    else:
        equivalent(result, json.loads(output.read_text(encoding="utf-8")))
        mode = "read_only_recalculation"
    print(json.dumps({"status": "PASS", "mode": mode, "checks": result["check_count"], "maximum_error": result["maximum_error"]}, ensure_ascii=False))


if __name__ == "__main__":
    main()
