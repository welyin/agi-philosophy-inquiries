"""1326 author: finite driven reference, full low band and fixed-source ledger.

The occupation vertices extend the frozen 1325 definition. No old author or
independent implementation is imported. The full finite boson cutoff is a
model input; this program does not claim an infinite-field or Higgs identity.
Default/--check are read-only. --write-new creates results.json exclusively.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
ROOT = STAGE.parent.parent
MATTER_C = (0, 0, 1, 1)  # D,L0,L1,U
MATTER_K = (0, 1, 1, 2)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def norm(a):
    a = np.asarray(a)
    if a.ndim == 0:
        return float(abs(a))
    if a.ndim == 1:
        return float(np.linalg.norm(a))
    return float(np.linalg.svd(a, compute_uv=False)[0]) if a.size else 0.


def same_result(a, b, path="result"):
    if isinstance(a, dict):
        assert isinstance(b, dict) and a.keys() == b.keys(), path
        for key in a:
            same_result(a[key], b[key], path + "." + key)
    elif isinstance(a, list):
        assert isinstance(b, list) and len(a) == len(b), path
        for i, (x, y) in enumerate(zip(a, b)):
            same_result(x, y, f"{path}[{i}]")
    elif isinstance(a, float):
        assert isinstance(b, (int, float)) and not isinstance(b, bool), path
        assert np.isfinite(b) and np.isclose(a, b, rtol=3e-9, atol=3e-10), path
    else:
        assert type(a) is type(b) and a == b, path


def matrix_data(a):
    return {"real": np.asarray(a).real.tolist(), "imag": np.asarray(a).imag.tolist()}


def reference(B, f):
    """All N1 sectors of the exact finite reference, not Poisson truncation."""
    sectors = []
    for n1 in range(B + 1):
        size = B - n1 + 1
        R = np.diag(np.arange(size, dtype=float) + n1)
        lower = np.sqrt(np.arange(1, size, dtype=float))
        S = -np.diag(lower, 1) - np.diag(lower, -1)
        R += f * S
        w, q = np.linalg.eigh(R)
        sectors.append((w, q, R, S))
    ordered = sorted((float(w), n1, j) for n1, sec in enumerate(sectors)
                     for j, w in enumerate(sec[0]))
    e0, ground_sector, _ = ordered[0]
    assert ground_sector == 0
    chi = sectors[0][1][:, 0].copy()
    if chi[0] < 0:
        chi *= -1
    return dict(e0=e0, gap=ordered[1][0] - e0, chi=chi,
                source=float(chi @ sectors[0][3] @ chi),
                number=float(np.abs(chi)**2 @ np.arange(B + 1)),
                sector_bottoms=[float(sec[0][0]) for sec in sectors],
                first_reference_levels=ordered[:6])


class Parent:
    def __init__(self, B, p, gd, gu, normalized, axis=True):
        self.B, self.p = B, p
        self.gd, self.gu, self.normalized = gd, gu, normalized
        self.basis = [(m, n0, N - n0) for m in range(4)
                      for N in range(B + 1) for n0 in range(N + 1)
                      if not axis or N - n0 + MATTER_C[m] in (0, 1)]
        self.pos = {s: i for i, s in enumerate(self.basis)}
        d = len(self.basis)
        self.V = np.zeros((d, d), complex)
        self.S = np.zeros((d, d), complex)
        self.N = np.diag([n0 + n1 for m, n0, n1 in self.basis]).astype(complex)
        self.D = np.diag([p * (1 if m in (0, 3) else -1)
                          for m, n0, n1 in self.basis]).astype(complex)
        self.C = np.array([n1 + MATTER_C[m] for m, n0, n1 in self.basis], float)
        self.K = np.array([n0 + n1 + MATTER_K[m] for m, n0, n1 in self.basis], float)
        for j, (m, n0, n1) in enumerate(self.basis):
            if n0:
                i = self.pos[(m, n0 - 1, n1)]
                self.S[i, j] = self.S[j, i] = -np.sqrt(n0)
            N = n0 + n1
            if m == 0 and N:
                for alpha, occ in enumerate((n0, n1)):
                    if occ:
                        state = (alpha + 1, n0 - int(alpha == 0), n1 - int(alpha == 1))
                        i = self.pos[state]
                        value = gd * np.sqrt(occ / N if normalized else occ)
                        self.V[i, j] += value
                        self.V[j, i] += value.conjugate()
            elif m == 3 and N < B:
                for alpha, beta, sign in ((0, 1, 1), (1, 0, -1)):
                    state = (alpha + 1, n0 + int(beta == 0), n1 + int(beta == 1))
                    i = self.pos[state]
                    occ = (n0, n1)[beta]
                    value = sign * gu * np.sqrt((occ + 1)/(N + 2) if normalized else occ + 1)
                    self.V[i, j] += value
                    self.V[j, i] += value.conjugate()
        self.global_g = max(abs(gd), abs(gu)) * (1. if normalized else np.sqrt(B + 2))
        # The meter is a declared finite two-outcome source instrument.
        self.source_scale = 2 * np.sqrt(B)
        w, u = np.linalg.eigh(self.S)
        self.meters = []
        for sign in (1, -1):
            effects = (1 + sign*w/self.source_scale)/2
            assert min(effects) > -1e-12
            self.meters.append((u * np.sqrt(np.maximum(0, effects))) @ u.conj().T)

    def embed(self, ref):
        E = np.zeros((len(self.basis), 4), complex)
        for j, (m, n0, n1) in enumerate(self.basis):
            if n1 == 0:
                E[j, m] = ref["chi"][n0]
        return E

    def solve(self, f):
        ref = reference(self.B, f)
        E = self.embed(ref)
        HR = self.N + f*self.S
        H = HR + self.D + self.V
        eig, U = np.linalg.eigh(H)
        low = U[:, :4]
        projection_E = low @ (low.conj().T @ E)
        c, q = np.linalg.eigh(E.conj().T @ projection_E)
        assert min(c) > 0
        J = projection_E @ (q @ np.diag(1/np.sqrt(c)) @ q.conj().T)
        hd = J.conj().T @ H @ J
        A = E.conj().T @ H @ E
        off = self.V @ E - E @ (E.conj().T @ self.V @ E)
        beta = norm(off)
        delta = ref["gap"] - 2*(abs(self.p) + self.global_g)
        assert delta > beta > 0
        return dict(ref=ref, E=E, HR=HR, H=H, eig=eig, U=U, low=low,
                    J=J, hd=hd, A=A, beta=beta, delta=delta,
                    energy_error=beta*beta/delta)


def compute():
    checks = []

    def close(name, a, b=0., tolerance=2e-9):
        residual = norm(np.asarray(a) - np.asarray(b))
        assert residual <= tolerance, (name, residual, tolerance)
        checks.append(dict(name=name, kind="equality", residual=residual, tolerance=tolerance))

    def bound(name, value, ceiling, tolerance=2e-10):
        value, ceiling = float(value), float(ceiling)
        assert value <= ceiling + tolerance, (name, value, ceiling)
        checks.append(dict(name=name, kind="upper_bound", value=value, ceiling=ceiling, tolerance=tolerance))

    B, p = 66, .005
    normalized = Parent(B, p, .013, .011, True)
    ref1 = reference(B, 1.)
    E1 = normalized.embed(ref1)
    bare_unit = Parent(B, p, 1., 1., False)
    mean_n = E1.conj().T @ normalized.V @ E1
    mean_b = E1.conj().T @ bare_unit.V @ E1
    gd = float((mean_n[1, 0]/mean_b[1, 0]).real)
    gu = float((mean_n[2, 3]/mean_b[2, 3]).real)
    bare = Parent(B, p, gd, gu, False)
    close("calibration.all_first_low_matrix", E1.conj().T @ bare.V @ E1, mean_n)
    close("axis.dimension", len(normalized.basis), 6*B + 4)
    cases = []
    stage_solutions = {}
    for label, model in (("normalized", normalized), ("bare_linear", bare)):
        prefix = label + "."
        close(prefix+"V_Hermitian", model.V, model.V.conj().T)
        close(prefix+"C_V_invariant", (model.C[:, None]-model.C[None, :])*model.V)
        close(prefix+"C_drive_invariant", (model.C[:, None]-model.C[None, :])*model.S)
        close(prefix+"K_V_invariant", (model.K[:, None]-model.K[None, :])*model.V)
        # The old K conservation is deliberately not claimed with drive on.
        bound(prefix+"drive_breaks_old_K", 1., norm((model.K[:, None]-model.K[None, :])*model.S))
        bound(prefix+"V_compression_bound", norm(model.V), model.global_g)
        bound(prefix+"source_bound", norm(model.S), model.source_scale)
        close(prefix+"source_instrument_complete",
              sum(M.conj().T @ M for M in model.meters), np.eye(len(model.basis)))
        close(prefix+"source_instrument_effect_difference",
              model.source_scale*(model.meters[0]@model.meters[0]-model.meters[1]@model.meters[1]), model.S)
        for f in (1., 2.):
            x = model.solve(f)
            stage_solutions[label, f] = x
            ref, E, H, low, J, hd = (x[k] for k in ("ref", "E", "H", "low", "J", "hd"))
            tag = prefix + "f" + str(int(f)) + "."
            a = np.linalg.eigvalsh(x["A"])
            exact = x["eig"][:4]
            delta, beta, error = x["delta"], x["beta"], x["energy_error"]
            effective = E.conj().T@model.V@E
            first_order_magnitudes = [float(abs(effective[1, 0])), float(abs(effective[2, 3]))]
            energies = np.sqrt(p*p+np.array(first_order_magnitudes)**2)
            close(tag+"inherited_first_order_Dirac_dictionary", a-ref["e0"],
                  np.sort(np.concatenate((-energies, energies))))
            for j in range(4):
                bound(tag+f"ordered_upper_{j}", exact[j], a[j])
                bound(tag+f"ordered_lower_{j}", a[j]-error, exact[j])
            bound(tag+"next_band_gap", delta, x["eig"][4]-exact[3])
            bound(tag+"global_other_C_threshold", exact[3], ref["e0"]+ref["gap"]-abs(p)-model.global_g)
            close(tag+"product_embedding", E.conj().T@E, np.eye(4))
            close(tag+"dressed_isometry", J.conj().T@J, np.eye(4))
            close(tag+"exact_dressed_intertwining", H@J, J@hd)
            y = low-E@(E.conj().T@low)
            rho = beta/delta
            jbound = np.sqrt(2-2*np.sqrt(1-rho*rho))
            bound(tag+"whole_subspace_angle", norm(y), rho)
            bound(tag+"canonical_dressing", norm(J-E), jbound)
            time = 20.
            evolved = x["U"] @ (np.exp(-1j*time*x["eig"])[:, None]*(x["U"].conj().T@E))
            eh, qh = np.linalg.eigh(hd)
            ue = (qh*np.exp(-1j*time*eh))@qh.conj().T
            bound(tag+"all_input_product_to_dressed_dynamics", norm(evolved-E@ue), 2*jbound)
            # Derivatives are fixed-parent diagnostics, not differentiated
            # error bars. All four energies and the entire source are checked.
            df = 1e-4
            plus, minus = model.solve(f+df), model.solve(f-df)
            hprime = (plus["hd"]-minus["hd"])/(2*df)
            jprime = (plus["J"]-minus["J"])/(2*df)
            connection = J.conj().T@jprime
            physical_source = J.conj().T@model.S@J
            close(tag+"dressed_source_connection", hprime-(hd@connection-connection@hd), physical_source, 3e-8)
            hf = np.diag(low.conj().T@model.S@low).real-ref["source"]
            fd = ((plus["eig"][:4]-plus["ref"]["e0"])
                  -(minus["eig"][:4]-minus["ref"]["e0"]))/ (2*df)
            close(tag+"all_four_HF_finite_difference", fd, hf, 3e-8)
            close(tag+"reference_HF_finite_difference",
                  (plus["ref"]["e0"]-minus["ref"]["e0"])/(2*df), ref["source"], 3e-8)
            number = np.diag(low.conj().T@model.N@low).real
            reference_energy = np.diag(low.conj().T@x["HR"]@low).real-ref["e0"]
            kinetic = np.diag(low.conj().T@model.D@low).real
            vertex = np.diag(low.conj().T@model.V@low).real
            close(tag+"all_four_full_energy_ledger", exact-ref["e0"], reference_energy+kinetic+vertex)
            bound(tag+"reference_excess_nonnegative", 0., min(reference_energy))
            # A real finite source instrument, including both outcomes,
            # disturbance and a nontrivial unknown partner test. The pointer
            # record energies are declared degenerate; no reset is performed.
            KJ = [M@J for M in model.meters]
            close(tag+"record_stinespring_isometry", sum(Z.conj().T@Z for Z in KJ), np.eye(4))
            close(tag+"source_readout_all_input",
                  model.source_scale*(KJ[0].conj().T@KJ[0]-KJ[1].conj().T@KJ[1]), physical_source)
            coeff = np.array([[1.,1j,0.], [.2j,.3,.4], [.1,-.2j,.6], [.5,.2,.1j]], complex)
            coeff /= np.linalg.norm(coeff)
            before = J@coeff
            after = [Z@coeff for Z in KJ]
            close(tag+"record_partner_marginal", sum(Z.conj().T@Z for Z in after), before.conj().T@before)
            read_energy = sum(Z.conj().T@H@Z for Z in KJ)-hd
            component_changes = []
            for component in (model.N, f*model.S, model.D, model.V):
                component_changes.append(sum(Z.conj().T@component@Z for Z in KJ)-J.conj().T@component@J)
            close(tag+"all_input_readout_work_ledger", read_energy, sum(component_changes))
            probabilities = [float(np.vdot(Z, Z).real) for Z in after]
            work = float(np.trace(coeff.conj().T@read_energy@coeff).real)
            source_read = model.source_scale*(probabilities[0]-probabilities[1])
            close(tag+"actual_record_source", source_read, np.vdot(before, model.S@before).real)
            # The exact instantaneous drive-switch work is a parent identity,
            # not a claim of an adiabatically implemented state preparation.
            df_switch = .03
            close(tag+"fixed_state_drive_switch_work",
                  np.vdot(before, (df_switch*model.S)@before).real, df_switch*source_read)
            # The classical effective drive also supplies the old number/charge
            # current. It is not an autonomous, already-quantized reservoir.
            K = np.diag(model.K)
            old_charge = .23*np.eye(len(model.basis))+.37*K
            K_rate = 1j*(H@K-K@H)
            drive_K_rate = 1j*f*(model.S@K-K@model.S)
            close(tag+"complete_K_current_is_drive_current", K_rate, drive_K_rate)
            close(tag+"complete_old_charge_current", 1j*(H@old_charge-old_charge@H), .37*K_rate)
            dt = 1e-4
            transformed = x["U"].conj().T@before
            future = x["U"]@(np.exp(-1j*dt*x["eig"])[:,None]*transformed)
            past = x["U"]@(np.exp(1j*dt*x["eig"])[:,None]*transformed)
            finite_K_rate = (np.vdot(future, K@future)-np.vdot(past, K@past)).real/(2*dt)
            current = float(np.vdot(before, K_rate@before).real)
            close(tag+"actual_finite_K_Ehrenfest", finite_K_rate, current, 3e-8)
            read_K_change = float((sum(np.vdot(Z, K@Z) for Z in after)-np.vdot(before,K@before)).real)
            cases.append(dict(menu=label, f=f, reference_ground=ref["e0"],
                reference_global_gap=ref["gap"], reference_sector_bottoms=ref["sector_bottoms"],
                reference_mean_number=ref["number"],
                global_vertex_norm_bound=model.global_g, actual_PQ_vertex_norm=beta,
                separation_denominator=delta, ordered_energy_error=error,
                actual_next_band_gap=float(x["eig"][4]-exact[3]),
                low_excitation_energies=(exact-ref["e0"]).tolist(),
                first_order_low_excitation_energies=(a-ref["e0"]).tolist(),
                first_order_conversion_magnitudes=first_order_magnitudes,
                physical_source_excitation=hf.tolist(), source_finite_difference=fd.tolist(),
                reference_excess_energy=reference_energy.tolist(),
                reference_mean_number_in_each_low_state=number.tolist(),
                kinetic_energy=kinetic.tolist(), vertex_energy=vertex.tolist(),
                dressed_source=matrix_data(physical_source),
                dressing_error=norm(J-E), dressing_bound=float(jbound),
                all_input_dynamics_time=time, all_input_dynamics_error=norm(evolved-E@ue),
                instrument=dict(scale=float(model.source_scale), probabilities=probabilities,
                    source_mean=source_read, mean_work=work,
                    component_work=[float(np.trace(coeff.conj().T@C@coeff).real) for C in component_changes],
                    record_energy=0., resets_performed=0, old_K_change=read_K_change,
                    old_charge_change=.37*read_K_change),
                drive_account=dict(old_K_current=current, old_charge_current=.37*current,
                                   switch_step=df_switch, instantaneous_switch_work=df_switch*source_read)))

    # One complete finite 760-state kernel checks the exact C0/1 restriction.
    smallB, smallf = 18, 1.25
    small_full = Parent(smallB, p, .013, .011, True, axis=False)
    small_axis = Parent(smallB, p, .013, .011, True, axis=True)
    Q = np.zeros((len(small_full.basis), len(small_axis.basis)), complex)
    for j, state in enumerate(small_axis.basis):
        Q[small_full.pos[state], j] = 1
    fullH = small_full.N+smallf*small_full.S+small_full.D+small_full.V
    axisH = small_axis.N+smallf*small_axis.S+small_axis.D+small_axis.V
    close("full760.axis_isometry", Q.conj().T@Q, np.eye(len(small_axis.basis)))
    close("full760.exact_invariant_subspace", fullH@Q, Q@axisH)
    full_eig = np.linalg.eigvalsh(fullH)
    axis_eig = np.linalg.eigvalsh(axisH)
    close("full760.all_four_low_energies", full_eig[:4], axis_eig[:4])
    close("full760.entire_C_selection", (small_full.C[:,None]-small_full.C[None,:])*fullH)
    bound("full760.global_normalized_norm", norm(small_full.V), .013)

    contrasts = []
    for label in ("normalized", "bare_linear"):
        one, two = stage_solutions[label, 1.], stage_solutions[label, 2.]
        first = np.linalg.eigvalsh(one["A"])-one["ref"]["e0"]
        second = np.linalg.eigvalsh(two["A"])-two["ref"]["e0"]
        estimate = second-first
        interval = np.column_stack((estimate-two["energy_error"], estimate+one["energy_error"]))
        actual = ((two["eig"][:4]-two["ref"]["e0"])
                  -(one["eig"][:4]-one["ref"]["e0"]))
        for j in range(4):
            bound(label+f".finite_response_interval_lower_{j}", interval[j,0], actual[j])
            bound(label+f".finite_response_interval_upper_{j}", actual[j], interval[j,1])
        contrasts.append(dict(menu=label, first_order_difference=estimate.tolist(),
                              exact_difference=actual.tolist(), intervals=interval.tolist()))
    separation = contrasts[1]["intervals"][3][0]-contrasts[0]["intervals"][3][1]
    bound("held_out_highest_low_level_intervals_disjoint", 0., separation)
    return dict(round=1326, status="PASS", count=len(checks), checks=checks,
        max_equality_residual=max(c["residual"] for c in checks if c["kind"]=="equality"),
        model=dict(boson_cutoff=B, full_dimension=4*(B+1)*(B+2)//2,
                   exact_axis_dimension=len(normalized.basis), p=p, omega=1.,
                   normalized_couplings=[.013,.011], once_calibrated_bare_couplings=[gd,gu],
                   calibration="E(f=1)^dagger V E(f=1); all parent coefficients fixed thereafter",
                   source="dH/df = -(a0+a0^dagger); explicit finite drive, not Higgs identity"),
        cases=cases, finite_response=contrasts, held_out_interval_gap=float(separation),
        complete_small_kernel=dict(B=smallB, f=smallf, full_dimension=len(small_full.basis),
                                   axis_dimension=len(small_axis.basis)),
        limitations=[
            "exact finite boson cutoff is a model input, not a demonstrated physical minimum scale",
            "one inherited helicity momentum kernel; no full Higgs, gauge field or QFT identity",
            "driven low-band setting and old resonant handoff are different settings",
            "dressed isometry is a dictionary, not a free preparation controller",
            "source instrument, records, reference preparation and control are explicitly declared finite permissions",
            "bare menu is a balance-contract deletion comparison, not a second model satisfying all old balance conditions",
            "finite response intervals do not arise by differentiating an uncontrolled energy error",
            "no new axioms adopted and no realistic mass, family count or UV prediction claimed",
        ])


def run():
    result = compute()
    paths = [HERE/"proof.md", HERE/"check.py", STAGE/"_admission/round1326_admission.md",
             STAGE/"1325/proof.md", STAGE/"1325/check.py", STAGE/"1325/results.json",
             STAGE.parent/"archive_1313_1322/1316/proof.md",
             STAGE.parent/"archive_1313_1322/1322/proof.md"]
    result["source_sha256"] = {p.relative_to(ROOT).as_posix(): digest(p) for p in paths}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--check", action="store_true")
    mode.add_argument("--write-new", action="store_true")
    args = parser.parse_args()
    result = run()
    target = HERE/"results.json"
    if args.write_new:
        with target.open("x", encoding="utf-8") as output:
            json.dump(result, output, ensure_ascii=False, indent=2)
            output.write("\n")
    else:
        same_result(result, json.loads(target.read_text(encoding="utf-8")))
    print(json.dumps(dict(round=1326, mode="write-new" if args.write_new else "read-only",
        passed=True, count=result["count"], max_equality_residual=result["max_equality_residual"],
        held_out_interval_gap=result["held_out_interval_gap"])))


if __name__ == "__main__":
    main()
