"""1274: inherited autonomous carrier, full trajectory and boundary entropy role.

One author group. Default is read-only replay. --write exclusively creates
results.json. Imports only the frozen 1237 carrier, with bytecode disabled.
Finite samples are diagnostics; continuous identities and strict bounds are
proved in proof.md. No new material, readout, actual temperature or horizon.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True
import numpy as np

HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parents[1]
OLD = HERE.parent / "1237" / "check.py"
TOL = 3e-10


def need(ok, msg):
    if not bool(ok):
        raise AssertionError(msg)


def norm(a):
    return float(np.linalg.norm(a))


def ent(vals):
    v = np.asarray(vals, dtype=float)
    need(v.min() > -TOL, "positive spectrum")
    v = v[v > 1e-13]
    return float(-np.sum(v * np.log(v)))


def rel_diag(values, references):
    total = 0.
    for v, r in zip(values, references):
        if v > 1e-13:
            need(r > 0, "relative entropy support")
            total += float(v * np.log(v / r))
    return total


def frac(x):
    return dict(numerator=x.numerator, denominator=x.denominator)


def inherited():
    spec = importlib.util.spec_from_file_location("round1237_readonly", OLD)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.carrier()


def strict_arithmetic():
    # ln 2 = 2 atanh(1/3) > 2/3; ratio > 4 proves D_* > 1/15.
    numerator, denominator = 2**33, 5**5 * 3**12
    need(numerator > 4 * denominator, "integer entropy gap")
    lower = F(1, 15)
    delta = F(1, 1000)
    entropy_each = delta * (2 + 7 + 1)
    entropy_remaining = lower - 2 * entropy_each
    amplitude_remaining = F(1, 2) - 4 * delta
    # Positive Taylor terms prove exp(2)>6 and exp(7)>1000.
    def exp_partial(x, n):
        term = F(1)
        value = term
        for j in range(1, n + 1):
            term *= F(x, j)
            value += term
        return value
    e2 = exp_partial(2, 3)
    e7 = exp_partial(7, 12)
    need(e2 > 6 and e7 > 1000, "log upper bounds")
    need(entropy_remaining == F(7, 150), "robust entropy")
    need(amplitude_remaining == F(62, 125), "robust area and energy")
    return dict(ratio=frac(F(numerator, denominator)),
                integer_gap_over_four=numerator - 4 * denominator,
                strict_entropy_drop_lower=frac(lower),
                state_half_trace_error_each=frac(delta),
                entropy_error_each_strict_upper=frac(entropy_each),
                robust_entropy_drop_strict_lower=frac(entropy_remaining),
                robust_area_and_energy_increase_lower=frac(amplitude_remaining),
                exp2_partial=frac(e2), exp7_partial=frac(e7),
                budget_is_mathematical_not_instrument=True)


def calculate():
    carrier_check, w, perms, tau = inherited()
    hb = np.array([4., 4., 6., 6., 6., 6., 6.])
    hw = np.array([0., 2.])
    bare = (hb[:, None] + hw).ravel()
    source, target = np.zeros((14, 2), complex), np.zeros((14, 2), complex)
    source[1, 0], source[3, 1] = 1., 1.
    target[4::2] = w
    u1 = (np.eye(14) - source @ source.conj().T - target @ target.conj().T
          + target @ source.conj().T + source @ target.conj().T)
    u2 = np.zeros((720, 7, 7), complex)
    u2[:, :2, :2] = np.eye(2)
    u2[:, 2:, 2:] = perms

    def gate1(x):
        return np.einsum("ij,jfk->ifk", u1, x.reshape(14, 720, -1),
                         optimize=True).reshape(x.shape)

    def gate2(x, adjoint=False):
        u = u2.conj().transpose(0, 2, 1) if adjoint else u2
        return np.einsum("fij,jwfk->iwfk", u, x, optimize=True)

    def interaction(x):
        out = np.zeros_like(x)
        out[0] = gate1(x[1])
        out[1] = gate1(x[0]) + gate2(x[2], True)
        out[2] = gate2(x[1])
        return out

    initial = np.zeros((7, 2, 720, 2), complex)
    initial[0, 1, :, 0] = 1 / np.sqrt(720)
    initial[1, 1, :, 1] = 1 / np.sqrt(720)
    states = [initial, gate1(initial)]
    states.append(gate2(states[1]))
    path = np.array([[0., 1., 0.], [1., 0., 1.], [0., 1., 0.]])
    ev, vec = np.linalg.eigh(path)
    tstar = np.pi / np.sqrt(2)
    beta_star = np.log(2.5) / 2
    Ldiag = np.log(np.array([2., 2., 5., 5., 5., 5., 5.]))
    reference_star = np.diag([.25, .25, .1, .1, .1, .1, .1])
    expected_star_spectrum = [.25, .25, .1, .1, .1, .1, .1]
    residuals = dict(inherited_operator_max=max(carrier_check["residuals"].values()))
    samples = []
    for fraction in (0., .25, .5, .75, 1., 1.5, 2.):
        time = fraction * tstar
        amplitudes = vec @ (np.exp(-1j * ev * time) * vec[0, :].conj())
        psi = np.stack([amplitudes[j] * states[j] for j in range(3)]) * np.exp(-6j*time)
        flat = psi.reshape(-1, 2)
        code_gram = flat.conj().T @ flat
        rho_b = np.einsum("cbwfr,cdwfr->bd", psi, psi.conj(),
                          optimize=True) / 2
        c = float(np.cos(np.pi * fraction))
        p0, p1, p2 = (1+c)**2/4, (1-c*c)/2, (1-c)**2/4
        expected = np.zeros((7, 7), complex)
        expected[:2, :2] = p0 * np.eye(2) / 2
        expected[2:, 2:] = p1 * tau + p2 * np.eye(5) / 5
        analytic_eigenvalues = np.array([p0/2]*2 + [p1/2+p2/5]*2 + [p2/5]*3)
        actual_eigenvalues = np.linalg.eigvalsh(rho_b)
        label = str(fraction)
        residuals["code_isometry_"+label] = norm(code_gram-np.eye(2))
        residuals["partial_trace_"+label] = norm(rho_b-expected)
        residuals["spectrum_"+label] = norm(np.sort(analytic_eigenvalues)-actual_eigenvalues)
        for i in range(2):
            for j in range(2):
                reduced = np.einsum("cbwf,cdwf->bd", psi[..., i], psi[..., j].conj(),
                                    optimize=True)
                unit = np.zeros((2, 2), complex)
                unit[i, j] = 1
                channel = np.zeros((7, 7), complex)
                channel[:2, :2] = p0 * unit
                channel[2:, 2:] = p1 * w @ unit @ w.conj().T
                if i == j:
                    channel[2:, 2:] += p2 * np.eye(5) / 5
                residuals["matrix_unit_%s_%d%d"%(label, i, j)] = norm(reduced-channel)

        sb = ent(actual_eigenvalues)
        area = float(np.trace(rho_b @ np.diag(hb)).real)
        lmean = float(np.trace(rho_b @ np.diag(Ldiag)).real)
        sector_mix = ent([p0, 1-p0])
        moving_reference_spectrum = [p0/2]*2+[(1-p0)/5]*5
        d_moving = rel_diag(analytic_eigenvalues, moving_reference_spectrum)
        d_fixed = rel_diag(analytic_eigenvalues, expected_star_spectrum)
        residuals["geometry_"+label] = abs(area-(6-2*p0))
        residuals["structural_affine_"+label] = abs(lmean-np.log(2)-beta_star*(area-4))
        residuals["structural_split_"+label] = abs(sb-sector_mix-lmean+d_moving)
        residuals["fixed_modular_split_"+label] = abs(d_fixed-lmean-np.log(2)+sb)

        interaction_state = interaction(psi)
        free_state = bare.reshape(1, 7, 2, 1, 1) * psi
        total_energy = float(np.vdot(psi, free_state+interaction_state).real/2)
        work_energy = float(np.vdot(psi, hw.reshape(1, 1, 2, 1, 1)*psi).real/2)
        coupling_mean = float(np.vdot(psi, interaction_state).real/2)
        coupling_second = float(np.vdot(interaction_state, interaction_state).real/2)
        residuals["total_energy_"+label] = abs(total_energy-6)
        residuals["work_energy_"+label] = abs(work_energy-2*p0)
        residuals["interaction_mean_"+label] = abs(coupling_mean)
        residuals["interaction_second_"+label] = abs(coupling_second-1)

        # Purification R of I2/2: E=CWF has the same nonzero spectrum as BR.
        br_columns = psi.transpose(1, 4, 0, 2, 3).reshape(14, -1)/np.sqrt(2)
        rho_br = br_columns @ br_columns.conj().T
        environment_spectrum = np.linalg.eigvalsh(rho_br)
        se = ent(environment_spectrum)
        s_be = ent(np.linalg.eigvalsh(code_gram/2))
        mutual = sb+se-s_be
        residuals["full_entropy_"+label] = abs(s_be-np.log(2))
        residuals["full_resource_balance_"+label] = abs(sb-np.log(2)-mutual+se)
        samples.append(dict(t_over_tstar=fraction, cosine=c,
                            program_probabilities=[p0,p1,p2],
                            actual_spectrum=actual_eigenvalues.tolist(),
                            analytic_spectrum=analytic_eigenvalues.tolist(),
                            boundary_entropy=sb,area=area,boundary_energy=area,
                            work_energy=work_energy,total_energy=total_energy,
                            interaction_mean=coupling_mean,
                            interaction_second_moment=coupling_second,
                            structural_mean=lmean,sector_mixing_entropy=sector_mix,
                            relative_entropy_to_moving_reference=d_moving,
                            relative_entropy_to_full_Gibbs_marginal=d_fixed,
                            environment_entropy=se,BE_entropy=s_be,
                            B_environment_mutual_information=mutual))

    # Every flag is conserved by the original Hamiltonian. Construct all 720
    # actual 42-dimensional (C,B,W) blocks and exponentiate the full generator.
    u2bw = np.asarray([np.kron(u, np.eye(2)) for u in u2])
    full_h = np.zeros((720, 42, 42), complex)
    for cindex in range(3):
        sl = slice(14*cindex, 14*(cindex+1))
        full_h[:, sl, sl] = np.diag(bare)
    full_h[:, :14, 14:28] = u1
    full_h[:, 14:28, :14] = u1
    full_h[:, 28:42, 14:28] = u2bw
    full_h[:, 14:28, 28:42] = u2bw.conj().transpose(0,2,1)
    energies, vectors = np.linalg.eigh(full_h)
    expected_energies = np.sort((bare[:,None]+ev[None,:]).ravel())
    residuals["full_Gibbs_block_spectra"] = float(np.max(np.abs(energies-expected_energies)))
    gibbs_checks = []
    for beta in (.37, beta_star):
        weights = np.exp(-beta*(energies-energies.min()))
        weights /= weights.sum()
        blocks = (vectors*weights[:,None,:]) @ vectors.conj().transpose(0,2,1)
        tensor = blocks.reshape(720,3,7,2,3,7,2)
        marginal = np.einsum("fcbwcdw->bd", tensor, optimize=True)
        bare_reference = np.exp(-beta*(hb-hb.min()))
        bare_reference /= bare_reference.sum()
        error = norm(marginal-np.diag(bare_reference))
        residuals["full_Gibbs_marginal_"+str(beta)] = error
        gibbs_checks.append(dict(beta=float(beta),
                                boundary_diagonal=marginal.diagonal().real.tolist(),
                                full_generator_marginal_error=error))
        if beta == beta_star:
            residuals["selected_reference"] = norm(marginal-reference_star)
            residuals["selected_modular_identity"] = norm(
                -np.log(marginal.diagonal().real)-Ldiag-np.log(2))

    mid = next(s for s in samples if s["t_over_tstar"]==.5)
    end = next(s for s in samples if s["t_over_tstar"]==1.)
    drop = mid["boundary_entropy"]-end["boundary_entropy"]
    exact_log_drop = (33*np.log(2)-5*np.log(5)-12*np.log(3))/20
    delta_s = -drop
    delta_h = end["boundary_energy"]-mid["boundary_energy"]
    delta_l = end["structural_mean"]-mid["structural_mean"]
    delta_d = (end["relative_entropy_to_full_Gibbs_marginal"]-
               mid["relative_entropy_to_full_Gibbs_marginal"])
    residuals["finite_drop"] = abs(drop-exact_log_drop)
    residuals["finite_modular_balance"] = abs(delta_s-beta_star*delta_h+delta_d)
    residuals["finite_structure_energy"] = abs(delta_l-beta_star*delta_h)
    need(drop>1/15, "diagnostic agrees with strict entropy gap")
    need(max(residuals.values())<TOL, "all inherited-operator trajectory checks")

    return dict(schema="research_round_1274_v1", status="PASS", author_groups=1,
                source_construction=dict(kinematic_dimension=144,
                    physical_dimension=7,old_input_dimension=2,
                    permutations=720,program_dimension=3,work_dimension=2,
                    direct_Gibbs_blocks=720,direct_Gibbs_block_dimension=42,
                    inherited_carrier_recounted_as_new_group=False),
                strict_certificate=strict_arithmetic(),
                samples=samples,full_Gibbs_checks=gibbs_checks,
                finite_interval=dict(t_start_over_tstar=.5,t_end_over_tstar=1.,
                    entropy_drop=drop,log_expression=exact_log_drop,
                    area_increase=delta_h,boundary_energy_increase=delta_h,
                    work_energy_change=end["work_energy"]-mid["work_energy"],
                    structural_entropy_increase=delta_l,
                    full_Gibbs_relative_entropy_increase=delta_d,
                    beta_star=float(beta_star),
                    entropy_time_derivative_at_mid_has_positive_sign=True),
                residuals=residuals,
                scope=dict(fixed_original_reference_trajectory=True,
                    all_code_channel_identity=True,
                    trajectory_entropy_claim_for_arbitrary_input=False,
                    independent_Casimir_area=True,
                    full_interacting_Gibbs_marginal_exact=True,
                    original_process_initial_state_is_Gibbs=False,
                    actual_reference_preparation_certified=False,
                    physical_Lorentz_area=False,actual_observer_temperature=False,
                    all_boundary_energy_is_heat=False,
                    user_full_conjecture_refuted=False,
                    Einstein_generation=False,new_devices=False))


def compare(current, saved, path="result"):
    if isinstance(current, dict):
        need(isinstance(saved,dict) and current.keys()==saved.keys(),path+" keys")
        for k in current:
            compare(current[k],saved[k],path+"."+k)
    elif isinstance(current,list):
        need(isinstance(saved,list) and len(current)==len(saved),path+" length")
        for j,(a,b) in enumerate(zip(current,saved)):
            compare(a,b,path+"["+str(j)+"]")
    elif isinstance(current,float):
        need(isinstance(saved,(int,float)) and np.isclose(current,saved,
             rtol=1e-10,atol=1e-10),path+" float")
    else:
        need(current==saved,path+" exact")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write",action="store_true",help="exclusive first result creation")
    args = parser.parse_args()
    result = calculate()
    bindings = [HERE/"admission.md",HERE/"proof.md",HERE/"check.py",
                OLD,OLD.with_name("proof.md"),OLD.with_name("results.json"),
                HERE.parent/"_admission/thermal_admission_independence_audit_after1238.md",
                HERE.parent/"_admission/thermal_mechanism_admission_after1270/geometric_thermal_identity_candidate_after1270.md"]
    result["source_sha256"] = {
        p.relative_to(RESEARCH).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in bindings}
    result_path=HERE/"results.json"
    if args.write:
        with result_path.open("x",encoding="utf-8",newline="\n") as stream:
            json.dump(result,stream,ensure_ascii=False,indent=2)
            stream.write("\n")
        mode="CREATED_EXCLUSIVE"
    else:
        compare(result,json.loads(result_path.read_text(encoding="utf-8")))
        mode="READ_ONLY_MATCH"
    print(json.dumps(dict(status="PASS",mode=mode,author_groups=1,
        entropy_drop=result["finite_interval"]["entropy_drop"],
        structural_increase=result["finite_interval"]["structural_entropy_increase"],
        beta_star=result["finite_interval"]["beta_star"],
        maximum_residual=max(result["residuals"].values())),ensure_ascii=False))


if __name__=="__main__":
    main()

