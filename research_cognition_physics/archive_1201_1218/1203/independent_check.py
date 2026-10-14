"""Independent 1203 audit. No imports of the author's 1203 implementation.
Default: recompute, compare stored JSON, print it; --write creates it once.
NumPy checks the finite ideal comparator. The infinite-dimensional clock is
covered by explicit analytical bounds, not a scattering discretization.
"""
from __future__ import annotations
import argparse
from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, localcontext
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def load(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))

def dec(x):
    if isinstance(x, F):
        return Decimal(x.numerator) / Decimal(x.denominator)
    return Decimal(str(x))

def out(x, places=15, upper=True):
    with localcontext() as ctx:
        ctx.prec = 80
        return str(dec(x).quantize(Decimal(1).scaleb(-places),
                   rounding=ROUND_CEILING if upper else ROUND_FLOOR))

def exp_h(h, t):
    e, u = np.linalg.eigh(h)
    return (u * np.exp(-1j * t * e)) @ u.conj().T

def entropy(rho):
    vals = np.linalg.eigvalsh((rho + rho.conj().T) / 2)
    assert vals.min() > -2e-12
    vals = vals[vals > 1e-14]
    return float(-np.dot(vals, np.log(vals)))

def classical_entropy(p):
    p = np.asarray(p)
    p = p[p > 0]
    return float(-np.dot(p, np.log(p)))

def clock_entropy_upper(epsilon):
    with localcontext() as ctx:
        ctx.prec = 80
        x = dec(epsilon) ** 2
        return -(x*x.ln()) - (1-x)*(1-x).ln() + x*Decimal(255).ln()

def compare(a, b, where="root"):
    if isinstance(a, dict):
        assert set(a) == set(b), where
        for k in a:
            compare(a[k], b[k], where+"."+k)
    elif isinstance(a, (list, tuple)):
        assert len(a) == len(b), where
        for i, (x,y) in enumerate(zip(a,b)):
            compare(x,y,f"{where}[{i}]")
    elif isinstance(a, float):
        assert math.isclose(a, b, rel_tol=1e-8, abs_tol=2e-12), (where,a,b)
    else:
        assert a == b, (where,a,b)

def run():
    source_path = STAGE/"1202"/"independent_results.json"
    old_path = STAGE/"1201"/"independent_results_1201.json"
    source, old = load(source_path), load(old_path)
    r = np.array(source["source_marginals"]["one"])
    pair = np.array(source["source_marginals"]["two"])
    assert np.max(abs(pair.sum(axis=1)-r)) < 2e-15
    assert abs(pair.sum()-1) < 2e-15
    assert pair.min() > 0

    I = np.eye(4)
    h = np.array([[0,-.2,0,-.05*np.sqrt(3)/4],
                  [-.2,1,0,0],[0,0,1,0],
                  [-.05*np.sqrt(3)/4,0,0,-.025]])
    q = np.zeros((4,4)); q[1,2] = q[2,1] = 1
    energy, basis = np.linalg.eigh(h)
    c = I[:,2]
    assert np.linalg.norm(h@c-c) < 1e-15
    H0 = np.kron(h,I)+np.kron(I,h)
    V = np.kron(q,q)
    HQ = H0+.1*V
    H064 = np.kron(H0,I)
    HQ64 = np.kron(HQ,I)

    G = np.zeros((64,64))
    copy = np.zeros((64,64))
    for j in range(4):
        S = I.copy()
        if j:
            S[[0,j],:] = S[[j,0],:]
        Pj = np.outer(basis[:,j],basis[:,j])
        G += np.kron(np.kron(Pj,I),math.pi/2*(I-S))
        copy += np.kron(np.kron(Pj,I),S)
    Uread = exp_h(G,1)
    g_eig = np.linalg.eigvalsh(G)
    read_error = float(np.linalg.norm(Uread-copy,2))
    comm_error = float(np.linalg.norm(G@H064-H064@G,2))
    assert g_eig.min() > -3e-14
    assert abs(g_eig.max()-math.pi) < 3e-14
    assert read_error < 2e-13 and comm_error < 2e-13
    copy_inputs = np.column_stack([
        np.kron(np.kron(basis[:,j],basis[:,k]),I[:,0])
        for j in range(4) for k in range(4)])
    copy_targets = np.column_stack([
        np.kron(np.kron(basis[:,j],basis[:,k]),I[:,j])
        for j in range(4) for k in range(4)])
    copy_isometry_error = float(np.linalg.norm(Uread@copy_inputs-copy_targets,2))
    assert copy_isometry_error < 2e-13

    T, M, sigma = F(606), F(10**24), F(1,10**10)
    delta, g, pi_upper, width = F(1,1000), F(1,10), F(22,7), F(547,500)
    assert energy[0] > -.055
    assert energy[-1]-energy[0] < float(width)
    spread = 1 + T*T/(8*M*M*sigma**4)
    bT = sigma*spread
    # q(u)=3u^2-2u^3 has integral 1/2 and max derivative 3/2.
    integral_q = 3*F(1,3)-2*F(1,4)
    assert integral_q == F(1,2)
    lipF, lipJ = F(3,2)/delta, F(3,2)
    assert g*lipF == 150
    assert 2*integral_q == 1
    assert 2*delta*(1-integral_q) == delta

    reports = {}
    ideal = {}
    epsilons = []
    entropies = []
    all_kernels = []
    for v in (1,2):
        # All branches are evaluated at the same physical final time 606.
        # The v=2 branch is NOT shut off at 303.
        lv = g*lipF+v*pi_upper*lipJ
        eps_clock = lv*T*bT
        eps_shape = g*delta/v
        epsilon = eps_clock+eps_shape
        epsilons.append(epsilon)
        sent = clock_entropy_upper(epsilon)
        entropies.append(sent)
        position = bT+lv*T*T/(2*M)
        tail0, tailT = lv*sigma, lv*position
        energy_error = 2*width*epsilon+tail0+tailT
        hc_mean = M*v*v/2+1/(8*M*sigma*sigma)
        hc_var = F(v*v,4)/sigma**2+1/(32*M*M*sigma**4)

        duration = 600/v
        Uid = exp_h(H064,606-601/v) @ Uread @ exp_h(HQ64,duration) @ exp_h(H064,1/v)
        inputs = np.column_stack([
            np.kron(np.kron(c,basis[:,j]),I[:,0]) for j in range(4)])
        states = Uid@inputs
        assert np.linalg.norm(states.conj().T@states-I,2) < 3e-11
        kernel = np.zeros((4,4))
        de_a = np.zeros(4); de_b = np.zeros(4)
        rho_ak, rho_b, rho_k = [],[],[]
        for j in range(4):
            wave = states[:,j].reshape(4,4,4)
            rhoA = np.einsum("abk,cbk->ac",wave,wave.conj())
            rhoB = np.einsum("abk,ack->bc",wave,wave.conj())
            rhoK = np.einsum("abk,abl->kl",wave,wave.conj())
            rhoAK = np.einsum("abk,cbl->akcl",wave,wave.conj()).reshape(16,16)
            kernel[:,j] = np.real(np.diag(rhoK))
            de_a[j] = np.trace(h@rhoA).real-1
            de_b[j] = np.trace(h@rhoB).real-energy[j]
            rho_ak.append(rhoAK); rho_b.append(rhoB); rho_k.append(rhoK)
        old_kernel = np.array(old["population_readout"]["matrix_rows_A_output_columns_B_population"]
                         if v == 1 else old["predeclared_holdout"]["matrix"])
        kernel_difference = float(np.max(abs(kernel-old_kernel)))
        assert kernel_difference < 2e-11
        assert np.max(abs(kernel.sum(axis=0)-1)) < 2e-11
        W = float((de_a+de_b)@r)
        mean_a, mean_b = float(de_a@r),float(de_b@r)
        marg_ak = sum(r[j]*rho_ak[j] for j in range(4))
        marg_b = sum(r[j]*rho_b[j] for j in range(4))
        marg_k = sum(r[j]*rho_k[j] for j in range(4))
        # AK entropy is the joint quantum entropy, not merely H(pointer labels).
        reports[str(v)] = {
            "contact_duration": duration, "common_final_time":606,
            "L_upper_exact":str(lv),"L_upper":float(lv),
            "clock_isometry_bound_exact":str(eps_clock),
            "shape_isometry_bound_exact":str(eps_shape),
            "total_isometry_bound_exact":str(epsilon),
            "total_isometry_bound":float(epsilon),
            "total_isometry_bound_up_12dp":out(epsilon,12),
            "free_clock_width_factor_upper_exact":str(spread),
            "position_rms_bound":float(position),
            "initial_interaction_tail_bound":float(tail0),
            "final_interaction_tail_bound":float(tailT),
            "initial_clock_mean_exact":str(hc_mean),
            "initial_clock_variance_exact":str(hc_var),
            "record_kernel":kernel.tolist(),
            "record_kernel_vs_1201_max_error":kernel_difference,
            "pointer_probabilities":(kernel@r).tolist(),
            "energy_coefficients":{"delta_A":de_a.tolist(),"delta_B":de_b.tolist(),
                                   "bare_AB_gain":(de_a+de_b).tolist()},
            "ideal_delta_A":mean_a,"ideal_delta_B":mean_b,"ideal_bare_AB_gain":W,
            "clock_energy_conservation_error_bound":float(energy_error),
            "clock_energy_interval_including_1e_minus_10_numeric_allowance":
                [-W-float(energy_error)-1e-10,-W+float(energy_error)+1e-10],
            "ideal_joint_AK_entropy":entropy(marg_ak),
            "ideal_B_entropy":entropy(marg_b),
            "ideal_pointer_entropy":entropy(marg_k),
            "clock_rank_upper":256,
            "clock_orthogonal_probability_bound":float(epsilon*epsilon),
            "clock_entropy_bound_decimal80_diagnostic":str(sent),
            "clock_entropy_upper_outward_15dp":out(sent+Decimal("1e-60"),15)
        }
        ideal[v]={"rho_ak":rho_ak,"rho_b":rho_b,"rho_k":rho_k,
                  "delta_b":mean_b}
        all_kernels.append(kernel)

    joint = all_kernels[0]@pair@all_kernels[1].T
    wrong_independent = np.outer(all_kernels[0]@r, all_kernels[1]@r)
    old_joint=np.array(source["original_1201_connection"]["joint_600_300"])
    joint_delta=float(np.max(abs(joint-old_joint)))
    assert joint_delta < 2e-11
    tv_model=float(np.sum(abs(joint-wrong_independent))/2)
    eps_sum=sum(epsilons)
    # This allowance is inherited from the accepted previous finite tables;
    # NumPy agreement alone is not a new enclosure of those tables.
    old_table_allowance=F(24,10**10)
    pointer_budget=eps_sum+old_table_allowance

    beta=1.2 # Fixed accepted comparison beta=6/5, not fitted beta.
    tau=np.exp(-beta*energy); tau/=tau.sum()
    Hpair=classical_entropy(pair.ravel())
    Dpair=float(np.sum(pair*np.log(pair/np.outer(tau,tau))))
    Dsingle=float(np.sum(r*np.log(r/tau)))
    source_mi=float(np.sum(pair*np.log(pair/np.outer(r,r))))
    assert abs(Dpair-2*Dsingle-source_mi) < 1e-13
    assert Dpair < .000912

    rhoAKpair=np.zeros((256,256),dtype=complex)
    rhoBpair=np.zeros((16,16),dtype=complex)
    for j in range(4):
        for k in range(4):
            rhoAKpair+=pair[j,k]*np.kron(ideal[1]["rho_ak"][j],ideal[2]["rho_ak"][k])
            rhoBpair+=pair[j,k]*np.kron(ideal[1]["rho_b"][j],ideal[2]["rho_b"][k])
    Sak,Sb=entropy(rhoAKpair),entropy(rhoBpair)
    eb_initial=2*float(energy@r)
    eb_final=float(np.trace(H0@rhoBpair).real)
    delta_eb=eb_final-eb_initial
    logZ=math.log(float(np.sum(np.exp(-beta*energy))))
    Dout=-Sb+beta*eb_final+2*logZ
    Iout=Sak+Sb-Hpair
    lhs=Sak+beta*delta_eb
    rhs=Iout+Dout-Dpair
    assert abs(lhs-rhs) < 2e-10
    with localcontext() as ctx:
        ctx.prec=80
        sent_sum=sum(entropies)
        actual_AK_lower=-Decimal("0.000912")-sent_sum
    inverse_bound=F(314,100)
    # Retain the previously accepted .015701 population budget in full;
    # only the new autonomous branch-1 error is transported through 3.14.
    population_upper=F(15701,10**6)+inverse_bound*epsilons[0]
    assert float(population_upper)<.016045

    result={
      "status":"independent 1203 finite comparator and analytic-bound audit; no clock PDE simulation",
      "all_checks_passed":True,
      "inputs":{
        "source_1202_independent_result_sha256":sha(source_path),
        "source_1201_independent_result_sha256":sha(old_path),
        "A_preparation_native_index":2,"A_preparation_energy":1,
        "B_pool_preparation_before_factory":"e0, not A's c state",
        "fixed_beta":1.2,"h":h.tolist(),"Q":q.tolist(),"energies":energy.tolist(),
        "source_one":r.tolist(),"source_two":pair.tolist(),
        "M_exact":str(M),"sigma_exact":str(sigma),"T_exact":str(T)
      },
      "finite_read_gate":{
        "dimension_ABK":64,"minimum_eigenvalue_numeric":float(g_eig.min()),
        "maximum_eigenvalue_numeric":float(g_eig.max()),
        "commutator_with_H0_norm":comm_error,
        "exponential_vs_controlled_swap_norm":read_error,
        "full_prepared_pointer_copy_isometry_error":copy_isometry_error,
        "program_compensation_deleted_v2_nonzero_label_copy_success_exact":"1/2",
        "initial_c_pointer0_G_expectation":math.pi/2
      },
      "gate_exact":{
        "smoothstep_integral_0_1":str(integral_q),
        "F_lipschitz":str(lipF),"J_lipschitz":str(lipJ),
        "J_integral":str(2*integral_q),
        "F_rectangular_L1_difference":str(delta),
        "semibounded_total_H_lower_greater_than":"-21/100",
        "constant_energy_shift_for_positive_reference":"1"
      },
      "branches":reports,
      "joint":{
        "pointer_table":joint.tolist(),"old_pointer_table_max_difference":joint_delta,
        "isometry_bound_exact":str(eps_sum),"isometry_bound":float(eps_sum),
        "isometry_bound_up_12dp":out(eps_sum,12),
        "inherited_old_table_allowance_exact":str(old_table_allowance),
        "pointer_TV_budget_up_12dp":out(pointer_budget,12),
        "wrong_independent_source_model_TV":tv_model,
        "model_separation_after_two_channel_error_budgets":tv_model-2*float(eps_sum),
        "separation_smaller_than_original_readout_budget":True,
        "branch1_population_infinity_budget":float(population_upper)
      },
      "entropy":{
        "fixed_beta":beta,"Gibbs_control_only":tau.tolist(),
        "D_source_single_numeric":Dsingle,"D_source_pair_numeric":Dpair,
        "source_pair_mutual_information":source_mi,
        "D_source_pair_inherited_certified_upper":.000912,
        "ideal_joint_AK_entropy":Sak,"ideal_joint_B_entropy":Sb,
        "ideal_source_delta_E":delta_eb,"ideal_mutual_information_AK_B":Iout,
        "ideal_D_B_out":Dout,"ideal_identity_residual":lhs-rhs,
        "clock_entropy_sum_decimal80_diagnostic":str(sent_sum),
        "clock_entropy_sum_upper_15dp":out(sent_sum+Decimal("1e-60"),15),
        "actual_AK_plus_beta_delta_EB_lower_decimal80_diagnostic":str(actual_AK_lower),
        "actual_AK_plus_beta_delta_EB_lower_outward_15dp":out(actual_AK_lower-Decimal("1e-60"),15,False),
        "clock_rank_proof_uses_A_K_C_pure_and_B_marginal_rank_at_most":4,
        "joint_AKC_rank_upper":256,
        "clock_projection_uses_squared_isometry_error":True
      },
      "scope":{
        "author_1203_code_imported":False,
        "unknown_input_reference_uniform_isometry_bound_analytical":True,
        "continuous_clock_evolution_numerically_solved":False,
        "two_branches_common_final_time":606,
        "actual_Gaussian_tails_discarded":False,
        "initial_read_interaction_claimed_exactly_zero":False,
        "R_omitted_or_pool_sources_assumed_independent":False,
        "physical_A_energy_pointer_constructed_in_declared_H":True,
        "B_energy_pointer_or_precision_clock_energy_instrument_constructed":False,
        "1202_factory_autonomous_controller_completed":False,
        "clock_or_pointer_reset_free":False,
        "internal_X_identified_with_space_or_metric":False,
        "full_final_admission_or_navigation_changed":False
      }
    }
    return result

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--write",action="store_true")
    args=p.parse_args()
    result=run()
    dest=HERE/"independent_results.json"
    if args.write:
        if dest.exists():
            raise SystemExit("--write refuses to overwrite existing results")
        dest.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    else:
        compare(result,load(dest))
    print(json.dumps(result,ensure_ascii=False,indent=2))

if __name__=="__main__":
    main()
