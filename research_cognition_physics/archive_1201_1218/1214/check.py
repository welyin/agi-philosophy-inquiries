"""Read-only author replay for round 1214; --write creates results once.
All sources are relative to research_cognition_physics. No earlier program imports.
Numerical small models diagnose analytic identities; they are not general proofs
or simulations of the continuous clock.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import itertools
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parents[1]
STAGE = HERE.parent
TOL = 3e-10
SOURCES = [
    "archive_370_428/research_note_421.md",
    "archive_429_466/research_note_457.md",
    "archive_1009_1043/research_note_1020.md",
    "archive_1009_1043/research_note_1037.md",
    "archive_1201_/_shared/after1208_reference_contract_admission.md",
    "archive_1201_/1202/results.json",
    "archive_1201_/1203/proof.md",
    "archive_1201_/1203/results.json",
    "archive_1201_/1204/proof.md",
    "archive_1201_/1204/results.json",
    "archive_1201_/1205/proof.md",
    "archive_1201_/1209/proof.md",
    "archive_1201_/1209/results.json",
    "archive_1201_/1210/proof.md",
    "archive_1201_/1211/proof.md",
    "archive_1201_/1212/proof.md",
    "archive_1201_/1213/proof.md",
]

def load_round(n):
    return json.loads((STAGE / str(n) / "results.json").read_text(encoding="utf-8-sig"))

def require(condition, label):
    if not bool(condition):
        raise AssertionError(label)

def near(a, b, label, tol=TOL):
    err = float(np.max(np.abs(np.asarray(a) - np.asarray(b))))
    require(err < tol, f"{label}: {err}")
    return err

def scalar(x):
    return round(float(x), 14)

def entropy(rho):
    p = np.linalg.eigvalsh((rho + rho.conj().T) / 2)
    require(float(p.min()) > -TOL, "density positivity")
    p = p[p > 1e-14]
    return -float(np.dot(p, np.log(p)))

def distance(a, b):
    return float(np.sum(np.abs(np.linalg.eigvalsh((a-b + (a-b).conj().T)/2))) / 2)

def partial(rho, ds, da, keep):
    t = rho.reshape(ds, da, ds, da)
    return np.trace(t, axis1=1, axis2=3) if keep == 0 else np.trace(t, axis1=0, axis2=2)

def apply(ks, x):
    return sum((k @ x @ k.conj().T for k in ks), np.zeros_like(x, dtype=complex))

def kraus_from_unitary(w, ds, da, env_prob):
    a = w.reshape(ds, da, ds, da)
    return [math.sqrt(env_prob[k]) * a[:, j, :, k] for j in range(da) for k in range(da)]

def unital_errors(ks, d):
    return {
        "trace_preservation_error": scalar(near(sum(k.conj().T @ k for k in ks), np.eye(d), "TP")),
        "unital_error": scalar(near(sum(k @ k.conj().T for k in ks), np.eye(d), "unital")),
    }

def exp_h(h, time):
    e, u = np.linalg.eigh(h)
    return (u * np.exp(-1j*time*e)) @ u.conj().T

def random_density(rng, d):
    z = rng.normal(size=(d,d)) + 1j*rng.normal(size=(d,d))
    a = z @ z.conj().T
    return a / np.trace(a)

def thermal(energies, beta):
    p = np.exp(-beta*np.asarray(energies))
    return p / p.sum()

def rel_to_diag(rho, q):
    return -entropy(rho) - float(np.dot(np.diag(rho).real, np.log(q)))

def exact_native(cold, hot, clock):
    rc = [F(x) for x in cold["one_marginal_rational"]]
    rh = [F(x) for x in hot["new_reference"]["one_marginal_rational"]]
    require(sum(rc) == sum(rh) == 1 and min(rc+rh) > 0, "exact populations")
    lam = min(rc)*min(rh)
    expected = F(
        "420812791211495342505037876598058009513148297805273888857249161790607/"
        "18746287226936270671051661441770688017997733890644189549902447495884150")
    require(lam == expected, "exact native spectral floor")
    t, mass, sigma = F(606), F(10**24), F(1,10**10)
    bt = sigma*(1+t*t/(8*mass*mass*sigma**4))
    rows = []
    for v, lower in [(1, F("0.269334240")), (2, F("0.269352096"))]:
        saved = clock["bounds_by_program"][str(v)]
        lip = F(150)+F(3,2)*v*F(22,7)
        ec = lip*t*bt
        es = F(1,10000*v)
        exact = ec+es
        upper = F(saved["isometry_error_upper"])
        require(ec == F(saved["clock_error_exact_upper"]), "clock error reconstruction")
        require(es == F(saved["edge_error_exact"]), "shape term")
        require(exact == F(saved["isometry_error_exact_upper"]) and exact < upper, "total error")
        require(upper < F(3,4), "positive calibration regime")
        strong, direct = lam*(12-16*upper), 12*lam-upper
        require(strong > lower and strong >= direct > 0, "strict native failure")
        rows.append({
            "v": v, "QQ_time": 600//v,
            "clock_error_exact": str(ec), "shape_error_exact": str(es),
            "half_diamond_exact_upper": str(exact), "outward_error_upper": str(upper),
            "normalised_unital_calibration_upper": str(upper),
            "strong_failure_exact_lower_using_outward_error": str(strong),
            "strong_failure_decimal": scalar(strong),
            "strict_decimal_lower": str(lower),
            "direct_half_diamond_failure_lower": str(direct),
            "strong_lower_using_exact_clock_formula": str(lam*(12-16*exact)),
            "same_bound_for_actual_K_nonzero_event_separately": True,
        })
    return rc, rh, {
        "name": "exact_native_population_and_clock_calibration", "passed": True,
        "cold_minimum_index": rc.index(min(rc)), "hot_minimum_index": rh.index(min(rh)),
        "joint_lambda_exact": str(lam), "joint_lambda_decimal": scalar(lam),
        "dimension": 16, "bad_rank": 12, "independent_CK_required": True,
        "source_state_diagonal_required_for_population_floor": True,
        "actual_1205_correlated_preparation_channel_claimed": False, "bounds": rows,
    }

def first_hit_lift():
    # Complete basis, including both ladder ends and both pointer inputs.
    words = list(itertools.product(range(3), repeat=2))
    index = {s:i for i,s in enumerate(words)}
    d = len(words)
    material_energy = np.array([sum(s) for s in words], dtype=int)
    perm = list(range(d))
    event = [0]*d
    for i, s in enumerate(words):
        hit = next((j for j,z in enumerate(s) if z in (0,1)), None)
        if hit is not None:
            t = list(s); t[hit] = 1-t[hit]
            perm[i] = index[tuple(t)]
            event[i] = int(s[hit] == 0)
    near(np.array(perm)[perm], np.arange(d), "material involution")
    out = []
    for length in (4, 7):
        ladder = length+2
        da = ladder*2
        joint = d*da
        def pos(s,k,q): return (s*ladder+k)*2+q
        lift = np.arange(joint)
        pairs = 0
        for s in range(d):
            u = perm[s]
            if s != u and material_energy[s] < material_energy[u]:
                require(material_energy[u]-material_energy[s] == 1, "ladder gap")
                for k in range(1, length+2):
                    for q in range(2):
                        a,b = pos(s,k,q),pos(u,k-1,q)
                        lift[a],lift[b] = b,a
                        pairs += 1
        full = np.empty(joint, dtype=int)
        energies = np.empty(joint, dtype=int)
        for i in range(joint):
            s,kq = divmod(i,da); k,q = divmod(kq,2)
            energies[i] = material_energy[s]+k
            z = int(lift[i]); so,kqo = divmod(z,da); ko,qo = divmod(kqo,2)
            full[i] = pos(so,ko,qo ^ event[so])
        require(len(set(full.tolist())) == joint, "complete W permutation")
        require(np.array_equal(energies[full], energies), "exact full energy conservation")
        v = np.zeros((joint,d), dtype=complex)
        for s in range(d):
            for k in range(1,length+1):
                v[full[pos(s,k,0)],s] = 1/math.sqrt(length)
        near(v.conj().T @ v, np.eye(d), "prepared isometry")
        aux = np.zeros((d,da), dtype=complex)
        for s in range(d):
            block = v[:,s].reshape(d,da)
            aux[s] = block[perm[s]]
            near(np.sum(np.abs(block)**2,axis=1), np.eye(d)[perm[s]], "full support material T")
        gram = aux @ aux.conj().T
        near(np.diag(gram), np.ones(d), "Gram diagonal")
        require(np.linalg.eigvalsh(gram).min() > -TOL, "Gram PSD")
        ks = [v.reshape(d,da,d)[:,j,:] for j in range(da)]
        errs = unital_errors(ks,d)
        tm = np.eye(d)[:,perm]
        # Every matrix unit: tests all coherent inputs, not just populations.
        gram_err = 0.
        for i in range(d):
            for j in range(d):
                x = np.zeros((d,d), complex); x[i,j] = 1
                gram_err = max(gram_err, near(apply(ks,x),tm@(gram*x)@tm.T, "Gram all matrix units"))
        p = np.array([.6,.25,.15])
        px = np.array([p[s[0]]*p[s[1]] for s in words])
        omega = (v*px)@v.conj().T
        material = partial(omega,d,da,0)
        auxiliary = partial(omega,d,da,1)
        near(material,tm@np.diag(px)@tm.T, "exact reference population permutation")
        near(entropy(material),entropy(np.diag(px)), "material entropy unchanged")
        rank = int(np.sum(np.linalg.eigvalsh(auxiliary)>1e-10))
        require(rank <= 3, "three pure auxiliary branches")
        info = entropy(material)+entropy(auxiliary)-entropy(omega)
        near(info,entropy(auxiliary),"I(X:A)=S(A)")
        out.append({
            "L":length,"joint_dimension":joint,"complete_swapped_pairs":pairs,
            "both_ladder_boundary_levels_included":True,"all_pointer_inputs_included":True,
            "max_matrix_unit_Gram_error":scalar(gram_err),"auxiliary_rank":rank,
            "material_entropy_change":scalar(entropy(material)-entropy(np.diag(px))),
            "I_X_A":scalar(info),"auxiliary_entropy":scalar(entropy(auxiliary)),**errs
        })
    return {"name":"complete_ladder_copy_all_input_Gram_channel","passed":True,"small_models":out,
            "general_classification_proved_analytically_not_by_finite_models":True}

def gibbs_and_zeroth():
    rng = np.random.default_rng(121403)
    d=3; energies=np.arange(d); total=(energies[:,None]+energies[None,:]).ravel()
    w=np.zeros((d*d,d*d),complex)
    for e in np.unique(total):
        idx=np.flatnonzero(total==e)
        a=rng.normal(size=(len(idx),len(idx)))+1j*rng.normal(size=(len(idx),len(idx)))
        block=exp_h((a+a.conj().T)/2,.73)
        w[np.ix_(idx,idx)]=block
    near(w.conj().T@w,np.eye(d*d),"joint W unitary")
    near(w@np.diag(total),np.diag(total)@w,"joint W exact block energy")
    aux0=np.eye(d)/d
    ks=kraus_from_unitary(w,d,d,[1/d]*d)
    errors=unital_errors(ks,d)
    rows=[]
    for beta in (-1.2,0.,1/3,1.2):
        p=thermal(energies,beta); initial=np.kron(np.diag(p),aux0)
        final=w@initial@w.conj().T
        sx=partial(final,d,d,0); ay=partial(final,d,d,1)
        de=float(np.dot(energies,np.diag(ay).real-1/d))
        ds=entropy(sx)-entropy(np.diag(p))
        dr=rel_to_diag(sx,p)
        near(dr,-beta*de-ds,"strengthened Gibbs identity")
        require(ds > -TOL and dr <= abs(beta)*abs(de)+TOL,"strengthened Gibbs bound")
        trace=distance(sx,np.diag(p))
        require(trace*trace <= abs(beta)*abs(de)/2+TOL,"Pinsker bound")
        if beta == 0:
            near(sx,np.eye(d)/d,"beta zero exactly fixed")
        rows.append({"beta":beta,"Delta_E_aux":scalar(de),"Delta_S_material":scalar(ds),
                     "relative_entropy":scalar(dr),"abs_beta_energy_bound":scalar(abs(beta*de)),
                     "material_distance":scalar(trace)})
    e1=np.array([0.,1.]); e2=np.array([0.,.4,1.7])
    z=[]
    for beta in (-.8,0.,.7):
        p=np.kron(thermal(e1,beta),thermal(e2,beta))
        q=thermal((e1[:,None]+e2[None,:]).ravel(),beta)
        z.append({"beta":beta,"common_beta_factorization_error":scalar(near(p,q,"common beta"))})
    unequal=np.kron(thermal(e1,.4),thermal(e1,.9))
    require(abs(unequal[1]-unequal[2])>0.05,"unequal beta degenerate populations")
    return {"name":"strengthened_Gibbs_and_common_beta_reuse","passed":True,
            "all_input_unital":errors,"Gibbs_diagnostics":rows,
            "common_beta_diagnostics":z,"positive_beta_not_selected_by_unitality":True}

def clausius():
    rng=np.random.default_rng(121404)
    d=3; e=np.arange(d); beta=.7
    swap=np.zeros((d*d,d*d),complex)
    for i in range(d):
        for j in range(d): swap[j*d+i,i*d+j]=1
    near(swap@np.diag((e[:,None]+e[None,:]).ravel()),
         np.diag((e[:,None]+e[None,:]).ravel())@swap,"SWAP exact energy")
    rho=random_density(rng,d); tau=np.diag(thermal(e,beta))
    initial=np.kron(rho,tau)
    rows=[]
    for p in (0.,.37,1.):
        # A classical zero-energy control of I/SWAP gives this unital joint channel.
        final=p*initial+(1-p)*swap@initial@swap.T
        sf,bf=partial(final,d,d,0),partial(final,d,d,1)
        delta_joint=entropy(final)-entropy(initial)
        info=entropy(sf)+entropy(bf)-entropy(final)
        db=rel_to_diag(bf,np.diag(tau).real)
        lhs=entropy(sf)-entropy(rho)+beta*float(np.dot(e,np.diag(bf-tau).real))
        near(lhs,delta_joint+info+db,"Clausius exact balance")
        require(lhs>=-TOL and delta_joint>=-TOL,"joint unital Clausius")
        rows.append({"identity_weight":p,"lhs":scalar(lhs),"Delta_S_SB":scalar(delta_joint),
                     "mutual_information":scalar(info),"bath_relative_entropy":scalar(db)})
    reduced_of_maxmix=partial(swap@np.kron(np.eye(d)/d,tau)@swap.T,d,d,0)
    near(reduced_of_maxmix,tau,"local replacement")
    nonunital=distance(reduced_of_maxmix,np.eye(d)/d)
    require(nonunital>.1,"SB unital does not imply S unital")
    return {"name":"joint_domain_Clausius_and_local_nonunital_deletion","passed":True,
            "balances":rows,"local_normalized_unital_deviation":scalar(nonunital),
            "energy_conserving_I_SWAP_control_exists":True}

def calibration_and_deletions():
    d=4; lam=.04; rho=np.diag([.5,.3,.16,.04]); bad=np.diag([0.,0.,1.,1.])
    reset=np.diag([1.,0.,0.,0.])
    rows=[]
    for t in (0.,.1,.5,.8,1.):
        def phi(x): return (1-t)*x+t*np.trace(x)*reset
        output=phi(rho); cal=distance(phi(np.eye(d)/d),np.eye(d)/d)
        floor=lam*d*max(0.,.5-cal)
        actual=float(np.trace(bad@output).real)
        require(actual+TOL>=floor,"positive-map calibration bound")
        near(phi(rho)-lam*phi(np.eye(d)),phi(rho-lam*np.eye(d)),"positive map decomposition")
        require(np.linalg.eigvalsh(output-lam*phi(np.eye(d))).min()>-TOL,"order preservation")
        rows.append({"reset_weight":t,"calibration_half_distance":scalar(cal),
                     "calibration_failure_floor":scalar(floor),"actual_failure":scalar(actual)})
    reset_cal=distance(reset,np.eye(4)/4)
    near(reset_cal,.75,"1213 reset threshold")
    # Old correlated example: all one-body marginals faithful, joint no-ground weight zero.
    corr=np.zeros((16,16))
    for j in range(1,4):
        corr[j,j]=1/6
        corr[4*j,4*j]=1/6
    marginal=partial(corr,4,4,0)
    pno=np.diag([float(i>0 and j>0) for i in range(4) for j in range(4)])
    near(np.trace(corr@pno),0,"old correlated pool no-ground projection")
    require(np.linalg.eigvalsh(marginal).min()>0 and np.linalg.matrix_rank(corr)<16,
            "marginal faithful not joint faithful")
    route=[]
    for n in (1,2,3):
        ground=1-(F(3,4))**n
        require(ground == F(1,4) if n==1 else ground>F(1,4),"pure record routing not unital")
        route.append({"n":n,"maxmixed_routing_ground_probability":str(ground)})
    return {"name":"finite_faithfulness_calibration_and_old_deletions","passed":True,
            "positive_map_diagnostics":rows,"reset_calibration_deviation":scalar(reset_cal),
            "correlated_joint_rank":int(np.linalg.matrix_rank(corr)),
            "correlated_marginal_min_eigenvalue":scalar(np.linalg.eigvalsh(marginal).min()),
            "old_pure_record_routing":route,
            "universal_finite_time_unattainability_claimed":False}

def native_channels(hot, rc, rh):
    h=np.array([[0,-1/5,0,-math.sqrt(3)/80],
                [-1/5,1,0,0],[0,0,1,0],[-math.sqrt(3)/80,0,0,-1/40]],float)
    q=np.zeros((4,4)); q[1,2]=q[2,1]=1
    e,basis=np.linalg.eigh(h)
    qe=basis.T@q@basis
    h0=np.diag((e[:,None]+e[None,:]).ravel())
    qq=np.kron(qe,qe)
    rho=np.diag(np.kron([float(x) for x in rc],[float(x) for x in rh]))
    projections=[np.diag([float(a==j) for a in range(4) for b in range(4)]) for j in range(4)]
    rows=[]
    for v in (1,2):
        time=600/v
        u=exp_h(h0+.1*qq,time)
        saved=hot["full_native_contacts"][str(int(time))]
        transition=np.array([[float(F(x)) for x in row] for row in saved["full_transition_P"]])
        transition_error=near(np.abs(u)**2,transition,"old complete QQ response",tol=1e-9)
        uv=exp_h(h0,606-601/v)@u@exp_h(h0,1/v)
        ks=[p@uv for p in projections]
        errs=unital_errors(ks,16)
        near(apply(ks,np.eye(16)/16),np.eye(16)/16,"native maxmixed calibration")
        vv=np.zeros((64,16),complex)
        for a in range(4):
            for b in range(4):
                vv[(4*a+b)*4+a,:]=uv[4*a+b,:]
        near(vv.conj().T@vv,np.eye(16),"actual ideal copy isometry")
        output=apply(ks,rho)
        full=(vv@rho)@vv.conj().T
        near(partial(full,16,4,0),output,"copy full-to-material reduction")
        pointer=partial(full,16,4,1)
        kernel=np.array([[float(F(x)) for x in row] for row in saved["pointer_K"]])
        predicted=kernel@np.diag(rho).real
        pointer_error=near(np.diag(pointer).real,predicted,"old pointer distribution",tol=1e-9)
        # Matrix-unit identity with the full phases retained.
        rng=np.random.default_rng(121406+v)
        z=rng.normal(size=(16,5))+1j*rng.normal(size=(16,5))
        z/=np.linalg.norm(z)
        full_reference=vv@z
        near(np.linalg.norm(full_reference),1,"passive unknown reference normalization")
        pexc=float(np.trace((np.eye(16)-projections[0])@output).real)
        entropy_change=entropy(output)-entropy(rho)
        require(entropy_change>=-TOL and pexc>12*min(float(x) for x in rc)*min(float(x) for x in rh)-TOL,
                "native ideal unital spectrum floor")
        rows.append({"v":v,"QQ_time":int(time),"ideal_A_excitation":scalar(pexc),
                     "ideal_K_nonzero":scalar(1-pointer[0,0].real),
                     "material_entropy_change":scalar(entropy_change),
                     "old_transition_max_error":scalar(transition_error),
                     "old_pointer_max_error":scalar(pointer_error),**errs})
    return {"name":"native_full_QQ_free_phases_and_copy_channel","passed":True,
            "energy_eigenvalues":[scalar(x) for x in e],"diagnostics":rows,
            "continuous_clock_numerically_simulated":False,
            "physical_actual_half_diamond_from_1203_analytic_proof":True}

def native_record_certificate(cold, old):
    r=[F(x) for x in cold["one_marginal_rational"]]
    cert=old["native_record_inclusive_certificate"]
    nx,ny=cert["x_occupation"],cert["y_occupation"]
    a=math.prod(p**n for p,n in zip(r,nx))
    b=math.prod(p**n for p,n in zip(r,ny))
    require(a>b>0,"native first-hit probabilities")
    q=a+b
    blocks=(q.denominator+q.numerator-1)//q.numerator
    require(str(blocks)==cert["blocks"],"frozen block count")
    length=int(cert["ladder_L"])
    infidelity=F(2,length)+b/q
    require(infidelity<F(cert["auxiliary_infidelity_upper"]),"native auxiliary infidelity")
    require(infidelity<F(cert["auxiliary_trace_distance_upper"])**2,"native trace budget")
    # e > 163/60 by its positive series, and Mq>=1: hit >103/163.
    require(sum(F(1,math.factorial(j)) for j in range(6))==F(163,60),"e-series rational")
    drift=F(103,163)*(a-b)/q
    require(drift>F(cert["material_distance_strict_lower"]),"native drift lower")
    emax=F(cert["auxiliary_infidelity_upper"])
    exp_lower=lambda x,n:sum((x**k/F(math.factorial(k)) for k in range(n+1)),F(0))
    require(exp_lower(F("7.71"),60)>1/emax,"log e certification")
    require(exp_lower(F(".694"),20)>2,"log2 certification")
    sbound=emax*(F("7.71")+1+F(".694"))
    require(sbound==F(cert["entropy_upper_exact"]) and sbound<F(".00423"),"entropy bound")
    gaps=[F(x) for x in cert["energy_gap_interval"]]
    require(max(abs(x) for x in gaps)<F(cert["energy_change_absolute_upper"]),"native energy budget")
    return {"name":"old_native_record_certificate_in_unital_menu","passed":True,
            "block_length":cert["block_length"],"blocks":str(blocks),
            "ladder_L":length,"record_dimension":2,
            "probability_ratio_a_over_b":scalar(a/b),
            "auxiliary_infidelity_decimal":scalar(infidelity),
            "auxiliary_trace_upper":cert["auxiliary_trace_distance_upper"],
            "auxiliary_entropy_rational_upper":str(sbound),
            "material_drift_certified_lower":cert["material_distance_strict_lower"],
            "same_all_input_Gram_proof_applies":True,
            "giant_operator_materialized":False,"original_QQ_implementation_claimed":False}

def compute():
    cold,clock,hot,old=load_round(1202),load_round(1203),load_round(1204),load_round(1209)
    rc,rh,g1=exact_native(cold,hot,clock)
    groups=[g1,first_hit_lift(),gibbs_and_zeroth(),clausius(),
            calibration_and_deletions(),native_channels(hot,rc,rh),native_record_certificate(cold,old)]
    require(len(groups)==7 and all(g["passed"] for g in groups),"all groups")
    return {
        "schema":"1214_unital_common_role_v1","round":1214,
        "status":"AUTHOR_CHECKS_PASS_NOT_ACCEPTANCE","all_checks_passed":True,
        "check_groups":7,"groups":groups,
        "source_path_base":"research_cognition_physics",
        "source_sha256":{s:hashlib.sha256((RESEARCH/s).read_bytes()).hexdigest() for s in SOURCES},
        "scope":{
            "new_adopted_axioms":False,
            "new_claim":"1209 inverse witnesses stay in unital common menu; native independent contact calibrated",
            "all_beta_reference_classification":True,
            "positive_temperature_branch_is_separate":True,
            "all_auxiliary_joint_budgets_rechecked_per_complete_process":True,
            "actual_independent_ABCK_only":True,
            "actual_1204_in_exact_1209_energy_menu_claimed":False,
            "actual_1205_prepared_correlated_cut_has_fixed_AB_channel_claimed":False,
            "general_finite_time_third_law_claimed":False,
            "clock_free_energy_equals_extracted_mechanical_work_claimed":False,
        },
    }

def compare(a,b,path="root"):
    require(type(a)==type(b),f"saved type mismatch {path}")
    if isinstance(a,dict):
        require(a.keys()==b.keys(),f"saved keys mismatch {path}")
        for k in a: compare(a[k],b[k],path+"."+str(k))
    elif isinstance(a,list):
        require(len(a)==len(b),f"saved length mismatch {path}")
        for i,(x,y) in enumerate(zip(a,b)): compare(x,y,f"{path}[{i}]")
    elif isinstance(a,float):
        require(math.isclose(a,b,rel_tol=2e-11,abs_tol=2e-11),f"saved float mismatch {path}")
    else:
        require(a==b,f"saved value mismatch {path}")

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--write",action="store_true",help="create results.json; refuse replacement")
    args=parser.parse_args()
    target=HERE/"results.json"
    if args.write and target.exists():
        raise SystemExit("--write refuses to overwrite existing results.json")
    result=compute()
    if args.write:
        with target.open("x",encoding="utf-8") as stream:
            json.dump(result,stream,ensure_ascii=False,indent=2,allow_nan=False)
            stream.write("\n")
        print(json.dumps({"status":"AUTHOR_RESULTS_WRITTEN_NOT_ACCEPTANCE","all_checks_passed":True,
                          "check_groups":7,"path":str(target)},ensure_ascii=False))
    else:
        if not target.exists():
            raise SystemExit("Missing results.json: first author creation requires --write")
        compare(result,json.loads(target.read_text(encoding="utf-8")))
        print(json.dumps({"status":"PASS","all_checks_passed":True,"results_match":True,
                          "check_groups":7,"mode":"read_only_replay"},ensure_ascii=False))

if __name__=="__main__":
    main()
