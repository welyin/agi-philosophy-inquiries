"""Independent 1274: original spin permutations and autonomous full-path maps.

No author checker is imported. Four/six-leg singlets are built from noncrossing
valence bonds, not a diagonalization of an author's Casimir construction.
Default execution compares a saved result without writing; --write creates it.
Floating diagnostics support, but do not replace, the analytic all-time proof.
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

def maxabs(x):
    return float(np.max(np.abs(x)))

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def pairings(labels):
    if not labels:
        yield []
        return
    a = labels[0]
    for j in range(1, len(labels), 2):
        b = labels[j]
        for left in pairings(labels[1:j]):
            for right in pairings(labels[j+1:]):
                yield [(a, b)] + left + right

def bond_vector(n, pairs):
    out = np.zeros(2**n)
    for choices in itertools.product((0, 1), repeat=len(pairs)):
        bits = [0]*n
        sign = 1
        for (i, j), bit in zip(pairs, choices):
            bits[i], bits[j] = bit, 1-bit
            sign *= (-1)**bit
        index = sum(bit << (n-1-i) for i, bit in enumerate(bits))
        out[index] += sign / math.sqrt(2**len(pairs))
    return out

def singlet_frame(n):
    vectors = np.column_stack([bond_vector(n, p) for p in pairings(list(range(n)))])
    q, r = np.linalg.qr(vectors)
    assert np.min(np.abs(np.diag(r))) > 0.1
    # Independent SU(2) check: Jz=0 and total raising annihilate every column.
    jz = np.array([sum(0.5 if not (i >> (n-1-j)) & 1 else -0.5
                       for j in range(n)) for i in range(2**n)])
    raising = np.zeros((2**n, 2**n))
    for i in range(2**n):
        for j in range(n):
            bit = 1 << j
            if i & bit:
                raising[i-bit, i] += 1
    assert maxabs(jz[:, None]*q) < 1e-13
    assert maxabs(raising @ q) < 1e-13
    return q

def permute_rows(q, pi):
    n = len(pi)
    index = np.arange(2**n)
    bits = (index[:, None] >> np.arange(n-1, -1, -1)) & 1
    destination = bits[:, pi] @ (2**np.arange(n-1, -1, -1))
    out = np.zeros_like(q)
    out[destination, :] = q
    return out

def entropy(matrix):
    values = np.linalg.eigvalsh((matrix+matrix.conj().T)/2)
    assert values.min() > -2e-12
    values = values[values > 1e-13]
    return float(-np.dot(values, np.log(values)))

def binary_entropy(p):
    return sum(-x*math.log(x) for x in (p, 1-p) if x > 1e-14)

def build_original_maps():
    q4, q6 = singlet_frame(4), singlet_frame(6)
    assert q4.shape == (16, 2) and q6.shape == (64, 5)
    w = q6.T @ np.kron(q4, np.array([0, 1, -1, 0])[:, None]/math.sqrt(2))
    assert maxabs(w.T@w-np.eye(2)) < 1e-13
    permutations = list(itertools.permutations(range(6)))
    reps = np.array([q6.T @ permute_rows(q6, p) for p in permutations])
    assert maxabs(np.einsum("fji,fjk->fik", reps, reps)-np.eye(5)) < 1e-12
    # Original B7 x W2 data space. All 720 flag values are kept in the isometry.
    s = np.zeros((14, 2), complex)
    t = np.zeros((14, 2), complex)
    s[1, 0] = s[3, 1] = 1
    t[4::2] = w
    u1 = np.eye(14)-s@s.conj().T-t@t.conj().T+t@s.conj().T+s@t.conj().T
    u2 = []
    for rep in reps:
        rb = np.zeros((7, 7))
        rb[:2, :2] = np.eye(2)
        rb[2:, 2:] = rep
        u2.append(np.kron(rb, np.eye(2)))
    u2 = np.array(u2)
    hb = np.diag([4.0]*2+[6.0]*5)
    hw = np.diag([0.0, 2.0])
    h0 = np.kron(hb, np.eye(2))+np.kron(np.eye(7), hw)
    assert maxabs(u1.conj().T@u1-np.eye(14)) < 1e-12
    assert maxabs(u1@s-t) < 1e-13
    assert maxabs(u1@h0-h0@u1) < 1e-12
    assert maxabs(np.einsum("fij,jk->fik",u2,h0)
                  -np.einsum("ij,fjk->fik",h0,u2)) < 1e-12
    paths = np.zeros((3, 14, 720, 2), complex)
    paths[0] = s[:, None, :]/math.sqrt(720)
    paths[1] = t[:, None, :]/math.sqrt(720)
    paths[2] = np.einsum("fij,jk->ifk",u2,t)/math.sqrt(720)
    return w, paths, u1, u2, h0


def full_thermal_reference(u1,u2,h0):
    # Compute exp(-beta H_total) from the actual 42D matrix separately for
    # every retained flag. This does not assume the claimed reduced Gibbs form.
    beta=math.log(2.5)/2
    target_bw=np.diag(np.exp(-beta*np.diag(h0)))
    target_bw/=np.trace(target_bw)
    target_b=np.diag([0.25]*2+[0.1]*5)
    reference_b=np.zeros((7,7),complex)
    max_marginal_error=0.0
    max_partition_difference=0.0
    first_partition=None
    for r in u2:
        h=np.kron(np.eye(3),h0).astype(complex)
        h[14:28,:14]=u1
        h[:14,14:28]=u1.conj().T
        h[28:,14:28]=r
        h[14:28,28:]=r.conj().T
        vals,vecs=np.linalg.eigh(h)
        weights=np.exp(-beta*vals)
        partition=float(weights.sum())
        if first_partition is None:first_partition=partition
        max_partition_difference=max(max_partition_difference,abs(partition-first_partition))
        density=(vecs*(weights/partition))@vecs.conj().T
        marginal=sum(density[j*14:(j+1)*14,j*14:(j+1)*14] for j in range(3))
        max_marginal_error=max(max_marginal_error,maxabs(marginal-target_bw))
        reduced=np.einsum("bwcw->bc",marginal.reshape(7,2,7,2))
        reference_b+=reduced/720
    ell=np.diag([math.log(2)]*2+[math.log(5)]*5)
    modular=-np.diag(np.log(np.diag(target_b)))
    modular_error=maxabs(modular-ell-math.log(2)*np.eye(7))
    assert max_marginal_error<2e-12
    assert max_partition_difference<2e-11
    assert maxabs(reference_b-target_b)<2e-12
    assert modular_error<1e-14
    return {"beta_star_in_inverse_e_units":beta,
            "direct_full_H_exponentiations":720,
            "each_flag_block_dimension":42,
            "max_BW_marginal_error":max_marginal_error,
            "all_flag_partition_difference":max_partition_difference,
            "B_reference_error":maxabs(reference_b-target_b),
            "B_reference_spectrum":[0.25]*2+[0.1]*5,
            "modular_operator_error":modular_error,
            "normal_full_H_reference_exists":True,
            "actual_reference_preparation_or_observer_temperature_certified":False}


def exact_bounds():
    ratio = F(2**33, 3**12*5**5)
    assert ratio > 5
    # ln r > 2(r-1)/(r+1), from integrating the positive atanh series.
    log_gap_lower = (ratio-1)/(10*(ratio+1))
    assert log_gap_lower > F(1, 15)
    # Each endpoint has half-trace error delta. ln6<2 and ln1000<7.
    # Both inequalities are independently bounded by finite positive exp sums.
    exp2_lower = sum((F(2)**j/F(math.factorial(j)) for j in range(5)), F())
    exp7_lower = sum((F(7)**j/F(math.factorial(j)) for j in range(19)), F())
    assert exp2_lower > 6 and exp7_lower > 1000
    delta = F(1,1000)
    each_entropy_error_upper = delta*(2+7+1)
    final_entropy_gap = F(1,15)-2*each_entropy_error_upper
    final_energy_gap = F(1,2)-4*delta
    assert final_entropy_gap == F(7,150)
    assert final_energy_gap == F(62,125)
    return {
        "ratio":str(ratio),
        "ratio_exceeds_five_exact":True,
        "gap_strict_lower_from_log_integral":str(log_gap_lower),
        "gap_exceeds":str(F(1,15)),
        "endpoint_half_trace_error":str(delta),
        "each_entropy_error_strict_upper":str(each_entropy_error_upper),
        "entropy_drop_strict_lower_with_both_endpoint_errors":str(final_entropy_gap),
        "area_and_energy_increase_lower_in_units_a_or_e":str(final_energy_gap),
        "normalization_of_error":"half trace distance on boundary B, dimension seven",
        "actual_clock_or_measurement_error_certified":False
    }

def compute():
    w, paths, u1, u2, h0 = build_original_maps()
    chain = np.array([[0,1,0],[1,0,1],[0,1,0]], float)
    eigenvalues, eigenvectors = np.linalg.eigh(chain)
    area = np.diag([4.0]*2+[6.0]*5)
    sigma4 = np.zeros((7,7)); sigma4[:2,:2] = np.eye(2)/2
    tau = np.zeros((7,7)); tau[2:,2:] = w@w.T/2
    sigma6 = np.zeros((7,7)); sigma6[2:,2:] = np.eye(5)/5
    station = np.diag([0.0,1.0,2.0])
    worst = dict(isometry=0.0, channel_matrix_units=0.0, spectrum=0.0,
                 station_probabilities=0.0, area_energy=0.0,
                 Schrodinger_full_path=0.0, full_resource_energy=0.0,
                 structure_account=0.0, full_reference_balance=0.0)
    rows = []
    endpoints = {}
    angles = (0, math.pi/6, math.pi/3, math.pi/2, 2*math.pi/3,
              3*math.pi/4, math.pi, 4*math.pi/3, 3*math.pi/2, 2*math.pi)
    for angle in angles:
        time = angle/math.sqrt(2)
        coeff = (eigenvectors*np.exp(-1j*eigenvalues*time))@eigenvectors.T[:,0]
        coeff_derivative = -1j*chain@coeff
        phase = np.exp(-6j*time)
        v = phase*coeff[:,None,None,None]*paths
        v_derivative = phase*(coeff_derivative-6j*coeff)[:,None,None,None]*paths
        # Actual full Hamiltonian action on all kept program/flag/data indices.
        hv = np.einsum("ij,cjfk->cifk",h0,v)
        hv[0] += np.einsum("ij,jfk->ifk",u1.conj().T,v[1])
        hv[1] += np.einsum("ij,jfk->ifk",u1,v[0])
        hv[1] += np.einsum("fji,jfk->ifk",u2.conj(),v[2])
        hv[2] += np.einsum("fij,jfk->ifk",u2,v[1])
        worst["Schrodinger_full_path"] = max(worst["Schrodinger_full_path"],maxabs(1j*v_derivative-hv))
        allmap = v.reshape(-1,2)
        worst["isometry"] = max(worst["isometry"],maxabs(allmap.conj().T@allmap-np.eye(2)))
        vbeq = v.reshape(3,7,2,720,2).transpose(1,0,2,3,4).reshape(7,4320,2)
        rho_b = np.einsum("bei,cei->bc",vbeq,vbeq.conj())/2
        c = math.cos(angle)
        p0,p1,p2 = (1+c)**2/4,(1-c*c)/2,(1-c)**2/4
        probs = np.sum(np.abs(v)**2,axis=(1,2,3))/2
        worst["station_probabilities"] = max(worst["station_probabilities"],maxabs(probs-np.array([p0,p1,p2])))
        # All four input matrix units verify the reduced map, not just sigma4.
        for i in range(2):
            for j in range(2):
                raw = np.einsum("be,ce->bc",vbeq[:,:,i],vbeq[:,:,j].conj())
                expected = np.zeros((7,7))
                expected[i,j] = p0
                expected[2:,2:] += p1*np.outer(w[:,i],w[:,j])
                if i == j: expected[2:,2:] += p2*np.eye(5)/5
                worst["channel_matrix_units"] = max(worst["channel_matrix_units"],maxabs(raw-expected))
        expected_spectrum = sorted([p0/2]*2+[p1/2+p2/5]*2+[p2/5]*3)
        actual_spectrum = np.linalg.eigvalsh(rho_b)
        worst["spectrum"] = max(worst["spectrum"],maxabs(actual_spectrum-expected_spectrum))
        area_mean = float(np.trace(area@rho_b).real)
        expected_area = 4+2*(1-p0)
        inferred_p0 = (6-area_mean)/2
        inferred_c = 2*math.sqrt(max(0,inferred_p0))-1
        assert abs(inferred_c-c) < 1e-6
        worst["area_energy"] = max(worst["area_energy"],abs(area_mean-expected_area))
        work_mean = float(np.sum(np.abs(v.reshape(3,7,2,720,2)[:,:,1])**2))
        total_mean = float(np.vdot(v,hv).real/2)
        # work energy is 2 times probability, with sigma4 normalization 1/2.
        worst["full_resource_energy"] = max(worst["full_resource_energy"],
                                            abs(work_mean-2*p0),abs(total_mean-6),
                                            abs(total_mean-area_mean-work_mean))
        sb = entropy(rho_b)
        env_columns = vbeq.transpose(1,0,2).reshape(4320,14)/math.sqrt(2)
        se = entropy(env_columns.conj().T@env_columns)
        sbe = entropy(allmap.conj().T@allmap/2)
        mutual = sb+se-sbe
        assert abs((sb-math.log(2))-(mutual-se)) < 2e-12
        # Fixed-reference sector bookkeeping, not a new thermal-entropy identity.
        p0 = min(1.0,max(0.0,p0))
        log_ref_entropy = p0*math.log(2)+(1-p0)*math.log(5)
        omega = p0*sigma4+(1-p0)*sigma6
        omega_diag = np.diag(omega)
        cross = -sum(float(rho_b[i,i].real)*math.log(x)
                     for i,x in enumerate(omega_diag) if x>1e-14)
        relative_entropy = cross-sb
        worst["structure_account"] = max(worst["structure_account"],
             abs(sb-(binary_entropy(p0)+log_ref_entropy-relative_entropy)))
        fixed_reference=np.array([0.25]*2+[0.1]*5)
        fixed_relative_entropy=-sb-float(np.dot(np.diag(rho_b).real,np.log(fixed_reference)))
        worst["full_reference_balance"]=max(worst["full_reference_balance"],
            abs(fixed_relative_entropy-(relative_entropy+math.log(2)-binary_entropy(p0))),
            abs(sb-(log_ref_entropy+math.log(2)-fixed_relative_entropy)))
        row={"angle_over_pi":float(angle/math.pi),
             "program_probabilities":probs.tolist(),
             "area_in_a_units":area_mean,"boundary_energy_in_e_units":area_mean,
             "work_energy_in_e_units":work_mean,"S_B":sb,"S_full_environment":se,
             "S_B_environment":sbe,"I_B_environment":mutual,
             "sector_log_dimension_expectation":log_ref_entropy,
             "relative_entropy_to_same_weight_sector_reference":relative_entropy,
             "relative_entropy_to_fixed_full_H_B_marginal":fixed_relative_entropy,
             "spectrum":actual_spectrum.tolist()}
        rows.append(row)
        if abs(angle-math.pi/2)<1e-12:endpoints["middle"]=row
        if abs(angle-math.pi)<1e-12:endpoints["terminal"]=row
    assert max(worst.values()) < 3e-11
    gap=endpoints["middle"]["S_B"]-endpoints["terminal"]["S_B"]
    direct_gap=math.log(2**33/(3**12*5**5))/20
    assert abs(gap-direct_gap)<2e-12 and gap>1/15
    slope=math.log(F(1152,3125))/10
    assert slope<0
    # Same c at pi/2 and 3pi/2 has the same B state despite a changed full path.
    assert abs(rows[3]["S_B"]-rows[8]["S_B"])<1e-12
    return {
        "round":1274,"verdict":"PASS_INDEPENDENT_SCOPED","independent_groups":1,
        "author_code_imported":False,
        "method":"noncrossing valence-bond frames; all 720 physical permutations; full C-B-W-F path isometry",
        "physical_block_dimensions":[2,5],"full_output_dimension":3*7*2*720,
        "matrix_units_per_time":4,"time_points":len(rows),
        "max_errors":worst,"timeline":rows,"entropy_drop":gap,
        "entropy_drop_independent_log_formula":direct_gap,
        "dS_dc_at_c_zero":slope,
        "initial_entropy_change_as_c_decreases_is_positive":True,
        "exact_finite_bounds":exact_bounds(),
        "full_H_thermal_reference":full_thermal_reference(u1,u2,h0),
        "unknown_R_scope":"full input isometry and all matrix units; stated spectrum only when marginal input equals sigma4",
        "geometry_scope":"independent Casimir area expectation on this fixed trajectory, not a Lorentz-area instrument",
        "thermal_scope":"rules out only simultaneous raw S_B, all dH_B as heat, positive T, and zero correction",
        "actual_TR_temperature_or_heat_flow_certified":False,
        "same_complete_initial_state_counterexample_to_user_conjecture":False,
        "code_sha256":sha(Path(__file__))
    }

if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write","--write-once",dest="write",action="store_true")
    args=parser.parse_args()
    result=compute()
    target=HERE/"independent_results.json"
    if args.write:
        with target.open("x",encoding="utf-8",newline="\n") as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
        mode="FIRST_WRITE"
    else:
        if json.loads(target.read_text(encoding="utf-8")) != result:
            raise AssertionError("independent saved result differs")
        mode="READ_ONLY_MATCH"
    print(json.dumps({"round":1274,"verdict":result["verdict"],"mode":mode,
        "independent_groups":1,"max_diagnostic_error":max(result["max_errors"].values()),
        "entropy_drop":result["entropy_drop"],
        "robust_entropy_gap":result["exact_finite_bounds"]["entropy_drop_strict_lower_with_both_endpoint_errors"]
        },ensure_ascii=False))
