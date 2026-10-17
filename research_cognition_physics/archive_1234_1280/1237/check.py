"""1237: one finite SU(2) reference conversion with autonomous full output.

Default is read-only replay; --write creates results.json once.
The proof carries exact/all-code statements. NumPy checks this fixed carrier,
not a scan, a physical horizon, a thermal equation of state, or GR.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parents[1]
TOL = 2e-10


def need(ok, label):
    if not bool(ok):
        raise AssertionError(label)


def norm(x):
    return float(np.linalg.norm(x))


def entropy(eigenvalues):
    values = np.asarray(eigenvalues, dtype=float)
    need(float(values.min()) > -TOL, "entropy positive spectrum")
    values = values[values > 1e-12]
    return float(-np.sum(values * np.log(values)))


def kron_all(mats):
    out = np.ones((1, 1), dtype=complex)
    for mat in mats:
        out = np.kron(out, mat)
    return out


SPIN = [np.array([[0, 1], [1, 0]], dtype=complex) / 2,
        np.array([[0, -1j], [1j, 0]], dtype=complex) / 2,
        np.diag([.5, -.5]).astype(complex)]


def local_spins(dimensions):
    sites = []
    for site, dimension in enumerate(dimensions):
        components = []
        for generator in SPIN:
            local = np.zeros((dimension, dimension), dtype=complex)
            local[-2:, -2:] = generator
            factors = [np.eye(d) for d in dimensions]
            factors[site] = local
            components.append(kron_all(factors))
        sites.append(components)
    return sites


def singlet_basis(n):
    sites = local_spins([2] * n)
    totals = [sum(site[c] for site in sites) for c in range(3)]
    casimir = sum(j @ j for j in totals)
    eigenvalues, vectors = np.linalg.eigh(casimir)
    return vectors[:, np.abs(eigenvalues) < 1e-9]


def carrier():
    v4, v6 = singlet_basis(4), singlet_basis(6)
    need(v4.shape == (16, 2) and v6.shape == (64, 5), "singlet dimensions")
    singlet = np.array([0., 1., -1., 0.]) / np.sqrt(2)
    w = v6.conj().T @ np.kron(v4, singlet[:, None])
    embed4 = np.zeros((144, 16), complex)
    embed6 = np.zeros((144, 64), complex)
    for old in range(16):
        embed4[9 * old, old] = 1
        for x, y in itertools.product(range(2), repeat=2):
            embed6[9 * old + 3 * (x + 1) + y + 1, 4 * old + 2 * x + y] = 1
    basis = np.column_stack((embed4 @ v4, embed6 @ v6))
    sites = local_spins([2, 2, 2, 2, 3, 3])
    totals = [sum(site[c] for site in sites) for c in range(3)]
    total_casimir = sum(j @ j for j in totals)
    full_eigenvalues = np.linalg.eigvalsh(total_casimir)
    need(int(np.count_nonzero(np.abs(full_eigenvalues) < 1e-8)) == 7,
         "complete 144-dimensional Gauss support")
    area = np.zeros((144, 144), complex)
    for site in sites:
        casimir = sum(j @ j for j in site)
        vals, vecs = np.linalg.eigh(casimir)
        area += (vecs * np.sqrt(np.maximum(vals, 0))) @ vecs.conj().T
    # a=1, a_star=2/sqrt(3); no Newton or Planck-scale input.
    area *= 2 / np.sqrt(3)
    area7 = basis.conj().T @ area @ basis
    q56 = basis[:, 2:].conj().T @ sum(sites[4][c] @ sites[5][c]
                                     for c in range(3)) @ basis[:, 2:]
    perms = list(itertools.permutations(range(6)))
    perm6 = []
    for perm in perms:
        indices = np.arange(64).reshape((2,) * 6).transpose(perm).reshape(64)
        perm6.append(v6.conj().T @ v6[indices, :])
    perm6 = np.asarray(perm6)
    generators = []
    for i in range(5):
        p = list(range(6))
        p[i], p[i + 1] = p[i + 1], p[i]
        generators.append(perm6[perms.index(tuple(p))])
    # Explicit all-matrix commutant map X -> ([X,U_i])_i.
    columns = []
    for k in range(25):
        matrix = np.zeros((5, 5), complex)
        matrix.flat[k] = 1
        columns.append(np.concatenate([(matrix @ u - u @ matrix).ravel()
                                       for u in generators]))
    comm_singular = np.linalg.svd(np.column_stack(columns), compute_uv=False)
    comm_rank = int(np.count_nonzero(comm_singular > 1e-9))
    need(comm_rank == 24, "full S6 has scalar commutant")
    tau = w @ w.conj().T / 2
    sigma6 = np.eye(5) / 5
    shape_tau = float(np.trace(tau @ q56).real)
    shape_sigma = float(np.trace(sigma6 @ q56).real)
    singlet56 = np.eye(5) / 4 - q56
    weak = []
    for p4 in itertools.permutations(range(4)):
        for p2 in ((4, 5), (5, 4)):
            weak.append(perm6[perms.index(p4 + p2)])
    weak_error = max(norm(u @ tau @ u.conj().T - tau) for u in weak)
    residuals = dict(orthonormal=norm(basis.conj().T @ basis - np.eye(7)),
                     Gauss=max(norm(j @ basis) for j in totals),
                     area=norm(area7 - np.diag([4., 4., 6., 6., 6., 6., 6.])),
                     W_isometry=norm(w.conj().T @ w - np.eye(2)),
                     permutation_unitarity=max(norm(u.conj().T @ u - np.eye(5))
                                               for u in perm6),
                     weak_reference_invariance=weak_error,
                     q_tau=abs(shape_tau + .75), q_sigma=abs(shape_sigma + .15))
    need(max(residuals.values()) < TOL, "group 1 operator identities")
    half_trace = float(np.sum(np.abs(np.linalg.eigvalsh(tau - sigma6))) / 2)
    need(abs(half_trace - .6) < TOL, "same area/energy shape deletion")
    result = dict(status="PASS",kinematic_dimension=144,physical_dimension=7,
                  sector_dimensions=[2,5], area_eigenvalues=[4,6],a_star=2/np.sqrt(3),
                  permutation_count=720,commutant_rank=comm_rank,
                  commutant_nullity=25-comm_rank,
                  minimum_nonzero_commutant_singular=float(comm_singular[-2]),
                  weak_group_size=len(weak), residuals=residuals,
                  q56_before=shape_tau,q56_full_reference=shape_sigma,
                  singlet56_before=float(np.trace(tau @ singlet56).real),
                  singlet56_full_reference=float(np.trace(sigma6 @ singlet56).real),
                  same_area_energy_half_trace_distance=half_trace,
                  reference_role_is_new_input=True)
    return result, w, perm6, tau


def autonomous(w, perm6):
    hb = np.array([4., 4., 6., 6., 6., 6., 6.])
    hw = np.array([0., 2.])
    bare = (hb[:, None] + hw).ravel()
    source, target = np.zeros((14, 2), complex), np.zeros((14, 2), complex)
    source[1, 0], source[3, 1] = 1, 1
    target[4::2] = w
    u1 = (np.eye(14) - source @ source.conj().T - target @ target.conj().T
          + target @ source.conj().T + source @ target.conj().T)
    u2 = np.zeros((720, 7, 7), complex)
    u2[:, :2, :2] = np.eye(2)
    u2[:, 2:, 2:] = perm6

    def gate1(x):
        return np.einsum("ij,jfk->ifk",u1,x.reshape(14,720,-1),optimize=True).reshape(x.shape)

    def gate2(x, adjoint=False):
        operators = u2.conj().transpose(0,2,1) if adjoint else u2
        return np.einsum("fij,jwfk->iwfk",operators,x,optimize=True)

    def program_action(x):
        out = np.zeros_like(x)
        out[0] = gate1(x[1])
        out[1] = gate1(x[0]) + gate2(x[2],True)
        out[2] = gate2(x[1])
        return out  # g=1 in this single fixed case.

    initial = np.zeros((7,2,720,2),complex)
    initial[0,1,:,0] = 1/np.sqrt(720)
    initial[1,1,:,1] = 1/np.sqrt(720)
    after1 = gate1(initial)
    after2 = gate2(after1)
    recovered = gate1(gate2(after2,True))
    initial_free = bare.reshape(7,2,1,1) * initial
    need(norm(initial_free - 6*initial)<TOL, "fixed total bare source energy")
    residuals = dict(gate1_unitarity=norm(u1.conj().T@u1-np.eye(14)),
                     gate1_energy=norm(u1*bare[None,:]-bare[:,None]*u1),
                     gate1_pairactivation=norm(u1@source-target),
                     gate2_unitarity=max(norm(u.conj().T@u-np.eye(7)) for u in u2),
                     gate2_energy=max(norm(u*hb[None,:]-hb[:,None]*u) for u in u2),
                     full_code_recovery=norm(recovered-initial),
                     full_code_isometry=norm(after2.reshape(-1,2).conj().T@
                                             after2.reshape(-1,2)-np.eye(2)))
    path = np.array([[0.,1.,0.],[1.,0.,1.],[0.,1.,0.]])
    ev, vec = np.linalg.eigh(path)
    tstar = np.pi/np.sqrt(2)
    states = [initial, after1, after2]
    samples = []
    for fraction in (0.,.5,1.):
        t = fraction*tstar
        amplitudes = vec@(np.exp(-1j*ev*t)*vec[0,:].conj())
        psi = np.stack([amplitudes[k]*states[k] for k in range(3)]) * np.exp(-6j*t)
        interaction = program_action(psi)
        free = bare.reshape(1,7,2,1,1)*psi
        energy = np.vdot(psi,free+interaction).real/2
        eint = np.vdot(psi,interaction).real/2
        eint2 = np.vdot(interaction,interaction).real/2
        area_now = np.vdot(psi,hb.reshape(1,7,1,1,1)*psi).real/2
        work_now = np.vdot(psi,hw.reshape(1,1,2,1,1)*psi).real/2
        p0 = ((1+np.cos(np.sqrt(2)*t))/2)**2
        residuals["energy_%s"%fraction] = abs(float(energy)-6)
        residuals["interaction_mean_%s"%fraction] = abs(float(eint))
        residuals["interaction_second_%s"%fraction] = abs(float(eint2)-1)
        residuals["dynamic_area_%s"%fraction] = abs(float(area_now)-(4+2*(1-p0)))
        residuals["dynamic_work_%s"%fraction] = abs(float(work_now)-2*p0)
        # Checks the exact path dressing action on all two input columns.
        expected = np.stack([(path@amplitudes)[k]*states[k] for k in range(3)])*np.exp(-6j*t)
        residuals["path_dressing_%s"%fraction] = norm(interaction-expected)
        samples.append(dict(t_over_tstar=fraction,total_energy=float(energy),
                            interaction_mean=float(eint),interaction_second_moment=float(eint2),
                            area=float(area_now),work_energy=float(work_now),start_probability=float(p0)))
        if fraction == 1:
            residuals["end_to_end"] = norm(psi[2]+np.exp(-6j*t)*after2)+norm(psi[:2])
    return_amplitudes = vec@(np.exp(-1j*ev*2*tstar)*vec[0,:].conj())
    residuals["forward_return"] = norm(return_amplitudes-np.array([1.,0.,0.]))
    delta_t=1e-4
    shifted=vec@(np.exp(-1j*ev*(tstar+delta_t))*vec[0,:].conj())
    timing_error=norm(shifted-np.array([0.,0.,-1.]))
    timing_bound=np.sqrt(2)*delta_t
    need(timing_error<=timing_bound+TOL,"finite terminal time error bound")
    need(max(residuals.values())<TOL,"group 2 all-code autonomous identities")
    min_total = float(bare.min()+ev.min())
    need(min_total>0, "positive total Hamiltonian by unitary path dressing")
    before_B = np.vdot(initial,hb.reshape(7,1,1,1)*initial).real/2
    after_B = np.vdot(after2,hb.reshape(7,1,1,1)*after2).real/2
    before_W = np.vdot(initial,hw.reshape(1,2,1,1)*initial).real/2
    after_W = np.vdot(after2,hw.reshape(1,2,1,1)*after2).real/2
    result = dict(status="PASS",e=1.,g=1.,E4=4.,E6=6.,work_gap=2.,
                  effective_data_dimension=7*2*720,
                  effective_program_dimension=3*7*2*720,
                  kinematic_program_dimension=3*144*2*720,
                  tstar=float(tstar),total_positive_lower_bound=min_total,
                  boundary_energy_before=float(before_B),boundary_energy_after=float(after_B),
                  work_energy_before=float(before_W),work_energy_after=float(after_W),
                  interaction_not_discarded=True,samples=samples,residuals=residuals,
                  terminal_time_error=delta_t,phase_removed_isometry_error=timing_error,
                  timing_error_bound=float(timing_bound),forward_return_time=float(2*tstar),
                  permanently_stable_reference=False,
                  free_pulses_used=False,physical_Lorentz_clock_supplied=False)
    return result,after2[2:,0,:,:]


def information(final, tau):
    sigma = np.eye(5)/5
    residuals = {}
    for i,j in itertools.product(range(2),repeat=2):
        out = final[:,:,i]@final[:,:,j].conj().T
        target = sigma if i==j else np.zeros((5,5))
        residuals["matrix_unit_%d%d"%(i,j)] = norm(out-target)
    reference_B = sum(final[:,:,i]@final[:,:,i].conj().T for i in range(2))/2
    # R purifies the initial I2/2. This 10x720 matrix has precisely the
    # nonzero flag spectrum and BR spectrum; no huge density matrix needed.
    br_amplitudes = final.transpose(0,2,1).reshape(10,720)/np.sqrt(2)
    rho_BR = br_amplitudes@br_amplitudes.conj().T
    flag_spectrum = np.linalg.eigvalsh(rho_BR)
    residuals["reference_B"] = norm(reference_B-sigma)
    residuals["BR_product"] = norm(rho_BR-np.eye(10)/10)
    residuals["flag_spectrum"] = norm(flag_spectrum-np.ones(10)/10)
    sb = entropy(np.linalg.eigvalsh(reference_B))
    sf = entropy(flag_spectrum)
    # For sigma4 input and pure F, the BF spectrum is unchanged from sigma4.
    whole_columns = final.reshape(-1,2)/np.sqrt(2)
    sbf = entropy(np.linalg.eigvalsh(whole_columns.conj().T@whole_columns))
    initial_entropy = np.log(2.)
    mutual = sb+sf-sbf
    residuals["entropy_balance"] = abs(sb-initial_entropy-mutual+sf)
    residuals["flag_entropy_exact"] = abs(sf-np.log(10))
    residuals["mutual_information_exact"] = abs(mutual-2*np.log(5))
    # Sharp boundary-only reversible comparison: rho -> rho tensor I2/2,
    # embedded in four of five dimensions. It preserves the entire input.
    best_reversible_reference = np.diag([.25,.25,.25,.25,0.])
    rank_distance = float(np.sum(np.abs(np.linalg.eigvalsh(best_reversible_reference-sigma)))/2)
    need(abs(rank_distance-.2)<TOL,"sharp exact-recovery rank bound")
    # Its concrete CPTP decoder traces the auxiliary factor on support4,
    # and sends unused basis5 to |0><0|. Verify all four matrix units.
    decoder = []
    for auxiliary in range(2):
        k = np.zeros((2,5),complex)
        for logical in range(2):
            k[logical,2*logical+auxiliary]=1
        decoder.append(k)
    last=np.zeros((2,5),complex); last[0,4]=1; decoder.append(last)
    residuals["decoder_TP"] = norm(sum(k.conj().T@k for k in decoder)-np.eye(5))
    for i,j in itertools.product(range(2),repeat=2):
        unit=np.zeros((2,2),complex);unit[i,j]=1
        encoded=np.zeros((5,5),complex);encoded[:4,:4]=np.kron(unit,np.eye(2)/2)
        decoded=sum(k@encoded@k.conj().T for k in decoder)
        residuals["reversible_%d%d"%(i,j)]=norm(decoded-unit)
    need(max(residuals.values())<TOL,"group 3 exact channel and entropy identities")
    return dict(status="PASS",flag_initial="pure_coherent_uniform_720",
                boundary_entropy_before=float(initial_entropy),boundary_entropy_after=sb,
                boundary_entropy_change=sb-initial_entropy,
                area_change=2.,flag_nonzero_eigenvalues=flag_spectrum.tolist(),
                flag_rank=10,flag_entropy=sf,BF_entropy=sbf,BF_mutual_information=mutual,
                reference_R_dimension=2,BR_product=True,
                stopped_twirl_entropy=entropy(np.linalg.eigvalsh(tau)),
                boundary_only_reversible_rank_max=4,
                boundary_only_reference_distance_lower_bound=.2,
                bound_attained_by_reversible_channel=rank_distance,
                residuals=residuals)


def calculate():
    one,w,perms,tau=carrier()
    two,final=autonomous(w,perms)
    three=information(final,tau)
    return dict(schema="research_round_1237_v1",status="PASS",science_groups=3,
                group_1=one,group_2=two,group_3=three,
                scope=dict(new_reference_role_adopted=False,
                           specified_reference_conversion=True,
                           autonomy_in_finite_effective_Hamiltonian=True,
                           all_unknown_input_R_retained_in_full_output=True,
                           physical_spacetime_area=False,
                           complete_Htotal_Gibbs_reference=False,
                           thermal_temperature_selected=False,
                           geometric_thermal_balance=False,Clausius=False,
                           Einstein_dynamics=False,complete_GR=False))


def compare(current,saved,path="result"):
    if isinstance(current,dict):
        need(isinstance(saved,dict) and current.keys()==saved.keys(),path+" keys")
        for key in current:compare(current[key],saved[key],path+"."+key)
    elif isinstance(current,list):
        need(isinstance(saved,list) and len(current)==len(saved),path+" length")
        for i,(a,b) in enumerate(zip(current,saved)):compare(a,b,path+"[%s]"%i)
    elif isinstance(current,float):
        need(isinstance(saved,(int,float)) and np.isclose(current,saved,rtol=1e-10,atol=1e-10),
             path+" float replay")
    else:
        need(current==saved,path+" exact replay")


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--write",action="store_true",help="create results.json only if absent")
    args=parser.parse_args()
    result=calculate()
    paths=[HERE/name for name in ("admission.md","proof.md","check.py","NEXT.md")]
    paths += [HERE.parent/"_admission/quantum_boundary_carrier_screen.md"]
    result["source_sha256"]={p.relative_to(RESEARCH).as_posix():
                            hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    output=HERE/"results.json"
    if args.write:
        with output.open("x",encoding="utf-8",newline="\n") as stream:
            json.dump(result,stream,ensure_ascii=False,indent=2)
            stream.write("\n")
        mode="CREATED"
    else:
        compare(result,json.loads(output.read_text(encoding="utf-8")))
        mode="READ_ONLY_MATCH"
    print(json.dumps(dict(status="PASS",mode=mode,science_groups=3,
                          delta_area=2.,delta_entropy=result["group_3"]["boundary_entropy_change"],
                          flag_rank=10,physical_spacetime_area=False),ensure_ascii=False))


if __name__=="__main__":
    main()
