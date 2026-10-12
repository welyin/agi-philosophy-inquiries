"""Read-only finite diagnostic for after1140 recoil qualification.
Not a simulation of the original Dirac sea, g, Y, or the full collection network.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve()
RESULT = HERE.with_name("after1140_recoil_results.json")

def norm(x):
    return float(np.linalg.norm(x))

def expect(psi, a):
    return float(np.vdot(psi, a @ psi).real)

def compare(a, b, path="root"):
    if isinstance(a, dict):
        assert set(a) == set(b), (path, "keys")
        for k in a:
            compare(a[k], b[k], path+"."+k)
    elif isinstance(a, list):
        assert len(a) == len(b), (path, "length")
        for k, (x, y) in enumerate(zip(a, b)):
            compare(x, y, path+f"[{k}]")
    elif isinstance(a, float):
        assert np.isclose(a, b, rtol=1e-8, atol=2e-10), (path, a, b)
    else:
        assert a == b, (path, a, b)

def run():
    sites, rails, modes = 3, 2, 6
    nf = 2**modes
    eye_f = np.eye(nf, dtype=complex)
    creators = []
    for j in range(modes):
        c = np.zeros((nf, nf), dtype=complex)
        for n in range(nf):
            if not ((n >> j) & 1):
                c[n | (1 << j), n] = (-1)**((n & ((1 << j)-1)).bit_count())
        creators.append(c)
    annih = [c.conj().T for c in creators]
    bilinear = [[creators[j] @ annih[k] for k in range(modes)] for j in range(modes)]
    def q(a):
        return sum((a[j,k]*bilinear[j][k] for j in range(modes) for k in range(modes)),
                   np.zeros((nf, nf), dtype=complex))
    car_error = max(norm(annih[j] @ creators[k] + creators[k] @ annih[j] -
                         (eye_f if j == k else 0))
                    for j in range(modes) for k in range(modes))
    shift = np.roll(np.eye(sites), 1, axis=0).astype(complex)
    t1 = np.kron(shift, np.eye(rails))
    phi = .37
    hs = np.eye(sites) - 1.2*(np.exp(1j*phi)*shift + np.exp(-1j*phi)*shift.conj().T)
    h = np.kron(hs, np.eye(rails))
    eps, basis = np.linalg.eigh(h)
    sea = basis[:, eps < 0]
    S = sea @ sea.conj().T
    omega = np.zeros(nf, complex)
    omega[0] = 1.
    for v in sea.T:
        cv = sum((v[j]*creators[j] for j in range(modes)), np.zeros((nf,nf),complex))
        omega = cv @ omega
    omega /= norm(omega)
    hf = q(h) - float(np.trace(S @ h).real)*eye_f
    vacuum_error = norm(hf @ omega)
    p = np.diag([1., 1., 1., 1., 0., 0.])
    a_hold = h @ p + p @ h - p @ h @ p
    bs = np.zeros((sites,sites),complex)
    bs[0,1], bs[1,0] = .4j, -.4j
    b = np.kron(bs, np.eye(rails))
    energies_y = 2*np.pi*np.fft.fftfreq(sites)
    fourier = np.exp(2j*np.pi*np.arange(sites)[:,None]*np.arange(sites)[None,:]/sites)/np.sqrt(sites)
    py = fourier @ np.diag(energies_y) @ fourier.conj().T
    ty = py @ py / 4
    pc = py.copy()
    hc = pc @ pc / 3
    chi = np.diag([.2, 1., .4])
    ey = np.eye(sites)/2 + (shift+shift.conj().T)/4
    E = np.kron(np.kron(ey, eye_f), np.eye(3))
    H0 = (np.kron(np.kron(np.eye(sites), hf), np.eye(3)) +
          np.kron(np.kron(ty, eye_f), np.eye(3)) +
          np.kron(np.kron(np.eye(sites), eye_f), hc))
    V = np.zeros_like(H0)
    holds, gates = [], []
    for y in range(sites):
        tt = np.linalg.matrix_power(t1, y)
        ay, by = tt @ a_hold @ tt.conj().T, tt @ b @ tt.conj().T
        Q, G = q(ay), q(by)
        holds.append(Q)
        gates.append(G)
        yy = np.diag([float(j == y) for j in range(sites)])
        H0 -= np.kron(np.kron(yy, Q), np.eye(3))
        V += np.kron(np.kron(yy, G), chi)
    psi_y = np.exp(-np.array([0.,1.,1.])/4).astype(complex)
    psi_y /= norm(psi_y)
    clock = np.array([.6, 1., .5],complex)*np.exp(1j*np.array([.1,.4,-.2]))
    clock /= norm(clock)
    psi = np.kron(np.kron(psi_y, omega), clock)
    # Actual finite Fock lift of the cyclic one-particle permutation.
    uf = np.zeros((nf,nf),complex)
    for n in range(nf):
        vector = np.zeros(nf,complex)
        vector[0] = 1.
        occupied = [j for j in range(modes) if (n >> j) & 1]
        for j in reversed(occupied):
            v = t1[:,j]
            cv = sum((v[k]*creators[k] for k in range(modes)),
                     np.zeros((nf,nf),complex))
            vector = cv @ vector
        uf[:,n] = vector
    joint_shift = np.kron(np.kron(shift, uf), np.eye(3))
    covariance_error = max(norm(joint_shift @ op - op @ joint_shift) for op in (H0,V))
    chi2 = expect(clock, chi @ chi)
    overlap = float(np.vdot(psi_y, shift @ psi_y).real)
    dsea = (gates[0] - gates[1]) @ omega
    predicted = -overlap*chi2*norm(dsea)**2/4
    vp, vvp = V @ psi, V @ (V @ psi)
    comm_coefficient = expect(vp, E)-float(np.vdot(psi, E @ vvp).real)
    rows, second, third_bounds, eigens = [], {}, {}, {}
    for z in (-1,0,1):
        H = H0 + z*V
        hp, h2p = H @ psi, H @ (H @ psi)
        h3p = H @ h2p
        m2 = expect(hp,E)-float(np.vdot(psi,E @ h2p).real)
        bound = norm(h3p)/3 + norm(hp)*norm(h2p)
        vals, vecs = np.linalg.eigh(H)
        coeff = vecs.conj().T @ psi
        eigens[z] = (vals,vecs,coeff)
        second[z], third_bounds[z] = m2,bound
        rows.append(dict(program=z, second_taylor_coefficient=m2,
                         H_norm=norm(hp), H2_norm=norm(h2p), H3_norm=norm(h3p),
                         cubic_remainder_bound=bound,
                         selfadjoint_error=norm(H-H.conj().T)))
    direct = .5*(second[1]+second[-1])-second[0]
    full_bound = .5*(third_bounds[1]+third_bounds[-1])+third_bounds[0]
    # Fixed auxiliary background offset: a nonzero mixed term must cancel too.
    # This is an algebra diagnostic, not a change of the actual continuum H.
    offset_second = {}
    offset = 1/3
    for z in (-1,0,1):
        hz = H0 + (offset+z)*V
        hp = hz @ psi
        h2p = hz @ hp
        offset_second[z] = expect(hp,E)-float(np.vdot(psi,E @ h2p).real)
    offset_contrast = .5*(offset_second[1]+offset_second[-1])-offset_second[0]
    offset_mixed = .5*(offset_second[1]-offset_second[-1])
    assert abs(offset_mixed) > 1e-5
    assert abs(offset_contrast-predicted) < 1e-11
    assert abs(offset_mixed-2*offset*predicted) < 1e-11
    tests = []
    for t in (.05,.025):
        probabilities = {}
        for z in (-1,0,1):
            vals,vecs,coeff = eigens[z]
            evolved = vecs @ (np.exp(-1j*t*vals)*coeff)
            probabilities[z] = expect(evolved,E)
        contrast = .5*(probabilities[1]+probabilities[-1])-probabilities[0]
        remainder = abs(contrast-predicted*t*t)
        bound = full_bound*t**3
        assert remainder <= bound + 1e-12
        tests.append(dict(time=t, p_minus=probabilities[-1], p_zero=probabilities[0],
                          p_plus=probabilities[1], contrast=contrast,
                          quadratic_prediction=predicted*t*t,
                          absolute_remainder=remainder, certified_cubic_bound=bound,
                          within_bound=True))
    assert car_error < 1e-12 and vacuum_error < 1e-12
    assert covariance_error < 1e-10
    assert predicted < -1e-5
    assert max(abs(direct-predicted),abs(comm_coefficient-predicted)) < 1e-11
    return dict(
        status="PASS",
        scope="Finite six-mode CAR, three-site periodic support and three-state clock only; not original Dirac sea/g/Y/17m network.",
        formal_round_opened=False,
        original_fixed_input_passed=False,
        original_fixed_precision_passed=False,
        new_calibration_and_terminal_read_permission=True,
        dimension=int(H0.shape[0]),
        car_error=car_error,
        normal_sea_occupied_modes=int(sea.shape[1]),
        normal_sea_free_vacuum_error=vacuum_error,
        common_cyclic_translation_error=covariance_error,
        cyclic_body_overlap=overlap,
        clock_chi_squared=chi2,
        sea_translation_difference_norm_squared=norm(dsea)**2,
        predicted_quadratic_coefficient=predicted,
        direct_symmetric_H_coefficient=direct,
        double_commutator_coefficient=comm_coefficient,
        background_second_coefficient=second[0],
        linear_H0_V_coefficient=(second[1]-second[-1])/2,
        auxiliary_mixed_term_diagnostic=dict(background_offset=offset,
            nonzero_linear_mixed_coefficient=offset_mixed,
            predicted_linear_mixed_coefficient=2*offset*predicted,
            symmetric_coefficient_after_cancellation=offset_contrast,
            not_an_actual_continuum_model_change=True),
        coefficient_agreement=max(abs(direct-predicted),abs(comm_coefficient-predicted)),
        vector_domain_diagnostics=rows,
        exact_finite_evolution=tests)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    data=run()
    if args.write_results:
        RESULT.write_text(json.dumps(data,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
        print("PASS: wrote new qualification results; no formal round opened.")
    else:
        compare(data,json.loads(RESULT.read_text(encoding="utf-8")))
        print("PASS: read-only replay matches saved qualification results.")
    print(json.dumps({k:data[k] for k in ("dimension","predicted_quadratic_coefficient",
          "direct_symmetric_H_coefficient","linear_H0_V_coefficient","coefficient_agreement")},
          ensure_ascii=False))
if __name__=="__main__":
    main()
