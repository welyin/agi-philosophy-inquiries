"""Independent small-matrix and native-instrument audit for round 1210.

This does not import the author's round-1210 code. Floating computations check
the algebra; the all-size classification is inherited from round 1209.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import numpy as np

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def gibbs(e, beta):
    x = -beta * np.asarray(e)
    x -= x.max()
    p = np.exp(x)
    return p / p.sum()


def fit_beta(e, p):
    e, y = np.asarray(e), np.log(p)
    ec, yc = e - e.mean(), y - y.mean()
    if np.linalg.norm(ec) < 1e-13:
        return None, float(np.max(np.abs(yc)))
    beta = -np.dot(ec, yc) / np.dot(ec, ec)
    return float(beta), float(np.max(np.abs(yc + beta * ec)))


def read(p):
    return json.loads(p.read_text(encoding='utf-8'))


def evaluate():
    h = np.array([[0, -.2, 0, -np.sqrt(3)/80],
                  [-.2, 1, 0, 0], [0, 0, 1, 0],
                  [-np.sqrt(3)/80, 0, 0, -.025]])
    q = np.zeros((4, 4)); q[1, 2] = q[2, 1] = 1
    e, r = np.linalg.eigh(h)
    specs = [np.array([0., 1.]), np.array([-.2, .1, 1.4]), e]
    cases = []
    for ea in specs:
        for eb in specs:
            for beta in [-1.3, 0., .4, 1.2]:
                energy = (ea[:, None] + eb).ravel()
                p = np.kron(gibbs(ea, beta), gibbs(eb, beta))
                b, residual = fit_beta(energy, p)
                assert abs(b-beta) < 1e-12 and residual < 1e-12
                variance = np.dot(gibbs(ea, beta), ea**2) - np.dot(gibbs(ea, beta), ea)**2
                assert variance > 0
                cases.append(dict(dim_a=len(ea), dim_b=len(eb), beta=beta,
                                  fitted_beta=b, log_residual=residual))
    # Non-identifying middle object: it belongs to every Gibbs-beta family.
    ea = np.array([0., 1.]); ez = np.zeros(2)
    a, c, z = gibbs(ea, np.log(2)), gibbs(ea, np.log(4)), np.ones(2)/2
    ab = fit_beta((ea[:, None]+ez).ravel(), np.kron(a, z))
    bc = fit_beta((ez[:, None]+ea).ravel(), np.kron(z, c))
    ac = fit_beta((ea[:, None]+ea).ravel(), np.kron(a, c))
    assert ab[1] < 1e-12 and bc[1] < 1e-12 and ac[1] > .3

    # A non-scalar middle object with partner-dependent, disjoint ports.
    hot, cold = gibbs(e, 1/3), gibbs(e, 6/5)
    left = np.einsum('i,j,k->ijk', cold, cold, hot)
    right = np.einsum('i,j,k->ijk', cold, hot, hot)
    assert np.max(abs(left-left.swapaxes(0, 1))) < 1e-15
    assert np.max(abs(right-right.swapaxes(1, 2))) < 1e-15
    ac_product = np.outer(cold, hot)
    td = float(np.abs(ac_product-ac_product.T).sum()/2)
    assert td > .18
    beta_middle, middle_residual = fit_beta((e[:, None]+e).ravel(), np.kron(cold, hot))
    assert middle_residual > .4

    # Reconstruct the actual finite QQ thermometry map by independent eigh.
    ident = np.eye(4)
    hh = np.kron(h, ident)+np.kron(ident, h)+.1*np.kron(q, q)
    eigen, v = np.linalg.eigh(hh)
    u = (v*np.exp(-600j*eigen)) @ v.conj().T
    basis = np.kron(r, r)
    ue = basis.T @ u @ basis
    matrix = np.empty((4, 4))
    for b in range(4):
        matrix[:, b] = np.abs(ue[:, 8+b].reshape(4, 4))**2 @ np.ones(4)
    saved = read(STAGE/'1201/results.json')
    contact = next(c for c in saved['contacts'].values() if c['t'] == 600)
    old_matrix = np.array(contact['K'], float)
    matrix_error = float(np.max(abs(matrix-old_matrix)))
    inverse_norm = float(np.linalg.norm(np.linalg.inv(matrix), np.inf))
    assert matrix_error < 1e-10 and inverse_norm < 3.14
    eta = .015701
    # Any data within eta of the true population produce the 2*eta envelope.
    envelopes = []
    width = e[-1]-e[0]
    for beta in [1/3, 6/5]:
        p = gibbs(e, beta)
        assert min(p) > 2*eta
        lo = np.log((p[0]-2*eta)/(p[-1]+2*eta))/width
        hi = np.log((p[0]+2*eta)/(p[-1]-2*eta))/width
        envelopes.append(dict(beta=beta, guaranteed_interval_envelope=[float(lo), float(hi)]))
        # Random normalized populations emulate data inside the same box.
        rng = np.random.default_rng(int(beta*3000))
        for _ in range(100):
            d = rng.uniform(-1, 1, 4); d -= d.mean()
            d *= eta / max(1., np.max(abs(d)))
            phat = p+d
            inner_lo = np.log((phat[0]-eta)/(phat[-1]+eta))/width
            inner_hi = np.log((phat[0]+eta)/(phat[-1]-eta))/width
            assert lo-1e-12 <= inner_lo <= beta <= inner_hi <= hi+1e-12
    assert envelopes[0]['guaranteed_interval_envelope'][1] < envelopes[1]['guaranteed_interval_envelope'][0]
    # A fixed pairwise tolerance does not produce transitivity.
    bs, tolerance = [0., .75, 1.5], 1.
    assert abs(bs[0]-bs[1]) <= tolerance and abs(bs[1]-bs[2]) <= tolerance
    assert abs(bs[0]-bs[2]) > tolerance
    # Energy origins do not affect beta. A scale change does.
    bshift, errshift = fit_beta(e+9, cold)
    bscale, errscale = fit_beta(3*e, cold)
    assert abs(bshift-1.2) < 1e-12 and abs(bscale-.4) < 1e-12
    sources = ['1201/results.json', '1201/acceptance.json', '1209/acceptance.json']
    return dict(status='PASS', all_checks_passed=True,
        source_sha256={x: sha(STAGE/x) for x in sources},
        independent_code_sha256=sha(Path(__file__)),
        shared_beta_cases=cases, scalar_middle=dict(AB=list(ab), BC=list(bc), AC=list(ac)),
        disjoint_port_example=dict(pair_AB_invariant=True, pair_BC_invariant=True,
            AC_joint_trace_distance=td, middle_not_single_beta_residual=middle_residual),
        native_readout=dict(matrix_max_error=matrix_error, inverse_inf_norm=inverse_norm,
            population_error=eta, thermal_family_interval_envelopes=envelopes),
        qualifications=dict(native_QQ_all_scale_permission_proved=False,
            finite_data_certify_exact_Gibbs_role=False, exact_beta_equality_transitive=True,
            fixed_epsilon_relation_transitive=False, negative_beta_not_excluded=True,
            new_axiom_adopted=False, stage_complete=False))


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--write', action='store_true')
    args = parser.parse_args(); result = evaluate(); output = HERE/'independent_results.json'
    if args.write:
        assert not (HERE/'acceptance.json').exists(), 'Do not overwrite accepted results'
        output.write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    else:
        assert read(output) == result, 'Saved independent result differs'
    print(json.dumps(dict(status='PASS', all_checks_passed=True,
        shared_beta_cases=len(result['shared_beta_cases']), results_match=True)))


if __name__ == '__main__':
    main()
