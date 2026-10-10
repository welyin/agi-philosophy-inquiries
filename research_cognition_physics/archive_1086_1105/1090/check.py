"""Round 1090: finite quantum calibration of one compatible typed family.

Default: recompute and compare to results.json without writing.
--write: create results.json once; refuse to replace any existing result.
No finite computation certifies the infinite quantifiers or grants primitive rights.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np

SEED = 109020261010
TOL = 2e-11
I2 = np.eye(2, dtype=complex)
SIGMA = np.array([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]]], dtype=complex)
ROOT = Path(__file__).resolve().parent


def norm(a):
    return float(np.linalg.norm(a, ord=2))


def maxabs(a):
    return float(np.max(np.abs(a))) if np.size(a) else 0.0


def sqrt_psd(a):
    vals, vecs = np.linalg.eigh((a + a.conj().T) / 2)
    if vals.min() < -TOL:
        raise AssertionError('nonpositive matrix')
    return (vecs * np.sqrt(np.maximum(vals, 0))) @ vecs.conj().T


def density(rng, d):
    a = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    a = a @ a.conj().T
    return a / np.trace(a)


def partial_second(rho, da, db):
    return np.trace(rho.reshape(da, db, da, db), axis1=0, axis2=2)


def effect(z):
    z = np.asarray(z, dtype=float)
    f = z / np.sqrt(1 + np.dot(z, z))
    return (I2 + 0.5 * np.einsum('i,ijk->jk', f, SIGMA)) / 2


def tent(z):
    return max(0.0, 1.0 - float(np.linalg.norm(z)))


def keys(points):
    return [tuple(np.round(np.asarray(x, dtype=float), 12)) for x in points]


def union_points(groups):
    return [np.array(k) for k in sorted(set(k for g in groups for k in keys(g)))]


def embedding(small, large):
    old, new = keys(small), keys(large)
    if len(set(old)) != len(old) or len(set(new)) != len(new):
        raise AssertionError('directory has duplicates')
    j = np.zeros((len(new), len(old)), dtype=complex)
    for col, x in enumerate(old):
        j[new.index(x), col] = 1
    return j


def index_embedding(indices, target_size):
    j = np.zeros((target_size, len(indices)), dtype=complex)
    for k, idx in enumerate(indices):
        j[idx, k] = 1
    return j


def block_diag(blocks):
    n = sum(b.shape[0] for b in blocks)
    ans = np.zeros((n, n), dtype=complex)
    at = 0
    for b in blocks:
        d = b.shape[0]
        ans[at:at+d, at:at+d] = b
        at += d
    return ans


def kraus_catalogue(s, anchors, frames, reader):
    """Order: coherent frame, coherent anchor, coherent endpoint, payload."""
    effects = []
    for r in frames:
        for a in anchors:
            for x in s:
                z = r.T @ (np.asarray(x) - np.asarray(a))
                effects.append(tent(z) * I2 if reader == 'tent' else effect(z))
    return [block_diag([sqrt_psd(e) for e in effects]),
            block_diag([sqrt_psd(I2 - e) for e in effects])]


def instrument_blocks(ks, rho, refdim):
    return [(np.kron(k, np.eye(refdim)) @ rho @ np.kron(k.conj().T, np.eye(refdim))) for k in ks]


def controlled_transport(s, maps):
    images = [[fn(np.asarray(x)) for x in s] for fn in maps]
    out = union_points(images)
    n, m, p = len(s), len(out), len(maps)
    v = np.zeros((p*m, p*n), dtype=complex)
    outkeys = keys(out)
    for j, img in enumerate(images):
        for k, y in enumerate(keys(img)):
            v[j*m + outkeys.index(y), j*n + k] = 1
    return v, out


def run():
    rng = np.random.default_rng(SEED)
    cases = []

    def check(name, residual, count=1, detail=None):
        residual = float(residual)
        if not np.isfinite(residual) or residual > TOL:
            raise AssertionError(f'{name}: residual {residual} exceeds {TOL}')
        cases.append({'name': name, 'instances': int(count), 'max_residual': residual,
                      'status': 'PASS', 'detail': detail or ''})

    rz = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
    rx = np.array([[1., 0., 0.], [0., 0., -1.], [0., 1., 0.]])
    s = [np.array(x) for x in [(0., 0., 0.), (.5, 0., 0.), (0., .5, .25)]]
    sbig = [np.array(x) for x in [(0., -.5, 0.), (0., .5, .25), (0., 0., 0.), (.75, .25, 0.), (.5, 0., 0.)]]
    anchors = [np.array([0., 0., 0.]), np.array([.25, -.25, 0.])]
    anchorsbig = [anchors[1], np.array([-.25, .5, 0.]), anchors[0]]
    frames = [np.eye(3), rz]
    framesbig = [rz, rx, np.eye(3)]
    j = np.kron(index_embedding([2, 0], 3), np.kron(embedding(anchors, anchorsbig), np.kron(embedding(s, sbig), I2)))
    d = 2 * len(s) * len(anchors) * len(frames)
    refdim = 3
    rho = density(rng, d * refdim)
    jr = np.kron(j, np.eye(refdim))
    rhobig = jr @ rho @ jr.conj().T
    for reader in ['tent', 'qubit']:
        ks = kraus_catalogue(s, anchors, frames, reader)
        kb = kraus_catalogue(sbig, anchorsbig, framesbig, reader)
        completeness = maxabs(sum(k.conj().T @ k for k in ks) - np.eye(d))
        natural = max(maxabs(b @ j - j @ a) for a, b in zip(ks, kb))
        smallout = instrument_blocks(ks, rho, refdim)
        largeout = instrument_blocks(kb, rhobig, refdim)
        unknown = max(maxabs(lb - jr @ sa @ jr.conj().T) for sa, lb in zip(smallout, largeout))
        marginal = maxabs(sum(partial_second(b, d, refdim) for b in smallout) - partial_second(rho, d, refdim))
        check(f'{reader}_full_instrument_naturality_with_coherent_controls_and_reference',
              max(completeness, natural, unknown, marginal), 2,
              f'Input dimension {d}; larger directory dimension {j.shape[0]}; reference {refdim}; both outcomes and poststates retained.')

    # Transport source, anchor and read-frame together. Output order is changed.
    g = rx @ rz
    b = np.array([.25, -.5, .125])
    gs = [g @ x + b for x in s]
    ga = [g @ a + b for a in anchors]
    gr = [g @ r for r in frames]
    ts, ta, tr = list(reversed(gs)), list(reversed(ga)), list(reversed(gr))
    transport = np.kron(index_embedding([1, 0], 2), np.kron(embedding(ga, ta), np.kron(embedding(gs, ts), I2)))
    for reader in ['tent', 'qubit']:
        kin = kraus_catalogue(s, anchors, frames, reader)
        kout = kraus_catalogue(ts, ta, tr, reader)
        residual = max(maxabs(ko @ transport - transport @ ki) for ki, ko in zip(kin, kout))
        check(f'CO2_{reader}_complete_Kraus_intertwining', residual, 2,
              'Frame R transports to gR and anchor a to ga+b; the probe is not silently replaced.')

    rigid = [lambda x: rz @ x + np.array([.25, 0., 0.]),
             lambda x: rx @ x + np.array([0., -.25, .5])]
    rigidbig = [rigid[1], lambda x: x - np.array([.5, .5, 0.]), rigid[0]]
    scaling = [lambda x: x, lambda x: .5*x]
    scalingbig = [scaling[1], lambda x: .75*x, scaling[0]]
    for name, maps, mapsbig in [('rigid', rigid, rigidbig), ('coordination', scaling, scalingbig)]:
        v, out = controlled_transport(s, maps)
        vb, outbig = controlled_transport(sbig, mapsbig)
        jp = index_embedding([2, 0], 3)
        jin = np.kron(jp, embedding(s, sbig))
        jout = np.kron(jp, embedding(out, outbig))
        iso = maxabs(v.conj().T @ v - np.eye(v.shape[1]))
        natural = maxabs(vb @ jin - jout @ v)
        source = density(rng, v.shape[1]*6)
        vr = np.kron(v, np.eye(6))
        transported = vr @ source @ vr.conj().T
        ref = maxabs(partial_second(transported, v.shape[0], 6) - partial_second(source, v.shape[1], 6))
        check(f'controlled_{name}_isometry_naturality_and_reference', max(iso, natural, ref), 3,
              'Control is retained; the output directory is the union of controlled images, not a falsely bijective common image.')
    for t in [0., .125, .5, .875, 1.]:
        f = lambda x, t=t: (1-t/2)*x
        v, _ = controlled_transport(s, [f])
        check(f'CO3_B_t_{t:g}_isometry', maxabs(v.conj().T @ v - np.eye(len(s))))
    nonzero = [np.array([.2, .1, -.3]), np.array([.9, 0., 0.])]
    if not all(np.linalg.norm(.5*x) < np.linalg.norm(x) for x in nonzero):
        raise AssertionError('coordination budget does not decrease')
    check('CO3_declared_budget_strict_decrease_samples', 0., len(nonzero), 'The budget is Euclidean norm, not energy or elapsed time.')

    # Distinct pairs can have the same midpoint; old records make the map injective.
    mids = union_points([[(x+y)/2 for x in s for y in s]])
    midkeys = keys(mids)
    n = len(s)
    v = np.zeros((n*n*len(mids), n*n), dtype=complex)
    raw_midpoints = []
    for i, x in enumerate(s):
        for k, y in enumerate(s):
            key = keys([(x+y)/2])[0]
            raw_midpoints.append(key)
            v[(i*n+k)*len(mids)+midkeys.index(key), i*n+k] = 1
    if len(set(raw_midpoints)) == len(raw_midpoints):
        raise AssertionError('test failed to include colliding midpoint values')
    source = density(rng, n*n*3)
    vr = np.kron(v, np.eye(3))
    out = vr @ source @ vr.conj().T
    residual = max(maxabs(v.conj().T @ v - np.eye(n*n)),
                   maxabs(partial_second(out, v.shape[0], 3)-partial_second(source, n*n, 3)))
    check('CO4_midpoint_with_both_inputs_retained', residual, n*n,
          f'{n*n} input pairs, {len(mids)} distinct midpoint labels; joint isometry preserves information and reference, but does not preserve the reduced position state after the new record is discarded.')

    cat = np.zeros(n*n, dtype=complex)
    cat[0] = cat[n+1] = 1/np.sqrt(2)
    cat_initial = np.outer(cat, cat.conj())
    cat_output = v @ cat_initial @ v.conj().T
    reduced_old = np.trace(cat_output.reshape(n*n,len(mids),n*n,len(mids)), axis1=1, axis2=3)
    if abs(cat_initial[0,n+1]) < .49:
        raise AssertionError('missing initial position coherence')
    check('CO4_discarding_new_midpoint_dephases_old_position', abs(reduced_old[0,n+1]), 1,
          'The joint isometry retains both original basis labels and quantum information; tracing the midpoint removes this nonzero coherence. No nondisturbing unknown-position read is claimed.')

    # Analytic Lipschitz/tent claims are only calibrated here, not proved by samples.
    lip, distance_witness, anti = 0., 0., 0.
    for _ in range(80):
        x, y = rng.normal(size=3), rng.normal(size=3)
        r = float(np.linalg.norm(x-y))
        lip = max(lip, max(0., norm(effect(x)-effect(y)) - .25*r))
        distance_witness = max(distance_witness, abs(abs(tent(x-x)-tent(y-x))-min(1., r)))
        z = .1*rng.normal(size=3)
        expected = np.linalg.norm(z)/(2*np.sqrt(1+np.dot(z,z)))
        anti = max(anti, abs(norm(effect(z)-effect(-z))-expected))
    check('declared_Q_tent_witness_and_qubit_Lipschitz_bound', max(lip, distance_witness), 80)
    check('CO5a_antipodal_traceless_effect_gap', anti, 80,
          'Gap equals |z|/(2 sqrt(1+|z|^2)); a finite sample does not certify every nonzero endpoint.')
    x = np.array([0., 0., .25])
    ab, ba = rz @ rx @ x, rx @ rz @ x
    signal = abs(tent(ab-ab)-tent(ba-ab))
    if signal <= .3:
        raise AssertionError('noncommuting endpoint witness absent')
    check('CO5b_actual_endpoint_order_witness', abs(signal-np.sqrt(2)/4), 1,
          'Rotations about the same anchor produce a tent-record gap; these are stipulated relation operations, not boosts.')

    # Finite mother-family U calibration using one state and a nonorthogonal ensemble.
    dim, nensemble = 3, 5
    states = []
    for _ in range(nensemble):
        u = rng.normal(size=dim) + 1j*rng.normal(size=dim)
        u /= np.linalg.norm(u)
        states.append(np.outer(u, u.conj()))
    weights = rng.random(nensemble)
    weights /= weights.sum()
    taus = [w*q for w, q in zip(weights, states)]
    target = sum(taus)
    vals, vecs = np.linalg.eigh(target)
    invsqrt = (vecs / np.sqrt(vals)) @ vecs.conj().T
    es = [(invsqrt @ tau @ invsqrt).T for tau in taus]
    psi = np.kron(np.eye(dim), sqrt_psd(target)) @ np.eye(dim).reshape(-1)
    resource = np.outer(psi, psi.conj())
    steering = maxabs(sum(es)-np.eye(dim))
    for e, tau in zip(es, taus):
        k = np.kron(sqrt_psd(e), np.eye(dim))
        steering = max(steering, maxabs(partial_second(k @ resource @ k.conj().T, dim, dim)-tau))
    check('U_one_fixed_resource_all_sampled_ensemble_members', steering, nensemble,
          'Same dimension-3 auxiliary and same purification for this finite ensemble; no unknown-state purifier is claimed.')

    def basis_from_first(vec):
        basis = [vec/np.linalg.norm(vec)]
        for e in np.eye(len(vec), dtype=complex):
            z = e.copy()
            for q in basis:
                z -= q*np.vdot(q,z)
            if np.linalg.norm(z) > 1e-10:
                basis.append(z/np.linalg.norm(z))
            if len(basis) == len(vec):
                break
        return np.column_stack(basis)
    aa = rng.normal(size=4) + 1j*rng.normal(size=4)
    bb = rng.normal(size=4) + 1j*rng.normal(size=4)
    aa, bb = aa/np.linalg.norm(aa), bb/np.linalg.norm(bb)
    ua, ub = basis_from_first(aa), basis_from_first(bb)
    unitary = ub @ ua.conj().T
    evals, evecs = np.linalg.eig(unitary)
    inverse = np.linalg.inv(evecs)
    cres = max(maxabs(unitary.conj().T @ unitary-np.eye(4)), maxabs(unitary @ aa-bb))
    for t in [0., .25, .5, .75, 1.]:
        ut = (evecs*np.exp(1j*t*np.angle(evals))) @ inverse
        cres = max(cres, maxabs(ut.conj().T @ ut-np.eye(4)))
    check('C_finite_pure_state_transport_connected_unitary_path', cres, 5)

    # P: independently mixed Kraus presentation of the same complete channel.
    p = .37
    ks = [np.sqrt(p)*I2] + [np.sqrt((1-p)/3)*z for z in SIGMA]
    fourier = np.exp(2j*np.pi*np.outer(np.arange(4), np.arange(4))/4)/2
    alt = [sum(fourier[i,j]*ks[j] for j in range(4)) for i in range(4)]
    local = 0.
    for i in range(2):
        for j0 in range(2):
            m = np.zeros((2,2), dtype=complex)
            m[i,j0] = 1
            local = max(local, maxabs(sum(k@m@k.conj().T for k in ks)-sum(k@m@k.conj().T for k in alt)))
    qr = density(rng, 6)
    first = sum(np.kron(k,np.eye(3))@qr@np.kron(k.conj().T,np.eye(3)) for k in ks)
    second = sum(np.kron(k,np.eye(3))@qr@np.kron(k.conj().T,np.eye(3)) for k in alt)
    pres = max(local, maxabs(first-second), maxabs(partial_second(first,2,3)-partial_second(qr,2,3)))
    check('P_complete_channel_Kraus_equivalence_and_reference', pres, 5,
          'Local equality checked on a matrix basis; correlated-reference extension and nonsignalling marginal checked.')

    return {
        'schema': 'research_round1090_typed_family_calibration_v1',
        'round': 1090,
        'status': 'PASS',
        'seed': SEED,
        'numpy_version': np.__version__,
        'tolerance': TOL,
        'claim_scope': 'Finite calibration of a declared exact compatible typed-family model; analytic universal statements and actual primitive permissions require the separate proof.',
        'infinite_inductive_family_is_an_F_object': False,
        'fixed_finite_programmer_for_all_continuous_inputs_claimed': False,
        'Lorentz_or_natural_motion_derived': False,
        'empirical_realization_claimed': False,
        'check_groups': len(cases),
        'max_residual': max(c['max_residual'] for c in cases),
        'checks': cases,
    }


def compare(actual, saved, location='root'):
    if isinstance(actual, dict):
        if actual.keys() != saved.keys():
            raise AssertionError(f'{location}: keys differ')
        for key in actual:
            compare(actual[key], saved[key], location+'.'+key)
    elif isinstance(actual, list):
        if len(actual) != len(saved):
            raise AssertionError(f'{location}: lengths differ')
        for i,(a,b) in enumerate(zip(actual,saved)):
            compare(a,b,f'{location}[{i}]')
    elif isinstance(actual, float):
        if not np.isfinite(saved) or abs(actual-saved) > 5e-13:
            raise AssertionError(f'{location}: saved numeric value differs')
    elif actual != saved:
        raise AssertionError(f'{location}: values differ')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true', help='Create results.json once; never overwrite.')
    args = parser.parse_args()
    target = ROOT/'results.json'
    if args.write and target.exists():
        raise SystemExit('Refusing to overwrite existing results.json')
    result = run()
    if args.write:
        with target.open('x', encoding='utf-8', newline='\n') as handle:
            json.dump(result, handle, ensure_ascii=False, indent=2, allow_nan=False)
            handle.write('\n')
        mode = 'created'
    else:
        if not target.exists():
            raise SystemExit('Missing results.json; use --write for first creation only')
        compare(result, json.loads(target.read_text(encoding='utf-8')))
        mode = 'read_only_recompute_match'
    print(json.dumps({'status':'PASS','round':1090,'mode':mode,'check_groups':result['check_groups'],
                      'max_residual':result['max_residual'],'Lorentz_or_natural_motion_derived':False}, ensure_ascii=False))


if __name__ == '__main__':
    main()
