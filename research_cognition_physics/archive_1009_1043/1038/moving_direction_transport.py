"""1038: canonical-spin comparison across adopted Lorentz frames.

Finite checks calibrate the analytic statements; no detector or full QFT claim.
Default execution is read-only. --write exclusively creates the first result.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
OUT = HERE / 'moving_direction_transport_results.json'
HISTORY = [
    'archive_1009_/1009/input_dependency_ledger_v0_1.md',
    'archive_956_989/957/drafts/unified_operation_hypotheses_v0_2.md',
    'archive_1009_/research_note_1035.md',
    'archive_370_428/research_note_372.md',
    'archive_370_428/research_note_382.md',
    'archive_370_428/research_note_383.md',
    'archive_370_428/research_note_384.md',
    'archive_370_428/research_note_386.md',
    'archive_370_428/research_note_425.md',
    'archive_467_530/research_note_522.md',
    'archive_467_530/research_note_523.md',
    'archive_702_741/research_note_720.md',
    'archive_923_934/research_note_923.md',
    'archive_342_369/research_note_369.md',
    'archive_956_989/research_note_963.md',
]
I = np.eye(2, dtype=complex)
SIGMA = np.array([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]],
                  [[1, 0], [0, -1]]], dtype=complex)


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def dot_sigma(v):
    return np.einsum('...i,ijk->...jk', v, SIGMA)


def opnorm(a):
    return float(np.linalg.norm(a, 2))


def energy(p, mass=1.):
    return np.sqrt(mass * mass + np.sum(np.asarray(p)**2, axis=-1))


def standard_boost(p, mass=1.):
    p = np.asarray(p, float)
    e = energy(p, mass)
    return ((e + mass) * I + dot_sigma(p)) / np.sqrt(2 * mass * (e + mass))


def sl_boost(n, rapidity):
    n = np.asarray(n, float)
    assert np.isclose(np.linalg.norm(n), 1.)
    return np.cosh(rapidity / 2) * I + np.sinh(rapidity / 2) * dot_sigma(n)


def transform_momentum(a, p, mass=1.):
    x = a @ (energy(p, mass) * I + dot_sigma(p)) @ a.conj().T
    return np.array([np.trace(x @ s).real / 2 for s in SIGMA])


def wigner_matrix(a, p, mass=1.):
    pp = transform_momentum(a, p, mass)
    return np.linalg.solve(standard_boost(pp, mass), a @ standard_boost(p, mass))


def wigner_closed(p, n, rapidity, mass=1.):
    p, n = np.asarray(p, float), np.asarray(n, float)
    v = p / (energy(p, mass)[..., None] + mass)
    t = np.tanh(rapidity / 2)
    aa = 1 + t * np.einsum('...i,i->...', v, n)
    bb = t * np.cross(n, v)
    norm = np.sqrt(aa * aa + np.sum(bb * bb, axis=-1))
    return (aa[..., None, None] * I + 1j * dot_sigma(bb)) / norm[..., None, None]


def representation_checks():
    rng = np.random.default_rng(1038)
    err = dict(closed_vs_induced=0., unitarity=0., determinant=0.,
               mass_shell=0., composition=0., rotation_independent=0.,
               infinitesimal=0., velocity_recovery=0., simultaneous_dictionary=0.)
    for _ in range(36):
        p = rng.normal(size=3) * 1.7
        n = rng.normal(size=3); n /= np.linalg.norm(n)
        r = rng.uniform(-1.3, 1.3)
        a = sl_boost(n, r)
        w = wigner_matrix(a, p)
        err['closed_vs_induced'] = max(err['closed_vs_induced'], opnorm(w - wigner_closed(p, n, r)))
        err['unitarity'] = max(err['unitarity'], opnorm(w.conj().T @ w - I))
        err['determinant'] = max(err['determinant'], float(abs(np.linalg.det(w) - 1)))
        pp = transform_momentum(a, p)
        ep = np.trace(a @ (energy(p) * I + dot_sigma(p)) @ a.conj().T).real / 2
        err['mass_shell'] = max(err['mass_shell'], float(abs(ep * ep - pp @ pp - 1)))
        b = sl_boost(np.array([0., 1., 0.]), -.43)
        err['composition'] = max(err['composition'], opnorm(wigner_matrix(b @ a, p) - wigner_matrix(b, pp) @ w))
        rotation = np.cos(.31) * I - 1j * np.sin(.31) * dot_sigma(n)
        err['rotation_independent'] = max(err['rotation_independent'], opnorm(wigner_matrix(rotation, p) - rotation))
        v = p / (energy(p) + 1)
        h = 1e-5
        for axis in range(3):
            direction = np.eye(3)[axis]
            wp = wigner_closed(p, direction, h)
            wm = wigner_closed(p, direction, -h)
            deriv = (wp.conj().T @ SIGMA[axis] @ wp - wm.conj().T @ SIGMA[axis] @ wm) / (2*h)
            target = dot_sigma(v) - v[axis] * SIGMA[axis]
            err['infinitesimal'] = max(err['infinitesimal'], opnorm(deriv - target))
            for other in range(3):
                if other != axis:
                    recovered = (SIGMA[other] @ deriv + deriv @ SIGMA[other]) / 2
                    err['velocity_recovery'] = max(err['velocity_recovery'], opnorm(recovered - v[other]*I))
        ket = rng.normal(size=2) + 1j*rng.normal(size=2); ket /= np.linalg.norm(ket)
        rho = np.outer(ket, ket.conj())
        effect = (I + dot_sigma(n)) / 2
        changed = w @ effect @ w.conj().T
        err['simultaneous_dictionary'] = max(err['simultaneous_dictionary'], float(abs(np.trace(rho@effect) - np.trace(w@rho@w.conj().T@changed))))
    assert max(v for k, v in err.items() if k not in ('infinitesimal', 'velocity_recovery')) < 2e-12
    assert err['infinitesimal'] < 1e-9 and err['velocity_recovery'] < 1e-9
    return dict(cases=36, maximum_residuals=err)


def packet_quadrature(radial, polar, azimuth, width=.01):
    """Integrate a C-infinity compact momentum bump in invariant d^3p/(2E)."""
    xr, wr = np.polynomial.legendre.leggauss(radial)
    z, wz = np.polynomial.legendre.leggauss(polar)
    rr = (xr + 1) / 2
    phi = 2*np.pi*np.arange(azimuth)/azimuth
    r, zz, ph = np.meshgrid(rr, z, phi, indexing='ij')
    sin = np.sqrt(1-zz*zz)
    offsets = width*np.stack([r*sin*np.cos(ph), r*sin*np.sin(ph), r*zz], axis=-1).reshape(-1,3)
    weights = (wr[:,None,None]/2 * wz[None,:,None] * np.ones((1,1,azimuth)) * 2*np.pi/azimuth)
    weights = (weights * r*r * np.exp(-2/(1-r*r))).reshape(-1)
    plus = offsets + np.array([0.,0.,4/3])
    minus = -plus
    weights /= 2*energy(plus)
    weights /= weights.sum()
    rapidity = math.log(3)
    wp, wm = [wigner_closed(pp, np.array([1.,0.,0.]), rapidity) for pp in (plus, minus)]
    pz = [(I+SIGMA[2])/2, (I-SIGMA[2])/2]
    def mean_state(w, rho):
        return np.einsum('n,nij,jk,nlk->il', weights, w, rho, w.conj())
    out_a = (mean_state(wp,pz[0])+mean_state(wm,pz[1]))/2
    out_b = (mean_state(wp,pz[1])+mean_state(wm,pz[0]))/2
    px=(I+SIGMA[0])/2
    pa,pb=[float(np.trace(o@px).real) for o in (out_a,out_b)]
    return dict(points_per_packet=len(weights), probability_a=pa, probability_b=pb,
                probability_gap=pa-pb, min_output_eigenvalue=float(min(np.linalg.eigvalsh(out_a).min(),np.linalg.eigvalsh(out_b).min())),
                max_energy=float(energy(plus).max()), trace_error=float(max(abs(np.trace(out_a)-1),abs(np.trace(out_b)-1))))


def exact_certificates():
    gap=F(8,17); width=F(1,100); lipschitz=F(1,2)
    lower=gap-2*lipschitz*width
    assert lower==F(783,1700)
    return dict(central_probability_a=str(F(25,34)), central_probability_b=str(F(9,34)),
                central_gap=str(gap), packet_radius=str(width), mass='1',
                rapidity='log(3)', velocity_coordinate_at_centers=['0','0','+/-1/2'],
                wigner_lipschitz=str(lipschitz), packet_gap_lower=str(lower),
                any_marginals_only_prediction_error_lower=str(lower/2),
                all_momenta_bounded_by=str(F(4,3)+width),
                total_prediction_errors_example=['1/100','1/100'],
                surviving_gap_if_those_errors_certified=str(lower-F(1,50)),
                physical_detector_errors_certified=False)


def finite_algebra_check():
    # Exact rational velocities, with distinct x values, generate four selectors.
    vals=[(F(0),F(0),F(0)),(F(1,10),F(1,5),F(0)),
          (F(-1,5),F(0),F(3,10)),(F(3,10),F(-1,10),F(1,5))]
    x=[v[0] for v in vals]
    selectors=[]
    for k in range(4):
        row=[]
        for xi in x:
            product=F(1)
            for j,xj in enumerate(x):
                if j!=k: product *= (xi-xj)/(x[k]-xj)
            row.append(product)
        selectors.append(row)
    assert selectors==[[F(int(k==j)) for j in range(4)] for k in range(4)]
    p=[]
    for v in vals:
        sq=sum(a*a for a in v)
        pp=[2*a/(1-sq) for a in v]
        ee=(1+sq)/(1-sq)
        assert ee*ee-sum(a*a for a in pp)==1
        p.append([str(a) for a in pp])
    basis=[]
    for k in range(4):
        for s in [I,*SIGMA]:
            a=np.zeros((8,8),complex);a[2*k:2*k+2,2*k:2*k+2]=s
            basis.append(a.reshape(-1))
    rank=int(np.linalg.matrix_rank(np.array(basis)))
    assert rank==16
    return dict(blocks=4, exact_selector_identity=True, algebra_dimension=rank,
                rational_mass_shell_momenta=p,
                infinite_closure_proved_by='analytic derivative and spectral-calculus proof, not this rank test',
                algebra_closure_is_an_extra_task_contract=True,
                single_shot_linear_span_classified=False)


def reference_and_binning_check():
    rng=np.random.default_rng(9038)
    centers=np.array([[0.,0.,4/3],[0.,0.,-4/3]])
    mom=[]; bin_index=[]
    for k,c in enumerate(centers):
        for _ in range(5):
            r=rng.normal(size=3);r*=.009*rng.random()/np.linalg.norm(r)
            mom.append(c+r);bin_index.append(k)
    w=wigner_closed(np.array(mom), np.array([1.,0.,0.]), math.log(3))
    wc=wigner_closed(centers, np.array([1.,0.,0.]), math.log(3))
    # Same unknown-reference input, including momentum coherences; no postselection.
    z=rng.normal(size=(60,4))+1j*rng.normal(size=(60,4))
    rho=z@z.conj().T;rho/=np.trace(rho)
    rr=rho.reshape(10,6,10,6)
    exact=np.zeros((6,6),complex);approx=np.zeros((6,6),complex)
    blocks=[np.zeros((6,6),complex) for _ in centers]
    err=0.
    for k in range(10):
        u=np.kron(w[k],np.eye(3));v=np.kron(wc[bin_index[k]],np.eye(3))
        b=rr[k,:,k,:]
        exact+=u@b@u.conj().T
        blocks[bin_index[k]]+=b
        err=max(err,opnorm(w[k]-wc[bin_index[k]]))
    for k,b in enumerate(blocks):
        v=np.kron(wc[k],np.eye(3));approx+=v@b@v.conj().T
    dist=float(np.abs(np.linalg.eigvalsh(exact-approx)).sum()/2)
    radius=float(max(np.linalg.norm(p-centers[k]) for p,k in zip(mom,bin_index)))
    analytic=radius/2
    assert dist<=err+1e-13 and err<=analytic+1e-13
    assert abs(np.trace(exact)-1)<1e-13 and min(np.linalg.eigvalsh(b).min() for b in blocks)>-1e-13
    return dict(momentum_bins=2, quadrature_labels=10, reference_dimension=3,
                input_dimension=60, trace_distance=dist, largest_unitary_deviation=err,
                max_bin_radius=radius, analytic_reference_uniform_bound=analytic,
                tail_weight=0., postselection=False, finite_block_dictionary_positive=True)


def run():
    coarse=packet_quadrature(18,8,16)
    fine=packet_quadrature(26,12,24)
    exact=exact_certificates()
    lower=float(F(exact['packet_gap_lower']))
    assert fine['probability_gap']>lower
    return dict(round=1038, date='2026-10-08', new_calibration_groups=1,
                new_adopted_cognitive_axioms=0, cumulative_count_assigned_by_root=True,
                code_sha256=sha(Path(__file__)),
                historical_source_sha256={p:sha(BASE/p) for p in HISTORY},
                representation=representation_checks(), exact_certificates=exact,
                normal_packet={'coarse':coarse,'fine':fine,
                               'quadrature_difference':abs(coarse['probability_gap']-fine['probability_gap']),
                               'quadrature_is_not_rigorous_error_bound':True,
                               'same_complete_momentum_marginal':True,'same_spin_marginal':'I/2',
                               'finite_energy_and_all_momentum_moments':True},
                algebra_closure=finite_algebra_check(), finite_reference_transport=reference_and_binning_check(),
                scope={'massive_free_one_particle':True,'canonical_spin_adopted':True,
                       'lorentz_and_dimension_generated':False,'actual_boost_device_certified':False,
                       'actual_stern_gerlach_instrument_certified':False,'full_parent_quantum_process_certified':False,
                       'continuum_dephasing_on_original_hilbert_space_claimed':False,
                       'single_shot_statistics_imply_algebra_closure':False},
                goal_completed=False)


def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys(),(a.keys(),b.keys())
        for k in a: compare(a[k],b[k])
    elif isinstance(a,(list,tuple)):
        assert len(a)==len(b)
        for x,y in zip(a,b): compare(x,y)
    elif isinstance(a,float):
        assert math.isclose(a,b,rel_tol=1e-9,abs_tol=2e-12),(a,b)
    else: assert a==b,(a,b)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true')
    args=parser.parse_args();result=run()
    if args.write:
        with OUT.open('x',encoding='utf8') as stream:
            json.dump(result,stream,ensure_ascii=False,indent=2,allow_nan=False);stream.write('\n')
    else: compare(result,json.loads(OUT.read_text('utf8')))
    print(json.dumps(dict(round=1038,packet_gap=result['normal_packet']['fine']['probability_gap'],
                          rigorous_gap=result['exact_certificates']['packet_gap_lower'],
                          max_representation_residual=max(result['representation']['maximum_residuals'].values()),
                          finite_reference_error=result['finite_reference_transport']['trace_distance']),ensure_ascii=False))
