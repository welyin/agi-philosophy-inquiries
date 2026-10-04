"""771 working entry: density/operator contact algebra, not a loop calculation.

The original 111-field principal block and scalar fiber metric are reused.
Calibration W=0 is only a homogeneous matrix bisolution, NOT a physical state.
The analytic argument retains the inherited physical Hadamard state.
"""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ARCHIVE = HERE.parent
sys.path.insert(0, str(ARCHIVE))
import joint_brst_relative_source as old
import joint_covariant_gauge_complex as cov


def st(a):
    return float(np.real(np.trace(a[:79, :79])-np.trace(a[79:, 79:])))


def scalar_metric_derivative(phi, direction):
    f = 2-phi@phi/6
    df = -phi@direction/3
    return (-df*np.eye(5)/f**2
            +(np.outer(phi, direction)+np.outer(direction, phi))/(6*f**2)
            -df*np.outer(phi, phi)/(3*f**3))


def contact_algebra():
    h, w, k, ks, p, d0, d, _, _, m, _ = old.setup([1.3, .2, -.3, .4])
    x1 = np.linalg.solve(h, np.diag(np.linspace(.03, .17, 63)))
    x0 = np.linalg.solve(-w, np.diag(np.linspace(.021, .043, 16)))
    hp = old.block_matrix([[x1, k@x0, None, None], [x0@ks, None, None, None],
                           [None, None, None, x0], [None, None, x0, None]],
                          [63, 16, 16, 16])
    cp = -hp
    phi = np.array([.13, .64, -.11, .08, .37])
    direction = np.array([.27, -.14, .31, .18, -.09])
    dm = np.zeros_like(m)
    dm[58:63, 58:63] = scalar_metric_derivative(phi, direction)
    theta = np.linalg.solve(m, dm)
    lh = m@d
    rng = np.random.default_rng(771)
    dl = np.zeros_like(m)
    u = np.zeros_like(m)
    for lo, hi in ((0, 79), (79, 111)):
        x = .003*rng.standard_normal((hi-lo, hi-lo))
        dl[lo:hi, lo:hi] = x+x.T
        u[lo:hi, lo:hi] = .03*rng.standard_normal((hi-lo, hi-lo))
    dd = np.linalg.solve(m, dl)-theta@d
    jl = .5*st(cp@np.linalg.solve(m, dl))
    jd = .5*st(cp@dd)
    contact = .5*st(d@hp@theta)
    # A general even change of frame calibrates congruence -> similarity.
    # It is not asserted to be a particular physical background gauge change.
    dmu = -u.conj().T@m-m@u
    dlu = -u.conj().T@lh-lh@u
    thetau = np.linalg.solve(m, dmu)
    ddu = np.linalg.solve(m, dlu)-thetau@d
    jlu = .5*st(cp@np.linalg.solve(m, dlu))
    jdu = .5*st(cp@ddu)
    contactu = .5*st(d@hp@thetau)
    remainder = -.5*st((d@hp-hp@d)@u)
    step = 1e-5
    fd = (cov.scalar_metric(phi+step*direction)
          -cov.scalar_metric(phi-step*direction))/(2*step)
    errors = dict(
        scalar_pairing_derivative=float(np.max(np.abs(fd-dm[58:63, 58:63]))),
        density_to_operator=abs(jd-jl-contact),
        congruence_to_similarity=float(np.max(np.abs(ddu-(u@d-d@u)))),
        gauge_contact_identity=abs(jdu-jlu-contactu),
        commutator_remainder=abs(jdu-remainder),
    )
    assert errors['scalar_pairing_derivative'] < 1e-10
    assert max(v for name, v in errors.items() if name != 'scalar_pairing_derivative') < 1e-12
    assert abs(contact) > 1e-6 and abs(remainder) > 1e-5
    return dict(errors=errors, density_insertion=jl, operator_insertion=jd,
                scalar_pairing_contact=contact, frame_change_contact=contactu,
                residual_commutator=remainder,
                scope='Finite algebra with original principal block and exact scalar pairing derivative. No Hadamard coefficient, anomaly, continuum state or permitted counterterm is computed.')


def scope_calibration():
    # Moretti's c_D and k_D share all factors except (D+2)/2^(D-1)
    # versus D/2^D; the ratio fixes the scalar EOM contact only.
    ratios = []
    for dim in (2, 4, 6, 8):
        cd = Fraction(dim+2, 2**(dim-1))
        kd = Fraction(dim, 2**dim)
        ratio = kd/cd
        assert ratio == Fraction(dim, 2*(dim+2))
        ratios.append(dict(dimension=dim, scalar_eta=str(ratio)))
    # Only operator differential orders, not background-jet orders or
    # derivatives of the singular kernels separately at the diagonal.
    remainder_orders = dict(r1=2, r0=2, q=1, Kq=2)
    required_orders = {name: order+1 for name, order in remainder_orders.items()}
    assert max(required_orders.values()) == 3
    return dict(scalar_literature_ratio=ratios,
                upper_bounds_after_one_gauge_derivative=required_orders,
                not_background_jet_order_bound=True,
                scope='Exact rational and differential-order bookkeeping only. Neither a matrix-field generalization of the scalar theorem nor the actual local Ward repair.')


def run():
    deps = ('research_note_734.md', 'research_note_735.md', 'research_note_765.md',
            'research_note_770.md', 'joint_brst_relative_source.py',
            'joint_covariant_gauge_complex.py')
    return dict(working_round=771, latest_completed_round=770,
                new_scientific_tests=0, entry_calibrations=2,
                contact_algebra=contact_algebra(), scope_calibration=scope_calibration(),
                dependency_hashes={n: hashlib.sha256((ARCHIVE/n).read_bytes()).hexdigest() for n in deps},
                actual_local_Ward_repair_proven=False,
                original_physical_state_not_replaced=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    result = run()
    target = HERE/'local_contact_entry_results.json'
    if args.write_results:
        with target.open('x', encoding='utf8') as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
    else:
        assert result == json.loads(target.read_text('utf8'))
    print(json.dumps({k: result[k] for k in ('working_round', 'latest_completed_round', 'new_scientific_tests', 'entry_calibrations')}))
