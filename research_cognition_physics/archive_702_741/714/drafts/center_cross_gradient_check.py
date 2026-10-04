"""714 entry supplement: original 589 cross gradients, with no diagonal truncation."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
HERE = Path(__file__).resolve().parent
ARCHIVE = HERE.parent
sys.path.insert(0, str(ARCHIVE))
import joint_full_spatial_metric as metric
import center_reference_entry as entry
TARGET = HERE / 'center_cross_gradient_results.json'


def run():
    rng = np.random.default_rng(7142)
    errors = []
    off_diagonals = []
    for _ in range(24):
        phi = rng.normal(size=(4, 5)) * .3
        links = [entry.old.sample(rng) for _ in range(3)]
        A = rng.normal(size=(3, 3))
        shape = np.eye(3) + A @ A.T
        shape /= np.linalg.det(shape) ** (1/3)
        inverse = np.linalg.inv(shape)
        off_diagonals.append(float(np.linalg.norm(inverse - np.diag(np.diag(inverse)))))
        edge_factors = np.exp(rng.normal(size=3) * .1)
        def average(twist):
            values = []
            for sign in (False, True):
                endpoints = phi[1:].copy()
                if sign:
                    endpoints[0] = entry.scalar_action(entry.CENTER, endpoints[0])
                gs = links.copy()
                if twist:
                    gs[0] = entry.old.product(entry.CENTER, gs[0])
                logs = np.stack([metric.target_log(phi[0], entry.scalar_action(g, p))
                                 for g, p in zip(gs, endpoints)])
                u = edge_factors[:, None] * logs
                gram = u @ metric.original.metric(phi[0]) @ u.T
                values.append(float(np.sum(inverse * gram)))
            return sum(values) / 2
        errors.append(abs(average(False) - average(True)))
    assert max(errors) < 1e-12 and min(off_diagonals) > .1
    return dict(samples=24,full_three_direction_cross_terms_retained=True,
                maximum_endpoint_average_error=max(errors),
                minimum_inverse_shape_off_diagonal_norm=min(off_diagonals),
                dependencies={name: hashlib.sha256((ARCHIVE/name).read_bytes()).hexdigest()
                              for name in ('research_note_589.md', 'joint_full_spatial_metric.py')})


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--write-results', action='store_true')
    result = run()
    if p.parse_args().write_results:
        with TARGET.open('x', encoding='utf8') as f:
            json.dump(result, f, ensure_ascii=False, indent=2)
    else:
        assert json.loads(TARGET.read_text('utf8')) == result
    print(json.dumps(result, ensure_ascii=False, indent=2))
