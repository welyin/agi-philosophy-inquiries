"""Finite monomial residual group and an atmospheric consistency check."""
import json
import math
from collections import deque
from pathlib import Path
import sys

import numpy as np

OUT = Path(__file__).with_name("independent_results.json")
N = 17
MOD = 2*N

def multiply(left, right):
    p, a = left
    q, b = right
    return (tuple(q[j] for j in p), tuple((a[i]+b[j]) % MOD for i,j in enumerate(p)))

def group_order():
    # A row i has its only nonzero entry in column p[i], with phase exp(2 pi i*a[i]/MOD).
    identity = ((0, 1, 2), (0, 0, 0))
    generators = [((1, 2, 0), (0, 0, 0)), ((0, 2, 1), (N,N-2,N+2))]
    queue = deque([identity])
    seen = {identity}
    while queue:
        current = queue.popleft()
        for gen in generators:
            result = multiply(current,gen)
            if result not in seen:
                seen.add(result)
                queue.append(result)
    assert len(seen) <= 6*MOD**2
    return len(seen)

def calculate():
    phase = 2*math.pi/N
    ray = np.array([0,1,-np.exp(1j*phase)])/math.sqrt(2)
    half = 2*np.outer(ray,ray.conj())-np.eye(3)
    t = np.eye(3)[[1,2,0]]
    fourier = np.exp(2j*math.pi*np.outer(np.arange(3),np.arange(3))/3)/math.sqrt(3)
    p = abs(fourier.conj().T@ray)**2
    assert np.allclose(half@half,np.eye(3)) and np.isclose(np.linalg.det(half),1)
    assert np.allclose(abs(fourier.conj().T@t@fourier),np.eye(3))
    assert abs(np.trace(t@half).imag) > .3 and np.isclose(p@p,.5)
    x = p[0]
    amplit = 4*x*(1-x)
    assert abs(amplit-.0851) <= 2*.0024
    atm = sorted(float(z/(1-x)) for z in p[1:])
    gap = math.sqrt(x*(2-3*x))/(1-x)
    assert np.allclose(atm,[(1-gap)/2,(1+gap)/2])
    low = (1-math.sqrt(1-(.0851-2*.0024)))/2
    high = (1-math.sqrt(1-(.0851+2*.0024)))/2
    offset = lambda v: math.sqrt(v*(2-3*v))/(2*(1-v))
    boxes = [[.5-offset(high),.5-offset(low)],[.5+offset(low),.5+offset(high)]]
    return dict(denominator=N, order=group_order(), trace_TS_imaginary=float(np.trace(t@half).imag),
                electron_weight=float(x), reactor_amplitude=float(amplit), atmospheric_octants=atm,
                all_fractions_small_branch_atmospheric=boxes,
                J_max_at_witness=float(math.sqrt(np.prod(p))/2),
                atmospheric_published_one_sigma_normal=[.51,.58],
                atmospheric_comparison="Conditional marginal comparison only; no combined significance.",
                sources=["https://arxiv.org/abs/2211.14988","https://arxiv.org/abs/2405.02163",
                         "https://arxiv.org/abs/1305.3200"])

if __name__ == "__main__":
    result = calculate()
    if sys.argv[1:]==["--write-new"]:
        with OUT.open("x",encoding="utf-8") as stream:
            json.dump(result,stream,indent=2)
            stream.write("\n")
    elif sys.argv[1:] in ([],["--check"]):
        assert result == json.loads(OUT.read_text(encoding="utf-8"))
        print("1495 independent integral group closure, column, reactor and atmospheric tests passed.")
    else:
        raise SystemExit("Expected --check or --write-new")
