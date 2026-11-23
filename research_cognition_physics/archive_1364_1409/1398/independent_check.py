"""Phase-independent two-sided light-current witness."""
import argparse
import json
from math import sqrt
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE/name).read_text(encoding="utf-8"))


def calculate():
    weak = read("../1381/inputs.json")
    flavor = read("../1386/inputs.json")["yukawa"]
    results = {}
    for n in (0, 1, 2):
        def edge(pair, sign):
            return pair[0]+sign*n*sum(pair[1:])
        bmax = edge(weak["K_product"], 1)/edge(weak["K_form"], -1)
        matrix = np.array([[edge(weak["ud_beta"], -1), -bmax],
                           [-edge(weak["charm_d"], 1), edge(weak["charm_s"], -1)]])
        singular = float(np.linalg.svd(matrix, compute_uv=False)[-1])
        assert np.min(np.linalg.inv(matrix)) >= 0
        assert np.isclose(singular, 1/np.linalg.norm(np.linalg.inv(matrix), 2))
        ratios = [(edge(flavor[sector][1], 1)/edge(flavor[sector][2], -1))**2
                  for sector in ("up","down")]
        f = sum(ratios)
        results[str(n)] = {"comparison_matrix": matrix.tolist(), "min_light_sigma": singular,
                           "max_visible_squared_on_some_input": f,
                           "max_inert_sigma": sqrt(f),
                           "lower_extra_norm_even_if_coherent": (singular-sqrt(f))/sqrt(1-f),
                           "min_apparent_gain_for_inert": singular/sqrt(f)}
    return results


def verify(a, b):
    if isinstance(a, dict):
        assert a.keys() == b.keys()
        for k in a:
            verify(a[k], b[k])
    else:
        assert np.allclose(a, b, atol=2e-12, rtol=2e-10)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--check", action="store_true")
    result = calculate()
    if p.parse_args().check:
        verify(result, read("independent_results.json"))
        print("1398 two-by-two worst-phase singular value and leakage verified")
    else:
        print(json.dumps(result, indent=2, allow_nan=False))
