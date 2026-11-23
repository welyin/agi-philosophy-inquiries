"""Integer lattice and conditional three-lepton spectral overlap. No source amplitudes fitted."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

BASE=Path(__file__).resolve().parent


def read(name):
    return json.loads((BASE/name).read_text(encoding="utf-8"))


def det(a,b):
    return a[0]*b[1]-a[1]*b[0]


def order(rows, anomaly=(3,3)):
    area=math.gcd(*(abs(det(a,b)) for a,b in itertools.combinations(rows,2)))
    mixed=math.gcd(*(abs(det(a,anomaly)) for a in rows))
    if area:
        return area//math.gcd(area,mixed)
    if mixed or not any(any(r) for r in rows):
        return None
    r=next(r for r in rows if any(r))
    i=0 if r[0] else 1
    stride=math.gcd(*(abs(a[i]) for a in rows))
    return stride//math.gcd(stride,anomaly[i])


def propagation(p,offset,over):
    # q phases for L/E in km/GeV, mass² in eV². Equal zero-mode flavour weights.
    phases=2*1.2669328*over*offset
    w=np.sqrt(p)*np.exp(1j*np.array([.23,-.7,.4]))
    matrix=np.column_stack([w*np.exp(-1j*k*phases) for k in (0,1,2)])
    gram=matrix.conj().T@matrix
    exact=64*np.prod(p)*math.prod(math.sin((a-b)/2)**2 for a,b in itertools.combinations(phases,2))
    assert np.isclose(abs(np.linalg.det(matrix))**2,exact,rtol=1e-11)
    assert np.isclose(np.linalg.det(gram).real,exact,atol=1e-15)
    assert np.linalg.matrix_rank(matrix)==3
    assert np.linalg.matrix_rank(np.repeat(matrix[:,1:2],3,axis=1))==1
    pwr=np.column_stack([w*offset**k for k in (0,1,2)])
    vand=math.prod(b-a for a,b in itertools.combinations(offset,2))*np.prod(w)
    assert np.isclose(np.linalg.det(pwr),vand,rtol=1e-10,atol=1e-20)
    return float(exact)


def audit():
    source=read("../1485/inputs.json")
    s=source["solar_sin2"][0]
    amp=read("../1476/inputs.json")["angles"]["reactor_amplitude"][0]
    r=(1-math.sqrt(1-amp))/2
    p=np.array([(1-r)*(1-s),(1-r)*s,r])
    spl=read("../1446/inputs.json")["juno_207_days"]
    dm=spl["solar_21"][0]
    values={}
    for name,sign,key in [("NO",1,"normal_atmospheric_31"),("IO",-1,"inverted_abs_31")]:
        spectrum=np.array([0.,dm,sign*spl[key][0]])
        values[name]=propagation(p,spectrum,500.)
    cases={"QQQL": [(1,1)],"L2":[(0,2)],"B2":[(2,0)],
           "B2_and_L2":[(2,0),(0,2)],"B4_and_L2":[(4,0),(0,2)],"B1_L3":[(1,3)]}
    harmonic={name:order(rows) for name,rows in cases.items()}
    assert list(harmonic.values())==[1,None,None,2,4,None]
    return {"integer_first_harmonic":harmonic,"oscillation_electron_probabilities":p.tolist(),
            "conditional_flavor_antisymmetric_projection_abs_det_squared_L_over_E_500":values,
            "for_normalized_wedge_projector_divide_each_by":6,
            "all_three_same_filter_rank":1,"minimum_total_polynomial_power":3}


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    m = p.add_mutually_exclusive_group()
    m.add_argument("--write-new", action="store_true")
    m.add_argument("--check", action="store_true")
    args = p.parse_args()
    path = BASE/"independent_results.json"
    out = json.loads(json.dumps(audit(), default=float))
    if args.write_new:
        with path.open("x", encoding="utf-8") as f:
            f.write(json.dumps(out, indent=2)+"\n")
    elif args.check:
        assert out == json.loads(path.read_text(encoding="utf-8"))
        print("Independent tests passed.")
    else:
        print(json.dumps(out, indent=2))
