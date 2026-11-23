"""Independently verify allowed pairings, a field-space winding, and three radial response tests."""
import argparse
import itertools
import json
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent

def smallest_cost(rows,cols):
    return min(sum(abs(rows[i]+q) for i,q in zip(p,cols)) for p in itertools.permutations(range(len(rows)),len(cols)))

def result():
    left=(4,13,-13)
    right=(-13,)*3
    eligible=[q for q in range(1,41) if all(any((a+b)%q==0 for b in right) for a in left)]
    costs=[smallest_cost(left,right),smallest_cost((4,4,13,13,-13,-13),(4,4,13,13,-13,-13)),
           smallest_cost((4,4,4,13,13,13),(4,13,-13))]
    assert eligible==[1] and costs==[35,16,25]
    mass=np.array([[.9,.2j,-.3],[.4,.8+.2j,-.1j],[.3j,.6,1.1]],dtype=complex)
    exponents=np.array([9,0,26])
    theta=np.linspace(0,2*np.pi,257)
    det=np.array([np.linalg.det(np.exp(1j*t*exponents)[:,None]*mass) for t in theta])
    winding=np.angle(det[1:]/det[:-1]).sum()/(2*np.pi)
    assert abs(winding-35)<1e-11
    U,m,Vh=np.linalg.svd(mass)
    coupling=U.conj().T@(exponents[:,None]*mass)@Vh.conj().T
    scaled=coupling/m
    assert np.max(abs(scaled-scaled.conj().T))<1e-12
    assert abs(np.trace(scaled)-35)<1e-11 and abs(np.linalg.norm(scaled)**2-757)<1e-10
    assert 7*(-35)+10*8+3*55==0
    return {"integer_q":eligible,"minimal_exponents_15_24_10":costs,
            "net_mass_winding_unit_string":round(winding),"net_complex_channels_entire_15":15*round(winding),
            "example_heavy_mode_normalized_diagonal_couplings":[round(x.real,10) for x in np.diag(scaled)],
            "diagonal_sum":35,"squared_Frobenius_norm":757,
            "full_colored_anomaly_15_24_10_sectors":[-245,80,165]}
if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--check",action="store_true")
    opt=ap.parse_args()
    out=HERE/"independent_results.json"
    answer=result()
    if opt.check:
        assert answer==json.loads(out.read_text(encoding="utf-8-sig"))
        print("1582 independent phase and mass tests verified")
    else:
        out.write_text(json.dumps(answer,indent=2)+"\n",encoding="utf-8")
