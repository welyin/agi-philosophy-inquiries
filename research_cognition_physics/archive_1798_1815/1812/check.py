"""Weak-CP budget for two split five-pairs with a common real adjoint."""
import argparse
import hashlib
import json
import math
import runpy
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
IN = json.loads((HERE / "inputs.json").read_text(encoding="utf-8"))


def load_record(field):
    return json.loads((HERE / IN[field]).read_text(encoding="utf-8"))


def pair_sum(edges):
    a, b, c = edges
    return a*b + b*c + c*a


def endpoint(data, band, gains, raw):
    y = np.array(data["high_yukawas"])
    *sines, phase = data["high_ckm_sines_delta"]
    mixing = runpy.run_path(str(HERE / IN["ckm_parameterization"]))["mixing"]
    v = mixing(np.arcsin(sines), phase)
    down = (y[3:6]/y[5])**2
    n = down[:2].sum() + np.sum((y[6:]/y[5])**2)
    h = (v*down) @ v.conj().T
    z = h[[0,1,2], [1,2,0]]
    gap = (1-down[0])*(1-down[1])*(down[1]-down[0])
    edge_budget = pair_sum(abs(z))
    j = abs((v[0,0]*v[1,1]*v[0,1].conj()*v[1,0].conj()).imag)
    np.testing.assert_allclose(abs(np.prod(z).imag), j*gap, rtol=1e-10, atol=0)
    need = j*gap/(n*edge_budget)
    d_pole,t_pole = data["D_T_poles_GeV"]
    gauge = gains["gauge_only_MS_mass_gains_from_each_pole_to_matching"]
    d_q, t_q = d_pole*gauge["D"], t_pole*gauge["T"]
    r_pole_max = math.exp(-band["original_window"]["B_log_min_T_over_D"][0])
    weak_j_cap = j * (d_q/t_q) / need
    return dict(physical_J_at_Q=j, down_mass_squared_over_b_squared=down.tolist(),
                trace_YdaggerY_over_b_squared=n, edge_moduli_over_b_squared=abs(z).tolist(),
                gap_product=gap, cyclic_modulus_pair_sum=edge_budget,
                rQ_necessary_lower=need, r_pole_whole_band_max=r_pole_max,
                necessary_rQ_to_rpole_gain_whole_band_fixed_flavor=need/r_pole_max,
                central_example_pole_ratio=d_pole/t_pole, gauge_only_central_Q_ratio=d_q/t_q,
                gauge_only_ratio_gain=gauge["D"]/gauge["T"],
                gauge_only_Q_ImA_necessary_GeV=max(0, t_q*need-d_q)/5,
                unit_mass_gain_ImA_necessary_GeV=max(0, t_pole*need-d_pole)/5,
                gauge_only_J_upper_at_measured_edge_moduli_and_masses=weak_j_cap,
                maximum_J_high_over_low_if_raw_B_minimum={k:weak_j_cap/raw[k]["KK_quartet_imaginary_abs_lower"] for k in ("1","2")})


def broad_box():
    b = IN["conditional_high_flavor_box"]
    a0,a1=b["s12"]
    b0,b1=b["s13"]
    c0,c1=b["s23"]
    y0,y1=b["ys_over_yb"]
    d=b["yd_over_ys_max"]**2
    j_min=a0*b0*c0*math.sqrt(1-a1*a1)*math.sqrt(1-c1*c1)*(1-b1*b1)*b["abs_sin_delta_min"]
    gap_min=(1-y1*y1)*(1-d*y1*y1)*(1-d)*y0*y0
    # Off diagonals may subtract y_d² I exactly, by CKM unitarity.
    edges=[b1*c1+y1*y1*a1,c1+y1*y1,b1+y1*y1*a1]
    return dict(J_lower=j_min, gap_lower=gap_min, edge_upper=edges,
                rQ_necessary_lower=j_min*gap_min/(b["trace_over_yb2_max"]*pair_sum(edges)))


def nonzero_cp_example(r):
    # Each 2-plane supplies one orthogonal T,D pair with identical imaginary row.
    t,d = np.zeros((2,5),complex),np.zeros((2,5),complex)
    for row in range(2):
        t[row,2*row:2*row+2]=[1,1j*r]
        d[row,2*row:2*row+2]=[-r*r,1j*r]
    def projector(m):
        return m.conj().T @ np.linalg.solve(m @ m.conj().T, m)
    pt,pd=projector(t),projector(d)
    y0=np.array([[1,.4,.2,1.2,0],[.2,.7,1.1,-.4,0],[0,0,0,0,3.]])
    # A real rotation relative to the diagonal nondegenerate up spectrum.
    q,_=np.linalg.qr(np.array([[1.,2,3],[2,1,4],[3,1,2]]))
    y=q@ y0
    hd=y@(np.eye(5)-pt)@y.T
    he=y@(np.eye(5)-pd)@y.T
    m2,v=np.linalg.eigh(hd)
    j=abs((v[0,0]*v[1,1]*v[0,1].conj()*v[1,0].conj()).imag)
    np.testing.assert_allclose(t@d.conj().T,0,atol=1e-15)
    np.testing.assert_allclose((t-d).imag,0,atol=1e-15)
    np.testing.assert_allclose(np.trace(hd+he).real-m2[-1],np.linalg.norm(y)**2)
    np.testing.assert_allclose(np.linalg.norm(pt.imag,2),r/(1+r*r))
    assert j > 1e-8
    return dict(rQ=r, projector_imag_norm=np.linalg.norm(pt.imag,2),
                weak_J=j, shared_bottom_squared=m2[-1])


def compute():
    fits,band,gains=map(load_record,("matching","gauge_band","central_mass_transport"))
    raw=load_record("direct_B")["bound_tables"]
    result={k:endpoint(fits[k],band[k],gains[k],raw) for k in ("common_MZ","independent_lattice")}
    box=broad_box()
    box["rQ_to_rpole_gain_necessary_over_gauge_union"]=box["rQ_necessary_lower"]/max(p["r_pole_whole_band_max"] for p in result.values())
    assert box["rQ_to_rpole_gain_necessary_over_gauge_union"]>1
    for name in result:
        assert result[name]["necessary_rQ_to_rpole_gain_whole_band_fixed_flavor"]>160
    result.update(conditional_broad_flavor_box=box,finite_soft_CP_examples=[nonzero_cp_example(r) for r in (.01,.001)])
    return result


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    checking=p.parse_args().check
    text=json.dumps(compute(),indent=2,ensure_ascii=False,allow_nan=False)+"\n"
    file=HERE/"results.json"
    if checking:
        assert json.loads(file.read_text(encoding="utf-8"))==json.loads(text)
    else:
        file.write_text(text,encoding="utf-8")
    print("1812",hashlib.sha256(text.encode()).hexdigest()[:12])
