"""Finite dihedral rays and the measured electron row (all masses/fields are synthetic in the witness)."""
import argparse
import cmath
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE / name).read_text(encoding="utf-8-sig"))


def row_bounds(solar, amplitude, widths):
    v, err = solar[:2]
    A, err_A = amplitude[:2]
    corners = []
    for a in (A - widths * err_A, A + widths * err_A):
        q = (1 - math.sqrt(1 - a)) / 2
        for t in (v - widths * err, v + widths * err):
            corners.append([(1 - t) * (1 - q), t * (1 - q), q])
    return [(min(col), max(col)) for col in zip(*corners)]


def angle(p):
    return math.degrees(math.acos(math.sqrt(max(0., min(1., p)))))


def certificate(bounds, grid, atm):
    cols = []
    for lo, hi in bounds:
        gaps = [max(angle(hi) - angle(p), angle(p) - angle(lo), 0) for p in grid]
        j = min(range(len(grid)), key=gaps.__getitem__)
        cols.append(dict(probability=[lo, hi], closest_grid=grid[j], fs_gap_deg=gaps[j]))
    # For columns 1 and 2, unitarity bounds each non-electron entry by their sum.
    cap = max(bounds[0][1], bounds[1][1], 1-bounds[0][0], 1-bounds[1][0],
              bounds[2][1], (1-bounds[2][0])*max(atm[1], 1-atm[0]))
    return dict(columns=cols, minimum_fs_deg=min(x["fs_gap_deg"] for x in cols),
                all_entry_probability_ceiling=cap,
                same_representation_same_reflections_gap_deg=angle(cap)-180/14)


def grid_from_representations(n):
    # The odd-m+even-m character matters: the singleton is the opposite
    # eigenvector to the sign of the singlet character on this reflection.
    vals=set()
    for j in range(1, n//2):
        for chi_r in [1, -1]:
            for chi_s in [1, -1]:
                rays = []
                for m in range(n):
                    a=j*m*math.pi/n
                    if chi_s*chi_r**m==1:
                        a+=math.pi/2
                    rays.append(np.array([math.cos(a), math.sin(a)]))
                vals.update(round(float(x@y)**2, 12) for x in rays for y in rays)
    target=sorted(math.cos(k*math.pi/n)**2 for k in range(n//2+1))
    assert len(vals)==len(target) and np.allclose(sorted(vals),target,atol=1.e-12,rtol=0)
    return target


def competitor(inputs, osc, other):
    a=osc["DayaBay_sin2_2theta13"][0]
    p13=(1-math.sqrt(1-a))/2
    p12=osc["JUNO_sin2_theta12_NO"][0]
    p23=other["NO"]["atmospheric_sine_squared"][0]
    s12,s13,s23=map(math.sqrt,[p12,p13,p23])
    c12,c13,c23=map(math.sqrt,[1-p12,1-p13,1-p23])
    q=math.cos(5*math.pi/14)**2
    cs=(q-(s12*c23)**2-(c12*s23*s13)**2)/(2*s12*c23*c12*s23*s13)
    delta=inputs["delta_sign_witness"]*math.acos(cs)
    z=cmath.exp(1j*delta)
    U=np.array([[c12*c13,s12*c13,s13*z.conjugate()],
                [-s12*c23-c12*s23*s13*z,c12*c23-s12*s23*s13*z,s23*c13],
                [s12*s23-c12*c23*s13*z,-c12*s23-s12*c23*s13*z,c23*c13]])
    I=np.eye(3)
    assert np.allclose(U.conj().T@U,I)
    e_mu=I[:,1]
    nu=U[:,0]
    # On the singlet both reflections act as -1; their fixed lines are mu and nu_1.
    a=2*np.outer(e_mu,e_mu)-I
    b=2*np.outer(nu,nu.conj())-I
    r=b@a
    assert np.allclose(np.linalg.matrix_power(r,14),I)
    assert abs(np.trace(r)-(1+2*math.cos(10*math.pi/14)))<1.e-12
    assert all(np.linalg.norm(np.linalg.matrix_power(r,k)-I)>1.e-4 for k in (1,2,7))
    assert np.allclose(a@r@a,r.conj().T)
    charged=np.diag([.00051,.106,1.777])
    mass=np.diag([.005,.010,.050])
    left=U@mass**2@U.conj().T
    majorana=U.conj()@mass@U.conj().T
    assert np.allclose(a@charged,charged@a)
    assert np.allclose(b@left,left@b)
    assert np.allclose(b.T@majorana@b,majorana)
    assert abs(abs(U[1,0])**2-q)<1.e-14
    arc=read(inputs["phase_check"])["phase_credible_arcs_pi"]["NO_99.73"]
    assert arc[0] < delta/math.pi < arc[1]
    return dict(doublet_j=5, mu_1=q,delta_deg=math.degrees(delta),
                all_moduli_squared=(abs(U)**2).tolist(),
                product_trace=float(np.trace(r).real),
                J=float((U[0,0]*U[1,1]*np.conj(U[0,1]*U[1,0])).imag))


def run():
    cfg=read("inputs.json")
    osc=read(cfg["reactor"])
    other=read(cfg["muon_disappearance"])
    old=read(cfg["older_reactor"])["juno_primary"]
    amp=osc["DayaBay_sin2_2theta13"]
    width=cfg["half_widths"]
    axes=grid_from_representations(cfg["n"])
    boxes={}
    for key, solar, w, order in [
            ("central",osc["JUNO_sin2_theta12_NO"],0,"NO"),
            ("normal_two_width",osc["JUNO_sin2_theta12_NO"],width,"NO"),
            ("historical_JUNO_59d_two_width",[old["sin_sq_theta12"],old["one_sigma"]],width,"NO"),
            ("inverted_historical_solar_prior_two_width",other["IO"]["solar_sine_squared"],width,"IO")]:
        atm=other[order]["atmospheric_sine_squared"]
        atm_bounds=[atm[0]-w*atm[1],atm[0]+w*atm[-1]]
        boxes[key]=certificate(row_bounds(solar,amp,w),axes,atm_bounds)
    return dict(grid_probabilities=axes,certificates=boxes,
                faithful_other_irrep_witness=competitor(cfg,osc,other))


def equal(a,b):
    if isinstance(a,dict):
        return isinstance(b,dict) and a.keys()==b.keys() and all(equal(v,b[k]) for k,v in a.items())
    if isinstance(a,list):
        return isinstance(b,list) and len(a)==len(b) and all(equal(x,y) for x,y in zip(a,b))
    return math.isclose(a,b,rel_tol=1.e-11,abs_tol=2.e-12)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    flags=parser.add_mutually_exclusive_group(required=True)
    flags.add_argument("--write-new",action="store_true")
    flags.add_argument("--check",action="store_true")
    args=parser.parse_args()
    data=run()
    if args.write_new:
        with (HERE/"results.json").open("x",encoding="utf8") as out:
            json.dump(data,out,ensure_ascii=False,allow_nan=False,indent=2)
            out.write("\n")
    elif not equal(data,read("results.json")):
        raise SystemExit("1616 saved results differ")
    print("1616 finite grid, complete electron row, and faithful counterexample verified")
