"""Vacuum-locked, two-path FX phase. Cabibbo and charged-B angles are holdouts."""

import argparse
import importlib.util
import json
import math
import sys
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
sys.dont_write_bytecode = True
spec = importlib.util.spec_from_file_location("historical_fx", HERE.parent / "1412" / "check.py")
fx = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fx)


def read(name):
    return json.loads((HERE/name).read_text(encoding="utf-8"))


def box(parts, n):
    return parts[0]-n*sum(parts[1:]), parts[0]+n*sum(parts[1:])


def angles(u, d, w):
    su, cu = fx.trig(u)
    sd, cd = fx.trig(d)
    s = w / cu
    c = math.sqrt(1 - s*s)
    k = math.hypot(su*cd*c, cu*sd)
    gamma = math.atan2(d*c*(1+u*u), u*(1-d*d*c*c))
    return k, gamma


def ranges(n, cfg):
    tree, belle, masses = read(cfg["tree"]), read(cfg["belle"]), read(cfg["masses"])["yukawa"]
    u, w = box(tree["ub_over_cb"], n), box(belle["BelleII_Vcb_BGL"], n)
    yd, ys = box(masses["down"][0], n), box(masses["down"][1], n)
    d = [math.sqrt(yd[0]/ys[1]), math.sqrt(yd[1]/ys[0])]
    c = [math.sqrt(1-w[1]**2*(1+u[1]**2)), math.sqrt(1-w[0]**2*(1+u[0]**2))]
    assert (1-d[1]**2)/(1+u[1]**2)**2 > w[1]**2
    kvar = [angles(u[0],d[0],w[1])[0], angles(u[1],d[1],w[0])[0]]
    ylo, yhi = d[0]*c[0]*(1+u[0]**2), d[1]*c[1]*(1+u[1]**2)
    xlo, xhi = u[0]*(1-d[1]**2*c[1]**2), u[1]*(1-d[0]**2*c[0]**2)
    gamma = [math.degrees(math.atan2(ylo,xhi)), math.degrees(math.atan2(yhi,xlo))]
    kf, form = box(tree["K_product"], n), box(tree["K_form"], n)
    observed = [kf[0]/form[1], kf[1]/form[0]]
    measured_gamma = read(cfg["charged_gamma"])["chargedB"]
    grace = measured_gamma["additional_degrees"]
    gobs = [measured_gamma["central_degrees"]-n*measured_gamma["lower_error"]-grace,
            measured_gamma["central_degrees"]+n*measured_gamma["upper_error"]+grace]
    cos_min = max(0., (kvar[0]**2-observed[1]**2)/(2*u[1]*d[1]))
    return {"inputs_u_d_w":[list(u),d,list(w)], "Vus_fixed_phase_range":kvar,
            "Vus_measured_range":observed, "gamma_positive_outer_deg":gamma,
            "gamma_conjugate_mod_180_outer_deg":[180-gamma[1],180-gamma[0]],
            "charged_gamma_data_mod_180_deg":gobs,
            "minimum_phase_departure_from_quadrature_deg":math.degrees(math.asin(cos_min)),
            "minimum_real_over_imaginary_source_if_only_phase_moves":cos_min/math.sqrt(1-cos_min*cos_min)}


def witness(cfg, sign):
    charges = cfg["z3"]
    for family, expected in [("Q1", ["Hu","Hd"]), ("Q23",["H0","H0"])]:
        for conj, target in zip([True,False], expected):
            allowed = [h for h in ["H0","Hu","Hd"]
                       if (-charges[family] + (-1 if conj else 1)*charges[h]) % 3 == 0]
            assert allowed == [target]
    widths = cfg["vev_example"]
    wu, wd = widths["wu_GeV"], widths["wd_GeV"]
    w0 = math.sqrt(widths["v_GeV"]**2/2 - wu**2-wd**2)
    phase_u, phase_d = 0.0, -sign*math.pi/2
    xu = widths["xu_GeV"]*np.exp(1j*phase_u)
    xd = widths["xd_GeV"]*np.exp(1j*phase_d)
    a = -sign*widths["abs_A_GeV"]
    assert np.isclose(xu**2, abs(xu)*xu.conjugate())
    assert np.isclose(xu*xd, 1j*abs(xu*xd/a)*a)
    assert np.isclose(w0*wu*np.exp(1j*phase_u), (w0*wu/abs(xu))*xu)
    assert np.isclose(w0*wd*np.exp(1j*phase_d), (w0*wd/abs(xd))*xd)
    point = cfg["witness"]
    u = point["up_tangent"]
    d = math.sqrt(point["down_yukawa"]/point["strange_yukawa"])
    w = point["Vcb"]
    su, cu = fx.trig(u)
    sd, cd = fx.trig(d)
    s, c = w/cu, math.sqrt(1-(w/cu)**2)
    rotation23 = np.array([[1.,0.,0.],[0.,c,s],[0.,-s,c]])
    ou, od = fx.rotation(u), rotation23 @ fx.rotation(d)
    up, down = read(cfg["masses"])["yukawa"]["up"], read(cfg["masses"])["yukawa"]["down"]
    yu = ou @ np.diag([-up[0][0], up[1][0], up[2][0]]) @ ou.T
    yd = od @ np.diag([-point["down_yukawa"], point["strange_yukawa"], down[2][0]]) @ od.T
    assert abs(yd[0,0]) < 1e-14
    pu = np.diag([np.exp(-1j*phase_u),1,1])
    pd = np.diag([np.exp(1j*phase_d),1,1])
    # Real row Yukawas rescaled by their own vevs generate these precise mass matrices.
    norm = widths["v_GeV"]/math.sqrt(2)
    for numeric, ratios, phases, target in [
        (yu, [wu,w0,w0], pu, pu@ou), (yd,[wd,w0,w0],pd,pd@od)
    ]:
        bare = norm * numeric / np.array(ratios)[:,None]
        mass = phases @ np.diag(ratios) @ bare
        diag = target.conjugate().T @ mass @ mass.conjugate().T @ target
        assert np.linalg.norm(diag-np.diag(np.diag(diag))) < 1e-8
    v = (pu@ou).conjugate().T @ (pd@od)
    quartet = (v[0,0]*v[1,1]*v[0,1].conjugate()*v[1,0].conjugate()).imag
    assert np.isclose(quartet, sign*su*cu*sd*cd*s*s*c, atol=1e-15)
    g = np.degrees(np.angle(-v[0,0]*v[0,2].conjugate()/(v[1,0]*v[1,2].conjugate())))
    return {"Vus":abs(v[0,1]), "gamma_DK_mod180_deg":float(g%180), "J":float(quartet)}


def same(a,b):
    if isinstance(a,dict):
        return a.keys()==b.keys() and all(same(x,b[k]) for k,x in a.items())
    if isinstance(a,list):
        return len(a)==len(b) and all(same(x,y) for x,y in zip(a,b))
    return math.isclose(a,b,rel_tol=1e-10,abs_tol=1e-12)


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--write",action="store_true")
    mode.add_argument("--check",action="store_true")
    args = parser.parse_args()
    config = read("inputs.json")
    result = {str(n): ranges(n,config) for n in [1,2]}
    result["witness_positive"] = witness(config,1)
    result["witness_conjugate"] = witness(config,-1)
    positive, two = result["witness_positive"], result["2"]
    assert two["Vus_measured_range"][0] < positive["Vus"] < two["Vus_measured_range"][1]
    assert two["charged_gamma_data_mod_180_deg"][0] < positive["gamma_DK_mod180_deg"] < two["charged_gamma_data_mod_180_deg"][1]
    output = HERE/"results.json"
    if args.write:
        with output.open("x",encoding="utf-8") as out:
            json.dump(result,out,indent=2)
            out.write("\n")
    elif not same(result,read("results.json")):
        raise ValueError("Data or phase certificate differs")
    print(json.dumps(result,indent=2))
