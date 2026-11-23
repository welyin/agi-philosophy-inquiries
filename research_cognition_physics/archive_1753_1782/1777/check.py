"""Coherent J=0 Yukawa amplitude, two-energy subtraction and finite weak spectra."""
import argparse
import json
import math
import runpy
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE / path).read_bytes())


def unitary(rng, n):
    z = rng.normal(size=(n,n))+1j*rng.normal(size=(n,n))
    return np.linalg.qr(z)[0]


def channel_verify(rng, nc, rows, columns):
    # Independent colour contraction before constructing the normalized flavour states.
    y = rng.normal(size=(rows,columns))+1j*rng.normal(size=(rows,columns))
    vector = y.ravel()
    projector = np.outer(vector, vector.conj())
    amplitude = -sum(projector for _ in range(nc*nc))/(16*math.pi*nc)
    strength = nc*np.linalg.norm(y, "fro")**2/(16*math.pi)
    bright = vector/np.linalg.norm(vector)
    probe = rng.normal(size=vector.size)+1j*rng.normal(size=vector.size)
    dark = probe-bright*np.vdot(bright,probe)
    dark /= np.linalg.norm(dark)
    incoming, outgoing = (bright+dark)/math.sqrt(2), (bright-dark)/math.sqrt(2)
    np.testing.assert_allclose(np.linalg.eigvalsh(amplitude)[0], -strength, rtol=1e-12)
    eye = np.eye(vector.size)
    for c in (-50., .015, 3.):
        measured = np.vdot(outgoing,(amplitude+c*eye)@incoming)
        np.testing.assert_allclose(measured,-strength/2,atol=5e-13)
    # A completely general Hermitian contact coefficient between all the same external states.
    raw = rng.normal(size=eye.shape)+1j*rng.normal(size=eye.shape)
    contact=(raw+raw.conj().T)/2
    vals = [np.vdot(outgoing,(amplitude+c*eye+s*contact)@incoming)
            for s,c in [(1.,-12.),(4.,28.)]]
    np.testing.assert_allclose((4*vals[0]-vals[1])/3,-strength/2,atol=8e-13)
    rotated = unitary(rng,rows).conj().T@y@unitary(rng,columns)
    np.testing.assert_allclose(np.linalg.norm(rotated, "fro")**2, strength*16*math.pi/nc)
    return float(strength)


def source_masses(cfg):
    old=read(cfg["mz_vev"])
    up=read(cfg["mz_up"])["MZ"]
    light=read(cfg["mz_down"])
    c=light["points"]["MZ"]["central"]
    h=read(cfg["mz_heavy"])
    v=old["v_MZ_GeV"]
    return v, dict(up=np.array([up["up"][0],up["charm"][0],h["yt"][0]])*v/math.sqrt(2),
                down=np.array([c[0]*light["units"][0],c[1]*light["units"][1],h["yb"][0]])*v/math.sqrt(2))


def hidden_schur():
    # Row labels (Q_L, Q'_L, S'_L), column labels (s_R, S'_R, Q'_R).
    # Entries quoted in common GeV units; the two hard vector masses exceed vev-induced entries.
    direct_a, direct_b, reverse=250.,330.,280.
    mq,ms=2500.,1700.
    m=np.array([[0,direct_a,0],[direct_b,0,mq],[0,ms,reverse]])
    eliminated=m[:1,:1]-m[:1,1:]@np.linalg.solve(m[1:,1:],m[1:,:1])
    np.testing.assert_allclose(eliminated[0,0],direct_a*direct_b*reverse/(mq*ms))
    return dict(open_zero_vev_Yukawa=0.,static_light_GeV=eliminated[0,0],
                physical_light_GeV=min(np.linalg.svd(m,compute_uv=False)))


def calculate():
    cfg=read("inputs.json")
    old=read(cfg["nearby"])
    previous=read(cfg["historical_mass_witness"])
    optimize=runpy.run_path(str(HERE/cfg["old_spectral"]))["optimized_floor"]
    v=read(old["old_width_and_vev"])["v_GeV"]
    nc=cfg["colour"]
    sratio=cfg["energy_ratio_Q2_over_Q1"]**2
    intercept_cost=(sratio+1)/(sratio-1)
    budget2={"Yukawa_block_only":8*math.pi/nc,
             "universal_gauge_transition":16*math.pi/nc,
             "universal_gauge_plus_arbitrary_dim6_two_energies":16*math.pi*intercept_cost/nc}
    v_mz, known=source_masses(cfg)
    results={}
    for name,y2 in budget2.items():
        delta=v*math.sqrt(y2/2)
        charged={}
        for label,masses in known.items():
            available=v_mz**2*y2/2-np.sum(masses**2)
            charged[label]={}
            for p in cfg["reference_family_counts"]:
                value=math.sqrt(available/(p-3))
                assert value>max(masses)
                charged[label][str(p)]=value
        results[name]=dict(maximum_Frobenius=math.sqrt(y2),bound_weak_mass_GeV=delta,
                          MZ_lightest_extra_GeV_if_transport_is_identity=charged,
                          weakest_new_channel_virtuality_GeV=math.sqrt(8)*delta,
                          normalized_extra_weak_spectra={str(q):optimize(q,delta,1,nc)["relative_rate"]
                                      for q in cfg["weak_readout_GeV"]})
    strong_y=previous["required_projected_Y_if_only_three_below"]["7171.07758"]["minimum_Y"]
    coherent_y=previous["entrywise_budget_counterexample"]["operator_Yukawa"]
    large={}
    for name,y in [("all_new_above_7171_GeV",strong_y),("six_by_six_real_example",coherent_y)]:
        amp=nc*y*y/(32*math.pi)
        large[name]=dict(orthogonal_flavor_partial_wave=amp,
                         one_energy_minimum_nonuniversal_cancellation=max(0.,amp-.5),
                         minimum_energy_ratio_for_two_energy_obstruction=math.sqrt((amp+.5)/(amp-.5)))
    eta=cfg["response_error_example"]
    column_tail=cfg["unseen_right_Frobenius_example"]
    drift=cfg["transport_Frobenius_example"]
    y_with_tail=math.sqrt(budget2["universal_gauge_transition"]*(1+2*eta)+column_tail**2)+drift
    compare=cfg["energy_ratio_Q2_over_Q1"]
    # Scalar contact attaining the two-energy conservative norm bound at the endpoints.
    intercept=.5*intercept_cost
    slope=-1/(sratio-1)
    np.testing.assert_allclose([intercept+slope,intercept+sratio*slope],[.5,-.5])
    # Degree two already recovers arbitrary intercept at BOTH fixed positive s.
    x=51.
    poly=x*np.poly1d([1.,-1.])*np.poly1d([1./sratio,-1.])
    np.testing.assert_allclose([poly(0),poly(1.),poly(sratio)],[x,0,0],atol=1e-12)
    rng=np.random.default_rng(cfg["tests_seed"])
    tests=[]
    for i in range(48):
        tests.append(channel_verify(rng,nc,2+i%4,2+i%3))
    return dict(one_doublet_bounds=results,known_retained_MZ_masses_GeV=known,
                required_escape=large,
                residual_and_partial_coverage_example=dict(absolute_partial_wave_error=eta,
                    unprobed_right_tail_F=column_tail,scale_drift_F=drift,
                    projected_F_upper=y_with_tail,tree_mass_upper_GeV=v*y_with_tail/math.sqrt(2)),
                dim6_optimal_endpoints=dict(amplitudes=[.5,-.5],ratio_Q=compare,
                    full_rational_intercept_bound=intercept,dim8_two_zero_control=poly.c),
                positive_closed_threshold_control=hidden_schur(),
                algebra_checks=dict(random_colour_and_family_blocks=len(tests)))


def clean(x):
    if isinstance(x,dict):
        return {k:clean(v) for k,v in x.items()}
    if isinstance(x,(list,np.ndarray)):
        return [clean(v) for v in x]
    if isinstance(x,int):
        return x
    return float(f"{x:.10g}")


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    check=parser.parse_args().check
    data=clean(calculate())
    output=HERE/"results.json"
    if check:
        assert data==read("results.json")
        print("1777 OK: colour singlet flavour amplitudes, two-energy contact subtraction,"
              " thresholds and time-like weak spectral bounds")
    else:
        output.write_text(json.dumps(data,indent=2,allow_nan=False)+"\n",encoding="utf-8")

