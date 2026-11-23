"""Finite electroweak Majorana self-energy on the actual four-pole branch."""
import argparse
import json
import math
from pathlib import Path
import runpy

import numpy as np


HERE = Path(__file__).resolve().parent
ARCH = HERE.parent


def read(path):
    return json.loads((HERE/path).read_bytes())


def finite(d, c, h, z_mass, higgs, vev, order):
    """Schur resolvent for the NN block, with exact lepton-number protection."""
    i2 = np.eye(2)
    nodes, weights = np.polynomial.legendre.leggauss(order)
    out = np.zeros((2, 2), complex)
    scale2 = np.linalg.norm(c, 2)**2
    for z, count in [(z_mass**2, 3), (higgs**2, 1)]:
        logarithm = math.log1p(scale2/z)
        for x, w in zip(nodes, weights):
            v = logarithm*(x+1)/2
            u = math.expm1(v)*z/scale2
            weight = logarithm/2*math.exp(v)*z/scale2*w
            a = (1-u)*z*i2+u*c.conj()@c.T
            off = u*c.conj()@h
            b = (1-u)*z*i2+u*(c.conj().T@c+h.conj().T@h)
            left = np.linalg.solve(a,off)
            core = -left @ np.linalg.solve(b-off.conj().T@left,c.conj().T)
            out += weight*count*z*core
    return d@out@d.T/(16*math.pi**2*vev**2)


def observables(mass, e):
    # Takagi vectors are eigenvectors of m†m for symmetric m, phases cancel here.
    val, vectors = np.linalg.eigh((mass.conj().T@mass)*1e18)
    u = e.conj().T@vectors
    p = abs(u)**2
    j = np.imag(u[0,0]*u[1,1]*np.conj(u[0,1]*u[1,0]))
    a13 = 4*p[0,2]*(1-p[0,2])
    return dict(nonzero_masses_eV=np.sqrt(val[1:]), R=math.sqrt(val[2]/val[1]),
                A13_solar_atmos=np.array([a13,p[0,1]/(1-p[0,2]),p[1,2]/(1-p[0,2])]),
                massless_flavor_probabilities=p[:,0], J=j)


def unitary(random, n):
    return np.linalg.qr(random.normal(size=(n,n))+1j*random.normal(size=(n,n)))[0]


def compute():
    common = read("../1698/results.json")["shared_witness"]
    inputs = read("../1698/inputs.json")
    old = runpy.run_path(str(ARCH/"1698/check.py"))
    d, c, h = old["neutral_graph"](inputs)
    t,b = common["neutral_t"],common["neutral_b"]
    d,c,h = d*math.sqrt(b),c*(b*math.sqrt(t)),h*(b*t)
    w = read("../1678/independent_results.json")["joint_witness"]
    weak = read("../1397/inputs.json")["high_scale"]
    qs, margins = np.array(weak["down"]+weak["angles_theta12_theta13_theta23_delta"][:3]).T
    port = runpy.run_path(str(ARCH/"1678/independent_check.py"))
    raw = qs+margins*w["quark_pulls"]
    s, l,_,_ = port["get_matrices"](raw,w["quark_gamma_high"])
    ye = (s+w["common_abs_t"]*np.exp(-1j*np.pi/6)*l).T*w["k_yb"]
    e = np.linalg.svd(ye)[0][:,::-1]
    # Remove the holomorphic 1689 family row phase, insert the previously calibrated unitary port.
    transport = np.exp(1j*np.r_[0.,port["ICE"]["phase_legs"]("third_leg_3_to_2")])
    angle = math.pi*read("../1685/inputs.json")["scalar_angle_over_pi"]
    family = read("../1688/inputs.json")["su5_generation_charges_5bar"]
    d = np.diag(transport*np.exp(-1j*angle*np.array(family)))@d
    mh = read("../1372/inputs.json")["atlas"]["mass_GeV"][0]
    mz = read("../1369/inputs.json")["mu_ref_GeV"]
    vev = read("../1492/inputs.json")["v246_GeV"]
    calc = lambda d,c,h,n=32: finite(d,c,h,mz,mh,vev,n)
    z = np.linalg.solve(c,d.T).T
    tree = z@h@z.T
    obs0 = observables(tree,e)
    np.testing.assert_allclose(obs0["A13_solar_atmos"],w["three_neutrino_angles"],rtol=3e-9,atol=0)
    correction = calc(d,c,h)
    ratio_error = np.linalg.norm(correction-calc(d,c,h,64))/np.linalg.norm(correction)
    assert ratio_error<1e-11
    assert np.count_nonzero(calc(d,c,0*h))==0
    np.testing.assert_allclose(calc(d,c,-h),-correction,rtol=2e-14,atol=1e-26)
    cooled = calc(d,c,h*1e-5)
    np.testing.assert_allclose(cooled*1e5,correction,rtol=2e-14,atol=1e-26)
    random = np.random.default_rng(1703)
    active, nr, hr = (unitary(random,n) for n in (3,2,2))
    rotated = calc(active@d@nr,nr.T@c@hr,hr.T@h@hr)
    np.testing.assert_allclose(rotated,active@correction@active.T,rtol=2e-13,atol=1e-26)
    full = tree+correction
    obs1 = observables(full,e)
    obsrot = observables(active@tree@active.T+rotated,active.conj()@e)
    np.testing.assert_allclose(obsrot["A13_solar_atmos"],obs1["A13_solar_atmos"],atol=4e-14)
    # Independently evaluate the spectral function directly from M†M eigenvectors.
    mass = np.block([[np.zeros_like(c), c], [c.T, h]])
    dm = np.column_stack([d,d*0])
    a, modes = np.linalg.eigh(mass.conj().T@mass)
    weights = sum(k*zz*a*np.log(a/zz)/(a-zz) for zz,k in [(mz**2,3),(mh**2,1)])
    rawloop = dm@(modes*weights)@modes.conj().T@np.linalg.solve(mass,dm.T)/(16*np.pi**2*vev**2)
    raw_err = np.linalg.norm(rawloop-correction)/np.linalg.norm(correction)
    assert raw_err<2e-5
    zplus = np.linalg.pinv(z)
    whiten = np.diag(1/np.sqrt(np.diag(h)))
    source = whiten@zplus@correction@zplus.T@whiten
    np.testing.assert_allclose(source,source.T,atol=1e-16)
    # Profile two independent real H-coefficients, retaining the same Yukawa/charged port.
    def profile(logweights):
        hs = h@np.diag(np.exp(logweights))
        return observables(z@hs@z.T+calc(d,c,hs),e)
    shifts = np.zeros(2)
    ref = obs0["nonzero_masses_eV"]
    for _ in range(5):
        y = np.log(profile(shifts)["nonzero_masses_eV"]/ref)
        deriv = np.column_stack([(np.log(profile(shifts+v*1e-5)["nonzero_masses_eV"]/ref)-
                                  np.log(profile(shifts-v*1e-5)["nonzero_masses_eV"]/ref))/2e-5 for v in np.eye(2)])
        shifts -= np.linalg.solve(deriv,y)
    obs2=profile(shifts)
    np.testing.assert_allclose(obs2["nonzero_masses_eV"],ref,rtol=1e-13)
    for rec in [obs1,obs2]:
        np.testing.assert_allclose(rec["massless_flavor_probabilities"],obs0["massless_flavor_probabilities"],atol=3e-14)
    data = read("../1668/inputs.json")
    central,errlo,errhi = data["juno_NO_dm31"]
    factor = ref[0]/obs1["nonzero_masses_eV"][0]
    obs_single = observables(z@(factor*h)@z.T+calc(d,c,factor*h),e)
    dm31 = obs_single["nonzero_masses_eV"][1]**2
    assert abs(obs_single["nonzero_masses_eV"][0]/ref[0]-1)<1e-14
    assert central-errlo<dm31<central+errhi
    angular_box = np.array(read("../1659/results.json")["joint_windows"])
    for rec in [obs1,obs_single,obs2]:
        ang = rec["A13_solar_atmos"]
        assert all((angular_box[:,0]<ang)&(ang<angular_box[:,1]))
    common_fit = np.vdot(tree,correction)/np.vdot(tree,tree)
    # Separate frozen-direction leading-log control. A physical common kinetic
    # operator must be canonically removed in both Yukawa and Majorana sectors.
    strength = math.log(common["common_flavor_cutoff_GeV"]/common["soft_GeV"])/(16*np.pi**2)
    nport = d.conj()*math.sqrt(2)/vev
    metric = np.eye(3)+strength*(nport@nport.conj().T)
    vals, vec = np.linalg.eigh(metric)
    undo = (vec/np.sqrt(vals))@vec.conj().T
    ek = np.linalg.svd(undo@ye)[0][:,::-1]
    with_metric = observables(undo.T@tree@undo,ek)
    null = np.linalg.eigh(tree.conj().T@tree)[1][:,0]
    np.testing.assert_allclose(undo@null,null,rtol=0,atol=2e-13)
    return dict(tree=obs0, finite_neutral_EW_self_energy_only=obs1,
                one_real_common_source_scale_fixing_dm21=factor, held_out_dm31_after_only_common_scale_eV2=dm31,
                fitted_real_H_factors=np.exp(shifts), after_two_real_source_mass_calibrations=obs2,
                relative_real_source_recalibration=np.exp(shifts[1]-shifts[0]),
                normalized_diagonal_source_corrections=np.diag(source).real,
                normalized_mixed_source_correction_real_imag=[source[0,1].real,source[0,1].imag],
                best_common_complex_rescaling_of_tree=[common_fit.real,common_fit.imag],
                remaining_fraction_of_tree_after_best_common_rescaling=np.linalg.norm(correction-common_fit*tree)/np.linalg.norm(tree),
                independent_kinetic_control_gamma=strength,
                independent_kinetic_control_largest_eigenvalue_shift=vals[-1]-1,
                independent_kinetic_control_after_common_canonical_normalization=with_metric,
                raw_spectral_subtraction_relative_error=raw_err, quadrature_relative_change=ratio_error)


def tidy(x):
    if isinstance(x,dict):
        return {k:tidy(v) for k,v in x.items()}
    if isinstance(x,(list,np.ndarray)):
        return [tidy(v) for v in x]
    return float(f"{x:.11g}") if isinstance(x,(float,np.floating)) else x


if __name__=="__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    out = tidy(compute())
    if p.parse_args().check:
        assert out==read("independent_results.json")
        print("1703 independent: finite self-energy, quasi-Dirac cancellation, covariant ports and source refit verified")
    else:
        (HERE/"independent_results.json").write_text(json.dumps(out,indent=2)+"\n",encoding="utf-8")
