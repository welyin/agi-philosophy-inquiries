"""A signed, gauge-covariant matter exchange on the R1712 tagged family frame."""
import argparse
import itertools
import json
import math
from pathlib import Path
import runpy

import numpy as np


HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE / name).read_bytes())


CFG = read("inputs.json")


def actual_target():
    old = runpy.run_path(str(HERE/CFG["mass_shooting"]))
    tag = runpy.run_path(str(HERE/CFG["flavor_transport"]))
    args = read(CFG["mass_case"])["points"][CFG["seed_label"]]
    fl = old["FL"]
    solve = fl["solve"]
    cache = {}
    def capture(*args, **kwargs):
        cache["match"] = solve(*args, **kwargs)
        return cache["match"]
    fl["solve"] = capture
    out = old["run_point"](**args)
    down, up, _ = tag["tagged_matrices"](fl, cache["match"], .20)
    return down, up, out


def norm(a):
    return np.linalg.norm(a, 2)


def polar(a):
    l, _, rh = np.linalg.svd(a)
    return l @ rh


def generated_alignment(y):
    h = y @ y.conj().T
    scale = np.trace(h).real
    h = h / scale  # positive rescaling leaves polar factors and relative gaps unchanged
    eigen = np.linalg.eigvalsh(h)
    e2 = sum(eigen[i]*eigen[j] for i,j in itertools.combinations(range(3),2))
    # The characteristic derivative uses no eigenvector choice or fitted sign table.
    poly = 2*h - 3*h@h - e2*np.eye(3)
    source = y.conj().T @ poly
    l, masses, rh = np.linalg.svd(y)
    expected = (rh.conj().T * [-1.,1.,-1.]) @ l.conj().T
    u = polar(source)
    np.testing.assert_allclose(u, expected, atol=2e-8)
    singular = np.linalg.svd(source, compute_uv=False) * scale**2
    slopes = np.array([np.prod([value-other for k,other in enumerate(eigen) if k!=i])
                       for i,value in enumerate(eigen)])
    eta_window = [float(slopes[1]),float(min(slopes[0],slopes[2]))]
    assert eta_window[0] < 0 < eta_window[1]
    return expected, {"dimensionless_polynomial_source_singular_values_desc": singular.tolist(),
                      "relative_weakest_alignment": singular[-1]/singular[0],
                      "same_degree_trace_squared_admixture_open_interval":eta_window,
                      "polar_error": norm(u-expected),
                      "positive_plain_mass_source_min_return": min(masses)}


def dilation(t, gain):
    c, basis = np.linalg.eigh(t / gain)
    c = np.clip(c, -1., 1.)
    w0 = []
    signs = []
    for index, value in enumerate(c):
        for parity in (1, -1):
            coefficient = (1+parity*value)/2
            if coefficient > 1e-13:
                w0.append(np.eye(3)[index]*math.sqrt(coefficient))
                signs.append(parity)
    w = np.array(w0) @ basis.conj().T
    a = np.diag(signs)
    np.testing.assert_allclose(w.conj().T @ w,np.eye(3),atol=2e-14)
    np.testing.assert_allclose(w.conj().T @ a @ w,t/gain,atol=2e-13)
    return a,w


def hermitian_invroot(a):
    val, u = np.linalg.eigh(a)
    return (u/np.sqrt(val)) @ u.conj().T


def exact_yukawa(a,w,link,kappa,g,split):
    # The singlet mixing has an overall minus sign, chosen in the actual vertex.
    signs = np.diag(a)
    masses = a @ np.diag(1+split*signs)
    right = -kappa * w @ link.conj().T
    heavy_component = -np.linalg.solve(masses, right)
    raw = g*w.conj().T @ heavy_component
    kin = np.eye(3) + heavy_component.conj().T @ heavy_component
    y = raw @ hermitian_invroot(kin)
    light = np.vstack([np.eye(3), heavy_component]) @ hermitian_invroot(kin)
    full = np.hstack([right,masses])
    np.testing.assert_allclose(full @ light,0,atol=1e-14)
    np.testing.assert_allclose(light.conj().T@light,np.eye(3),atol=3e-14)
    np.testing.assert_allclose(g*w.conj().T @ light[3:,:],y,atol=2e-14)
    return y,np.linalg.svd(full,compute_uv=False)


def compute():
    y,u,_ = actual_target()
    link,alignment = generated_alignment(y)
    ll,_,rrh=np.linalg.svd(y)
    surrogate=(ll*np.array([.17,.083,.047]))@rrh
    surrogate_link,_=generated_alignment(surrogate)
    np.testing.assert_allclose(surrogate_link,link,atol=5e-14)
    alignment["same_axes_arbitrary_ordered_reference_error"]=norm(surrogate_link-link)
    t = y @ link
    np.testing.assert_allclose(t,t.conj().T,atol=1e-15)
    mass = np.linalg.svd(t,compute_uv=False)
    gain = mass[0]
    a, w = dilation(t, gain)
    u /= norm(u)
    source = w @ u
    target = a @ source
    assert abs(np.vdot(source,target))<2e-13
    np.testing.assert_allclose(a @ target, source, atol=1e-15)
    # Equal bare masses make the Euclidean propagator proportional to A at every momentum.
    momentum = .23
    flip = np.linalg.solve(momentum**2*np.eye(len(w))+a@a, a)
    np.testing.assert_allclose(w.conj().T @ flip @ w,
                               t/(gain*(1+momentum**2)), atol=3e-13)
    assert abs(np.vdot(u,t@u))<2e-15 and norm(t@u)>0
    defect = np.eye(3)-(t/gain)@(t/gain)
    assert np.linalg.matrix_rank(defect,tol=1e-10)==2 and len(w)==5
    a_extra, _ = dilation(t,1.1*gain)
    assert len(a_extra)==6
    g = CFG["left_messenger_strength"]
    kappa = gain / math.sqrt(g*g-gain*gain)
    match, heavy_masses = exact_yukawa(a,w,link,kappa,g,0.)
    np.testing.assert_allclose(match,y,atol=2e-15)
    # The determinant in the same specified chiral color basis, including all heavy rows.
    # The dimensionless Higgs insertion is arbitrary and cancels from the ratio.
    vev = .1
    chiral = np.block([[np.zeros((3,3)),vev*g*w.conj().T],
                       [-kappa*w@link.conj().T,a]])
    chiral_flipped=chiral.copy()
    chiral_flipped[3:] *= -1.
    light_det = np.linalg.det(vev*y)
    full_ratios = np.array([np.linalg.det(m)/light_det for m in (chiral,chiral_flipped)])
    np.testing.assert_allclose(full_ratios, np.array([-1,1])*(1+kappa*kappa)**1.5,
                               atol=2e-11)
    def leakage(delta):
        output,_ = exact_yukawa(a,w,link,kappa,g,delta)
        return float(np.vdot(u,output@link@u).real/mass[-1])
    roots = []
    for sign in (-1,1):
        lo,hi=0.,.001
        for _ in range(50):
            mid=(lo+hi)/2
            if abs(leakage(sign*mid))>CFG["max_direct_return_over_yd"]:
                hi=mid
            else:
                lo=mid
        roots.append(sign*(lo+hi)/2)
    positive = (1+np.diag(a))/2
    pos_weight = norm(positive*source)**2
    inv_alpha_delta = -2/(3*math.pi*2) * np.log(
        CFG["flavor_boundary_GeV"]/(CFG["mediator_matching_GeV"]*heavy_masses)).sum()
    charge = read(CFG["family_charges"])
    ur = np.mod(-np.array(charge["su5_generation_charges_10"]),12)
    dr = np.mod(-np.array(charge["su5_generation_charges_5bar"]),12)
    residue_ranks=[sum(min(np.count_nonzero((ur+shift)%12==i),
                           np.count_nonzero(dr==i)) for i in range(12)) for shift in range(12)]
    assert max(residue_ranks)<3
    # Explicit negative example: arbitrarily rotate a family link that both
    # generates Y and is inserted into the would-be alignment energy.
    alt = polar(np.array([[1j,2,0],[0,1,1],[1,0,2j]],complex))
    def self_energy(v):
        data=t@v.conj().T
        h=data@data.conj().T
        e2=(np.trace(h)**2-np.trace(h@h))/2
        candidate=data.conj().T@(2*np.trace(h)*h-3*h@h-e2*np.eye(3))
        return float(np.trace(v.conj().T@candidate).real)
    np.testing.assert_allclose(self_energy(alt),self_energy(link),rtol=1e-12,atol=2e-25)
    return {"tagged_yukawa_values_desc": mass.tolist(),"alignment":alignment,
            "canonical_matched_yukawa_error":norm(match-y),"inertia_of_signed_yukawa":[1,2,0],
            "u_positive_weight":pos_weight,"u_total_return_rel_yd":leakage(0),
            "allowed_transition_from_u_norm":norm(t@u),"candidate_nonzero_total_spectrum":True,
            "isometry_dimension_minimum_at_gain_yb":len(w),"minimum_at_gain_1p1yb":len(a_extra),
            "kernel_compression_defect_rank":int(np.linalg.matrix_rank(defect,tol=1e-10)),
            "kappa":kappa,"physical_heavy_singular_masses_over_M":heavy_masses.tolist(),
            "full_down_determinant_over_light_for_two_chiral_color_boundaries":full_ratios.real.tolist(),
            "number_of_negative_bare_color_masses": int(np.sum(np.diag(a)<0)),
            "split_at_one_percent_return_near_zero":roots,"leakage_yd_at_signed_split_1e_minus4":leakage(1e-4),
            "extra_dirac_triplets_one_loop_inverse_alpha3_10PeV_shift":inv_alpha_delta,
            "same_charge_RH_current_max_rank_per_shift":list(map(int,residue_ranks)),
            "self_sourced_alignment_energy_two_links":[self_energy(link),self_energy(alt)]}


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args=parser.parse_args()
    dest=HERE/"results.json"
    results=compute()
    if args.check:
        assert json.loads(dest.read_bytes())==results
        print("1713 signed alignment, canonical Schur exchange and source-degeneracy checks passed")
    else:
        dest.write_text(json.dumps(results,indent=2,allow_nan=False)+"\n",encoding="utf-8")
