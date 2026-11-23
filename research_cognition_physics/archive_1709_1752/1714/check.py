"""Tree-level doublet transmission of a Cabibbo amplitude through two dynamical tags."""
from pathlib import Path
import argparse
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE/path).read_bytes())


def chain(u, s, masses, shortcut):
    row = np.zeros((4,5),complex)
    row[:3,0] = -s
    row[:3,1:4] = np.eye(3)
    row[3,0] = shortcut
    row[3,1:4] = -u.conj()
    row[3,4] = 1.
    row *= np.array([masses[0]]*3+[masses[1]])[:,None]
    _, sm, v = np.linalg.svd(row,full_matrices=True)
    light = v[-1].conj().copy()
    light *= light[0].conjugate()/abs(light[0])
    z = np.vdot(u,s)-shortcut
    expected = np.r_[1.,s,z]/math.sqrt(2+abs(z)**2)
    np.testing.assert_allclose(light,expected,atol=1e-13)
    gram = row @ row.conj().T
    prod = np.prod(sm)**2/(masses[0]**6*masses[1]**2)
    np.testing.assert_allclose(prod,2+abs(z)**2,rtol=3e-14)
    shape = 2*abs(z)**2
    k = (4-16*shape+121*shape**2)/(6*shape)
    return {"normalized_H_amplitude":abs(light[0]), "normalized_B_amplitude":abs(light[-1]),
            "B_over_H":[z.real,z.imag],"Dirac_ratio_squared":shape,
            "tree_mass_ratio":(k+math.sqrt(k*k-4))/2,
            "heavy_masses_GeV":sorted(sm.tolist()),"product_ratio":prod,
            "spectrum_only_p_equal_mass_zero_contact":6**4*math.prod(sm**2)/np.trace(gram).real**4-2}


def main():
    cfg = read("inputs.json")
    point = read(cfg["common_point"])["cases"][cfg["point"]]
    # The primary point supplies the physical down singular values and first up-row weights.
    weights = np.array(point["upper_probabilities"])
    d,s,b = np.array(point["upper_yukawas"])**2
    c = (d+s+b)/3
    sq = (np.array([d,s,b])-c)**2
    gaps = (s-d)*(2*b-s-d)/3, (b-s)*(s+b-2*d)/3
    np.testing.assert_allclose(sq[[0,2]]-sq[1],gaps,rtol=2e-11)
    v = np.sqrt(weights/weights.sum())
    # In the down mass basis, it is legitimate to make the up first row positive.
    target = np.eye(3)[1]
    z = np.dot(v,target)
    np.testing.assert_allclose(z,point["quark_neutrino"]["Vus_at_upper_neutral_cluster"],rtol=1e-10)
    k_data = read("../1702/inputs.json")["direct_semileptonic_moduli"]["Vus"]
    width = np.linalg.norm(k_data[1:])
    length = np.linalg.norm(v[:2])
    angle = math.asin((z+width)/length)-math.asin(z/length)
    rot = np.eye(3)
    rot[:2,:2] = [[math.cos(angle),math.sin(angle)],[-math.sin(angle),math.cos(angle)]]
    bias = np.linalg.norm(rot@np.diag(sq)@rot.T-np.diag(sq),2)
    np.testing.assert_allclose(v @ (rot@target),z+width,atol=1e-13)
    np.testing.assert_allclose(bias,gaps[0]*math.sin(angle),rtol=2e-10)
    masses = cfg["example_heavy_doublet_mass_units_GeV"]
    baseline = chain(v,target,masses,0.)
    control = chain(v,target,masses,cfg["example_direct_mixing"])
    q,_ = np.linalg.qr(np.array([[1+1j,2,3],[2j,-1,1],[1,1,1j]]))
    rotated = chain(q@v,q@target,masses,0.)
    np.testing.assert_allclose(baseline["heavy_masses_GeV"],rotated["heavy_masses_GeV"],rtol=1e-14)
    np.testing.assert_allclose(baseline["spectrum_only_p_equal_mass_zero_contact"],z*z,atol=2e-14)
    # Independent hierarchy check: universal trace position picks the middle state.
    count=0
    for i in range(1,8):
        for j in range(i+1,14):
            for k in range(j+1,20):
                spectrum=np.array([i,j,k])**2
                assert np.argmin((spectrum-spectrum.mean())**2)==1
                count+=1
    align = read(cfg["sterile_alignment"])
    couplings = read(cfg["sterile_couplings"])
    vac = read(cfg["flavor_vacuum"])
    a,bv = np.array(align["source_a"]),np.array(align["source_b"])
    rho = vac["scalar_magnitude_over_cutoff"]
    flav = rho*np.exp(1j*np.pi*vac["scalar_angle_over_pi"])
    za = np.array([0,flav**2*vac["real_X_over_cutoff"]/rho**2])
    zb = np.array([1,couplings["coupling_to_s2"]/couplings["m_N0H0"]*flav**2])
    A = np.outer(a,za)/math.sqrt(np.dot(a,a))
    B = np.outer(bv,zb)
    cross = np.vdot(B,A)
    assert abs(cross.imag)<1e-14 and cross.real>0
    old = read(cfg["old_gauge_and_soft_spectrum"])["shared_witness"]
    gauge_log = sum(math.log(old["U_GeV"]/m) for m in baseline["heavy_masses_GeV"])/(2*math.pi)
    inv = old["inverse_alpha_U"]
    return {"Vus_static":z,"normalized_real_anchor_value":math.sqrt(z),
            "squared_singular_values_down_at_7171GeV":[d,s,b],"traceless_shift":c,
            "minimum_alignment_gaps_down":gaps,"enumerated_nondegenerate_spectra":count,
            "one_K_width_alignment_rotation_rad":angle,"one_K_width_alignment_operator_perturbation":bias,
            "same_1712_point":baseline,"allowed_direct_contact":control,
            "cross_spurion_in_units_of_shared_bare_Yukawa_squared":cross.real,
            "normalized_two_arm_overlap":abs(cross)/(np.linalg.norm(A)*np.linalg.norm(B)),
            "relative_alpha_inverse_changes_at_old_unification":[-.6*gauge_log,-gauge_log,0.],
            "old_unification_GeV":old["U_GeV"],"old_inverse_alpha":inv,
            "inverse_alpha_if_completed_with_equal_mass_color_partners":inv-gauge_log,
            "separate_nonholomorphic_shape_and_matching":True}


def tidy(x):
    if isinstance(x,dict): return {k:tidy(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)): return [tidy(v) for v in x]
    if isinstance(x,(float,np.floating)): return float(f"{x:.11g}")
    return x


if __name__=="__main__":
    args=argparse.ArgumentParser()
    args.add_argument("--check",action="store_true")
    out=tidy(main())
    if args.parse_args().check:
        assert read("results.json")==out
        print("1714 verified")
    else:
        (HERE/"results.json").write_text(json.dumps(out,indent=2)+"\n",encoding="utf8")
