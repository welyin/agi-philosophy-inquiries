"""Round 1550: norm, physical Z transfer, and a three-family shared-matter witness."""
import argparse
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent


def read(p):
    return json.loads(p.read_text(encoding="utf-8-sig"))


def mixing(angles, phase):
    s, t, u = np.sin(angles)
    c, v, w = np.cos(angles)
    z = np.exp(1j*phase)
    return np.array([[c*v,s*v,t/z],[-s*w-c*u*t*z,c*w-s*u*t*z,u*v],
                     [s*u-c*w*t*z,-c*u-s*w*t*z,w*v]])


def bounds(cfg, charged):
    top, dt = cfg["up_at_1TeV"][-1]
    yb, db = charged[2]
    tau, de = charged[5]
    ratio = 1-(yb+2*db+tau+2*de)/(top-2*dt)
    v = cfg["illustrative_common_v_GeV"]
    direct = cfg["CMS_high_bin"]
    top_low = direct["mt_GeV"]-2*sum(direct["down_variations_GeV"])
    primary = v/math.sqrt(2)*(top-2*dt)*ratio
    secondary = direct["min_total_transfer"]*top_low*ratio
    mz2 = cfg["Zmass_GeV"]/2
    nlo = cfg["LEP_Nnu"][0]-2*cfg["LEP_Nnu"][1]
    mlo = cfg["light_mass_max_eV"]*1e-9
    rows = {}
    for xi in cfg["N_to_extracted_upper_factors"]:
        delta = 1-math.sqrt(nlo/xi-2)
        rows[str(xi)] = {"largest_loss_upper":delta,
            "heavy_max_lower_standard_matching_GeV":math.sqrt(mlo**2+(primary**2-mlo**2)/delta),
            "heavy_max_lower_CMS_matching_GeV":math.sqrt(mlo**2+(secondary**2-mlo**2)/delta)}
    return {"Dirac_over_top_lower":ratio,"table_singular_lower_GeV":primary,
            "CMS_top_working_lower_GeV":top_low,"CMS_singular_lower_GeV":secondary,
            "CMS_combined_transfer_critical_for_Z_threshold":mz2/(top_low*ratio),
            "LEP_two_width_lower":nlo,"necessary_transfer_for_only_two_species":nlo/2,
            "Majorana_conditional":rows}


def construction(cfg, charged):
    scenario = cfg["neutral_witness"]
    v = cfg["illustrative_common_v_GeV"]/math.sqrt(2)
    up = v*np.array(cfg["up_at_1TeV"])[:,0]
    down, lepton = (v*np.array(charged)[:,0]).reshape(2,3)
    ckm = cfg["ckm_theta12_theta13_theta23_delta"]
    V = mixing(ckm[:3],ckm[3])
    A = V.conj().T@np.diag(up)@V + np.diag(down)
    eps = scenario["auxiliary_epsilon_GeV"]
    eig, vec = np.linalg.eigh(A-eps*np.eye(3))
    assert lepton[0]<eig[0]<lepton[1]<eig[1]<lepton[2]<eig[2]
    residues = [np.prod(a-lepton)/np.prod(a-np.delete(eig,i)) for i,a in enumerate(eig)]
    t = vec@np.sqrt(residues)
    D = np.outer(t,t.conj())+eps*np.eye(3)
    spectrum, lep_axes = np.linalg.eigh(A-D)
    assert np.allclose(spectrum,lepton,atol=1e-10)
    # Rephase the charged-lepton eigenbasis and use a real singular basis for D.
    z = abs(lep_axes.conj().T@t)/np.linalg.norm(t)
    axis = z-np.eye(3)[:,2]
    O = np.eye(3)-2*np.outer(axis,axis)/(axis@axis)
    assert np.allclose(O[:,2],z)
    diag_d = np.diag([eps,eps,np.vdot(t,t).real+eps])

    osc = read(HERE/cfg["neutrino"])["NO"]
    s23 = read(HERE/cfg["theta23"])["NO"]["atmospheric_sine_squared"][0]
    m0 = scenario["lightest_eV"]
    masses = 1e-9*np.sqrt(m0**2+np.array([0,osc["solar"][0],osc["atm_abs"][0]]))
    U = mixing(np.arcsin(np.sqrt([osc["s12"][0],osc["s13"][0],s23])),
               scenario["Dirac_phase_rad"])
    U = np.diag(np.exp(1j*np.array(scenario["row_phase_rad"])))@U
    p = (z@U)**2/masses
    l1, l2, l3 = abs(p)
    assert 2*max(abs(p))<sum(abs(p))
    theta = math.acos((l3**2-l1**2-l2**2)/(2*l1*l2))
    target = np.array([1,np.exp(1j*theta),-(l1+l2*np.exp(1j*theta))/l3])
    phases = target/(p/abs(p))
    assert abs(p@phases) < 1e-12*max(abs(p))
    K = O.T@((U*(phases/masses))@U.T)@O
    assert abs(K[2,2]) < 1e-12*np.linalg.norm(K)
    K[2,2] = 0  # exact triangle identity, avoid amplified double-precision remainder
    K = (K+K.T)/2
    R = -diag_d@K@diag_d
    big = np.sort(np.linalg.svd(R,compute_uv=False))
    eff = np.linalg.inv(K)
    tiny = np.sort(np.linalg.svd(eff,compute_uv=False))
    assert np.allclose(tiny,masses,rtol=1e-12,atol=1e-25)
    eta = np.linalg.norm(np.linalg.solve(R,diag_d),2)
    relative = (1+eta)**2-1
    upper_m = max(masses)*(1+relative)
    gap = min(np.diff(masses))-upper_m*relative
    assert gap>0 and big[0]>cfg["Zmass_GeV"] and big[2]<1e9
    # Takagi residual projected to the light sector, then a Hermitian dilation gap.
    vector_error = upper_m*eta**2/((1-upper_m/big[0])*gap)
    heavy_loss_upper = 3*(eta/(1-upper_m/big[0]))**2
    return {"charged_lepton_singular_GeV":spectrum.tolist(),"interlacing_GeV":eig.tolist(),
            "D_singular_GeV":np.diag(diag_d).tolist(),"dominant_axis_abs_charged_basis":z.tolist(),
            "triangle_sides_inverse_eV":(abs(p)*1e-9).tolist(),
            "free_column_phase_angles_rad":np.angle(phases).tolist(),
            "target_light_eV":(masses*1e9).tolist(),"right_block_singular_GeV":big.tolist(),
            "all_physical_singular_relative_error_bound":relative,
            "light_takagi_projective_direction_error_bound":vector_error,
            "largest_loss_upper_from_singular_vectors":heavy_loss_upper,
            "block_mixing_norm":eta,
            "largest_mass_without_phase_cancellation_scale_GeV":diag_d[2,2]**2/max(masses)}


def stable(v):
    if isinstance(v,dict):
        return {k:stable(x) for k,x in v.items()}
    if isinstance(v,list):
        return [stable(x) for x in v]
    return float(format(float(v),".9g"))


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument("--write-new",action="store_true")
    mode.add_argument("--check",action="store_true")
    args=parser.parse_args()
    cfg=read(HERE/"inputs.json")
    charged=read(HERE/cfg["down_lepton"])["light_scales"]["1000"]
    out=stable({"bounds":bounds(cfg,charged),"joint_tree_witness":construction(cfg,charged)})
    path=HERE/"results.json"
    if args.write_new:
        with path.open("x",encoding="utf-8") as stream:
            json.dump(out,stream,indent=2)
            stream.write("\n")
    elif args.check:
        assert read(path)==out, "Stored result differs"
        print("1550: inclusive Dirac tests, heavy-mass product and same-matter phase witness verified.")
    else:
        print(json.dumps(out,indent=2))
