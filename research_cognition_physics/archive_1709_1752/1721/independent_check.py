"""Test the identification of light d^c and L through the original family ports."""
import argparse
import json
import math
from pathlib import Path
import runpy

import numpy as np


HERE = Path(__file__).resolve().parent


def read(relative):
    return json.loads((HERE/relative).read_bytes())


def calculate():
    prev = read("../1720/results.json")
    model = prev["central_conditional_combined_point"]
    R = model["R"]
    a0 = 1/R
    kappa = math.sqrt(R**2-1)
    color = runpy.run_path(str(HERE/"../1713/check.py"))
    y, _, _ = color["actual_target"]()
    link, _ = color["generated_alignment"](y)
    vals = np.linalg.svd(y,compute_uv=False)
    a,w = color["dilation"](y@link, vals[0])
    embed = a@w@link.conj().T
    d = np.vstack((np.eye(3), kappa*embed))/R
    hrow = np.column_stack((-kappa*w@link.conj().T, a))
    np.testing.assert_allclose(hrow@d, 0, atol=1e-14)
    np.testing.assert_allclose(d.conj().T@d, np.eye(3), atol=1e-14)
    bare = model["five_color_poles_GeV"][0]
    vev = read("../1492/inputs.json")["v246_GeV"]/math.sqrt(2)
    x = model["g"]*vev
    mat = np.block([[np.zeros((3,3)),x*w.conj().T],[bare*hrow[:,:3],bare*a]])
    _, singular, rh = np.linalg.svd(mat)
    right = rh.conj().T[:,::-1][:,:3]
    projections = np.linalg.norm(right[:3],axis=0)
    k2 = (vals[::-1]/vals[0])**2
    q = (singular[::-1][:3]/bare)**2
    x2 = (x/bare)**2
    predicted = 1/np.sqrt(1+kappa**2*((1-k2)/(1-q)**2+k2/(1+x2-q)**2))
    np.testing.assert_allclose(projections, predicted, rtol=1e-10, atol=1e-14)
    # Explicit optimal lepton embeddings for all four resource bounds.
    capacities = []
    for n in range(4):
        L = np.vstack((np.eye(3),np.zeros((5,3),complex)))
        L[:,:n] = d[:,:n]
        np.testing.assert_allclose(L.conj().T@L, np.eye(3), atol=1e-14)
        port = L.conj().T@d
        probabilities = np.linalg.svd(port,compute_uv=False)**2
        np.testing.assert_allclose(probabilities, [1.]*n+[a0**2]*(3-n), atol=1e-14)
        capacities.append(float(sum(probabilities)))
    # Tree mass [−κM I, M I] on the weak partner of each aligned light line.
    weak_mass = bare*np.column_stack((-kappa*np.eye(3),np.eye(3)))
    np.testing.assert_allclose(np.linalg.svd(weak_mass,compute_uv=False), bare*R)
    U = prev["baseline_unification_GeV"]
    extra_Bw = 2*math.log(U/(bare*R))
    weak_new = np.array(model["four_weak_poles_GeV"])*math.exp(extra_Bw/4)
    assert max(weak_new) < U
    np.testing.assert_allclose(np.log(U/weak_new).sum()+extra_Bw,
                               np.log(U/np.array(model["four_weak_poles_GeV"])).sum(), atol=5e-14)
    tau = read("../1715/independent_results.json")["charged_yukawas_high"][-1]
    G = math.sqrt(4*math.pi)
    min_original_tau = tau/G
    max_tau_current = a0*min_original_tau+math.sqrt((1-a0**2)*(1-min_original_tau**2))
    return {"R":R, "squared_original_port_on_each_light_right_line":a0**2,
            "finite_vev_original_amplitudes":projections.tolist(),
            "sum_current_strength_max_0_1_2_3_new_weak_matter_directions":capacities,
            "isolated_X_exchange_neutrino_mode_rate_relative_to_unsplit_1698":
                (read("../1698/results.json")["shared_witness"]["inverse_alpha_U"] /
                 model["common_inverse_coupling"]/R)**2,
            "three_matter_weak_Dirac_poles_each_GeV":bare*R,
            "Higgs_transfer_doublet_masses_after_admitting_3_weak_matter_pairs_GeV":weak_new.tolist(),
            "original_only_tau_yukawa_needed_for_exact_common_port":R*tau,
            "original_only_tau_yukawa_at_minimum_R_of_fixed_weak_shape":
                prev["at_fixed_baseline_R_min_from_lightest_four_weak_poles"]*tau,
            "original_only_yukawa_budget":G,"minimum_original_tau_amplitude_at_that_budget":min_original_tau,
            "max_current_amplitude_on_tau_for_original_only_at_that_budget":max_tau_current}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    data=json.loads(json.dumps(calculate()),parse_float=lambda x:float(f"{float(x):.10g}"))
    out=HERE/"independent_results.json"
    if args.check:
        assert data==json.loads(out.read_bytes())
        print("1721 verified actual right-family projection, sharp ranks, matter matching and tau constraint")
    else:
        with out.open("x",encoding="utf8") as stream:
            stream.write(json.dumps(data,indent=2,allow_nan=False)+"\n")
