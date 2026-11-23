"""Test the charged-lepton / rank-two Majorana sum rule at the same matching surface."""
import argparse
import json
import math
from pathlib import Path
import runpy

import numpy as np

HERE = Path(__file__).resolve().parent
EYE = np.eye(3, dtype=complex)


def parent(step):
    prior = runpy.run_path(str(HERE.parent/"1739"/"check.py"))
    params = prior["CFG"]["nearby_measured_inputs_in_marginal_widths"]["absolute_and_oscillations"]
    args = (params[k] for k in ("weak", "up", "g3", "top", "charm"))
    _, pt, nu = prior["charged_endpoint"](2, step, *args)
    flow = pt["flow"]
    y = prior["initial_state"](pt)
    transport = EYE.copy()
    poles = flow["fields"]["neutral_poles_GeV"]
    for start, susy, dt, grid in flow["segments"](step):
        if start >= poles[2]:
            break
        assert not susy
        for j, m in enumerate(poles):
            if start == m:
                y[3,:,j] = flow["physics"]["Y0"][:,j]
        y, transport = prior["advance"](y, transport, dt, grid, flow)
    to_low = np.linalg.inv(transport @ pt["low_E"])
    neutral_low = pt["C"] + nu["FROZEN"][0]*np.outer(pt["b"],pt["b"])
    neutral = to_low.T @ neutral_low @ to_low
    c = abs(np.vdot(np.linalg.svd(y[0,:,:3])[0][:,-1], np.linalg.svd(y[1,:,:3])[0][:,-1]))
    return y[2,:,:3], neutral, np.array([c, 1., math.sqrt((1+c)/2)]), float(poles[2])


def observable(charged, neutral):
    v, lepton, _ = np.linalg.svd(charged)
    _, mass, vh = np.linalg.svd(neutral)
    u0 = vh.conj().T[:,-1]
    assert mass[2] < 2e-12 * mass[1]
    prob = abs(v.conj().T @ u0)**2
    moment = np.dot(lepton**2, prob)
    np.testing.assert_allclose(moment, np.linalg.norm(charged.conj().T @ u0)**2,rtol=1e-12)
    charged_volume = np.prod(lepton)
    area = mass[0]*mass[1]
    inverse = np.linalg.inv(charged.conj().T)
    pulled_back = inverse.T @ neutral @ inverse
    singular = np.linalg.svd(pulled_back,compute_uv=False)
    return dict(area=area, moment=moment, volume=charged_volume,
                null=u0, joint=area*moment/charged_volume**2,
                shape=singular[0]/singular[1], probabilities=prob[::-1])


def compare_row(charged, neutral, f, original):
    new = observable(f@charged, f.T@neutral@f)
    norm = np.linalg.norm(np.linalg.solve(f,original["null"]))**2
    det = abs(np.linalg.det(f))
    ratios = {k: new[k]/original[k] for k in ("area","moment","volume","joint","shape")}
    np.testing.assert_allclose([ratios["area"],ratios["moment"],ratios["volume"],ratios["joint"],ratios["shape"]],
                               [det**2*norm,1/norm,det,1.,1.],rtol=3e-9)
    return dict(**ratios, null_inverse_norm_squared=norm,
                PMNS_massless_column_probabilities=new["probabilities"].tolist())


def majorana_test(charged, neutral, inherited_shape):
    """Erase the unmeasured relative Majorana phase; reconstruct its two possible readings."""
    left, lepton, _ = np.linalg.svd(charged)
    left, lepton = left[:,::-1], lepton[::-1]
    n = left.T @ neutral @ left
    nu, mass, _ = np.linalg.svd(n)
    v = nu[:,[1,0]]
    v = v*(v[0].conj()/abs(v[0]))
    cols = v/lepton[:,None]
    gram = cols.conj().T @ cols
    area = sum(abs(np.linalg.det(cols[[i,j],:]))**2 for i,j in ((0,1),(0,2),(1,2)))
    mass_ratio = mass[0]/mass[1]
    constant = gram[0,0].real**2 + (mass_ratio*gram[1,1].real)**2
    interference = 2*mass_ratio*abs(gram[0,1])**2
    def shape_from_trace(trace):
        h = trace/(mass_ratio*area)
        return (h+math.sqrt(h*h-4))/2
    extremes = [shape_from_trace(constant+x*interference) for x in (-1,1)]
    cosval = ((inherited_shape+1/inherited_shape)*mass_ratio*area - constant)/interference
    assert -1-1e-10 < cosval < 1+1e-10
    angles = np.array([-1,1])*math.acos(np.clip(cosval,-1,1))-2*np.angle(gram[0,1])
    beta = abs(abs(v[0,0])**2 + mass_ratio*abs(v[0,1])**2*np.exp(1j*angles))
    actual = abs(n[0,0])/mass[1]
    assert min(abs(beta-actual)) < 2e-8
    return dict(shape_range_with_free_Majorana_phase=extremes, cosine=cosval,
                m_beta_beta_over_m2_candidates=sorted(beta.tolist()), realized=actual)


def gauge_matching(pole):
    folder = HERE.parent/"1739"
    old = json.loads((folder/"gauge_results.json").read_bytes())["R1698_two_new_colours_and_soft_scalars_at_g3_minus_0p67"]
    soft = json.loads((folder/"higgs_results.json").read_bytes())["absolute_and_oscillations"]["soft_GeV"]
    unification, needed = old["U"], old["required_weak_log_load"]
    def weight(mass):
        return 2*(2*math.log(unification/mass)+math.log(unification/math.hypot(mass,soft)))/3
    residual = needed-weight(pole)
    assert residual < 0
    lo,hi = pole,2*pole
    for _ in range(64):
        mid = (lo+hi)/2
        if weight(mid) > needed:
            lo=mid
        else:
            hi=mid
    assert abs(weight(hi)-needed)<2e-13
    return dict(weak_required_log_load=needed, common_two_pair_log_load=weight(pole),
                required_additional_log_load=residual, unification_scale_GeV=unification,
                inverse_gauge_spreads_at_unification=[.6*residual/(2*math.pi),residual/(2*math.pi),0.],
                weak_dirac_mass_repair_GeV=hi, weak_scalar_mass_repair_GeV=math.hypot(hi,soft))


def run():
    charged, neutral, legs, pole = parent(.25)
    obs = observable(charged,neutral)
    unitary = np.linalg.svd(obs["null"].reshape((3,1)), full_matrices=True)[0]
    poles = {}
    for i, label in enumerate(("max_neutral_area","min_neutral_area","intermediate")):
        oriented = np.roll(legs,-i)
        filt = (unitary*oriented) @ unitary.conj().T
        poles[label] = compare_row(charged, neutral, filt, obs)
    lower, upper = np.prod(legs)**2, max(legs[0],legs[2])**2
    np.testing.assert_allclose([poles["min_neutral_area"]["area"],poles["max_neutral_area"]["area"]],
                               [lower,upper],rtol=2e-12)
    rng = np.random.default_rng(1743)
    max_error = 0.
    for _ in range(48):
        unitary = np.linalg.qr(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))[0]
        filt = (unitary*legs) @ unitary.conj().T
        row = compare_row(charged, neutral, filt, obs)
        assert lower <= row["area"] <= upper
        max_error = max(max_error, abs(row["joint"]-1))
    gain_e, gain_n = .987, 1.011
    scaled = observable(gain_e*filt@charged, gain_n*filt.T@neutral@filt)
    np.testing.assert_allclose(scaled["joint"]/obs["joint"],gain_n**2/gain_e**4,rtol=2e-12)
    np.testing.assert_allclose(scaled["shape"]/obs["shape"],1.,rtol=3e-9)
    # A split operator on just the neutral bilinear gives a test of the joint condition.
    split = observable(charged, filt.T@neutral@filt)
    assert abs(split["joint"]/obs["joint"]-1)>1e-3
    return dict(pole_GeV=pole, leg_eigenvalues=legs.tolist(),
                complete_multiplets_on_previous_gauge_witness=gauge_matching(pole),
                original_massless_probabilities_e_mu_tau=obs["probabilities"].tolist(),
                original_gain_independent_pulled_back_Takagi_ratio=obs["shape"],
                parent_phase_blind_test=majorana_test(charged,neutral,obs["shape"]),
                transported_phase_blind_test=majorana_test(filt@charged,filt.T@neutral@filt,obs["shape"]),
                normalized_neutral_area_range=[lower,upper],
                normalized_null_charged_moment_range=[legs[0]**2,1.],
                charged_volume_ratio=np.prod(legs),
                max_error_joint_ratio_48_complex_frames=max_error,extrema=poles,
                neutral_only_example=dict(joint=split["joint"]/obs["joint"],shape=split["shape"]/obs["shape"]),
                independent_charged_and_mass_gains=dict(charged=gain_e,neutral=gain_n,joint=scaled["joint"]/obs["joint"]))


def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:
            compare(a[k],b[k])
    else:
        np.testing.assert_allclose(a,b,rtol=3e-8,atol=3e-12)


if __name__ == "__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--check",action="store_true")
    args=ap.parse_args()
    out=HERE/"independent_results.json"
    result=run()
    if args.check:
        compare(result,json.loads(out.read_bytes()))
        print("1743 independent charged/neutral cofactor identity, shape and attainable endpoints verified")
    else:
        with out.open("x",encoding="utf8") as f:
            json.dump(result,f,indent=2)
            f.write("\n")
