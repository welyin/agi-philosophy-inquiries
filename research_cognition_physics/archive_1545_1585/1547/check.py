"""Round 1547: mass spectra and complete/one-missing real octet pole sums."""
import argparse
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
NORMALIZATION = 3/(32*math.pi)


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def mass_box(rows, widths):
    values = np.asarray(rows, dtype=float)
    center, half = values.T
    low, high = center-widths*half, center+widths*half
    assert np.all(low > 0)
    for ordinary in (slice(0,3),slice(3,6)):
        assert np.all(high[ordinary][:-1]<low[ordinary][1:])
    gaps = np.maximum(0,np.maximum(low[:3]-high[3:], low[3:]-high[:3]))
    bb = max(0., low[2]-high[5])
    return {"minimal_spectral_gaps":gaps.tolist(),
            "normalized_all_down_sum_lower":NORMALIZATION*float(gaps@gaps),
            "normalized_bb_sum_lower":NORMALIZATION*bb**2,
            "bb_lower_if_one_real_pole_missing":NORMALIZATION*bb**2/2,
            "minimum_unrelated_bb_yukawa_shift_to_erase_bound":bb}


def basis_and_mixing_tests(central):
    # Yukawas are dimensionless; electroweak vevs divided by their common norm.
    d = np.diag(central[:3])
    e = np.diag(central[3:])
    n = 3
    omega = np.exp(2j*math.pi/n)
    fourier = omega**np.outer(np.arange(n), np.arange(n))/math.sqrt(n)
    e = fourier @ e @ np.diag([1j, 1., -1.]) @ fourier
    p = np.array([.36, .48j])
    p2 = np.vdot(p,p).real
    tensors = np.array([(math.sqrt(3)/2)*z.conjugate()*(d-e)/p2 for z in p])
    assert np.allclose(2/math.sqrt(3)*np.einsum("i,ijk->jk", p, tensors), d-e)
    # The orthogonal vev combination is physically free.
    arbitrary = np.array([[.001j, .002, -.003], [0., .001, .002j], [0., 0., .004]])
    tensors += np.array([p[1]*arbitrary, -p[0]*arbitrary])
    assert np.allclose(2/math.sqrt(3)*np.einsum("i,ijk->jk", p, tensors), d-e)
    real_fields = np.stack([tensors[0],1j*tensors[0],tensors[1],-1j*tensors[1]])
    # A fixed real rotation mixes all neutral CP assignments and both doublets.
    trial = np.array([[1,3,4,2],[3,-2,1,-3],[5,4,-2,1],[-1,2,3,4]],float)
    orthogonal, _ = np.linalg.qr(trial)
    mixed = np.einsum("pk,kij->pij",orthogonal,real_fields)
    for block in (mixed, mixed[:,2:3,2:3]):
        rates = np.sum(abs(block)**2,axis=(1,2))/(16*math.pi)
        assert max(rates) <= sum(rates)/2 + 1e-15
    rates = np.sum(abs(mixed)**2,axis=(1,2))/(16*math.pi)
    assert math.isclose(sum(rates),np.linalg.norm(tensors)**2/(8*math.pi))
    bb = np.abs(mixed[:,2,2])**2/(16*math.pi)
    assert sum(rates) >= NORMALIZATION/p2*np.sum((central[:3]-central[3:])**2)
    assert sum(bb) >= NORMALIZATION/p2*max(0.,central[2]-central[5])**2
    # Exact down-sector saturation at finite v15/v=.6; directions may be diagonal.
    aligned_y = (math.sqrt(3)/2)*(np.diag(central[:3]-central[3:]))/.6
    exact = np.linalg.norm(aligned_y)**2/(8*math.pi)
    assert math.isclose(exact,NORMALIZATION/.6**2*np.sum((central[:3]-central[3:])**2))
    # Two source-free real poles can exhaust a visible two-pole subspace.
    hidden_pair = np.array([aligned_y,1j*aligned_y,0*aligned_y,0*aligned_y])
    assert np.linalg.norm(hidden_pair[2:])==0 and np.linalg.norm(hidden_pair[:2])>0
    t15 = np.diag([1.,1.,1.,-3.])/(2*math.sqrt(6))
    t8 = np.diag([1.,1.,-2.])/(2*math.sqrt(3))
    assert np.isclose(np.trace(t15 @ t15),.5) and np.isclose(np.trace(t8 @ t8),.5)
    assert np.isclose(math.sqrt(2)*(t15[0,0]-t15[3,3]),2/math.sqrt(3))
    return {"random_orientation_all":sum(rates), "random_orientation_bb":sum(bb),
            "largest_real_pole_share_bb":max(bb)/sum(bb),
            "aligned_saturating_sum_at_vev_fraction_0p6":exact,
            "visible_fraction_two_silent_poles":0.}


def study():
    cfg = read(HERE/"inputs.json")
    low = cfg["light_scales"]
    high = read(HERE/cfg["high_scale"])["scales"]["1e16"]
    data = dict(low)
    data["1e16"] = high
    scales = {key:{"central":mass_box(rows,0),"two_width":mass_box(rows,cfg["widths"])}
              for key, rows in data.items()}
    visible = cfg["visible_example"]
    beta3 = (1-4*(visible["max_kinematic_down_mass_GeV"]/visible["min_scalar_GeV"])**2)**1.5
    b = scales["1000"]["two_width"]["normalized_bb_sum_lower"]
    # A left-handed lepton rotation can cancel bb alone when tau exceeds bottom.
    c = high[2][0]/high[5][0]
    s = math.sqrt(1-c*c)
    u = np.array([[c,0.,s],[0.,1.,0.],[-s,0.,c]])
    residual = np.diag(np.array(high)[:3,0])-u@np.diag(np.array(high)[3:,0])
    assert abs(residual[2,2])<1e-15 and np.linalg.norm(residual)>0
    # Traces include the spectator color and right-isospin components.
    base_ps = np.array([22/3, 3., -13/3])
    extra_dynkin = np.array([4*2*2, .5*15*2, .5*15*2])
    complex_beta = base_ps-extra_dynkin/3
    real_beta = base_ps-extra_dynkin/6
    assert complex_beta[1] == -2 and real_beta[1] == .5
    running = {"full_complex_b4_L_R":complex_beta.tolist(), "real_b4_L_R":real_beta.tolist(),
               "critical_log_fraction_full_complex_for_positive_bL":base_ps[1]/(extra_dynkin[1]/3)}
    return {"running":running,"scales":scales,"finite_phase_space_floor":beta3,
            "one_TeV_bb_all_open_width_ratio":b*beta3,
            "one_TeV_bb_one_missing_width_ratio":b*beta3/2,
            "one_TeV_four_real_poles_at_least_one_width_GeV_min_mass_1TeV":1000*b*beta3/4,
            "high_scale_cos_left_rotation_for_exact_bb_zero":c,
            "source_checks":basis_and_mixing_tests(np.asarray(low["1000"])[:,0])}


def stable(v):
    if isinstance(v,dict):
        return {k:stable(x) for k,x in v.items()}
    if isinstance(v,list):
        return [stable(x) for x in v]
    return float(format(float(v),".13g"))


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument("--write-new",action="store_true")
    mode.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=stable(study())
    target=HERE/"results.json"
    if args.write_new:
        with target.open("x",encoding="utf-8") as output:
            json.dump(result,output,indent=2)
            output.write("\n")
    elif args.check:
        assert read(target)==result, "Stored output differs"
        print("1547: color normalization, mass boxes, neutral pole completeness and b readout verified.")
    else:
        print(json.dumps(result,indent=2))
