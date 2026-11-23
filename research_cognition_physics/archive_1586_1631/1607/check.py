"""Right-flavor source recognition and a conditional inclusive Z-pole test."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def load(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def mixing(sines, delta):
    a,b,c = sines
    A,B,C = np.sqrt(1-sines**2)
    z = np.exp(1j*delta)
    return np.array([[A*B,a*B,b/z],
                     [-a*C-A*c*b*z,A*C-a*c*b*z,c*B],
                     [a*c-A*C*b*z,-A*c-a*C*b*z,C*B]])


def ranges(table, widths):
    middle, error = np.array(table).T
    return middle - widths*error, middle + widths*error


def dominance_floors(slo, shi, dlo, dhi):
    a,b,c = slo
    A,B,C = np.sqrt(1-shi**2)
    x,y,z = shi
    # Bounds on |V_id| from below, |V_is| and |V_ib| from above.
    first = np.array([A*B, a*C-z*y, a*c-y])
    other = np.array([x, 1, z+x*y])
    last = np.array([y, z, 1])
    assert min(first) > 0
    return 1 / (1+(other*dhi[0]/(first*dlo[1]))**2+
                (last*dhi[0]/(first*dlo[2]))**2)


def singlet_mass_witness(v, masses, enlarged):
    """Two vectorlike down singlets, a single Higgs insertion and hard heavy rows."""
    # Columns give the three tagged ordinary gauge states in the full mass basis.
    _, _, vh = np.linalg.svd(enlarged.conj().T, full_matrices=True)
    mass_from_gauge = np.column_stack((enlarged, vh.conj().T[:,3:]))
    left = np.eye(5, dtype=complex)
    left[:3,:3] = v
    physical = np.diag([*masses, 3., 7.])
    matrix = left @ physical @ mass_from_gauge
    higgs_rows = np.diag([1,1,1,0,0]) @ matrix
    back = mass_from_gauge.conj().T
    assert np.allclose(mass_from_gauge @ back, np.eye(5))
    assert np.allclose(left.conj().T @ matrix @ back, physical)
    assert np.allclose(left.conj().T @ higgs_rows @ back, np.diag([*masses,0.,0.]))
    assert np.allclose(left.conj().T @ np.diag([1,1,1,0,0]) @ left, np.diag([1,1,1,0,0]))
    assert np.allclose(left[:3,3:], 0) and np.allclose(left[:3,:3],v)


def tag_test(sines, phase, masses, floors):
    v = mixing(sines, phase)
    inv = np.diag(1/masses) @ v.conj().T
    scales = 1/np.linalg.norm(inv, axis=0)
    tags = inv * scales
    gram = tags.conj().T @ tags
    assert np.allclose(v@np.diag(masses)@tags, np.diag(scales))
    eig, vec = np.linalg.eigh(gram)
    f = 1/eig[-1]
    # A single orthonormal triple in 3 light + 2 extra dimensions.
    extra = np.diag(np.sqrt(1-f*eig[:2])) @ vec[:,:2].conj().T
    enlarged = np.vstack((np.sqrt(f)*tags, extra))
    assert np.allclose(enlarged.conj().T @ enlarged, np.eye(3))
    singlet_mass_witness(v, masses, enlarged)
    assert np.all(abs(tags[0])**2>=floors) and min(floors)>0.5
    j = (v[0,0]*v[1,1]*v[0,1].conjugate()*v[1,0].conjugate()).imag
    cyclic = (gram[0,1]*gram[1,2]*gram[2,0]).imag
    inv2 = 1/masses**2
    product = np.prod([inv2[i]-inv2[l] for i,l in itertools.combinations(range(3),2)])
    cp_formula = j*np.prod(scales**2)*product
    assert np.isclose(abs(cp_formula), abs(cyclic), rtol=3e-4, atol=1e-16)
    pairs = {}
    for i,l in itertools.combinations(range(3),2):
        q = abs(gram[i,l])
        bound = math.sqrt(floors[i]*floors[l])-math.sqrt((1-floors[i])*(1-floors[l]))
        assert 0 < bound < q < 1
        light = np.column_stack((tags[:,i],tags[:,l])) / np.sqrt(1+q)
        hidden = np.sqrt(q/(1+q))*np.array([1, -gram[i,l]/q])
        pair = np.vstack((light, hidden))
        assert np.allclose(pair.conj().T@pair, np.eye(2))
        pairs["uct"[i]+"uct"[l]] = {"overlap":q, "overlap_floor":bound,
            "unambiguous_equal_success_max":1-q, "hidden_weight_floor":2*bound/(1+bound)}
    return {"pairs":pairs,"d_probability_floors":floors.tolist(),"central_d_probabilities":(abs(tags[0])**2).tolist(),
            "common_light_fraction_max":f,"imag_cyclic":cyclic,
            "minimum_number_hidden_for_three":2, "aligned_singlet_mass_and_Higgs_witness":True}


def width_floor(overlap, s_window):
    """Analytic conservative lower enclosure over sin^2 theta in this interval."""
    lo,hi = s_window
    eta = 2*overlap/(1+overlap)
    # Five standard left currents. Three right d currents plus one irreducible right u.
    base = 2*(.5-2*hi/3)**2+3*(.5-hi/3)**2+3*(lo/3)**2+(2*lo/3)**2
    down_gain = eta**2/4-hi*eta/3
    assert eta/2>hi/3
    lepton_max = max((.5-s)**2+s*s for s in [lo,hi])
    return 3*(base+down_gain)/lepton_max


def calculate():
    inputs = load(HERE/"inputs.json")
    weak = load(HERE/inputs["weak_inputs"])
    high = load(HERE/inputs["high_inputs"])["high_scale"]
    scales = {}
    for name,table,angle_table,is_angle in [
        ("MZ",weak["yukawa"]["down"],weak["quark_sines_and_phase"],False),
        ("1e16",high["down"],high["angles_theta12_theta13_theta23_delta"],True)]:
        m = np.asarray(table)[:,0]
        dlo,dhi = ranges(table, inputs["widths"])
        angle = np.asarray(angle_table)[:,0]
        slo,shi = ranges(angle_table[:3],inputs["widths"])
        middle = angle[:3]
        if is_angle:
            slo,shi,middle = np.sin(slo),np.sin(shi),np.sin(middle)
        floors = dominance_floors(slo,shi,dlo,dhi)
        scales[name] = tag_test(middle,angle[3],m,floors)
    q = min(p["overlap_floor"] for p in scales["MZ"]["pairs"].values())
    lower = width_floor(q,inputs["leptonic_sin2_window"])
    obs,err = inputs["pole_ratio"]
    allowed = obs+inputs["widths"]*err
    factor = inputs["lower_pole_transfer"]
    assert lower*factor>allowed
    # Invert the analytic bound on the increasing physical branch.
    small, large = 0.25,1.
    for _ in range(70):
        mid = (large+small)/2
        if width_floor(mid,inputs["leptonic_sin2_window"])*factor>allowed:
            large = mid
        else:
            small = mid
    return {"scales":scales, "weak_all_pairs_floor":q, "tree_ratio_floor":lower,
            "working_pole_ratio_floor":factor*lower, "measured_two_width_upper":allowed,
            "largest_pole_transfer_to_escape_this_bound":allowed/lower,
            "overlap_above_which_bound_excludes":large}


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    args = p.parse_args()
    file = HERE/"results.json"
    result = calculate()
    if args.check:
        assert result == load(file)
    else:
        with file.open("x", encoding="utf-8") as stream:
            stream.write(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print("1607: simultaneous sharp tags, projection cost and Z ratio verified")
