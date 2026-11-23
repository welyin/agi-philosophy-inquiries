"""Test the entire up-aligned diagonal-congruence orbit using a physical row."""
import argparse
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def load(path):
    return json.loads(path.read_text(encoding="utf-8"))


def left_axes(y):
    u, values, _ = np.linalg.svd(y)
    return values[::-1], u[:,::-1]


def tensors(point):
    t = point["tensors"]
    h = np.diag(t["H_diag"])
    f = np.zeros((3,3),dtype=complex)
    i,j = np.triu_indices(3)
    f[i,j] = np.array(t["F_real_upper"]) + 1j*np.array(t["F_imag_upper"])
    f[j,i] = f[i,j]
    # Symmetric two-component masses: conjugate to use the left-column convention.
    return np.conj(t["r"]*(h+complex(t["s_real"],t["s_imag"])*f)), np.conj(h+f)


def charged_bounds(inp, widths, fixed_first_mass=None):
    def box(x):
        return np.array([x[0]-widths*x[1],x[0]+widths*x[1]])
    down = np.array([box(x) for x in inp["down"]])
    if fixed_first_mass is not None:
        down[0] = fixed_first_mass
    t12,t13 = [box(x) for x in inp["angles_theta12_theta13_theta23_delta"][:2]]
    p1 = np.cos(t12[::-1])**2*np.cos(t13[::-1])**2
    p2 = np.sin(t12)**2*np.cos(t13[::-1])**2
    p3 = np.sin(t13)**2
    p = np.vstack((p1,p2,p3))
    # Every positive term in the pair expansion is monotone with the mass ratio.
    def endpoint(side):
        w = 1.
        for i in range(3):
            for j in range(i+1,3):
                ratio = down[j,side]/down[i,1-side]
                assert ratio > 1
                w += p[i,side]*p[j,side]*(ratio-1/ratio)**2
        return float(w)
    return [endpoint(0),endpoint(1)]


def invariants(gram):
    diag = np.diag(gram).real
    corr = gram / np.sqrt(diag[:,None]*diag[None,:])
    return (np.abs(corr[np.triu_indices(3,1)]),
            corr[0,1]*corr[1,2]*corr[2,0],
            diag*np.diag(np.linalg.inv(gram)).real)


def evaluate():
    cfg = load(HERE/"inputs.json")
    src = load(HERE/cfg["shared_tensors"])["candidates"]
    high = load(HERE/cfg["charged_target"])["high_scale"]
    window = charged_bounds(high,cfg["widths"])
    result = {"target_Wu_outer_box":window,"points":{}}
    for key in cfg["keys"]:
        up,down = tensors(src[key])
        _,u = left_axes(up)
        masses,v = left_axes(down)
        ckm = u.conj().T@v
        basis = u.conj().T@down
        gram = basis@basis.conj().T
        rho,phase,ws = invariants(gram)
        row = abs(ckm[0])**2
        direct = (row@masses**2)*(row@(1/masses**2))
        np.testing.assert_allclose(direct,ws[0],rtol=2e-10)
        for gains in [[0.4+0.5j,-1.2j,2.4],[2.,0.7,1.1j]]:
            s = np.diag(gains)
            test = invariants(s@gram@s.conj().T)
            for a,b in zip((rho,phase,ws),test):
                np.testing.assert_allclose(a,b,rtol=2e-10)
        assert direct > window[1]
        cond_min = np.sqrt(direct/window[1])
        low = masses[0]/2
        high_bound = high["down"][1][0]/2
        def maximum_W(first):
            return charged_bounds(high,cfg["widths"],first)[1]
        assert maximum_W(low) > direct > maximum_W(high_bound)
        for _ in range(72):
            middle = (low+high_bound)/2
            if maximum_W(middle) >= direct:
                low = middle
            else:
                high_bound = middle
        a,b,c = (masses/masses[-1])**2
        power = high["down"][0][0]/masses[0]
        # In the existing right down-mass basis, one extra spectral response suffices.
        _,_,right_h = np.linalg.svd(down)
        right_axes = right_h.conj().T[:,::-1]
        projector = np.outer(right_axes[:,0],right_axes[:,0].conj())
        corrected = down@(np.eye(3)+(power-1)*projector)
        new_mass,new_axes = left_axes(corrected)
        np.testing.assert_allclose(new_mass, [high["down"][0][0],*masses[1:]],rtol=2e-11)
        np.testing.assert_allclose(abs(u.conj().T@new_axes),abs(ckm),rtol=3e-10)
        assert 0<a<b<c and power>1
        result["points"][key] = {"down_masses":masses.tolist(),"up_row_probabilities":row.tolist(),
            "Wu":float(direct),"additional_right_condition_number_lower":float(cond_min),
            "first_mass_yukawa_upper_with_other_inputs_in_box":float(high_bound),
            "right_spectral_gain_witness":float(power),
            "necessary_slope_of_normalised_right_gain":float((power-1)/(b-a)),
            "complex_triangle":[float(phase.real),float(phase.imag)]}
    return result


def compare(x,y):
    if isinstance(x,dict):
        assert x.keys() == y.keys()
        for k in x:
            compare(x[k],y[k])
    else:
        np.testing.assert_allclose(x,y,rtol=5e-10,atol=1e-15)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    check = parser.parse_args().check
    result = evaluate()
    dest = HERE/"results.json"
    if check:
        compare(result,load(dest))
        print("1561 full diagonal-congruence obstruction and independent down direction verified")
    else:
        dest.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
