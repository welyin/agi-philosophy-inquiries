"""Coherent propagation with reversed actual layer order, using NumPy only."""
import argparse
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def rotation(i, j, sine):
    matrix = np.eye(3)
    matrix[i, i] = matrix[j, j] = np.sqrt(1-sine*sine)
    matrix[i, j], matrix[j, i] = sine, -sine
    return matrix


def transport(mass, energies, lengths, densities, phase_unit, matter):
    result = np.broadcast_to(np.eye(3), (len(energies),3,3)).copy()
    for length, density in zip(lengths, densities):
        h = mass / energies[:,None,None]
        h[:,0,0] += matter*density
        roots, frames = np.linalg.eigh(h)
        step = (frames*np.exp(-1j*phase_unit*length*roots)[:,None,:]) @ frames.conj().transpose(0,2,1)
        result = step @ result
    return result


def run(d):
    pos = np.asarray(d["position_km"])
    lengths = np.diff(np.r_[pos[0], (pos[1:]+pos[:-1])/2, pos[-1]])
    energies = np.linspace(*d["energy_GeV"], d["intervals"]+1)
    center = len(energies)//2
    center_energy = energies[center:center+1]
    phase = 1e-13/(2*d["hbarc_GeV_cm"])
    matter = 2e18*np.sqrt(2)*d["GF"]*d["NA"]*d["hbarc_GeV_cm"]**3*d["electron_fraction"]
    s13 = np.sqrt((1-np.sqrt(1-d["theta13_amplitude"]))/2)
    alpha, alo, ahi = d["alpha_rad"]
    phases = {"minus_a":-alpha, "plus_a":alpha, "pi_minus_a":np.pi-alpha, "pi_plus_a":np.pi+alpha}
    result = {}
    for profile, densities in d["density_g_cm3"].items():
        result[profile] = {}
        for order, (s12sq,s23sq,solar,atm) in d["mass_and_angles"].items():
            r = rotation(0,2,s13) @ rotation(0,1,np.sqrt(s12sq))
            mass = (r*[0,solar,atm]) @ r.T
            sf = transport(mass, energies, lengths, densities, phase, matter)
            sr = transport(mass, center_energy, lengths[::-1], densities[::-1], phase, matter)
            assert np.allclose(sf[center].T, sr[0], atol=1e-13)
            x, y = sf[:,1,0], sf[:,2,0]
            product = np.conj(x)*y
            loss = abs(x)**2+abs(y)**2
            assert np.allclose(loss, 1-abs(sf[:,0,0])**2)
            # Each of Re/Im(a* b) and the survival probability has this Lipschitz bound.
            lip = phase*sum(lengths)*np.ptp([0,solar,atm]) / energies[0]**2
            margin = lip*(energies[-1]-energies[0])/(2*d["intervals"])
            def signed_abs_bound(values):
                assert np.all(values*values[center]>0)
                bound = np.min(np.abs(values))-margin
                assert bound>0
                return bound
            real_lower = signed_abs_bound(product.real)*np.cos(ahi)
            imag_lower = signed_abs_bound(product.imag)*np.sin(alo)
            loss_upper = np.max(loss)+margin
            s23, c23 = np.sqrt(s23sq), np.sqrt(1-s23sq)
            factor = 2*s23*c23
            # Sign of the centered even signal without knowing the propagation-basis split.
            cos_floor = factor*real_lower - abs(s23sq-.5)*loss_upper
            assert cos_floor>0
            span = {}
            for label, delta in phases.items():
                k = rotation(1,2,s23) @ np.diag([1,1,np.exp(1j*delta)])
                central_f = k @ sf[center] @ k.conj().T
                central_r = k @ sr[0] @ k.conj().T
                f,rev = abs(central_f[1,0])**2, abs(central_r[0,1])**2
                assert np.isclose((f-rev)/2, -factor*product[center].imag*np.sin(delta))
                assert np.isclose((f+rev)/2-loss[center]/2,
                    (.5-s23sq)*(abs(x[center])**2-abs(y[center])**2)
                    +factor*product[center].real*np.cos(delta))
                span[label] = dict(forward_e_mu=f, reversed_mu_e=rev, forward_e_survive=1-loss[center],
                    half_difference=(f-rev)/2, centered_half_sum=(f+rev-loss[center])/2)
            zero_f = rotation(1,2,s23) @ sf[center] @ rotation(1,2,s23).T
            result[profile][order] = dict(mean_density=float(np.average(densities,weights=lengths)),
                central=span, grid_margin=margin, min_abs_odd=factor*imag_lower,
                min_abs_centered_even=cos_floor, loss_upper=loss_upper,
                sign_Re=int(np.sign(product[center].real)), sign_Im=int(np.sign(product[center].imag)),
                sufficient_max_deviation_s23sq_from_half=real_lower/np.hypot(loss_upper,2*real_lower),
                wrong_order_same_direction_half_difference_at_zero_delta=(abs(zero_f[1,0])**2-abs(zero_f[0,1])**2)/2)
    return result


def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for key in a:
            compare(a[key],b[key])
    else:
        assert np.isclose(a,b,rtol=1e-9,atol=1e-12)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    mode=p.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check",action="store_true")
    mode.add_argument("--write-new",action="store_true")
    args=p.parse_args()
    result=run(json.loads((HERE/"inputs.json").read_text(encoding="utf-8-sig")))
    target=HERE/"results.json"
    if args.write_new:
        with target.open("x",encoding="utf-8") as stream:
            json.dump(result,stream,indent=2)
            stream.write("\n")
    else:
        compare(result,json.loads(target.read_text(encoding="utf-8")))
        print("1423 reciprocal physical profiles verified")
