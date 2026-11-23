"""Independent two-pole flavor and CP test for the Dirac charged-singlet loop."""
import argparse
from itertools import product
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def outer(readout):
    return readout[0]-2*readout[1], readout[0]+2*readout[-1]


def rotation(i, j, sin2, phase=0):
    s, c = math.sqrt(sin2), math.sqrt(1-sin2)
    u = np.eye(3, dtype=complex)
    u[i,i] = u[j,j] = c
    u[i,j] = s * np.exp(-1j*phase)
    u[j,i] = -u[i,j].conj()
    return u


def weights(g):
    return np.sum(abs(g)**2, axis=1)


def make_joint_example(cfg, osc):
    a = osc["DayaBay_amplitude"][0]
    s13 = (1-math.sqrt(1-a))/2
    par = osc["NO"]
    s12 = par["solar_sine_squared"][0]
    s23 = par["atmospheric_sine_squared"][0]
    delta = math.pi/3
    U = rotation(1,2,s23) @ rotation(0,2,s13,delta) @ rotation(0,1,s12)
    n = U[:,0].conj()
    # Cross product n with the three coordinate directions, giving an antisymmetric f.
    f = np.cross(n, np.eye(3)).T
    f *= 0.001
    light = np.array(cfg["charged_masses_GeV"])
    masses = np.array(cfg["benchmark"]["charged_poles_GeV"])
    sine = cfg["benchmark"]["sin_theta"]
    cosine = math.sqrt(1-sine**2)
    ratio = (light[:,None]/masses[None,:])**2
    d = light*np.diff(-np.log(ratio)/(1-ratio), axis=1)[:,0]
    a_matrix = (2*sine*cosine/(16*np.pi**2))*f@np.diag(d)
    nu = 1e-9*np.sqrt([0, par["solar_gap"][0],
                       par["solar_gap"][0]+par["atmospheric_gap"][0]])
    target = U @ np.diag(nu)
    # Fix R_e=0. The remaining two rows are then determined, with no PMNS input
    # from the proposed scalar readout. Equal total widths are just a witness.
    g = np.zeros((3,3), dtype=complex)
    g[1:], _, rank, _ = np.linalg.lstsq(a_matrix[:,1:], target, rcond=None)
    assert rank == 2
    scaling = (weights(g).sum()/(4*weights(f).sum()))**.25
    f, g = scaling*f, g/scaling
    assert np.allclose(scaling*a_matrix@g, target, atol=1e-23)
    left, right = 4*weights(f), weights(g)
    assert np.isclose(left.sum(), right.sum(), rtol=1e-12, atol=0)
    truth = left/left.sum()
    assert np.allclose(truth, (1-abs(n)**2)/2)
    rows = []
    all_rates = [cosine**2*left+sine**2*right,sine**2*left+cosine**2*right]
    for m,rate in zip(masses, all_rates):
        phase = (1-(light/m)**2)**2
        detected = rate*phase
        corrected = detected/phase
        rows.append(corrected/corrected.sum())
    # Independently convert the sharp one-flavor maximum to a single Hermitian
    # eigenproblem on the actual affine space (include its null components).
    A = scaling*a_matrix
    _,_,Vh = np.linalg.svd(A)
    p = Vh.conj().T[:,-1]
    g0 = g-np.outer(p,p.conj()@g)
    g0 /= np.linalg.norm(g0)
    frame = np.column_stack([g0.ravel()] + [np.outer(p,v).ravel() for v in np.eye(3)])
    assert np.allclose(frame.conj().T@frame,np.eye(4))
    cut = frame[6:9]
    sharp_tau = np.linalg.eigvalsh(cut.conj().T@cut)[-1]
    pred = abs(p[2])**2+weights(g0)[2]
    assert abs(sharp_tau-pred) < 1e-12 and sharp_tau<.006701
    result = dict(phase_degrees=60, largest_abs_f=float(np.max(abs(f))),
                  largest_abs_g=float(np.max(abs(g))), two_pole_corrected_fractions=np.array(rows).tolist(),
                  left=truth.tolist(), right=(right/right.sum()).tolist(),
                  single_right_tau_supremum=float(sharp_tau))
    # Parent benchmark's direct dipole formula, used only to audit this new witness.
    f12 = abs((f@f.conj().T)[1,0])
    g12 = abs((g@g.conj().T)[1,0])
    gf = cfg["benchmark"]["fermi_GeV_minus2"]
    alpha = 1/cfg["benchmark"]["inverse_alpha"]
    br = alpha/(768*np.pi*gf**2)*(16*f12**2*(cosine**2/masses[0]**2+sine**2/masses[1]**2)**2
                                   +g12**2*(sine**2/masses[0]**2+cosine**2/masses[1]**2)**2)
    result["illustrative_dipole_branch"] = float(br)
    return result


def invert(rows, e_left, uncertainty):
    # Each variable enters linearly or as a linear fraction with fixed
    # denominator sign, so each extreme over the box occurs at a corner.
    a, b = np.asarray(rows)[:,:2]
    assert abs(a[0]-b[0])>2*uncertainty
    values = []
    for x1, y1, x2, y2 in product(*[(u-uncertainty,u+uncertainty)
                                       for u in (a[0],a[1],b[0],b[1])]):
        for x in e_left:
            values.append(y1+(x-x1)*(y2-y1)/(x2-x1))
    return min(values), max(values)


def analyze(cfg, osc):
    result = make_joint_example(cfg,osc)
    rows = result["two_pole_corrected_fractions"]
    no, io = osc["NO"], osc["IO"]
    s13_range = [(1-math.sqrt(1-v))/2 for v in outer(osc["DayaBay_amplitude"])]
    s13 = (1-math.sqrt(1-osc["DayaBay_amplitude"][0]))/2
    s12,s23 = no["solar_sine_squared"][0],no["atmospheric_sine_squared"][0]
    error = .001
    e0 = (1-(1-s12)*(1-s13))/2
    mu = invert(rows,[e0],error)
    constant = s12*(1-s23)+(1-s12)*s23*s13
    coefficient = 2*math.sqrt(s12*(1-s12)*s23*(1-s23)*s13)
    fixed_cos = [(1-2*v-constant)/coefficient for v in mu[::-1]]
    assert 0<fixed_cos[0]<.5<fixed_cos[1]<1
    # NO with existing oscillation errors: outer intervals, including all phases.
    s12lo,s12hi = outer(no["solar_sine_squared"])
    s23lo,s23hi = outer(no["atmospheric_sine_squared"])
    low, hi = s13_range
    e_range = ((1-(1-s12lo)*(1-low))/2,(1-(1-s12hi)*(1-hi))/2)
    mu_box = invert(rows,e_range,error)
    const_min = s12lo*(1-s23hi)+(1-s12hi)*s23lo*low
    const_max = s12hi*(1-s23lo)+(1-s12lo)*s23hi*hi
    cden_min = 2*math.sqrt(s12lo*(1-s12hi)*s23lo*(1-s23hi)*low)
    cden_max = 2*math.sqrt(s12hi*(1-s12lo)*s23hi*(1-s23lo)*hi)
    cos_box = [n/d for n,d in product([1-2*mu_box[1]-const_max,
                                       1-2*mu_box[0]-const_min],[cden_min,cden_max])]
    imulo,imuhi = outer(io["atmospheric_sine_squared"])
    io_e = [(1-hi)/2,(1-low)/2]
    io_mu = [(1-imuhi*(1-low))/2,(1-imulo*(1-hi))/2]
    intersection = invert(rows,io_e,error)
    assert intersection[1] < io_mu[0]
    result.update(face_absolute_error=error, fixed_angle_cos_delta=fixed_cos,
                  NO_cos_delta_outer=[max(-1,min(cos_box)),min(1,max(cos_box))],
                  IO_required_mu=io_mu, IO_line_extrapolated_mu=intersection,
                  common_phase_conjugation_same=True)
    return result


def close(x,y):
    if isinstance(x,dict):
        return x.keys()==y.keys() and all(close(x[k],y[k]) for k in x)
    return np.allclose(x,y,rtol=3e-9,atol=1e-25)


if __name__=="__main__":
    pars = argparse.ArgumentParser()
    mode = pars.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write",action="store_true")
    mode.add_argument("--check",action="store_true")
    opts = pars.parse_args()
    cfg = read(HERE/"inputs.json")
    result = analyze(cfg,read(HERE/cfg["oscillation_input"]))
    path = HERE/"independent_results.json"
    if opts.write:
        with path.open("x",encoding="utf-8") as handle:
            json.dump(result,handle,indent=2,allow_nan=False)
            handle.write("\n")
    elif not close(read(path),result):
        raise SystemExit("Independent mismatch")
    else:
        print("1538: joint two-pole flavor inversion, ordering and direct row checks passed")
