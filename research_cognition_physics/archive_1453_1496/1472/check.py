"""Common-scale, one-loop test of the charged matter used in the flavor witness."""
import argparse
from fractions import Fraction
import itertools
import json
import math
from pathlib import Path
import runpy

import numpy as np

BASE = Path(__file__).resolve().parent


def read(name):
    return json.loads((BASE / name).read_text(encoding="utf-8-sig"))


def screening(multiplets):
    hypercharge = Fraction(0)
    weak = Fraction(0)
    for twice_spin, charge in multiplets:
        spin = Fraction(twice_spin, 2)
        dim = twice_spin + 1
        hypercharge += dim * Fraction(charge) ** 2
        weak += dim * spin * (spin + 1) / 3
        assert all((Fraction(charge) + i - spin).denominator == 1 for i in range(dim))
    return hypercharge, weak


def crossing(g, db, z_mass):
    inverse = 4 * math.pi / np.square(g)
    b = np.array([41/10, -19/6, -7.])
    db = np.array(db)
    # L is log(U/MZ), t is log(U/M). Both independent differences vanish.
    L, t = np.linalg.solve(np.column_stack([b[:2]-b[2], db[:2]-db[2]]),
                           2 * math.pi * (inverse[:2]-inverse[2]))
    final = inverse - (b*L+db*t)/(2*math.pi)
    assert np.ptp(final) < 1.e-10
    return dict(L=L, t=t, U_GeV=z_mass*math.exp(L), M_GeV=z_mass*math.exp(L-t),
                common_inverse=float(final[0]), physical=bool(0 < t < L and final[0] > 0))


def match_below_threshold(low, g, z_mass, target, charges, epsilon, steps):
    # At M^- the low-energy theory of this witness has the SM fields only.
    sines_phase = [v[0] for v in low["quark_sines_and_phase"]]
    angles = list(np.arcsin(sines_phase[:3])) + [sines_phase[3]]
    mixing = runpy.run_path(str(BASE/"../1450/check.py"))["ckm_from_angles"](angles)
    yu,yd,ye = [np.diag([y for y,dy in low["yukawa"][name]])
                for name in ("up","down","charged_leptons")]
    state = np.array([yu,mixing@yd,ye],dtype=complex)
    gauge_c = np.array([[17/20,9/4,8],[1/4,9/4,8],[9/4,9/4,0]])
    beta = np.array([41/10,-19/6,-7])

    def rhs(log_scale, state):
        g_squared = 1/(1/np.square(g)-beta*log_scale/(8*math.pi**2))
        hu,hd,he = state@state.conj().transpose(0,2,1)
        trace = np.trace(3*hu+3*hd+he).real
        gauge_part = (trace-gauge_c@g_squared)[:,None,None]*state
        flavor_part = 1.5*np.array([hu-hd,hd-hu,he])@state
        return (flavor_part+gauge_part)/(16*math.pi**2)

    h = math.log(target/z_mass)/steps
    for i in range(steps):
        t = i*h
        a = rhs(t,state)
        b = rhs(t+h/2,state+h*a/2)
        c = rhs(t+h/2,state+h*b/2)
        d = rhs(t+h,state+h*c)
        state += h*(a+2*b+2*c+d)/6
    # Right-handed flavor bases are free because all left charges vanish.
    norms=[]
    for mat, exponents in zip(state,charges):
        u, values, vh = np.linalg.svd(mat)
        coefficient = u[:,::-1]@np.diag(values[::-1]/epsilon**np.array(exponents))
        reconstructed = coefficient@np.diag(epsilon**np.array(exponents))@vh[::-1]
        assert np.max(abs(reconstructed-mat)) < 1.e-13
        norms.extend(np.linalg.svd(coefficient,compute_uv=False))
    return [float(min(norms)),float(max(norms))]


def calculate():
    spec = read("inputs.json")
    flavor = read(spec["flavor"])["high_scale_witnesses"][0]
    qcd = read(spec["qcd"])["qcd"]
    ew = read(spec["gauge"])
    D = flavor["D"]
    original = Fraction(13, 5)
    assert D == 30 and abs(float(original)-flavor["E_over_N"]) < 1.e-12
    c = Fraction(str(qcd["photon_subtraction"]))
    dc = Fraction(str(qcd["photon_subtraction_width"]))
    integers = list(range(math.ceil((c-dc-original)*D/2), math.floor((c+dc-original)*D/2)+1))
    assert integers == [-10]
    k0 = qcd["m_times_f_GeV2"]
    dk = sum(qcd["m_times_f_errors_GeV2"])
    mass = spec["axion_example_micro_eV"] * 1.e-15
    U0 = spec["high_scale_GeV"]
    cap = spec["yukawa_cap"]
    max_mass = cap*D*(k0+dk)/(math.sqrt(2)*mass)
    cost = abs(integers[0])
    g = [ew[key] for key in ("g1_GUT", "g2", "g3")]
    rows = {}
    for name, spectrum in spec["spectra"].items():
        Y, W = screening(spectrum)
        assert Y+W == cost
        db = [float(Fraction(4,5)*Y),float(Fraction(4,3)*W),0.]
        row = crossing([v[0] for v in g], db, ew["mu_ref_GeV"])
        envelope = [crossing(values, db, ew["mu_ref_GeV"]) for values in
                    itertools.product(*[(v-2*dv,v+2*dv) for v,dv in g])]
        assert all(r["physical"] == row["physical"] for r in envelope)
        row.update(beta=db, two_width_ranges={key:[min(r[key] for r in envelope),
                    max(r[key] for r in envelope)] for key in ("L","t","M_GeV","U_GeV","common_inverse")})
        if row["physical"]:
            row["required_y_central"] = math.sqrt(2)*row["M_GeV"]*mass/(D*k0)
            low,high=row["two_width_ranges"]["M_GeV"]
            row["required_y_range"]=[math.sqrt(2)*low*mass/(D*(k0+dk)),
                                     math.sqrt(2)*high*mass/(D*(k0-dk))]
            row["axion_mass_ceiling_micro_eV_central_M"] = cap*D*k0/(math.sqrt(2)*row["M_GeV"])*1.e15
        rows[name]=row
    assert not rows["integer_singlets"]["physical"] and rows["weak_quintet"]["physical"]
    threshold = rows["weak_quintet"]["M_GeV"]
    unify = rows["weak_quintet"]["U_GeV"]
    epsilon = read(spec["old_flavor_inputs"])["witness_epsilon"]
    g0 = [v[0] for v in g]
    extrema = [match_below_threshold(read(spec["low_flavor"]),g0,ew["mu_ref_GeV"],threshold,
               flavor["right_charges"],epsilon,steps) for steps in (100,200)]
    assert np.max(np.abs(np.array(extrema[0])-extrema[1])) < 1.e-8
    assert 1/3 < extrema[1][0] < extrema[1][1] < 3
    F = D*k0/(spec["joint_example_micro_eV"]*1.e-15)
    cutoff = F/(math.sqrt(2)*epsilon)
    assert threshold < F < unify < cutoff
    joint = dict(mass_micro_eV=spec["joint_example_micro_eV"],F_GeV=F,
                 effective_flavor_scale_GeV=cutoff,y=math.sqrt(2)*threshold/F,
                 core_extrema_at_M_minus=extrema[1],
                 mass_ceiling_for_flavor_scale_above_U_micro_eV=D*k0/(math.sqrt(2)*epsilon*unify)*1.e15)
    return dict(allowed_delta_E=integers, corrected_E_over_N=str(original+2*Fraction(integers[0],D)),
                minimum_sum_Q_squared=cost, min_delta_bY_plus_b2=4*cost/3,
                latest_new_threshold_GeV_at_axion_point=max_mass,
                minimum_inverse_Y_plus_2_shift_at_high_scale=4*cost/3/(2*math.pi)*max(0,math.log(U0/max_mass)),
                axion_mass_ceiling_if_all_above_high_scale_micro_eV=cap*D*(k0+dk)/(math.sqrt(2)*U0)*1.e15,
                crossings=rows,joint_broken_phase_witness=joint)


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check",action="store_true")
    mode.add_argument("--write-new",action="store_true")
    args=parser.parse_args()
    data=json.loads(json.dumps(calculate()),parse_float=lambda x:float(format(float(x),".11g")))
    if args.write_new:
        with (BASE/"results.json").open("x",encoding="utf-8") as fp:
            json.dump(data,fp,ensure_ascii=False,indent=2)
            fp.write("\n")
    else:
        assert data == read("results.json")
        print("1472: anomaly cost, threshold bound and all eight gauge input corners verified")
