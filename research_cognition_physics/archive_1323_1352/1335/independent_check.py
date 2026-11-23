"""Round 1335: independent finite checks for a formed common left code.

No imports from the author. By default compare to saved JSON; pass --write for
the initial result. All observables here use a declared constant Higgs field
or a separately controlled finite neutral measurement protocol.
"""
import argparse
import json
from fractions import Fraction
from pathlib import Path

import numpy as np


def dagger(matrix):
    return matrix.conj().T


def spin_pulse(j, theta):
    n = int(2 * j + 1)
    m = np.arange(n) - j
    raising = np.diag(np.sqrt((j - m[:-1]) * (j + m[:-1] + 1)), -1)
    jz = np.diag(m)
    generator = np.block([[jz, raising.T], [raising, -jz]])
    eig, basis = np.linalg.eigh(generator)
    return (basis * np.exp(-1j * theta * eig)) @ basis.T


def port(chi, unitary):
    n = len(unitary) // 2
    # chi is a ket. The portal acting on the left is a row.
    bra = np.kron(np.conj(chi)[None, :], np.eye(n))
    return bra @ unitary[:, :n], bra @ unitary[:, n:]


def masses_and_left_vectors(unitary, chi, t, heavy_mass):
    n = len(unitary) // 2
    # Right rows: ordinary singlet, heavy doublet. All left columns are
    # doublets. The bottom row measures the orthogonal complement of W.
    matrix = np.vstack((t * np.kron(np.conj(chi)[None, :], np.eye(n)),
                        heavy_mass * dagger(unitary[:, n:])))
    right, mass, left_dag = np.linalg.svd(matrix)
    order = np.argsort(mass)
    return mass[order], dagger(left_dag)[:, order], right[:, order]


def reduced_masses(singular, t, heavy_mass):
    total = t*t + heavy_mass*heavy_mass
    prod = t*t * heavy_mass*heavy_mass * singular**2
    # Form the small root without cancellation.
    large = (total + np.sqrt(np.maximum(total*total - 4*prod, 0))) / 2
    small = prod / large
    return np.sqrt(np.sort(np.concatenate((small, large))))


def tr(matrix):
    return np.trace(matrix)


class TestLog:
    def __init__(self):
        self.residuals = {}
        self.counts = {}

    def close(self, name, actual, target, tolerance=3e-11):
        difference = np.asarray(actual) - np.asarray(target)
        residual = float(np.max(np.abs(difference)) /
                         max(1., np.max(np.abs(target))))
        if not np.isfinite(residual) or residual > tolerance:
            raise AssertionError((name, residual, tolerance))
        self.counts[name] = self.counts.get(name, 0) + 1
        self.residuals[name] = max(self.residuals.get(name, 0.), residual)


def exact_effect_numerator(a, b_real, b_imag):
    # theta=pi/3, omit global phase. Basis order m=-1,0,1.
    # Bra spinor = (a,b)/d. L0=diag(-1,1,3)/3 and L1=2 sqrt(2) S/3.
    # Similarity with diag(1,sqrt(2),1) makes L^*L Gaussian rational.
    # Use an integer realification [[Re,-Im],[Im,Re]] multiplied by 9 d^2.
    diag_short = a*a + 8*(b_real*b_real + b_imag*b_imag)
    real = np.array([[diag_short, 4*a*b_real, 0],
                     [2*a*b_real, diag_short, 6*a*b_real],
                     [0, 12*a*b_real, 9*a*a]], dtype=object)
    imag = np.array([[0, -4*a*b_imag, 0], [2*a*b_imag, 0, -6*a*b_imag],
                     [0, 12*a*b_imag, 0]], dtype=object)
    return np.block([[real, -imag], [imag, real]])


def exact_observables(flip):
    d_u, d_d = 5, 65
    u = exact_effect_numerator(3, 4, 0)
    d = exact_effect_numerator(25, flip*36, -flip*48)
    den_u, den_d = 9*d_u**2, 9*d_d**2
    comm = u@d - d@u
    cube = comm@comm@comm
    cp = Fraction(int(tr(cube[3:, :3])), (den_u*den_d)**3)
    three_yes = Fraction(int(tr((u@u@d)[:3, :3])), 3*den_u**2*den_d)
    return cp, three_yes


def neutral_history(unitary, chi_u, chi_d):
    # Multiply actual ordered projections and pulse matrices in the
    # 6-dimensional space. A fresh record bit can store each yes/no flag.
    # The squared norm here is unconditional on selection; failure branches
    # remain in the instrument. Ramified history bits are accounted below.
    n = len(unitary)//2
    identity = np.eye(2*n)
    p0 = np.diag(np.r_[np.ones(n), np.zeros(n)])
    projector = lambda c: np.kron(np.outer(c, np.conj(c)), np.eye(n))
    code = unitary[:, :n]
    history = {}
    for i, p in enumerate((identity-projector(chi_u), projector(chi_u))):
        for j, q in enumerate((identity-p0, p0)):
            for k, r in enumerate((identity-projector(chi_d), projector(chi_d))):
                block = r @ unitary @ q @ dagger(unitary) @ p @ code
                history[f"{i}{j}{k}"] = float(np.linalg.norm(block, "fro")**2 / n)
    return history


def rational_witness(log):
    unitary = spin_pulse(1, np.pi/3)
    chi_u = np.array([3, 4], dtype=complex)/5
    answers = {}
    up_mass, up_left, _ = masses_and_left_vectors(unitary, chi_u, .2, 7)
    lu, _ = port(chi_u, unitary)
    eu = dagger(lu) @ lu
    for label, flip in (("shared_frame", 1), ("independent_down_pi_offset", -1)):
        chi_d = np.array([25, flip*(36+48j)])/65
        ld, _ = port(chi_d, unitary)
        ed = dagger(ld) @ ld
        comm = eu@ed-ed@eu
        cp, yes = exact_observables(flip)
        log.close("exact_cp_witness", tr(comm@comm@comm), 1j*float(cp))
        history = neutral_history(unitary, chi_u, chi_d)
        log.close("neutral_record_three_yes", history["111"], float(yes))
        log.close("neutral_history_normalized", sum(history.values()), 1.)
        down_mass, down_left, _ = masses_and_left_vectors(unitary, chi_d, .3, 7)
        current = dagger(up_left) @ down_left
        light = current[:3, :3]
        quartet = (light[0,0]*light[1,1]*light[0,1].conj()*light[1,0].conj()).imag
        log.close("complete_left_weak_current", current @ dagger(current), np.eye(6))
        if abs(quartet) < 1e-6:
            raise AssertionError("degenerate finite quartet witness")
        answers[label] = dict(cp_trace_imag=str(cp), joint_yes=str(yes),
                              joint_yes_decimal=float(yes), history=history,
                              masses_up=up_mass.tolist(), masses_down=down_mass.tolist(),
                              charged_current_abs_light=np.abs(light).tolist(),
                              cp_quartet_finite=float(quartet))
    baseline, competing = answers.values()
    log.close("separate_frame_spectra_equal", baseline["masses_down"],
              competing["masses_down"])
    if abs(baseline["joint_yes_decimal"] - competing["joint_yes_decimal"]) < .05:
        raise AssertionError("record discriminator missing")
    return answers


def full_checks(log):
    low_min, formed_min = 1., 1.
    samples = 0
    probe = np.array([1., 0., 0.])
    for theta in (0., .17, .63, np.pi/3, 1.85, 2*np.pi/3):
        u = spin_pulse(1, theta)
        w = u[:, :3]
        log.close("formed_code_isometry", dagger(w) @ w, np.eye(3))
        # The maximal weight state is m=+1 (last index).
        for angle in (0., .21, .57, 1.14, np.pi/2):
            for phase in (.0, .44, 1.23):
                chi = np.array([np.cos(angle), np.sin(angle)*np.exp(1j*phase)])
                other = np.array([-np.conj(chi[1]), np.conj(chi[0])])
                l, b = port(chi, u)
                ell, _ = port(other, u)
                log.close("row_isometry",l@dagger(l)+b@dagger(b),np.eye(3))
                e, ec = dagger(l)@l, dagger(ell)@ell
                log.close("complementary_formed_effects",e+ec,np.eye(3))
                log.close("complementary_commuting_leading_H",e@ec,ec@e)
                p = abs(chi[0])**2
                log.close("extremal_probe",np.real(e[2,2]),p)
                k0, _ = port((1.,0.),u)
                log.close("triangular_determinant",np.linalg.det(l),
                          np.conj(chi[0])**3*np.linalg.det(k0))
                s = np.linalg.svd(l,compute_uv=False)
                for t, heavy in ((.2,7.), (.88,1.), (1.,1.32)):
                    mass, _, _ = masses_and_left_vectors(u,chi,t,heavy)
                    log.close("complete_mass_pairing",mass,reduced_masses(s,t,heavy))
                    inferred = mass[:3]**2*(t*t+heavy*heavy-mass[:3]**2)/(t*t*heavy*heavy)
                    log.close("inverted_light_masses",inferred,np.sort(s*s))
                    if s[0] > 1e-12:
                        formed_min=min(formed_min,float(s[1]/s[0]))
                        low_min=min(low_min,float(mass[1]/mass[2]))
                        if s[1]/s[0]<1/3-1e-11 or mass[1]/mass[2]<1/(3*np.sqrt(2))-1e-11:
                            raise AssertionError("violated analytic hierarchy bound")
                    samples += 1
    return dict(cases=samples,formed_ratio_min=formed_min,finite_ratio_min=low_min)


def skew_checks(log):
    rng = np.random.default_rng(1335)
    for _ in range(32):
        a, x, y, z = rng.normal(size=4)+1j*rng.normal(size=4)
        v=np.array([x,y,z])
        c=np.array([[0.,-z,y],[z,0.,-x],[-y,x,0.]])
        log.close("cross_gram",dagger(c)@c,np.vdot(v,v)*np.eye(3)-np.outer(v,v.conj()))
        s=np.linalg.svd(a*np.eye(3)+c,compute_uv=False)
        if s[0]>3*s[1]+1e-12:
            raise AssertionError("complex spin-one pair bound")


def calculate():
    log = TestLog()
    witness = rational_witness(log)
    spectra = full_checks(log)
    skew_checks(log)
    return dict(status="PASS",round=1335,scope="finite spin-one common code; tree-level mass; switched neutral protocol",
                witness=witness,spectral_samples=spectra,counts=log.counts,
                max_residuals=log.residuals)


def compare(actual, reference):
    if isinstance(reference,dict):
        if actual.keys()!=reference.keys():
            raise AssertionError("JSON keys differ")
        for key in actual:
            compare(actual[key],reference[key])
    elif isinstance(reference,list):
        np.testing.assert_allclose(actual,reference,rtol=2e-9,atol=2e-11)
    elif isinstance(reference,float):
        np.testing.assert_allclose(actual,reference,rtol=2e-9,atol=2e-11)
    elif actual!=reference:
        raise AssertionError((actual,reference))


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--write",action="store_true")
    args=parser.parse_args()
    outfile=Path(__file__).with_name("independent_results.json")
    result=calculate()
    if args.write:
        outfile.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    else:
        compare(result,json.loads(outfile.read_text(encoding="utf-8")))
    print(json.dumps(dict(status=result["status"],counts=result["counts"],spectrum=result["spectral_samples"],
                          joint_yes=[w["joint_yes"] for w in result["witness"].values()],
                          finite_quartet=[w["cp_quartet_finite"] for w in result["witness"].values()])))
