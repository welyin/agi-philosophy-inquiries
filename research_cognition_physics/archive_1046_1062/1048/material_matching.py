"""1048: conditional microscopic dipole-to-material task matching.

TRK and gauge truncation are established tools. The new project bridge is a
uniform finite-history/instrument bound, including the old full charge poststate.
No microscopic Hubbard embedding or unbounded photon/source error is certified.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import importlib.util
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULT = HERE / "material_matching_results.json"
HISTORY = [
    "archive_956_989/956/native_material_interface.py",
    "archive_956_989/research_note_956.md",
    "archive_956_989/research_note_960.md",
    "archive_956_989/research_note_965.md",
    "archive_702_741/research_note_712.md",
]


def norm(a):
    return float(np.linalg.norm(a, 2))


def evolution(h, t):
    e, v = np.linalg.eigh(h)
    return (v * np.exp(-1j * t * e)) @ v.conj().T


def comm(a, b):
    return a @ b - b @ a


def load_material():
    path = ROOT / HISTORY[0]
    spec = importlib.util.spec_from_file_location("material956_1048", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def exact_virtual_identity():
    """General algebra test, not a canonical NR microscopic parent."""
    energies = [F(0), F(1), F(2), F(10), F(11)]
    x = [[F(0) if i == j else F((i + j) % 3 + 1, 7)
          for j in range(5)] for i in range(5)]
    errors = []
    offdiagonal_heavy = []
    for m in range(3):
        for n in range(3):
            total = sum((2 * energies[k] - energies[m] - energies[n])
                        * x[m][k] * x[k][n] for k in range(5))
            retained = sum((2 * energies[k] - energies[m] - energies[n])
                           * x[m][k] * x[k][n] for k in range(3))
            # j_mh j_hn = (E_h-E_m)(E_h-E_n) x_mh x_hn.
            virtual = sum(
                (energies[k] - energies[m]) * (energies[k] - energies[n])
                * x[m][k] * x[k][n]
                * (1 / (energies[m] - energies[k])
                   + 1 / (energies[n] - energies[k]))
                for k in range(3, 5))
            errors.append(total + virtual - retained)
            if m != n:
                offdiagonal_heavy.append(virtual)
    assert all(e == 0 for e in errors)
    assert all(e != 0 for e in offdiagonal_heavy)
    return {"exact_matrix_entries": len(errors),
            "nonzero_offdiagonal_virtual_entries": len(offdiagonal_heavy)}


def displacement_entry(n, m, beta):
    total = sum((-beta.conjugate()) ** (m-k) * beta ** (n-k)
                / (math.factorial(m-k) * math.factorial(n-k) * math.factorial(k))
                for k in range(min(n, m)+1))
    return (math.exp(-abs(beta)**2/2)
            * math.sqrt(math.factorial(m)*math.factorial(n))*total)


def exact_oscillator_checks():
    """Exact infinite-oscillator evolution compressed to three input levels.

Only the 3x3 cross-Gram matrix is needed to compute the full isometry distance;
there is no Fock cutoff used to approximate evolution.
    """
    d = 3
    h = np.diag(np.arange(d)+.5)
    x = np.zeros((d, d))
    for n in range(d-1):
        x[n, n+1] = x[n+1, n] = math.sqrt((n+1)/2)
    k = comm(x, comm(h, x))
    residual = np.eye(d)-k
    assert norm(residual-np.diag([0, 0, 3])) < 2e-15
    b = math.sqrt(1.5)
    rows = []
    for force, time in ((.03, .1), (.05, .5), (.04, 1.0)):
        alpha = force/math.sqrt(2)
        beta = alpha*(np.exp(-1j*time)-1)
        phase = np.exp(1j*alpha**2*(time-math.sin(time)))
        full_pp = np.array([[displacement_entry(n, m, beta)
                             for m in range(d)] for n in range(d)])
        full_pp = phase*full_pp@evolution(h, time)
        effective = evolution(h+force*x, time)
        difference_gram = 2*np.eye(d)-full_pp.conj().T@effective-effective.conj().T@full_pp
        actual = math.sqrt(max(0., float(np.linalg.eigvalsh(difference_gram)[-1])))
        upper = b*abs(force)*time
        assert actual <= upper+2e-11
        annihilator = np.zeros((d,d),complex)
        for n in range(1,d):
            annihilator[n-1,n] = math.sqrt(n)
        rotation = evolution(h,time)
        energy_operator = h+rotation.conj().T@(
            beta*annihilator.conj().T+beta.conjugate()*annihilator
            )@rotation+abs(beta)**2*np.eye(d)
        max_bare_energy = float(np.linalg.eigvalsh(energy_operator)[-1])
        energy_upper = (math.sqrt(2.5)+abs(force)*time/math.sqrt(2))**2
        assert max_bare_energy <= energy_upper+1e-13
        rows.append({"force":force,"time":time,"exact_infinite_isometry_error":actual,
                     "trk_duhamel_upper":upper,
                     "exact_max_bare_energy":max_bare_energy,
                     "conditional_energy_upper":energy_upper})
    return {"rank":d,"C":1,"heavy_gap":1,
            "off_block_norm_squared":1.5,"trk_bound_squared":1.5,
            "finite_history_checks":rows}


def run():
    model = load_material()
    material = model.material(1., .1)
    h, D = material["H"], material["D"]
    cs = [model.annihilate(i) for i in range(4)]
    ns = [c.T@c for c in cs]
    sector = np.eye(16)[:, [s for s in range(16) if s.bit_count() == 2]]
    x = sector.T@((ns[0]+ns[1]-ns[2]-ns[3])/2)@sector
    assert norm(x@x-D) < 1e-15
    j = 1j*comm(h, x)
    k = comm(x, comm(h, x))
    assert norm(k-.1*material["hop"]) < 1e-15
    assert abs(np.trace(k)) < 1e-15
    e, v = np.linalg.eigh(h)
    f_s = float(v[:,0]@k@v[:,0])
    Delta = 1+material["J"]
    formula = 2*Delta*material["occupancy"]
    assert abs(f_s-formula) < 2e-15
    C = 2.
    residual_diagonal = C-np.diag(v.T@k@v)
    assert residual_diagonal.min() > 0
    assert abs(residual_diagonal.sum()-6*C) < 2e-14
    triplet_dark = norm(x@material["W"]@(np.eye(4)-material["ps"]))
    assert triplet_dark < 1e-14
    gauge = []
    for A in (.1, .5, 1.):
        G = evolution(x, -A)
        right = G@h@G.conj().T
        wrong = h-A*j+.5*C*A*A*np.eye(6)
        # A^2 scalar cannot repair the changed singlet-triplet energy gap.
        wrong_gap = (math.sqrt(1+.16*(1+A*A))-1)/2
        wrong_spectrum = np.linalg.eigvalsh(wrong)-.5*C*A*A
        assert abs(wrong_spectrum[0]+wrong_gap) < 2e-14
        gauge.append({"A":A,"correct_spectrum_max_error":float(np.max(abs(np.linalg.eigvalsh(right)-e))),
                      "bare_contact_false_gap_change":wrong_gap-material["J"],
                      "same_source_transport_error":norm(-G@j@G.conj().T-1j*comm(x,right))})
    sign = evolution(x, math.pi)
    K0, K1 = (np.eye(6)+sign)/2, (np.eye(6)-sign)/2
    assert norm(K0-(np.eye(6)-D)) < 1e-14 and norm(K1-D) < 1e-14
    target = np.vstack([np.eye(6)-D,D])
    center = (e[-1]+e[0])/2
    radius = (e[-1]-e[0])/2
    meter = []
    for tau in (.001, .01, .1):
        branch0 = evolution(h-center*np.eye(6),tau)
        branch1 = evolution(h-center*np.eye(6)+math.pi*x/tau,tau)
        actual = np.vstack([(branch0+branch1)/2,(branch0-branch1)/2])
        error = norm(actual-target)
        assert norm(actual.conj().T@actual-np.eye(6)) < 3e-14
        assert error <= tau*radius+2e-14
        meter.append({"duration":tau,"complete_instrument_isometry_error":error,
                      "internal_drift_upper":tau*radius})
    W = material["W"]
    outside = np.eye(6)-W@W.T
    leakage = W.T@(D@outside@D+(np.eye(6)-D)@outside@(np.eye(6)-D))@W
    assert norm(leakage-2*material["occupancy"]*(1-material["occupancy"])*material["ps"]) < 2e-14
    # This is a sufficient hypothetical certificate; the old material has no
    # established microscopic spectral/dipole embedding or heavy gap.
    gap = 1e8
    beta = math.sqrt(6*C/(2*gap))
    certificate = beta*(math.pi+.01)+.001*radius
    assert certificate < .00132
    frequency_rows = []
    for gap in (10.,100.):
        for omega in (0.,.5):
            r = C-f_s
            # A single positive residual line at the lower allowed gap
            # saturates this scalar spectral bound; not a whole QED completion.
            heavy_response = r/(gap*gap-omega*omega)
            frequency_rows.append({"heavy_gap":gap,"frequency":omega,
                                   "singlet_remainder_upper":heavy_response,
                                   "triplet_remainder_upper":C/(gap*gap-omega*omega)})
    result = {
        "round":1048,"new_scientific_groups":1,"new_cognitive_axioms":0,
        "all_scientific_checks_passed":True,
        "old_material":{"U":1.,"v":.1,"rank":6,"electron_count":2,
                        "C_diagnostic_units":C,"singlet_retained_strength":f_s,
                        "singlet_strength_formula":formula,
                        "triplet_retained_strength":0.,
                        "required_heavy_strength_diagonal":residual_diagonal.tolist(),
                        "total_required_heavy_strength":float(residual_diagonal.sum()),
                        "four_dimensional_charge_poststate_leakage":2*material["occupancy"]*(1-material["occupancy"])},
        "exact_virtual_matrix_check":exact_virtual_identity(),
        "canonical_oscillator_certificate":exact_oscillator_checks(),
        "gauge_checks":gauge,"finite_meter_checks":meter,
        "conditional_microscopic_certificate_example":{
            "heavy_gap":1e8,"extra_field_area":.01,"meter_duration":.001,
            "isometry_error_upper":certificate,
            "microscopic_material_embedding_certified":False},
        "subgap_linear_response_bounds":frequency_rows,
        "scope":{"parent":"Nonrelativistic local-potential canonical electrons with stated domains and spectral gap",
                 "uniform_tasks":"Finite low-input histories with bounded auxiliary couplings to one dipole component",
                 "unbounded_cavity_coupling_certified":False,
                 "unbounded_stress_error_certified":False,
                 "full_P981_or_GR_or_roadmap_completed":False},
        "historical_sha256":{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in HISTORY},
    }
    return result


def compare(a,b,path="root"):
    if isinstance(a,dict):
        assert a.keys()==b.keys(),path
        for k in a: compare(a[k],b[k],path+"."+k)
    elif isinstance(a,list):
        assert len(a)==len(b),path
        for i,(x,y) in enumerate(zip(a,b)): compare(x,y,path+"."+str(i))
    elif isinstance(a,float):
        assert math.isclose(a,b,rel_tol=2e-11,abs_tol=2e-11),(path,a,b)
    else:
        assert a==b,(path,a,b)


def main():
    p=argparse.ArgumentParser();p.add_argument("--write",action="store_true");args=p.parse_args()
    result=run()
    if args.write:
        with RESULT.open("x",encoding="utf-8") as f:
            json.dump(result,f,ensure_ascii=False,indent=2);f.write("\n")
    else:
        compare(result,json.loads(RESULT.read_text(encoding="utf-8")))
    print(json.dumps({"round":1048,"all_scientific_checks_passed":True,
                      "mode":"exclusive_write" if args.write else "read_only_compare",
                      "singlet_retained_strength":result["old_material"]["singlet_retained_strength"],
                      "meter_max_error":max(x["complete_instrument_isometry_error"] for x in result["finite_meter_checks"]),
                      "canonical_checks":result["canonical_oscillator_certificate"]["finite_history_checks"],
                      "scope":result["scope"]},ensure_ascii=False))


if __name__=="__main__": main()
