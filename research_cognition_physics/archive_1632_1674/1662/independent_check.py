"""Independent interval, rectangular-matrix and explicit threshold cross-checks."""
import argparse
from fractions import Fraction
import importlib.util
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("outward1486", HERE.parent / "1486/check.py")
iv = importlib.util.module_from_spec(spec)
spec.loader.exec_module(iv)
iv.SCALE = 10**40
I = iv.I
qpi = Fraction("3.1415926535897932384626433832795028841971")
pi = I(I.point(qpi).lo, I.point(qpi + Fraction(1, 10**40)).hi)


def load(name):
    return json.loads((HERE.parent / name).read_bytes())


def trig(angle, sine):
    def endpoints(n):
        x = Fraction(n, iv.SCALE)
        order = int(sine)
        lower = sum((-1)**j*x**(order+2*j)/math.factorial(order+2*j) for j in range(14))
        return I.point(lower).lo, I.point(lower+x**(order+28)/math.factorial(order+28)).hi
    assert 0 <= angle.lo <= angle.hi < pi.lo // 2
    at_lo, at_hi = endpoints(angle.lo), endpoints(angle.hi)
    return I(at_lo[0], at_hi[1]) if sine else I(at_hi[0], at_lo[1])


def flavour():
    down, strange, bottom, electron, muon, tau = map(iv.width, load("1402/inputs.json")["scales"]["1e16"])
    angles = list(map(iv.width, load("1397/inputs.json")["high_scale"]["angles_theta12_theta13_theta23_delta"]))
    s12, s13, s23, _ = [trig(x, True) for x in angles]
    c12, c13, c23, cd = [trig(x, False) for x in angles]
    k_min = (electron*muon/(down*strange)*c13*c23).sqrt()
    y2d, y2s, y2b = [y.square() for y in (down, strange, bottom)]
    x = c13*c23*s12*c12*(y2s-y2d)
    y = c13*s23*s13*(y2b-y2d*c12.square()-y2s*s12.square())
    assert min(x.lo, y.lo, cd.lo) > 0
    h12 = (x.square()+y.square()+2*x*y*cd).sqrt()
    a_upper = (y2d+y2s*s12.square()+y2b*s13.square()).sqrt()
    k_max = tau*a_upper/(3*h12)
    k_limits = [k_min.values()[0], k_max.values()[1]]
    assert k_limits[0] > 1.30003 and k_limits[1] < 15.701
    return k_limits


def gauge(ks):
    src = load("1369/inputs.json")
    a = [4*pi/iv.width(src[key]).square() for key in ("g1_GUT", "g2", "g3")]
    b = [Fraction(41,10),Fraction(-19,6),Fraction(-7)]
    triplet = [Fraction(4,15),Fraction(0),Fraction(2,3)]
    doublet = [Fraction(2,5),Fraction(2,3),Fraction(0)]
    projection = [5,-3,-2]
    assert all(sum(x*y for x,y in zip(projection,col)) == expected for col,expected in
               ((b,44),(triplet,0),(doublet,0),([1,1,1],0)))
    r = (b[2]-b[1])/(b[0]-b[1])
    factor = 1-Fraction(2,5)*r
    u = a[2]-a[1]-I.point(r)*(a[0]-a[1])
    log_mu = (sum(w*ai for w,ai in zip(projection,a))*pi*2/44).values()
    log_k = (-pi*u/factor).values()
    mu = [src["mu_ref_GeV"]*math.exp(x) for x in log_mu]
    k_needed = [math.exp(x) for x in log_k]
    u_model = sorted(-float(factor)*math.log(k)/math.pi for k in ks)
    gap = u_model[0]-u.values()[1]
    assert gap > 4.375 and k_needed[0] > 1_337_100
    assert 4.75e13 < mu[0] and mu[1] < 5.26e13
    limits = {"U_observed":u.values(),"U_allowed_by_flavour":u_model,"projection_gap_at_least":gap,
              "gauge_required_heavy_ratio_per_family":k_needed, "threshold_free_scale_GeV":mu}
    point = load("1657/independent_results.json")["extended_box_GJ_test"]
    k = point["k_times_yb"]/point["marginal_coordinates"][2]
    central = np.array([4*math.pi/src[key][0]**2 for key in ("g1_GUT","g2","g3")])
    u0 = central[2]-central[1]-float(r)*(central[0]-central[1])
    heavy_ratio = math.exp(-math.pi*u0/float(factor))
    m_l, mz = 1e5, src["mu_ref_GeV"]
    m_t = m_l*heavy_ratio
    scale = mz*math.exp(2*math.pi*np.dot(projection,central)/44)
    assert mz < m_l < m_t < scale
    b, triplet, doublet = [np.asarray(col,dtype=float) for col in (b,triplet,doublet)]
    unified = central - (b*math.log(scale/mz)+3*triplet*math.log(scale/m_t)
                         +3*doublet*math.log(scale/m_l))/(2*math.pi)
    mismatch = float(np.ptp(unified))
    assert mismatch < 1e-11 and min(unified) > 0
    return limits, {"k_flavour":k,"bare_mass_doublet_GeV":m_l,"bare_mass_triplet_GeV":m_t/k,
                    "physical_doublet_GeV":m_l,"physical_triplet_GeV":m_t,
                    "unification_GeV":scale,"common_inverse_alpha":float(np.mean(unified)),
                    "max_gauge_splitting":mismatch}


def rect_test():
    rng = np.random.default_rng(1662)
    def draws(rows, cols):
        return rng.normal(size=(rows,cols))+1j*rng.normal(size=(rows,cols))
    residuals = []
    for n in (1,2,3,4,5):
        bare = 4*np.eye(n)+draws(n,n)
        for _ in range(3):
            mix = draws(n,3)
            a = np.linalg.solve(bare,mix)
            eig, vec = np.linalg.eigh(np.eye(3)+a.conj().T@a)
            overlap = (vec/np.sqrt(eig))@vec.conj().T
            embed = np.vstack([np.eye(3),-a])@overlap
            mass = np.column_stack([mix,bare])
            prod = np.prod(np.linalg.svd(mass,compute_uv=False))
            checks = [np.linalg.norm(mass@embed), np.linalg.norm(embed.conj().T@embed-np.eye(3)),
                      abs(prod*np.linalg.det(overlap)/abs(np.linalg.det(bare))-1)]
            if n < 3:
                checks.append(abs(np.max(np.linalg.eigvalsh(overlap))-1))
            residuals.extend(map(float,checks))
    assert max(residuals) < 1e-11
    return {"seed":1662,"dimensions":[1,2,3,4,5],"samples_per_dimension":3,"checks":len(residuals),
            "matrix_tests_pass":True}


def report():
    ks = flavour()
    bounds, counter = gauge(ks)
    return {"two_width_flavour_gain":ks, "one_loop":bounds,
            "independent_hard_masses_example":counter,"rectangular_complex_matrix_checks":rect_test()}


def tidy(value):
    if isinstance(value,dict):
        return {k:tidy(v) for k,v in value.items()}
    if isinstance(value,list):
        return list(map(tidy, value))
    if isinstance(value,float):
        return 0.0 if abs(value) < 1e-12 else float(f"{value:.11g}")
    return value


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    opts = parser.add_mutually_exclusive_group(required=True)
    opts.add_argument("--check",action="store_true")
    opts.add_argument("--write-new",action="store_true")
    args = parser.parse_args()
    path = HERE/"independent_results.json"
    result = tidy(report())
    if args.write_new:
        with path.open("x",encoding="utf-8") as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    else:
        assert result == json.loads(path.read_bytes())
        print("PASS 1662 independent interval, rank and threshold checks")
