"""One-loop threshold geometry for the two charged-lepton radiative multiplets."""
import argparse
import itertools
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
TAU = 2*math.pi


def read(name):
    return json.loads((HERE / name).read_text(encoding="utf-8-sig"))


def evaluate(couplings, spec, mz):
    A = [4*math.pi/g**2 for g in couplings]
    x, y = A[0]-A[1], A[1]-A[2]
    b1, b2, b3 = spec["SM_b"]
    B12, B23 = b1-b2, b2-b3
    ratio = x/y
    bare = B12-ratio*B23
    efficiency = 1/15+ratio/6
    r0 = B12/B23
    C = 1/15 + r0/6
    gap = r0*y - x
    W0 = gap/C
    T0 = (y-W0/6)/B23
    def mass(log_window):
        return mz*math.exp(TAU*log_window)
    minimum_spans = {}
    for n in spec["ordinary_copies"]:
        r = (B12-n/15)/(B23+n/6)
        minimum_spans[str(n)] = max(0., (r*y-x)/(1+r))
    n = spec["example"]["doublets"]
    earliest_t = y/(B23+n/6)
    max_S = C*n*earliest_t-gap
    return {"rho":ratio, "split_continuous_count":bare/efficiency,
            "intact_positive_continuous_count":bare/(efficiency-.2),
            "negative_intact_efficiency":efficiency-.8,
            "span":minimum_spans, "W0":W0, "T0":T0,
            "seven_scale_min_GeV":mass(earliest_t),"seven_scale_max_GeV":mass(T0),
            "seven_doublet_geomean_cap_GeV":mass(T0-W0/n),
            "seven_positive_singlet_geomean_floor_GeV":mass(earliest_t-max_S/(.2*n))}


def positive_witness(gauges, cfg, mz):
    A = [4*math.pi/g**2 for g in gauges]
    x, y = A[0]-A[1], A[1]-A[2]
    b1, b2, b3 = cfg["SM_b"]
    D12, D23 = b1-b2, b2-b3
    n = cfg["example"]["doublets"]
    md = cfg["example"]["common_unbroken_doublet_GeV"]
    ell = math.log(md/mz)/TAU
    t = (y+n*ell/6)/(D23+n/6)
    w = n*(t-ell)
    S = x-D12*t+w/15
    s_single = S/(.2*n)
    grand = mz*math.exp(TAU*t)
    ms = grand*math.exp(-TAU*s_single)
    assert mz<md<ms<grand and 0<s_single<t
    invU = A[2]-b3*t
    rebuilt = [invU+b1*t+.1*w+.2*n*s_single,
               invU+b2*t+w/6,invU+b3*t]
    assert all(abs(x-y)<2e-12 for x,y in zip(A,rebuilt))
    gauge_rows = []
    bd = cfg["gauge_split"]["doublet"]
    bs = cfg["gauge_split"]["positive_charged_scalar"]
    for point in [mz,md,ms,grand]:
        def length(threshold):
            return max(0.,math.log(point/threshold)/TAU)
        inv = [a-b*length(mz)-n*d*length(md)-n*s*length(ms)
               for a,b,d,s in zip(A,cfg["SM_b"],bd,bs)]
        assert min(inv)>0
        gauge_rows.append({"scale_GeV":point,"inverse_alphas":inv})
    # Exact two-state mass eigenvalues of an explicitly nonnegative gauge-invariant potential.
    loop = read(cfg["muon_loop"])
    v = loop["v_GeV"]
    trilinear = cfg["example"]["trilinear_over_singlet"]*ms
    off = trilinear*v/math.sqrt(2)
    hd, hs = md**2+off**2/ms**2, ms**2
    split = math.hypot(hs-hd,2*off)
    large = (hs+hd+split)/2
    small = md**2*ms**2/large
    M = math.sqrt(large)
    xx = small/large
    sin2 = 2*off/split
    diff = 1-xx*math.log(xx)/(xx-1)
    # Q=1: Baker et al Eq.7, second scalar x2=1 derivative limits.
    loop_q1 = ((3*xx-1)/(1-xx)**2+2*xx**2*math.log(xx)/(1-xx)**3+
               2*(1/(1-xx)+xx*math.log(xx)/(1-xx)**2))
    kernel = (loop_q1-1/3)/(diff*M**2)
    flavours = {}
    for name, mass in [("electron",cfg["electron_mass_GeV"]),("muon",loop["muon_mass_GeV"])]:
        prod = 32*math.pi**2*mass/(M*sin2*diff)
        assert prod<1
        flavours[name] = {"Yukawa_product":prod,"leading_magnetic_moment":mass**2*kernel}
    # Three flavours, five active labels: the tau mass adds coherently from three tags.
    tau_mass = cfg["tau_mass_GeV"]
    taur = 32*math.pi**2*tau_mass/(3*M*sin2*diff)
    assert 0 < taur < 4*math.pi
    three = {"electron_labels":1,"muon_labels":1,"tau_labels":3,"spectator_labels":2,
             "tau_product_per_tag":taur,"equal_tau_Yukawas_per_tag":math.sqrt(taur),
             "tau_mass_check_GeV":3*taur*M*sin2*diff/(32*math.pi**2),
             "light_relative_mass_error_bound":(trilinear*v/math.sqrt(2)/ms**2)**2}
    return {"three_flavour_option":three,"gauge_running":gauge_rows,"neutral_fermion_mass_GeV":M,
            "physical_charged_scalar_masses_GeV":[math.sqrt(small),M],
            "physical_sin_2theta":sin2,"trilinear_GeV":trilinear,
            "diagonal_flavour_inputs_fitted":flavours}


def calculate():
    cfg = read("inputs.json")
    g = read(cfg["gauge"])
    gs = [g[k] for k in ("g1_GUT","g2","g3")]
    mz = g["mu_ref_GeV"]
    middle = [p[0] for p in gs]
    corners = [evaluate([v+e*cfg["working_widths"]*s for (v,e),s in zip(gs,signs)],cfg,mz)
               for signs in itertools.product((-1,1),repeat=3)]
    center = evaluate(middle,cfg,mz)
    def extent(vals):
        return [min(vals),max(vals)]
    boxes = {key:extent([r[key] for r in corners]) for key in center if key!="span"}
    boxes["span"] = {str(n):extent([r["span"][str(n)] for r in corners])
                     for n in cfg["ordinary_copies"]}
    assert 6<boxes["split_continuous_count"][0]<=boxes["split_continuous_count"][1]<7
    assert boxes["negative_intact_efficiency"][1]<0
    witness = positive_witness(middle,cfg,mz)
    return {"central":center,"two_width_working_box":boxes,"seven_positive_witness":witness}


def match(a,b):
    if isinstance(a,dict):
        return set(a)==set(b) and all(match(v,b[k]) for k,v in a.items())
    if isinstance(a,list):
        return len(a)==len(b) and all(match(v,w) for v,w in zip(a,b))
    return math.isclose(a,b,rel_tol=5e-11,abs_tol=1e-34)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=calculate()
    if args.check:
        if not match(result,read("results.json")):
            raise SystemExit("1592 threshold results differ")
        print("1592: split multiplicities, finite matching costs, common mass witness all verified")
    else:
        (HERE/"results.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
