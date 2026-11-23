"""Same physical spectra for canonical qque matching, gauge running, and proton data."""
import argparse
import json
import math
from pathlib import Path
import runpy

import numpy as np

HERE = Path(__file__).resolve().parent
def read(path):
    return json.loads((HERE/path).read_bytes())

CONFIG = read("inputs.json")
BASE = read(CONFIG["baseline"])["shared_witness"]
REFERENCE = read(CONFIG["gauge_result"])[CONFIG["spectrum_case"]]
GAUGE = runpy.run_path(str(HERE/"../1696/check.py"))
SM, SUSY = GAUGE["SM"], GAUGE["MSSM"]
COLOR, WEAK = GAUGE["COLOR"], GAUGE["WEAK"]
CONFIG_1725 = read(CONFIG["superk_and_masses"])
DATA = read(CONFIG["nucleon"])


def flavour(step):
    original = runpy.run_path(str(HERE/CONFIG["frozen_flavour"]))
    v = original["CFG"]["nearby_measured_inputs_in_marginal_widths"][CONFIG["case"]]
    _, pt, _ = original["charged_endpoint"](2, step, v["weak"], v["up"],
                                              v["g3"], v["top"], v["charm"])
    # E is the high-family Ye transpose in the down-texture frame; physical e^c lies in the 10.
    df = np.linalg.svd(pt["D"])[0][:,-1]
    ef = np.linalg.svd(pt["E"])[0][:,-1]
    vec = np.array([2,1,1])*df
    assert np.allclose(np.linalg.svd(pt["E"],compute_uv=False)[::-1],
                       pt["rep"]["Y_charged_at_10PeV"],rtol=1e-7)
    terms = ef.conj()*vec
    amps = abs(terms)
    hi = sum(amps)
    lo = max(0.,2*max(amps)-hi)
    phases = np.where(np.arange(3)==np.argmax(amps),1.,-1.)*np.exp(-1j*np.angle(terms))
    assert np.isclose(abs(terms@phases),lo)
    return dict(terms=amps.tolist(), tagged_amplitude=abs(sum(terms)), min_amplitude=lo,
                max_amplitude=hi, norm_enlarged_d=np.linalg.norm(vec))


def anomalous():
    # Positive coefficients multiply alpha_i dlog(mu)/(2pi) in downward running.
    # Ordinary SM hypercharge, then convert g_Y^2=(3/5)*g_1^2 and halve the beta.
    yq,yu,ye=1/6,2/3,-1
    g1_sm=-6*(yq*yq+yu*ye)*(3/5)/2
    g1_susy=(23/9)*(3/5)/2
    elec=np.array([[g1_sm,9/4,2],[g1_susy,3/2,4/3]])
    nu=np.array([[11/20,9/4,2],[(11/9)*(3/5)/2,3/2,4/3]])
    # Recover Eq.80's alpha1 powers in ordinary minimal SM / MSSM running.
    np.testing.assert_allclose(elec[:,0]/np.array([SM[0],SUSY[0]]), [23/82,23/198])
    np.testing.assert_allclose(nu[:,0]/np.array([SM[0],SUSY[0]]),[11/82,1/18])
    return elec,nu


def thresholds(weak_poles):
    rules=[]
    for label, charge in (("colored",COLOR),("weak",WEAK)):
        for species,weight in (("fermions",2/3),("scalars",1/3)):
            rules.extend((m,weight*charge) for m in BASE[f"{label}_{species}_GeV"])
    new=read("../1706/results.json")["cluster_reference_scales_GeV"][1]
    rules+= [(new,4*COLOR/3),(math.hypot(new,BASE["soft_GeV"]),2*COLOR/3)]
    rules.extend((m,WEAK) for m in weak_poles)
    return rules


def transport(weak_poles):
    params=read(CONFIG["gauge_inputs"])
    mz=params["mu_ref_GeV"]
    start=np.array([params[k] for k in ("g1_GUT","g2","g3")])
    inv=4*math.pi/(start[:,0]+start[:,1]*REFERENCE["signs"])**2
    soft=BASE["soft_GeV"]
    limit=REFERENCE["U"]
    rules=thresholds(weak_poles)
    grid=sorted({mz,soft,limit,*(m for m,_ in rules)})
    gamma,other=anomalous()
    contributions=np.zeros((2,3))
    quadrature=np.zeros_like(contributions)
    nodes,weights=np.polynomial.legendre.leggauss(24)
    for l,h in zip(grid,grid[1:]):
        regime=int(l>=soft)
        if l==soft:
            inv-=GAUGE["ADJ"]/(12*math.pi)
        b=(SUSY if regime else SM)+sum((w for m,w in rules if m<=l),np.zeros(3))
        dt=math.log(h/l)/(2*math.pi)
        after=inv-dt*b
        assert min(after)>1
        integrals=dt/inv
        nz=abs(b)>1e-12
        integrals[nz]=np.log1p(dt*b[nz]/after[nz])/b[nz]
        contributions[regime]+=integrals
        quadrature[regime]+=dt/2*(weights@(1/(inv-(nodes[:,None]+1)*dt*b/2)))
        inv=after
    np.testing.assert_allclose(inv,REFERENCE["new_inverse_colour"],rtol=5e-9)
    np.testing.assert_allclose(contributions,quadrature,rtol=2e-13,atol=2e-15)
    gains=np.exp(np.sum(gamma*contributions,axis=1))
    gains_other=np.exp(np.sum(other*contributions,axis=1))
    return dict(weak_poles_GeV=weak_poles,final_inverse_gauge=inv.tolist(),
                integrals_SM_SUSY=contributions.tolist(),A_qque_SM_and_SUSY=gains.tolist(),
                passes_direct_gauss_quadrature=True,
                total_qque=CONFIG["long_distance_QCD"]*np.prod(gains),
                total_duql_for_comparison=CONFIG["long_distance_QCD"]*np.prod(gains_other))


def constraint(tag,transfer,lower,finite):
    mp=CONFIG_1725["proton_mass_GeV"]
    mpi=CONFIG_1725["neutral_pion_mass_GeV"]
    years=CONFIG_1725["electron_mode_partial_lifetime_lower_years"]
    width=DATA["hbar_GeV_s"]/(DATA["year_seconds"]*years)
    w=DATA["pion_W0_GeV2"]
    w0=(abs(w[0])-lower*sum(w[1:]))/math.sqrt(2)
    if finite:
        old=runpy.run_path(str(HERE/"../1725/independent_check.py"))
        c2_old=old["coefficients"]("pi_e",lower)
    else:
        pref=mp/(32*math.pi)*(1-(mpi/mp)**2)**2
        c2_old=math.sqrt(width/pref)/w0
    c_match_max=c2_old/(tag["min_amplitude"]*transfer)
    canon=4*math.pi/(REFERENCE["new_inverse_colour"]*REFERENCE["U"]**2)
    rate=(canon/c_match_max)**2
    cap=tag["min_amplitude"]/math.sqrt(rate)
    n=tag["norm_enlarged_d"]
    rotate=math.asin(tag["min_amplitude"]/n)-math.asin(min(1,cap/n))
    return dict(W0_abs_GeV2=w0,low_tagged_coefficient_upper_GeV_minus2=c2_old,
                U_flavour_stripped_coefficient_max_GeV_minus2=c_match_max,
                predicted_rate_min_over_experimental_upper=rate,
                maximum_partial_lifetime_at_canonical_U_year=years/rate,
                maximum_relative_amplitude_at_U=1/math.sqrt(rate),
                fixed_coefficient_equivalent_scale_over_U=rate**.25,
                additional_same_helicity_cancellation_fraction_min=max(0.,1-1/math.sqrt(rate)),
                additional_electron_line_angle_min_deg=max(0.,math.degrees(rotate)))


def calculate():
    ep=flavour(CONFIG["flavour_step"])
    check=flavour(2*CONFIG["flavour_step"])
    assert abs(ep["min_amplitude"]-check["min_amplitude"])<1e-8
    avg=REFERENCE["three_equal_complete_weak_poles"]
    lo=CONFIG["weak_mass_floor_GeV"]
    num=CONFIG["weak_pair_count"]
    high=lo*(avg/lo)**num
    assert high<REFERENCE["U"]
    eq=transport([avg]*num)
    spl=transport([lo]*(num-1)+[high])
    test=transport([lo,avg,avg**2/lo])
    assert eq["total_qque"]<=test["total_qque"]<=spl["total_qque"]
    opts={}
    for lower,label in [(0,"central_matrix_element"),(1,"each_lattice_error_lower")]:
        opts[label]={str(k):constraint(ep,a,lower,finite) for k,a,finite in [
            ("unit_transfer_for_diagnosis",1.,False),("full_massless",eq["total_qque"],False),
            ("full_finite_electron_bound",eq["total_qque"],True)]}
    assert opts["each_lattice_error_lower"]["full_finite_electron_bound"][
        "predicted_rate_min_over_experimental_upper"]>3.30
    return dict(case=CONFIG["case"],U_GeV=REFERENCE["U"],inverse_alpha_U=REFERENCE["new_inverse_colour"],
                canonical_coefficient_GeV_minus2=4*math.pi/(REFERENCE["U"]**2*REFERENCE["new_inverse_colour"]),
                charged_family=ep,qque_gamma_and_duql_check=[a.tolist() for a in anomalous()],
                degenerate_weak=eq,maximally_split_at_fixed_product=spl,
                constraints=opts)


def tidy(o):
    if isinstance(o, dict):return {k:tidy(v) for k,v in o.items()}
    if isinstance(o, list):return [tidy(v) for v in o]
    if isinstance(o, float):return float(f"{o:.10g}")
    return o


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    dest=HERE/"results.json"
    out=tidy(calculate())
    if args.check:
        assert out==json.loads(dest.read_bytes())
        print("1742 OK: two gauge regimes, three weak thresholds, charged flavour and e+ pion bound.")
    else:
        dest.write_text(json.dumps(out,indent=2)+"\n",encoding="utf-8")

