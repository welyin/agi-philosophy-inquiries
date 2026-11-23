"""Independent neutral-kaon CP test of unitary charged-sector witnesses."""
import argparse
import cmath
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
COEFFICIENTS = np.array([32900., 20.55])
ERRORS = np.array([1500., .81])
CORRELATION = .34
COV = np.outer(ERRORS, ERRORS) * [[1., CORRELATION], [CORRELATION, 1.]]
EPS_EXP, EPS_ERROR = .002228, .000011


def read(path):
    return json.loads((HERE / path).read_bytes())


def amplitudes(moduli, gamma_deg):
    us, ub, cb = (moduli[k] for k in ("Vus", "Vub", "Vcb"))
    c13 = math.sqrt(1 - ub**2)
    s12, s23 = us / c13, cb / c13
    c12, c23 = math.sqrt(1 - s12**2), math.sqrt(1 - s23**2)
    gamma = math.radians(gamma_deg)
    delta = gamma + math.asin(c12*s23*ub*math.sin(gamma)/(s12*c23))
    phase = cmath.exp(1j*delta)
    ud, cd = c12*c13, -s12*c23-c12*s23*ub*phase
    td, ts = s12*s23-c12*c23*ub*phase, -c12*s23-s12*c23*ub*phase
    tb = c23*c13
    assert abs(cmath.phase(-ud*ub*phase/(cd*cb))-gamma) < 2e-15
    lam_u, lam_t = us*ud, ts.conjugate()*td
    basis = np.array([-(lam_t**2).imag / 2, (lam_u*lam_t).imag])
    assert min(basis)>0
    beta = cmath.phase(-cd*cb/(td*tb))
    approximate = cb*ub*math.sin(gamma)*np.array([
        cb**2*s12*(1-s12**2/2) - cb*ub*(1-2*s12**2)*math.cos(gamma)-ub**2*s12,
        s12*(1-s12**2/2)])
    return basis, approximate, beta, np.abs(np.array([td,ts])*tb)


def loop(mt, mc, mw):
    x, c = (mt/mw)**2, (mc/mw)**2
    t = (4*x-11*x*x+x**3)/(4*(1-x)**2) - 1.5*x**3*math.log(x)/(1-x)**3
    u = c*(1-math.log(x/c) + .75*x/(1-x) + .75*x*x*math.log(x)/(1-x)**2)
    return np.array([t, -u])


def new_central_factors():
    # Change only the specified inputs relative to Altmannshofer-Lewis.
    mass_shift = loop(162.92, 1.27, 80.379) / loop(163.48, 1.27, 80.379)
    return 1.01*(.7627/.7625)*mass_shift*np.array([1.,1.005])


def test(moduli, gamma_deg):
    v, approx, beta, products = amplitudes(moduli, gamma_deg)
    pred = float(v@COEFFICIENTS)
    sig = math.sqrt(v@COV@v)
    approx_shift = float(approx@COEFFICIENTS/pred-1.)
    assert abs(approx_shift)<.005
    corrected = float(v@(COEFFICIENTS*new_central_factors()))
    # Distance in the joint Gaussian metric of (a,b,experimental epsilon).
    residual = EPS_EXP-pred
    denom = sig**2 + EPS_ERROR**2
    coeff_move = COV@v*residual/denom
    exp_move = -EPS_ERROR**2*residual/denom
    assert abs(v@(COEFFICIENTS+coeff_move)-(EPS_EXP+exp_move))<1e-16
    mixing_coeff = np.array([6770.,10180.])
    clock_cov = np.outer([260.,370.],[260.,370.])*[[1.,.76],[.76,1.]]
    prob = products**2
    clocks = read("../1736/independent_results.json")["data"]
    observed = np.array([clocks["delta_md"],clocks["delta_ms"]])
    clock_residual = prob*mixing_coeff-observed[:,0]
    covariance = np.outer(prob,prob)*clock_cov+np.diag(observed[:,1]**2)
    clock_distance = math.sqrt(clock_residual@np.linalg.solve(covariance,clock_residual))
    radiative = products[1]*.3624
    return dict(moduli=moduli, gamma_deg=gamma_deg, sin2beta=math.sin(2*beta),
                absolute_B_clocks_ps_inv=(prob*mixing_coeff).tolist(),
                B_clocks_joint_distance=clock_distance,
                radiative_C7_times_lambda=radiative,
                radiative_marginal_pull=(radiative-.01477)/math.hypot(products[1]*.0151,.00078),
                coefficient_weights=v.tolist(), tt_ut_parts=(v*COEFFICIENTS).tolist(),
                epsilon=pred, correlated_error_at_fixed_CKM=sig,
                signed_pull_at_fixed_CKM=(pred-EPS_EXP)/math.sqrt(denom),
                relative_eq14_difference=approx_shift, alternate_central_epsilon=corrected,
                nearest_coefficient_shifts=coeff_move.tolist(),
                nearest_experimental_shift=exp_move)


def light_ratio(down):
    # FNAL/MILC Table IV; common arbitrary scale and any movement of b,u are profiled.
    s,d = .08487,.004291
    es,ed,rho = .00018,.000039,.43609054
    h = 2.
    roots = np.roots([s*s-h*h*es*es,-2*s*d+2*h*h*rho*es*ed,d*d-h*h*ed*ed])
    lower,upper = sorted(roots.tolist())
    ratio = down[0]/down[1]
    pull = (ratio*s-d)/math.sqrt(ed**2+ratio**2*es**2-2*ratio*rho*ed*es)
    return dict(high_ds=ratio,conditional_pull_Kds_one=pull,
                low_ds_interval_at_radius_two=[lower,upper],
                allowed_Kds=[ratio/upper,ratio/lower])


def intersect_two_CP_inputs():
    """Solve the two CP centres on the charged family and propagate its neutral sector."""
    import runpy
    family = runpy.run_path(str(HERE / "../1736/check.py"))
    ref = read("../1736/results.json")["witnesses"]["old_ub_cb_and_central_up"]
    us_shift = ref["weak_offsets"]["us"]
    target_B, B_error = .710, .011

    def evaluate(x):
        weak = dict(us=us_shift, ub=x[0], cb=x[1])
        row = family["shoot"](weak, 0., .25)
        result = test(row["moduli_MZ"],row["gamma"])
        residual = np.array([(result["epsilon"]-EPS_EXP)/EPS_ERROR,
                             (result["sin2beta"]-target_B)/B_error])
        return residual, row, result

    offsets = np.array([-.15,1.8])
    for _ in range(6):
        residual, row, result = evaluate(offsets)
        if np.linalg.norm(residual, np.inf) < 1e-7:
            break
        step = .002
        jac = np.column_stack([(evaluate(offsets+d)[0]-residual)/step for d in step*np.eye(2)])
        offsets -= np.linalg.solve(jac,residual)
    assert max(abs(residual))<1e-7
    raw = family["HAD"]["observables"](row["up_MZ"],row["down_MZ"])
    budget = family["hadron_budgets"](raw,2.)["both_two"]
    return dict(weak_offsets=row["weak_offsets"],up_MZ=row["up_MZ"],down_MZ=row["down_MZ"].tolist(),
                gamma=row["gamma"],CP_centre_result=result,largest_1736_pull=row["largest_width"],
                independent_light_ratio=light_ratio(row["down_MZ"]),
                signed_widths=row["signed_widths"],solar_eV2=row["solar_eV2"],
                neutral_FC68_90=row["observables"]["SK_T2K_FC68_90"],
                hadrons_at_old_centres=row["hadrons_at_central_matching"],
                minimum_old_relative_matching_at_two=budget["min_budget"],
                old_two_width_matching_example=budget["K_over_K_reference"].tolist())


def check_same_charge_profile():
    """Independent SVD and linear-contrast version of the Table IV ellipsoid distance."""
    old = read("../1736/results.json")
    kernel = read("../1734/results.json")["MZ_over_QCD_ratio_ds_ud_bs_reference"]
    primary = read("results.json")["independent_profiled_points"]
    source = read("inputs.json")["table_IV"]
    mu = np.array(source["means"][:3])
    sig = np.array(source["sigmas"][:3])
    cov = np.outer(sig,sig)*np.array(source["correlation"])[:3,:3]
    ck = read("../1702/inputs.json")["direct_semileptonic_moduli"]
    rows = {}
    for name, row in primary.items():
        moduli = {k: val[0]+row["weak"][k[1:]]*math.hypot(*val[1:]) for k,val in ck.items()}
        us,ub,cb = (moduli[k] for k in ("Vus","Vub","Vcb"))
        c13 = math.sqrt(1-ub*ub)
        s12,s23 = us/c13,cb/c13
        c12,c23 = math.sqrt(1-s12*s12),math.sqrt(1-s23*s23)
        gamma = math.radians(row["gamma_deg"])
        z = cmath.exp(1j*(gamma+math.asin(c12*s23*ub*math.sin(gamma)/(s12*c23))))
        top = np.array([s12*s23-c12*c23*ub*z,-c12*s23-s12*c23*ub*z,c23*c13])
        matrix = (np.eye(3)+math.expm1(row["profile_t"])*np.outer(top.conj(),top))@np.diag(
            [row["D_mz"]/kernel[0],1.,row["B_mz"]/kernel[2]])
        b,s,d = np.linalg.svd(matrix,compute_uv=False)
        contrast = np.array([[1.,-b/s,0.],[0.,-d/s,1.]])
        residual = contrast@mu
        rho = math.sqrt(residual@np.linalg.solve(contrast@cov@contrast.T,residual))
        assert abs(rho-row["profile_witness_radius"])<2e-7
        rows[name] = dict(same_charge_profile_radius=rho,CP_and_heldouts=test(moduli,row["gamma_deg"]))
        if name in old["witnesses"]:
            assert max(abs(moduli[k]-old["witnesses"][name]["moduli_MZ"][k]) for k in moduli)<1e-10
    return rows


def transfer_certificate():
    """Uniform covariance half-space with an interval for QCD-to-SM matching."""
    # These weak projections are analytic interval enclosures from R1736.
    projections = read("../1736/results.json")["two_width_Vus_projection"]
    lattice = read("inputs.json")["table_IV"]
    mu = np.array(lattice["means"][:3])
    sig = np.array(lattice["sigmas"][:3])
    cov = np.outer(sig,sig)*np.array(lattice["correlation"])[:3,:3]
    sem = read("../1702/inputs.json")["direct_semileptonic_moduli"]
    ub_max,cb_max = (sem[k][0]+2*math.hypot(*sem[k][1:]) for k in ("Vub","Vcb"))
    pb,eta = ub_max**2,ub_max**2+cb_max**2
    y = np.array(read("../1386/inputs.json")["yukawa"]["down"])
    # B0=b/s divided by Kbs, D0=d/s divided by Kds.
    klo,khi = 1.00195,1.00307
    Dcap,Bmin,Bmax = .056,51.,56.
    assert (y[2,0]-2*y[2,1])/(y[1,0]+2*y[1,1])/khi>Bmin
    assert (y[2,0]+2*y[2,1])/(y[1,0]-2*y[1,1])/klo<Bmax
    assert (y[0,0]+2*y[0,1])/(y[1,0]-2*y[1,1])/.99 < Dcap
    radius95 = math.sqrt(-2*math.log(.05))
    def positive_on_ellipsoid(a):
        a=np.array(a)
        return a@mu - radius95*math.sqrt(a@cov@a)>0
    assert all(positive_on_ellipsoid(a) for a in ((1,-53.5,0),(-1,54.5,0),(0,.055,-1)))
    tdmin = math.prod(sem[k][0]-2*math.hypot(*sem[k][1:]) for k in ("Vus","Vcb"))-ub_max
    def small_tail(x):
        return (Bmax*(eta+x)+(1+Dcap**2)/Bmin)/math.sqrt(x+(1-x)*tdmin**2)
    assert max(small_tail(0.),small_tail(.25))<53.5 and 2*Bmin*(1-eta)>54.5
    def leak(limit):
        gap=math.exp(-2*limit)*(1-eta)
        perturb=math.exp(2*limit)/Bmin**2
        return (math.sqrt(eta/gap)+perturb/(gap-perturb))**2
    tmax=math.log(2.)
    blog=max(abs(math.log(53.5/Bmax)),abs(math.log(54.5/Bmin)))
    for _ in range(2):
        tmax=blog/(1-2*leak(tmax))
    deltaD=.055*math.expm1(leak(tmax)*tmax)
    deltaB=max(Bmax-53.5,54.5-Bmin)
    # Rounded outwards compared to the bound on the full registered one-loop flow.
    flow_bound=.000081623
    assert read("results.json")["flow_bound"]["upper_F_mz_if_zero_at_upper"]<flow_bound
    direct = read("../1737/results.json")["thresholds"]
    assert abs(direct["top_minus_light_one_loop_in_decoupling"])*eta < .000006
    assert direct["HZ_bs_print_rounding_range"][0] - 3*eta*.00257 > klo
    assert direct["HZ_bs_print_rounding_range"][1] + 3*eta*.00257 < khi
    out={}
    for label,box in projections.items():
        if box is None:
            continue
        # Round the stored lower weak bound downwards before squaring.
        ps=(box[0]-1e-10)**2
        pd=1-ps-pb
        for radius in (2.,radius95):
            def support_value(Kds):
                drift=pd*Kds*deltaD+pb*khi*deltaB+flow_bound
                a=np.array([pb*khi,drift-ps,pd*Kds])
                return a@mu + radius*math.sqrt(a@cov@a)
            lo,hi=.99,1.05
            assert support_value(lo)<0<support_value(hi)
            for _ in range(60):
                mid=(lo+hi)/2.
                if support_value(mid)<0: lo=mid
                else: hi=mid
            out[f"{label}_R{radius:.5f}"]=dict(necessary_Kds_at_least=lo-1e-12,
                support_at_max_registered_top_only_matching=support_value(1.000006))
    return dict(Kbs_window=[klo,khi], Kds_window_considered=[.99,1.05],
                low_b_s_bounds=[53.5,54.5],abs_t_bound=tmax,
                delta_D_upper=deltaD, delta_B_upper=deltaB, lower_bounds=out)


def main():
    rows = {key:test(row["moduli_MZ"],row["gamma"])
            for key,row in read("../1736/results.json")["witnesses"].items()}
    return {"Altmannshofer_Lewis_2022_eq10_15":{"coefficients":COEFFICIENTS.tolist(),
            "standard_errors":ERRORS.tolist(), "correlation":CORRELATION},
            "experimental_abs_epsilon": [EPS_EXP, EPS_ERROR],
            "selected_2022_2025_updates_relative_to_2022":new_central_factors().tolist(),
            "endpoints":rows, "intersection_of_B_and_K_CP_centres":intersect_two_CP_inputs(),
            "main_profile_crosscheck": check_same_charge_profile(),
            "direct_threshold_windows":transfer_certificate()}


if __name__ == "__main__":
    cli = argparse.ArgumentParser()
    cli.add_argument("--check",action="store_true")
    args = cli.parse_args()
    result = main()
    path = HERE / "independent_results.json"
    if args.check:
        assert result == json.loads(path.read_bytes())
        print("1738 independent correlated neutral-kaon CP predictions verified")
    else:
        path.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf-8")
