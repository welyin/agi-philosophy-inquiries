"""Invert d/s and the reactor angle on a single conditional, historically selected matter branch."""
import argparse
import functools
import itertools
import json
import math
from pathlib import Path
import runpy
from types import SimpleNamespace, FunctionType

import numpy as np


HERE=Path(__file__).resolve().parent
def read(path):return json.loads((HERE/path).read_bytes())
CFG=read("inputs.json")
JOINT=runpy.run_path(str(HERE/CFG["pipeline"]))
OLD=read(CFG["qcdsf_and_inclusive"])
R,SIG_R=OLD["QCDSF_s_d"]
T=read(CFG["reactor"])["DayaBay_sin2_2theta13"]
CB,SIG_CB=OLD["inclusive_cb"]
DIRECT=read(CFG["exclusive"])["direct_semileptonic_moduli"]
CP=read(CFG["CP"])
FREQ=read(CFG["angles_and_freq"])
HIST=sorted(r["raw_real_b_coefficient"] for r in read(CFG["relative_gain_history"])["original_basis_verification"]["direct_forward_both_roots"])
SCALE_X=np.array([.005,.0002])
SCALE_F=np.array([SIG_R/R**2,T[1]])
DEFAULT_SHAPE=np.array([.2242,.00374])


@functools.lru_cache(maxsize=64)
def point(s,u,c,step=.25,gain=1.):
    # Run the exact archived point, adding a recorder at the existing daughter call.
    # The only intervention in sensitivity tests is the relative Majorana source gain.
    saved={}
    def loader(path):
        mod=runpy.run_path(str(path))
        if Path(path).parent.name=="1748":
            original=mod["daughter"]
            def with_gain(bg,fil):
                if not np.allclose(bg["nu"]["FROZEN"],HIST,rtol=1e-12,atol=1e-16):
                    raise ValueError("Unexpected neutral source boundary")
                new=dict(bg)
                new["nu"]=dict(bg["nu"],FROZEN=[HIST[0]*gain,HIST[1]])
                a,row=original(new,fil)
                saved.update(S_Bd=row["S_Bd"],epsilonK_pull=row["epsilon_K_pull"],
                             angles_in_box=row["neutrinos"]["angles_in_box"])
                return a,row
            mod["daughter"]=with_gain
        return mod
    f=JOINT["point"]
    patched=FunctionType(f.__code__,dict(f.__globals__,runpy=SimpleNamespace(run_path=loader)),
                         f.__name__,f.__defaults__)
    row=patched(HERE,s,u,c,step)
    row.update(saved)
    return row


def outputs(row):return np.array([row["down_over_s"],row["angle_readings"][0]])


def jac(x,c,step):
    eps=np.array([1.e-5,1.e-6])
    directions=np.eye(2)*eps
    return np.column_stack([(outputs(point(*(x+dx),c,step))-outputs(point(*(x-dx),c,step)))/(2*e)
                            for dx,e in zip(directions,eps)])


def solve(r,react,c,step,guess,J,gain=1.):
    x=np.array(guess,dtype=float)
    target=np.array([1/r,react])
    # Broyden in dimensionless coordinates. Recreate the physical material on every evaluation.
    A=J*SCALE_X[None,:]/SCALE_F[:,None]
    row=point(*x,c,step,gain)
    residual=(outputs(row)-target)/SCALE_F
    for _ in range(10):
        if max(abs(residual))<4e-7:break
        dx=np.linalg.solve(A,-residual)
        xnew=x+dx*SCALE_X
        assert .20<xnew[0]<.25 and .0027<xnew[1]<.0055
        row=point(*xnew,c,step,gain)
        new=(outputs(row)-target)/SCALE_F
        A+=np.outer(new-residual-A@dx,dx)/(dx@dx)
        x,residual=xnew,new
    assert max(abs(residual))<4e-7
    return x,row


def measured(row):
    us,ub,cb=[row["low_moduli"][k] for k in ("Vus","Vub","Vcb")]
    ud=math.sqrt(1-us**2-ub**2)
    def pull(x,data):return (x-data[0])/math.hypot(*data[1:])
    gp=CP["time_integrated_gamma_deg"]
    def fp(a,meaning):
        m,low,high=FREQ[meaning]
        return (a-m)/(low if a<m else high)
    return dict(us=us,ub=ub,cb=cb,z=row["z"],e_mu_tau_common_gain=row["tau_calibration"],
                calibrator_rebuilt_s_over_d=1/row["down_over_s"],angles=row["angle_readings"],
                angle_boxes=row["angles_in_box"],gamma_deg=row["gamma_deg"],
                time_integrated_gamma_in_68=gp["68"][0]<=row["gamma_deg"]<=gp["68"][1],
                S_Bd=row["S_Bd"],S_Bd_pull=pull(row["S_Bd"],CP["S_d"]),
                b_s=row["bottom_over_s"],epsilonK_pull=row["epsilonK_pull"],
                solar_calibrated_dm31=row["sun_calibrated_dm31_eV2"],solar_calibrated_dm32=row["sun_calibrated_dm32_eV2"],
                atmo_clock_pulls=[fp(row["sun_calibrated_dm31_eV2"],"juno_NO_dm31"),
                                 fp(row["sun_calibrated_dm32_eV2"],"accelerator_NO_dm32")],
                leptonic_delta_deg=row["delta_lepton_deg"],joint_FC68_90=row["SK_T2K_FC68_90"],
                mbb_meV=row["effective_double_beta_meV"],new_data_comparisons=dict(
                    exclusive_ub_pull=pull(ub,DIRECT["Vub"]),inclusive_ub_pull=pull(ub,OLD["inclusive_ub"]),
                    kaon_vector_us_pull=pull(us,OLD["holdout_Kl3"]),
                    kaon_axial_pull=pull(us/ud,OLD["holdout_axial"]),
                    nuclear_ud_pull=pull(ud,OLD["holdout_nuclear"])))


def scanning():
    # Locate changes of sign while re-solving d/s for each trial; do not clip to measured ub.
    records=[]
    for trial in CFG["scan_ub"]:
        row=JOINT["shoot"](HERE,1/R,trial,CB,.4)
        records.append([trial,row["trial_us"],row["angle_readings"][0]])
    assert all(x[2]<y[2] for x,y in zip(records,records[1:]))
    brackets=[[x[0],y[0]] for x,y in zip(records,records[1:]) if x[2]<=T[0]<=y[2]]
    assert brackets==[[.0036,.0038]]
    return dict(points=records,observed_sign_change_brackets=brackets)


def calculate():
    origin=DEFAULT_SHAPE
    J0=jac(origin,CB,.25)
    x,row=solve(R,T[0],CB,.25,origin,J0)
    J=jac(x,CB,.25)
    ec=2e-5
    dFc=(outputs(point(*x,CB+ec,.25))-outputs(point(*x,CB-ec,.25)))/(2*ec)
    data=np.column_stack((np.array([-SIG_R/R**2,0]),np.array([0,T[1]]),-dFc*SIG_CB))
    loads=np.linalg.solve(J,data)
    cov=loads@loads.T
    cg=cov/np.sqrt(np.outer(np.diag(cov),np.diag(cov)))
    corners={}
    for sr,st,sc in itertools.product((-1,1),repeat=3):
        offsets=np.array([sr,st,sc])
        _,candidate=solve(R+sr*SIG_R,T[0]+st*T[1],CB+sc*SIG_CB,.25,x+loads@offsets,J)
        corners[" ".join(map(str,offsets))]=measured(candidate)
    derivatives=J*SCALE_X[None,:]/SCALE_F[:,None]
    dJ=jac_at_smaller(x,CB)
    rho={}
    for d in CFG["gain_variations"]:
        _,candidate=solve(R,T[0],CB,.25,x,J,1+d)
        rho[str(d)]=measured(candidate)
    fine,yfine=solve(R,T[0],CB,.125,x,J)
    _,reactor_upper=solve(R,T[0]+T[1],CB,.25,x+loads[:,1],J)
    uvals=[p["ub"] for p in corners.values()]
    svals=[p["us"] for p in corners.values()]
    return dict(historical_raw_gain=HIST,calibrator_center=dict(strong_s_d=R,reactor=T[0],external_cb=CB),
                central=measured(row),central_trial_coordinates=x,finite_differences=dict(
                    J_rows_D_reactor_columns_us_ub=J,normalized_J=derivatives,normalized_determinant=np.linalg.det(derivatives),
                    normalized_condition_2=np.linalg.cond(derivatives),
                    maximum_relative_change_on_halving_probe=np.max(abs((dJ-J)/J)),
                    reactor_slope_along_fixed_ds=J[1,1]-J[1,0]*J[0,1]/J[0,0]),
                local_error=dict(input_order=["s/d","reactor","cb"],one_input_width_response=loads,
                                 linear_marginal_std=np.sqrt(np.diag(cov)),correlation=cg[0,1]),
                eight_one_width_vertices=corners,attained_us_range_vertices=[min(svals),max(svals)],
                attained_ub_range_vertices=[min(uvals),max(uvals)],relative_gain_offsets=rho,
                step_halving_trial_us_ub_delta=fine-x,central_fine=measured(yfine),
                only_reactor_at_upper_one_width=measured(reactor_upper),scan=scanning())


def jac_at_smaller(x,c):
    eps=[5e-6,5e-7]
    cols=[]
    for ax,h in enumerate(eps):
        shift=np.zeros(2)
        shift[ax]=h
        cols.append((outputs(point(*(x+shift),c,.25))-outputs(point(*(x-shift),c,.25)))/(2*h))
    return np.column_stack(cols)


def simple(x):
    if isinstance(x,dict):return {k:simple(v) for k,v in x.items()}
    if isinstance(x,(tuple,list,np.ndarray)):return [simple(v) for v in x]
    if isinstance(x,(bool,np.bool_)):return bool(x)
    if isinstance(x,(float,np.floating)):return float(f"{x:.9g}") if abs(x)>1e-15 else 0.
    return x


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    cli.add_argument("--check",action="store_true")
    args=cli.parse_args()
    result=simple(calculate())
    out=HERE/"results.json"
    if args.check:
        assert result==json.loads(out.read_bytes())
        print("1757 OK: two measured inputs inverted, signed branch, source sensitivity, eight vertices and held-out data.")
    else:
        out.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf8")

