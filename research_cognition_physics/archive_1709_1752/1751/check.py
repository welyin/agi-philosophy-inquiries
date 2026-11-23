"""Finite radial vacuum response and positive physical threshold matching on the 1748 point."""
import argparse
import json
import itertools
from pathlib import Path
import runpy

import numpy as np


HERE = Path(__file__).resolve().parent
CFG = json.loads((HERE/"inputs.json").read_bytes())
BASE = runpy.run_path(str(HERE/CFG["mother"]))
SPECTATOR = runpy.run_path(str(HERE/CFG["source_check"]))


def inverse_sqrt(H):
    spectrum, L = np.linalg.eigh(H)
    assert min(spectrum) > 0
    return (L/np.sqrt(spectrum)) @ L.conj().T


def source(bg):
    Y = bg["yW"][2,:,:3]
    H = Y @ Y.conj().T
    h,L = np.linalg.eigh(H)
    h /= h[-1]
    j = np.array([np.prod([h[i]-h[k] for k in range(3) if k!=i]) for i in range(3)])
    j /= max(abs(j))
    np.testing.assert_allclose(SPECTATOR["spectator_test"](H)["normalized_eigenvalues"],h)
    return L,j


def vacuum(eta,j):
    field = []
    for x in j:
        lo,hi = 1.,1+eta*abs(x)/2
        for _ in range(60):
            mid = (lo+hi)/2
            if mid*(mid-1)*(mid+1)<eta*abs(x):
                lo=mid
            else:
                hi=mid
        field.append(np.sign(x)*(lo+hi)/2)
    field=np.array(field)
    np.testing.assert_allclose(field*(field-1)*(field+1),eta*j,atol=2e-15)
    radial=3*field**2-1
    # Divided difference of the stationarity polynomial. Two real modes per pair.
    pair=[]
    for i,k in ((0,1),(0,2),(1,2)):
        curve=field[i]**2+field[i]*field[k]+field[k]**2-1
        np.testing.assert_allclose(curve,eta*(j[i]-j[k])/(field[i]-field[k]),rtol=2e-8,atol=1e-14)
        pair.append(curve)
    assert min(pair)>0
    return field,radial,np.array(pair)


def loop_comparison(j,eta=.001):
    roots=[np.sort(np.roots([1,0,-1,-eta*q]).real) for q in j]
    assert max(abs(eta*j))<2/np.sqrt(27)
    rows={}
    for x in itertools.product(*roots):
        x=np.array(x)
        radial=3*x*x-1
        pairs=[(i,k) for i in range(3) for k in range(i)]
        angular=np.array([x[i]**2+x[i]*x[k]+x[k]**2-1 for i,k in pairs])
        curvature=np.r_[radial,np.repeat(angular,2)]
        if min(curvature)<=0:
            continue
        base=np.sum((x*x-1)**2/4-eta*j*x)
        cw=np.sum(curvature**2*(np.log(curvature/2)-1.5))/(64*np.pi**2)
        def deriv(q):
            return q*(np.log(q/2)-1)/(32*np.pi**2)
        grad=deriv(radial)*6*x
        for q,(i,k) in zip(angular,pairs):
            grad[i]+=2*deriv(q)*(2*x[i]+x[k])
            grad[k]+=2*deriv(q)*(2*x[k]+x[i])
        label="".join("+" if t>0 else "-" for t in x)
        rows[label]=dict(tree_over_lambda_f4=base,scalar_one_loop_over_lambda_squared_f4=cw,
                         radial_shift_at_example_lambda=-CFG["scalar_lambda"]*grad/radial)
    assert len(rows)==4
    ref,comp=rows["+-+"],rows["+++"]
    gap=comp["tree_over_lambda_f4"]-ref["tree_over_lambda_f4"]
    slope=comp["scalar_one_loop_over_lambda_squared_f4"]-ref["scalar_one_loop_over_lambda_squared_f4"]
    lam=-gap/slope
    def winner(c):
        return min(rows,key=lambda k:rows[k]["tree_over_lambda_f4"]+c*rows[k]["scalar_one_loop_over_lambda_squared_f4"])
    assert 0<lam<.001 and winner(CFG["scalar_lambda"])=="+-+" and winner(1.)=="+++"
    return dict(sectors=rows,all_positive_minus_alternating_tree=gap,
                all_positive_minus_alternating_CW=slope,
                crossing_lambda_at_fixed_MSbar_boundary=lam,
                lowest_of_classical_minima_to_one_loop_at_example_lambda=winner(CFG["scalar_lambda"]),
                lowest_of_classical_minima_to_one_loop_at_lambda_one=winner(1.))


def observables(bg,F):
    gain,row=BASE["daughter"](bg,F)
    bounds=bg["model"]["JOINT"]["LEPTONS"]
    pulls=(np.array(row["MZ_leptons"])-bounds[:,0])/bounds[:,1]
    nu=row["neutrinos"]
    clocks=np.array([nu["dm31_eV2"],nu["dm32_eV2"]])
    cp=[]
    for k,x in zip(("juno_NO_dm31","accelerator_NO_dm32"),clocks):
        center,low,high=bg["nu"]["FREQ"][k]
        cp.append((x-center)/(low if x<center else high))
    return dict(tau_scale=gain,charged_pulls=pulls,angles=nu["angle_readings"],
                angle_passes=nu["angles_in_box"],solar_center_clocks=clocks,
                clock_pulls=cp,FC68_FC90=nu["SK_T2K_FC68_90"],delta_deg=nu["delta_deg"],
                quark_mass_radius=row["absolute_four_ratio_radius"],epsilonK_pull=row["epsilon_K_pull"],
                e_mu_pass=bool(max(abs(pulls[:2]))<=2),
                passed=bool(max(abs(pulls))<=2 and all(nu["angles_in_box"]) and
                            max(abs(np.array(cp)))<=2 and nu["SK_T2K_FC68_90"][1]))


def response(bg,X,L,j,eta,bias=0):
    field,radial,pairs=vacuum(eta,j+bias)
    G=(L*field)@L.conj().T
    Xt=X@G
    F=inverse_sqrt(np.eye(3)+Xt.conj().T@Xt)
    norm=inverse_sqrt(np.eye(2)+Xt@Xt.conj().T)
    norm0=inverse_sqrt(np.eye(2)+X@X.conj().T)
    M=bg["Mweak"]
    blocked=np.hstack([Xt,np.eye(2)])
    np.testing.assert_allclose(np.linalg.svd(norm@blocked,compute_uv=False),1,atol=2e-14)
    frozen_bare_poles=M*np.linalg.svd(norm0@blocked,compute_uv=False)[::-1]
    # At the vacuum the full bilinear is M norm [X Phi/f, I], with constant norm.
    scalar_yuk=M*np.linalg.norm(norm@X,2)/CFG["example_f_GeV"]
    # Two Dirac SU(2) doublets (eight real fermionic degrees for each singular pole).
    doublet_gradient=2*(M/CFG["example_f_GeV"])**4*np.linalg.norm(norm@X,2)**2*max(abs(field))
    doublet_gradient/=4*np.pi**2*CFG["scalar_lambda"]
    return dict(eta=eta,bias=bias,phi_over_f=field,radial_m2_over_lambda_f2=radial,
                pair_labels=["e_mu","e_tau","mu_tau"],pair_m2_over_lambda_f2=pairs,
                min_to_max_mass_ratio=np.sqrt(min(*radial,*pairs)/max(*radial,*pairs)),
                scalar_mass_min_GeV=CFG["example_f_GeV"]*np.sqrt(CFG["scalar_lambda"]*min(*radial,*pairs)),
                scalar_heavy_light_yukawa_opnorm=scalar_yuk,
                log10_vector_doublet_gradient_bound_in_units_lambda_f4_at_Q_equal_M=np.log10(doublet_gradient),
                fixed_bare_poles_GeV=frozen_bare_poles,
                fixed_bare_pole_shifts_GeV=frozen_bare_poles-M,
                required_bare_block_relative_norm_shift=np.linalg.norm(norm-norm0,2)/np.linalg.norm(norm0,2),
                contraction_spectrum=np.linalg.eigvalsh(F),low=observables(bg,F))


def calculate(step):
    bg=BASE["OLD"]["mother"](step)
    X,F0,_,_=BASE["terms"](bg,CFG["root_zero_based"],(1,1,1))
    L,j=source(bg)
    rows={str(eta):response(bg,X,L,j,eta) for eta in CFG["eta_values"]}
    assert rows["0.001"]["low"]["passed"]
    assert not rows["1"]["low"]["passed"]
    # With pre-existing charged gains recalibrated by the inverse radial vevs,
    # the original rows, heavy masses, and complete endpoint are recovered exactly.
    field,_,_=vacuum(.5,j)
    compensated=X@((L/abs(field))@L.conj().T)@((L*field)@L.conj().T)
    expected=BASE["terms"](bg,0,(1,-1,1))[0]
    np.testing.assert_allclose(compensated,expected,atol=3e-14)
    flip=response(bg,X,L,j,.001,-1.01*j[1])
    assert not flip["low"]["passed"] and all(flip["phi_over_f"]>0)
    loop=loop_comparison(j)
    x=vacuum(.001,j)[0]+loop["sectors"]["+-+"]["radial_shift_at_example_lambda"]
    Xshift=X@((L*x)@L.conj().T)
    flow_with_shift=observables(bg,inverse_sqrt(np.eye(3)+Xshift.conj().T@Xshift))
    assert flow_with_shift["passed"]
    loop["endpoint_with_first_order_singlet_vev"]=flow_with_shift
    return dict(charged_source=j,heavy_weak_pole_GeV=bg["Mweak"],points=rows,
                singlet_effective_potential_at_eta_0p001=loop,
                competing_universal_bias=flip,exact_charged_gain_compensation_at_eta_half=1/abs(field),
                compensated_mass_row_norm_error=np.linalg.norm(compensated-expected,2))


def tidy(x):
    if isinstance(x,dict):return {k:tidy(v) for k,v in x.items()}
    if isinstance(x,(tuple,list,np.ndarray)):return [tidy(v) for v in x]
    if isinstance(x,(bool,np.bool_)):return bool(x)
    if isinstance(x,(float,np.floating)):return float(f"{x:.10g}") if abs(x)>1e-13 else 0.
    return x


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=calculate(.125)
    other=calculate(.25)
    delta=abs(result["points"]["0.001"]["low"]["solar_center_clocks"][0]-
              other["points"]["0.001"]["low"]["solar_center_clocks"][0])
    assert delta<1e-9
    result["step_halving_dm31_shift_at_eta_0p001"]=delta
    result=tidy(result)
    output=HERE/"results.json"
    if args.check:
        assert result==json.loads(output.read_bytes())
        print("1751 OK: scalar vacuum, Hessian, pole matching and transported flavour data.")
    else:
        with output.open("w",encoding="utf8") as f:
            json.dump(result,f,indent=2,allow_nan=False)
            f.write("\n")

