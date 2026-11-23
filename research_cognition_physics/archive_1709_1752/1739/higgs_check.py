"""Propagate the parent and physical daughter Higgs quartic through the same pole."""
import argparse
import json
import math
from pathlib import Path
import runpy

import numpy as np


HERE = Path(__file__).resolve().parent
MAIN = runpy.run_path(str(HERE / "check.py"))
COLOR = 3


def advance(y, quartic, dt, grid, flow):
    """Integrate all charged and active-neutral Yukawas with one or more quartics."""
    def beta(z, lam, gauge):
        squares = z @ z.conj().transpose(0, 2, 1)
        weights = np.array([3, 3, 1, 1])
        trace = weights @ np.trace(squares, axis1=1, axis2=2).real
        trace4 = weights @ np.trace(squares @ squares, axis1=1, axis2=2).real
        gy2, gw2 = .6 * gauge[0], gauge[1]
        dl = (24 * lam**2 + (4*trace - 9*gw2 - 3*gy2)*lam - 2*trace4
              + 3/8*(2*gw2**2 + (gy2+gw2)**2))
        return flow["beta"](z, gauge, False), dl

    for i in range(0, len(grid)-1, 2):
        k1, a1 = beta(y, quartic, grid[i])
        k2, a2 = beta(y + dt*k1/2, quartic + dt*a1/2, grid[i+1])
        k3, a3 = beta(y + dt*k2/2, quartic + dt*a2/2, grid[i+1])
        k4, a4 = beta(y + dt*k3, quartic + dt*a3, grid[i+2])
        y = y + dt*(k1+2*k2+2*k3+k4)/6
        quartic = quartic + dt*(a1+2*a2+2*a3+a4)/6
    return y, quartic


def upper_running(y, quartic, segments, poles, columns, flow):
    y = y.copy()
    for start, susy, dt, grid in segments:
        assert not susy
        for col, pole in enumerate(poles):
            if start == pole:
                y[3,:,col] = columns[:,col]
        y, quartic = advance(y, quartic, dt, grid, flow)
    return y, quartic


def cw_terms(mass, soft, yukawas, fractions):
    r = 1-fractions**2
    fermion_slope = COLOR*mass**2/(8*math.pi**2) * r*yukawas**2
    fermion_quartic = COLOR/(8*math.pi**2) * r*(1-r)*yukawas**4
    # Heavy squared mass excess, normalized by M²; this form avoids subtraction.
    z = 1e-5
    shift = 2*r*z/(np.sqrt((1-z)**2+4*r*z)+1-z)
    np.testing.assert_allclose((shift-r*z)/z**2,r*(1-r),rtol=2e-5)
    k = (soft/mass)**2
    soft_ratio = math.log1p(k)+k*math.log1p(1/k)
    return dict(positive_potential_x_per_arm_GeV2=fermion_slope.tolist(),
                potential_x_at_pole_GeV2=sum(fermion_slope),
                finite_fermion_potential_quartic=sum(fermion_quartic),
                uniform_soft_complete_minus_aligned_x_GeV2=soft_ratio*sum(fermion_slope),
                uniform_soft_ratio_to_fermion_at_pole=soft_ratio)


def solve(name, step):
    cases = MAIN["CFG"]["nearby_measured_inputs_in_marginal_widths"]
    data = cases.get(name)
    args = () if data is None else tuple(data[k] for k in ("weak","up","g3","top","charm"))
    _, pt, _ = MAIN["charged_endpoint"](2, step, *args)
    f = pt["flow"]
    poles = f["fields"]["neutral_poles_GeV"]
    mass = poles[2]
    segments = [s for s in f["segments"](step) if s[0] < mass]
    bounds = f["physics"]["quartic"]
    central, error = np.array(bounds["lambda_paper"])/bounds["lambda_paper_factor"]
    y_src, parent = upper_running(MAIN["initial_state"](pt),central,segments,poles,f["physics"]["Y0"],f)
    up = f["basis"](y_src[0])[0][:,0]
    left, vals, vh = np.linalg.svd(y_src[1,:,:3])
    c = abs(up.conj()@left[:,-1])
    fb = math.sqrt((1+c)/2)
    factors = np.array([fb, 1., c])
    cw = cw_terms(mass,f["MS"],vals[[2,0]],factors[[2,0]])
    registered = MAIN["read"]("results.json")["rows"][name]
    assert abs(registered["c"]-c) < 2e-7
    np.testing.assert_allclose(cw["positive_potential_x_per_arm_GeV2"],
          [x["colored_Higgs_Yukawa_moment_GeV2"] for x in registered["finite_vev_d_and_b"]],rtol=1e-6)
    y = y_src.copy()
    y[1,:,:3] = y[1,:,:3] @ (vh.conj().T@np.diag(factors)@vh)
    daughter_at_threshold = y.copy()
    daughter_sources = np.zeros_like(f["physics"]["Y0"])
    lam = parent + np.array([0., cw["finite_fermion_potential_quartic"]])
    for start, _, dt, g in reversed(segments):
        y, lam = advance(y,lam,-dt,g[::-1],f)
        for col, pole in enumerate(poles):
            if pole == start:
                daughter_sources[:,col] = y[3,:,col]
                y[3,:,col] = 0
    # Solve the inverse question: which finite matching at M yields the published
    # MZ quartic or either edge? Carry the actual daughter source on the return leg.
    starting = np.array([central-error, central, central+error, lam[0]])
    replay_y, predicted = upper_running(y,starting,segments,poles,daughter_sources,f)
    assert np.max(abs(replay_y-daughter_at_threshold)) < 4e-9
    assert abs(predicted[-1]-parent)<2e-10
    correction = predicted[:3]-parent
    assert correction[0]<0<correction[2] and abs(correction[1])<1e-5
    return dict(M_GeV=mass,soft_GeV=f["MS"],c=c,fb=fb,
                lambda_MZ_input=central,lambda_MZ_marginal_width=error,
                parent_quartic_at_M=parent,delta_lambda_MZ_without_local_matching=lam[0]-central,
                delta_lambda_MZ_with_fermion_potential_matching=lam[1]-central,
                finite_matching_at_M_for_MZ_lower_center_upper=correction.tolist(),
                fraction_of_marginal_quartic_width=(lam[0]-central)/error,
                down_yukawas_at_M=vals[::-1].tolist(),**cw)


def main():
    names = ["2_d_and_b_bisector","absolute_and_oscillations"]
    answer = {}
    for name in names:
        normal = solve(name, .25)
        fine = solve(name, .125)
        shift = abs(normal["delta_lambda_MZ_without_local_matching"]-fine["delta_lambda_MZ_without_local_matching"])
        assert shift<2e-9
        fine["halfstep_absolute_change_in_delta_lambda_MZ"]=shift
        answer[name]=fine
    return answer


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    verifying = parser.parse_args().check
    answer = main()
    path = HERE/"higgs_results.json"
    if verifying:
        assert answer == json.loads(path.read_bytes())
        print("1739 physical daughter quartic, local matching window and colored curvature verified")
    else:
        path.write_text(json.dumps(answer,indent=2,allow_nan=False)+"\n",encoding="utf8")
