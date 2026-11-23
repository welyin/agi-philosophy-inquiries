"""Right-handed flavour resolution and one-loop transport on the registered two-port point."""
import argparse
import json
import math
from pathlib import Path
import runpy

import numpy as np

HERE = Path(__file__).resolve().parent
OLD = runpy.run_path(str(HERE.parent/"1739"/"check.py"))
TAG = "absolute_and_oscillations"
REG = json.loads((HERE.parent/"1739"/"results.json").read_bytes())["rows"][TAG]
LOOP = 16*math.pi**2


def right_port(y):
    left, sd, vh = np.linalg.svd(y[1, :, :3])
    left, sd, right = left[:, ::-1], sd[::-1], vh.conj().T[:, ::-1]
    up = np.linalg.svd(y[0, :, :3])[0][:, -1]
    c = abs(up.conj() @ left[:, 0])
    f = np.array([c, 1, math.sqrt((1+c)/2)])
    return sd, right, f


def sharp_cost(y, factors):
    q = 1-factors**2
    h = (y/y[-1])**2
    pairs = [(i, j) for i in range(3) for j in range(i)]
    value, i, j = max((abs(q[i]-q[j])/(h[i]+h[j]), i, j) for i,j in pairs)
    scalar = (q[i]*h[j]+q[j]*h[i])/(h[i]+h[j])
    diag = (q-scalar)/h
    assert np.isclose(max(abs(diag)), value)
    # For each pair the best uniform absolute loss error is the excess gap/2.
    def error_at_budget(b):
        return max(0., max((abs(q[i]-q[j])-b*(h[i]+h[j]))/2 for i,j in pairs))
    gains = np.sqrt(1/factors**2-1)
    pos = gains / np.sqrt(h)
    alternative_bottom = q[0]-q[2]
    return dict(positive_loss_dsb=q.tolist(), norm_with_unscaled_Y=value/y[-1]**2,
                exact_min_F_norm_with_Y_over_yb=value, scalar_at_min=scalar,
                sharp_F_diagonal=diag.tolist(), bottleneck_pair_dsb=[j,i],
                min_absolute_probability_error_norm_le_one=error_at_budget(1),
                min_absolute_probability_error_norm_le_4pi=error_at_budget(4*math.pi),
                converting_bisector_to_two_direct_with_single_DdaggerD={
                    "coefficient_for_normalized_Y":alternative_bottom,
                    "induced_d_s_loss_shifts":(alternative_bottom*h[:2]).tolist()},
                min_direct_heavy_map_norm=np.max(pos),
                min_direct_heavy_map_norm_with_unscaled_Y=np.max(pos)/y[-1])


def rates(y, x, flow, gauge, susy):
    dy = flow["beta"](y, gauge, susy)
    d = y[1, :, :3]
    # Flavour-universal factors in the heavy mass subblock cancel in X=M_h^{-1} M_ord.
    dx = (2 if susy else 1)*x@(d.conj().T@d)
    return dy, dx


def diagnostics(y, x, initial_right):
    sd, right, f = right_port(y)
    gram = x.conj().T@x
    loss = np.linalg.solve(np.eye(3)+gram, gram)
    physical = right.conj().T @ loss @ right
    target = np.diag(1-f*f)
    offdiag = physical - np.diag(np.diag(physical))
    res = physical-target
    probability = abs(initial_right.conj().T @ right)**2
    return dict(y_down_ascending=sd.tolist(), f_inst=f.tolist(),
                right_line_mixing_amplitudes_ds_db_sb=[math.sqrt(probability[i,j])
                                                       for i,j in [(0,1),(0,2),(1,2)]],
                diagonal_loss_drift_from_instantaneous=np.diag(res).real.tolist(),
                max_offdiagonal_loss=float(np.max(abs(offdiag))),
                distance_from_instantaneous_loss=float(np.linalg.norm(res, 2)))


def old_discrete_port(tagged_down, f):
    # Three inequivalent 5bar charges in the registered SU(5) basis. Diagonal
    # family rephasings used for the real texture leave all magnitudes unchanged.
    _, _, vh = np.linalg.svd(tagged_down)
    right = vh.conj().T[:, ::-1]
    target = right @ np.diag(1-f*f) @ right.conj().T
    off = target - np.diag(np.diag(target))
    charges = np.array([2,0,1])
    vac = json.loads((HERE.parent/"1685"/"inputs.json").read_bytes())["scalar_magnitude_over_cutoff"]
    # A2=S²; either the field or its adjoint can supply the Hermitian charge transfer.
    pairs = [(0,1),(0,2),(1,2)]
    return dict(abs_tagged_offdiagonal_12_13_23=[abs(off[i,j]) for i,j in pairs],
                any_charge_only_operator_error_norm_lower=np.max(abs(off)),
                best_diagonal_Frobenius_error=np.linalg.norm(off),
                K_diagonal_in_charge_basis=np.diag(target).real.tolist(),
                magnitudes_if_each_kinetic_tensor_has_own_coefficient=[
                    abs(off[i,j])/vac**abs(charges[i]-charges[j]) for i,j in pairs])


def solve(step):
    args = OLD["CFG"]["nearby_measured_inputs_in_marginal_widths"][TAG]
    summary, pt, _ = OLD["charged_endpoint"](2, step, args["weak"], args["up"],
                                              args["g3"], args["top"], args["charm"])
    y = OLD["initial_state"](pt)
    flow = pt["flow"]
    mass = flow["fields"]["neutral_poles_GeV"][2]
    x = np.zeros((2,3), complex)
    rec = {}
    for start, susy, dt, grid in flow["segments"](step):
        if start == mass:
            sd, basis, factors = right_port(y)
            angles = np.sqrt(1/factors**2 - 1)
            x = np.diag(angles[[0,2]]) @ basis[:,[0,2]].conj().T
            rec["sharp_spurion_cost"] = sharp_cost(sd, factors)
            rec["M_threshold"] = diagnostics(y, x, basis)
            # Derivative of the target spectral projections in the instantaneous mass frame.
            d, u = y[1, :, :3], y[0, :, :3]
            dHd_noncommuting = -3*basis.conj().T @ (d.conj().T@u@u.conj().T@d) @ basis
            lossvals = 1-factors**2
            speed = []
            for i,j in ((0,1),(0,2),(1,2)):
                speed.append(abs((lossvals[j]-lossvals[i])*dHd_noncommuting[i,j]/
                                 (sd[j]**2-sd[i]**2))/LOOP)
            rec["threshold_target_offdiagonal_speed_per_log_mu_ds_db_sb"] = speed
            assert min(speed)>0 and max(speed)>1e-8
            assert np.allclose(sd, REG["down_above"], rtol=1e-7, atol=0)
        if start == flow["MS"]:
            rec["soft_scale_before_matching"] = diagnostics(y, x, basis)
            t = summary["positive_soft_tan_beta"]
            cosb = 1/math.hypot(1,t)
            y /= (cosb*np.array([t,1,1,t]))[:,None,None]
        for j, pole in enumerate(flow["fields"]["neutral_poles_GeV"]):
            if start == pole:
                y[3,:,j] = flow["physics"]["Y0"][:,j]
        for k in range(0, len(grid)-1, 2):
            a,b = rates(y,x,flow,grid[k],susy)
            c,d = rates(y+dt*a/2,x+dt*b/2,flow,grid[k+1],susy)
            e,f = rates(y+dt*c/2,x+dt*d/2,flow,grid[k+1],susy)
            g,h = rates(y+dt*e,x+dt*f,flow,grid[k+2],susy)
            y += dt*(a+2*c+2*e+g)/6
            x += dt*(b+2*d+2*f+h)/6
    rec["10PeV"] = diagnostics(y,x,basis)
    rec["10PeV_old_unbroken_flavour_charge_only"] = old_discrete_port(pt["D"],right_port(y)[2])
    assert np.allclose(rec["10PeV"]["y_down_ascending"], REG["high_yukawas"], rtol=2e-7, atol=0)
    return rec


def tidy(obj):
    if isinstance(obj, dict):
        return {k:tidy(v) for k,v in obj.items()}
    if isinstance(obj, list):
        return [tidy(v) for v in obj]
    if isinstance(obj, float):
        # Remove machine-zero diagonalization noise at the exactly aligned initial condition.
        return 0 if abs(obj)<1e-14 else float(f"{obj:.8g}")
    return obj


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    coarse, fine=solve(.25),solve(.125)
    diff=abs(coarse["10PeV"]["distance_from_instantaneous_loss"]-
             fine["10PeV"]["distance_from_instantaneous_loss"])
    assert diff<2e-10
    out=tidy(dict(reference=TAG,main=fine,step_half_distance_difference=diff))
    dest=HERE/"independent_results.json"
    if args.check:
        assert out==json.loads(dest.read_bytes())
        print("1740 independent check OK: sharp spurion cost, three nonzero drift terms, 10 PeV transfer.")
    else:
        dest.write_text(json.dumps(out,indent=2)+"\n",encoding="utf-8")

