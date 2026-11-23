"""Independent Majorana-plane, electron projector, beta and oscillation checks."""
import argparse
import json
from math import acos, asin, atan2, cos, pi, sin, sqrt
from pathlib import Path

import numpy as np


CENTRAL_A, HALF_A = .0851, .0024
INITIAL_COLUMN = np.array([2.,-1.,1.])/sqrt(6)
E_WEIGHT = INITIAL_COLUMN[0]**2


def small_s13(amp):
    # Rationalized solution avoids cancellation on the small root of 4x(1-x)=amp.
    return amp/(2*(1+sqrt(1-amp)))


def solar(x):
    return 1-E_WEIGHT/(1-x)


def ray_angle(p):
    return atan2(sqrt(1-p),sqrt(p))


def minimal_tilt(p_box):
    nearest = min(max(E_WEIGHT, p_box[0]), p_box[1])
    return abs(ray_angle(nearest)-ray_angle(E_WEIGHT))


def solar_run(s, half, xs):
    xmin,x,xmax = xs
    target = [(1-xmax)*(1-s-half),(1-xmin)*(1-s+half)]
    z = (solar(x)-s)/half
    curve_band = [solar(xmax),solar(xmin)]
    first_box_gap = max(0., curve_band[0]-s-half)
    return {"input_solar": [s,half], "central_gap_in_published_marginal_widths": z,
            "gap_between_one_width_rectangles": first_box_gap,
            "electron_probability_box":target,
            "central_minimum_plane_angle":minimal_tilt([(1-x)*(1-s)]*2),
            "full_box_minimum_plane_angle":minimal_tilt(target)}


def beta_readouts(x):
    m2,m3 = sqrt(.0000749),sqrt(.002513)
    p2 = 1/3-x
    return np.array([sqrt(p2*m2*m2+x*m3*m3), abs(p2*m2-x*m3), p2*m2+x*m3])*1000


def plane_reconstruction(x):
    # Build the Takagi matrix from physical mass eigenpairs. The last two columns
    # admit arbitrary mixing and Majorana phases on the common orthogonal plane.
    v2 = np.array([1.,1.,-1.])/sqrt(3)
    v3 = np.array([0.,1.,1.])/sqrt(2)
    theta = asin(sqrt(3*x))
    phase = .73
    columns = np.column_stack((INITIAL_COLUMN,cos(theta)*v2-sin(theta)*np.exp(-1j*phase)*v3,
                               np.exp(1j*.58)*(sin(theta)*np.exp(1j*phase)*v2+cos(theta)*v3)))
    m = np.array([0.,sqrt(.0000749),sqrt(.002513)])
    target = columns.conj() @ np.diag(m) @ columns.conj().T
    assert np.allclose(columns.T @ target @ columns, np.diag(m))
    assert np.isclose(abs(columns[0,2])**2,x)
    # Different complex n and an additional complex invertible column change.
    for n in [-2.,0.,1.+2.5j,4j]:
        sources = np.array([[0,1],[1,n],[1,n-2]],complex)
        for change in (np.eye(2),np.array([[1j,.35],[.1-.3j,2]])):
            cols = sources @ change
            dual = np.linalg.pinv(cols)
            symmetric_kernel = dual @ target @ dual.T
            assert np.allclose(symmetric_kernel, symmetric_kernel.T)
            assert abs(np.linalg.det(symmetric_kernel)) > 1e-7
            realised = cols @ symmetric_kernel @ cols.T
            assert np.allclose(realised,target) and np.allclose(realised@INITIAL_COLUMN,0)
    return {"mass_singular_values_eV":m.tolist(),
            "actual_beta_meV":1000*sqrt(np.sum(abs(columns[0])**2*m**2)),
            "one_Majorana_phase_example_mbb_meV":1000*abs(target[0,0])}


def compute():
    xs=[small_s13(v) for v in (CENTRAL_A-HALF_A,CENTRAL_A,CENTRAL_A+HALF_A)]
    xmin,x,xmax=xs
    first=solar_run(.3092,.0087,xs)
    update=solar_run(.3036,.0064,xs)
    io_low=minimal_tilt([xmin,xmax])
    io_high=minimal_tilt([1-xmax,1-xmin])
    assert io_low > io_high > 0
    all_beta = np.array([beta_readouts(v) for v in xs])
    reconstruction=plane_reconstruction(x)
    assert np.isclose(reconstruction["actual_beta_meV"],all_beta[1,0])
    assert all_beta[1,1]<reconstruction["one_Majorana_phase_example_mbb_meV"]<all_beta[1,2]
    assert 3*x < 1 and first["gap_between_one_width_rectangles"]>0
    return {
        "inputs": {"Daya_Bay":"https://arxiv.org/pdf/2211.14988",
                   "JUNO_59_days":"https://doi.org/10.1038/s41586-026-10538-z",
                   "JUNO_207_days":"https://arxiv.org/abs/2609.37895",
                   "mass_splittings":"../1370/inputs.json"},
        "small_x13_low_center_high":xs, "NO_Tm1_solar_low_center_high":[solar(v) for v in xs[::-1]],
        "first_JUNO_NO":first,"overlapping_update_JUNO_NO":update,
        "IO_predicted_oscillation_amplitude":4*E_WEIGHT*(1-E_WEIGHT),
        "IO_minimum_plane_tilt_radians_low_octant":io_low,
        "IO_minimum_plane_tilt_radians_both_octants":min(io_low,io_high),
        "minimum_IO_matrix_change_over_observed_m2nd_low_octant":sin(io_low),
        "minimum_IO_matrix_change_over_observed_m2nd_both_octants":sin(min(io_low,io_high)),
        "NO_beta_mbb_min_max_meV_central":all_beta[1].tolist(),
        "NO_beta_meV_full_reactor_box":[all_beta[0,0],all_beta[-1,0]],
        "NO_mbb_meV_full_reactor_box":[all_beta[-1,1],all_beta[-1,2]],
        "plane_reconstruction":reconstruction
    }


def rounded(obj):
    if isinstance(obj, dict):
        return {k: rounded(v) for k,v in obj.items()}
    if isinstance(obj, list):
        return [rounded(v) for v in obj]
    return round(obj,9) if isinstance(obj,float) else obj


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write",action="store_true")
    parser.add_argument("--check",action="store_true",help="check saved data (default)")
    args=parser.parse_args()
    file=Path(__file__).with_name("independent_results.json")
    answer=rounded(compute())
    if args.write:
        with file.open("x",encoding="utf-8") as out:
            json.dump(answer,out,indent=2,ensure_ascii=False)
            out.write("\n")
    else:
        assert answer==json.loads(file.read_text(encoding="utf-8"))
    print("1371 independent Takagi, projector, oscillation and joint beta checks passed.")
