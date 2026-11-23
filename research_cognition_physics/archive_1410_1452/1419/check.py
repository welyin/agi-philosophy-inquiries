"""First-generation semileptonic vector coefficients in a real, CKM-aligned contact model."""
import argparse
import itertools
import json
from math import hypot, pi, sqrt
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
SIGNS = np.array([1., -1.])


def read(path):
    return json.loads((HERE / path).read_text(encoding="utf-8"))


def expand(wilson):
    singlet, triplet, ru, rd, eq, eu, ed = wilson
    ll = singlet-triplet*SIGNS
    lr, rr = np.array([ru, rd]), np.array([eu, ed])
    nu = singlet + triplet*SIGNS + lr
    return nu, ll+lr-eq-rr, ll-lr+eq-rr, 2*triplet


def signed_pv_box(pv, width):
    proton, cs = pv["proton"], pv["caesium"]
    ap = proton["pv_electron_scattering_only"][0]-proton["SM"][0]
    dp = hypot(proton["pv_electron_scattering_only"][1],proton["SM"][1])
    ac = cs["value_and_errors_experiment_atom"][0]-cs["SM"][0]
    dc = hypot(*cs["value_and_errors_experiment_atom"][1:],cs["SM"][1])
    charge_map = -2*np.array([[2,1],[2*cs["Z"]+cs["N"],cs["Z"]+2*cs["N"]]])
    corners = [np.linalg.solve(charge_map,[ap+i*width*dp,ac+j*width*dc])
               for i,j in itertools.product([-1,1],repeat=2)]
    return corners


def projection_floor(energies_squared, delta_cc, nu, v, l):
    """Amplitude floor at a single hard energy, or rigorously over the full hard interval."""
    low, high = energies_squared
    weak = l["hard_inputs"]
    w = weak["sin2_weak"]
    mz, gamma = weak["MZ_GeV"], weak["widthZ_GeV"]
    e2 = 4*pi/weak["alpha_em_inverse"]
    charge = np.array([2/3,-1/3])
    isospin = SIGNS/2
    zfactor = e2*(-.5+w)*(isospin-2*w*charge)/(w*(1-w))
    norm = 3*low/(96*pi)*l["normalization_pb_per_GeV_minus2"]
    if low == high and delta_cc[0] == delta_cc[1]:
        sm = -2*e2*charge/low + zfactor/(low-mz*mz+1j*mz*gamma)
        contact = -2*(nu-SIGNS*delta_cc[0])/v**2
        return norm*sum(abs(sm+contact)**2)
    assert low > mz*mz
    sm_ceiling = 2*e2*abs(charge)/low+abs(zfactor)/(low-mz*mz)
    ends = np.array([nu-SIGNS*d for d in delta_cc])
    contact_floor = np.where(ends[0]*ends[1] <= 0,0,np.min(abs(ends),axis=0))
    summed = np.maximum(2*contact_floor/v**2-sm_ceiling,0)
    return norm*sum(summed**2)


def compute():
    config = read("inputs.json")
    low, lep = config["low_energy"],config["LEP"]
    pv = read(low["PV"])
    pion = read(low["pion"])["PIENU"]
    beta = read(low["beta"])["ud_beta"]
    gf = read(low["GF"])["const_GeV"]["GF"]
    v = 1/sqrt(sqrt(2)*gf)
    width, angle = low["width"],low["sin2_low_illustrative_tree_1418"]
    nu_sm = SIGNS/2-2*angle*np.array([2/3,-1/3])
    flip = -2*nu_sm
    c1_corners = signed_pv_box(pv,width)
    iso = [sum(p) for p in c1_corners]
    exp,err,sys = pion["exp_in_units_1e_4"]
    pred,errsm = pion["same_paper_SM"]
    a,b = sqrt((exp-width*hypot(err,sys))/(pred+width*errsm)),sqrt((exp+width*hypot(err,sys))/(pred-width*errsm))
    cc_branches = [(a-1,b-1),(-1-b,-1-a)]
    beta_abs = [(beta[0]-width*beta[1])/b,(beta[0]+width*beta[1])/a]
    charged_ranges = []
    for branch in cc_branches:
        hold = [float(c1 @ SIGNS + 2*d - 2*(flip @ SIGNS)) for c1 in c1_corners for d in branch]
        floor = projection_floor([(lep["sqrt_s_GeV"]*lep["sqrt_sprime_over_s_cut"])**2,
                                  lep["sqrt_s_GeV"]**2],branch,flip,v,lep)
        e=lep["sigma_had_pb_and_errors"]
        charged_ranges.append({
            "delta_cc":list(branch), "axial_isovector_shift_if_flip":[min(hold),max(hold)],
            "hard_ud_pb_floor_whole_sprime_window":floor,
            "critical_positive_response_weight_at_zero_additive_error": (e[0]+width*hypot(*e[1:]))/floor,
        })
    # Real dimension-six gauge operators hiding a whole mirror from C1 and 2C2u-C2d.
    # Set all target shifts in the registered compressed partner data to zero.
    eq = 2*flip[0]-flip[1]
    coeffs = [0.,0.,*flip,eq,flip[0]-eq,flip[1]-eq]
    nu,c1,c2,cc = expand(coeffs)
    assert np.allclose(nu,flip) and np.allclose(c1,0) and abs(2*c2[0]-c2[1])<1e-13 and cc == 0
    assert nu_sm.sum()*-2 > max(iso)
    assert max(charged_ranges[1]["axial_isovector_shift_if_flip"])<0
    assert min(charged_ranges[0]["axial_isovector_shift_if_flip"])>0
    hard_points = [projection_floor([lep["sqrt_s_GeV"]**2]*2,[d,d],flip,v,lep) for d in [0,-2]]
    # Verify the flavour expansion and exactly the C1/C2 identity in the seven-operator space.
    for basis in np.eye(7):
        n,p1,p2,x=expand(basis)
        assert abs((n@SIGNS)-x-(p1-p2)@SIGNS/2)<1e-13
    mat=np.array([np.r_[expand(k)[1],expand(k)[3]] for k in np.eye(7)[:4]]).T
    assert np.linalg.matrix_rank(mat)==3 and np.allclose(mat@np.array([1.,0.,-1.,-1.]),0)
    return {
        "flip_nu_quark_shift":flip.tolist(), "flip_required_sum":sum(flip),
        "two_width_four_operator_allowed_sum":[min(iso),max(iso)],
        "minimal_extra_isoscalar_shift_for_whole_mirror":sum(flip)-max(iso),
        "two_width_absolute_Vud_from_beta_and_pion":beta_abs,
        "two_width_CC_sign_branches":charged_ranges,
        "seven_operator_mirror_witness":list(coeffs),"witness_C2_shifts":c2.tolist(),
        "PVDIS_shift_2u_minus_d":2*c2[0]-c2[1],
        "observed_DIS_ppm":low["PVDIS_observed_ppm"],"baseline_DIS_ppm":low["PVDIS_SM_asym_ppm"],
        "reference_at_200_GeV_ud_lower_pb_delta_cc_zero_or_minus_two":hard_points,
        "LEP_hadronic_upper_two_width_pb":e[0]+width*hypot(*e[1:]),
    }


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    group=parser.add_mutually_exclusive_group()
    group.add_argument("--write",action="store_true")
    group.add_argument("--check",action="store_true")
    args=parser.parse_args()
    report=json.dumps(compute(),indent=2,allow_nan=False)+"\n"
    output=HERE/"results.json"
    if args.write:
        with output.open("x",encoding="utf-8") as handle:
            handle.write(report)
    elif json.loads(output.read_text(encoding="utf-8"))!=json.loads(report):
        raise SystemExit("1419 results mismatch")
    print("1419 operator identity, two charge signs, PV mirror and positive LEP floor verified")
