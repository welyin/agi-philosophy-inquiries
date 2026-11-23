"""Vacuum orientation of two equal-mass flavor ports and a signed material probe."""
import argparse
import itertools
import json
import math
from pathlib import Path
import runpy

import numpy as np
from numpy.polynomial.legendre import leggauss


HERE = Path(__file__).resolve().parent
def read(path):
    return json.loads((HERE/path).read_bytes())

CFG = read("inputs.json")
ROW = read(CFG["mass_point"])[CFG["case"]]
Y = np.asarray(ROW["down_yukawas_at_M"])
M, SOFT = ROW["M_GeV"], ROW["soft_GeV"]
STAT = CFG["color_dimension"]/(8*math.pi**2)
C, FB = ROW["c"], ROW["fb"]
K = np.array([-Y[0], Y[1], -Y[2]])/Y[2]
LOSS = np.array([1-C**2, 0., 1-FB**2])
T = (SOFT/M)**2
XEW = CFG["higgs_v_GeV"]**2/2


def limits(x, order):
    """Critical signed couplings for all orientations, with endpoint derivatives."""
    z = x/M**2
    s, b = Y[1:]
    weight = LOSS[2]
    delta = K[1]-K[2]
    nodes, gauss = leggauss(order)
    cuts = sorted({-70., -32., -16., 0., math.log(T), 16., 36., 65.})
    result = np.zeros(3)
    for lo, hi in zip(cuts, cuts[1:]):
        q = np.exp((lo+hi)/2+(hi-lo)*nodes/2)
        ratio_over_z = (b*b-s*s)*weight*q/((1+q)*(q+z*s*s)*(q+z*b*b*(1-weight/(1+q))))
        h = z*ratio_over_z
        # Stable integrals for the onset, continuum interval width and secant.
        unit = STAT*M**2 * q*np.minimum(q,T)*ratio_over_z/delta
        lmean = np.divide(np.log1p(h),h,out=np.ones_like(h),where=(h!=0))
        pieces = np.array([unit,unit*h/(1+h),unit*lmean])
        result += (hi-lo)/2*(pieces@gauss)
    upper, gap, crossing = result
    return {"g_for_desired_from_below_GeV2":upper,
            "g_above_which_free_ds_loses_stability_GeV2":upper-gap,
            "width_of_intermediate_rotations_GeV2":gap,
            "g_equal_endpoint_energies_GeV2":crossing}


def daughter(parent, pt, scheme):
    flow = pt["flow"]
    poles = flow["fields"]["neutral_poles_GeV"]
    segments = [a for a in flow["segments"](.125) if a[0]<poles[2]]
    y = parent["initial_state"](pt)
    unit = np.eye(3,dtype=complex)
    for start, susy, dt, grid in segments:
        assert not susy
        for i,pole in enumerate(poles):
            if start==pole:
                y[3,:,i]=flow["physics"]["Y0"][:,i]
        y,_ = parent["advance"](y,unit,dt,grid,flow)
    left, val, vh = np.linalg.svd(y[1,:,:3])
    c = abs(flow["basis"](y[0])[0][:,0].conj()@left[:,-1])
    sb = math.sqrt((1+c)/2)
    f = np.array([sb,1.,c] if scheme=="db" else [1.,sb,c])
    y[1,:,:3] = y[1,:,:3]@(vh.conj().T@np.diag(f)@vh)
    for start, _, dt, grid in reversed(segments):
        y,_ = parent["advance"](y,unit,-dt,grid[::-1],flow)
        for i, pole in enumerate(poles):
            if start==pole:
                y[3,:,i]=0.
    pars, _, gamma, _, ups = flow["extraction"](y)
    mods = dict(zip(("Vus","Vub","Vcb"),np.sin(pars[3:])*[np.cos(pars[4]),1,np.cos(pars[4])]))
    width = {k[1:]:(v-parent["DIRECT"][k][0])/math.hypot(*parent["DIRECT"][k][1:]) for k,v in mods.items()}
    mk = parent["MATCH"]
    ks,_ = mk["thresholds"](dict(moduli_MZ=mods,weak_offsets=width,gamma=gamma),mk["source_thresholds"]())
    h = np.array([pars[2]/pars[1]/ks[2],1.,pars[0]/pars[1]/ks[0],ups[0]/pars[1]/(ks[0]*ks[1])])
    fit = parent["project_shared"](h)
    absolute = parent["absolute_table"](h,ups[0]/ups[1],fit)
    return dict(absolute_four_ratio_radius=absolute["radius_at_Kuc_one"],
                profile_four_ratio_radius=fit["radius"], b_over_s=h[0], d_over_s=h[2],
                MZ_down=pars[:3].tolist(),gamma_deg=gamma)


def solve():
    logfactor = math.log1p(T)+T*math.log1p(1/T)
    slope = STAT*M**2*logfactor
    key = ("d","s","b")
    table = {}
    for i,j in itertools.permutations(range(3),2):
        moment = slope*(Y[i]**2*LOSS[0]+Y[j]**2*LOSS[2])
        tag = K[i]+K[j]
        table[key[i]+key[j]] = {"cw_slope_GeV2":moment, "tag_slope_per_g":tag}
    best, alt = table["db"],table["ds"]
    critical = slope*LOSS[2]*Y[2]*(Y[2]-Y[1])
    crossing = (best["cw_slope_GeV2"]-alt["cw_slope_GeV2"])/(alt["tag_slope_per_g"]-best["tag_slope_per_g"])
    np.testing.assert_allclose(crossing,critical,rtol=1e-12)
    matched = -best["cw_slope_GeV2"]/best["tag_slope_per_g"]
    assert matched>critical
    high = limits(0.,64)
    assert abs(high["g_for_desired_from_below_GeV2"]-critical)<1e-9
    controls={}
    for label, x in [("electroweak",XEW)]+[(str(z),M*M*z) for z in CFG["x_over_M_squared_controls"]]:
        coarse, fine = (limits(x,order) for order in CFG["quadrature_orders"])
        assert max(abs(fine[k]-coarse[k]) for k in fine)<2e-7
        assert matched>fine["g_for_desired_from_below_GeV2"]>fine["g_above_which_free_ds_loses_stability_GeV2"]
        controls[label]={"x_GeV2":x, **fine}
    parent = runpy.run_path(str(HERE/CFG["parent"]))
    data = parent["CFG"]["nearby_measured_inputs_in_marginal_widths"][CFG["case"]]
    _, pt, _ = parent["charged_endpoint"](2,.125,*(data[k] for k in ("weak","up","g3","top","charm")))
    readouts = {mode:daughter(parent,pt,mode) for mode in ("db","ds")}
    original = parent["read"]("results.json")["rows"][CFG["case"]]
    assert abs(readouts["db"]["absolute_four_ratio_radius"]-original["absolute_table_IV"]["radius_at_Kuc_one"])<1e-5
    pole, error = CFG["comparison_Higgs_pole_GeV"]
    return dict(mother_yukawas=Y.tolist(),loss_d_and_half=[LOSS[0],LOSS[2]],signed_probe_eigenvalues=K.tolist(),
                orientation_slope_candidates=table, common_soft_coefficient_GeV2=slope,
                sharp_threshold_all_backgrounds_GeV2=critical,
                g_at_zero_added_slope_GeV2=matched,g_margin_GeV2=matched-critical,
                sufficient_slope_precision_to_exclude_untagged_branch_GeV2=alt["cw_slope_GeV2"],
                tree_single_doublet_mass_error_slope_equivalent_GeV2=pole*error,
                controls=controls,measured_points=readouts)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    output=solve()
    if args.check:
        assert output==read("results.json")
        print("1741 all-background sorting, finite signed portal and same-parent flavor tests verified")
    else:
        (HERE/"results.json").write_text(json.dumps(output,indent=2,allow_nan=False)+"\n",encoding="utf8")
