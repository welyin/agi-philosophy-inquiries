"""Inclusive Z width and K -> pi invisible for a measured first-row defect."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def read(p):
    return json.loads((HERE/p).read_bytes())


# The second SM value is an alternative input-parameter evaluation, not
# another Z measurement. Sources and qualifications are in independent_review.
EXPERIMENT = (20.767,.025)
BASELINES = {"de_Blas_conservative_individual":(20.7458,.0087),
             "CMS_EWPO_table4":(20.7578,.0086)}
SW2, AS = .2315,.118
DOWN_LEFT = .5 - SW2/3
LEPTON = (-.5+SW2)**2+SW2**2
KZ = 3*(1+AS/math.pi)/LEPTON


def width(D):
    return KZ*(DOWN_LEFT*np.trace(D)-np.trace(D@D)/4).real


def complete_row(row):
    first = np.asarray(row,dtype=float)
    first = first/np.linalg.norm(first)
    second = np.cross(first,[0,0,1.])
    second /= np.linalg.norm(second)
    return np.stack([first,second,np.cross(first,second)])


def check_projector(V):
    D = np.eye(3)-V.conj().T@V
    L = np.eye(3)-V@V.conj().T
    np.testing.assert_allclose(np.linalg.eigvalsh(D),np.linalg.eigvalsh(L),atol=1e-14)
    g = -DOWN_LEFT*np.eye(3)+D/2
    loss = KZ*(3*DOWN_LEFT**2-np.linalg.norm(g)**2)
    np.testing.assert_allclose(width(D),loss,rtol=2e-10)
    return {"D_dd":float(D[0,0].real),"abs_D_ds":abs(D[0,1]),"trace_D":np.trace(D).real,
            "first_row_defect":L[0,0].real,"R_l_loss":loss}


def data_ranges(shift):
    out = {}
    for name,(ref,dr) in BASELINES.items():
        obs,do = EXPERIMENT
        boxes = {}
        for n in (1,2):
            model=[ref-n*dr-1.1*shift,ref+n*dr-.9*shift]
            measured=[obs-n*do,obs+n*do]
            boxes[str(n)] = {"nominal_model_interval_with_response_allowance":model,
                             "experiment_interval":measured,
                             "intersect":bool(max(model[0],measured[0])<=min(model[1],measured[1]))}
        out[name]=boxes
    return out


def run():
    cfg=read("inputs.json")
    mod=read(cfg["quark"])["direct_semileptonic_moduli"]
    root=abs(complex(*read(cfg["branch"])["exact_and_conditional_implementations"]["second_complex_root"]["z"]))
    s,b=(mod[name][0] for name in ("Vus","Vub"))
    delta=1-root**2-s**2-b**2
    vo=np.array([root,s,b])
    # Extremizer has an up-basis rank-one loss. All missing weight is in u.
    A=np.diag([math.sqrt(1-delta),1,1])@complete_row(vo)
    V_unit=complete_row([math.sqrt(1-s**2-b**2),s,b])
    diag_d=delta/(1-s**2-b**2)
    B=V_unit@np.diag([math.sqrt(1-diag_d),1,1])
    for V in (A,B):
        np.testing.assert_allclose(V[0],vo,rtol=3e-14)
    minimum=check_projector(A)
    diagonal=check_projector(B)
    np.testing.assert_allclose(minimum["R_l_loss"],KZ*delta*(DOWN_LEFT-delta/4),rtol=3e-12)
    np.testing.assert_allclose(minimum["abs_D_ds"],delta*root*s/(1-delta),rtol=3e-12)
    assert diagonal["abs_D_ds"]<1e-14 and diagonal["R_l_loss"]>minimum["R_l_loss"]
    # Conservatively normalize the kaon upper bound using the total non-Z
    # amplitude norm and an explicit transfer floor for the rate.
    ke3=.04965-2*math.hypot(.00038,.00037)
    na=(13.+2*3.3)*1e-11
    other=2e-10
    sup=s+2*np.linalg.norm(mod["Vus"][1:])
    kaon_co=.8*1.5*ke3/sup**2
    ds_limit=(math.sqrt(na)+math.sqrt(other))/math.sqrt(kaon_co)
    floor=max(0.,math.sqrt(kaon_co)*minimum["abs_D_ds"]-math.sqrt(other))**2
    assert floor>na
    miss=minimum["abs_D_ds"]-ds_limit
    vmin=2*miss/(math.sqrt(delta+4*miss)+math.sqrt(delta))
    trace_min=delta+(1-delta)*vmin**2
    with_kaon=KZ*trace_min*(DOWN_LEFT-trace_min/4)
    # For trace >= 1/2, the general linear estimate already exceeds this floor.
    assert .5*KZ*(DOWN_LEFT-.25)>with_kaon
    assert minimum["R_l_loss"]<with_kaon<diagonal["R_l_loss"]
    actual=read("results.json")["coloured_five_tree_light_d_s_b_losses_R800"]
    dref=np.diag(actual)
    if min(np.linalg.eigvalsh(dref))<0:
        raise ValueError("Invalid frozen singlet fractions")
    # One- and two-width simultaneous boxes are deterministic envelopes.
    # Include the root (frozen) and independently vary the two other moduli.
    spreads=[np.linalg.norm(mod[key][1:]) for key in ("Vus","Vub")]
    flo=[]
    for n in (0,1,2):
        dd=1-root**2-(s+n*spreads[0])**2-(b+n*spreads[1])**2
        flo.append(KZ*max(0.,dd)*(DOWN_LEFT-max(0.,dd)/4))
    apv=read("../1462/inputs.json")
    qsm,err=apv["baseline_Cs"]
    ah={}
    for name,(measurement,stat,syst) in apv["atomic_branches"].items():
        ah[name]=abs(qsm+211*diag_d-measurement)/(stat+syst+err)
    return {"gamma_normalization":KZ, "target_physical_row":vo.tolist(),
            "frozen_R800_R_l_shift":width(dref),
            "least_Z_loss_row":minimum, "diagonal_d_competitor":diagonal,
            "loss_minimum_at_0_1_2_widths_in_other_quark_inputs":flo,
            "R_l_conditional_comparisons_for_central_light_quark_input":data_ranges(minimum["R_l_loss"]),
            "rare_kaon_two_widths_and_extra_budgets":{"KLOE_Ke3_lower":ke3,
                    "NA62_upper":na,"other_amplitude_solo_BR_ceiling":other,"relative_rate_coefficient_floor":.8,
                    "physical_Vus_upper":sup,"max_D_ds_modulus":ds_limit,
                    "minimum_total_BR_in_extremal_orientation":floor,
                    "extremal_over_allowed_D_ds":minimum["abs_D_ds"]/ds_limit,
                    "trace_floor_after_rare_K":trace_min,
                    "R_l_loss_floor_after_rare_K":with_kaon,
                    "non_FC_witness_R_l_loss":diagonal["R_l_loss"]},
            "diagonal_d_PV_shifts_at_fixed_phys_G_F":{"proton":diag_d,"Cs133":211*diag_d,
                    "distance_in_linear_marginal_widths_Cs_older_alternatives":ah}}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=run()
    target=HERE/"independent_results.json"
    if args.check:
        assert read("independent_results.json")==result
        print("1730 independent neutral-current width and kaon checks passed")
    else:
        target.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf-8")
