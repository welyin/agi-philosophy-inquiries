"""Compare the two adjoint-mass roots using the complete broken-generator current."""
import argparse
import json
import math
from pathlib import Path
import runpy

import numpy as np

HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE/path).read_bytes())


def mass_object(rows):
    """Two massive rows in five labeled anti-fundamentals. Return actual light columns."""
    mass = np.zeros((2,5))
    light = np.eye(5)[:,:3].copy()
    for j,i in enumerate((0,2)):
        a,b = rows[j]
        mass[j,[i,3+j]] = [a,b]
        light[[i,3+j],i] = [b,-a]/np.linalg.norm(rows[j])
    np.testing.assert_allclose(mass@light,0,atol=1e-10)
    return mass,light


def kernel(mass):
    return np.linalg.svd(mass,full_matrices=True)[2].conj().T[:,2:]


def bounds(color, weak, mc, mw):
    heavy_overlap = (color@weak.conj().T)/(mc*mw)
    u,cosine,vh = np.linalg.svd(heavy_overlap)
    cosine=np.minimum(1.,cosine)
    norm_min=math.sqrt((mw-mc)**2+2*mc*mw*(1-min(cosine)))/5
    frob_min=math.sqrt(2*(mw-mc)**2+2*mc*mw*sum(1-cosine))/5
    aligned=(color-u@vh@weak)/5
    np.testing.assert_allclose([np.linalg.norm(aligned,2),np.linalg.norm(aligned,"fro")],
                               [norm_min,frob_min],rtol=1e-12)
    return cosine,norm_min,frob_min


def tagged(mc,mw,c,sign):
    f = np.array([c,math.sqrt((1+c)/2)])
    color=np.column_stack((np.sqrt(1-f*f),f))*mc
    k = color.sum(axis=0)
    k/=np.linalg.norm(k)
    mean = np.dot(color[0],k)
    dist = np.linalg.norm(color[0]-mean*k)
    shift = sign*math.sqrt(mw**2-dist**2)-mean
    weak=color+shift*k
    rc,ec=mass_object(color)
    rw,ew=mass_object(weak)
    cross=ec.T@ew
    spurion=(rc-rw)/5
    np.testing.assert_allclose(rw@rw.T,mw**2*np.eye(2),rtol=3e-13)
    # Exact inclusive sums in the full five-plet current, including the partner copies.
    loss_c=25*(spurion@ec).T@(spurion@ec)/mw**2
    loss_l=25*(spurion@ew).T@(spurion@ew)/mc**2
    np.testing.assert_allclose(np.eye(3)-cross@cross.T,loss_c,atol=5e-16)
    np.testing.assert_allclose(np.eye(3)-cross.T@cross,loss_l,atol=5e-16)
    cosines,minimum,frob=bounds(rc,rw,mc,mw)
    np.testing.assert_allclose(np.linalg.svd(kernel(rc).T@kernel(rw),compute_uv=False),
                               np.r_[1., cosines],atol=8e-16)
    assert abs(abs(cross[0,0])-abs(cross[2,2])) < 1e-13
    projected_only=ec[:3].T@ew[:3]
    return dict(color_ordinary_light_amplitudes=f.tolist(),
                weak_ordinary_light_amplitudes=(weak[:,1]/mw).tolist(),
                normalized_inclusive_light_X_rates_d_s_b=np.diag(cross@cross.T).tolist(),
                normalized_rates_if_partner_current_dropped=np.diag(projected_only@projected_only.T).tolist(),
                adjoint_mass_operator_norm_GeV=np.linalg.norm(spurion,2),
                min_over_relative_heavy_labels_operator_GeV=minimum,
                min_over_relative_heavy_labels_Frobenius_GeV=frob,
                weak_light_to_heavy_color_adjoint_norm_GeV=np.linalg.norm(spurion@ew,2),
                color_light_to_heavy_weak_adjoint_norm_GeV=np.linalg.norm(spurion@ec,2))


def split_running(mw):
    model=runpy.run_path(str(HERE.parent/"1742"/"check.py"))
    # Keep the registered colour poles. Replace all three old extra weak poles by
    # the two actual weak fermion and scalar pairs on their own thresholds.
    rules=model["thresholds"]([])+[(mw,4*model["WEAK"]/3),
        (math.hypot(mw,model["BASE"]["soft_GeV"]),2*model["WEAK"]/3)]
    transfer=model["transport"]
    transfer.__globals__["thresholds"]=lambda unused: rules
    run=transfer([])
    ref=model["REFERENCE"]
    return run["total_duql_for_comparison"],4*math.pi/(ref["U"]**2*ref["new_inverse_colour"])


def mesons(light_d,ckm,gain,canonical):
    pion=read("../1520/inputs.json")
    kaon=read("../1483/inputs.json")
    mp=kaon["mass_GeV"]["p"]
    ud,us=ckm
    wpi=abs(pion["pion_W0_GeV2"][0])
    w1,w2=(abs(kaon["W0_GeV2"][k][0]) for k in ("ud_s","us_d"))
    def width(mass,amplitude_squared):
        return (gain*canonical)**2 * mp*(1-(mass/mp)**2)**2*amplitude_squared/(32*math.pi)
    gp=width(pion["pion_mass_GeV"],(ud*wpi)**2*light_d)
    # On the actual roots d and s rows remain orthogonal in the inclusive light space.
    gk=width(kaon["mass_GeV"]["K+"],(us*w1)**2*light_d+(ud*w2)**2)
    unit=pion["hbar_GeV_s"]/pion["year_seconds"]
    return dict(conditional_pion_lifetime_years=unit/gp,conditional_kaon_lifetime_years=unit/gk,
                charged_kaon_to_pion_width_ratio=gk/gp,
                ratios_to_separate_published_rate_limits=[pion["pion_lifetime_years"]*gp/unit,
                                                         kaon["proton"]["partial_lifetime_lower_yr"]*gk/unit])


def random_test(mc,mw):
    rng=np.random.default_rng(1745)
    worst=0.
    for _ in range(36):
        frames=[np.linalg.qr(rng.normal(size=(5,5))+1j*rng.normal(size=(5,5)))[0] for i in range(2)]
        rc=mc*frames[0][:,:2].conj().T
        rw=mw*frames[1][:,:2].conj().T
        ec,ew=(f[:,2:] for f in frames)
        cross=ec.conj().T@ew
        v=(rc-rw)/5
        q=v@ec
        err=np.linalg.norm(np.eye(3)-cross@cross.conj().T-25*q.conj().T@q/mw**2)
        worst=max(worst,err)
        cosine,minimum,_=bounds(rc,rw,mc,mw)
        np.testing.assert_allclose(np.linalg.svd(cross,compute_uv=False),np.r_[1.,cosine],atol=2e-14)
        assert minimum<=np.linalg.norm(v,2)+1e-10
    assert worst<4e-14
    return worst


def solve():
    model=runpy.run_path(str(HERE.parent/"1744"/"check.py"))
    b=model["mother"](.125)
    mc,mw=b["Mcolor"],b["Mweak"]
    cc=model["daughter"](b,np.eye(3))["color_scale_c"]
    gain,coupling=split_running(mw)
    # Use the physical 10 PeV labels of the mother; extrapolation to X is an explicit assumption.
    ckm=abs(np.linalg.svd(b["pt"]["D"])[0][0,::-1][:2])
    rows={}
    for sign,name in [(1,"near"),(-1,"reflected")]:
        row=tagged(mc,mw,cc,sign)
        weak_reference=tagged(mc,mw,b["cW"],sign)
        row["current_rate_change_if_c_fixed_at_weak_pole"]=(
            weak_reference["normalized_inclusive_light_X_rates_d_s_b"][0]
            - row["normalized_inclusive_light_X_rates_d_s_b"][0])
        row["mesons"]=mesons(row["normalized_inclusive_light_X_rates_d_s_b"][0],ckm,gain,coupling)
        rows[name]=row
    return dict(color_mass_GeV=mc,weak_mass_GeV=mw,reference_c_at_weak_threshold=b["cW"],
                reference_c_at_color_threshold=cc,
                CKM_10PeV_ud_us=ckm.tolist(),conditional_U_matching_GeV_minus2=coupling,
                duql_gain_to_2GeV=gain,cases=rows,
                complex_frame_max_identity_residual=random_test(mc,mw))


def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:
            compare(a[k],b[k])
    else:
        np.testing.assert_allclose(a,b,rtol=2e-9,atol=1e-13)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    output=HERE/"independent_results.json"
    result=solve()
    if args.check:
        compare(result,json.loads(output.read_bytes()))
        print("1745 complete broken current, adjoint mass identities and meson widths verified")
    else:
        with output.open("x",encoding="utf8") as fp:
            json.dump(result,fp,indent=2)
            fp.write("\n")
