"""Two paramagnetic EDM directions applied to the two shared Yukawa sources."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import runpy

import numpy as np


HERE = Path(__file__).resolve().parent
V = 246.22
VB = .95
VA = math.sqrt(1-VB**2)
ME = .00051099895
MN = .939  # isospin-symmetric overall normalization; the p,n scalar fractions differ
# D^equiv = d_e + alpha (Z/A C_p+(1-Z/A) C_n), positive-alpha convention.
MOLECULES = {"HfF": (72/180, .9e-20, 4.1e-30), "ThO": (90/232,1.5e-20,1.1e-29)}
SMALL = (.0139,.0253,.0116,.0302)  # BMW: fu_p,fd_p,fu_n,fd_n
STRANGE = .030/MN  # independent ETMC sigma_s, GeV
BANDS = ("A_vanishes_first_edge","A_vanishes_third_edge")


def sources():
    old = runpy.run_path(str(HERE.parent/"1795/check.py"))
    sig = runpy.run_path(str(old["OLD"]["HERE"]/old["OLD"]["CFG"]["charged_background"]))
    fl,_ = sig["charged_background"](.1)
    dm,em = (np.array(fl["central"](i))/math.sqrt(2) for i in ("down","charged_leptons"))
    out={}
    for name,(_,d,e) in old["sources"]()[0].items():
        if name in BANDS:
            n = (d-VB**2*np.diag(dm))/(VB*VA)
            out[name] = (np.diag(n)/dm,(e[0,0]/em[0]-VB**2)/(VB*VA),n,dm)
    return out


def hadronic(nq, light, fs):
    fu,fd = light
    return fd*nq[0]+fs*nq[1]+2*(1-fu-fd-fs)*nq[2]/27


def coeff(e, n, m, split):
    a=1/(m*m+V**2*split)
    b=1/(m*m)
    prod=(a-b)*e*n
    return ME*MN*np.array([(a+b)*(e*n.conjugate()).imag,prod.imag,prod.real])


def envelope(h):
    amp = math.hypot(h[1],h[2])
    return [float(h[0]-amp),float(h[0]+amp)]


def case(nq,e,m,fs,split=1):
    nuc = [hadronic(nq,SMALL[:2],fs),hadronic(nq,SMALL[2:],fs)]
    vals = [coeff(e,x,m,split) for x in nuc]
    obs = {}
    hm = []
    for name,(frac,k,_) in MOLECULES.items():
        c = frac*vals[0]+(1-frac)*vals[1]
        hm.append(k*c)
        obs[name] = dict(CS_all_phi=envelope(c), dipole_for_de_zero_all_phi_ecm=envelope(k*c))
    diff = hm[1]-hm[0]
    inv = vals[0]-vals[1]
    q = sum(r[2] for r in MOLECULES.values())
    ends = envelope(diff)
    lim = max(map(abs,ends))
    # Verify products with direct scalar and pseudoscalar vertices, independent of harmonic formula.
    for phi in [-.79, .0, .91]:
        phase = np.exp(1j*phi)
        h = np.array([1.,math.cos(2*phi),math.sin(2*phi)])
        for x,v in zip(nuc,vals):
            direct=2*ME*MN*((phase*e).imag*(phase*x).real/(m*m+V*V*split)
                            -(phase*e).real*(phase*x).imag/(m*m))
            np.testing.assert_allclose(direct,v@h,rtol=3e-13,atol=1e-26)
    return dict(MA_GeV=m, MH_GeV=math.sqrt(m*m+V**2*split),
                scalar_harmonics_Cp_Cn=[v.tolist() for v in vals],molecules=obs,
                ThO_minus_HfF_sine_cos_coeff_order_1_cos_sin_ecm=diff.tolist(),
                ThO_minus_HfF_all_phi_ecm=ends,phase_span_ecm=ends[1]-ends[0],
                fraction_of_joint_marginal_sum=lim/q,
                Cp_minus_Cn_all_phi=envelope(inv))


def run():
    main=runpy.run_path(str(HERE/"check.py"))
    out={}
    matter=sources()
    for label,(nq,e,nfull,dm) in matter.items():
        cases=[case(nq,e,m,STRANGE,split) for m,split in [(500,1),(5000,1),(24000,4),(30000,1)]]
        srange=[max(0.,(.030-.017)/MN),(.030+.017)/MN]
        for data in cases:
            m=data["MA_GeV"]
            split = 4 if m==24000 else 1
            c0,c1=[case(nq,e,m,f,split) for f in srange]
            env=[min(c0["ThO_minus_HfF_all_phi_ecm"][0],c1["ThO_minus_HfF_all_phi_ecm"][0]),
                 max(c0["ThO_minus_HfF_all_phi_ecm"][1],c1["ThO_minus_HfF_all_phi_ecm"][1])]
            data["delta_all_phi_one_strange_width_ecm"]=env
            data["delta_degenerate_poles_all_phi_ecm"]=case(nq,e,m,STRANGE,0)["ThO_minus_HfF_all_phi_ecm"]
            # Eliminate strange and universal heavy response through the p−n difference.
            np.testing.assert_allclose(c0["Cp_minus_Cn_all_phi"],c1["Cp_minus_Cn_all_phi"],
                                       rtol=2e-12,atol=1e-27)
            gain=main["mass_run"](m)
            k=main["full_kernel"](dm*V*gain,m,split*V**2)
            val=np.sum(k*(nfull*gain)*(nfull.T*gain))
            qrow=-np.array([val.imag,val.real])
            h=np.array(data["ThO_minus_HfF_sine_cos_coeff_order_1_cos_sin_ecm"])
            normalized=np.array([qrow/np.linalg.norm(qrow), h[1:]/np.linalg.norm(h[1:])])
            orth=np.array([qrow[1],-qrow[0]])/np.linalg.norm(qrow)
            gap=abs(h[1:]@orth)*2
            data["signed_neutral_strong_threshold_and_semileptonic"]={
                "strong_row_cos_sin":qrow.tolist(),"angle_rows_abs_sine":abs(float(np.linalg.det(normalized))),
                "normalized_singular_values":np.linalg.svd(normalized,compute_uv=False).tolist(),
                "two_molecular_values_when_neutral_strong_threshold_zero_ecm":
                    [h[0]-gap/2,h[0]+gap/2],"gap_ecm":gap}
        hist={}
        for fs,tag in [(STRANGE,"ETMC_central"),(.113,"BMW_central")]:
            hist[tag]=case(nq,e,5000,fs)["ThO_minus_HfF_all_phi_ecm"]
        # With degenerate poles and any common up-type real Yukawa shift t,
        # a third aggregate contributes t [fu+2fc], with 2 from c,t.
        def invar(t):
            ans=[]
            for frac,k,_ in MOLECULES.values():
                nn=[]
                for light in (SMALL[:2],SMALL[2:]):
                    fu,fd=light
                    extra=fu+4*(1-fu-fd-STRANGE)/27
                    nn.append(hadronic(nq,light,STRANGE)+extra*t)
                n=frac*nn[0]+(1-frac)*nn[1]
                ans.append(k*(e*n.conjugate()).imag)
            return ans[1]-ans[0]
        t=-invar(0)/(invar(1)-invar(0))
        assert abs(invar(t)) < 1e-33
        out[label]=dict(N_down_over_mq_by_v=[[z.real,z.imag] for z in nq],
                        N_e_over_me_by_v=[e.real,e.imag],real_common_up_return_blind_at_degeneracy=t,
                        archived_lattice_comparison_5TeV=hist,mass_cases=cases)
    branch_gaps={}
    for m,split in [(5000,1),(24000,4),(30000,1)]:
        # Even with independent phases, extrema over a common f_s interval occur at an end:
        # the difference center is affine and the sum of the two radii is convex.
        margins=[]
        for fs in srange:
            a,b=[case(*matter[label][:2],m,fs,split)["ThO_minus_HfF_all_phi_ecm"] for label in BANDS]
            margins.append(a[0]-b[1])
        branch_gaps[str(m)]=min(margins)
    return dict(convention="electronic and proton/neutron scalar operators only; isospin-weighted nuclear impulse",
                inputs=dict(v=V,vB_over_v=VB,me=ME,mN_iso=MN,fu_fd_pn=SMALL,sigma_s_GeV=[.030,.017],
                            molecule_Z_over_A_alpha_ecm_and_published_absolute_limits=MOLECULES),
                minimum_first_minus_third_at_common_strange_fraction_and_any_phases_ecm=branch_gaps,cases=out)


if __name__=="__main__":
    ap=argparse.ArgumentParser(); ap.add_argument("--check",action="store_true"); args=ap.parse_args()
    data=json.dumps(run(),indent=2,ensure_ascii=False,allow_nan=False)+"\n"
    path=HERE/"independent_results.json"
    if args.check:
        assert json.loads(path.read_text(encoding="utf-8"))==json.loads(data)
        print("1797 independent OK",hashlib.sha256(data.encode()).hexdigest()[:12])
    else:
        path.write_text(data,encoding="utf-8")
