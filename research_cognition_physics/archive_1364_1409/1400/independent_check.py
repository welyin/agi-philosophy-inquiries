"""Absolute on-shell leptonic width measured with CMS/ATLAS and with the LEP ensemble."""
import argparse
from itertools import product
import json
from math import pi, sqrt
from pathlib import Path


HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE/name).read_text(encoding="utf-8"))


def endpoints(m, w, n):
    return (m-n*w,m+n*w)


def collect(mass, width, b_leptonic, n, muon_constant, inv_alpha):
    z=(91.1876,.0021)
    pools=[endpoints(*x,n) for x in (mass,width,b_leptonic,z)]
    values=[]
    for m,g,b,m_z in product(*pools):
        partial=g*b
        angle=1-(m/m_z)**2
        squared_coupling=48*pi*partial/m
        gamma0=muon_constant*m**3/(6*pi*sqrt(2))
        k_mu=sqrt(2)*muon_constant*m*m*angle*inv_alpha/pi
        k_width=squared_coupling*angle*inv_alpha/(4*pi)
        assert abs(partial/gamma0-k_width/k_mu) < 1e-12
        values.append((partial,k_mu,k_width,partial/gamma0,sqrt(squared_coupling)))
    labels=["partial_width_GeV","ell_squared_R_mu_rho1","ell_squared_R_leptonic_rho1",
            "ratio_R_leptonic_over_R_mu","g_ell_sqrt_R_leptonic"]
    out={label:[min(column),max(column)] for label,column in zip(labels,zip(*values))}
    return out


def run():
    old=read("../1394/inputs.json")["formula"]
    gf, alpha_inv=old["Gmu_GeV_minus2"],old["alpha_inverse"]
    old_bounds=read("../1398/independent_results.json")
    trace_cap=read("../1398/results.json")["widths"]
    data={}
    for n in (0,1,2):
        atlas=collect((80.3548,.0161),(2.198,.049),(.1089,.0008),n,gf,alpha_inv)
        # LEP LFU closure: the hadronic fraction and the total width were combined
        # by LEP collaborations. The outer box does not assume statistical independence.
        lep=collect((80.376,.033),(2.195,.083),((1-.6741)/3,.0027/3),n,gf,alpha_inv)
        ratio_lower=1.984-n*.021
        visible_cap=trace_cap[str(n)]["visible_two_row_strength_upper"]
        ell_sq_sum=visible_cap/ratio_lower
        spectral=old_bounds[str(n)]
        ell_sq_four=(spectral["max_inert_sigma"]/spectral["min_light_sigma"])**2
        min_g_ell=atlas["g_ell_sqrt_R_leptonic"][0]/sqrt(1.1)
        min_g=min_g_ell/sqrt(ell_sq_four)
        assert ell_sq_four<ell_sq_sum<1 and min_g>20
        data[str(n)]={"ATLAS_width_CMS_B":atlas,"LEP_width_and_B":lep,
                      "ell_upper_hadronic":sqrt(ell_sq_sum),"ell_upper_four_direct_moduli":sqrt(ell_sq_four),
                      "R_width_min_rho1_hadronic":atlas["ell_squared_R_leptonic_rho1"][0]/ell_sq_sum,
                      "R_width_min_rho1_four":atlas["ell_squared_R_leptonic_rho1"][0]/ell_sq_four,
                      "gl_min_if_R_width_le_1p1":min_g_ell,
                      "g_min_if_four_channels_inert_and_R_width_le_1p1":min_g,
                      "g2_over_16pi2_min":(min_g/(4*pi))**2}
    return data


def verify(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:
            verify(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b):
            verify(x,y)
    else:
        assert abs(a-b)<1e-10*max(abs(a),1.)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    p.add_argument("--write",action="store_true")
    args=p.parse_args()
    result=run()
    out=HERE/"independent_results.json"
    if args.write:
        with out.open("x",encoding="utf-8") as handle:
            handle.write(json.dumps(result,indent=2,allow_nan=False)+"\n")
    else:
        verify(result,read(out.name))
        print("1400 independent leptonic partial-width and perturbative consistency bound verified")
