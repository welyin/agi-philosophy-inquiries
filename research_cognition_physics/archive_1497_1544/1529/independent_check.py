"""Independent mixed-meson cap and simultaneous symmetric Yukawa diagonalization."""
import argparse
import itertools
import json
import math
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE/path).read_text(encoding="utf-8-sig"))


def mixing(angles):
    a,b,c,d=angles
    s1,s2,s3=np.sin([a,b,c])
    c1,c2,c3=np.cos([a,b,c])
    z=np.exp(1j*d)
    return np.array([[c1*c2,s1*c2,s2/z],
                     [-s1*c3-c1*s3*s2*z,c1*c3-s1*s3*s2*z,s3*c2],
                     [s1*s3-c1*c3*s2*z,-c1*s3-s1*c3*s2*z,c3*c2]])


def evaluate():
    inputs=read("inputs.json")
    primary=read("results.json")
    kte=inputs["ktev"]
    widths=inputs["marginal_widths"]
    dm,w=kte["delta_m_1e6_per_second"]
    ts,ws=kte["tau_s_ps"]
    tl,*wl=inputs["kloe_tau_L_ns"]
    upper=0.5*inputs["hbar_GeV_s"]*math.hypot((dm+widths*w)*1e6,
        1/math.sqrt((ts-widths*ws)*1e-12*(tl-widths*sum(wl))*1e-9))
    np.testing.assert_allclose(upper,primary["two_width_total_M12_upper_GeV"],rtol=2e-13,atol=0)
    # Minimize LR +/− SLL, separately multi-affine in the independent interval variables.
    gauge=read(inputs["gauge"])["sm_dr_at_one_TeV"]["g3"]
    strong=[(gauge[0]+k*widths*gauge[1])**2/(16*math.pi**2) for k in (-1,1)]
    box=inputs["hadronic_relative_box"]
    budgets=[]
    for x,m,sgn in itertools.product(strong,inputs["mass_window_GeV"],(-1,1)):
        L=math.log((m/inputs["mu_GeV"])**2)
        coefficients=[-1.5*x,1-x,sgn*(1+x*(4.5-3*L)),sgn*x*(1/8-L/12)]
        products=np.array(coefficients)*inputs["table1_K_LR1_LR2_SLL1_SLL2_GeV3"]
        budgets.append(sum(products)-box*sum(abs(products)))
    gap=min(budgets)
    np.testing.assert_allclose(gap,primary["NLO_response_gap_GeV3"],rtol=2e-13)
    s=primary["regions"]["matched_rectangle"]["geometric_smallest_singular_value"]
    lower=gap*(s*inputs["frame_transfer_floor"])**2/(2*max(inputs["mass_window_GeV"])**2)
    np.testing.assert_allclose(lower,primary["regions"]["matched_rectangle"]["scalar_M12_abs_lower_GeV"],
                               rtol=2e-13,atol=0)
    angles=np.array(read("../1527/inputs.json")["one_TeV_DR"]["angles_order_12_13_23_delta"])
    e=2*np.array([.00082,.00008,.00043,.023])
    sl=np.sin(angles-e); ch=np.cos(angles+e)
    su=np.sin(angles+e)
    jfloor=sl.prod()*ch[0]*ch[1]**2*ch[2]
    tmax=3-(ch[0]*ch[1])**4-(ch[0]*ch[2]-su[0]*su[2]*su[1])**4-(ch[1]*ch[2])**4
    uniform=48*jfloor*jfloor/tmax**2
    v=mixing(angles)
    b=v*v
    p=abs(v)**2
    smax=np.linalg.svd(b,compute_uv=False)[0]
    J=abs((v[0,0]*v[1,1]*v[0,1].conjugate()*v[1,0].conjugate()).imag)
    det=np.linalg.det(np.eye(3)-b.conj().T@b).real
    floor=4*J**2*(9-np.sum(p*p))
    assert det>=floor*(1-2e-8) and 1-smax>=uniform>0
    # Verify smallest constrained quadratic form on the whole complex symmetric space.
    components=list(itertools.combinations_with_replacement(range(3),2))
    pairs=[(i,j) for i,j in components if i!=j]
    columns=[]
    for i,j in components:
        z=np.zeros((3,3),complex)
        z[i,j]=z[j,i]=1 if i==j else 1/math.sqrt(2)
        w=v.conj().T@z@v.conj()
        columns.append([math.sqrt(2)*h[i,j] for h in (z,w) for i,j in pairs])
    actual=np.linalg.svd(np.array(columns).T,compute_uv=False)[-1]**2
    np.testing.assert_allclose(actual,1-smax,rtol=1e-7)
    coeff=v[0]*v[1]
    d=np.array([coeff[1],-coeff[0],0])
    d/=np.linalg.norm(d)
    yu=v@np.diag(d)@v.T
    assert abs(yu[0,1])<1e-14 and np.linalg.norm(yu[:2,2])>0
    return dict(kaon_cap_GeV=upper,kaon_response_gap_GeV3=gap,additional_over_cap=lower/upper,
                central_J=J,central_defect=actual,determinant=det,four_cycle_floor=floor,
                box_J_floor=jfloor,box_trace_upper=tmax,box_defect_floor=uniform,
                one_absolute_offdiagonal_floor=math.sqrt(uniform/12),
                all_four_standard_mesons_blind_example_top_norm=float(np.linalg.norm(yu[:2,2])))


if __name__=="__main__":
    p=argparse.ArgumentParser()
    g=p.add_mutually_exclusive_group(required=True)
    for flag in ("check","write-new"):
        g.add_argument("--"+flag,action="store_true")
    opt=p.parse_args()
    result=evaluate()
    if opt.check:
        frozen=read("independent_results.json")
        assert frozen.keys()==result.keys()
        for k in frozen:
            np.testing.assert_allclose(frozen[k],result[k],atol=1e-25,rtol=2e-7)
        print("1529 independent: kaon cap, matrix-element box, CP obstruction and four-meson null sample.")
    else:
        with (HERE/"independent_results.json").open("x",encoding="utf8") as f:
            json.dump(result,f,indent=2)
            f.write("\n")
