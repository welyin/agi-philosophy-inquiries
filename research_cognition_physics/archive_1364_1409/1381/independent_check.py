"""Dimension of a measured weak-current block, from singular values and row polygons."""
import argparse
import json
from cmath import exp
from math import acos, pi, sqrt
from pathlib import Path

import numpy as np


A0,D0,R0,C0=.97373,.0351/.93,.083,.2238
K0,CS0=.2157/.9704,.9589
CW=sum([.0011,.0015,.0004,.0002,.0022])


def block(a,b,c,d,degree):
    return np.array([[a,b],[c,-d*exp(1j*degree*pi/180)]])


def transition(a,b,d,gamma):
    # The largest singular value crosses one as c increases in this rectangle.
    lo,hi=.18,.28
    for _ in range(55):
        mid=(lo+hi)/2
        if np.linalg.svd(block(a,b,mid,d,gamma),compute_uv=False)[0]<1:
            lo=mid
        else:
            hi=mid
    return (lo+hi)/2


def extend_rows(rows):
    dim=rows.shape[1]
    assert np.allclose(rows@rows.conj().T,np.eye(2),atol=5e-14,rtol=0)
    # Complete via the kernel of the measured rows. The top block stays fixed.
    _,_,vt=np.linalg.svd(rows,full_matrices=True)
    out=np.vstack([rows,vt[2:,:]])
    assert np.allclose(out@out.conj().T,np.eye(dim),atol=5e-14,rtol=0)
    return out


def fourway(gamma):
    v=block(A0,D0*R0,C0,D0,gamma)
    t=sqrt(1-v[0].real@v[0].real-K0**2)
    q=sqrt(1-np.vdot(v[1],v[1]).real-CS0**2)
    h=-(v[0]@v[1].conj())
    first,second=K0*CS0,t*q
    assert abs(first-second)<abs(h)<first+second
    deviation=acos((abs(h)**2+first**2-second**2)/(2*abs(h)*first))
    a=first*h/abs(h)*exp(1j*deviation)
    p=h-a
    assert abs(abs(p)-second)<1e-14
    rows=np.column_stack((v,[K0,np.conj(a)/K0],[t,np.conj(p)/t]))
    out=extend_rows(rows)
    assert np.allclose(abs(out[1,:]),[C0,D0,CS0,q])
    return {"residual_row_moduli":[t,q],"known_pair_sum":abs(h),
            "triangle_range":[abs(first-second),first+second],
            "full_rows_real":out.real.tolist(),"full_rows_imag":out.imag.tolist()}


def threeway(a,d,r,gamma):
    b=d*r
    c=transition(a,b,d,gamma)
    v=block(a,b,c,d,gamma)
    k=sqrt(1-a*a-b*b)
    h=-(v[0]@v[1].conj())
    u=extend_rows(np.column_stack([v,[k,h.conjugate()/k]]))
    return {"a":a,"d":d,"ratio":r,"angle":gamma,"c":c,"k":k,"cs":abs(u[1,2]),
            "full_rows_real":u.real.tolist(),"full_rows_imag":u.imag.tolist()}


def range_at_width(m):
    amin,amax=A0-m*.00031,A0+m*.00031
    # All historic errors expanded together, lower form factor kept fixed.
    dmax=(.0351+m*(.0019+.0018+.0008))/.8
    bmax=(R0+m*.008)*dmax
    pmin,pmax=1-amax**2-bmax**2,1-amin**2
    kl,kh=(.2157-m*.0006)/(.9704+m*.0032),(.2157+m*.0006)/(.9704-m*.0032)
    cmin=sqrt((1-amax**2-bmax**2)*(1-bmax**2-dmax**2))-amax*bmax*dmax
    cmax=(sqrt(1-amin**2)+amax*bmax*dmax)/(1-bmax**2)
    return {"three_c_outer":[cmin,cmax],"three_k_outer":[sqrt(pmin),sqrt(pmax)],
            "K_semileptonic_box":[kl,kh],"first_row_missing_prob_box":[pmin-kh**2,pmax-kl**2]}


def compute():
    bn=D0*R0
    central_v=block(A0,bn,C0,D0,68.7)
    central_eig=np.linalg.eigvalsh(np.eye(2)-central_v@central_v.conj().T)
    assert np.min(central_eig)>0
    allphase=[transition(A0,bn,D0,g) for g in [180,0]]
    assert C0<min(allphase)
    # Independent directional bound: the known row inner product is at most ac+bd.
    # Even the least possible determinant at the four central moduli is positive.
    p1,p2=1-A0**2-bn**2,1-C0**2-D0**2
    detmin=p1*p2-(A0*C0+bn*D0)**2
    assert detmin>0
    width1,width2=range_at_width(1),range_at_width(2)
    assert width1["first_row_missing_prob_box"][0]>0>width2["first_row_missing_prob_box"][0]
    c_only=threeway(A0,D0,R0,68.7)
    assert abs(c_only["c"]-C0)<CW
    # A 2-width witness for all six moduli and both alternative phase windows.
    two=[threeway(A0+2*.00031,.0495,R0+2*.008,ang) for ang in [68.7,248.7]]
    for sol in two:
        assert abs(sol["c"]-C0)<2*CW
        assert width2["K_semileptonic_box"][0]<sol["k"]<width2["K_semileptonic_box"][1]
        assert abs(sol["cs"]-CS0)<2*sum([.0023,.0040,.0015,.0005,.0095])
    kaon_high=width1["K_semileptonic_box"][1]
    k_three_low=width1["three_k_outer"][0]
    one_input_changes={"needed_max_form_factor_for_fixed_product_upper":.2163/k_three_low,
        "needed_min_first_entry_for_fixed_K_upper":sqrt(1-kaon_high**2-(.0495*.091)**2),
        "needed_min_product_for_fixed_form_factor_lower":k_three_low*.9672}
    return {"central_min_det_any_phase":detmin,"central_defect_eigenvalues_68deg":central_eig.tolist(),
            "one_width_change_one_input_thresholds":one_input_changes,
            "three_c_central_all_phases":allphase,"one_width":width1,"two_widths":width2,
            "three_dim_c_only_one_width":c_only,
            "three_dim_all_six_two_widths":two,
            "four_dim_all_six_central_matched_phase_68":fourway(68.7),
            "four_dim_all_six_central_matched_phase_248":fourway(248.7)}


def rounded(x):
    if isinstance(x,dict):
        return {k:rounded(v) for k,v in x.items()}
    if isinstance(x,list):
        return [rounded(v) for v in x]
    return float(format(x,".10g")) if isinstance(x,(float,np.floating)) else x


def check(actual, saved):
    # SVD nullspace bases can rotate across libraries. Compare the realized rows,
    # and verify all saved full matrices separately.
    if isinstance(actual,dict):
        assert actual.keys()==saved.keys()
        if "full_rows_real" in actual:
            u=np.array(saved["full_rows_real"])+1j*np.array(saved["full_rows_imag"])
            assert np.allclose(u@u.conj().T,np.eye(len(u)),rtol=0,atol=3e-9)
            v=np.array(actual["full_rows_real"])+1j*np.array(actual["full_rows_imag"])
            assert np.allclose(v[:2],u[:2],atol=2e-9,rtol=0)
        for k in actual.keys()-{"full_rows_real","full_rows_imag"}:
            check(actual[k], saved[k])
    elif isinstance(actual,list):
        assert len(actual)==len(saved)
        for a,b in zip(actual,saved):
            check(a,b)
    else:
        assert np.isclose(actual,saved,rtol=4e-9,atol=4e-12)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--write",action="store_true")
    p.add_argument("--check",action="store_true")
    args=p.parse_args()
    file=Path(__file__).with_name("independent_results.json")
    answer=rounded(compute())
    if args.write:
        with file.open("x",encoding="utf-8") as f:
            json.dump(answer,f,indent=2,ensure_ascii=False)
            f.write("\n")
    else:
        check(answer,json.loads(file.read_text(encoding="utf-8")))
    print("1381 phase-independent bounds, SVD transitions and six-modulus row completions verified.")
