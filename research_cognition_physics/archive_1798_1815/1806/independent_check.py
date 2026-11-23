"""Pure type-I polygon constraints on a singlet-generated Majorana portal."""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
V = 246.22


def read(name):
    return json.loads((HERE/name).read_bytes())


def edge_bound(l, n):
    # Normalised pair source (a,b), t=b/a. Its maximum on t in [l,1]
    # occurs at an endpoint, with all n-2 spectators at mass l.
    return max((1-l+l*l)/math.sqrt(1+(n-1)*l**4),
               1/math.sqrt(2+(n-2)*l**4))


def spectrum_bound(n, tolerance, target):
    lo, hi = 0., 1.
    for _ in range(75):
        mid=(lo+hi)/2
        f=edge_bound(mid,n)
        if f + tolerance*(1-f) >= target:
            lo=mid
        else:
            hi=mid
    return 1/hi


def canonical_trace(Y,B):
    # Takagi frame, with M_R=s B. Check all three quartics through the
    # particle/conjugate trace 2 Tr[(M^\dagger M)^2] for complex Y.
    left, right=Y.shape
    h,s=.67,-.41
    d=Y*h/math.sqrt(2)
    full=np.block([[np.zeros((left,left)),d],[d.T,s*B]])
    A=Y.conj().T@Y
    G=B.conj().T@B
    formula=np.trace(A@A).real*h**4+4*np.trace(A@G).real*h*h*s*s+2*np.trace(G@G).real*s**4
    dense=2*np.trace(np.linalg.matrix_power(full.conj().T@full,2)).real
    np.testing.assert_allclose(dense,formula,rtol=4e-14)


def two_neutrinos(k, central, mother_mass):
    # Two explicitly separated poles admit a physical inclusive-width comparison.
    omega, norm=central["omega"],central["S"]
    A0=central["DeltaA"]/omega
    ms=np.array([mother_mass/k,mother_mass])
    osc=read("../1370/inputs.json")["NO"]
    masses=np.sqrt([osc["solar"][0],osc["atm_abs"][0]])*1e-9
    ch=(V**2*A0-(ms[0]-ms[1])*(masses[0]-masses[1]))/(sum(ms)*sum(masses))
    t=math.acosh(ch)/2
    rot=np.array([[math.cosh(t),1j*math.sinh(t)],[-1j*math.sinh(t),math.cosh(t)]])
    phases=np.exp(2j*math.pi*np.outer(np.arange(3),[1,2])/3)/math.sqrt(3)
    Y=math.sqrt(2)/V * (phases*np.sqrt(masses)) @ rot @ np.diag(np.sqrt(ms))
    light=V**2/2*(Y/ms)@Y.T
    desired=(phases*masses)@phases.T
    assert np.linalg.norm(light-desired,2)/max(masses)<3e-10
    A=Y.conj().T@Y
    np.testing.assert_allclose(A.trace().real,A0,rtol=1e-12)
    w=math.sqrt(omega/norm*sum(ms*ms))
    B=np.diag(ms/w)
    canonical_trace(Y,B)
    G=B@B
    lh=(central["B_ch"]+omega*np.trace(A@A).real)/norm
    ls=2*omega*np.trace(G@G)/norm
    lp=2*omega*np.trace(A@G).real/norm
    shift=lp*lp/ls
    heavy2=2*ls*w*w
    # These quadratic terms are set separately to support (v,w).
    # Positive local curvature is verified for these examples.
    assert lh-shift>0
    hh,hs=2*lh*V*V,2*lp*V*w
    heavy=(hh+heavy2+math.hypot(heavy2-hh,2*hs))/2
    widths=ms*np.diag(A).real/(8*math.pi)
    normed=widths/ms**2
    tol=2*np.linalg.norm(light,2)/(8*math.pi*V**2)
    assert abs(normed[0]-normed[1])<=tol*(1+1e-8)
    return dict(M_N_GeV=ms.tolist(),same_vev_w_GeV=w,quartics_H_Hsigma_sigma=[lh,lp,ls],
                shift=shift,remaining_lambda=lh-shift,
                normalized_positive_overlap=np.trace(A@G).real/(A0*np.linalg.norm(G)),
                light_mass_target_eV=(masses*1e9).tolist(), C5_relative_residual_below_3e_minus10=True,
                tree_radial_pole_over_Mmax=math.sqrt(heavy)/max(ms),
                in_diagonal_basis_inclusive_LH_width_GeV=widths.tolist(),
                width_light_over_heavy=widths[0]/widths[1],
                ratio_of_widths_over_M_squared=normed[0]/normed[1])


def compute():
    inputs=read("../1805/results.json")
    seesaw=read("../1805/seesaw_rank_results.json")
    central,box=inputs["central"],inputs["outer"]
    e=seesaw["transport_and_low_energy_bounds"]
    eps=e["max_second_Y_squared"]
    amins=read("../1805/inputs.json")["gauge_ratio_outer"]
    omega_lo,omega_hi=amins
    amin=box["DeltaA_lower"]/omega_hi
    eta=eps/amin
    fraction=(1-2*eta)**2+2*eta**2
    # Quantity needed by the same trace to give the Higgs quartic, with all
    # extrema picked independently and in the conservative direction.
    norm_lo,norm_hi=box["S"]
    dmax=norm_hi-box["A_ch"][0]
    target=math.sqrt(.5*(fraction+omega_lo*(box["B_ch_lower"]-norm_lo*box["lambda_data"][1])/dmax**2))
    c=central
    A0=c["DeltaA"]/c["omega"]
    e0=eps/A0
    t0=math.sqrt(.5*((1-2*e0)**2+2*e0**2 + (c["B_ch"]-c["S"]*c["lambda_data"])/(c["omega"]*A0**2)))
    normed_width_tolerance=2*e["max_mZ_singular_eV"]*1e-9*e["max_SM_operator_growth"]/(8*math.pi*V**2)
    assert target<1<t0
    mother=read("../1369/results.json")["one_loop_central"]["mu_GeV"]
    d_equal=2*c["omega"]*A0*A0/(3*c["S"])
    return dict(conditional_width_triangle_additive_tolerance_inverse_GeV=normed_width_tolerance,
                pureI_width_triangle_dimensionless_tolerance=eta,
                pureI_neutral_B_over_A_squared_box_floor=fraction,
                overlap_required_box=target,overlap_required_at_central=t0,
                mass_spread_necessary={str(n):spectrum_bound(n,eta,target) for n in [2,3]},
                equal_heavy_3_threshold_at_central=d_equal,
                min_central_lambda_after_equal_three_threshold=seesaw["conditional_central_gauge_and_charged_inputs"]["lambda_lower_at_most_three_pure_I"]-d_equal,
                angular_bounds_by_Mmax_over_Mmin={
                    str(k):{str(n):edge_bound(1/k,n)+(1-edge_bound(1/k,n))*eta for n in [2,3]}
                    for k in [1,2,10,30,35,40]},
                two_state_common_matching_examples={str(k):two_neutrinos(k,c,mother) for k in [1,2,40]})


def clean(value):
    if isinstance(value,dict):
        return {key:clean(x) for key,x in value.items()}
    if isinstance(value,list):
        return [clean(x) for x in value]
    if isinstance(value,(float,np.floating)):
        return float(f"{value:.10g}")
    return value


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    args=p.parse_args()
    result=clean(compute())
    blob=json.dumps(result,indent=2,allow_nan=False)+"\n"
    out=HERE/"independent_results.json"
    if args.check:
        assert result==json.loads(out.read_bytes())
    else:
        out.write_bytes(blob.encode())
    print("1806 independent",hashlib.sha256(blob.encode()).hexdigest()[:12])
