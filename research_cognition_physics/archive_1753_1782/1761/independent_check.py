"""Invert the orthonormal two-port source with all unmeasured Majorana phases."""
import argparse
import json
import math
import runpy
from pathlib import Path

import numpy as np
from numpy.polynomial import Polynomial


HERE = Path(__file__).resolve().parent
OLD = runpy.run_path(str(HERE/"../1758/check.py"))
STURM = runpy.run_path(str(HERE/"../1757/independent_check.py"))["sturm_real_count"]


def scalar(value):
    return {"real": float(value.real), "imag": float(value.imag)}


def centered_roots(a, b):
    top = [a+b, 1j*(b-a)]
    aa, dd, bb = (Polynomial([entry[i, j] for entry in top])
                  for i, j in [(0, 0), (1, 1), (0, 1)])
    prod = aa*dd*Polynomial(bb.coef.conjugate())**2
    poly = Polynomial(prod.coef.imag)
    assert abs(poly.coef[-1]) > 1e-12*np.max(abs(poly.coef))
    zeros = poly.roots()
    real = sorted(z.real for z in zeros if abs(z.imag)<1e-8)
    assert STURM(poly.coef) == len(real)
    return real, poly.coef


def channels(W, target):
    u, masses, _ = np.linalg.svd(target)
    first = np.outer(u[:, 0], u[:, 0]) * (u[:, 0].conj()@target@u[:, 0].conj())
    second = target-first
    assert np.linalg.svd(second, compute_uv=False)[1] < masses[1]*1e-12
    reverse = np.linalg.pinv(W)
    return [(reverse@k@reverse.T) for k in [first,second]], first, second


def reconstruct(W, target, poles):
    scale = np.linalg.svd(target, compute_uv=False)[0]
    target = target/scale
    (aa, bb), first, second = channels(W, target)
    roots, coeffs = centered_roots(aa, bb)
    # A zero off-diagonal would add a diagonal solution, including a double root.
    offdiag_min = abs(abs(aa[0,1])-abs(bb[0,1]))
    assert offdiag_min>1e-7
    # Squeezing diagonal moduli can never replace a missing phase. A positive
    # bound on Re(a*d*b*²) excludes every Majorana phase in one step.
    def sign_function(phi):
        J=aa+np.exp(1j*phi)*bb
        return (J[0,0]*J[1,1]*J[0,1].conjugate()**2).real
    count=32768
    phi=np.arange(count)*(2*np.pi/count)
    aa0,aa1=aa[0,0],aa[1,1]
    bd0,bd1=bb[0,0],bb[1,1]
    ac,bc=aa[0,1],bb[0,1]
    caps=np.array([abs(aa0)+abs(bd0),abs(aa1)+abs(bd1),abs(ac)+abs(bc)])
    lip=(abs(bd0)*caps[1]+abs(bd1)*caps[0])*caps[2]**2+2*abs(bc)*caps.prod()
    grid = ((aa0+np.exp(1j*phi)*bd0)*(aa1+np.exp(1j*phi)*bd1)*
            (ac+np.exp(1j*phi)*bc).conjugate()**2).real
    lower=float(np.min(grid)-lip*math.pi/count-1e-12)
    assert abs(sign_function(0)-grid[0])<1e-13
    assert lower>0
    caps4=np.array([caps[0],caps[1],caps[2],caps[2]])
    lo,hi=0.,1.
    for _ in range(56):
        mid=(lo+hi)/2
        if np.prod(caps4+mid)-np.prod(caps4)<lower:
            lo=mid
        else:
            hi=mid
    smin=np.linalg.svd(W,compute_uv=False)[-1]
    low_distance=smin**2*lo
    rotation_floor=np.degrees(2*np.arcsin(low_distance/4))
    kernel=np.linalg.svd(W,full_matrices=True)[0][:,-1]
    assert min(abs(kernel))>.1
    points = []
    for t in roots:
        z = (1+1j*t)/(1-1j*t)
        J = aa+z*bb
        a, b, d = J[0,0], J[0,1], J[1,1]
        ratio2 = -(a*b.conjugate())/(b*d.conjugate())
        assert abs(ratio2.imag)<1e-10
        assert ratio2.real<0
        points.append({"relative_Majorana_rotation_deg": np.angle(z, deg=True),
                       "required_signed_square": ratio2.real})
    # Independent new physical coupling: on entrance and exit insert a
    # reciprocal (unitary and adjoint) rotation relative to nondegenerate poles.
    J=aa+bb
    U,sv,_=np.linalg.svd(J)
    diag=np.diag(U.conj().T@J@U.conjugate())
    takagi=U@np.diag(np.sqrt(diag/abs(diag)))
    np.testing.assert_allclose(takagi@np.diag(sv)@takagi.T,J,atol=1e-12)
    h=math.sqrt(sv[0]/sv[1])
    g=math.sqrt(h)
    x1,x2=poles**2
    msq=(x2-g*x1)/(g-1)
    assert msq>0
    R=takagi
    sorted_transfer=(msq+x2)/(msq+poles**2)
    side=R@np.diag(sorted_transfer**2)@R.conj().T
    P=R@R.T
    rebuilt=sv[1]*W@side@P@side.T@W.T
    np.testing.assert_allclose(rebuilt,target,atol=1e-12)
    rotated={"source_condition_number":sv[0]/sv[1], "leg_condition":h,"each_stage_ratio":g,
             "common_mediator_mass_GeV":math.sqrt(msq),
             "angle_deg_new_axis_to_original_first_pole":np.degrees(np.arccos(abs(R[0,0]))),
             "axis_projector_11":abs(R[0,0])**2, "axis_projector_12":scalar(R[0,0]*R[1,0].conjugate()),
             "max_light_matrix_residual_at_unit_m3":np.max(abs(rebuilt-target))}
    zero = np.linalg.svd(target, compute_uv=False)
    return {"mass_ratio": zero[0]/zero[1], "classifying_quartic": coeffs,
            "all_real_roots":len(roots),"off_diagonal_modulus_floor":offdiag_min,
            "min_Re_ad_bstar_squared_whole_Majorana_circle_bound":lower,
            "source_entry_distance_lower_in_units_low_m3":lo,
            "low_operator_distance_lower_in_units_low_m3":low_distance,
            "necessary_min_rotation_degrees_mod_axis_exchange":rotation_floor,
            "massless_line_moduli":abs(kernel), "rotated_positive_two_stage_counterexample":rotated,
            "branches":points, "observables":OLD["OBS"](target,np.zeros(3),0.)["angle_readings"]}


def analyze(ub, step):
    _,d=OLD["endpoint"](ub,step)
    g,vec=np.linalg.eigh(d["B"].conj().T@d["B"])
    W=d["L"]@(vec@np.diag(1/np.sqrt(g))@vec.conj().T)
    qbasis=np.linalg.qr(W,mode="reduced")[0]
    p0=np.linalg.norm(qbasis[0])**2
    pmin,pmax=(1-np.sqrt(1-np.array([.0803,.0899])))/2
    assert p0<.5
    c_amplitudes=abs(W[0])/np.linalg.norm(W,axis=0)
    p=np.clip(c_amplitudes**2,pmin,pmax)
    overlaps=(c_amplitudes*np.sqrt(p)/p0
              +np.sqrt((1-c_amplitudes**2/p0)*(1-p/p0)))**2
    rmin=np.sqrt(1+.00219/.0000773)
    required=rmin/(1+rmin)
    gap_angle=np.arccos(np.sqrt(max(overlaps)))-np.arccos(np.sqrt(required))
    axis_floor=np.degrees(np.arcsin(np.sin(gap_angle)/np.linalg.cond(W)))
    assert gap_angle>0
    # The expectation identity is basis independent; exercise generic complex P.
    rng=np.random.default_rng(1761)
    identity_error=[]
    for _ in range(9):
        V=np.linalg.qr(rng.normal(size=(2,2))+1j*rng.normal(size=(2,2)))[0]
        f=rng.uniform(.15,2,size=2)*np.exp(1j*rng.normal(size=2))
        X=W@np.diag(f)
        U,sing,_=np.linalg.svd(X@(V@V.T)@X.T)
        lhs=np.linalg.norm(X.conj().T@U[:,0])**2/np.linalg.norm(X,"fro")**2
        identity_error.append(abs(lhs-sing[0]/(sing[0]+sing[1])))
    assert max(identity_error)<1e-13
    interval_check=dict(plane_electron_projection_square=p0,max_column_overlaps_square=overlaps,
                        required_overlap_square=required, necessary_axis_rotation_deg=axis_floor,
                        max_random_projection_identity_error=max(identity_error))
    clock_ratio=math.sqrt(1+.0025/.0000754)
    roots=OLD["diagonal_clock_reactor"](d["L"],clock_ratio,.0851)
    out=[]
    for frozen in roots:
        t=frozen["physical_mu22_over_mu11"]
        z=complex(t["real"],t["imag"])
        mat=d["L"]@np.diag([1,z])@d["L"].T
        out.append(dict(original_t=scalar(z), original_angles=frozen["angles"],
                        accepted=all(frozen["angles_in_box"]) and frozen["SK_T2K_68_90"][-1],
                        inverse=reconstruct(W,mat,d["M"])))
    return {"M_GeV":d["M"],"whole_oscillation_window":interval_check,"diagonal_1758_sources":out}


def tidy(obj):
    if isinstance(obj, dict):
        return {k:tidy(v) for k,v in obj.items()}
    if isinstance(obj,(list,tuple,np.ndarray)):
        return [tidy(v) for v in obj]
    if isinstance(obj,(float,np.floating)):
        return float(f"{obj:.10g}")
    if isinstance(obj,np.bool_):
        return bool(obj)
    return obj


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    out={"primary_inclusive":analyze(OLD["WEAK"]["inclusive_ub"][0],.25),
         "profiled_comparison":analyze(OLD["PROFILE"]["Vub"],.25)}
    dest=HERE/"independent_results.json"
    if args.check:
        assert tidy(out)==json.loads(dest.read_bytes())
        print("1761 independent OK: all Majorana phases and exact source-filter inverse")
    else:
        dest.write_text(json.dumps(tidy(out),indent=2,allow_nan=False)+"\n",encoding="utf-8")
