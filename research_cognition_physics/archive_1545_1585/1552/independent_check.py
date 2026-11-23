"""Muon disappearance at all phases of the archived charged-only completion."""
import argparse
import importlib.util
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.dont_write_bytecode = True


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def charged_inputs():
    old = HERE.parent / "1454"
    spec = importlib.util.spec_from_file_location("round1454_charged",old/"check.py")
    previous = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(previous)
    cfg = read(old/"inputs.json")
    quarks = read(old/cfg["high_quarks"])["high_scale"]
    up,down = (np.array(quarks[k])[:,0] for k in ("up","down"))
    e = np.array(read(old/cfg["high_leptons"])["scales"]["1e16"])[3:,0]
    angles = np.array(quarks["angles_theta12_theta13_theta23_delta"])[:,0]
    U,D = previous.matrices(angles,up,down,cfg["charged_only_completion_relative_phases"])
    root,err = previous.charged_solution(U,D,e,cfg["charged_only_start_x_y"])
    assert err<2e-11
    scale = e[2]/np.linalg.norm(D+root*U,2)
    c, d = scale/down[2], scale*root/up[2]
    U,D = U*up[2],D*down[2]
    E = c*D + d*U
    np.testing.assert_allclose(np.linalg.svd(E,compute_uv=False)[::-1],e,rtol=2e-10)
    return U,D,E,c,d


def poly_coeffs(U,D,E,c,d):
    """adjugate by bilinear row products; no small inverse in the phase polynomial."""
    f = np.array([D, -E])
    adj = np.zeros((3,3,3),complex)
    for p,a in enumerate(f):
        for q,b in enumerate(f):
            for k in range(3):
                adj[p+q,:,k] += np.cross(a[(k+1)%3],b[(k+2)%3])
    det = np.zeros(4,complex)
    for p,a in enumerate(f):
        for q,b in enumerate(adj):
            det[p+q] += a[0]@b[:,0]
    y = np.array([3*D,d*U+c*D-3*E,-c*E])
    tensor = np.zeros((7,3,3),complex)
    for p,a in enumerate(y):
        for q,b in enumerate(adj):
            for r,g in enumerate(y):
                tensor[p+q+r] += a@b@g.T
    # Cross-check the group-independent commuting pencil at two regular points.
    for z in (-1.,1j):
        F = D-z*E
        Y = U + (3+z*c)/(z*d)*F
        K = Y@np.linalg.solve(F,Y.T)
        kernel = sum(t*z**j for j,t in enumerate(tensor))/(z*d)**2/np.linalg.det(F)
        np.testing.assert_allclose(kernel,K,rtol=1e-9,atol=1e-10)
        comm = F@np.linalg.solve(U,K)-K@np.linalg.solve(U,F)
        assert np.linalg.norm(comm,2)<1e-11*np.linalg.norm(K,2)*np.linalg.cond(U)
        mixed = K+(2+3j)*F
        common = D@np.linalg.solve(U,mixed)-mixed@np.linalg.solve(U,D)
        assert np.linalg.norm(common,2)<5e-11*np.linalg.norm(mixed,2)
    return tensor,det


def circle_test():
    U,D,E,c,d = charged_inputs()
    mu = np.linalg.svd(E)[0][:,1]
    t,det = poly_coeffs(U,D,E,c,d)
    norm = np.linalg.norm(t,ord=2,axis=(1,2)).sum()
    t /= norm
    k = np.arange(-3,4)
    tnorm = np.linalg.norm(t,ord=2,axis=(1,2))
    second = np.dot(k*k,tnorm)
    cells = 65536
    h = math.pi/cells
    qmin, observed, max_radius_fraction = 1., 1., 0.
    det_min = 1.
    for segment in np.array_split(np.arange(cells),32):
        angle = 2*math.pi*(segment+.5)/cells
        powers = np.exp(1j*angle[:,None]*k)
        mat = np.einsum("ij,jkl->ikl",powers,t)
        deriv = np.einsum("ij,jkl->ikl",1j*k*powers,t)
        delta = h*np.linalg.norm(deriv,ord=2,axis=(1,2)) + h*h/2*second + 1e-10
        left,s,_ = np.linalg.svd(mat)
        gaps = np.minimum(s[:,0]-s[:,1],s[:,1]-s[:,2])
        assert np.all(gaps>delta)
        loss = delta/(gaps-delta)
        projected = abs(left[:,:,1]@mu.conj())
        angle = np.arccos(np.clip(projected,0.,1.)) + 2*np.arcsin(loss)
        qmin = min(qmin,float(np.min(np.cos(np.minimum(angle,math.pi/2))**2)))
        observed = min(observed,float(np.min(projected**2)))
        max_radius_fraction = max(max_radius_fraction,float(max(loss)))
        ds = np.exp(1j*angle[:,None]*np.arange(4))@det
        det_min = min(det_min,float(min(abs(ds)))-h*sum(np.arange(4)*abs(det)))
    assert det_min>2.45e-10
    return {"cells":cells,"min_abs_det_F_floor":2.45e-10,"min_middle_mu_sample":observed,
            "min_middle_mu_certified":qmin,"max_sine_dilation_radius":max_radius_fraction}


def poisson_upper(count,tail):
    def cdf(v):
        return math.exp(-v)*sum(v**i/math.factorial(i) for i in range(count+1))
    lo,hi=float(count),2.*count+10.
    for _ in range(64):
        mid=(lo+hi)/2
        if cdf(mid)>tail:
            lo=mid
        else:
            hi=mid
    return hi


def counting(q):
    # A budget on half the integrated spectral range of the matter Hamiltonian.
    phase = 295.e3/(1.973269804e-7)*7.64e-14*3*.55/2
    vacuum_angle = 2*math.acos(math.sqrt(q))
    pmin = math.cos(vacuum_angle+phase)**2
    upper = poisson_upper(31,.05)
    usable = 104.-2*14.-2.
    extra = math.acos(math.sqrt(upper/usable))-vacuum_angle-phase
    return {"matter_half_phase":phase,"vacuum_mu_survival_min":math.cos(vacuum_angle)**2,
            "earth_mu_survival_min":pmin,"observed_count":31,"poisson_upper95":upper,
            "required_unoscillated_mu_count":upper/pmin,
            "minimum_retained_fraction_at_76":upper/(pmin*76.),
            "predicted_count_lower_with_two_event_removal":pmin*usable,
            "extra_relative_evolution_angle_needed_rad":extra,
            "extra_mu_middle_ray_angle_if_only_matching_deg":math.degrees(extra)/2}


def compact(x):
    if isinstance(x,dict):
        return {k:compact(v) for k,v in x.items()}
    return float(f"{x:.9g}") if isinstance(x,float) else x


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument("--check",action="store_true")
    mode.add_argument("--write-new",action="store_true")
    args=parser.parse_args()
    circle=circle_test()
    results=compact({"circle":circle,"T2K":counting(circle["min_middle_mu_certified"])})
    path=HERE/"independent_results.json"
    if args.write_new:
        with path.open("x",encoding="utf-8") as stream:
            json.dump(results,stream,indent=2)
            stream.write("\n")
    elif args.check:
        assert results==read(path)
        print("1552 independent: full phase circle, muon survival and counting window verified")
    else:
        print(json.dumps(results,indent=2))
