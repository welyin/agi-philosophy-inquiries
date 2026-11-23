"""Common positive propagation on the two polar-normalized neutral ports."""
import argparse
import json
import math
from pathlib import Path
import runpy

import numpy as np
from numpy.polynomial import Polynomial


HERE = Path(__file__).resolve().parent
ARCH = HERE.parent
CFG = json.loads((HERE/"inputs.json").read_bytes())
PARENT = runpy.run_path(str(ARCH/"1758/check.py"))
ANG = runpy.run_path(str(ARCH/"1731/check.py"))["ANG"]


def frame(data):
    g = data["B"].conj().T@data["B"]
    eig, v = np.linalg.eigh(g)
    invroot = (v/np.sqrt(eig))@v.conj().T
    return data["L"]@invroot


def hierarchy_cap(a,r):
    sv = np.linalg.svd(a@np.diag([r,1.]),compute_uv=False)
    return (sv[0]/sv[1])**2


def root_threshold(a,R):
    # tr(A F^2 A*)/sqrt(det on the two-plane) = r u + v/r.
    gram = a.conj().T@a
    den = math.sqrt(np.linalg.det(gram).real)
    c,d = np.diag(gram).real/den
    roots = np.roots([c,-(math.sqrt(R)+1/math.sqrt(R)),d])
    assert min(roots)<1<max(roots)
    r = max(roots)
    assert abs(hierarchy_cap(a,r)/R-1)<1e-12
    return r


def phase_solutions(a,r,R):
    # Singular frame converts symmetric-unitary congruence to a one-phase
    # circle of allowed heavy eigenlines after imposing the eigenmass ratio.
    U,s,vh = np.linalg.svd(a@np.diag([r,1.]),full_matrices=False)
    upper = (s[0]/s[1])**2
    if upper<R:
        return dict(max_ratio=upper)
    c2 = (R+1/R-2)/(upper+1/upper-2)
    prod = np.prod(s)**2
    h11 = s[0]**2*(s[1]**2+(s[0]**2-s[1]**2)*c2)
    p = (h11-prod/R)/(prod*(R-1/R))
    p = np.clip(p,0,1)
    v,w = U[:,0]*math.sqrt(p),U[:,1]*math.sqrt(1-p)
    off = np.conj(v)*w
    zero = 1-np.sum(abs(U)**2,axis=1)
    constant = abs(v)**2+abs(w)**2
    def prob(phi):
        return abs(v+w*np.exp(1j*phi))**2
    def readings(phi):
        e,m,_ = prob(phi)
        return np.array([4*e*(1-e),1-zero[0]/(1-e),m/(1-e)])
    # Put boundaries for every angle cut on the circle; signs are constant
    # in every interval. Both quadratic branches of the reactor angle retained.
    cuts=[0.,2*math.pi]
    def linear(base,z):
        if abs(z)<1e-14: return
        t=-base/(2*abs(z))
        if abs(t)<=1:
            for angle in (math.acos(t),-math.acos(t)):
                cuts.append((angle-np.angle(z))%(2*math.pi))
    for bound in ANG[0]:
        for e in ((1-math.sqrt(1-bound))/2,(1+math.sqrt(1-bound))/2):
            linear(constant[0]-e,off[0])
    for bound in ANG[1]:
        e=1-zero[0]/(1-bound)
        linear(constant[0]-e,off[0])
    for bound in ANG[2]:
        linear(constant[1]+bound*constant[0]-bound,off[1]+bound*off[0])
    cuts=sorted(set(cuts))
    phases=sorted(cuts+[.5*(l+h) for l,h in zip(cuts,cuts[1:])])
    ok=[]
    for phi in phases:
        vread=readings(phi)
        if np.all((ANG[:,0]-1e-12<=vread)&(vread<=ANG[:,1]+1e-12)):
            ok.append([phi,vread.tolist()])
    # Construct an actual symmetric unitary U' with the requested eigenline.
    # The Takagi singular-frame symmetric unitary has cos^2(chi)=c2.
    phi=ok[len(ok)//2][0] if ok else .37
    c=math.sqrt(c2)
    ss=math.sqrt(1-c2)
    j0=np.array([[c,1j*ss],[1j*ss,c]])
    ph=np.diag([1.,np.exp(1j*(phi-math.pi/2))])
    # The test Hermitian matrix's off-diagonal is -i * positive before rephasing.
    j=vh.conj().T@ph@j0@ph@vh.conj()
    physical=a@np.diag([r,1])@j@np.diag([r,1])@a.T
    left,sv,_=np.linalg.svd(physical)
    np.testing.assert_allclose(abs(left[:,0])**2,prob(phi),atol=5e-10)
    np.testing.assert_allclose(sv[0]/sv[1],R,rtol=1e-12)
    return dict(max_ratio=upper,heaviest_first_singular_frame_weight=p,
                electron_probability_range=[constant[0]-2*abs(off[0]),constant[0]+2*abs(off[0])],
                any_joint_angle_solution=bool(ok),example_angles=readings(phi).tolist() if ok else None)


def overlap_limit(a, Rlow):
    proj_e = np.linalg.norm((a@np.linalg.pinv(a))[0])**2
    # At most the measured two-width upper probability, including both branches.
    # Projection norm excludes the high branch at these physical endpoints.
    large = (1+math.sqrt(1-ANG[0,1]))/2
    assert proj_e < large
    high = (1-math.sqrt(1-ANG[0,1]))/(2*proj_e)
    low = (1-math.sqrt(1-ANG[0,0]))/(2*proj_e)
    r = []
    for col in a.T:
        u = abs(col[0])**2/(proj_e*np.vdot(col,col).real)
        p = np.clip(u,low,high)
        r.append((math.sqrt(u*p)+math.sqrt((1-u)*(1-p)))**2)
    gamma=Rlow/(1+Rlow)
    angle=max(0.,math.acos(math.sqrt(max(r)))-math.acos(math.sqrt(gamma)))
    return dict(ceiling_per_tag=r,required_common_overlap=gamma,
                necessary_tag_rotation_degrees=math.degrees(angle),
                all_positive_weights_impossible=bool(max(r)<Rlow/(1+Rlow)))


def actual_target(d,a,rho):
    inv=np.linalg.pinv(a)
    ab=inv@d["C0"]@inv.T
    col=inv@d["b0"]
    rankone=np.outer(col,col)
    # Scale real rho for better-conditioned roots; retain the rephasing invariant.
    coeff=[Polynomial([ab[i,j],rho*rankone[i,j]]) for i,j in ((0,0),(0,1),(1,1))]
    aa,bb,dd=coeff
    # ad-b^2 is exactly affine because the new source has rank one.
    det0=ab[0,0]*ab[1,1]-ab[0,1]**2
    det1=rho*(ab[0,0]*rankone[1,1]+ab[1,1]*rankone[0,0]-2*ab[0,1]*rankone[0,1])
    reduced=Polynomial([det0,det1])*Polynomial(bb.coef.conj())**2
    def invariant(x):return reduced(x)+abs(bb(x))**4
    phase=float(np.angle(invariant(1),deg=True))
    rr=Polynomial(reduced.coef.imag).roots()
    sturm=runpy.run_path(str(ARCH/"1757/independent_check.py"))["sturm_real_count"]
    assert sum(abs(z.imag)<1e-7 for z in rr)==sturm(reduced.coef.imag)
    records=[]
    for z in rr:
        if abs(z.imag)>1e-7:continue
        k=float(z.real)
        M=ab+k*rankone
        assert min(abs(M[0,0]),abs(M[1,1]),abs(M[1,0]))>1e-12
        allowed=bool(invariant(k).real<0)
        rec=dict(rho=k*rho,negative_invariant=allowed)
        if allowed:
            r=math.sqrt(abs(M[0,0]/M[1,1]))
            F=np.diag([1/r,1])
            check=np.linalg.svd(F@M@F,compute_uv=False)
            assert abs(check[0]/check[1]-1)<1e-7
            obs=PARENT["OBS"](d["C0"],d["b0"],k*rho)
            sv=np.linalg.svd(M,compute_uv=False)
            ls=np.linalg.svd(d["C0"]+k*rho*np.outer(d["b0"],d["b0"]),compute_uv=False)
            rec.update(required_filter_ratio=r,source_hierarchy=sv[0]/sv[1],mass_ratio=ls[0]/ls[1],
                       reactor=obs["angle_readings"][0],all_angles_in_box=all(obs["angles_in_box"]))
        records.append(rec)
    # If any entry vanished on the real affine line, separate special cases
    # could be missed by the generic invariant criterion.
    ratios=[-f.coef[0]/f.coef[1] for f in coeff]
    assert all(abs(z.imag)>1e-7 for z in ratios)
    return dict(registered_source_phase_invariant_degrees=phase,all_real_zero_invariant=records)


def calculate():
    Rlow=math.sqrt(1+.00219/.0000773)
    Rc=math.sqrt(1+.0025/.0000754)
    prof=json.loads((ARCH/"1757/independent_results.json").read_bytes())["combined_central_solution"]
    cases={}
    for name,u,rho in (("inclusive",.0041,None),("oscillation_profile",prof["Vub"],prof["rho"])):
        _,d=PARENT["endpoint"](u,.25)
        rho=d["R"] if rho is None else rho
        a=frame(d)
        m1,m2=d["M"]
        rmin=root_threshold(a,Rlow)
        ratio=m2**2/m1**2
        table={}
        for n,shift_tag,val in ((1,"zero",0.),(2,"zero",0.),(2,"lower",m1),(3,"lower",m1),(3,"upper",m2)):
            r=((m2**2+val**2)/(m1**2+val**2))**n
            table[f"{n}_{shift_tag}"]=dict(raw_ratio=r,central_frequency=phase_solutions(a,r,Rc))
        bounds={}
        for n in range(1,5):
            step=rmin**(1/n)
            cap=(m2**2-step*m1**2)/(step-1)
            bounds[str(n)]=math.sqrt(cap) if cap>0 else None
        cases[name]=dict(masses_GeV=[m1,m2],phase_reactor_solar_atmos_boxes=ANG.tolist(),
                         mass_only_filter_ratio_floor=rmin,largest_common_gap_GeV=bounds,known_spectra=table,
                         electron_overlap=overlap_limit(a,Rlow),full_light_matrix=actual_target(d,a,rho))
    return cases


def tidy(value):
    if isinstance(value,dict):return {k:tidy(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [tidy(v) for v in value]
    if isinstance(value,(float,np.floating)):return float(f"{float(value):.10g}")
    return value

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    args=p.parse_args()
    result=tidy(calculate())
    file=HERE/"results.json"
    if args.check:
        assert result==json.loads(file.read_bytes())
        print("1761: positive common propagation, arbitrary symmetric-unitary seed and phase witnesses verified")
    else:
        if (HERE/"acceptance.json").exists() and file.exists():
            raise SystemExit("Use --check to verify an archived result.")
        with file.open("w",encoding="utf8") as out:
            json.dump(result,out,indent=2,allow_nan=False)
            out.write("\n")
