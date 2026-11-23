"""Can a shared complex charged transition carry one LNV phase reference to the other pole?"""
import argparse
import json
import math
from fractions import Fraction
from pathlib import Path
import runpy

import numpy as np
from numpy.polynomial import Polynomial


HERE = Path(__file__).resolve().parent
ARCH = HERE.parent
ONE = runpy.run_path(str(ARCH/"1761/check.py"))
PARENT = ONE["PARENT"]
ANG = ONE["ANG"]
OBS = PARENT["OBS"]
RC = math.sqrt(1+.00250/.0000754)
RBOX = [math.sqrt(1+.00219/.0000773),math.sqrt(1+.00263/.0000736)]


def positive_count(polynomial):
    poly=[Fraction.from_float(float(c)) for c in polynomial.trim().coef]
    seq=[poly,[j*poly[j] for j in range(1,len(poly))]]
    assert poly[0]
    while len(seq[-1])>1:
        rem=seq[-2].copy()
        div=seq[-1]
        while len(rem)>=len(div):
            h=rem[-1]/div[-1]
            for j,v in enumerate(div):
                rem[len(rem)-len(div)+j]-=h*v
            while rem and not rem[-1]:
                rem.pop()
        assert rem
        # Positive rescaling keeps the Sturm signs and limits Fraction growth.
        seq.append([-v/abs(rem[-1]) for v in rem])
    def changes(seq):
        signs=[v>0 for v in seq if v]
        return sum(i!=j for i,j in zip(signs,signs[1:]))
    return changes([s[0] for s in seq])-changes([s[-1] for s in seq])


def profile(a, axis, ref, eta):
    transfer = ref[1,0]/abs(ref[1,0])
    seed = np.diag([1.,eta*transfer**2])
    _,v = np.linalg.eigh(axis)
    projection = np.outer(v[:,-1],v[:,-1].conj())
    other = np.eye(2)-projection
    # F = (1-P) + rP. Remove the common scale; either orientation can be enhanced.
    mats = [a@other@seed@other.T@a.T,
            a@(other@seed@projection.T+projection@seed@other.T)@a.T,
            a@projection@seed@projection.T@a.T]
    # Check the same physical source after giving both original Dirac tags
    # different phases, including transport of the reference on the first tag.
    d=np.diag(np.exp(1j*np.array([.37,-.64])))
    transformed=d.conj().T@ref@d
    t2=transformed[1,0]/abs(transformed[1,0])
    j2=np.conj(d[0,0])**2*np.diag([1.,eta*t2**2])
    f2=d.conj().T@(other+1.72*projection)@d
    np.testing.assert_allclose(a@d@f2@j2@f2.T@d.T@a.T,
                               sum(m*1.72**j for j,m in enumerate(mats)),atol=1e-12)
    mass_product = np.linalg.det(a.conj().T@a).real
    trace = Polynomial([sum(np.vdot(mats[k],mats[j-k]).real for k in range(3)
                         if 0<=j-k<3) for j in range(5)])
    r2=Polynomial([0,0,1])
    prod=mass_product*r2

    def real_roots(poly):
        roots=poly.roots()
        out=[float(t.real) for t in roots if abs(t.imag)<1e-6 and t.real>0]
        assert len(out)==positive_count(poly), (out,poly.coef)
        return sorted(out)
    def from_ratio(R):
        return real_roots(trace-(R+1/R)*prod)
    def observation(r):
        mass=sum(m*r**j for j,m in enumerate(mats))
        return OBS(mass,np.zeros(3),0.)
    def describe(r):
        obs=observation(r)
        return dict(fast_over_other=r,angles=obs["angle_readings"],cp_deg=obs["delta_deg"],
                    SK_T2K_68_90=obs["SK_T2K_FC68_90"],all_angles_in_window=all(obs["angles_in_box"]))
    # Electron's projection onto the neutral plane is fixed. For a target
    # heavy-line probability p solve (2E-mT)^2=(2p-m)^2(T^2-4P).
    e=Polynomial([sum(np.vdot(mats[k][0],mats[j-k][0]).real for k in range(3)
                         if 0<=j-k<3) for j in range(5)])
    pe=(a@np.linalg.pinv(a))[0,0].real
    eqs = [trace-(R+1/R)*prod for R in RBOX]
    def angle_poly(x,b,p):
        return (2*x-b*trace)**2-(2*p-b)**2*(trace**2-4*prod**2)
    for alpha in ANG[0]:
        for z in [(1-math.sqrt(1-alpha))/2,(1+math.sqrt(1-alpha))/2]:
            eqs.append(angle_poly(e,pe,z))
    for solar in ANG[1]:
        z=1-(1-pe)/(1-solar)
        eqs.append(angle_poly(e,pe,z))
    # P_mu3 + b P_e3 = b. This is another matrix element of a diagonal
    # weight with calibrated trace on the two-plane.
    emu=Polynomial([sum(np.vdot(mats[k][1],mats[j-k][1]).real for k in range(3)
                            if 0<=j-k<3) for j in range(5)])
    pmu=(a@np.linalg.pinv(a))[1,1].real
    for atm in ANG[2]:
        eqs.append(angle_poly(emu+atm*e,pmu+atm*pe,atm))
    knots=sorted({r for eq in eqs for r in real_roots(eq)})
    inner=[math.sqrt(lo*hi) for lo,hi in zip(knots,knots[1:])]
    tests=sorted([.5*knots[0],2*knots[-1],*inner,*knots])
    def valid(r):
        ratio=(trace(r)/prod(r))
        if ratio < RBOX[0]+1/RBOX[0]-1e-10 or ratio>RBOX[1]+1/RBOX[1]+1e-10:
            return False
        ang=np.array(observation(r)["angle_readings"])
        return bool(np.all((ANG[:,0]-1e-11<=ang)&(ang<=ANG[:,1]+1e-11)))
    witnesses=[r for r in tests if valid(r)]
    # Four positive poles per leg with the independently named upper mass.
    r4=((1+(7171.077576/5129.332691)**-2)/2)**4
    return dict(source_transfer_phase_deg=float(np.angle(transfer,deg=True)),relative_parity=eta,
          central_mass_ratio_roots=[describe(r) for r in from_ratio(RC)],
          frequency_and_angle_box_witness=(describe(witnesses[len(witnesses)//2]) if witnesses else None),
          positive_cut_counts=[positive_count(q) for q in eqs],
          four_step_known_scale_examples=[describe(r4),describe(1/r4)])


def calculation():
    # For F use the independently distinguishable electron; for the transported
    # common phase compare transfer by that same tag vs all charged states.
    weights=json.loads((ARCH/"1386/inputs.json").read_bytes())["yukawa"]["charged_leptons"]
    weights=np.array(weights)[:,0]**2
    source=json.loads((ARCH/"1757/independent_results.json").read_bytes())["combined_central_solution"]
    cases={}
    for name,u in (("semileptonic",.0041),("clocks",source["Vub"])):
        color=runpy.run_path(str(HERE/"same_color_axis.py"))
        d,bg,m,dr,U,a,target,chosen=color["endpoint"](u,.25)
        gate=U.conj().T@(dr*[bg["cW"],1.,bg["fW"]])@dr.conj().T@U
        down=U.conj().T@np.outer(dr[:,0],dr[:,0].conj())@U
        e=np.outer(a[0].conj(),a[0])
        full=a.conj().T@(weights[:,None]*a)
        cases[name]={}
        for tag,axis,h in (("same_e",e,e),("charged_spectrum",e,full),
                           ("electron_axis_color_phase",e,gate),("color_axis_same_phase",gate,gate),
                           ("color_axis_e_phase",gate,e),("down_axis_color_phase",down,gate),
                           ("down_axis_e_phase",down,e)):
            for parity in [-1,1]:
                cases[name][tag+str(parity)]=profile(a,axis,h,parity)
    return cases


def tidy(x):
    if isinstance(x,dict):return {k:tidy(v) for k,v in x.items()}
    if isinstance(x,(list,tuple,np.ndarray)):return [tidy(v) for v in x]
    if isinstance(x,(bool,np.bool_)):return bool(x)
    if isinstance(x,(float,np.floating)):return float(f"{float(x):.10g}")
    return x


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    opt=p.parse_args()
    result=tidy(calculation())
    target=HERE/"transported_phase_results.json"
    if opt.check:
        assert result==json.loads(target.read_bytes())
        print("1762 common LNV phase transported by measured charged tags: both parities verified")
    else:
        if (HERE/"acceptance.json").exists(): raise SystemExit("Use --check for an archived round.")
        target.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf8")
