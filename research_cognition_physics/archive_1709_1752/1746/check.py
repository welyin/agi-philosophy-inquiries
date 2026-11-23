"""Invert a common lepton contraction using only masses, a mixing matrix and its spectrum."""
import argparse
import json
import math
from pathlib import Path
import runpy

import numpy as np

HERE = Path(__file__).resolve().parent

def read(path):
    return json.loads((HERE/path).read_bytes())

CFG=read("inputs.json")


def takagi(matrix):
    u, s, _ = np.linalg.svd(matrix)
    d = np.diag(u.conj().T @ matrix @ u.conj())
    u[:,:2] *= np.exp(.5j*np.angle(d[:2]))
    assert np.linalg.norm((u*s)@u.T-matrix) < 2e-12*np.linalg.norm(matrix)
    assert s[0]>s[1] and s[1]>s[2]*1e8
    return u,s


def charged_frame(y,n):
    basis, ell, _ = np.linalg.svd(y)
    basis=basis[:,::-1]
    return ell[::-1], basis.T@n@basis, basis


def weights(ell,n):
    q = n / np.outer(ell,ell)
    return takagi(q)


def visible(n):
    v,m=takagi(n)
    v=v[:,[1,0]]
    assert min(abs(v[0])) > 1e-8
    return v*(v[0].conj()/abs(v[0])), m[0]/m[1]


def phases(ell,v,ratio,k):
    """Use exterior area to avoid subtracting nearly equal Gram entries."""
    z=v/ell[:,None]
    w=z.conj().T@z
    area=sum(abs(np.linalg.det(z[[i,j],:]))**2 for i,j in ((0,1),(0,2),(1,2)))
    mean=w[0,0].real**2+(ratio*w[1,1].real)**2
    width=2*ratio*abs(w[0,1])**2
    target=(k+1/k)*ratio*area
    if width < 1e-15*mean:
        raise ValueError("The registered nonorthogonal weighted mass columns are required.")
    q=(target-mean)/width
    if abs(q)>1+2e-10:
        return q,[]
    a=math.acos(np.clip(q,-1,1))
    return q,[p-2*np.angle(w[0,1]) for p in (-a,a)]


def pencil(b0,b1,ell0,ell1,fac,sign):
    part=np.diag([1,sign,0])
    null=np.diag([0,0,1])
    def direct(which):
        return fac*((b0.conj()@which@b1.T)*ell1[None,:])/ell0[:,None]
    def inverse(which):
        return ((b1.conj()@which@b0.T)*ell0[None,:])/(fac*ell1[:,None])
    return direct(part),direct(null),inverse(part),inverse(null)


def circle_rows(t0,t1,s0,s1,goal):
    """Represent trace and inverse trace as affine functions of (cos chi, sin chi)."""
    a,b = np.vdot(t0,t1),np.vdot(s0,s1)
    coeff=2*np.array([[a.real,-a.imag],[b.real,b.imag]])
    centers=np.array([np.linalg.norm(t0)**2+np.linalg.norm(t1)**2,
                      np.linalg.norm(s0)**2+np.linalg.norm(s1)**2])
    norms=np.linalg.norm(coeff,axis=1)
    j=np.argmax(norms)
    if norms[j]==0:
        return ([0.] if np.max(abs(goal-centers))<3e-9 else []), centers,coeff
    unit=coeff[j]/norms[j]
    q=(goal[j]-centers[j])/norms[j]
    if abs(q)>1+2e-9:
        return [],centers,coeff
    a=math.atan2(unit[1],unit[0])
    w=math.acos(np.clip(q,-1,1))
    return [a+w,a-w],centers,coeff


def inverse_case(e0,n0,e1,v,ratio,legs,recover=False):
    b0,s0=weights(e0,n0)
    cosine,ph=phases(e1,v,ratio,s0[0]/s0[1])
    goal=np.array([sum(legs**2),sum(legs**-2)])
    fac=(np.prod(legs)*np.prod(e0)/np.prod(e1))**(1/3)
    rows=[]
    for phi in ph:
        target=(v*np.array([1,ratio*np.exp(1j*phi)]))@v.T
        b1,s1=weights(e1,target)
        np.testing.assert_allclose(s0[0]/s0[1],s1[0]/s1[1],rtol=4e-10)
        for sign in (-1,1):
            p,q,ip,iq=pencil(b0,b1,e0,e1,fac,sign)
            angles,centers,harmonics=circle_rows(p,q,ip,iq,goal)
            # Exact lower bound on the first trace for all choices of the null-line phase.
            a=np.vdot(p,q)
            low_trace=float(np.linalg.norm(p-q*np.exp(-1j*np.angle(a)))**2)
            trials=[]
            for chi in angles:
                T=p+q*np.exp(1j*chi)
                left,sv,_=np.linalg.svd(T)
                uv=T.T@n0@T
                mass_gain=np.vdot(target,uv)/np.vdot(target,target)
                np.testing.assert_allclose(uv,mass_gain*target,rtol=2e-8,atol=1e-10*abs(mass_gain))
                ce=T.conj().T @ np.diag(e0**2) @ T
                np.testing.assert_allclose(ce,np.diag((fac*e1)**2),atol=2e-13*max(e1)**2*fac**2)
                diff=float(max(abs(sv-np.sort(legs)[::-1])))
                trial=dict(chi_rad=float(chi),max_singular_error=diff,match=diff<1e-7)
                if recover and trial["match"]:
                    filter_src=(left*sv)@left.conj().T
                    inv=np.linalg.inv(T)
                    trial.update(filter_src=filter_src,heavy_gram=inv@inv.conj().T-np.eye(3))
                trials.append(trial)
            rows.append(dict(phi_rad=float(phi),sign=sign,min_trace_minus_required=low_trace-goal[0],
                             harmonic_matrix=harmonics.tolist(),center_minus_required=(centers-goal).tolist(),
                             trials=trials))
    return dict(majorana_cosine=float(cosine),phase_sign_branches=rows,
                any_solution_at_tolerance=any(t["match"] for r in rows for t in r["trials"]))


def pmns(s12,s13,s23,delta):
    a,b,c=np.sqrt([s12,s13,s23])
    A,B,C=np.sqrt(1-np.array([s12,s13,s23]))
    z=np.exp(1j*delta)
    u=np.array([[A*B,a*B,b*z.conjugate()],
                [-a*C-A*c*b*z,A*C-a*c*b*z,c*B],
                [a*c-A*C*b*z,-A*c-a*C*b*z,C*B]])
    out=u[:,[1,2]].conj()
    return out*(out[0].conj()/abs(out[0]))


def solve():
    old=runpy.run_path(str(HERE/CFG["physical_mother"]))
    b=old["mother"](CFG["step"])
    pt=b["pt"]
    tr=np.linalg.inv(b["S"]@pt["low_E"])
    n=tr.T@(pt["C"]+b["nu"]["FROZEN"][0]*np.outer(pt["b"],pt["b"]))@tr
    y=b["yW"][2,:,:3]
    e0,n0,_=charged_frame(y,n)
    legs=np.array([b["cW"],1.,b["fW"]])
    cases={}
    # Controls include the actual named material gate and an independently oriented complex gate.
    rng=np.random.default_rng(1746)
    z=np.linalg.qr(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))[0]
    for name,f in (("named_hidden_gate",b["F"]),("complex_hidden_gate",(z*legs)@z.conj().T)):
        e1,n1,_=charged_frame(f@y,f.T@n@f)
        v,r=visible(n1)
        cases[name]=inverse_case(e0,n0,e1,v,r,legs)
        assert cases[name]["any_solution_at_tolerance"]
    fre=read("../1668/inputs.json")
    measured=read("../1484/inputs.json")
    atm=read("../1659/inputs.json")["atmospheric_NOvA_upper_octant"]
    s13=(1-math.sqrt(1-measured["DayaBay_sin2_2theta13"][0]))/2
    angle=[measured["JUNO_sin2_theta12_NO"][0],s13,atm["central"]]
    # These central/edge clocks and mixing patterns are *threshold benchmarks*.
    # The independent calculation instead compares genuine low-scale outputs to data.
    rcenter=math.sqrt(fre["juno_NO_dm31"][0]/fre["juno_NO_dm21"][0])
    redge=math.sqrt((fre["juno_NO_dm31"][0]-2*fre["juno_NO_dm31"][1])/
                    (fre["juno_NO_dm21"][0]+2*fre["juno_NO_dm21"][2]))
    vm,rm=visible(n0)
    for name,v in (("mother_columns",vm),("registered_central_angles",pmns(*angle, math.radians(CFG["probe_dirac_phase_deg"])))):
        for rname,r in (("central_ratio",rcenter),("two_window_edge",redge)):
            cases[name+"_"+rname]=inverse_case(e0,n0,e0,v,r,legs)
    return dict(pole_GeV=b["Mweak"], spectrum=legs.tolist(),determinant=float(np.prod(legs)),
                mother_mass_ratio=rm,central_reference_ratio=rcenter,two_window_reference_ratio=redge,
                cases=cases)


def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:
            compare(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for aa,bb in zip(a,b):
            compare(aa,bb)
    else:
        np.testing.assert_allclose(a,b,rtol=3e-8,atol=5e-8)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=solve()
    out=HERE/"results.json"
    if args.check:
        compare(result,json.loads(out.read_bytes()))
        print("1746 inverse common filter, Majorana branches and null circle verified")
    else:
        with out.open("x",encoding="utf8") as f:
            json.dump(result,f,indent=2,allow_nan=False)
            f.write("\n")
