"""Independent quark-side test of two symmetric SO(10) flavour tensors."""
import argparse
import json
import itertools
import math
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
def read(p):
    return json.loads(p.read_text(encoding="utf-8"))


def ckm(a):
    t12,t13,t23,delta=a
    s12,s13,s23=np.sin(a[:3]); c12,c13,c23=np.cos(a[:3])
    z=np.exp(1j*delta)
    return np.array([[c12*c13,s12*c13,s13/z],
        [-s12*c23-c12*s23*s13*z,c12*c23-s12*s23*s13*z,s23*c13],
        [s12*s23-c12*c23*s13*z,-c12*s23-s12*c23*s13*z,c23*c13]])


def worst_case(up,down,leptons,mix,widths):
    # mix contains sin(theta13) and sin(theta23) intervals.
    lu=(leptons[1,0]-widths*leptons[1,1])/(leptons[2,0]+widths*leptons[2,1])
    ru=(up[1,0]+widths*up[1,1])/(up[2,0]-widths*up[2,1])
    rd=(down[1,0]+widths*down[1,1])/(down[2,0]-widths*down[2,1])
    eta=1-(1-mix[0]**2)*(1-mix[1]**2)
    delta=abs(rd-ru)/2
    def perturb(g):
        return ((rd+ru)*g/2+delta)/(1-delta)
    def test(g):
        eps=perturb(g)
        return max(1-eps,0)*max(lu-eps,0)-eta*g*g/4
    lo,hi=0,20
    for _ in range(80):
        mid=(hi+lo)/2
        if test(mid)>0: lo=mid
        else: hi=mid
    ell_sq=1+sum(((leptons[i,0]+widths*leptons[i,1])/
                    (leptons[2,0]-widths*leptons[2,1]))**2 for i in [0,1])
    return {"sum_normalized_gain_floor":lo, "signed_Frobenius_interference_ceiling":2*ell_sq/lo**2-1,
        "each_normalized_gain_floor":(lo-1-perturb(lo))/2, "CKM_tail_squared_upper":eta}


def matrices(angles,up,down,p):
    u=up/up[-1]*np.exp(1j*np.array([p[0],p[1],0]))
    v=ckm(angles)
    d=down/down[-1]*np.exp(1j*np.array([p[2],p[3],0]))
    return np.diag(u), (v*d)@v.T


def charged_solution(u,d,leptons,start):
    def residual(z):
        sig=np.linalg.svd(d+(z[0]+1j*z[1])*u,compute_uv=False)[::-1]
        return np.log((sig[:2]/sig[-1])/(leptons[:2]/leptons[-1]))
    z=np.array(start,dtype=float)
    for _ in range(120):
        f=residual(z)
        if max(abs(f))<5e-12: break
        step=1e-6
        jac=np.column_stack([(residual(z+e*step)-residual(z-e*step))/(2*step) for e in np.eye(2)])
        dz=np.linalg.lstsq(jac,-f,rcond=None)[0]
        for line in range(18):
            trial=z+dz/(2**line)
            if np.linalg.norm(residual(trial))<np.linalg.norm(f):
                z=trial
                break
        else: break
    err=max(abs(residual(z)))
    return complex(*z),float(err)


def type_ii(d,e):
    ve,se,_=np.linalg.svd(e)
    # electron is the last column returned in descending SVD.
    v=ve[:,::-1]
    # Normalize the squared-mass Hamiltonian to have finite absolute Lipschitz error.
    s=np.linalg.norm(d,2)+np.linalg.norm(e,2)
    db=d/s; eb=e/s
    n=2048
    grid=2*np.pi*np.arange(n)/n
    nu=db[None,:,:]-np.exp(1j*grid)[:,None,None]*eb[None,:,:]
    h=nu@np.transpose(nu.conj(),(0,2,1))
    lam,axes=np.linalg.eigh(h)
    gap=lam[:,2]-lam[:,0]
    rho=(lam[:,1]-lam[:,0])/gap
    evec=v[:,0]
    he=h@evec
    orth=he-(he@evec.conj())[:,None]*evec[None,:]
    # H = DD* + EE* - exp(+i phi) ED* - exp(-i phi) DE*.
    eta=2*np.linalg.norm(db@eb.conj().T,2)*np.pi/n
    var=np.linalg.norm(orth,axis=1)
    cert=np.max((var+eta)/(gap-2*eta)) if np.all(gap>2*eta) else None
    no_min=np.min(np.maximum(0,(lam[:,1]-lam[:,0])-2*eta)/(gap+2*eta))
    io_min=np.min(np.maximum(0,(lam[:,2]-lam[:,1])-2*eta)/(gap+2*eta))
    out={"grid_cells":n,"operator_norm_radius":float(eta),"variance_over_span_certified_upper":None if cert is None else float(cert),
        "NO_gap_ratio_certified_floor":float(no_min),"IO_gap_ratio_certified_floor":float(io_min),
        "grid_variance_over_span_range":[float(np.min(var/gap)),float(np.max(var/gap))]}
    return out

def evaluate():
    cfg=read(HERE/"inputs.json")
    hi=read(HERE/cfg["high_quarks"])["high_scale"]
    qhi=[np.array(hi[k]) for k in ("up","down")]
    ehi=np.array(read(HERE/cfg["high_leptons"])["scales"]["1e16"])[3:]
    lo=read(HERE/cfg["low_spectra"])
    qlo=[np.array(lo["yukawa"][k]) for k in ("up","down")]
    elo=np.array(lo["yukawa"]["charged_leptons"])
    wid=cfg["working_mass_halfwidths"]
    hiangs=np.array(hi["angles_theta12_theta13_theta23_delta"])
    lowangs=np.array(lo["quark_sines_and_phase"])
    ret={}
    for name,quarks,leptons,angle,is_angle in [
        ("MZ",qlo,elo,lowangs,False),("1e16",qhi,ehi,hiangs,True)]:
        edge=angle[1:3,0]+wid*angle[1:3,1]
        ret[name]=worst_case(*quarks,leptons,np.sin(edge) if is_angle else edge,wid)
    up,down=(q[:,0] for q in qhi)
    charged=ehi[:,0]
    u,d=matrices(hiangs[:,0],up,down,cfg["charged_only_completion_relative_phases"])
    z,err=charged_solution(u,d,charged,cfg["charged_only_start_x_y"])
    assert err<2e-11
    tau=charged[2]; top=up[2]; btm=down[2]
    x=1/np.linalg.norm(d+z*u,2)
    yp=x*z
    e=x*(d+z*u)
    # UNIT: down and up are now divided by bottom and top, E by tau.
    c=x*tau/btm; gain_u=yp*tau/top
    jquark=np.imag(np.prod([ckm(hiangs[:,0])[0,0],ckm(hiangs[:,0])[1,1],
        ckm(hiangs[:,0])[0,1].conjugate(),ckm(hiangs[:,0])[1,0].conjugate()]))
    interference=np.real(np.vdot(x*d,yp*u))/(x*abs(yp))
    # Two sources can be reconstructed at each phase. Illustrate phi=0.
    yd=btm*d; ye=tau*e; yu=top*u
    H=(3*yd+ye)/4; F=(yd-ye)/4
    r,s=(1-c)/gain_u,(-3-c)/gain_u
    resid=np.linalg.norm(yu-r*H-s*F,2)
    assert resid<1e-12
    model=type_ii(yd,ye)
    criterion=(1-cfg["oscillation"]["smaller_gap_full_span_upper"])*math.sqrt(cfg["oscillation"]["sin2_2theta13_min"])/2
    assert model["variance_over_span_certified_upper"]<criterion
    ret["charged_only_completion"]={"real_imag_z":[z.real,z.imag],"largest_relative_log_mass_error":err,
        "actual_charged_yukawas":np.linalg.svd(ye,compute_uv=False)[::-1].tolist(),
        "quark_J":float(jquark),"dimensionless_down_and_up_gains":[x,abs(yp)],"normalized_trace_interference":float(interference),
        "coeff_c_d_at_phi0":[[c,0.],[gain_u.real,gain_u.imag]],"two_source_reconstruction_operator_error":float(resid)}
    ret["all_phi_of_this_completion"]=model
    ret["oscillation_variance_floor_either_order"]=criterion
    ret["minimum_joint_transfer_for_this_completion"]=criterion-model["variance_over_span_certified_upper"]
    def clean(a):
        if isinstance(a,dict): return {k:clean(v) for k,v in a.items()}
        if isinstance(a,list): return [clean(v) for v in a]
        if isinstance(a,(float,np.floating)):return round(float(a),10)
        return a
    return clean(ret)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    g=p.add_mutually_exclusive_group()
    g.add_argument("--write",action="store_true");g.add_argument("--check",action="store_true")
    args=p.parse_args()
    results=json.dumps(evaluate(),indent=2,allow_nan=False)+"\n"
    path=HERE/"results.json"
    if args.write:
        with path.open("x",encoding="utf-8") as f: f.write(results)
    elif results!=path.read_text(encoding="utf-8"):
        raise SystemExit("1454 mismatch")
    else: print("1454 OK")

