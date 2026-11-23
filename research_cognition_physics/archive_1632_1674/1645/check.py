"""Independent CP-even rotation and spectrum calibration followed by a signed neutral-B holdout."""
import argparse
import json
import math
import runpy
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent
def read(name):
    return json.loads((HERE/name).read_bytes())


def bounds(row,n=2):
    return row[0]-n*sum(row[1:]),row[0]+n*sum(row[1:])


def checkin(v,box):
    assert box[0]<=v<=box[1], (v,box)


def pack(z):
    a=np.asarray(z)
    return np.stack((a.real,a.imag),axis=-1).tolist()


def calibrate(cfg,locked):
    t=read(cfg["tree"]); semi=read(cfg["semileptonic"]);inc=read(cfg["inclusive"])
    masses=read(cfg["spectra"])["yukawa"]["down"]
    yd,ys=map(bounds,masses[:2])
    dl,dh=math.sqrt(yd[0]/ys[1]),math.sqrt(yd[1]/ys[0])
    ul,uh=bounds(t["ub_over_cb"]);wl,wh=bounds(semi["BelleII_Vcb"])
    nk,f=bounds(t["K_product"]),bounds(t["K_form"])
    kl,kh=nk[0]/f[1],nk[1]/f[0]
    nulo,nuhi=bounds(t["ud_beta"])
    required_ub=math.sqrt(1-nuhi**2-kh**2)
    # Exact inverse for the upper K boundary on the right-angle fixture.
    x=kh**2
    def largest_u(b,w):
        T=w*w
        disc=1-T-x*(1+b)
        numer=x*(1+b)-b
        return math.sqrt(2*numer/(disc+math.sqrt(disc*disc-4*T*numer)))
    def inferred_d(u,w):
        a=u*u
        return math.sqrt((x*(1+a)-a*(1-w*w*(1+a)))/(1-x*(1+a)))
    if locked:
        uh=min(uh,largest_u(dl*dl,wh))
        d_at_min=min(dh,inferred_d(ul,wl))
        assert dl<=d_at_min<=dh
        tlo=ul/(d_at_min*math.sqrt(1-wl**2*(1+ul**2)))
        thi=uh/(dl*math.sqrt(1-wh**2*(1+uh**2)))
        cos_bounds=[0.,0.]
    else:
        ilo,ihi=bounds(inc["inclusive_holdout"]["Vcb"])
        wl,wh=max(wl,ilo),min(wh,ihi)
        ul=max(ul,required_ub/wh)
        # Nuclear beta and maximum Ub also imply a tighter lower kaon modulus.
        kl=max(kl,math.sqrt(1-nuhi**2-(uh*wh)**2))
        cl=math.sqrt(1-wh**2*(1+uh**2));ch=math.sqrt(1-wl**2*(1+ul**2))
        numer_lo=ul**2*cl**2+dl**2-kh**2*(1+uh**2)*(1+dh**2)
        numer_hi=uh**2*ch**2+dh**2-kl**2*(1+ul**2)*(1+dl**2)
        lo=numer_lo/(2*uh*dh*ch);hi=numer_hi/(2*ul*dl*cl)
        assert 0<lo<hi<1
        cos_bounds=[lo,hi]
        tlo=ul*math.sqrt(1-hi*hi)/(dh*ch-ul*lo)
        thi=uh*math.sqrt(1-lo*lo)/(dl*cl-uh*hi)
    return {"u":[ul,uh],"d":[dl,dh],"w":[wl,wh],"Kus":[kl,kh],
            "cos_phi":cos_bounds,"beta_tangent":[tlo,thi],
            "ub_required_by_nuclear":required_ub,"ub_max":uh*wh,
            "positive_b_over_mu_a_min":cos_bounds[0]/math.sqrt(1-cos_bounds[0]**2)}


def short_kernel(cfg):
    h=read(cfg["decay"]);w=read(cfg["kinematics"])
    masses=h["up_charm_mass_GeV"]+[w["short_distance"]["mtMS"][0]]
    x=(np.array(masses)/w["const_GeV"]["MW"])**2
    one=runpy.run_path(str(HERE.parent/"1644"/"check.py"))["short_distance"]
    return np.array([[one(i,j) for j in x] for i in x])


def holdout_bounds(domain,loop):
    ul,uh=domain["u"];dl,dh=domain["d"];wl,wh=domain["w"];kl,kh=domain["Kus"]
    cmin=math.sqrt(1-wh*wh*(1+uh*uh))
    udmin=math.sqrt(1-kh*kh-(uh*wh)**2)
    # The phase of the favored B tree and the neutral-D tag compared with the charm axis.
    tdts=wh**2*(1+uh**2)*dh/(1+dl**2)
    tag=2*math.asin(tdts/(udmin*kl))
    # Divide the GIM box by the top term and bound both conjugate branches at once.
    lu=uh/(dl/math.sqrt(1+dl**2)*cmin)
    r=[lu,1+lu,1.]
    gsum=sum(abs(loop[i,j])*r[i]*r[j] for i in range(3) for j in range(3) if (i,j)!=(2,2))
    gcap=gsum/abs(loop[2,2])
    half=tag+math.asin(gcap)
    tlo,thi=domain["beta_tangent"]
    theta=[2*math.atan(tlo)-half,2*math.atan(thi)+half]
    assert 0<theta[0]<theta[1]<math.pi/2
    bounds_s=[math.sin(x) for x in theta]
    return {"beta_degrees":[math.degrees(math.atan(t)) for t in (tlo,thi)],
            "tag_phase_cap_rad":tag,"GIM_relative_cap":gcap,"S_reduced_positive_outer":bounds_s,
            "S_reduced_negative_outer":[-bounds_s[1],-bounds_s[0]]}


def evaluate(cfg,point,positive,locked,kernel):
    inp=read(cfg["fx"]);y=read(cfg["spectra"])["yukawa"]
    fx=runpy.run_path(str(HERE.parent/"1412"/"check.py"))
    u=point["up_tangent"];w=point["cb"];yd=point["yd"];ys=point["ys"]
    d=math.sqrt(yd/ys)
    s=w*math.sqrt(1+u*u);c=math.sqrt(1-s*s)
    cosphi=0. if locked else (u*u*c*c+d*d-point["us"]**2*(1+u*u)*(1+d*d))/(2*u*d*c)
    assert abs(cosphi)<1
    phi=(1 if positive else -1)*math.acos(cosphi)
    ru,rd=fx["rotation"](u),fx["rotation"](d)
    r23=np.array([[1.,0.,0.],[0.,c,s],[0.,-s,c]])
    ou,od=ru,r23@rd
    real_u=ou@np.diag([-y["up"][0][0],y["up"][1][0],y["up"][2][0]])@ou.T
    real_d=od@np.diag([-yd,ys,y["down"][2][0]])@od.T
    vu=inp["vev_example"];xu,xd,aa=vu["xu_GeV"],vu["xd_GeV"],vu["abs_A_GeV"]
    b=xu*xd*cosphi
    mu=xu*xd*abs(math.sin(phi))/aa
    avev=-math.copysign(aa,phi)
    np.testing.assert_allclose(xu*xd*np.exp(-1j*phi),b+1j*mu*avev,atol=1e-11)
    norm=vu["v_GeV"]/math.sqrt(2)
    w0=math.sqrt(norm*norm-vu["wu_GeV"]**2-vu["wd_GeV"]**2)
    hu=real_u/np.array([vu["wu_GeV"],w0,w0])[:,None]*norm
    hd=real_d/np.array([vu["wd_GeV"],w0,w0])[:,None]*norm
    muq=np.diag([vu["wu_GeV"],w0,w0])@hu/norm
    mdq=np.diag([vu["wd_GeV"]*np.exp(-1j*phi),w0,w0])@hd/norm
    assert abs(real_d[0,0])<1e-15
    U,su,VRu=np.linalg.svd(muq);D,sd,VRd=np.linalg.svd(mdq)
    U=U[:,np.argsort(su)];D=D[:,np.argsort(sd)]
    K=U.conj().T@D
    np.testing.assert_allclose(np.sort(su),[e[0] for e in y["up"]],atol=1e-14)
    np.testing.assert_allclose(np.sort(sd),[yd,ys,y["down"][2][0]],atol=1e-13)
    ph=np.diag([np.exp(-1j*phi),1.,1.])
    np.testing.assert_allclose(abs(K),abs(ou.T@ph@od),atol=1e-12)
    np.testing.assert_allclose(K.conj().T@K,np.eye(3),atol=1e-14)
    # A neutral scalar is a linear combination of P1 M and (I-P1) M.
    assert abs((D.conj().T@np.diag([1,0,0])@D)[0,2])<1e-12
    t=read(cfg["tree"]);semi=read(cfg["semileptonic"]);inc=read(cfg["inclusive"])
    checkin(abs(K[0,2]/K[1,2]),bounds(t["ub_over_cb"]));checkin(abs(K[1,2]),bounds(semi["BelleII_Vcb"]))
    checkin(yd,bounds(y["down"][0]));checkin(ys,bounds(y["down"][1]))
    checkin(point["fplus"],bounds(t["K_form"]));checkin(abs(K[0,1])*point["fplus"],bounds(t["K_product"]))
    if not locked:
        checkin(abs(K[0,0]),bounds(t["ud_beta"]))
        checkin(abs(K[1,2]),bounds(inc["inclusive_holdout"]["Vcb"]))
    lam=K[:,0].conjugate()*K[:,2]
    dk=K[1,2]*K[0,0].conjugate()*K[1,1].conjugate()*K[0,1]
    M=lam@kernel@lam
    phase=np.angle(M/dk**2)
    beta=np.angle(lam[2]**2/lam[1]**2)/2
    expected=u*math.sin(phi)/(d*c-u*math.cos(phi))
    np.testing.assert_allclose(math.tan(beta),expected,rtol=1e-11)
    return {"phi_deg":math.degrees(phi),"b_GeV2":b,"mu_A_GeV":mu,"A_vev_GeV":avev,
         "cos_phi":cosphi,"rotations_u_d_w":[u,d,w],"Y_u":pack(muq),"Y_d":pack(mdq),
         "weak_moduli":np.abs(K).tolist(),"Kus_times_fplus":abs(K[0,1])*point["fplus"],
         "J":(K[0,0]*K[1,1]*K[0,1].conj()*K[1,0].conj()).imag,
         "beta_deg":math.degrees(beta),"reduced_S":math.sin(phase),"direct_C":0.}


def norm(x):
    if isinstance(x,list):return [norm(v) for v in x]
    if isinstance(x,dict):return {k:norm(v) for k,v in x.items()}
    if isinstance(x,float):return 0. if abs(x)<1e-13 else float(f"{x:.11g}")
    return x


def main():
    p=argparse.ArgumentParser();p.add_argument("--check",action="store_true");args=p.parse_args()
    cfg=read("inputs.json");kernel=short_kernel(cfg);w=read(cfg["fx"])["witness"]
    fx={"up_tangent":w["up_tangent"],"cb":w["Vcb"],"yd":w["down_yukawa"],
        "ys":w["strange_yukawa"],"fplus":.965}
    results={}
    for name,point,locked in [("orthogonal",fx,True),("offset",cfg["calibrated_competitor"],False)]:
        d=calibrate(cfg,locked)
        dom=holdout_bounds(d,kernel)
        branches=[evaluate(cfg,point,sgn,locked,kernel) for sgn in (True,False)]
        np.testing.assert_allclose(branches[0]["reduced_S"],-branches[1]["reduced_S"],atol=2e-13)
        for i,k in enumerate(["S_reduced_positive_outer","S_reduced_negative_outer"]):
            checkin(branches[i]["reduced_S"],dom[k])
        obs=read(cfg["decay"])["BaBar_Belle_2015"]
        slo,shi=bounds(obs["minus_eta_S"])
        assert bounds(obs["C"])[0]<0<bounds(obs["C"])[1]
        assert slo<dom["S_reduced_positive_outer"][0]<=dom["S_reduced_positive_outer"][1]<shi
        assert dom["S_reduced_negative_outer"][1]<slo
        phase_budget=math.asin(slo)+math.asin(dom["S_reduced_positive_outer"][0])
        dom["negative_branch_min_extra_phase_deg"]=math.degrees(phase_budget)
        dom["mixing_only_min_relative_extra"]=math.sin(min(math.pi/2,phase_budget))
        results[name]={"calibration":d,"holdout":dom,"points":branches}
    assert results["orthogonal"]["calibration"]["ub_max"]<results["orthogonal"]["calibration"]["ub_required_by_nuclear"]
    out=norm(results)
    if args.check:assert out==read("results.json")
    else:(HERE/"results.json").write_text(json.dumps(out,indent=2)+"\n",encoding="utf-8")
    for k,e in out.items():
        print(k,e["calibration"]["cos_phi"],e["holdout"]["S_reduced_positive_outer"],
              e["points"][0]["reduced_S"])

if __name__=="__main__":
    main()
