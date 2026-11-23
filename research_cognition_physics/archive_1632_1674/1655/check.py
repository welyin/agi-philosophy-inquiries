"""Projector identities, explicit matched Yukawas, and interval bounds for the 1486 charged textures."""
import argparse
import importlib.util
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE/path).read_bytes())


spec = importlib.util.spec_from_file_location("r1486", HERE.parent/"1486"/"check.py")
old = importlib.util.module_from_spec(spec)
spec.loader.exec_module(old)


def svd(Y):
    U, d, Vt = np.linalg.svd(Y)
    return U[:, ::-1], d[::-1], Vt[::-1, :].T


def off(G):
    return np.abs([G[0,1], G[1,0]])


def calculate(cfg):
    t = read(cfg["texture"])
    form = read(cfg["charges"])
    hist = read(cfg["old_results"])["textures"]
    masses = read("../1480/inputs.json")
    mu = read(cfg["muon"])["mu_mass_GeV"]
    v = form["v_GeV"]
    xi = v*form["sin_theta"]/form["f_GeV"]
    KK = (np.asarray(form["left"])[:,None]+np.asarray(form["right"]))%form["N"]
    assert np.array_equal(KK[:2,:2], [[2,0],[2,0]])
    q = masses["lattice"]["md"][0]/masses["lattice"]["ms"][0]
    _,s,ds = next(row for row in read("../1476/inputs.json")["edges"] if row[0]=="cd")
    ratio = masses["muonium"]["pole_ratio"][0]
    I = old.I
    qb=old.width(masses["lattice"]["md"])/old.width(masses["lattice"]["ms"])
    err=sum(map(I.point,ds),I.point(0))
    sb=old.hull(((I.point(s)-2*err)/I.point("1.02"),(I.point(s)+2*err)/I.point(".98")))
    b=math.sqrt(q*q*(1-s*s)+s*s)
    a,c=q/b,(1-q*q)*s*math.sqrt(1-s*s)/b
    Yd=np.array([[0,b],[a,c]])
    P=np.diag([1.,0.])
    U,d,V=svd(Yd)
    np.testing.assert_allclose(d,[q,1])
    Qoff=off(U.T@(2*xi*Yd@P)@V)
    sinR=2*a*c/(1-q*q)
    np.testing.assert_allclose(Qoff,[q*xi*sinR,xi*sinR])
    rrun = read(cfg["scalar_competitors"])["two_width_solutions"]
    mh,W0=cfg["mh_GeV"],cfg["non_LFV_width_GeV"]
    def branch(g):
        w=mh*np.dot(g,g)/(8*math.pi)
        return w/(w+W0)
    limit=cfg["CMS_observed_95"]
    gmax=math.sqrt(8*math.pi*W0/mh*limit/(1-limit))
    out={}
    for name,(ca,cb,cc) in t["clebsch"].items():
        E=np.array([[0,ca*a],[cb*b,cc*c]])
        L,yp,R=svd(E)
        high=yp[1]/yp[0]
        np.testing.assert_allclose(high,hist[name]["real_minus_witness"]["lepton_mass_ratio"])
        sin2e=2*abs(L[0,0]*L[0,1])
        mass=np.array([mu/ratio,mu])
        factors=mass/yp
        target=L@np.diag(mass)@R.T
        TR=R@np.diag(factors)@R.T
        TL=L@np.diag(factors)@L.T
        mass_right=E@TR
        mass_left=TL@E
        np.testing.assert_allclose(mass_right,target)
        np.testing.assert_allclose(mass_left,target)
        def poles(matrix):
            return off(L.T@matrix@R/v)
        Vh=2*xi*P@E
        right=poles(Vh@TR)
        left=poles(TL@Vh)
        common=poles(Vh*mu/yp[1])
        np.testing.assert_allclose(right,xi*sin2e*mass[::-1]/v)
        np.testing.assert_allclose(common[0]/common[1],high)
        np.testing.assert_allclose(left,right*np.array([high/ratio,ratio/high]))
        double=(right@right)/(mu/v)**2/(Qoff@Qoff)
        np.testing.assert_allclose(double, (ca*cc*(1-q*q)/(yp[1]**2-yp[0]**2))**2*(1+ratio**-2)/(1+q*q))

        # Exact outer rational enclosures on the inherited two-width quark rectangle.
        boxes=[]
        for qi in qb.parts(20):
            for si in sb.parts(20):
                bi=(qi.square()*(1-si.square())+si.square()).sqrt()
                ai,ci=qi/bi,(1-qi.square())*si*(1-si.square()).sqrt()/bi
                se,lg=old.charged(qi,si,(ca,cb,cc))
                sin2=2*se*(1-se.square()).sqrt()
                sr=2*ai*ci/(1-qi.square())
                ratios=(sin2/sr).square()*(1+1/I.point(ratio).square())/(1+qi.square())
                boxes.append((sin2,ratios))
        ang=old.hull(pair[0] for pair in boxes).values()
        cross=old.hull(pair[1] for pair in boxes).values()
        ob0=xi*mu/v
        def br_angle(sine,z):
            return branch(z*mu/v*sine*np.array([1.,1/ratio]))
        capped=[gmax/(mu/v*x*math.sqrt(1+1/ratio**2)) for x in reversed(ang)]
        out[name]={"G_high":high,"sin_2theta_e":sin2e,"sin_2theta_e_interval":ang,
          "h_e_mu_and_mu_e_abs_right":right.tolist(),"h_e_mu_and_mu_e_abs_left":left.tolist(),
          "h_e_mu_and_mu_e_abs_common_vertex":common.tolist(),"BR_right_at_xi_ref":branch(right),
          "BR_right_range_cover": [br_angle(x,xi) for x in ang],
          "quark_lepton_normalized_width_double_ratio":double,"double_ratio_cover":cross,
          "CMS_reference_ceiling_abs_xi_central":gmax/(mu/v*sin2e*math.sqrt(1+1/ratio**2)),
          "CMS_reference_ceiling_abs_xi_cover":capped,
          "1654_xi_references":[{"xi":row["xi"],"central_BR":br_angle(sin2e,row["xi"]),
                                  "BR_cover":[br_angle(x,row["xi"]) for x in ang]} for row in rrun]}
    assert out["small_left"]["double_ratio_cover"][1]<out["large_left"]["double_ratio_cover"][0]
    # With the same charges: shared factor F(t)=2-t and t=(r/f)^2 turns the nonzero power from 2 to 0 at the vacuum.
    from numpy.polynomial import Polynomial
    r=Polynomial([1.,1.])
    dressed=r**2*(2-r**2)
    assert dressed(0)==1 and dressed.deriv()(0)==0
    replacement=np.array([0,0,1])
    changed=(np.array(form["left"])[:,None]+replacement)%4
    assert changed[1,0]==0 and np.array_equal(changed[:,1:],KK[:,1:])
    return {"xi_reference":xi,"down_quark_chiral_norm_by_ms_over_v":Qoff.tolist(),
            "sin_2theta_down_right":sinR,"quark_q_s": [q,s],
            "CMS_reference_norm_from_exact_BR_formula":gmax,"branches":out,
            "same_charge_dressing_normalized_value_and_slope":[dressed(0),dressed.deriv()(0)]}


def canonical(x):
    if isinstance(x,dict):
        return {k:canonical(v) for k,v in x.items()}
    if isinstance(x,list):
        return [canonical(z) for z in x]
    return float(f"{x:.12g}") if isinstance(x,(float,np.floating)) else x


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    cli.add_argument("--check",action="store_true")
    args=cli.parse_args()
    data=canonical(calculate(read("inputs.json")))
    if args.check:
        assert data==read("results.json")
    else:
        (HERE/"results.json").write_text(json.dumps(data,indent=2)+"\n",encoding="utf-8")
    print("1655 OK: two-sided correlated pole transfer, three Clebsches, shared quark bounds")
    for name,x in data["branches"].items():
        print(name,x["h_e_mu_and_mu_e_abs_right"],x["BR_right_at_xi_ref"],x["double_ratio_cover"],
              x["CMS_reference_ceiling_abs_xi_central"],x["BR_right_range_cover"])
