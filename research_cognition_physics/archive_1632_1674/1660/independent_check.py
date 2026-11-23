"""Majorana/Dirac equality of physical left spectral projectors; all-phase reactor bound."""
import argparse
import cmath
import json
import math
from pathlib import Path
import runpy

import numpy as np

HERE=Path(__file__).resolve().parent


def read(rel):
    return json.loads((HERE/rel).read_bytes())


def left_observables(H,Ue):
    lam,L=np.linalg.eigh(H)
    V=Ue.conj().T@L
    p=np.abs(V)**2
    J=(V[0,0]*V[1,1]*V[0,1].conjugate()*V[1,0].conjugate()).imag
    return [4*p[0,2]*(1-p[0,2]),p[0,1]/(1-p[0,2]),p[1,2]/(1-p[0,2]),J,
            math.sqrt(max(0,np.dot(p[0],lam)))]


def spectral_matches(Ue,old):
    W=np.array([[0,1],[1,3],[1,1]],dtype=float)
    gram=W.T@W
    ghalf=(gram+math.sqrt(6)*np.eye(2))/math.sqrt(13+2*math.sqrt(6))
    np.testing.assert_allclose(ghalf@ghalf,gram)
    rng=np.random.default_rng(67181)
    rows={}
    for name,scenario in old["scenarios"].items():
        rows[name]={}
        for eta,point in scenario["small_root_by_eta_deg"].items():
            z=point["r"]*cmath.exp(1j*math.radians(int(eta)))
            P=np.diag(np.exp(1j*np.array([0,*point["mass_P_alpha_beta_rad"]])))
            left=P.conj()@W
            m2=math.sqrt(scenario["dm21"])
            m3=math.sqrt(scenario["dm21"]+scenario["dm32"])
            m0=math.sqrt(m2*m3/(6*abs(z)))
            M=m0*P@W@np.diag([1,z])@W.T@P
            right=m0*np.diag([1,z.conjugate()])@ghalf
            D=left@right
            HM=M.conj().T@M
            HD=D@D.conj().T
            np.testing.assert_allclose(HD,HM,atol=1e-16)
            residue=np.linalg.norm(HD-HM)
            # Reconstruct the right Gram from a left inverse; verify all choices of a squared factor agree.
            Gdata=np.linalg.pinv(left)@HD@np.linalg.pinv(left).conj().T
            norms=np.sqrt(Gdata.diagonal().real)
            overlap=Gdata[0,1]/np.prod(norms)
            np.testing.assert_allclose(norms[1]/norms[0],abs(z)*math.sqrt(11/2))
            np.testing.assert_allclose(overlap,4/math.sqrt(22)*z/abs(z))
            unitary,_=np.linalg.qr(rng.normal(size=(2,2))+1j*rng.normal(size=(2,2)))
            X=right@unitary
            loc=np.linalg.solve(right,X)
            np.testing.assert_allclose(loc@loc.conj().T,np.eye(2),atol=2e-14)
            values=left_observables(HD,Ue)
            np.testing.assert_allclose(values,[*point["A13_solar_atmospheric"],point["J_dirac"],point["m_beta_eV"]],rtol=3e-8)
            mbb=abs((Ue.T@M@Ue)[0,0])
            np.testing.assert_allclose(mbb,point["m_beta_beta_eV"],rtol=2e-8)
            # The symmetric mass matrix of the five Weyl states has a zero active-active block.
            five=np.zeros((5,5),complex)
            five[:3,3:]=D
            five[3:,:3]=D.T
            sv=np.sort(np.linalg.svd(five,compute_uv=False))
            np.testing.assert_allclose(sv,[0,m2,m2,m3,m3],rtol=2e-9,atol=1e-14)
            rows[name][eta]={"left_gram_max_error_eV2":float(residue),
                            "width_angles_J_mbeta_eV":values,"majorana_mbb_eV":mbb,
                            "dirac_active_LNV_mass_block_norm":np.linalg.norm(five[:3,:3]),
                            "relative_right_row_norms":norms[1]/norms[0],
                            "right_normalized_overlap_modulus":abs(overlap),"right_overlap_phase_deg":np.angle(overlap,deg=True)}
    return rows


def overlap_exclusions(old):
    iv=runpy.run_path(str(HERE.parent/"1486/check.py"))
    I=iv["I"]
    source=read("../1484/inputs.json")
    R=(1+iv["width"](source["T2K_dm32_NO_eV2"])/iv["width"](source["KamLAND_dm21_eV2"])).sqrt()
    # Deliberately weaken the charged projector determined at the fixed 1659 GJ point.
    ue=I(997619000000,997620000000)
    assert ue.lo/iv["SCALE"]<old["abs_Ue"][0][0]<ue.hi/iv["SCALE"]
    be=(1-ue.square()).sqrt()
    target=iv["width"](source["DayaBay_sin2_2theta13"])
    outputs={}
    for label in ("0",".22",".23"):
        cap=I.point(label)
        k=I((6*(1-cap.square())).sqrt().lo,I.point(6).sqrt().hi)
        points=[[],[]]
        for Ri in R.parts(1024):
            base=k*(Ri+1/Ri)
            arg=I((base-8*cap).lo,(base+8*cap).hi)
            # z is the norm ratio of the two right row vectors, separate from the normalized overlap.
            h=arg+(arg.square()-88).sqrt()
            for j,z in enumerate((4/h,h/22)):
                pr=(z/k-1/(3*Ri))/(Ri-1/Ri)
                pr=I(max(0,pr.lo),pr.hi)
                points[j].append(pr.sqrt())
        report=[]
        for j,ls in enumerate(points):
            bare=iv["hull"](ls)
            low=I(bare.lo,bare.lo);high=I(bare.hi,bare.hi)
            lower=ue*low-be*(1-low.square()).sqrt()
            upper=ue*high+be*(1-high.square()).sqrt()
            assert bare.hi<ue.lo
            u_phys=I(max(0,lower.lo),upper.hi)
            w=u_phys.square()
            assert w.hi < iv["SCALE"]//2
            # A13=4 p(1-p) is increasing throughout this p range.
            w0=I(w.lo,w.lo);w1=I(w.hi,w.hi)
            amp=I((4*w0*(1-w0)).lo,(4*w1*(1-w1)).hi)
            if j==0: assert amp.hi < target.lo
            else: assert amp.lo > target.hi
            report.append({"bare_e3":bare.values(),"charged_e3":u_phys.values(),"A13":amp.values()})
        outputs[label]=report
    return outputs


def tidy(x):
    if isinstance(x,dict):return {k:tidy(v) for k,v in x.items()}
    if isinstance(x,list):return [tidy(v) for v in x]
    if isinstance(x,(float,np.floating)):return 0 if abs(x)<1e-16 else float(f"{x:.10g}")
    return x


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--check",action="store_true")
    args=ap.parse_args()
    old=read("../1659/results.json")
    oldcode=runpy.run_path(str(HERE.parent/"1659/check.py"))
    Ue,_=oldcode["charged_basis"](read("../1659/inputs.json"))
    data=tidy({"same_left_matrix_in_each_solution":spectral_matches(Ue,old),
               "excluded_overlap_below_or_equal_to":overlap_exclusions(old)})
    path=HERE/"independent_results.json"
    if args.check:assert data==read("independent_results.json")
    else:path.write_text(json.dumps(data,indent=2)+"\n",encoding="utf-8")
    print("1660 independent PASS: four full matrices, certified two-branch reactor exclusion")
    print(data["excluded_overlap_below_or_equal_to"][".23"])
