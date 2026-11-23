"""Independent charged reconstruction, exact-electron phase elimination and Takagi tests."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from numpy.polynomial import Polynomial

HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE/path).read_bytes())


def charged(point):
    rows = point["marginal_coordinates"]
    down = np.array(rows[:3])/rows[2]
    s12,s13,s23 = np.sin(rows[3:6])
    c12,c13,c23 = np.cos(rows[3:6])
    ph = np.exp(1j*rows[6])
    v = np.array([[c12*c13,s12*c13,s13*ph.conjugate()],
                  [-s12*c23-c12*s23*s13*ph,c12*c23-s12*s23*s13*ph,s23*c13],
                  [s12*s23-c12*c23*s13*ph,-c12*s23-s12*c23*s13*ph,c23*c13]])
    h = (v*down**2)@v.conj().T
    phase = np.r_[1,h[0,1:]/abs(h[0,1:])]
    h = phase[:,None]*h*phase.conj()
    a = np.sqrt(h[0,0].real)
    c,e = h[0,1:].real/a
    f0 = np.sqrt(h[2,2].real-e**2)
    b = -np.prod(down)/(a*f0)
    z = (h[1,2]-c*e)/f0
    t = point["t"]
    ye = np.array([[0,b,0],[a,3*c,t*e],[0,t*z,t*f0]])*point["k_times_yb"]
    left,y,_ = np.linalg.svd(ye)
    left,y = left[:,::-1],y[::-1]
    np.testing.assert_allclose(y,[point["predicted_e"],*rows[7:]],rtol=5e-9)
    return left,y


def mass_source(r, sign, dm21):
    a,b = np.array([0,1,1]),np.array([1,3,1])
    eta = np.exp(1j*sign*2*np.pi/3)
    m = np.outer(a,a)+r*eta*np.outer(b,b)
    w = np.linalg.eigh(m.conj().T@m)[1]
    for i in (1,2):
        w[:,i] *= np.exp(-.5j*np.angle(w[:,i].T@m@w[:,i]))
    diag = np.diag(w.T@m@w).real
    scale = math.sqrt(dm21)/diag[1]
    masses = np.array([0,diag[1]*scale,diag[2]*scale])
    np.testing.assert_allclose(w.T@(m*scale)@w,np.diag(masses),atol=1e-15)
    return w,m*scale,masses


def rates(ue, core, alpha, beta):
    w,m,masses = core
    p = np.exp(1j*np.array([0,alpha,beta]))
    mix = ue.conj().T@(p.conj()[:,None]*w)
    weights = abs(mix)**2
    e3 = weights[0,2]
    a = [4*e3*(1-e3),weights[0,1]/(1-e3),weights[1,2]/(1-e3)]
    j = (mix[0,0]*mix[1,1]*mix[0,1].conjugate()*mix[1,0].conjugate()).imag
    mb = math.sqrt(weights[0]@masses**2)
    mbb = abs(np.sum(masses*mix[0]**2))
    ep = p*ue[:,0]
    np.testing.assert_allclose(mbb,abs(ep.T@m@ep),atol=3e-16)
    terms = masses[1:]*weights[0,1:]
    return {"mass_P_rad":[alpha,beta],"angles":a,"J":j,"m_beta_eV":mb,"m_bb_eV":mbb,
            "m_bb_if_only_oscillation_H_is_fixed_eV":[abs(terms[0]-terms[1]),sum(terms)]}


def align(electron_coeff, p1, p3):
    # For B=e^{ib}=(1+it)^2/(1+t^2), solve linearly for cos(a),sin(a).
    den, bnum = np.array([1.,0.,1.]),np.array([1.,2j,-1.])
    equations=[]
    for col,target in ((0,p1),(2,p3)):
        p,q,r = electron_coeff[col]
        g = (p*den+r*bnum).conj()*q
        v = (target-abs(p)**2-abs(q)**2-abs(r)**2)*den-2*np.real(p.conjugate()*r*bnum)
        equations.append((Polynomial(2*g.real),Polynomial(-2*g.imag),Polynomial(v)))
    l,m,y = equations[0]
    L,M,Y = equations[1]
    det,nx,ny = l*M-m*L,y*M-m*Y,l*Y-y*L
    real = []
    for v in (nx*nx+ny*ny-det*det).roots():
        t = v.real
        if abs(v.imag)<1e-7 and abs(det(t))>1e-14:
            # Return phases of the symmetric mass, the opposite of the Takagi mode phases.
            real.append((-np.angle((nx(t)+1j*ny(t))/det(t)),-2*math.atan(t)))
    return real


def phase_enclosure(ue, w, p1, p3, n=1024):
    # Every pair on the torus is within pi/n in both angles from some center.
    ph = np.exp(2j*np.pi*np.arange(n)/n)
    def mesh(l,k):
        p,q,r = ue[:,l].conj()*w[:,k]
        center = abs(p+q*ph[:,None]+r*ph[None,:])
        width = (abs(q)+abs(r))*np.pi/n + 2e-9
        return np.maximum(center-width,0)**2,np.minimum(center+width,1)**2
    lo,hi = mesh(0,0)
    ok = (lo<=p1)&(p1<=hi)
    lo,hi = mesh(0,2)
    ok &= (lo<=p3)&(p3<=hi)
    if not np.any(ok):
        return None
    lo,hi = mesh(1,2)
    return [np.min(lo[ok])/(1-p3),np.max(hi[ok])/(1-p3)]


def roots(R):
    return np.sort(np.roots([121,-(6*(R+1/R)+16),4]))


def run():
    cfg = read("inputs.json")
    osc = read(cfg["oscillation"])
    point = read(cfg["charged"])["extended_box_GJ_test"]
    ue,y = charged(point)
    osc_center = [osc["DayaBay_sin2_2theta13"][0],osc["JUNO_sin2_theta12_NO"][0]]
    p3 = (1-math.sqrt(1-osc_center[0]))/2
    p1 = (1-p3)*(1-osc_center[1])
    dm21 = osc["KamLAND_dm21_eV2"][0]
    rmain = math.sqrt(1+osc["T2K_dm32_NO_eV2"][0]/dm21)
    working_R, alternative_ds = 5.675425,7.65e-5
    assembled = {}
    certificate = {}
    for label,R,ds in (("both_gaps_central",rmain,dm21),
                       ("both_electron_rates_central",working_R,alternative_ds)):
        solutions={}
        for root,r in enumerate(roots(R)):
            for sign in (-1,1):
                key=f"root{root}_eta{sign*120}"
                core = mass_source(r,sign,ds)
                amps = (ue[:,0].conj()[:,None]*core[0]).T
                vals = [rates(ue,core,a,b) for a,b in align(amps,p1,p3)]
                for val in vals:
                    np.testing.assert_allclose(val["angles"][:2],osc_center,atol=2e-12)
                solutions[key] = vals
                if label=="both_gaps_central":
                    bound=phase_enclosure(ue,core[0],p1,p3)
                    assert bound is None or bound[1]<.49
                    certificate[key]=bound
        assembled[label]={"R":R,"dm21":ds,"dm32":ds*(R**2-1),"phase_solutions":solutions}
    saved=read("results.json")
    points={}
    for label,scenario in saved["scenarios"].items():
        d21,d32=scenario["dm21"],scenario["dm32"]
        R=math.sqrt(1+d32/d21)
        for key,p in scenario["small_root_by_eta_deg"].items():
            core=mass_source(roots(R)[0],np.sign(float(key)),d21)
            audit=rates(ue,core,*p["mass_P_alpha_beta_rad"])
            np.testing.assert_allclose(audit["angles"],p["A13_solar_atmospheric"],atol=2e-9)
            for lhs,rhs in (("J","J_dirac"),("m_beta_eV","m_beta_eV"),("m_bb_eV","m_beta_beta_eV")):
                np.testing.assert_allclose(audit[lhs],p[rhs],atol=2e-10)
            points[f"{label}_eta{key}"]=audit
    preliminary=rates(ue,mass_source(roots(working_R)[0],1,alternative_ds),2.725027882,1.475333567)
    np.testing.assert_allclose(preliminary["angles"],[.0853851,.3040222,.5619923],atol=7e-7)
    fixed_R = 5.68
    fixed_gaps = [alternative_ds, alternative_ds*(fixed_R**2-1)]
    fixed_core = mass_source(roots(fixed_R)[0],1,alternative_ds)
    competitors = [rates(ue,fixed_core,*np.radians(pair)) for pair in ((156,85),(302,230))]
    for entry in competitors:
        assert all(lo<a<hi for (lo,hi),a in zip(saved["joint_windows"],entry["angles"]))
    assert competitors[0]["J"]>0>competitors[1]["J"]
    # Select the large positive root over the whole marginal frequency box. Formula for
    # |W_e3|^2 comes from the heavy projector (H-m2^2*(1-P0))/(m3^2-m2^2).
    def window(v): return v[0]-2*v[1],v[0]+2*v[2]
    slo,shi=window(osc["KamLAND_dm21_eV2"])
    alo,ahi=window(osc["T2K_dm32_NO_eV2"])
    borders=np.linspace(math.sqrt(1+alo/shi),math.sqrt(1+ahi/slo),97)
    amin,amax=1.,0.
    u=abs(ue[0,0]); orth=np.linalg.norm(ue[1:,0])
    for lo,hi in zip(borders[:-1],borders[1:]):
        zmin=math.sqrt((11*roots(lo)[1]-2/lo)/(6*(hi-1/hi)))
        zmax=math.sqrt((11*roots(hi)[1]-2/hi)/(6*(lo-1/lo)))
        emin=max(0,u*zmin-orth*math.sqrt(1-zmin**2))
        emax=min(1,u*zmax+orth*math.sqrt(1-zmin**2))
        amin=min(amin,emin**2)
        amax=max(amax,emax**2)
    assert amax<.5 and 4*amin*(1-amin)>.0899
    return {"charged_yukawas":y.tolist(),"abs_Ue":abs(ue).tolist(),"published_points_recomputed":points,
            "preliminary_point_requested":preliminary,"exact_electron_cases":assembled,
            "same_positive_eta_two_CP_signs":{"dm21_dm32":fixed_gaps,"R":fixed_R,"points":competitors},
            "central_all_phase_cell_atmosphere_covers":certificate,
            "large_root_all_frequencies":{"e3_weight":[amin,amax],"reactor_min":4*amin*(1-amin)}}


def stable(x):
    if isinstance(x,dict):return {k:stable(v) for k,v in x.items()}
    if isinstance(x,list):return [stable(v) for v in x]
    return float(format(x,".10g")) if isinstance(x,(float,np.floating)) else x


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    out=stable(run())
    if args.check:
        assert out==read("independent_results.json")
    else:
        with (HERE/"independent_results.json").open("x",encoding="utf8") as f:
            json.dump(out,f,indent=2)
            f.write("\n")
    print("1659 independent PASS: reconstructed charged matrix, common phases, mass form and torus enclosures")
