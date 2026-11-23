"""Reconstruct all real left frames for three fixed semileptonic moduli."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

BASE = Path(__file__).resolve().parent


def read(path):
    return json.loads((BASE / path).read_text(encoding="utf-8-sig"))


def parent_core(data):
    cfg = read(data["parent"])
    table = json.loads((BASE.parent / "1515" / cfg["matched_spectra"]).read_text(encoding="utf-8-sig"))
    d, e = np.split(np.array(table["scales"][cfg["scale"]])[:, 0], 2)
    s = np.r_[np.sqrt(np.maximum(d, e)[:2] * np.minimum(d, e)[1:]), 1.2 * max(d[-1], e[-1])]
    poles = np.r_[0., s*s]
    def amp(values):
        residues = [np.prod(p-values**2)/np.prod(p-np.delete(poles, i)) for i, p in enumerate(poles)]
        assert min(residues) > 0 and np.isclose(sum(residues), 1)
        return np.sqrt(residues)
    h = amp(d)
    magn_e = np.prod(d/e) * amp(e)
    gram = np.diag(s) @ (np.eye(3)-np.outer(h[1:], h[1:])) @ np.diag(s)
    values, real_frame = np.linalg.eigh(gram)
    np.testing.assert_allclose(values, d**2, rtol=1e-8)
    return s, d, e, h, magn_e, real_frame


def mixing_core(core, phi):
    s, d, e, hdabs, heabs, real_frame = core
    z = np.exp(1j * math.radians(phi))
    r2 = (heabs[1:]/hdabs[1:])**2
    a, q, c = 9-4*r2, (6+4*r2)*z.real, 1-r2
    disc = q*q-4*a*c
    assert disc.min() > 0 and q.min() > 0
    x = 2*c/(q+np.sqrt(disc))
    bare = hdabs[1:]/abs(1+2*x*z)
    hd = np.r_[hdabs[0], bare*(1+2*x*z)]
    he = np.r_[hdabs[0], bare*(1-3*x*z)]
    v = (hd[1:]/abs(hd[1:]))[:, None]*real_frame
    y0 = np.column_stack([np.zeros(3), np.diag(s)])
    for h, target in [(hd, d), (he, e)]:
        u = h/np.linalg.norm(h)
        gram = y0@(np.eye(4)-np.outer(u,u.conj()))@y0.T
        np.testing.assert_allclose(np.linalg.eigvalsh(gram), target**2, rtol=2e-8)
    det = np.linalg.det(np.vstack([y0,-hd.conj()]))
    assert det.real > 0 and abs(det.imag) < 1e-17
    return v, x, det.real


def polynomial_roots(coeffs):
    """Bracket on all intervals of monotonicity. No repeated or infinite roots at these inputs."""
    p = np.asarray(coeffs)/max(abs(np.asarray(coeffs)))
    assert abs(p[-1]) > 1e-10
    if len(p) == 2:
        return [-p[0]/p[1]]
    critical = polynomial_roots(np.arange(1,len(p))*p[1:])
    radius = 1+max(abs(p[:-1]/p[-1]))
    cuts = [-radius] + [v for v in critical if -radius < v < radius] + [radius]
    evaluate = lambda x: np.polynomial.polynomial.polyval(x,p)
    assert all(abs(evaluate(v)) > 2e-13 for v in cuts)
    found=[]
    for lo,hi in zip(cuts,cuts[1:]):
        flo=evaluate(lo)
        if flo*evaluate(hi) > 0:
            continue
        for _ in range(90):
            mid=(lo+hi)/2
            if flo*evaluate(mid) <= 0:
                hi=mid
            else:
                lo=mid
                flo=evaluate(lo)
        found.append((lo+hi)/2)
    return found


def all_frames(v, a, b, c):
    bm = np.outer(v[:,2],v[:,2].conj()).real
    sm = np.outer(v[:,1],v[:,1].conj()).real
    ev, e = np.linalg.eigh(bm-b*b*np.eye(3))
    assert np.min(np.abs(ev)) > 1e-10
    pos = np.count_nonzero(ev>0)
    if pos in (0,3):
        return [],0
    i = np.argmin(ev) if pos==2 else np.argmax(ev)
    j,k = [axis for axis in range(3) if axis!=i]
    u=e[:,i]
    u1=math.sqrt(abs(ev[i]/ev[j]))*e[:,j]
    u2=math.sqrt(abs(ev[i]/ev[k]))*e[:,k]
    poly_vec=np.array([u+u1,2*u2,u-u1])
    mat=sm-a*a*np.eye(3)
    poly=np.zeros(5)
    for i, vi in enumerate(poly_vec):
        for j, vj in enumerate(poly_vec):
            poly[i+j]+=vi@mat@vj
    roots=polynomial_roots(poly)
    solutions=[]
    for t in roots:
        x=poly_vec[0]+t*poly_vec[1]+t*t*poly_vec[2]
        x/=np.linalg.norm(x)
        if x[0]<0:
            x=-x
        plane=np.linalg.svd(x[None,:])[2][1:].T
        vals, axes=np.linalg.eigh(plane.T@bm@plane)
        assert vals[1]>vals[0]+1e-7
        f=(c*c-vals[0])/(vals[1]-vals[0])
        if not 0 <= f <= 1:
            continue
        for sign in (-1,1):
            y=plane@axes@np.array([math.sqrt(1-f), sign*math.sqrt(f)])
            if y[0]<0:
                y=-y
            o=np.array([x,y,np.cross(x,y)])
            np.testing.assert_allclose(o@o.T,np.eye(3),atol=1e-11)
            assert abs(np.linalg.det(o)-1)<1e-11
            ck=o@v
            np.testing.assert_allclose(abs(ck[[0,0,1],[1,2,2]]),[a,b,c],rtol=1e-9)
            gamma=float(np.angle(-ck[0,0]*ck[0,2].conj()/(ck[1,0]*ck[1,2].conj()),deg=True)%180)
            jinv=float(np.imag(ck[0,0]*ck[1,1]*ck[0,1].conj()*ck[1,0].conj()))
            solutions.append((gamma,jinv,o))
    return sorted(solutions,key=lambda s:s[0]),len(roots)


def continuous_certificate(core, phase_max, lower, upper, ratio_lo, cb_lo, cb_hi, scale, us_scale):
    _, x0, _=mixing_core(core,0)
    _, x1, _=mixing_core(core,phase_max)
    xmin=np.minimum(x0,x1)
    xmax=np.maximum(x0,x1)
    rad=math.radians(phase_max)
    qmin=1+2*np.minimum(xmin,xmin*math.cos(rad))
    assert min(qmin)>0
    eps=[]
    for weights in (core[-1]**2).T:
        total=0
        for i,j in itertools.combinations(range(3),2):
            diff=max(abs(xmin[i]-xmax[j]),abs(xmax[i]-xmin[j]))
            sin_max=min(1,2*math.sin(rad)*diff/(qmin[i]*qmin[j]))
            total+=weights[i]*weights[j]*sin_max**2
        assert 4*total<1
        eps.append(math.sqrt(2*total/(1+math.sqrt(1-4*total))))
    a_lo,a_hi=lower*us_scale[0],upper*us_scale[1]
    c_lo,c_hi=cb_lo*scale[0],cb_hi*scale[1]
    b_lo=c_lo*ratio_lo
    # |Vub| <= |Vcb| here because the upper baryonic ratio is below 1.
    ud=math.sqrt(1-a_hi*a_hi-c_hi*c_hi)
    cd=math.sqrt(a_lo*a_lo-c_hi*c_hi)
    bounds=[(eps[0],ud),(eps[2],b_lo),(eps[0],cd),(eps[2],c_lo)]
    degrees=sum(math.degrees(math.asin(min(1,err/amp))) for err,amp in bounds)
    return {"column_imag_norm_upper":eps,"angle_to_zero_mod_180_deg_max":degrees}


def calculate():
    data=read("inputs.json")
    core=parent_core(data)
    tree, sem, charge=read(data["historical_tree"]),read(data["belle"]),read(data["charged_gamma"])["chargedB"]
    a0,da=tree["K_product"]
    f,df=tree["K_form"]
    cc,*cerr=sem["BelleII_Vcb_BGL"]
    rb,*rerr=tree["ub_over_cb"]
    hi=read(data["high_mixing"])["high_scale"]["angles_theta12_theta13_theta23_delta"][2][0]
    lo=read(data["low_mixing"])["quark_sines_and_phase"][2][0]
    sm_like=math.sin(hi)/lo
    gamma_interval=[charge["central_degrees"]-charge["lower_error"]-charge["additional_degrees"],
                    charge["central_degrees"]+charge["upper_error"]+charge["additional_degrees"]]
    model_results={}
    example_phase=data["secondary_examples"]["additional_phase_after_transport_check"]
    for tag,scale in [("direct",1),("common_cb_ub_transfer",sm_like)]:
        points=list(data["fixed_phase_degrees"]) + [example_phase]
        fixed={}
        for deg in points:
            v, _, det=mixing_core(core,deg)
            sol,nroot=all_frames(v,a0/f,rb*cc*scale,cc*scale)
            vb=v[:,2]
            max_tb=(1+abs(vb.T@vb))/2
            assert abs(np.linalg.eigvalsh(np.outer(vb,vb.conj()).real)[-1]-max_tb)<1e-12
            passing=[s for s in sol if gamma_interval[0]<=s[0]<=gamma_interval[1]]
            fixed[str(deg)]={"first_ray_real_roots":nroot,"gamma_mod_180_each_frame":[s[0] for s in sol],
                             "max_abs_Vtb_squared":max_tb,"passing_gamma_deg":[s[0] for s in passing]}
            if passing:
                g,j,o=passing[0]
                chosen=o@v
                signed=float(np.angle(-chosen[0,0]*chosen[0,2].conj()/(chosen[1,0]*chosen[1,2].conj()),deg=True))
                fixed[str(deg)]["one_passing_frame"]={"O":o.tolist(),"gamma_deg":g,"quartet_J":j,
                                                      "colored_tree_determinant":det, "gamma_signed_deg":signed,
                                                      "abs_Vtd_over_Vts":float(abs(chosen[2,0]/chosen[2,1]))}
        model_results[tag]=fixed
    ex=data["secondary_examples"]
    assert rb+2*sum(rerr)<1
    cert=continuous_certificate(core,ex["lower_abs_phase_interval_deg"],(a0-2*da)/(f+2*df),
                                (a0+2*da)/(f-2*df),rb-2*sum(rerr),cc-2*sum(cerr),cc+2*sum(cerr),
                                ex["common_amplitude_transfer_box"],ex["us_transfer_box"])
    gamma_lo=charge["central_degrees"]-2*charge["lower_error"]-charge["additional_degrees"]
    gamma_hi=charge["central_degrees"]+2*charge["upper_error"]+charge["additional_degrees"]
    assert cert["angle_to_zero_mod_180_deg_max"]<min(gamma_lo,180-gamma_hi)
    td,*terr=read(data["charged_gamma"])["mixing"]["td_over_ts"]
    mix_bound=[td-2*sum(terr),td+2*sum(terr)]
    return {"parent_singular_sources":core[0].tolist(), "tree_centers_us_ub_cb":[a0/f,rb*cc,cc],
            "gamma_holdout_one_and_two_widths":[gamma_interval,[gamma_lo,gamma_hi]],"transfer_example":sm_like,
            "neutral_B_mixing_two_width":mix_bound,"continuous_workbox_certificate":cert,"solutions":model_results}


def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:
            compare(a[k],b[k])
    else:
        np.testing.assert_allclose(a,b,atol=3e-8,rtol=3e-8)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    ab=p.add_mutually_exclusive_group(required=True)
    ab.add_argument("--check",action="store_true")
    ab.add_argument("--write",action="store_true")
    opt=p.parse_args()
    result=calculate()
    path=BASE/"results.json"
    if opt.write:
        with path.open("x",encoding="utf8") as out:
            json.dump(result,out,indent=2,allow_nan=False)
            out.write("\n")
    else:
        compare(result,read("results.json"))
        print("1516 six spectra, calibrated real frames, phase and continuum certificate verified.")
