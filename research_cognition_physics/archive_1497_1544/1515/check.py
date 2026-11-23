"""One CP phase in a split vectorlike five: masses, the colored determinant and gauge running."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def load(name):
    return json.loads((HERE/name).read_text(encoding="utf-8"))


def turn(i, j, angle):
    r = np.eye(3)
    r[i,i] = r[j,j] = math.cos(angle)
    r[i,j], r[j,i] = math.sin(angle), -math.sin(angle)
    return r


def weights(poles, spectrum):
    roots = spectrum**2
    w = np.array([np.prod(x-roots)/np.prod(x-np.delete(poles,i)) for i,x in enumerate(poles)])
    assert np.all(w>0)
    np.testing.assert_allclose(sum(w),1.)
    return w


def witness(d,e,phi,angles):
    R = np.prod(d/e)
    s = np.r_[np.sqrt(np.maximum(d,e)[:2]*np.minimum(d,e)[1:]),1.2*max(d[2],e[2])]
    poles = np.r_[0.,s*s]
    d_abs = np.sqrt(weights(poles,d))
    e_abs = R*np.sqrt(weights(poles,e))
    np.testing.assert_allclose(d_abs[0],e_abs[0])
    m = [d_abs[0]]
    ell = [0.]
    for u,v in zip(d_abs[1:],e_abs[1:]):
        r2=(v/u)**2
        aa,bb,cc=9-4*r2,-(6+4*r2)*math.cos(phi),1-r2
        disc=bb**2-4*aa*cc
        assert disc>=0
        # Smaller-magnitude root, avoiding subtraction of nearly equal numbers.
        x=2*cc/(-bb+math.sqrt(disc))
        bare=u/abs(1+2*x*np.exp(1j*phi))
        m.append(bare)
        ell.append(bare*x)
    m,ell=np.array(m),np.array(ell)
    O=turn(0,1,angles[0])@turn(0,2,angles[1])@turn(1,2,angles[2])
    Y=O@np.column_stack((np.zeros(3),np.diag(s)))
    out={}
    for sign in (-1,0,1):
        if sign==0:
            # Preserve the same two spectra with a physically different, real set of inputs.
            h_d,h_e=d_abs.copy(),e_abs.copy()
        else:
            phase=np.exp(1j*sign*phi)
            h_d,h_e=m+2*phase*ell,m-3*phase*ell
        grams=[]
        for h,target in ((h_d,d),(h_e,e)):
            u=h/np.linalg.norm(h)
            H=Y@(np.eye(4)-np.outer(u,u.conj()))@Y.T
            np.testing.assert_allclose(np.linalg.eigvalsh(H),target**2,rtol=2e-9,atol=2e-18)
            grams.append(H)
        Hd,He=grams
        # Up flavors define a real nondegenerate diagonal reference.
        j=(Hd[0,1]*Hd[1,2]*Hd[2,0]).imag/np.prod([d[j]**2-d[i]**2 for i,j in [(0,1),(0,2),(1,2)]])
        np.testing.assert_allclose(He.imag,-1.5/R**2*Hd.imag,atol=2e-18)
        # Use one real orientation in the colored mass rows to make the tree determinant positive.
        mass=np.vstack((Y,h_d.conj()))
        mass[-1]*=-1
        det=np.linalg.det(mass)
        assert det.real>0 and abs(det.imag)<1e-15
        out[str(sign)]={"signed_weak_quartet_in_fixed_convention":float(j),
                        "tree_colored_determinant":float(det.real),
                        "heavy_lepton_to_color":float(np.linalg.norm(h_e)/np.linalg.norm(h_d))}
    assert abs(out["1"]["signed_weak_quartet_in_fixed_convention"])>1e-12
    assert out["0"]["signed_weak_quartet_in_fixed_convention"]==0
    np.testing.assert_allclose(out["1"]["signed_weak_quartet_in_fixed_convention"],
                               -out["-1"]["signed_weak_quartet_in_fixed_convention"])
    return out


def gauge(g):
    a1,a2,a3=4*math.pi/np.square(g)
    b1,b2,b3=41/10,-19/6,-7
    t=((a1-a3)-.6*(a2-a3))/((b1-b3)-.6*(b2-b3))
    delta=(a2-a3)-(b2-b3)*t
    length=2*math.pi*t
    required_log=2*math.pi*delta/(2/3)
    return [length,required_log,required_log-length,math.exp(-required_log)]


def calc():
    data=load("inputs.json")
    table=np.array(load(data["matched_spectra"])["scales"][data["scale"]])
    d,e=np.split(table[:,0],2)
    R=float(np.prod(d/e))
    corners=[table[:,0]+data["widths"]*table[:,1]*z for z in itertools.product((-1,1),repeat=6)]
    rv=[float(np.prod(row[:3]/row[3:])) for row in corners]
    rmax=max(rv)
    threshold=lambda r:2*math.sqrt((9-4*r*r)*(1-r*r))/(6+4*r*r)
    tau_max=max(float(np.prod(row[:3])/np.prod(row[3:5])) for row in corners)
    couplings=load(data["gauge_inputs"])
    g=np.array([couplings[x] for x in ("g1_GUT","g2","g3")])
    corner_g=[g[:,0]+data["widths"]*g[:,1]*z for z in itertools.product((-1,1),repeat=3)]
    # Logs and log differences are affine in the 3 inverse coupling values.
    ranges=np.array([gauge(z) for z in corner_g])
    cfg=data["parent_example"]
    return {"R_center":R,"R_two_width":[min(rv),rmax],
            "alignment_times_cosine_needed_center":threshold(R),
            "alignment_times_cosine_needed_two_width":threshold(rmax),
            "max_acute_phase_angle_two_width_deg":math.degrees(math.acos(threshold(rmax))),
            "five_data_quadrature_tau_ceiling":tau_max,
            "tau_observed_lower":float(table[5,0]-data["widths"]*table[5,1]),
            "independent_per_mass_fraction_threshold":(1-rmax**(1/3))/(1+rmax**(1/3)),
            "flavor_implied_dirac_pair_threshold_delta23":[-math.log(rmax)/(3*math.pi),-math.log(min(rv))/(3*math.pi)],
            "minimal_gauge_test_logGUT_over_MZ_logsplit_excess_Rrequired":{
                "center":gauge(g[:,0]),"two_width_min":ranges.min(axis=0).tolist(),
                "two_width_max":ranges.max(axis=0).tolist()},
            "all_six_central_spectra_examples":witness(d,e,math.radians(cfg["phi_degrees"]),cfg["orthogonal_angles_radians"])}


def verify(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:
            verify(a[k],b[k])
    else:
        np.testing.assert_allclose(a,b,rtol=2e-7,atol=1e-26)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write-new",action="store_true")
    mode.add_argument("--check",action="store_true")
    opt=parser.parse_args()
    result=calc()
    path=HERE/"results.json"
    if opt.write_new:
        with path.open("x",encoding="utf8") as f:
            json.dump(result,f,indent=2)
            f.write("\n")
    else:
        verify(result,load("results.json"))
    print("Verified light spectra, colored determinant, phase constraint and one-loop gauge test.")
