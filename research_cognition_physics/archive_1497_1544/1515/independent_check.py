"""Check a quark/right-lepton common CP relation using full projected resolvents."""
import argparse
import json
from pathlib import Path

import numpy as np

BASE=Path(__file__).resolve().parent


def cycle(h):
    return (h[0,1]*h[1,2]*h[2,0]).imag


def effective(y,h,z):
    v=y@h.conj()
    mass2=np.vdot(h,h).real
    return y@y.T-np.outer(v,v.conj())/(mass2-z)


def from_resolvent(full,z):
    squared=full@full.conj().T
    resolvent=np.linalg.inv(z*np.eye(4)-squared)[:3,:3]
    return z*np.eye(3)-np.linalg.inv(resolvent)


def audit():
    y=np.array([[1,1,0,0],[0,1,1,0],[1,0,2,0]],float)
    m=np.array([1,-1,2,4])
    ell=np.array([1,2,-1,0])
    phase=(3+4j)/5
    hu=np.diag([1.,4.,11.])
    rot=np.diag([1.,1j,-1j])
    output={}
    for scale in (1,10,100):
        hd,he=[scale*(m+c*phase*ell) for c in (2,-3)]
        md2,me2=[float(np.vdot(h,h).real) for h in (hd,he)]
        full_d,full_e=[np.vstack([y,h]) for h in (hd,he)]
        svd,sve=[np.linalg.svd(a,compute_uv=False) for a in (full_d,full_e)]
        dets=[np.linalg.det(a) for a in (full_d,full_e)]
        for det in dets:
            assert abs(det-12*scale)<1.e-8 and det.real>0
        r_phys=np.prod(svd[1:])/np.prod(sve[1:])
        assert abs(r_phys-sve[0]/svd[0])<1.e-11
        normy2=np.linalg.norm(y,ord=2)**2
        for bare,mass in ((md2,svd[0]),(me2,sve[0])):
            assert bare<=mass**2+2.e-10 and mass**2<=bare+normy2+2.e-10
        w_true=1.5*md2/me2
        w_poles=1.5/r_phys**2
        ratio_range=[1-normy2/svd[0]**2,1/(1-normy2/sve[0]**2)]
        assert ratio_range[0]<=w_true/w_poles<=ratio_range[1]
        d0,e0=[effective(y,h,0) for h in (hd,he)]
        r_light=np.sqrt(np.linalg.det(d0).real/np.linalg.det(e0).real)
        assert abs(1.5/r_light**2-w_true)<2e-13
        assert abs(cycle(d0))>1.e-3
        combo=e0+w_true*d0
        assert np.max(np.abs(combo.imag))<2.e-12
        cam=rot@e0@rot.conj().T+w_true*d0
        assert abs(cycle(cam))>1.e-3
        quasi=cam@hu-hu@cam
        assert abs(np.trace(quasi@quasi@quasi).imag)>1
        # Keep the left lepton-doublet matrix and all charged-lepton masses;
        # change only the physical right-lepton orientation relative to the 10 source.
        rot4=np.eye(4,dtype=complex)
        rot4[:3,:3]=rot
        modified=rot4@full_e
        assert np.allclose(modified.conj().T@modified,full_e.conj().T@full_e)
        assert abs(np.linalg.det(modified)-dets[1])<1e-8
        offshell=[]
        for zd,ze in ((-1.,-2.),(.2,-.5)):
            d,e=from_resolvent(full_d,zd),from_resolvent(full_e,ze)
            assert np.allclose(d,effective(y,hd,zd)) and np.allclose(e,effective(y,he,ze))
            wz=1.5*(md2-zd)/(me2-ze)
            residual=np.max(np.abs((e+wz*d).imag))
            assert residual<1e-10
            offshell.append({"z_down":zd,"z_lepton":ze,"weight":wz,"imaginary_residual":float(residual)})
        output[str(scale)]={"equal_full_det":float(dets[0].real),"full_light_product_ratio":float(r_phys),
                            "effective_ratio":r_light,"effective_weight":w_true,
                            "physical_mass_inferred_weight":w_poles,"weight_ratio_enclosure":ratio_range,
                            "down_CP_cycle":float(cycle(d0)), "combined_cycle":float(cycle(combo)),
                            "changed_right_lepton_cycle":float(cycle(cam)),
                            "changed_commutator_trace_imag":float(np.trace(quasi@quasi@quasi).imag),
                            "resolvents":offshell}
    return output


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    m = p.add_mutually_exclusive_group()
    m.add_argument("--write-new", action="store_true")
    m.add_argument("--check", action="store_true")
    args = p.parse_args()
    path = BASE/"independent_results.json"
    out = json.loads(json.dumps(audit(), default=float))
    if args.write_new:
        with path.open("x", encoding="utf-8") as f:
            f.write(json.dumps(out, indent=2)+"\n")
    elif args.check:
        assert out == json.loads(path.read_text(encoding="utf-8"))
        print("Independent tests passed.")
    else:
        print(json.dumps(out, indent=2))
