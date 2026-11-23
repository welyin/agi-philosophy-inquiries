"""A common exclusion on both real roots, the full two-width JUNO box and all d >= 1."""
import argparse
import json
import math
import runpy
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ARCH = HERE.parent


def read(name):
    return json.loads((HERE/name).read_text(encoding="utf-8"))


def compute():
    flavor = runpy.run_path(str(ARCH/"1678/independent_check.py"))
    data = read("../1397/inputs.json")["high_scale"]
    q, dq = np.array(data["down"] + data["angles_theta12_theta13_theta23_delta"][:3]).T
    phi = -np.pi/6
    angle_window = np.array(read("../1659/results.json")["joint_windows"])
    osc = read("../1668/inputs.json")
    l0,lm,lp = osc["juno_NO_dm21"]
    h0,hm,hp = osc["juno_NO_dm31"]
    Rbox = [math.sqrt((h0-2*hm)/(l0+2*lp)), math.sqrt((h0+2*hp)/(l0-2*lm))]
    assert 5.70 < Rbox[0] < Rbox[1] < 5.96
    # All constants in the certificate have deliberately rounded outer bounds.
    cmax = 7.38*(5.96+1/5.96) + 2*4.28**2/math.sqrt(13)
    rfloor = 2*2.23**2/(cmax+math.sqrt(cmax**2-4*(2.23*11.23)**2))
    phase_floor = -0.065/math.sqrt(13) - 0.006
    numerator = 11.23*1.21*rfloor + 2*4.28*phase_floor - 7.38*0.31/5.70
    pfloor = numerator/(7.38*(5.96-1/5.96))
    afloor = 4*pfloor*(1-pfloor)
    assert afloor > 0.099 and afloor > angle_window[0,1]
    rows = []
    for w in read("../1682/independent_results.json")["witnesses"]:
        sh,lg,_,_ = flavor["get_matrices"](q + dq*np.array(w["quark_marginal_pulls"]), w["gamma_matching"])
        ye = (sh+w["mag_t"]*np.exp(1j*phi)*lg).T * w["charged_normalization"]
        _,e = np.linalg.eigh(ye@ye.conj().T)
        h = w["neutral_modulus"]
        l = np.array([1, np.exp(2j*phi), h*np.exp(3j*phi)*((3+2*np.exp(4j*phi))/math.sqrt(7))])
        a = e.T@(l*[0,1,1])
        b = e.T@(l*[1,3,1])
        aa,bb,g = np.vdot(a,a).real,np.vdot(b,b).real,np.vdot(a,b).real
        G = aa*bb-g**2
        null = e.conj().T@(np.array([-2,1,-1])/l)
        null /= np.linalg.norm(null)
        pspan = 1-abs(null[0])**2
        z = a[0].conjugate()*b[0]
        assert 2.23 < aa and 11.23 < bb and g < 4.28 and 0 < G < 7.38
        assert 0 < z.real < .065 and abs(z.imag) < .006 and abs(b[0])**2 > 1.21 and pspan < .31
        cases = []
        for d in [0,1]:
            R = w["R"]
            phase = (d-2+2j*math.sqrt(3))/abs(d-2+2j*math.sqrt(3))
            c = G*(R+1/R)-2*g*g*phase.real
            rr = sorted(np.roots([bb*bb,-c,aa*aa]))
            for ratio in rr:
                mass = np.outer(a,a)+ratio*phase*np.outer(b,b)
                vals, modes = np.linalg.eigh(mass.conj().T@mass)
                assert np.linalg.norm(mass@null) < 2e-14
                p = abs(modes)**2
                angles = np.array([4*p[0,2]*(1-p[0,2]),p[0,1]/(1-p[0,2]),p[1,2]/(1-p[0,2])])
                inside = np.all((angle_window[:,0]<=angles)&(angles<=angle_window[:,1]))
                p_formula = ((aa*abs(a[0])**2/ratio+bb*abs(b[0])**2*ratio+2*g*(z*phase).real)
                             /G-pspan/R)/(R-1/R)
                np.testing.assert_allclose([vals[2]/vals[1],p[0,2]], [R*R,p_formula],rtol=3e-12)
                cases.append([d,round(float(ratio),9),round(float(angles[0]),9),bool(inside)])
        assert [x[-1] for x in cases] == [True,False,False,False]
        rows.append({"frozen_R":w["R"],"electron_kernel": [round(float(i),9) for i in
                     [aa,bb,g,G,pspan,abs(b[0])**2,z.real,z.imag]],"d_root_reactor_pass":cases})
    return {"JUNO_two_width_R_box": Rbox, "certificate_enclosing_R": [5.70,5.96],
            "all_d_at_least_1_positive_root_lower_bound": round(rfloor,10),
            "electron_probability_lower_bound": round(pfloor,10),
            "reactor_amplitude_lower_bound": round(afloor,10), "frozen_frame_checks":rows}


if __name__ == "__main__":
    cli=argparse.ArgumentParser()
    cli.add_argument("--check",action="store_true")
    args=cli.parse_args()
    out=compute()
    if args.check:
        assert out == read("independent_results.json")
        print("PASS: full two-width frequency ratio and all d >= 1, both frozen ports")
    else:
        print(json.dumps(out,indent=2,allow_nan=False))
