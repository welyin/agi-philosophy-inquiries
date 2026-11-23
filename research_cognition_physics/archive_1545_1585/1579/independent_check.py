"""Local phase test using six running masses, three angles and a finite logdet remainder."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def load(path):
    return json.loads((HERE/path).read_text(encoding="utf-8-sig"))


def box(entry):
    mean, width = entry
    return mean-2*width, mean+2*width


def ckm_matrix(s12,s13,s23,phase):
    c12,c13,c23 = (math.sqrt(1-x*x) for x in (s12,s13,s23))
    z = np.exp(1j*phase)
    return np.array([[c12*c13,s12*c13,s13/z],
                     [-s12*c23-c12*s23*s13*z,c12*c23-s12*s23*s13*z,s23*c13],
                     [s12*s23-c12*c23*s13*z,-c12*s23-s12*c23*s13*z,c23*c13]])


def run():
    raw = load("../1386/inputs.json")
    up, down = ([box(x) for x in raw["yukawa"][key]] for key in ("up","down"))
    angles = raw["quark_sines_and_phase"]
    a, b, c = [box(x)[1] for x in angles[:3]]
    # Independent outer bounds for any CKM phase, including unmeasured alternatives.
    vmax = [[1,a,b],[a+c*b,1,c],[a*c+b,c+a*b,1]]
    tb_abs2_floor = (1-c*c)*(1-b*b)
    dom = tb_abs2_floor*(up[2][0]/down[2][1]-down[2][1]/up[2][0])
    rest = 0
    for i,j in itertools.product(range(3),repeat=2):
        if i==j==2:
            continue
        ratios = [up[i][0]/down[j][1], up[i][1]/down[j][0]]
        rest += vmax[i][j]**2*max(abs(r-1/r) for r in ratios)
    # Round down the certified gap to leave room for floating arithmetic.
    gap = math.floor(100*(dom-rest))/100
    assert gap > 0
    t = 0.02
    slope = gap*t/(1-t*t)
    K_u = t*(down[-1][1]+t*up[-1][1])/((1-t*t)*up[0][0])
    K_d = t*(up[-1][1]+t*down[-1][1])/((1-t*t)*down[0][0])
    cost = 1.5*(K_u+K_d)+3*(K_u*K_u+K_d*K_d)
    radius = min(.5/max(K_u,K_d),slope/(2*cost))
    edm = load("../1389/inputs.json")["edm"]
    a0,err1,err2 = edm["abs_slope_e_fm"]
    theta = edm["limit_90_ecm"]/(edm["fm_in_cm"]*(a0-2*(err1+err2)))
    alpha_limit = 2*theta/slope
    assert alpha_limit < radius
    # Verify the local identity for all 64 signed spectra in actual complex flavour bases.
    v = ckm_matrix(*(x[0] for x in angles))
    u0,d0 = [np.array([x[0] for x in raw["yukawa"][label]]) for label in ("up","down")]
    signed_min = math.inf
    for signs in itertools.product((-1,1),repeat=6):
        u = u0*np.array(signs[:3])
        d = d0*np.array(signs[3:])
        U=np.diag(u)
        D=v@np.diag(d)@v.conj().T
        A=(U-t*D)/(1-t*t)
        B=(D-t*U)/(1-t*t)
        matrix_derivative = (t*np.trace(np.linalg.solve(D,A)-np.linalg.solve(U,B))).real
        generated = float(np.sum(abs(v)**2*(u[:,None]/d-d/u[:,None])))*t/(1-t*t)
        assert abs(matrix_derivative-generated) < 1e-9
        signed_min = min(signed_min,abs(generated)*(1-t*t)/t)
    return {"signed_branches_checked":64, "all_phase_tb_term_floor":dom,
            "other_eight_terms_absolute_ceiling":rest,"trace_abs_certified_floor":gap,
            "central_minimum_abs_trace_all_signs":signed_min,"example_t":t,
            "absolute_slope_min":slope, "certified_radius_rad":radius,
            "neutron_theta_upper_in_declared_matching":theta,
            "maximum_abs_vacuum_phase_inside_radius_rad":alpha_limit}


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--check",action="store_true")
    group.add_argument("--write-new",action="store_true")
    args=parser.parse_args()
    result=run()
    out=HERE/"independent_results.json"
    if args.write_new:
        with out.open("x",encoding="utf-8") as f:
            json.dump(result,f,indent=2)
            f.write("\n")
    else:
        assert result == json.loads(out.read_text())
    print(json.dumps(result))
