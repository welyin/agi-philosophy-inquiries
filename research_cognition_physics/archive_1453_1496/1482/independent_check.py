"""Check the common muon normalization in a second set of actual channels."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

BASE = Path(__file__).resolve().parent


def read(path):
    return json.loads((BASE/path).read_text(encoding="utf-8"))


def audit():
    cfg = read("inputs.json")
    osc = read(cfg["oscillations"])
    ck = read("../1381/inputs.json")
    prior = read("../1383/results.json")
    angle, spread = osc["DayaBay_amplitude"]
    zlo,zhi = [(1-math.sqrt(1-i))/2 for i in (angle-2*spread,angle+2*spread)]
    results = {}
    for order in ("NO", "IO"):
        solar, e = osc[order]["solar_sine_squared"]
        b, elo, *upper = osc[order]["atmospheric_sine_squared"]
        ehi = upper[0] if upper else elo
        if order == "NO":
            ze = (1-solar-2*e)*(1-zhi)
            zm = max(0,math.sqrt((solar-2*e)*(1-b-2*ehi))-
                       math.sqrt((1-solar+2*e)*(b+2*ehi)*zhi))**2
        else:
            ze,zm = zlo,(b-2*elo)*(1-zhi)
        tau_max = 1-ze-zm
        quadratic = (1-ze)*(1-zm)
        bound = (cfg["MEG_II"]["BR_90CL_upper"]*24*math.pi*
                 cfg["inverse_alpha_upper"]/cfg["matched_dipole_rate_lower_factor"])
        # Solve x_e x_mu * s² <= bound * D_mu using rigorous outer extrema.
        a = ze*zm/bound-quadratic
        strength = (tau_max+math.sqrt(tau_max**2+a))/a
        mu_rate = 1+2*tau_max*strength+quadratic*strength**2
        for chi in (0.,1.4,3.13):
            t13 = math.asin(math.sqrt((zlo+zhi)/2))
            t12,t23 = math.asin(math.sqrt(solar)),math.asin(math.sqrt(b))
            def rot(i,j,ang,phase=0.):
                r=np.eye(3,dtype=complex)
                r[i,i]=r[j,j]=math.cos(ang)
                r[i,j]=math.sin(ang)*np.exp(-1j*phase)
                r[j,i]=-r[i,j].conjugate()
                return r
            u=rot(1,2,t23)@rot(0,2,t13,chi)@rot(0,1,t12)
            null=u[:,0 if order=="NO" else 2]
            f=np.column_stack([np.cross(unit,null) for unit in np.eye(3)])
            np.testing.assert_allclose(f@null,0,atol=2e-15)
            M=f@np.diag([1.,3.,2.])@f.T
            np.testing.assert_allclose(M@null,0,atol=2e-15)
            x=abs(null)**2
            def decay(i,j,s):
                sm=np.zeros((3,3),complex)
                sm[i,j]=1
                amp=sm-s*np.outer(f[j].conj(),f[i])
                return np.linalg.norm(amp)**2
            np.testing.assert_allclose(decay(1,0,.7),1+1.4*x[2]+.49*(1-x[0])*(1-x[1]))
            assert (decay(2,1,.7)>decay(2,0,.7))==(order=="NO")
        results[order]={"contact_upper":strength,"Gmu_over_G0_upper":math.sqrt(mu_rate),
                        "semileptonic_row_lower":1/mu_rate}
    observed={}
    for width in (1,2):
        # Same upper/lower exclusive |Vub| conventions as round 1383; use its
        # archived complete working rectangle with the common K average.
        observed[str(width)]=prior["bounds"][str(width)]["official_mixed"]["row"]
    assert min(block["semileptonic_row_lower"] for block in results.values())>observed["1"][1]
    assert observed["2"][1]>=min(block["semileptonic_row_lower"] for block in results.values())
    a,ea=ck["ud_beta"]
    f,ef=ck["K_form"]
    v,ev=ck["K_product"]
    upper_2term=(a+ea)**2+((v+ev)/(f-ef))**2
    assert 0 < observed["1"][1]-upper_2term < .0001
    return {"model":results,"archived_inferred_row":observed,
            "minimum_additional_row_mismatch_one_width":
                {k:v["semileptonic_row_lower"]-observed["1"][1] for k,v in results.items()}}

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
