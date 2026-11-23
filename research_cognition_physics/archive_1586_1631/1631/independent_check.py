"""Independent two-triplet width and lifetime audit using complex Yukawa witnesses."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
OUT = HERE/"independent_results.json"


def light_masses(order, abs31, solar_ratio):
    if order == "NO":
        return np.sqrt([abs31*solar_ratio, abs31])
    return np.sqrt([abs31, abs31*(1+solar_ratio)])


def partial_widths_neutral(mass, m_eff_ev, v, mass_ports):
    # Equation (9) specifies per conjugate channel coefficients.
    norm_y_sq = 2 * mass * (m_eff_ev*1e-9) / v**2
    width_unbroken = norm_y_sq*mass/(8*math.pi)
    flavored = []
    for name, branching in (("W",1/4),("Z",1/8),("H",1/8)):
        r = (mass_ports[name]/mass)**2
        phase = (1-r)**2 * (1+2*r if name != "H" else 1)
        for _ in range(2):
            flavored.append(width_unbroken*branching*phase if r<1 else 0.)
    return math.fsum(flavored)


def check_witnesses(masses_ev, m_light, m_heavy, v):
    # Exercise both attainable endpoint assignments using the actual 2x2 seesaw.
    a,b = masses_ev
    cap = 1e9*v**2/(2*m_heavy)
    y = math.acosh((2*cap+b-a)/(a+b))/2
    for z, reference in ((0.,a),(math.pi/2+1j*y,cap+b-a)):
        rot=np.array([[np.cos(z),-np.sin(z)],[np.sin(z),np.cos(z)]])
        md=np.diag(np.sqrt(1e-9*masses_ev))@rot@np.diag(np.sqrt([m_light,m_heavy]))
        reconstructed=md@np.diag([1/m_light,1/m_heavy])@md.T
        assert np.allclose(reconstructed*1e9,np.diag(masses_ev),rtol=3e-9,atol=3e-11)
        norms=2*np.sum(np.abs(md)**2,axis=0)/v**2
        assert norms[1] <= 1+1e-12
        mt_ev=1e9*np.linalg.norm(md[:,0])**2/m_light
        assert math.isclose(mt_ev,reference,rel_tol=4e-13)


def calculate():
    boxes=json.loads((STAGE/"1625/inputs.json").read_bytes())["box"]
    p_bounds=json.loads((STAGE/"1611/results.json").read_bytes())["whole_triplet_geometric_mean_GeV"]
    ref=json.loads((STAGE/"1610/inputs.json").read_bytes())
    v,ports,hbarc=ref["v_GeV"],ref["kinematic_GeV"],ref["hbarc_GeV_mm"]
    results={}
    for order in ["NO","IO"]:
        masses=[light_masses(order,d,r) for d,r in
                itertools.product(boxes["abs31_eV2"],boxes["solar_over_abs31"])]
        floor=min(x[0] for x in masses)
        biggest=max(x[1] for x in masses)
        delta=max(x[1]-x[0] for x in masses)
        for ml in [200,1000]:
            for p,m_pair in itertools.product(p_bounds,masses):
                check_witnesses(m_pair,ml,p*p/ml,v)
            mh_min=p_bounds[0]**2/ml
            ceiling=v**2/(2*mh_min)*1e9+delta
            width=lambda mt: partial_widths_neutral(ml,mt,v,ports)
            pref=hbarc/width(1.0)
            # Extra competition: a real orthogonal R bounds each effective mass by m_b.
            # If a future proper lifetime is 1 micron, this much imaginary angle is necessary.
            sample_ev = pref/.001
            required_cosh = min((2*sample_ev-(b-a))/(a+b) for a,b in masses)
            imag_bound = math.acosh(max(1,required_cosh))/2
            results[f"{order}_{ml}"]={
                "mixing_norm_squared":[floor*1e-9/ml,ceiling*1e-9/ml],
                "width_GeV":[width(floor),width(ceiling)],
                "proper_length_mm":[pref/ceiling,pref/floor],
                "minimum_proper_length_real_R_mm":pref/biggest,
                "min_abs_imag_angle_for_1um_example":imag_bound}
    return results


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=calculate()
    if args.check:
        assert result==json.loads(OUT.read_bytes())
    else:
        with OUT.open("x",encoding="utf-8") as out:
            json.dump(result,out,indent=2)
            out.write("\n")
    print(json.dumps(result,indent=2))
