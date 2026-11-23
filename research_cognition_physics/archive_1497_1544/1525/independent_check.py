"""Does the tauonic B_c partial width itself exceed the measured total width?"""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
FERMI = 1.1663787e-5
HBAR = 6.582119569e-25
TAU = 1.77686
# Units GeV and fs. Central, followed by marginal widths added linearly here.
DATA = [(6.27447,.00027,.00017), (.434,.015), (513.4,11.,5.7), (4.578,.008)]
BANDS = {"R_D":(.44,.02), "R_Dstar":(.75,.04)}
MASS_FLOOR = 80.


def sm_fraction(mass, decay_constant, lifetime, vcb):
    phase = (1-(TAU/mass)**2)**2
    return (lifetime*1e-15/HBAR * FERMI**2*vcb**2*decay_constant**2
            *mass*TAU**2*phase/(8*math.pi))


def width(point, x, inv_h2):
    m,f,t,ratio,vcb = point
    b = sm_fraction(m,f,t,vcb)
    rho = ratio/(1+ratio)
    amplitude = 1-m*m*(rho*x*x-(1-rho)*inv_h2)/(1-m*m*inv_h2)
    return b*amplitude**2


def upper_x(point, inv_h2):
    m,f,t,ratio,vcb=point
    b=sm_fraction(m,f,t,vcb)
    rho=ratio/(1+ratio)
    return math.sqrt(((1+1/math.sqrt(b))*(1-m*m*inv_h2)/m**2
                      +(1-rho)*inv_h2)/rho)


def compute():
    vcb=json.loads((HERE/"../1412/inputs.json").read_text(encoding="utf8"))["BelleII_Vcb_BGL"]
    records=[*DATA,vcb]
    centre=[entry[0] for entry in records]
    box=[(v[0]-2*sum(v[1:]),v[0]+2*sum(v[1:])) for v in records]
    rmax=1/MASS_FLOOR**2
    x_limit=max(upper_x(p,a) for p in itertools.product(*box) for a in (0,rmax))
    cases={}
    for name,(x,err) in BANDS.items():
        z=(x-2*err,x+2*err)
        # These particular bands are wholly past the amplitude zero.
        assert z[0] > .30
        samples=[width(p,xx,a) for p in itertools.product(*box)
                 for xx in z for a in (0,rmax)]
        cases[name]={"central_SM_normalization_contact":width(centre,x,0),
                     "two_width_total_fraction_range":[min(samples),max(samples)],
                     "extra_amplitude_at_least_in_SM_units_if_rate_exceeds_total":
                         max(0,math.sqrt(min(samples))-1)/math.sqrt(sm_fraction(
                             box[0][1],box[1][1],box[2][1],box[4][1])),
                     "SM_normalization_multiplier_needed_at_most_if_excess":
                         min(1,1/math.sqrt(min(samples)))}
    assert cases["R_Dstar"]["two_width_total_fraction_range"][0]>1
    assert cases["R_D"]["two_width_total_fraction_range"][0]<1
    for a in (0,rmax):
        np.testing.assert_allclose(width(centre,upper_x(centre,a),a),1.,rtol=1e-13)
    return {"SM_Bc_to_tau_fraction":sm_fraction(*centre[:3],centre[4]),
            "x_ceiling_central_contact_GeV_inv":upper_x(centre,0),
            "x_ceiling_over_two_widths_GeV_inv":x_limit,
            "M_H_floor_for_this_test_GeV":MASS_FLOOR,"historical_branches":cases}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    flags=parser.add_mutually_exclusive_group(required=True)
    flags.add_argument("--write-new",action="store_true")
    flags.add_argument("--check",action="store_true")
    args=parser.parse_args()
    value=compute()
    dest=HERE/"independent_results.json"
    if args.write_new:
        with dest.open("x",encoding="utf8") as f:
            json.dump(value,f,indent=2)
            f.write("\n")
    else:
        expected=json.loads(dest.read_text(encoding="utf8"))
        def compare(a,b):
            if isinstance(a,dict):
                assert a.keys()==b.keys()
                for k in a:
                    compare(a[k],b[k])
            else:
                np.testing.assert_allclose(a,b,rtol=1e-11,atol=1e-15)
        compare(value,expected)
    print("Verified independent total-width bound and the two historical tau branches.")
