"""Independent majorization bound for the whole solar/atmospheric box."""
import argparse
import json
import math
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent


def check():
    read=lambda p: json.loads(p.read_text(encoding="utf-8"))
    high=read(HERE.parent/"1402"/"inputs.json")["scales"]["1e16"]
    masses=np.array(high)
    down=np.linalg.norm((masses[:,0]+2*masses[:,1])[:3])
    gap=masses[-1,0]-2*masses[-1,1]-(masses[-2,0]+2*masses[-2,1])
    # Loose, common enclosing intervals, deliberately avoid fitting delta or using correlations.
    A=(.28,.335); B=(.020,.024); C=(.398,.626)
    def fractions(a,b,c):
        lx=math.sqrt(a[0]*c[0]); ux=math.sqrt(a[1]*c[1])
        vy=math.sqrt((1-a[0])*(1-c[0])*b[1])
        return (max(0,lx-vy)**2,(ux+vy)**2)
    box=[fractions(A,B,C),fractions((1-A[1],1-A[0]),B,C),
         ((1-C[1])*(1-B[1]),(1-C[0])*(1-B[0]))]
    assert min(t[0] for t in box)>.05 and max(t[1] for t in box)<.6
    P=1-(down/gap)**2
    base=.6**2+.35**2+.05**2
    # increasing top, decreasing bottom by TV; middle fixed for the maximizing sorted candidate.
    t=(-1.1+math.sqrt(1.1**2+8*(P-base)))/4
    assert 0<t<.05
    return {"probability_independent_envelope":box,"max_tau_purity_coarse_majorization":base,
       "required_tau_purity_min":P,"required_extra_matrix_Frobenius_coarse":gap*math.sqrt(1-base)-down,
       "required_probability_transport_TV":t,"required_Takagi_row_angle_degrees":math.degrees(math.asin(t))}


if __name__=="__main__":
    p=argparse.ArgumentParser()
    c=p.add_mutually_exclusive_group()
    c.add_argument("--write",action="store_true");c.add_argument("--check",action="store_true")
    opt=p.parse_args()
    out=json.dumps(check(),indent=2)+"\n"; path=HERE/"independent_results.json"
    if opt.write:
        with path.open("x",encoding="utf-8") as f: f.write(out)
    elif out!=path.read_text(encoding="utf-8"): raise SystemExit("1453 independent mismatch")
    else:print("1453 independent OK")
