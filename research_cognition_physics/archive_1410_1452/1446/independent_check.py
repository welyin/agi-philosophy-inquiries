"""Takagi-family distance and Dirac alternatives, using published JUNO inputs."""

import argparse
import json
from math import sqrt
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def ratio(matrix):
    s=np.linalg.svd(matrix,compute_uv=False)
    return s[1]/s[0]


def results():
    t=sqrt(2)-1
    r=t*t
    # Original signed inverted Delta31 errors reverse when taking its magnitude.
    Slo,Shi=(7.388-2*.078)*1e-5,(7.388+2*.078)*1e-5
    Alo,Ahi=(2.509-2*.025)*1e-3,(2.509+2*.027)*1e-3
    io_Alo=(2.482-2*.025)*1e-3
    qlo,qhi=sqrt(Slo/Ahi),sqrt(Shi/Alo)
    def right_orthogonal_left_overlap(q):
        return sqrt(1-((t+1/t)/(q+1/q))**2)
    hcrit=sqrt(((qhi+1/qhi)-(t+1/t))/((qhi+1/qhi)+2))
    left=np.array([1.,0.,0.])
    def second(h):
        return np.array([h,sqrt(1-h*h),0.])
    both=np.outer(left,left)+t*np.outer(second(hcrit),second(hcrit))
    assert abs(ratio(both)-qhi)<1e-13
    for q in (qlo,r,qhi):
        h=right_orthogonal_left_overlap(q)
        k = np.column_stack([left,t*second(h),np.zeros(3)])
        assert abs(ratio(k)-q)<1e-13
    assert abs(right_orthogonal_left_overlap(r)-sqrt(7)/3)<1e-13
    low = sqrt(io_Alo)
    high = sqrt(io_Alo+Shi)
    distance=(low-r*high)/(1+r)
    target=np.diag([high,low,0])
    amplitude=(low+high)/(1+r)
    witness=np.diag([amplitude,r*amplitude,0])
    assert abs(np.linalg.norm(target-witness,ord=2)-distance)<1e-14
    # Check phase dependence for genuinely complex physical source columns.
    v=np.array([1,2j,1])/sqrt(6)
    w=np.array([1,-1j,.5j])/1.5
    inner=np.vdot(v,w)
    for phase in np.linspace(0,2*np.pi,20):
        c=np.exp(1j*phase)
        q=ratio(np.outer(v,v)+c*r*np.outer(w,w))
        result=(r+1/r+2*(c*inner**2).real)/(1-abs(inner)**2)
        assert abs(q+1/q-result)<1e-12
    return {"JUNO_NO_q_range":[round(qlo,10),round(qhi,10)],
            "distance_to_IO_meV_all_masses":round(1000*distance,8),
            "relative_to_heavy_at_nearest_IO":round(distance/high,9),
            "orthogonal_right_Dirac_needed_left_overlap_range":[round(right_orthogonal_left_overlap(qhi),9),
                                                                round(right_orthogonal_left_overlap(qlo),9)],
            "either_side_Dirac_sharp_required_max_overlap":round(hcrit,9),
            "Dirac_left_overlap_for_exact_squared_seed":round(sqrt(7)/3,9)}


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    mode=cli.add_mutually_exclusive_group()
    mode.add_argument("--write",action="store_true")
    mode.add_argument("--check",action="store_true")
    args=cli.parse_args()
    value=results()
    dest=HERE/"independent_results.json"
    if args.write:
        with dest.open("x",encoding="utf8") as stream:
            json.dump(value,stream,indent=2)
            stream.write("\n")
    else:
        assert value==json.loads(dest.read_text(encoding="utf8"))
    print(json.dumps(value))
