"""Enumerate all charged irreps that saturate the anomaly cost; use pairwise running slopes."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent


def active_reps():
    reps=[]
    for length in range(1,7):
        for first in range(-4,5):
            # A charge-conjugate spectrum has the same running. Retain nonnegative Y.
            charges=range(first,first+length)
            twice_y=2*first+length-1
            cost=sum(q*q for q in charges)
            if twice_y>=0 and 0<cost<=10:
                twice_index=length*(length*length-1)//6
                reps.append((length,twice_y,cost,twice_index))
    return reps


def enumerate_spectra(reps):
    groups={}
    def visit(start,remaining,index,spectrum):
        if remaining==0:
            groups.setdefault(index,[]).append(spectrum)
        for i in range(start,len(reps)):
            n,y,cost,twice=reps[i]
            if cost<=remaining:
                visit(i,remaining-cost,index+twice,spectrum+[(n,y)])
    visit(0,10,0,[])
    return groups


def crossing(g,two_index):
    a=[8*math.pi**2/v**2 for v in g]
    b12,b23=109/15,23/6
    index=two_index/2
    shift12=.8*(10-index)-(4/3)*index
    shift23=(4/3)*index
    det=b12*shift23-b23*shift12
    if abs(det)<1e-10:
        raise ArithmeticError("parallel slopes")
    L=((a[0]-a[1])*shift23-(a[1]-a[2])*shift12)/det
    run=(b12*(a[1]-a[2])-b23*(a[0]-a[1]))/det
    inv_common=(a[2]+7*L)/(2*math.pi)
    return L,run,inv_common


def hermitian_running(log_window,steps):
    low=json.loads((HERE/"../1386/inputs.json").read_text(encoding="utf-8"))
    s12,s13,s23,phase=[v[0] for v in low["quark_sines_and_phase"]]
    def rotation(i,j,s,phase=0.):
        transform=np.eye(3,dtype=complex)
        transform[i,i]=transform[j,j]=math.sqrt(1-s*s)
        transform[i,j]=s*np.exp(-1j*phase)
        transform[j,i]=-transform[i,j].conjugate()
        return transform
    v=rotation(1,2,s23)@rotation(0,2,s13,phase)@rotation(0,1,s12)
    state=np.array([np.diag([x[0]**2 for x in low["yukawa"][key]])
                    for key in ("up","down","charged_leptons")],dtype=complex)
    state[1]=v@state[1]@v.conj().T
    gauge0=np.array([.461228,.65096,1.2123])
    beta=np.array([4.1,-19/6,-7.])
    gauge_weights=np.array([[.85,2.25,8],[.25,2.25,8],[2.25,2.25,0.]])
    def derivative(t,h):
        g2=1/(1/gauge0**2-t*beta/(8*math.pi**2))
        tr=np.trace(3*h[0]+3*h[1]+h[2]).real
        run=(tr-gauge_weights@g2)[:,None,None]*np.eye(3)
        run=run+1.5*np.array([h[0]-h[1],h[1]-h[0],h[2]])
        return (run@h+h@run)/(16*math.pi**2)
    step=log_window/steps
    for i in range(steps):
        midpoint=state+.5*step*derivative(i*step,state)
        state+=step*derivative((i+.5)*step,midpoint)
    masses=np.sqrt(np.linalg.eigvalsh(state))
    coeff=masses/.22**np.array([[8,4,1],[8,6,3],[8,5,3]])
    return [float(coeff.min()),float(coeff.max())]


def calculate():
    reps=active_reps()
    assert all(n<=5 for n,_,_,_ in reps)
    groups=enumerate_spectra(reps)
    g0=(.461228,.65096,1.2123)
    e=(.000026,.00004,.0046)
    corners=[tuple(x+2*s*d for x,d,s in zip(g0,e,sign)) for sign in itertools.product((-1,1),repeat=3)]
    rows={}
    for twice,spectra in sorted(groups.items()):
        L,run,alpha=crossing(g0,twice)
        points=[crossing(g,twice) for g in corners]
        valid=[r>0 and l-r>0 and a>0 for l,r,a in points]
        # L-run, run and alpha are affine in all three inverse squared couplings.
        mins=[min(p[i] for p in points) for i in range(3)]
        maxs=[max(p[i] for p in points) for i in range(3)]
        gap=[min(l-r for l,r,a in points),max(l-r for l,r,a in points)]
        all_ok=mins[1]>0 and gap[0]>0 and mins[2]>0
        all_bad=maxs[1]<0 or gap[1]<0 or maxs[2]<0
        if not all_ok and not all_bad:
            assert any(valid) and not all(valid)
        assert all(valid) == all_ok
        rows[str(twice/2)]={"spectrum_count":len(spectra),"all_two_width_box_physical":all_ok,
                            "physical_corner_count":sum(valid),"log_M_over_MZ_full_range":gap}
        if 0<run<L and alpha>0:
            m=91.1876*math.exp(L-run)
            mx=91.1876*math.exp(L)
            rows[str(twice/2)].update({"M_central":m,"MX_central":mx,
                 "ma_max_micro_eV_central": math.sqrt(2*math.pi)*30*.00570e15/m,
                 "flavor_cutoff_above_MX_if_ma_below_micro_eV":30*.00570e15/(math.sqrt(2)*.22*mx)})
    log_window=math.log(rows["10.0"]["M_central"]/91.1876)
    coeff=hermitian_running(log_window,1600)
    assert abs(coeff[0]-.502804)<.00001 and abs(coeff[1]-2.195453)<.00003
    examples=[[(-2,-1,0,1,2)],[(-1,0,1)]*5]
    fourth=[sum(q**4 for slot in slots for q in slot) for slots in examples]
    assert fourth==[34,10]
    return {"active_irreps_length_twiceY_cost_twiceIndex":reps,"quintet_five_triplets_sum_Q4":fourth,
            "spectrum_count":sum(map(len,groups.values())),"T_table":rows,"joint_flavor_coef_RK2":coeff}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--write-new",action="store_true")
    parser.add_argument("--check",action="store_true")
    options=parser.parse_args()
    out=json.loads(json.dumps(calculate()))
    path=HERE/"independent_results.json"
    if options.write_new:
        with path.open("x",encoding="utf-8") as f:
            json.dump(out,f,indent=2)
            f.write("\n")
    else:
        assert out==json.loads(path.read_text(encoding="utf-8"))
        print("1472 independent OK")
