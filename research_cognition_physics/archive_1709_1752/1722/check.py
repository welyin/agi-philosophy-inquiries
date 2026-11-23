"""Integer certificates for the selection obstruction and a finite weak-path budget."""
import argparse
from itertools import product
import json
import math
from pathlib import Path
import runpy

import numpy as np


HERE=Path(__file__).resolve().parent
def read(f):
    return json.loads((HERE/f).read_bytes())
CFG=read("inputs.json")
FIELDS="s x a2 t1 t2 t3 f2 f3 hu hd n2 h0 h6 b l2 lx li w".split()

def word(*fields):
    return np.array([fields.count(field)-(field=="w") for field in FIELDS],dtype=int)

def identities():
    i=word("s","s","s","s","x")-word()
    R=word("t2","f3","hd","s")
    D=word("t2","f3","hd","s","s","s","s","s","x")
    v0=word("li")
    v1=word("li","s","s","s","s","x")
    A,E=word("t1","f2","hd"),word("t3","f2","hd","s")
    c,E0=word("t2","f2","hd"),word("t3","t3","hu","s","s")
    nu=word("n2","h6","x")-word("n2","h0","a2")
    term=word("b","h6","h6")-word("b","h0","h0")
    lock=word("lx","x","x","x")-word("lx")
    a2=word("l2","a2")-word("l2","s","s")
    cert={"two_named_weak_paths":D-R,"linear_I_vacuum_lock":v1-v0,
          "upper_u_top_and_A_E":word("t1","t1","hu")+i-E0-2*(A-E),
          "common_neutral_terminal_and_real_X":lock-2*nu+term-2*a2,
          "charm_top_and_C_E_twice":word("t2","t2","hu")+2*i-E0-2*(c-E)}
    for label,row in cert.items():
        np.testing.assert_array_equal(row,(2 if "twice" in label else 1)*i)
    return {name:{f:int(k) for f,k in zip(FIELDS,row) if k} for name,row in cert.items()}

def weak_data():
    old=runpy.run_path(str(HERE/CFG["charged_run"]))
    _,rep=old["charged_background"](.1)
    D=old["unpack"](rep["down_texture_10PeV"])
    f=read(CFG["flavor"])
    amp=f["scalar_magnitude_over_cutoff"]
    assert amp==f["real_X_over_cutoff"]==.5
    angle=np.pi*f["scalar_angle_over_pi"]
    omega=np.exp(4j*angle)
    cD=D[1,2].imag/omega.imag
    cR=D[1,2].real-cD*omega.real
    edges=np.array([D[0,1].real,D[1,0].real,D[1,1].real,D[2,1].real,D[2,2].real,cR,cD])
    base=np.array(read(CFG["down_paths"])["powers_A_B_C_E_F_R_D"])
    coeff=edges/amp**(base+[0,0,0,0,0,0,1])
    return base,coeff,angle

def example(ex,base,coeff,angle):
    N,s,x=ex["group"],ex["S"],ex["X"]
    t,f=np.array(ex["T"]),np.array(ex["F"])
    n=np.array(ex["n"])
    iq=(4*s+x)%N
    assert iq==1 and 3*x%N==0
    # The last weak channel has the same bare charge as R, with its first I still explicit.
    bare=base.copy()
    bare[-1]-=4
    m=n%3
    power=bare+4*n
    ports=[(0,1),(1,0),(1,1),(2,1),(2,2),(1,2),(1,2)]
    arms=np.array([t[i]+f[j]+ex["Hd"] for i,j in ports])
    assert np.all((arms+s*power+x*m)%N==0)
    assert np.all((2*t+ex["Hu"]+s*np.array([0,0,2])+np.array([1,2,0])*iq)%N==0)
    assert (2*t[0]+ex["Hu"])%N and (2*t[1]+ex["Hu"])%N and (2*t[1]+ex["Hu"]+iq)%N
    factor=2.**(power+m-base-np.array([0,0,0,0,0,0,1]))
    old=coeff*0.5**(base+[0,0,0,0,0,0,1])*np.exp(1j*base*angle)
    new=coeff*factor*.5**(power+m)*np.exp(1j*power*angle)
    np.testing.assert_allclose(new,old,atol=1e-14)
    direct=np.array([[0,new[0],0],[new[1],new[2],new[5]+new[6]],[0,new[3],new[4]]])
    fphase=np.angle(np.linalg.det(direct)*np.exp(1j*(2+12+6)*angle))
    # Down determinant -ABF with the raw real signs; the heavy determinant is +s^6.
    assert abs(fphase)<1e-12
    earlier=[list(range(int(v)%N,int(v),N)) for v in n]
    return {"cyclic_order":N,"edge_wilsons":(coeff*factor).tolist(),"maximum_increase":max(factor),
            "max_abs_Wilson":max(abs(coeff*factor)),"max_change_named_amplitudes":max(abs(new-old)),
            "earlier_dressings_in_same_charge_class":earlier,"tree_strong_angle_rad":fphase}

def budget(coeff):
    best=2**24
    possibilities=[]
    for a,r in product(range(0,7,3),range(1,8,3)):
        diff=abs(a-r)
        for order in range(3,diff+1):
            if diff%order or math.gcd(order,6)!=1:continue
            cost=max(2**(4*a+a%3),2**(4*r+r%3-5))
            possibilities.append([a,r,order,cost])
    assert min(v[-1] for v in possibilities)==best
    return {"weak_path_relative_floor":best,"actual_weak_path_absolute_floor":min(abs(coeff[5:]))*best,
            "minimal_cases":[v for v in possibilities if v[-1]==best]}

def calculate():
    base,coeff,angle=weak_data()
    rows=[example(e,base,coeff,angle) for e in CFG["examples"]]
    b=budget(coeff)
    # With maximal budget at or below the first exhibit, every individual
    # basic path has n <= 6, and every long path has n <= 4.
    ceiling=rows[0]["max_abs_Wilson"]
    nmax=[]
    for i,c in enumerate(abs(coeff)):
        start=1 if i==6 else 0
        feasible=[n for n in range(start,19,3) if c*2.**(4*n+(n%3)-5*start)<=ceiling*(1+1e-10)]
        nmax.append(feasible)
    lower=[]
    for E,C,R,F,D in product(nmax[3],nmax[2],nmax[5],nmax[4],nmax[6]):
        d=abs(D-R)
        for order in range(3,d+1):
            if math.gcd(order,6)!=1 or d%order:continue
            if (C-E-1)%order==0 and (R-F-1)%order==0:
                lower.append(max(abs(coeff[2])*2.**(4*C),abs(coeff[3])*2.**(4*E),
                                 abs(coeff[4])*2.**(4*F),abs(coeff[5])*2.**(4*R),
                                 abs(coeff[6])*2.**(4*D+D%3-5)))
    assert min(lower)>=ceiling*(1-1e-10)
    return {"integer_row_certificates":identities(),"original_named_Wilsons":coeff.tolist(),
            **b,"two_relaxed_weak_examples":rows,"optimal_max_abs_Wilson_relaxed_charged_graph":ceiling,
            "allowed_n_at_or_below_optimal_max":nmax}
if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=calculate()
    dest=HERE/"results.json"
    if args.check:
        assert json.loads(dest.read_bytes())==result
        print("1722 exact character certificates and optimal charged-path budget verified")
    else:
        dest.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf-8")
