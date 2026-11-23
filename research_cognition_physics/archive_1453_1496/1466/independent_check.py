"""WW fiducial transfer from analytic spin-one cap matrices and rigorous midpoint errors."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


ROOT=Path(__file__).resolve().parent
CUT=math.cos(math.radians(20))
EYE=np.eye(3)
JZ=np.diag([1.,0.,-1.])
JX=np.array([[0.,1.,0.],[1.,0.,1.],[0.,1.,0.]])/math.sqrt(2)
PAIRS=[(a,b) for a in (1,0,-1) for b in (1,0,-1) if abs(a-b)<=1]
N_MID=24000


def efficiency(theta,beta,which):
    """Decay matrix for lepton from W- (which=1) or W+ (-1)."""
    gamma=1/math.sqrt(1-beta**2)
    sign_helicity=-which
    eff=np.broadcast_to(EYE,(len(theta),3,3)).copy()
    for beam in (-1,1):
        z=beam*which*np.cos(theta)
        vx=-beam*which*np.sin(theta)
        vz=gamma*(z-beta*CUT)
        h=gamma*(CUT-beta*z)
        length=np.hypot(vx,vz)
        d=h/length
        a=(1-d)/2
        p_aligned=1-(1-a)**3
        p_opposite=a**3
        p_zero=3*a*a-2*a**3
        spin=(vx/length)[:,None,None]*JX+(vz/length)[:,None,None]*JZ
        spin2=spin@spin
        eff-=p_zero[:,None,None]*(EYE-spin2)
        eff-=((p_aligned+p_opposite)/2)[:,None,None]*spin2
        eff-=((p_aligned-p_opposite)*sign_helicity/2)[:,None,None]*spin
    return eff


def grams(zlo,zhi,beta):
    lo,hi=math.acos(zhi),math.acos(zlo)
    step=(hi-lo)/N_MID
    accum=np.zeros((2,7,7))
    for start in range(0,N_MID,512):
        theta=lo+(np.arange(start,min(start+512,N_MID))+.5)*step
        em,ep=efficiency(theta,beta,1),efficiency(theta,beta,-1)
        block=np.zeros((len(theta),7,7))
        for i,(a,b) in enumerate(PAIRS):
            for j,(aa,bb) in enumerate(PAIRS):
                if b==bb:
                    block[:,i,j]+=.5*em[:,1-a,1-aa]
                if a==aa:
                    block[:,i,j]+=.5*ep[:,1-b,1-bb]
        for k,m in enumerate((1,-1)):
            d=np.array([(1+m*(a-b)*np.cos(theta))/2 if a!=b
                        else -m*np.sin(theta)/math.sqrt(2) for a,b in PAIRS]).T
            accum[k]+=np.einsum("t,ti,tij,tj->ij",np.sin(theta),d,block,d)*step
    # Matrix-valued integrand Lipschitz <= 40 on 0<=beta<=.65 (review).
    remainder=40*(hi-lo)**2/(4*N_MID)+2e-12
    return accum,remainder


def quotient(F,C):
    chol=np.linalg.cholesky(C)
    S=np.linalg.solve(chol,F)
    S=np.linalg.solve(chol,S.T).T
    return np.linalg.eigvalsh((S+S.T)/2)[-1]


def check():
    rows={}
    for energy, front,middle in [(198.4,(6.98,.47,.04),(2.41,.25,.01)),
                                  (205.9,(7.23,.44,.05),(2.91,.26,.02))]:
        beta=math.sqrt(1-(2*80.379/energy)**2)
        front_m,df=grams(.8,1.,beta)
        mid_m,dc=grams(.4,.6,beta)
        nom=max(quotient(f,c) for f,c in zip(front_m,mid_m))
        upper=max(quotient(f+df*np.eye(7),c-dc*np.eye(7)) for f,c in zip(front_m,mid_m))
        lo=front[0]-2*(front[1]+front[2])
        hi=middle[0]+2*(middle[1]+middle[2])
        ratio=lo/hi
        amp=(max(0,math.sqrt(ratio)-math.sqrt(upper))**2/((1+upper)*(1+ratio)))
        rows[str(energy)]={"fixed_s_accepted_ratio":round(nom,9),
                           "guaranteed_upper_after_quadrature_error":round(upper,9)+1e-9,
                           "two_width_ratio_floor":round(ratio,9),
                           "extra_amplitude_squared_over_observed_two_bins_floor":round(amp,12)}
    # No angular cut: all three positive Wigner-square shapes; largest is (1+z)^2.
    primitive=lambda z:z+z*z+z**3/3
    bare=(primitive(1)-primitive(.8))/(primitive(.6)-primitive(.4))
    return {"on_shell_MW_GeV":80.379,"no_cut_ratio_upper":round(bare,10),"bins":rows}


if __name__=="__main__":
    parse=argparse.ArgumentParser()
    mode=parse.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check",action="store_true")
    mode.add_argument("--write-new",action="store_true")
    args=parse.parse_args()
    res=check()
    if args.write_new:
        with (ROOT/"independent_results.json").open("x",encoding="utf-8") as f:
            json.dump(res,f,indent=2)
            f.write("\n")
    else:
        assert res==json.loads((ROOT/"independent_results.json").read_text(encoding="utf-8"))
    print(json.dumps(res))
