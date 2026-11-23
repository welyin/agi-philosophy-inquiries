"""Basis-independent absolute-flavour CP source and a finite splitting necessity bound."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

BASE=Path(__file__).resolve().parent


def loops(ratio):
    x=ratio**2
    f=ratio*(1/(1-x)+1-(1+x)*math.log1p(1/x))
    return f,1/(1-x)


def source(y,i):
    j=1-i
    norms=y.conj().T@y
    hij=norms[i,j]
    hii,hjj=norms[i,i].real,norms[j,j].real
    return hij,hii,hjj,np.outer(y[:,j],y[:,i].conj())


def sharp(y,mass,i):
    j=1-i
    f,g=loops(mass[j]/mass[i])
    hij,hii,hjj,outer=source(y,i)
    z=f*hij+g*hij.conjugate()
    op=z*outer
    op=(op-op.conj().T)/(2j)
    direct=np.linalg.eigvalsh(op)
    val=math.sqrt((z*hij).imag**2+abs(z)**2*(hii*hjj-abs(hij)**2))
    assert np.isclose(sum(abs(direct)),val,rtol=1e-10,atol=1e-24)
    norm=8*math.pi*hii
    return val/norm,abs(np.trace(op))/norm,sum(abs(np.diag(op)))/norm


def audit():
    cfg=json.loads((BASE.parent/"1484"/"inputs.json").read_text(encoding="utf-8"))
    sol=cfg["KamLAND_dm21_eV2"]
    atm=cfg["T2K_dm32_NO_eV2"]
    m2=math.sqrt(sol[0])
    m3=math.sqrt(sol[0]+atm[0])
    sth=math.sqrt(m2*m2+m3*m3+10*m2*m3)
    product=11*m2*m3/3
    small=(sth-math.sqrt(sth*sth-4*product))/2
    big=sth-small
    u=246.22/math.sqrt(2)
    conversion=7.04*(28/79)*135*1.202056903159594/(4*math.pi**4*106.75)
    rho=4/math.sqrt(22)
    sample={}
    for ta,tb in [(small,big),(big,small)]:
        m0=(8*math.pi**2*u*u*125.11**2/((ta+27*tb)*1.e-9))**(1/3)
        M=np.array([m0,3*m0])
        y=np.column_stack([np.array([0,1,1])*math.sqrt(ta*1.e-9*m0/2)/u,
                           np.array([1,3,1])*math.sqrt(tb*1.e-9*3*m0/11)*np.exp(1j*math.pi/3)/u])
        eps=[sharp(y,M,i) for i in (0,1)]
        sample["t_a_low" if ta==small else "t_a_high"]={
            "mass_low_GeV":m0,"mass_high_GeV":3*m0,"sum_individually_sharp_eta_cap":
                conversion*sum(e[0] for e in eps),"single_source":[list(map(float,e)) for e in eps]}
    upper_sol=sol[0]+2*sol[2]
    m2=math.sqrt(upper_sol)
    m3=math.sqrt(upper_sol+atm[0]+2*atm[2])
    t_upper=2e-9*math.sqrt(m2*m2+m3*m3+10*m2*m3)
    # Each positive threshold transfer >=1/2, independently broadened t_i<=2 t_i0.
    # Use deliberately rounded-up Higgs upper 126 GeV; kinematic mass scheme after matching.
    weighted=(16*math.pi**2*u*u*126**2)**(1/3)*t_upper**(2/3)
    pref=conversion*rho*weighted/(8*math.pi*u*u)
    target=5.8e-10
    def cap(r):
        return pref*(r/(r-1)+2/3)
    crit=1+1/(target/pref-5/3)
    assert cap(3)<target and np.isclose(cap(crit),target,rtol=1e-11,atol=1e-20)
    for r in np.geomspace(1.001,30,90):
        f,g=loops(r)
        assert r*(abs(f)+abs(g))<=r/(r-1)+.5
        f,g=loops(1/r)
        assert abs(f)+abs(g)<=r/(r-1)+2/3
    return {"max_light_times_sum_t_GeV2":weighted,"yield_to_photon":conversion,
            "conservative_eta_cap_r_ge_3":cap(3),"target_eta":target,
            "necessary_ratio_upper":crit,"central_r3_examples":sample}


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
