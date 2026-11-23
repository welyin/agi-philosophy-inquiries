"""CMS Table 3, five-flavour transport and conditional one-triplet discrimination."""
import argparse
import hashlib
import json
import math
from pathlib import Path


HERE=Path(__file__).resolve().parent
MZ=91.1876
# (bin edges, representative Q, alpha at Z, alpha at Q, total alpha(Q) errors -/+)
DATA=[(74,220,103.06,.1182,.1160,.0012,.0014),
      (220,395,266.63,.1184,.1019,.0009,.0008),
      (395,638,464.31,.1179,.0947,.0008,.0009),
      (638,1410,753.66,.1184,.0898,.0007,.0008),
      (1410,3103,1600.5,.1170,.0821,.0008,.0010)]


def rhs(inv):
    # Conventional coefficients of beta(a), a=alpha/(4*pi), MSbar at nf=5.
    nf=5
    z3,z4,z5=1.2020569031595943,math.pi**4/90,1.0369277551433699
    b0=11-2*nf/3
    b1=102-38*nf/3
    b2=2857/2-5033*nf/18+325*nf**2/54
    b3=(149753/6+3564*z3+nf*(-1078361/162-6508*z3/27)
        +nf**2*(50065/162+6472*z3/81)+1093*nf**3/729)
    b4=(8157455/16+621885*z3/2-88209*z4/2-288090*z5
        +nf*(-336460813/1944-4811164*z3/81+33935*z4/6+1358995*z5/27)
        +nf**2*(25960913/1944+698531*z3/81-10526*z4/9-381760*z5/81)
        +nf**3*(-630559/5832-48722*z3/243+1618*z4/27+460*z5/9)
        +nf**4*(1205/2916-152*z3/81))
    a=1/(4*math.pi*inv)
    return (b0+a*(b1+a*(b2+a*(b3+a*b4))))/(2*math.pi)


def evolve(inv, log_ratio, steps=8):
    h=log_ratio/steps
    for _ in range(steps):
        k1=rhs(inv)
        k2=rhs(inv+h*k1/2)
        k3=rhs(inv+h*k2/2)
        k4=rhs(inv+h*k3)
        inv+=h*(k1+2*k2+2*k3+k4)/6
    return inv


def boxes(z,k,mass,binned,start,transport_error=0.):
    out=[]
    for lo,hi,q,_,alpha,down,up in DATA[start:]:
        central=max(0.,math.log(q/mass))
        t=[central+float(binned)*(max(0.,math.log(energy/mass))-central) for energy in (lo,hi)]
        l,h=sorted(z*v/(2*math.pi) for v in t)
        ylo=1/(alpha+k*up)-transport_error+l
        yhi=1/(alpha-k*down)+transport_error+h
        out.append((evolve(ylo,math.log(MZ/q)),evolve(yhi,math.log(MZ/q))))
    return out


def profile(z,k,mass,binned,start,err=0.):
    intervals=boxes(z,k,mass,binned,start,err)
    lower=max(x[0] for x in intervals)
    upper=min(x[1] for x in intervals)
    return dict(common_inverse_aZ_interval=[lower,upper], feasible=lower<=upper)


def transition(f,lo,hi):
    before=f(lo)
    assert before!=f(hi)
    for _ in range(38):
        mid=(lo+hi)/2
        if f(mid)==before:
            lo=mid
        else:
            hi=mid
    return (lo+hi)/2


def compute():
    reproduction=[1/evolve(1/aZ,math.log(q/MZ))-aQ for _,_,q,aZ,aQ,_,_ in DATA]
    assert max(map(abs,reproduction))<.00011
    step=max(abs(evolve(1/aQ,math.log(MZ/q),8)-evolve(1/aQ,math.log(MZ/q),16))
             for _,_,q,_,aQ,_,_ in DATA)
    assert step<1e-6
    scans={}
    for start,name in [(0,"all_five"),(1,"four_above_top")]:
        mass=DATA[start][0]
        for binned in [False,True]:
            p=lambda z,k,err=0.: profile(z,k,mass,binned,start,err)
            kcrit=transition(lambda k:p(2/3,k)["feasible"],0.,5.)
            rows={}
            for k in [1,2]:
                assert p(0,k)["feasible"]
                cap=transition(lambda z:p(z,k)["feasible"],0.,6.)
                test=p(2/3,k)
                need=0. if test["feasible"] else transition(lambda err:p(2/3,k,err)["feasible"],0.,2.)
                rows[str(k)]=dict(upper_z=cap,one_triplet=test,extra_inverse_alpha_room_to_restore=need)
            scans[name+("_bin_support" if binned else "_representative")]=dict(minimum_k_for_triplet=kcrit,widths=rows)
    kappa=transition(lambda t:profile(2/3,1,220,t,1)["feasible"],0.,1.)
    z=2/3
    low,high=profile(z,1,220,True,1)["common_inverse_aZ_interval"]
    anchor=(low+high)/2
    witness=[]
    for lo,hi,q,_,alpha,down,up in DATA[1:]:
        base=evolve(anchor,math.log(q/MZ))
        target=max(math.log(lo/220),min(math.log(hi/220),(base-1/alpha)*2*math.pi/z))
        effective=220*math.exp(target)
        observable=1/(base-z*target/(2*math.pi))
        assert lo-1e-10<=effective<=hi+1e-10 and alpha-down-1e-10<=observable<=alpha+up+1e-10
        witness.append(dict(Q=effective, alpha=observable))
    return dict(CMS_Q_GeV=[d[2] for d in DATA],reproduction_difference_in_alpha_Q=reproduction,
                max_inverse_alpha_step_doubling_difference=step, profiles=scans,
                high_four_fractional_log_bin_freedom_for_reentry=kappa,
                high_four_support_witness=dict(alphaZ=1/anchor, effective_response=witness),
                threshold_examples_high_four={str(m): profile(2/3,1,m,False,1) for m in [220,400,600,1000,2000]},
                highest_pair_expected_inverse_alpha_shift=(2/3)*math.log(1600.5/266.63)/(2*math.pi))


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    check=p.parse_args().check
    result=json.dumps(compute(),indent=2,allow_nan=False)+"\n"
    dest=HERE/"independent_results.json"
    if check:
        assert json.loads(dest.read_bytes())==json.loads(result)
    else:
        dest.write_text(result,encoding="utf-8")
    print("1783 CMS NNLO jet slope discrimination",hashlib.sha256(result.encode()).hexdigest()[:12])
