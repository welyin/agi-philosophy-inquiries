"""Center characters and a chirality-resolved finite-mass Pauli form factor."""
import argparse
import itertools
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
READ = lambda path: json.loads((HERE/path).read_text(encoding="utf-8"))
cfg = READ("inputs.json")
parent = READ(cfg["weak_benchmarks"])["by_E_over_L"]
pole = READ(cfg["earlier_family"])["pole_MeV"]
masses = [pole["electron"][0]/1000, pole["muon"][0]/1000,
          READ(cfg["earlier_family_result"])["roots_MeV"][1]/1000]


def center_kernel():
    # H and Phi first force the hypercharge angle to be a sixth root,
    # and the central U(3)_F action to be trivial, respectively.
    fields = cfg["left_weyl_table_triality_doublet_y6_family_fundamental_charge"]
    return [[a,b,n] for a,b,n in itertools.product(range(3),range(2),range(6))
            if all((2*a*t+3*b*s+n*y)%6 == 0 for t,s,y,f in fields.values())]


def dipole(i, base, g2, step=0.125):
    mi = masses[i]
    components = []
    for j,mj in enumerate(masses):
        vec2 = base**2*(mi+mj)/(2*masses[0])
        lam2, u = mi**2/vec2, mj/mi
        def integrand(t):
            gap = math.exp(-t)
            x = -math.expm1(-t)
            den = gap*(1-lam2*x) + lam2*u*u*x
            if i == j:
                numerator = 2*x*gap*(x-4)-4*lam2*x**3
            else:
                numerator = 2*x*gap*(x-2)+lam2*x*x*(gap-(1+x)*u*u)
            return gap*numerator/den
        count = round(72/step)
        integral = step/3*math.fsum((1 if n in (0,count) else 4 if n%2 else 2)
                                    *integrand(n*step) for n in range(count+1))
        # diag axial g^2/2; off diagonal g^2/4 for each, averaged in numerator
        components.append(g2*lam2*integral/(16*math.pi**2))
    return components


def invert_electron(target, g2):
    assert target < 0
    a,b = 1.,2000.
    value = lambda m: sum(dipole(0,m,g2))
    assert value(a)<target<value(b)
    for _ in range(40):
        mid=(a+b)/2
        if value(mid)<target:
            a=mid
        else:
            b=mid
    return (a+b)/2


def rectangular(offset,width,n):
    trap,trap_error=163,15
    lo=offset-trap-n*(width+trap_error)
    hi=offset-trap+n*(width+trap_error)
    vals=[t*s*1e-9 for t in (lo,hi) for s in (.0000084,.00000852)]
    return [min(vals),max(vals)]


def calculate():
    ratio=math.sqrt(masses[2]/masses[0])
    modern={k:[rectangular(offset,width,n) for n in (1,2)]
            for k,(offset,width) in {"Cs":(45,27),"Rb":(205.2,9.7)}.items()}
    earlier=READ(cfg["recoil"])["recoil_ae_rectangles_one_two"]
    scenarios={}
    for k,b in parent.items():
        m=b["scale_GeV"]/ratio
        g2=b["g2_squared"]
        e=dipole(0,m,g2,.03125)
        mu=dipole(1,m,g2,.03125)
        for i,fine in enumerate((e,mu)):
            coarse=dipole(i,m,g2,.0625)
            assert all(abs(x/y-1)<5e-6 for x,y in zip(coarse,fine))
        evalue=sum(e)
        assert max(e)<0 and max(mu)<0
        scenarios[k]={"M11_ceiling_GeV":m,"gF_squared":g2,
                      "electron_terms":e,"muon_terms":mu,"electron_sum":evalue,
                      "muon_sum":sum(mu),"muon_over_electron":sum(mu)/evalue,
                      "Rb_two_width_min_extra_positive":modern["Rb"][1][0]-evalue}
    g2=parent["1"]["g2_squared"]
    windows={}
    for period,data in (("2023",earlier),("2025",modern)):
        windows[period]={}
        for atom,groups in data.items():
            windows[period][atom]={}
            for n,(lo,hi) in enumerate(groups,1):
                bounds=([invert_electron(lo,g2),invert_electron(hi,g2) if hi<0 else None]
                        if lo<0 else [])
                windows[period][atom][str(n)]=bounds
    kernel=center_kernel()
    assert kernel==[[0,0,0],[0,1,3],[1,0,4],[1,1,1],[2,0,2],[2,1,5]]
    return {"center_kernel_a_b_n":kernel,"allowed_subgroup_orders":[1,2,3,6],
            "modern_recoil_ae_rectangles":modern,"light_pole_upper_edge":scenarios,
            "M11_windows_GeV_middle_coupling":windows}


def rounded(obj):
    if isinstance(obj,dict):return {key:rounded(v) for key,v in obj.items()}
    if isinstance(obj,list):return [rounded(v) for v in obj]
    return float(format(obj,".9g")) if isinstance(obj,float) else obj


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    group=parser.add_mutually_exclusive_group()
    group.add_argument("--check",action="store_true")
    group.add_argument("--write",action="store_true")
    args=parser.parse_args()
    path=HERE/"results.json"
    result=rounded(calculate())
    if args.write:
        with path.open("x",encoding="utf-8") as f:
            json.dump(result,f,indent=2)
            f.write("\n")
    else:
        assert result==json.loads(path.read_bytes())
        print("PASS 1624: central subgroup, exact paired vector loops and recoil windows")
