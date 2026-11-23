"""Round 1459: W polarization sums, full-width envelope, and a two-tensor control."""
import argparse
import json
import math
from pathlib import Path


ROOT=Path(__file__).resolve().parent
MW=80.379
FERMI=1.1663787e-5
HBARC=1.973269804e-13


def helicity_data(mass):
    n=mass**2/(2*MW**2)-1
    if n<=1:
        raise ValueError("two real on-shell W required")
    beta=math.sqrt(1-4*(MW/mass)**2)
    g2=4*math.sqrt(2)*FERMI*MW**2
    # For unit v_delta, tree metric tensor coefficient in this convention.
    a=math.sqrt(2)*g2
    kww=beta*(2+n*n)*a*a/(32*math.pi*mass)
    k=kww*mass/(16*math.pi)
    inv=(n*n+1)/(n*n-1)
    b=-2*n/(n*n-1)
    T,L=inv+n*b,n*inv+b
    assert math.isclose(T,-1,abs_tol=5e-14) and math.isclose(L,n,rel_tol=2e-14)
    return k,kww,{"F_longitudinal_pair":n*n/(2+n*n),
                   "max_metric_amplitude_factor_at_identical_helicity_yields":inv,
                   "relative_second_tensor_at_reflected_branch":b/inv}


def compute():
    out={}
    for mass in [200,500,1000]:
        k,kww,pol=helicity_data(mass)
        kll=mass/(16*math.pi)
        rows={}
        for name,atm in [("NO",(2.509-.050)*1e-3),("IO",2*(2.482-.050)*1e-3)]:
            trace=atm+(7.388-.078*2.05)*1e-5
            optimal_vev=(kll*trace*1e-18/kww)**.25
            ll=kll*trace*1e-18/optimal_vev**2
            ww=kww*optimal_vev**2
            assert abs(ll-ww)<ll*1e-13
            for factor in [.2,.5,1,2,5]:
                assert ll/factor**2+ww*factor**2 >= (ll+ww)*(1-1e-13)
            rows[name]={"longest_mm_at_zero_other_widths":HBARC/(ll+ww),
                        "vev_GeV_at_lifetime_envelope":optimal_vev}
        out[str(mass)]={"K_from_explicit_polarizations":k,"polarization":pol,
                        "orders":rows}
    return out


def equivalent(a,b):
    if isinstance(a,dict) and isinstance(b,dict):
        return a.keys()==b.keys() and all(equivalent(a[k],b[k]) for k in a)
    return math.isclose(a,b,rel_tol=2e-12,abs_tol=1e-18)


if __name__=="__main__":
    cmd=argparse.ArgumentParser()
    mode=cmd.add_mutually_exclusive_group()
    mode.add_argument("--write",action="store_true")
    mode.add_argument("--check",action="store_true")
    args=cmd.parse_args()
    path=ROOT/"independent_results.json"
    output=compute()
    if args.write:
        with path.open("x",encoding="utf-8") as stream:
            stream.write(json.dumps(output,indent=2,allow_nan=False)+"\n")
    else:
        assert equivalent(output,json.loads(path.read_text(encoding="utf-8")))
    print("PASS identical-particle gauge width, lifetime envelope and polarized competitor")
