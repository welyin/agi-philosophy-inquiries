"""Independent integer-determinant and integrated screening tests for round 1451."""
import argparse
from itertools import product
import json
import math
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent


def read(p):
    return json.loads(p.read_text(encoding="utf-8"))


def quark_log_ranges(data):
    # Ordered MZ or GUT input; no reference to author's log table.
    out=[]
    for key in ("up","down"):
        center,err=np.array(data[key]).T
        out.append([-sum(np.log(center+2*err)), -sum(np.log(center-2*err))])
    return np.array(out)


def joint(data,budget,e):
    bands=quark_log_ranges(data)+3*math.log(budget)*np.array([-1.,1.])
    entries=[]
    for a,b in product(range(1,10),repeat=2):
        if a+b>10: continue
        width=bands/np.array([[a],[b]])
        if max(width[:,0])<=min(width[:,1]):
            entries.append([math.exp(-max(width[:,0])),a,b])
    best=max(entries)
    ceilings=np.ceil(bands[:,0]/math.log(1/e)).astype(int)
    assert np.all(ceilings*math.log(1/e)<=bands[:,1])
    return {"largest_epsilon_mass_volumes_for_D_at_most_10":best[0],
        "optimizing_indices":best[1:], "minimum_indices_at_epsilon":ceilings.tolist(),
        "sum_min":int(sum(ceilings))}


def compare():
    inp=read(HERE/"inputs.json")
    low=read(HERE.parent/"1386"/"inputs.json")["yukawa"]
    high=read(HERE.parent/"1397"/"inputs.json")["high_scale"]
    eps=inp["illustrative_epsilon"]
    data={label:joint(side,inp["coefficient_budget_B"],eps) for label,side in [("MZ",low),("1e16",high)]}
    # m_0 is bare at Phi=0, two diagonal Phi insertions, the reverse step is forbidden.
    m0=np.array([[0.,.5],[0.,0.]])
    coupling=.5*np.eye(2)
    masses=np.linalg.svd(m0+eps*coupling,compute_uv=False)
    assert np.linalg.matrix_rank(m0)==1
    assert np.isclose(np.prod(masses),(eps/2)**2) and np.max(masses)<1
    logweight=(2/3)*sum(-np.log(masses))
    assert logweight>(2/3)*2*math.log(1/eps)
    assert np.allclose(np.linalg.svd(eps*np.eye(2),compute_uv=False),[eps,eps])
    charges=inp["quark_charge_control"]
    # Verify telescoping at the level of the complete multiset of signed U(1) charges.
    for x in sum(charges.values(),[]):
        states=[x]+[k-1 for k in range(1,x+1)]
        conjugate=[k for k in range(1,x+1)]
        assert sorted(states)==[0]+conjugate
    net=sum(2*q+u+d for q,u,d in zip(charges["Q"],charges["uc"],charges["dc"]))
    q=inp["cms_R32"]
    a,b=q["alpha_centers"]; da,db=q["working_symmetric_errors"]
    l,h=q["bins_GeV"]
    perfect=(7-2*net/3)*math.log(h/l)/(2*math.pi)
    radii=[]
    for n in [1,2]:
        candidates=[1/(b+i*n*db)-1/(a+j*n*da) for i,j in product([-1,1],repeat=2)]
        radii.append(max(0,min(candidates)-perfect,perfect-max(candidates)))
    data.update({"test_block_masses_over_cutoff":masses.tolist(),
        "test_block_integrated_Dirac_screening":logweight,"saturated_bound_for_two_units":(4/3)*math.log(1/eps),
        "D27_finite_window_inverse_alpha_shortfall":net*math.log(1/eps)/(3*math.pi),
        "CMS_transfer_required_if_all_partners_active_across_two_bins":radii})
    return data


if __name__=="__main__":
    p=argparse.ArgumentParser()
    grp=p.add_mutually_exclusive_group()
    grp.add_argument("--write",action="store_true"); grp.add_argument("--check",action="store_true")
    cmd=p.parse_args()
    result=json.dumps(compare(),indent=2,allow_nan=False)+"\n"
    path=HERE/"independent_results.json"
    if cmd.write:
        with path.open("x",encoding="utf-8") as f: f.write(result)
    elif not path.read_text(encoding="utf-8")==result:
        raise SystemExit("1451 independent mismatch")
    else: print("1451 independent OK")
