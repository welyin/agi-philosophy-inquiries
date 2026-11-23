"""An independent exact one-channel test for the rank-one splitting."""
import argparse
import json
import math
from itertools import product
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def overlap_floor(s,b,mu,tau):
    assert s < mu < b < tau
    return (mu**2-s**2)*(tau**2-b**2)/((b**2-s**2)*(tau**2-mu**2))


def complete_right(d,e,w,k):
    """Make two rows identical by a right isometry, using the left neutral plane."""
    first = np.diag(d)
    second = (w*e[None,:])/k
    vals, cols = np.linalg.eigh(first@first.conj().T-second@second.conj().T)
    scale = max(np.linalg.norm(first)**2,np.linalg.norm(second)**2)
    assert vals[0]<0<vals[-1] and abs(vals[1])<2e-11*scale
    neutral = (cols[:,0]*math.sqrt(vals[-1])+cols[:,-1]*math.sqrt(-vals[0]))/math.sqrt(vals[-1]-vals[0])
    common = np.stack((neutral.conj(),cols[:,1].conj()))
    first_rows, second_rows = common@first, common@second
    left, _, right = np.linalg.svd(second_rows.conj().T@first_rows)
    residual = first-second@(left@right)
    singular = np.linalg.svd(residual,compute_uv=False)
    assert np.linalg.norm(common@residual)<1e-8*math.sqrt(scale)
    assert singular[1]<1e-8*math.sqrt(scale)
    return singular


def run():
    inputs = json.loads((HERE.parent/"1402"/"inputs.json").read_text(encoding="utf-8-sig"))
    values = np.array(inputs["scales"]["1e16"])
    d,e = np.split(values[:,0],2)
    theta_floor = overlap_floor(d[1],d[2],e[1],e[2])
    rotation = np.eye(3)
    c,s = math.sqrt(1-theta_floor),math.sqrt(theta_floor)
    rotation[1:,1:] = [[c,-s],[s,c]]
    witness = complete_right(d,e,rotation,1.)
    np.testing.assert_allclose(abs(rotation[2,1])**2,theta_floor,rtol=1e-14)
    median = np.linalg.svd(rotation*e[None,:]/d[:,None],compute_uv=False)[1]
    assert abs(median-1.)<1e-12
    floors=[]
    for signs in product([-1.,1.],repeat=4):
        sch,b,mu,tau = (values[[1,2,4,5],0]+2*values[[1,2,4,5],1]*signs)
        floors.append(overlap_floor(sch,b,mu,tau))
    # A cross-check with unrelated spectra and complex W, not a proof by sampling.
    rng=np.random.default_rng(143301)
    for _ in range(32):
        w,_=np.linalg.qr(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))
        point=np.sort(np.exp(rng.uniform(-1,1,6)).reshape(2,3),axis=1)
        first,second=point
        central=np.linalg.svd(w*second[None,:]/first[:,None],compute_uv=False)[1]
        complete_right(first,second,w,central)
    min_prob=min(floors)
    n=math.ceil(math.log(.05)/math.log1p(-min_prob))
    return {"physical_match":"kappa=1, common left current, increasing rank order",
            "central_b_mu_sharp_floor":theta_floor,"sharp_floor_range_two_widths":[min_prob,max(floors)],
            "attaining_weak_basis_probabilities":np.round(abs(rotation)**2,11).tolist(),
            "single_splitting_singular_value":float(witness[0]),
            "null_mu_probability_at_minimum_ideal_trials_95":(1-min_prob)**n,
            "previous_null_mu_probability":(1-min_prob)**(n-1),"required_ideal_b_tags":n}


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group()
    group.add_argument("--check",action="store_true")
    group.add_argument("--write-new",action="store_true")
    args=parser.parse_args()
    result=run()
    target=HERE/"independent_results.json"
    if args.write_new:
        with target.open("x",encoding="utf-8") as handle:
            json.dump(result,handle,indent=2)
            handle.write("\n")
    elif args.check:
        recorded=json.loads(target.read_text(encoding="utf-8"))
        for key,value in result.items():
            if isinstance(value,str):
                assert value==recorded[key]
            else:
                np.testing.assert_allclose(value,recorded[key],rtol=2e-10,atol=1e-16)
        print("1433 independent reconstruction and sharp tagged probability verified")
    else:
        print(json.dumps(result,indent=2))
