"""Independent auditor of the center/electric identity and its finite-window bound."""
import argparse
from itertools import product
import json
from pathlib import Path

import numpy as np


def norm(a):
    return float(np.linalg.norm(a,ord=2))


def blind_window():
    # Charge in thirds. Equal charges from different tensor words and integer
    # attachments occur with multiplicity. Only the input/output charge is visible.
    charge = np.array([-4,-3,-3,-2,-1,0,0,1,2,2,3,4])
    triality = (-charge) % 3
    z = np.diag(np.exp(2j*np.pi*triality/3))
    allowed = np.equal.outer(charge,charge)
    rotations = np.equal.outer(triality,triality)
    assert np.all(rotations[allowed])
    rng = np.random.default_rng(1356)
    mix=rng.normal(size=allowed.shape)+1j*rng.normal(size=allowed.shape)
    mix *= allowed
    mix /= norm(mix)
    np.testing.assert_allclose(z @ mix, mix @ z, atol=2e-15)
    np.testing.assert_allclose(np.diag(z),np.exp(-2j*np.pi*charge/3),atol=3e-15)
    # A raw two-leg contraction for the vector representation of SO(3).
    d=np.eye(3).flatten()/np.sqrt(3)
    rotation=np.array([[0,-1,0],[1,0,0],[0,0,1]])
    np.testing.assert_allclose(d @ np.kron(rotation,rotation),d)
    # Independently tagged, strictly colorless endpoint, Q=4/3.
    return {"electric_blocks":int(allowed.sum()),"also_center_invariant":int((allowed&rotations).sum()),
            "ideal_pair_delivery_given_delta":float(abs(np.vdot(d,d)/np.vdot(d,d))**2),
            "joint_success_two_independent_mixed_legs":float(np.vdot(d,d).real/(9*np.vdot(d,d).real)),
            "center_mismatch_of_pair_vertex":float(abs(1-np.exp(4j*np.pi/3)))}


def finite_error():
    # Simultaneous imperfections in electric conservation and both endpoint dictionaries.
    rng=np.random.default_rng(51356)
    worst=0.
    for _ in range(96):
        sectors=rng.integers(-7,8,size=(2,5))
        z=[np.diag(np.exp(-2j*np.pi*row/3)) for row in sectors]
        h=[np.diag(row/3+rng.uniform(-.03,.03,5)) for row in sectors]
        phases=[np.diag(np.exp(-2j*np.pi*np.diag(q))) for q in h]
        k=rng.normal(size=(5,5))+1j*rng.normal(size=(5,5)); k/=norm(k)
        err=2*np.pi*norm(h[1]@k-k@h[0])+sum(norm(zk-pk) for zk,pk in zip(z,phases))
        assert norm(z[1]@k-k@z[0]) <= err+1e-13
        for a in range(3):
            for b in range(3):
                if a==b: continue
                proj_in=np.diag(((-sectors[0]) % 3)==a)
                proj_out=np.diag(((-sectors[1]) % 3)==b)
                actual=norm(proj_out @ k @ proj_in)
                assert actual <= err/np.sqrt(3)+1e-13
                worst=max(worst,actual*np.sqrt(3)/err)
    return {"cases":96,"largest_observed_fraction_of_bound":worst}


def graph_quotient():
    # A square plus diagonal 0--2 gives b1=2. Exhaust all discrete lifts of a fixed w.
    edges=[(0,1),(1,2),(2,3),(3,0),(0,2)]
    flat=list(product(range(3),repeat=len(edges)))
    unseen=set(flat)
    loops=set()
    orbit_sizes=[]
    while unseen:
        origin=next(iter(unseen))
        orbit=set()
        for tail in product(range(3),repeat=3):
            gauge=(0,*tail)
            orbit.add(tuple((e+gauge[b]-gauge[a])%3 for e,(a,b) in zip(origin,edges)))
        invariants={((f[0]+f[1]-f[4])%3,(f[4]+f[2]+f[3])%3) for f in orbit}
        assert len(invariants)==1 and not (invariants & loops)
        loops.update(invariants)
        unseen-=orbit
        orbit_sizes.append(len(orbit))
    assert sorted(orbit_sizes)==[27]*9
    # Exchanging w-marginal weights changes visible average, even with all probes locked.
    return {"edge_lifts":len(flat),"vertex_orbit_sizes":sorted(orbit_sizes),"independent_holonomies":len(loops),
            "different_effective_holonomies_locked_bright_ports":
            [float((1+np.cos(t))/2) for t in (0,2*np.pi/3)]}


def optical_holdouts():
    zeta=np.exp(2j*np.pi/3)
    lifts=((1,zeta),(zeta**-1,1))
    def ports(n,t,u,z):
        return [(1+(np.exp(1j*d)*u**n*z**t).real)/2 for d in (0,np.pi/2)]
    for n,t in [(2,1),(-1,1),(3,0)]:
        np.testing.assert_allclose(ports(n,t,*lifts[0]),ports(n,t,*lifts[1]),atol=4e-15)
    a=ports(1,0,*lifts[0]); b=ports(1,0,*lifts[1])
    # |cos(theta+eta)-cos(theta)| <= |eta|; loss of contrast <= .04.
    radius=(.02+.04)/2+.005
    gaps=np.abs(np.array(a)-np.array(b))-2*radius
    assert all(gaps>0)
    return {"colorless_N1_lift_A":a,"colorless_N1_lift_B":b,"conservative_halfwidth":radius,
            "disjoint_interval_gaps":gaps.tolist()}


def sample():
    return {"round":1356,"restricted_ports":blind_window(),"perturbations":finite_error(),
            "gauge_orbits":graph_quotient(),"holdouts":optical_holdouts()}


def compare(actual,wanted):
    if isinstance(actual,dict):
        assert actual.keys()==wanted.keys()
        for k in actual:compare(actual[k],wanted[k])
    else:
        np.testing.assert_allclose(actual,wanted,rtol=1e-12,atol=2e-13)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--check",action="store_true")
    group.add_argument("--write-new",action="store_true")
    args=parser.parse_args()
    path=Path(__file__).with_name("independent_results.json")
    result=sample()
    if args.write_new:
        with path.open("x",encoding="utf-8") as stream:
            json.dump(result,stream,indent=2,allow_nan=False)
            stream.write("\n")
    else:
        compare(result,json.loads(path.read_bytes()))
    print(json.dumps(result))

