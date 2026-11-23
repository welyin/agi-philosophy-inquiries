"""Round 1345 independent charge pairing checks, with all partner words retained."""
import argparse
import itertools
import json
from pathlib import Path
import numpy as np


def charges(kind):
    q=[]
    if kind == "SU3":
        for i,j in itertools.combinations(range(3),2):
            x=np.zeros((3,3),complex); x[i,j]=.5
            q.extend([x+x.T,1j*(x.T-x)])
        q.extend([np.diag([1.,-1,0])/2,np.diag([1.,1,-2])/np.sqrt(12)])
    else:
        for i,j in [(1,2),(2,0),(0,1)]:
            x=np.zeros((3,3),complex); x[i,j]=1j
            q.append(x+x.conj().T)
    C=np.trace(sum(x@x for x in q)).real/3
    return np.array(q)/np.sqrt(C)


def rms(Z):
    return float(np.sqrt(np.sum(abs(Z)**2)/Z.shape[-1]))


def pair(a):
    return np.array([np.kron(x,np.eye(3))+np.kron(np.eye(3),x) for x in a])


def connected(q):
    partial=q.reshape(-1,3,3,3,3).trace(axis1=2,axis2=4)/3
    mean=partial.trace(axis1=1,axis2=2)/3
    return q-pair(partial)+mean[:,None,None]*np.eye(9)


def embed(q,i,j):
    order=[i,j]+[k for k in range(4) if k not in (i,j)]
    perm=np.arange(81).reshape([3]*4).transpose(order).ravel()
    out=np.empty((len(q),81,81),complex)
    for k,m in enumerate(q):
        out[k][np.ix_(perm,perm)]=np.kron(m,np.eye(9))
    return out


def readout(blocks,t):
    # Taylor propagation: cosine and sinc of positive V^dagger V.
    V=blocks.reshape(-1,81)
    K=V.conj().T@V
    unit=np.eye(81)
    term=unit.copy()
    ground=unit.astype(complex)
    drive=unit.astype(complex)
    for h in range(1,17):
        term=term@K*(-t**2)/((2*h)*(2*h-1))
        ground+=term
        drive+=term/(2*h+1)
    excited=-1j*t*V@drive
    assert np.allclose(ground.conj().T@ground+excited.conj().T@excited,unit,atol=1e-12)
    return ground,excited


def sym4_projector():
    words=np.array(list(itertools.product(range(3),repeat=4)))
    labels=np.array([np.bincount(w,minlength=3) for w in words])
    groups=np.unique(labels,axis=0)
    W=np.column_stack([np.all(labels==g,axis=1) for g in groups]).astype(float)
    W/=np.linalg.norm(W,axis=0)
    return W@W.T


def run():
    swap=np.eye(9).reshape(3,3,3,3).transpose(1,0,2,3).reshape(9,9)
    anti=(np.eye(9)-swap)/2
    volume=np.zeros((9,3))
    for k,(a,b) in enumerate([(1,2),(2,0),(0,1)]):
        volume[3*a+b,k],volume[3*b+a,k]=1/np.sqrt(2),-1/np.sqrt(2)
    assert np.allclose(volume@volume.T,anti)
    P=sym4_projector()
    assert np.isclose(np.trace(P),15)
    rho=P.reshape([3]*8)/15
    for _ in range(3):
        rho=rho.trace(axis1=1,axis2=rho.ndim//2+1)
    assert np.allclose(rho,np.eye(3)/3)
    rng=np.random.default_rng(1345)
    general_wedge=0.
    for _ in range(12):
        b=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
        residual=volume.T@pair([b])[0]@volume-(np.trace(b)*np.eye(3)-b.T)
        general_wedge=max(general_wedge,np.max(abs(residual)))
    assert general_wedge<1e-12
    checks={}
    for group in ["SU3","SO3"]:
        T=charges(group)
        total=pair(T)
        q=total@anti
        z=rms(connected(q))
        assert abs(z**2-5/18)<1e-12
        assert np.allclose(volume.T@q@volume,-T.transpose(0,2,1))
        V0=embed(q,0,1)+embed(q,2,3)
        V1=embed(q,0,2)+embed(q,1,3)
        assert abs(rms(V0-V1)-2*z)<1e-12
        assert rms(V0@P)<1e-12 and rms(V1@P)<1e-12
        u0,v0=readout(V0,.25)
        u1,v1=readout(V1,.25)
        records={"minus_exc":rms((v0-v1)/2)**2,"plus_exc":rms((v0+v1)/2)**2,
                 "minus_ground":rms((u0-u1)/2)**2,"plus_ground":rms((u0+u1)/2)**2}
        assert np.isclose(sum(records.values()),1.)
        y=rms(q-total)
        assert abs(y-np.sqrt(6)*z)<1e-12
        checks[group]=dict(connected=z,source_error=y,sharp_ratio=y/z,
                           symmetric_correlated_dark_max=float(max(rms(V0@P),rms(V1@P))),
                           records=records)
    largest=0.
    eye=np.eye(9)
    ps=(eye+swap)/2
    for _ in range(12):
        raw=rng.normal(size=(9,9))+1j*rng.normal(size=(9,9))
        yf=ps@(raw+raw.conj().T)@ps
        yf-=np.trace(yf)*ps/6
        Y=yf[None,:,:]
        z=rms(connected(Y))
        assert rms(Y)<=np.sqrt(6)*z+1e-12
        corrected=Y+pair(np.array([np.diag([.1,-.3,.2])]))
        eta=rms(volume.T@corrected@volume)
        assert rms(corrected)<=np.sqrt(6)*rms(connected(corrected))+np.sqrt(2)*eta+1e-12
        largest=max(largest, rms(Y)/z)
    return dict(groups=checks, general_wedge_identity_error_max=float(general_wedge),
                arbitrary_symmetric_sector_checks=12, max_random_ratio=largest,
                optimal_calibrated_RMS_coefficient=float(np.sqrt(6)))


def compare(actual,saved):
    if isinstance(actual,dict):
        assert actual.keys()==saved.keys()
        for k in actual:
            compare(actual[k],saved[k])
    else:
        assert np.isclose(actual,saved,atol=2e-12,rtol=2e-12),(actual,saved)


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    cli.add_argument("--check",action="store_true")
    args=cli.parse_args()
    data=run()
    output=Path(__file__).with_name("independent_results.json")
    if args.check:
        compare(data,json.loads(output.read_text(encoding="utf-8")))
    else:
        with output.open("x",encoding="utf-8") as f:
            json.dump(data,f,indent=2,allow_nan=False);f.write("\n")
    print("PASS", data["groups"]["SU3"]["records"]["minus_exc"],data["groups"]["SO3"]["records"]["minus_exc"])
