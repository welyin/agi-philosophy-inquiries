"""Independent physical-kernel reconstruction for inversion-even orbifold modes."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PAIRS = list(itertools.combinations_with_replacement(range(3),2))


def read(rel):
    return json.loads((HERE/rel).read_text(encoding="utf-8-sig"))


def kernel(degree, odd):
    # Direct values of holomorphic sections. This does not use the theta product matrix.
    basis = np.eye(degree)
    folded = []
    for i in range(degree//2+1):
        v = basis[i]+(-1 if odd else 1)*basis[-i % degree]
        if np.linalg.norm(v)>0:
            folded.append(v/np.linalg.norm(v))
    assert len(folded)==3
    n = np.arange(-36,37)
    def products(z):
        fields = np.array([np.exp(-np.pi*n[n%degree==k,None]**2/degree+
                            2j*np.pi*n[n%degree==k,None]*z).sum(axis=0) for k in range(degree)])
        vec = np.array(folded)@fields
        return np.array([vec[i]*vec[j]*(1 if i==j else math.sqrt(2)) for i,j in PAIRS]).T
    points = np.linspace(.001,.483,30)+.073j
    sample = products(points)
    _,s,vh = np.linalg.svd(sample,full_matrices=False)
    assert s[-2]>.001 and s[-1]<1e-12
    relation = vh[-1].conj()
    np.testing.assert_allclose(products(np.linspace(.02,.71,23)-.103j)@relation,0,atol=1e-12)
    relation *= np.exp(-1j*np.angle(relation[np.argmax(abs(relation))]))
    assert np.max(abs(relation.imag))<1e-12
    q = np.zeros((3,3))
    for c,(i,j) in zip(relation.real,PAIRS):
        q[i,j]=q[j,i]=c/(1 if i==j else math.sqrt(2))
    return q


def physical_data():
    raw=read("../1527/inputs.json")["one_TeV_DR"]
    up=np.diag(np.asarray(raw["up"])[:,0])
    ds=np.asarray(read("../1527/"+raw["down_from"])["down_yukawa_dr_1TeV"])[:,0]
    a12,a13,a23,phase=raw["angles_order_12_13_23_delta"]
    def rotation(i,j,a,delta=0):
        m=np.eye(3,dtype=complex)
        m[i,i]=m[j,j]=math.cos(a)
        m[i,j]=math.sin(a)*np.exp(-1j*delta)
        m[j,i]=-m[i,j].conjugate()
        return m
    v=rotation(1,2,a23)@rotation(0,2,a13,phase)@rotation(0,1,a12)
    yd=v@np.diag(ds)@v.T
    np.testing.assert_allclose(np.linalg.svd(yd,compute_uv=False)[::-1],ds)
    return up,yd,v,ds


def takagi(a):
    # Both spectra of the actual witness are simple.
    val,u=np.linalg.eigh(a@a.conj().T)
    u=u[:,::-1]
    s=np.sqrt(val[::-1])
    diag=np.diag(u.conj().T@a@u.conj())
    u=u*np.exp(.5j*np.angle(diag))
    np.testing.assert_allclose((u*s)@u.T,a,rtol=1e-9,atol=1e-12)
    return u,s


def compute():
    qvalues={}
    for n,odd in [(4,False),(5,False),(7,True),(8,True)]:
        q=kernel(n,odd)
        qvalues[str(n)]=np.linalg.svd(q,compute_uv=False).tolist()
    q=kernel(5,False)
    u,y,v,ds=physical_data()
    k=np.zeros((3,3),complex)
    k[0,1]=k[1,0]=k[0,2]=k[2,0]=1
    k[1,2]=k[2,1]=-((y[0,1]+y[0,2])/y[1,2]).conjugate()
    for mat in (u,y):
        assert abs(np.vdot(k,mat))<1e-14
    assert abs(np.linalg.det(k))>1e-3
    ug,sg=takagi(q)
    up,sp=takagi(k)
    coords=read("inputs.json")["common_frame_witness"]
    gen=np.zeros((3,3),complex)
    for (i,j),xy in zip(itertools.combinations(range(3),2),zip(coords[::2],coords[1::2])):
        gen[i,j]=complex(*xy)
        gen[j,i]=-gen[i,j].conjugate()
    identity=np.eye(3)
    frame=2*np.linalg.inv(identity-gen)-identity
    margins=[]
    for basis,masses in [(frame,u.diagonal()),(frame@v,ds)]:
        sizes=np.abs(np.einsum("i,ij,ij,j->j",sg,basis,basis,masses))
        margins.append(float((sizes.sum()-2*sizes.max())/masses[-1]))
    assert min(margins)>1e-7
    ratios=sp/sg
    r=(up*np.sqrt(ratios))@ug.conj().T
    np.testing.assert_allclose(r@q@r.T,k,atol=1e-12)
    c=np.linalg.inv(r).conjugate()
    invc=np.linalg.inv(c)
    metric=invc.conj().T@invc
    np.testing.assert_allclose(c.conj().T@metric@c,np.eye(3),atol=2e-12)
    for yp in (u,y):
        yg=invc.T@yp@invc
        assert abs(np.vdot(q,yg))<2e-12
        np.testing.assert_allclose(c.T@yg@c,yp,atol=1e-13)
    rho=math.hypot(abs(y[0,1]),abs(y[0,2]))/abs(y[1,2])
    # All complex hollow missing tensors that annihilate the same measured down witness.
    family_bound=max((sg[1]+sg[2])/sg[0],sg[2]/(sg[0]*rho))
    achieved=ratios.max()/ratios.min()
    np.testing.assert_allclose(np.linalg.cond(metric),achieved,rtol=1e-10)
    assert 1<family_bound<=achieved
    return dict(direct_theta_takagi_spectra=qvalues, down_offdiag_abs=[float(abs(y[i,j]))
                for i,j in ((0,1),(0,2),(1,2))], physical_K23=[k[1,2].real,k[1,2].imag],
                physical_hollow_takagi_spectrum=sp.tolist(), geometric_spectrum=sg.tolist(),
                hollow_witness_optimal_kinetic_condition=achieved,
                all_hollow_blind_witnesses_kinetic_condition_lower=family_bound,
                ratio_control=rho, hollow_trace_identity=sp[0]-sp[1]-sp[2],
                fixed_bulk_kinetic_author_triangle_margins=margins,
                canonical_map_condition=float(np.linalg.cond(c)),
                joint_zero_trace_error=max(float(abs(np.vdot(k,x))) for x in (u,y)))


def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for key in a:
            compare(a[key],b[key])
    else:
        np.testing.assert_allclose(a,b,rtol=2e-9,atol=4e-12)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write-new",action="store_true")
    mode.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=compute()
    if args.write_new:
        with (HERE/"independent_results.json").open("x",encoding="utf8") as f:
            json.dump(result,f,indent=2)
            f.write("\n")
    else:
        compare(result,read("independent_results.json"))
        print("1531 independent parity modes, dual tensor, metric and CKM witness verified.")
