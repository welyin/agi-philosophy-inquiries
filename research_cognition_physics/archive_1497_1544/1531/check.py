"""Parity-projected torus products and a common-frame quark test. Python and NumPy."""
import argparse
import itertools
import json
import math
from decimal import Decimal, localcontext
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def parity_basis(n, parity):
    columns = []
    for j in range(n//2+1):
        partner = -j % n
        if j == partner and parity == -1:
            continue
        v = np.zeros(n)
        v[j] += 1
        v[partner] += parity
        columns.append(v/np.linalg.norm(v))
    return np.column_stack(columns)


def multiplication(n, parity, cutoff):
    basis = parity_basis(n, parity)
    out = parity_basis(2*n, 1)
    product = np.zeros((2*n,n,n))
    for k in range(2*n):
        for i in range(n):
            j = (k-i) % n
            momenta = i+n*np.arange(-cutoff,cutoff+1)
            product[k,i,j] = np.exp(-math.pi*(2*momenta-k)**2/(2*n)).sum()
    symmetric = []
    for i,j in itertools.combinations_with_replacement(range(basis.shape[1]),2):
        value = np.einsum("oij,i,j->o",product,basis[:,i],basis[:,j])
        symmetric.append((1 if i == j else math.sqrt(2))*value)
    matrix = out.T@np.column_stack(symmetric)
    return n**.25 * matrix


def quadratic(vector):
    result = np.zeros((3,3))
    for element,(i,j) in zip(vector,itertools.combinations_with_replacement(range(3),2)):
        result[i,j] = result[j,i] = element/(1 if i == j else math.sqrt(2))
    return result


def high_precision_kernel():
    """Bracket the three absolute eigenvalues without using singular-vector routines."""
    with localcontext() as ctx:
        ctx.prec=52
        z=Decimal(0)
        root2=Decimal(2).sqrt()
        pi=Decimal("3.14159265358979323846264338327950288419716939937510")
        def basis(n,j):
            jj=-j % n
            return [(j,Decimal(1))] if j==jj else [(j,1/root2),(jj,1/root2)]
        matrix=[]
        for k in range(5):          # the first five output rows are independent
            row=[]
            for i,j in itertools.combinations_with_replacement(range(3),2):
                value=z
                for a,wa in basis(5,i):
                    for b,wb in basis(5,j):
                        for c,wc in basis(10,k):
                            if (c-a-b)%5:
                                continue
                            exponent=[-pi*(2*(a+5*t)-c)**2/10 for t in range(-9,10)]
                            value+=wa*wb*wc*sum(v.exp() for v in exponent)
                row.append(value*(1 if i==j else root2))
            matrix.append(row)
        def determinant(rows):
            n=len(rows)
            result=z
            for perm in itertools.permutations(range(n)):
                odd=sum(perm[i]>perm[j] for i,j in itertools.combinations(range(n),2))
                result+=(-1)**odd*math.prod(rows[i][perm[i]] for i in range(n))
            return result
        co=[(-1)**j*determinant([row[:j]+row[j+1:] for row in matrix]) for j in range(6)]
        norm=sum(x*x for x in co).sqrt()
        assert norm>Decimal(".4")
        a,b,c,d,e,f=[x/norm for x in co]
        b,c,e=b/root2,c/root2,e/root2
        centers=[Decimal("-.5590169944"),Decimal("-.4478219619"),Decimal(".6978219619")]
        for center in centers:
            vals=[]
            for t in [center-Decimal("1e-9"),center+Decimal("1e-9")]:
                vals.append(determinant([[a-t,b,c],[b,d-t,e],[c,e,f-t]]))
            assert min(abs(x) for x in vals)>Decimal("1e-12") and vals[0]*vals[1]<0


def data():
    cfg = read(HERE/"inputs.json")
    parent = (HERE/cfg["same_curve"]).resolve()
    cfg = read(parent)["one_TeV_DR"]
    u = np.array(cfg["up"])
    d = np.array(read(parent.parent/cfg["down_from"])["down_yukawa_dr_1TeV"])
    a,b,c,phase = cfg["angles_order_12_13_23_delta"]
    def plane(i,j,theta,delta):
        rot=np.eye(3,dtype=complex)
        rot[i,i] = rot[j,j] = math.cos(theta)
        rot[i,j] = np.exp(-1j*delta)*math.sin(theta)
        rot[j,i] = -rot[i,j].conjugate()
        return rot
    return u,d,plane(1,2,c,0)@plane(0,2,b,phase)@plane(0,1,a,0)


def common_frame(parameters):
    g = np.zeros((3,3),dtype=complex)
    for p,(i,j) in enumerate([(0,1),(0,2),(1,2)]):
        g[i,j]=complex(*parameters[2*p:2*p+2])
        g[j,i]=-g[i,j].conjugate()
    return np.linalg.solve(np.eye(3)-g,np.eye(3)+g)


def weights(q,frame,ys):
    return ys * np.diag(frame.T @ np.diag(q) @ frame)


def find_frame(q,u,d,v):
    # Deterministic search merely locates a witness. Its *existence* is tested
    # by strict triangle margins with the archived coordinates during --check.
    def residual(x):
        w=common_frame(x)
        terms=[]
        for frame,y in [(w,u),(w@v,d)]:
            z=weights(q,frame,y)
            val=(z[2]-1j*z[1])/y[2]
            terms.extend([val.real,val.imag])
        return np.asarray(terms)
    rng=np.random.default_rng(1531)
    for _ in range(8):
        x=rng.uniform(-.8,.8,6)
        for _ in range(100):
            error=residual(x)
            if np.linalg.norm(error)<2e-13:
                return x.tolist()
            h=1e-6
            jac=np.column_stack([(residual(x+h*e)-residual(x-h*e))/(2*h) for e in np.eye(6)])
            step=np.linalg.lstsq(jac,-error,rcond=1e-12)[0]
            for n in range(15):
                trial=x+step/(2**n)
                if np.linalg.norm(residual(trial))<np.linalg.norm(error):
                    x=trial
                    break
            else:
                break
    raise ArithmeticError("No bracketed triangle witness located")


def certificate():
    high_precision_kernel()
    results={}
    for n,parity in [(4,1),(5,1),(7,-1),(8,-1)]:
        m=multiplication(n,parity,9)
        np.testing.assert_allclose(m,multiplication(n,parity,13),atol=1e-14,rtol=0)
        _,s,vh=np.linalg.svd(m,full_matrices=True)
        assert np.linalg.matrix_rank(m,tol=1e-11)==5
        q=np.linalg.svd(quadratic(vh[-1]),compute_uv=False)
        results[str(n)]=dict(dimension_out=m.shape[0],rank=5,kernel=1,cokernel=m.shape[0]-5,
                             quadratic_takagi=q.tolist(),min_nonzero_singular=s[4])
    q=np.asarray(results["5"]["quadratic_takagi"])
    u,d,v=data()
    parameters=read(HERE/"inputs.json")["common_frame_witness"]
    w=common_frame(parameters)
    np.testing.assert_allclose(w.conj().T@w,np.eye(3),atol=1e-14)
    cases={}
    for name,frame,ys in [("up",w,u),("down",w@v,d)]:
        z=weights(q,frame,ys[:,0])
        lens=abs(z)
        margin=(sum(lens)-2*max(lens))/ys[2,0]
        assert margin>1e-7
        cos_angle=(lens[2]**2-lens[0]**2-lens[1]**2)/(2*lens[0]*lens[1])
        phi=math.acos(float(cos_angle))
        tri=np.array([lens[0],lens[1]*np.exp(1j*phi),-lens[0]-lens[1]*np.exp(1j*phi)])
        phases=tri/z
        np.testing.assert_allclose(abs(phases),np.ones(3),atol=2e-12)
        y=frame@np.diag(ys[:,0]*phases)@frame.T
        assert abs(np.trace(np.diag(q)@y))<1e-12*ys[2,0]
        np.testing.assert_allclose(y@y.conj().T,frame@np.diag(ys[:,0]**2)@frame.conj().T,atol=2e-14)
        hi,lo=ys[:,0]+2*ys[:,1],ys[:,0]-2*ys[:,1]
        err=q[-1]-q[0]*(hi[0]+hi[1])/lo[2]
        near_axis=math.asin(math.sqrt(err/(q[-1]+q[0])))
        cases[name]=dict(weight_over_heavy=(lens/ys[2,0]).tolist(),triangle_margin_over_heavy=margin,
                         aligned_takagi_axis_min_error_over_heavy=err,
                         minimum_deg_to_any_identified_quadratic_axis=math.degrees(near_axis))
    return {"parity_cases":results, "N5_common_frame":cases,
            "ckm_abs_J":float(abs(np.imag(v[0,0]*v[1,1]*v[0,1].conj()*v[1,0].conj())))}


def compare(a,b):
    if isinstance(a,dict):
        return a.keys()==b.keys() and all(compare(x,b[k]) for k,x in a.items())
    return np.allclose(a,b,atol=2e-13,rtol=5e-9)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    flags=parser.add_mutually_exclusive_group(required=True)
    for name in ("--check","--write","--find-frame"):
        flags.add_argument(name,action="store_true")
    args=parser.parse_args()
    if args.find_frame:
        _,_,vh=np.linalg.svd(multiplication(5,1,9))
        q=np.linalg.svd(quadratic(vh[-1]),compute_uv=False)
        u,d,v=data()
        print(json.dumps(find_frame(q,u[:,0],d[:,0],v)))
    else:
        result=certificate()
        out=HERE/"results.json"
        if args.write:
            with out.open("x",encoding="utf-8") as file:
                file.write(json.dumps(result,indent=2,allow_nan=False)+"\n")
        else:
            assert compare(result,read(out))
        print("1531 parity ranks, physical null tensor and common-frame checks passed")
