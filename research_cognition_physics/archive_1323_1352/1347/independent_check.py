"""Independent round 1347: amplitude-level, transverse polarization double-soft currents."""
import argparse
import itertools
import json
from pathlib import Path
import numpy as np


def generators(name):
    out=[]
    for a,b in [(0,1),(1,2),(2,0)]:
        E=np.zeros((3,3),complex)
        E[a,b]=1.
        if name=="SO3":
            out.append(1j*(E.T-E))
        else:
            out.extend([(E+E.T)/2,1j*(E.T-E)/2])
    if name=="SU3":
        out.extend([np.diag([1.,-1,0])/2,np.diag([1.,1,-2])/np.sqrt(12)])
    return np.array(out)


def m_dot(a,b):
    return a[0]*b[0]-a[1:]@b[1:]


def polarizations(q):
    axis=q[1:]/q[0]
    seed=np.eye(3)[np.argmin(abs(axis))]
    first=np.cross(axis,seed)
    first/=np.linalg.norm(first)
    return [np.r_[0.,first],np.r_[0.,np.cross(axis,first)]]


def currents(p,legs,q,k,e,f,omega):
    # Catani-Grazzini Eq. (102). legs here are generators of one two-leg singlet.
    j1=sum(m_dot(P,e)/m_dot(P,q)*T for P,T in zip(p,legs))
    j2=sum(m_dot(P,f)/m_dot(P,k)*T for P,T in zip(p,legs))
    u=m_dot(q,k)
    gamma=[]
    for P in p:
        a,b=m_dot(P,q),m_dot(P,k)
        first=(m_dot(P,e)*m_dot(q,f)-m_dot(P,f)*m_dot(k,e))/(u*(a+b))
        second=(a-b)/(2*(a+b))*(m_dot(P,e)*m_dot(P,f)/(a*b)+m_dot(e,f)/u)
        gamma.append(first-second)
    G=sum(t*L for t,L in zip(gamma,legs))
    # omega[a,b,c] T_c = [T_a,T_b]; it already includes i.
    anticom=(np.einsum("aij,bjk->abik",j1,j2)+np.einsum("bij,ajk->abik",j2,j1))/2
    double=anticom+np.einsum("abc,cij->abij",omega,G)
    return j1,j2,double


def piece_norms(name,hard,direction1,direction2,energies):
    T=generators(name); d=len(T)
    C=sum(np.trace(x@x).real for x in T)/3
    ind=np.trace(T[0]@T[0]).real
    comm=np.einsum("aij,bjk->abik",T,T)-np.einsum("bij,ajk->abik",T,T)
    omega=np.einsum("abij,cji->abc",comm,T)/ind
    legs=[np.array([np.kron(x,np.eye(3)) for x in T]),
          np.array([np.kron(np.eye(3),-x.T) for x in T])]
    singlet=np.eye(3).reshape(-1)/np.sqrt(3)
    P=[np.r_[1.,n] for n in hard]
    q,k=[t*np.r_[1.,n] for t,n in zip(energies,[direction1,direction2])]
    blocks=np.zeros((4,4))
    check_bose=0.
    for e,f in itertools.product(polarizations(q),polarizations(k)):
        A,a,A2=currents(P[:2],legs,q,k,e,f,omega)
        B,b,B2=currents(P[2:],legs,q,k,e,f,omega)
        _,_,flipped=currents(P[:2],legs,k,q,f,e,omega)
        check_bose=max(check_bose,float(np.max(abs(A2-flipped.transpose(1,0,2,3)))))
        _,_,ward_q=currents(P[:2],legs,q,k,q,f,omega)
        _,_,ward_k=currents(P[2:],legs,q,k,e,k,omega)
        assert np.max(abs(ward_q@singlet))<1e-9
        assert np.max(abs(ward_k@singlet))<1e-9
        a1,a2,b1,b2=A@singlet,a@singlet,B@singlet,b@singlet
        every=np.array([
            np.einsum("abi,j->abij",A2@singlet,singlet),
            np.einsum("i,abj->abij",singlet,B2@singlet),
            np.einsum("ai,bj->abij",a1,b2),
            np.einsum("bi,aj->abij",a2,b1)
        ]).reshape(4,-1)
        blocks+=(every.conj()@every.T).real
    def w(i,j,l):
        return m_dot(P[i],P[j])/(m_dot(P[i],l)*m_dot(P[j],l))
    def h(l):
        return w(0,2,l)-w(0,3,l)-w(1,2,l)+w(1,3,l)
    expected=4*C*C/d*h(q)*h(k)
    isolated=np.trace(blocks)
    independent=blocks.sum()-isolated
    assert np.allclose(independent,expected,rtol=2e-10,atol=2e-10)
    assert np.allclose(2*blocks[0,1],expected/2,rtol=2e-10,atol=2e-10)
    assert np.allclose(2*blocks[2,3],expected/2,rtol=2e-10,atol=2e-10)
    assert np.max(abs(blocks[:2,2:]))<1e-9
    assert check_bose<1e-9
    norm=4*C*C*w(0,1,q)*w(2,3,k)
    return dict(sources_sum=float(blocks.sum()), source_label_overlay=float(isolated),
                difference=float(independent), normalized=float(independent/norm),
                same_source_branch_cross=float(2*blocks[0,1]),
                interchanged_emission_cross=float(2*blocks[2,3]), bose_residual=check_bose)


def angular(angles,opening):
    cosnode, weights=np.polynomial.legendre.leggauss(32)
    low=np.cos(np.deg2rad(opening))
    mu=(1+low+(1-low)*cosnode[:,None])/2
    phi=2*np.pi*(np.arange(192)+.5)/192
    r=np.sqrt(1-mu**2)
    answer=[]
    for i,j in [(0,1),(2,3),(0,2),(0,3),(1,2),(1,3)]:
        delta=angles[i]-angles[j]
        denom=(1-r*np.cos(phi-angles[i]))*(1-r*np.cos(phi-angles[j]))
        answer.append(weights@np.mean((1-np.cos(delta))/denom,axis=1))
    contrast=answer[2]-answer[3]-answer[4]+answer[5]
    return contrast**2/(answer[0]*answer[1])


def run():
    angles=np.deg2rad([0.,180.,60.,240.])
    hard=np.column_stack([np.cos(angles),np.sin(angles),np.zeros(4)])
    rng=np.random.default_rng(134705)
    data={}
    for name,d in [("SU3",8),("SO3",3)]:
        normal=piece_norms(name,hard,np.array([0.,0.,1.]),np.array([0.,0.,-1.]),[.91,.63])
        assert abs(normal["normalized"]-1/d)<1e-12
        flips=[]
        for _ in range(3):
            dirs=rng.normal(size=(2,3))
            dirs/=np.linalg.norm(dirs,axis=1)[:,None]
            flips.append(piece_norms(name,hard,*dirs,[.5,.72])["normalized"])
        windows={str(beta):float(angular(angles,beta)/d) for beta in (20,30)}
        null_angles=np.deg2rad([0.,180.,90.,270.])
        assert max(abs(angular(null_angles, beta)) for beta in (20,30))<1e-22
        data[name]=dict(normal=normal,opposite_cones=windows,
                        unrelated_soft_directions=flips)
    return data


def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:
            compare(a[k],b[k])
    else:
        assert np.allclose(a,b,rtol=5e-10,atol=5e-10), (a,b)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    args=p.parse_args()
    data=run()
    out=Path(__file__).with_name("independent_results.json")
    if args.check:
        compare(data,json.loads(out.read_text(encoding="utf-8")))
    else:
        with out.open("x",encoding="utf-8") as f:
            json.dump(data,f,indent=2,allow_nan=False);f.write("\n")
    print("PASS",data["SU3"]["opposite_cones"],data["SO3"]["opposite_cones"])
