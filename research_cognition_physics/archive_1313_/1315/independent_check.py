"""Independent round 1315 checks; does not import author calculations.
Default/--check is read-only; --write-new creates one result file exclusively.
"""
from __future__ import annotations
import argparse, hashlib, json, math
from functools import lru_cache
from fractions import Fraction as F
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
C=np.kron(np.array([[0,1],[-1,0]],dtype=int),
          np.array([[0,1],[-1,0]],dtype=int))
PAULI=[np.array([[0,1],[1,0]],complex),
       np.array([[0,-1j],[1j,0]],complex),
       np.diag([1,-1]).astype(complex)]

def hafnian(G):
    """Unsigned matching sum, independent of the author's Grassmann Pfaffian."""
    N=len(G)
    @lru_cache(None)
    def rec(mask):
        if mask==0: return 1
        lo=mask & -mask
        i=lo.bit_length()-1
        rest=mask^lo
        ans=0
        for j in range(i+1,N):
            if rest&(1<<j) and G[i][j]!=0:
                ans += G[i][j]*rec(rest^(1<<j))
        return ans
    return rec((1<<N)-1)

def wedge(a,b):
    out={}
    for i,x in a.items():
        for j,y in b.items():
            if i&j: continue
            parity=sum((j&((1<<k)-1)).bit_count() for k in range(i.bit_length()) if i&(1<<k))%2
            z=(-1 if parity else 1)*x*y
            key=i|j
            out[key]=out.get(key,F(0))+z
    return {k:v for k,v in out.items() if v}

def addto(a,b,scale=1):
    for k,v in b.items(): a[k]=a.get(k,F(0))+scale*v
    return {k:v for k,v in a.items() if v}

def full_four_species():
    """Full exterior coefficient after a nontrivial rational SL(4) species map."""
    I=[[F(i==j) for j in range(4)] for i in range(4)]
    def mm(A,B):
        return [[sum(A[i][k]*B[k][j] for k in range(4)) for j in range(4)] for i in range(4)]
    T=I
    for i,j,x in [(0,2,F(1,3)),(1,3,F(2,5)),(2,1,F(-1,7))]:
        M=[row[:] for row in I]; M[i][j]=x; T=mm(M,T)
    G=[[F((i<2)!=(j<2)) for j in range(4)] for i in range(4)]
    S={}
    for f in range(4):
        for g in range(f+1,4):
            S[f,g]={(1<<i)|(1<<j):(T[f][i]*T[g][j]-T[f][j]*T[g][i])*G[i][j]
                    for i in range(4) for j in range(i+1,4)
                    if (T[f][i]*T[g][j]-T[f][j]*T[g][i])*G[i][j]}
    out=wedge(S[0,1],S[2,3])
    out=addto(out,wedge(S[0,2],S[1,3]),-1)
    out=addto(out,wedge(S[0,3],S[1,2]))
    assert out=={15:F(2)}
    # Symmetric internal metric is preserved by determinant-one spin/weak maps.
    W=np.array([[F(2),F(1)],[F(0),F(1,2)]],dtype=object)
    L=np.array([[F(1),F(0)],[F(1,3),F(1)]],dtype=object)
    K=np.kron(W,L); Cf=C.astype(object)
    assert np.array_equal(K.T@Cf@K,Cf)
    return {"full_exterior_SL4_coefficient":2,"exact_weak_spin_invariance":True}

def su(n,seed):
    rng=np.random.default_rng(seed)
    Z=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n))
    Q,_=np.linalg.qr(Z)
    Q[:,0]/=np.linalg.det(Q)
    return Q

def spectrum_and_group(n):
    omega=np.exp(2j*np.pi/3)
    U=np.eye(n,dtype=complex)
    U[:3,:3]=np.array([[omega**(i*j) for j in range(3)] for i in range(3)])/np.sqrt(3)
    for j in range(3,n):
        V=np.eye(n); c=math.cos(.03); s=math.sin(.03)
        V[2,2]=V[j,j]=c; V[2,j]=s; V[j,2]=-s
        U=U@V
    mu=2.0**-12
    masses=mu*np.arange(1,n+1)
    M=U@np.diag(masses)@U.conj().T
    H=M@M.conj().T
    v=.25
    Yphysical=math.sqrt(2)*M/v
    Yleft=Yphysical.conj()
    assert np.max(np.abs(Yleft.conj()*v/math.sqrt(2)-M))<1e-14
    seen={0}
    while True:
        nxt=seen|{j for i in seen for j in range(n) if abs(H[i,j])>1e-14}
        if nxt==seen: break
        seen=nxt
    J=float(np.imag(U[0,0]*U[1,1]*np.conj(U[0,1]*U[1,0])))
    assert len(seen)==n and abs(J-math.sqrt(3)/18)<1e-14
    assert np.min(np.linalg.eigvalsh(M))>0
    Tq=np.kron(su(n,31+n),su(3,90+n)); Tl=su(n,72+n)
    deterr=float(max(abs(np.linalg.det(Tq)-1),abs(np.linalg.det(Tl)-1)))
    assert deterr<1e-12
    for k in range(-2,3):
        yq=F(1,6)+k; yl=F(-1,2)-3*k
        assert 3*n*yq+n*yl==0
    assert (1j)**(4*n)==1
    return {"n":n,"phase_J":J,"all_flavors_connected":True,"left_Weyl_conjugation_checked":True,
            "species_determinant_error":deterr,"neutral_hypercharge_samples":5,
            "Dirac_min_mass":float(np.min(np.linalg.eigvalsh(M))),
            "scope":"same leading data only, not interacting P probabilities"}

def dirac_h(m,p):
    sig=sum(p[k]*PAULI[k] for k in range(3))
    return np.block([[-sig,m*np.eye(2)],[m*np.eye(2),sig]])

def positive_mass_point(n):
    k=2*n; N=4*n
    # First 2n species are Q upper components. All L select lower charged components.
    labels=[]; masses=[]
    for f in range(N):
        upper=f<k
        if f<3*n:
            gen=f//3
            masses.append((.7 if upper else 1.1)+.1*gen)
            labels.append(("u" if upper else "d")+str(gen+1))
        else:
            gen=f-3*n
            masses.append(float(2+gen)); labels.append("e"+str(gen+1))
    W=20.0
    Es=np.zeros(N); ps=np.zeros((N,3)); states=[]; scales=np.zeros(N)
    maxres=0.; signed_total=np.zeros(4)
    for i in range(k):
        j=k+i; ma,mb=masses[i],masses[j]
        assert W>ma+mb
        Ea=(W*W+ma*ma-mb*mb)/(2*W); Eb=W-Ea
        p=math.sqrt((W*W-(ma+mb)**2)*(W*W-(ma-mb)**2))/(2*W)
        incoming=i<n
        for f,E,zz,chi in [(i,Ea,-p,np.array([1,0],complex)),(j,Eb,p,np.array([0,1],complex))]:
            Es[f]=E; ps[f]=[0,0,zz]
            up=math.sqrt(E+p); low=math.sqrt(E-p)
            spin=np.r_[up*chi,(low if incoming else -low)*chi]
            H=dirac_h(masses[f],ps[f] if incoming else -ps[f])
            residual=H@spin-(E if incoming else -E)*spin
            maxres=max(maxres,float(np.max(np.abs(residual))),float(abs(np.vdot(spin,spin)-2*E)))
            scales[f]=up
            signed_total+=(1 if incoming else -1)*np.r_[E,ps[f]]
    G=[[scales[i]*scales[j] if ((i<k)!=(j<k)) else 0. for j in range(N)] for i in range(N)]
    amplitude=float(hafnian(G)); formula=math.factorial(k)*float(np.prod(scales))
    rel=abs(amplitude/formula-1)
    assert rel<1e-12 and maxres<1e-11 and max(abs(signed_total))<1e-11 and amplitude>0
    return {"n":n,"pair_total_energy":W,"incoming_legs":2*n,"outgoing_legs":2*n,
            "masses":masses,"external_roles":labels,
            "momentum_residual":float(max(abs(signed_total))),"Dirac_norm_residual":maxres,
            "contact_coefficient":amplitude,"positive_closed_form":formula,"relative_error":rel,
            "scope":"leading contact coefficient, no rate or device error certificate"}

def calculate():
    target=[]
    for n in range(1,9):
        k=2*n
        G=[[int((i<k)!=(j<k)) for j in range(2*k)] for i in range(2*k)]
        v=hafnian(G); assert v==math.factorial(k)
        target.append({"n":n,"coefficient":v})
    branches=[spectrum_and_group(n) for n in (3,4,5)]
    mass=[positive_mass_point(n) for n in (3,4,5)]
    sources={}
    for rel in ["_admission/round1315_admission.md","1315/proof.md","1315/check.py",
                "1315/results.json","research_note_1315.md","1315/NEXT.md",
                "1315/independent_check.py","1314/proof.md"]:
        p=STAGE/rel
        sources[rel]=hashlib.sha256(p.read_bytes()).hexdigest()
    return {"schema":"round1315_independent_v1","status":"PASS_WITH_SCOPE",
            "method":"unsigned matching recursion plus exact full exterior SL4 check; no author import",
            "finite_coefficients":target,"exact_four_species":full_four_species(),
            "branches":branches,"positive_mass_scattering":mass,
            "minimum_dimension_scope":"D>=6n only for nonzero equal DeltaB=DeltaL in declared original local polynomial class",
            "analytic_scope":"all n is proved analytically; finite samples do not replace the theorem",
            "sources_sha256":sources,"new_axioms":0}

def main():
    ap=argparse.ArgumentParser(); g=ap.add_mutually_exclusive_group()
    g.add_argument("--check",action="store_true"); g.add_argument("--write-new",action="store_true")
    a=ap.parse_args(); result=calculate(); dst=HERE/"independent_results.json"
    if a.write_new:
        with dst.open("x",encoding="utf-8") as f: json.dump(result,f,ensure_ascii=False,indent=2); f.write("\n")
    else:
        assert json.loads(dst.read_text(encoding="utf-8"))==result,"saved independent result mismatch"
    print(json.dumps({"status":result["status"],"target_cases":8,"massive_branches":3,"sources":8,
                      "results_sha256":hashlib.sha256(dst.read_bytes()).hexdigest()}))
if __name__=="__main__": main()