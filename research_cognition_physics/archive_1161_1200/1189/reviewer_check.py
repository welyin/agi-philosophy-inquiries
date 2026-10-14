"""Independent 1189 review: reconstruct 940 Dirac/CAR sources, no author imports.
The untruncated theorem is analytic; the total-number-4 evolution is diagnostic only.
Default: recompute and compare this review's saved results. --write: exclusive creation.
"""
import argparse
import json
import math
from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent
TOL = 3e-10

def dgamma(A, masks):
    out = np.zeros((len(masks), len(masks)), complex)
    where = {v:k for k,v in enumerate(masks)}
    for col, mask in enumerate(masks):
        for j in range(4):
            if not mask & (1 << j):
                continue
            after = mask ^ (1 << j)
            s1 = (-1)**((mask & ((1 << j)-1)).bit_count())
            for i in range(4):
                if after & (1 << i):
                    continue
                s2 = (-1)**((after & ((1 << i)-1)).bit_count())
                out[where[after | (1 << i)], col] += s1*s2*A[i,j]
    return out

def reconstruct():
    eye = np.eye(2)
    zero = np.zeros((2,2))
    pauli = [np.array([[0,1],[1,0]],complex),
             np.array([[0,-1j],[1j,0]],complex),
             np.diag([1.,-1.])]
    gamma = [np.block([[eye,zero],[zero,-eye]])]
    gamma += [np.block([[zero,s],[-s,zero]]) for s in pauli]
    p, pp = np.array([1.,0.,0.]), np.array([1.,0.,1.])
    def spinors(momentum):
        e = math.sqrt(1 + momentum@momentum)
        return math.sqrt((e+1)/(2*e))*np.vstack(
            (eye, sum(x*s for x,s in zip(momentum,pauli))/(e+1)))
    ep, epp = math.sqrt(2), math.sqrt(3)
    up, upp = spinors(p), spinors(pp)
    four_sum = np.r_[ep+epp,p+pp]
    tensor = np.empty((4,4,2,2),complex)
    for mu,nu in product(range(4),repeat=2):
        tensor[mu,nu] = (upp.conj().T@gamma[0]@
            (gamma[mu]*four_sum[nu]+gamma[nu]*four_sum[mu])@up/4)
    qcov = np.r_[epp-ep,-(pp-p)]
    ward = float(np.max(np.abs(np.einsum("m,mnij->nij",qcov,tensor))))
    masks = sorted(sum(1 << j for j in indices)
                   for indices in combinations(range(4),2))
    def real_pair(block):
        transfer = np.zeros((4,4),complex)
        transfer[2:, :2] = block
        return ((transfer+transfer.conj().T)/math.sqrt(2),
                (transfer-transfer.conj().T)/(1j*math.sqrt(2)))
    one_sources = []
    for block in ((tensor[1,1]-tensor[2,2])/math.sqrt(2),
                  (tensor[1,2]+tensor[2,1])/math.sqrt(2)):
        one_sources.extend(real_pair(block))
    sources = [dgamma(A,masks) for A in one_sources]
    hm = dgamma(np.diag([ep,ep,epp,epp]),masks)
    def normal(A,B):
        aa,bb = dgamma(A,masks), dgamma(B,masks)
        return (aa@bb+bb@aa-dgamma(A@B+B@A,masks))/2
    pairs = {(mu,nu):real_pair(tensor[mu,nu])
             for mu,nu in product(range(4),repeat=2)}
    v0 = np.zeros((6,6),complex)
    for phase in range(2):
        rho = pairs[0,0][phase]
        sigma = sum(pairs[i,i][phase] for i in range(1,4))
        flow = [pairs[0,i][phase] for i in range(1,4)]
        cg = normal(rho,rho)/2 + normal(rho,sigma)
        cg -= 2*(normal(flow[0],flow[0])+normal(flow[1],flow[1]))
        cg -= F(3,2).__float__()*normal(flow[2],flow[2])
        v0 -= cg/8
    return masks, hm, sources, v0, ward

def exact_bounds():
    eps = F(16,5)*F(3,5)**33 + F(3,10**12)
    delta_abs = 2*eps + eps*eps
    # Four TT modes at omega0=1; J<=1, V0<=3.
    d, r = F(5)*delta_abs, F(6)*delta_abs
    # sqrt(2)<=3/2, b<=2/5, u<=3/5, w<=4.
    eta = d*(F(2)+F(3,10)+F(1,75))
    z0 = F(17,10)
    assert F(0) < r/2 < 1
    zp = z0/(1-r/2)  # exp(x)<=1/(1-x), 0<=x<1
    source_bound = eta*(F(6,5)+4*(z0+zp))
    # Direct one-particle gamma bounds using Ep+Ep'<=16/5.
    rho, sigma = F(8,5), F(3,2)
    vx,vy,vz = F(13,10),F(4,5),F(21,20)
    v_bound = (rho*rho/2+rho*sigma+2*(vx*vx+vy*vy)
               +F(3,2)*vz*vz)/4
    assert v_bound == F(1599,640) < 3
    assert eta < F(4,10**6)
    assert source_bound < F(6,10**5)
    # Preserve 1180 n=32 coefficients at s=1, including positive base.
    S = sum(F(math.comb(2*j,j),4**j)*F(3,5)**j for j in range(33))
    actual_delta = F(2,5)*S*S + F(1,10**12) - 1
    assert abs(actual_delta) <= delta_abs
    result = {}
    for name,value in (("epsilon",eps),("delta_abs",delta_abs),("d",d),
                       ("r",r),("eta",eta),("source_error",source_bound),
                       ("V0_bound",v_bound),("actual_delta",actual_delta)):
        result[name] = {"exact":str(value),"display":float(value)}
    result["initial_code_dimension"] = 6*5
    result["per_mode_cap_one_dimension_not_covered_by_current_m"] = 6*16
    return result

def bosons(cap=4, modes=4):
    basis = sorted((x for x in product(range(cap+1),repeat=modes)
                    if sum(x)<=cap), key=lambda x:(sum(x),x))
    where = {x:i for i,x in enumerate(basis)}
    Q,Q2 = [],[]
    for k in range(modes):
        q = np.zeros((len(basis),len(basis)))
        q2 = np.zeros_like(q)
        for col, occ in enumerate(basis):
            n = occ[k]
            q2[col,col] = n+.5
            for jump, weight in ((1,math.sqrt(n+1)/math.sqrt(2)),
                                 (-1,math.sqrt(n)/math.sqrt(2))):
                new = list(occ);new[k] += jump;new=tuple(new)
                if new in where:q[where[new],col] = weight
            for jump, weight in ((2,math.sqrt((n+1)*(n+2))/2),
                                 (-2,math.sqrt(max(n*(n-1),0))/2)):
                new = list(occ);new[k] += jump;new=tuple(new)
                if new in where:q2[where[new],col] = weight
        Q.append(q);Q2.append(q2)
    number = np.array([sum(v) for v in basis],float)
    return basis,Q,Q2,number

def diagnostic(hm,sources,v0,masks,delta):
    basis,Q,Q2,number = bosons()
    nb = len(basis);eye_b=np.eye(nb);eye_m=np.eye(6);kap=.1
    driven = sum(np.kron(j,q) for j,q in zip(sources,Q))
    second = sum(Q2)
    def ham(k):
        return np.kron(hm+k*k*v0,eye_b)+np.kron(eye_m,np.diag(number+2))-k*driven/2
    h0=ham(kap);hp=h0+float(delta)*np.kron(eye_m,second)/2
    src=np.kron(2*kap*v0,eye_b)-driven/2
    psi=np.zeros(6*nb,complex)
    for mask in (3,5):psi[masks.index(mask)*nb] = 1/math.sqrt(2)
    def evolved(h):
        vals,vecs=np.linalg.eigh(h)
        return vecs@(np.exp(-1j*vals)*(vecs.conj().T@psi))
    v,z=evolved(h0),evolved(hp)
    mean_num=float(np.vdot(v,np.tile(number,6)*v).real)
    out_prob=float(np.sum(abs(v[np.tile(number>1,6)])**2))
    gap=float(np.linalg.norm(z-v))
    srcgap=float(abs(np.vdot(z,src@z)-np.vdot(v,src@v)))
    fd=(ham(kap+1e-5)-ham(kap-1e-5))/(2e-5)
    deriv=float(np.max(abs(fd-src)))
    energy=float(abs(np.vdot(v,h0@v)-np.vdot(psi,h0@psi)))
    inside=np.flatnonzero(np.tile(number<=1,6))
    outside=np.flatnonzero(np.tile(number>1,6))
    leakage=float(np.linalg.norm(h0[np.ix_(outside,inside)],2))
    assert leakage>.04
    assert mean_num>1e-5 and out_prob>1e-8
    assert gap<4e-6 and srcgap<6e-5
    assert deriv<1e-9 and energy<1e-11
    # Q2 is the compression of full Q^2, not the square of compressed Q.
    boundary=float(np.linalg.norm(Q2[0]-Q[0]@Q[0],2))
    assert boundary > 1
    return dict(dimension=len(psi),cap=4,driven_mean_number=mean_num,
                outside_initial_code=out_prob,one_input_vector_gap=gap,
                source_expectation_gap=srcgap,
                source_derivative_max_residual=deriv,energy_residual=energy,
                exact_initial_code_generator_leakage=leakage,
                true_Q2_vs_compressed_Q_squared_norm=boundary,
                full_Fock_certificate=False)

def run():
    masks,hm,sources,v0,ward = reconstruct()
    snorm=[float(np.linalg.norm(j,2)) for j in sources]
    comm_h=max(float(np.linalg.norm(hm@j-j@hm,2)) for j in sources)
    comm_j=max(float(np.linalg.norm(a@b-b@a,2)) for a,b in product(sources,repeat=2))
    herm=max(float(np.max(abs(a-a.conj().T))) for a in [hm,v0]+sources)
    assert ward<1e-12 and herm<1e-12
    assert max(snorm)<=1+1e-12 and np.linalg.norm(v0,2)<3
    assert comm_h>.2 and comm_j>.49
    bounds=exact_bounds()
    diag=diagnostic(hm,sources,v0,masks,F(bounds["actual_delta"]["exact"]))
    # Independent exact comparisons with author output; never import its code.
    authored=json.loads((ROOT/"results.json").read_text(encoding="utf-8-sig"))
    ab=authored["groups"]["2_exact_full_Fock_and_source_bounds"]
    keys={"epsilon":"epsilon","d":"d","r":"r",
          "eta":"full_Fock_prepared_distance",
          "source_error":"same_kappa_source_insertion_error",
          "V0_bound":"constraint_V0_norm_bound"}
    for own,other in keys.items():
        assert F(bounds[own]["exact"]) == F(ab[other]["exact"])
    return dict(round=1189,review_status="PASS_WITH_EXPLICIT_SCOPE",
                independent_author_code_imports=False,
                domain_argument="analytic shared D(N+1); number compression only for proof limit",
                sources=dict(Ward=ward,hermiticity=herm,norms=snorm,
                             V0_norm=float(np.linalg.norm(v0,2)),
                             Hm_commutator=comm_h,pair_commutator=comm_j),
                exact_bounds=bounds,diagnostic=diag,
                limits=dict(full_GR=False,actual_internal_clock=False,
                            unknown_input="matter6 times total oscillator N<=1, arbitrary R",
                            output="untruncated four-mode Fock",
                            source="same fixed-parameter kappa insertion only"))

def compare(a,b,path=""):
    if isinstance(a,dict):
        assert a.keys()==b.keys(),path
        for k in a:compare(a[k],b[k],path+"/"+k)
    elif isinstance(a,list):
        assert len(a)==len(b),path
        for k,(x,y) in enumerate(zip(a,b)):compare(x,y,path+f"/{k}")
    elif isinstance(a,float):
        assert abs(a-b)<=TOL*max(1,abs(a),abs(b)),(path,a,b)
    else:assert a==b,(path,a,b)

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--write",action="store_true")
    args=parser.parse_args()
    result=run()
    target=ROOT/"review_checks.json"
    if args.write:
        with target.open("x",encoding="utf-8") as f:
            json.dump(result,f,ensure_ascii=False,indent=2);f.write("\n")
    else:
        compare(result,json.loads(target.read_text(encoding="utf-8-sig")))
    print(json.dumps(dict(status=result["review_status"],
                          eta=result["exact_bounds"]["eta"]["display"],
                          source_error=result["exact_bounds"]["source_error"]["display"],
                          diagnostic=result["diagnostic"]),ensure_ascii=False))
