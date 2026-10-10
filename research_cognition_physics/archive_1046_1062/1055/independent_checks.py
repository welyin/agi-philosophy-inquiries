"""1055 independent calibration: no imports from author modules."""
import json, math, argparse
from fractions import Fraction as F
from pathlib import Path
import numpy as np

def maxabs(x):
    return float(np.max(np.abs(x)))

def run():
    t=F(1,100); m=10**20; G1=10**6; G2=10**12
    A0=F(5,4)*G2+3*G1+2+4*t*G1+3*t+3*t*t+4*t
    A1=F(5,4)*G2+11*G1+19
    C=F(1401,42875000)
    Jerror=3*8*t*t*10**12/m+18*C/10**6
    prob=F(4,10**10)+24*t*t/(49*10**6)
    kappa=F(589541,39711168000000)
    margin=kappa-6*C/10**5-Jerror
    Adir=F(5,4)*G2+3*G1+2+4*G1+3+3+4
    edir=2*(F(2,10**8)+F(2,7)*F(13,10**6))
    gap=F(6889,4000000)-F(3,200000)
    vals=dict(A0=A0,A1=A1,Jacobian_error=Jerror,probability_error=prob,
              common_inverse_slope=margin,A0_direction=Adir,
              direction_effect_error=edir,direction_gap=gap)
    assert A0<2*10**12 and A1<2*10**12 and Adir<2*10**12
    assert Jerror<F(1,10**9) and prob<F(1,10**9) and margin>F(1,10**8)
    assert edir<F(3,400000) and gap>F(17,10000)
    assert 6*(F(1,2)-F(1,3))==1
    assert 18*(F(1,3)-F(1,2)+F(1,5))==F(3,5)

    # Independent direct PDE residual for the normal driven Gaussian.
    mu=1.7; sigma=.4; d=np.array([.21,-.13,.08])
    xyz=[np.array(x) for x in [(0,0,0),(.1,.3,-.2),(-.7,.4,.6)]]
    residuals=[]
    for s in [.07,.31,.72,1.0]:
        force=6*mu*d*(1-2*s)
        impulse=6*mu*d*(s-s*s)
        b=d*(3*s*s-2*s**3)
        bdot=impulse/mu
        phidot=float(impulse@impulse)/(2*mu)
        zeta=1+1j*s/(2*mu*sigma*sigma)
        zdot=1j/(2*mu*sigma*sigma)
        for r in xyz:
            z=r-b
            grad=-z/(2*sigma*sigma*zeta)
            free_dt=-1.5*zdot/zeta+float(z@z)*zdot/(4*sigma*sigma*zeta*zeta)
            log_dt=1j*float(force@r)-1j*phidot+free_dt-float(bdot@z)*(-1/(2*sigma*sigma*zeta))
            lap_ratio=np.sum((1j*impulse+grad)**2)-3/(2*sigma*sigma*zeta)
            lhs=1j*log_dt
            rhs=-lap_ratio/(2*mu)-float(force@r)
            residuals.append(float(abs(lhs-rhs)))
    assert max(residuals)<2e-13

    # Physical momenta versus reference representation: no symbolic K reused.
    rng=np.random.default_rng(1055)
    recoil=[]
    for _ in range(30):
        P,p,k=rng.normal(size=(3,3))
        k=k/np.linalg.norm(k)*rng.uniform(1,2)
        px=P/2+p; py=P/2-p
        laboratory=(px@px+(py-k)@(py-k))/2
        source=((P-k)@(P-k))/4+(p+k/2)@(p+k/2)
        recoil.append(float(abs(laboratory-source)))
    assert max(recoil)<2e-14

    # Paired angular quadrature solely checks block signs and source normalization.
    # It is not a proof of continuum rotational covariance or of the error bound.
    pauli=np.array([[[0,1],[1,0]],[[0,-1j],[1j,0]],[[1,0],[0,-1]]],complex)
    modes=[]; pair=[]
    for radius in [1.2,1.7]:
        for axis in range(3):
            eps=[np.eye(3)[j] for j in range(3) if j!=axis]
            first=len(modes)
            for sign in [1,-1]:
                k=sign*radius*np.eye(3)[axis]
                for e in eps:
                    modes.append((k.copy(),e.copy()))
            pair.extend([first+2,first+3,first,first+1])
    n=len(modes); ng=2*n; dim=ng+2
    coeff=1/math.sqrt(12)
    J=np.zeros((3*dim,2),complex)
    for j in range(3):
        for q in range(2):
            for l,(k,e) in enumerate(modes):
                J[j*dim+q*n+l,q]=coeff*e[j]/math.sqrt(2)
    PF=np.zeros((n,n),complex)
    for l,j in enumerate(pair):PF[l,j]=-1
    parity=np.zeros((dim,dim),complex)
    parity[:ng,:ng]=np.kron(np.eye(2),PF)
    parity[ng:,ng:]=-np.eye(2)
    def optical(r):
        A=np.zeros((2,ng),complex)
        for q in range(2):
            for l,(k,e) in enumerate(modes):
                se=np.tensordot(e,pauli,axes=1)
                A[:,q*n+l]=coeff*np.exp(1j*(k@r))*se[:,q]
        H=np.diag(np.r_[np.tile([np.linalg.norm(k) for k,e in modes],2),[1.5,1.5]]).astype(complex)
        H[ng:,:ng]=.1*A;H[:ng,ng:]=.1*A.conj().T
        return H
    def kinetic(P,p):
        ks=np.array([k for k,e in modes])
        diag=np.sum((P-ks)**2,axis=1)/4+np.sum((p+ks/2)**2,axis=1)
        return np.diag(np.r_[np.tile(diag,2),[P@P/4+p@p]*2])
    r=np.array([.14,-.21,.08]);P=np.array([.2,-.1,.3]);p=np.array([.1,.3,-.4])
    parity_H=maxabs(parity@optical(r)@parity-optical(-r))
    parity_K=maxabs(parity@kinetic(P,p)@parity-kinetic(-P,-p))
    inject=maxabs(J.conj().T@J-np.eye(2))
    source_parity=maxabs(np.kron(np.eye(3),parity)@J+J)
    def effect(r,P,p):
        H=optical(r)+kinetic(P,p)/3.4
        w,V=np.linalg.eigh(H)
        U=(V*np.exp(-.17j*w))@V.conj().T
        out=np.kron(np.eye(3),U)@J
        idx=np.concatenate([np.arange(j*dim+ng,(j+1)*dim) for j in range(3)])
        return out[idx].conj().T@out[idx]
    parity_E=maxabs(effect(r,P,p)-effect(-r,-P,-p))
    parity_errors=dict(H_block=parity_H,K_block=parity_K,
                       injection_normalization=inject,
                       source_parity=source_parity,effect=parity_E)
    assert max(parity_errors.values())<2e-13
    return {
      "scope":"Independent calibration only; continuum domain, uniform reference quantifiers and derivative bounds require the analytic review.",
      "author_modules_imported":False,
      "fraction_bounds":{k:{"exact":str(v),"float":float(v)} for k,v in vals.items()},
      "gaussian_driven_PDE":{"points":len(residuals),"max_residual":max(residuals)},
      "source_recoil":{"samples":len(recoil),"max_residual":max(recoil)},
      "paired_mode_parity":{"mode_count":n,"optical_dimension":dim,"residuals":parity_errors},
      "all_checks_passed":True
    }

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--write",action="store_true",help="exclusively create independent result")
    args=parser.parse_args()
    result=run()
    output=Path(__file__).with_name("independent_checks_results.json") if "__file__" in globals() and __file__!="<stdin>" else None
    if args.write:
        if output is None:raise RuntimeError("save script first")
        with output.open("x",encoding="utf-8") as f:
            json.dump(result,f,ensure_ascii=False,indent=2);f.write("\n")
    elif output is not None and output.exists():
        old=json.loads(output.read_text(encoding="utf-8"))
        assert old["fraction_bounds"]==result["fraction_bounds"]
        assert old["all_checks_passed"] is True
    print(json.dumps(result,ensure_ascii=False,indent=2))
