"""1191 internal independent implementation. Does not import author check.py.
Default is read-only; --write exclusively creates review_checks.json.
Uses q,p construction plus an augmented ODE exponential, not author kernels.
"""
import argparse
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import math
import numpy as np

ROOT = Path(__file__).resolve().parent
AUTHOR_FILES = ("proof.md","sources.md","admission.md","check.py","results.json","baseline.json")

def exp_small(A):
    """Scaling/squaring Taylor exponential; handles defective resonance matrices."""
    size = A.shape[0]
    nrm = float(np.linalg.norm(A, 1))
    power = max(0, math.ceil(math.log2(nrm))) if nrm else 0
    B = A / 2**power
    result = np.eye(size, dtype=complex)
    term = result.copy()
    for j in range(1,90):
        term = term @ B / j
        result += term
        if np.linalg.norm(term, 1) < 2e-17:
            break
    else:
        raise AssertionError("Exponential series did not converge")
    for _ in range(power):
        result = result @ result
    return result

def forced(mu, nu, t):
    # y''+mu^2 y = exp(i nu t), y(0)=y'(0)=0.
    G = np.array([[0,1,0],[-mu*mu,0,1],[0,0,1j*nu]], complex)
    state = exp_small(t*G) @ np.array([0,0,1],complex)
    return state[0],state[1]

def source_from_qp(N,t=0):
    # At least two unused output levels; then restriction to initial P_N is exact.
    dim = N+5
    a = np.zeros((dim,dim),complex)
    for j in range(1,dim):
        a[j-1,j] = math.sqrt(j)
    q = (a+a.conj().T)/math.sqrt(2)
    p = (a-a.conj().T)/(1j*math.sqrt(2))
    S = .5*(q@q)-(p@p)
    phases = np.exp(1j*t*(np.arange(dim)+.5))
    evolved = phases[:,None]*S*phases.conj()[None,:]
    assert np.linalg.norm(evolved[N+3:,:N+1]) < 1e-13
    return evolved[:N+3,:N+1]

def Cbound(N):
    lower = math.sqrt(N*(N-1)) if N>=2 else 0.0
    return (N+.5)/2+.75*(lower+math.sqrt((N+1)*(N+2)))

def response_ode(N,K,u,t,k=0.):
    source = source_from_qp(N)
    mu = math.sqrt(k*k+u/(2*K))
    response = np.zeros_like(source)
    velocity = np.zeros_like(source)
    cache = {}
    for row in range(source.shape[0]):
        for col in range(source.shape[1]):
            if abs(source[row,col])>1e-15:
                nu = row-col
                if nu not in cache:
                    cache[nu] = forced(mu,nu,t)
                integ,deriv = cache[nu]
                response[row,col] = -source[row,col]*integ/(2*K)
                velocity[row,col] = -source[row,col]*deriv/(2*K)
    return response,velocity,mu

def calculate():
    author = json.loads((ROOT/"results.json").read_text(encoding="utf-8"))
    parameter_rows = []
    for b,v in ((F(1),F(1)),(F(5,7),F(4,3)),(F(3),F(2,5))):
        for z in (F(1,7),F(2),F(17)):
            K=3*v*z*z/(4*b)
            upp=288*b/v**3
            kap=2*b*v
            cosm=24/v**2
            mass2=upp/(2*K)
            alpha2=1/(8*v*v*K)
            assert alpha2/mass2 == 1/(24*kap*cosm)
            parameter_rows.append({"b":str(b),"v":str(v),"z":str(z),"ratio":str(alpha2/mass2)})
    # Exact independent derivative checks of the cold energy.
    for row in author["cold_exact"]:
        x=F(row["x"]);r=F(row["density_ratio"])
        yn=math.isqrt(x.numerator);yd=math.isqrt(x.denominator)
        assert yn*yn==x.numerator and yd*yd==x.denominator
        y=F(yn,yd)
        derivative=F(3,4)*(1-x**-4)-r/(2*y**3)
        second=3*x**-5+F(3,4)*r*y**-5
        assert derivative==0 and second>0
        energy=(x**-3+3*x)/4+r/y
        assert energy==F(row["normalized_energy"])
    cold_residual=0.
    for row in author["cold_small_density"]:
        x=row["x"];r=row["r"]
        residual=abs(.75*(1-x**-4)-.5*r*x**-1.5)
        cold_residual=max(cold_residual,residual)
        assert residual<2e-13
    # Augmented matrix exponential includes zero frequency and resonance.
    kernel_error=0.
    for mu in (0.,.1,2.,7.):
        for nu in (-2.,0.,2.):
            for t in (.2,1.,2.3):
                y,dy=forced(mu,nu,t)
                assert abs(y)<=t*t/2+2e-12
                assert abs(dy)<=t+2e-12
                if nu==0:
                    ref=t*t/2 if mu==0 else 2*math.sin(mu*t/2)**2/mu**2
                    kernel_error=max(kernel_error,abs(y-ref))
                if mu==2 and abs(nu)==2:
                    sign=1 if nu>0 else -1
                    ref=t*math.sin(mu*t)/(2*mu)+1j*sign*(math.sin(mu*t)/(2*mu*mu)-t*math.cos(mu*t)/(2*mu))
                    kernel_error=max(kernel_error,abs(y-ref))
    assert kernel_error<2e-12
    source_rows=[]
    for row in author["full_Fock_source"]:
        N=row["initial_max_occupation"]
        source=source_from_qp(N)
        norm=float(np.linalg.norm(source,2))
        assert abs(norm-row["exact_source_norm"])<2e-12
        assert norm<=Cbound(N)+2e-12
        for t in (.17,.73,2.4):
            assert abs(np.linalg.norm(source_from_qp(N,t),2)-norm)<2e-12
        # Numerical tensor tests supplement, but do not prove, arbitrary-R identity.
        for rdim in (2,3):
            assert abs(np.linalg.norm(np.kron(source,np.eye(rdim)),2)-norm)<3e-12
        source_rows.append({"N":N,"norm_from_qp":norm,"bound":Cbound(N)})
    s0=source_from_qp(0)
    assert abs(np.linalg.norm(s0)-math.sqrt(19)/4)<1e-14
    assert abs(s0[2,0]-.75*math.sqrt(2))<1e-14
    time_rows=[]
    max_author_response_error=0.
    max_author_energy_error=0.
    for row in author["finite_time"]:
        K=row["K"];t=row["T"];N=2
        D,V,mu=response_ode(N,K,288.,t)
        response_norm=float(np.linalg.norm(D,2))
        energy_matrix=K*(V.conj().T@V+mu*mu*(D.conj().T@D))
        energy=float(np.linalg.eigvalsh(energy_matrix)[-1])
        response_error=abs(response_norm-row["full_source_response_norm"])
        energy_error=abs(energy-row["computed_energy_norm"])
        max_author_response_error=max(max_author_response_error,response_error)
        max_author_energy_error=max(max_author_energy_error,energy_error)
        assert response_error<2e-12 and energy_error<2e-12
        assert response_norm<=Cbound(N)*t*t/(4*K)+2e-12
        assert energy<=Cbound(N)**2*t*t/(4*K)+2e-12
        time_rows.append({"z":row["z"],"response_ode":response_norm,"energy_ode":energy})
    # Genuine resonance in the response: mu=sqrt(288/(2K))=2.
    resonance=[]
    for N in (0,2):
        for t in (.5,1.,2.3):
            K=36.
            D,V,mu=response_ode(N,K,288.,t)
            assert abs(mu-2)<1e-14
            value=float(np.linalg.norm(D,2))
            upper=Cbound(N)*t*t/(4*K)
            assert value<=upper+2e-12
            if N==0:
                f2,_=forced(mu,2.,t)
                assert abs(abs(D[2,0])-3*math.sqrt(2)*abs(f2)/(8*K))<2e-14
            resonance.append({"N":N,"t":t,"norm":value,"bound":upper})
    # Multiple spatial eigenmodes check the L2 and joint energy bounds.
    spatial=[]
    weights=np.array([.4,-.7,1.1])
    for K in (36.,7500.):
        blocks=[];vblocks=[];eblocks=[]
        for weight,k in zip(weights,(0.,1.3,7.)):
            D,V,mu=response_ode(2,K,288.,1.2,k)
            blocks.append(weight*D)
            vblocks.append(weight*V)
            eblocks.append(weight*mu*D)
        total=np.vstack(blocks)
        joint=np.vstack(vblocks+eblocks)
        M=float(np.linalg.norm(weights))*Cbound(2)
        assert np.linalg.norm(total,2)<=M*1.2**2/(4*K)+2e-12
        energy=K*np.linalg.norm(joint,2)**2
        assert energy<=M*M*1.2**2/(4*K)+2e-12
        spatial.append({"K":K,"response_norm":float(np.linalg.norm(total,2)),"energy":float(energy)})
    hashes={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in AUTHOR_FILES}
    return {
        "status":"PASS_INDEPENDENT_IMPLEMENTATION_WITH_SCOPE",
        "implementation":"q,p quadratic source and 3x3 augmented ODE matrix exponential; no import of author check.py",
        "author_files_sha256":hashes,
        "parameter_identity":parameter_rows,
        "cold_max_stationary_residual":cold_residual,
        "kernel_anchor_max_error":float(kernel_error),
        "source_norms":source_rows,
        "finite_time":time_rows,
        "max_author_response_error":max_author_response_error,
        "max_author_energy_error":max_author_energy_error,
        "resonant_response":resonance,
        "multiple_spatial_modes":spatial,
        "scope":{
            "full_Fock_output":"single quadratic source/linear response exact from P_N to P_(N+2)",
            "arbitrary_R":"proved tensor norm identity; finite examples supplementary",
            "no_coupled_material_geometry_error_certificate":True,
            "no_actual_clock_or_GR_certificate":True,
            "not_external_peer_review":True
        }
    }

def compatible(a,b):
    if isinstance(a,float):
        return isinstance(b,(float,int)) and math.isclose(a,b,rel_tol=3e-9,abs_tol=3e-12)
    if isinstance(a,dict):
        return a.keys()==b.keys() and all(compatible(a[k],b[k]) for k in a)
    if isinstance(a,list):
        return len(a)==len(b) and all(compatible(x,y) for x,y in zip(a,b))
    return a==b

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--write",action="store_true")
    args=parser.parse_args()
    result=calculate()
    out=ROOT/"review_checks.json"
    if args.write:
        with out.open("x",encoding="utf-8") as f:
            json.dump(result,f,ensure_ascii=False,indent=2)
            f.write("\n")
    else:
        assert compatible(result,json.loads(out.read_text(encoding="utf-8"))),"Review result or author hashes changed"
    print(json.dumps({"status":result["status"],"response_error":result["max_author_response_error"],"energy_error":result["max_author_energy_error"],"coupled_GR_certificate":False}))
if __name__=="__main__":
    main()
