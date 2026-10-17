"""1269 material source closure; one author group, exact finite algebra.
Default is read-only. --write-results exclusively creates the result file.
No PDE evolution, quantum instrument, or full causal implementation is simulated.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import permutations
import json
from math import isqrt
from pathlib import Path

HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
RESULT=HERE/"material_results.json"

def fs(x):
    x=F(x)
    return str(x.numerator) if x.denominator==1 else f"{x.numerator}/{x.denominator}"

class D:
    def __init__(self,v=0,d=0):
        if isinstance(v,D):
            self.v,self.d=v.v,v.d
        else:
            self.v,self.d=F(v),F(d)
    def __add__(self,o):
        o=D(o); return D(self.v+o.v,self.d+o.d)
    __radd__=__add__
    def __neg__(self): return D(-self.v,-self.d)
    def __sub__(self,o): return self+-D(o)
    def __rsub__(self,o): return D(o)+-self
    def __mul__(self,o):
        o=D(o); return D(self.v*o.v,self.d*o.v+self.v*o.d)
    __rmul__=__mul__
    def __truediv__(self,o):
        o=D(o); return D(self.v/o.v,(self.d*o.v-self.v*o.d)/o.v**2)
    def __rtruediv__(self,o): return D(o)/self
    def __pow__(self,k):
        assert isinstance(k,int)
        if k<0: return D(1)/(self**(-k))
        out=D(1)
        for _ in range(k): out=out*self
        return out
    def sqrt(self):
        assert self.v>0
        a,b=isqrt(self.v.numerator),isqrt(self.v.denominator)
        assert a*a==self.v.numerator and b*b==self.v.denominator
        root=F(a,b)
        return D(root,self.d/(2*root))

def inverse_and_det(mat):
    n=len(mat)
    a=[[D(x) for x in row]+[D(int(i==j)) for j in range(n)]
       for i,row in enumerate(mat)]
    det=D(1)
    for col in range(n):
        pivot=next(i for i in range(col,n) if a[i][col].v)
        if pivot!=col:
            a[col],a[pivot]=a[pivot],a[col]
            det=-det
        p=a[col][col]; det=det*p
        a[col]=[x/p for x in a[col]]
        for i in range(n):
            if i!=col:
                m=a[i][col]
                a[i]=[a[i][j]-m*a[col][j] for j in range(2*n)]
    return [row[n:] for row in a],det

def source_metric_checks():
    signs=[-1,1,1,1]
    U=[F(5,3),F(4,3),F(0),F(0)]
    n=F(2,3); j=[n*x for x in U]
    Ucov=[signs[i]*U[i] for i in range(4)]
    M0=F(11,10); d=F(-1,2000)
    gradients=[
        [F(1,7),F(2,3),F(0),F(0)],
        [F(0),F(1,11),F(3,4),F(0)],
        [F(1,13),F(0),F(0),F(4,5)]]
    S=[[sum(q[i]*q[jj] for q in gradients) for jj in range(4)] for i in range(4)]
    B=sum(signs[i]*S[i][i] for i in range(4))
    assert B>0
    T=[[n*M0*Ucov[i]*Ucov[jj]+2*d*S[i][jj]
        -d*B*(signs[i] if i==jj else 0) for jj in range(4)] for i in range(4)]
    out=[]
    for i in range(4):
        for jidx in range(i,4):
            gi=[[D(signs[a] if a==b else 0,
                    int((a==i and b==jidx) or (a==jidx and b==i)))
                 for b in range(4)] for a in range(4)]
            cov,det=inverse_and_det(gi)
            sqrtg=D(1)/(-det).sqrt()
            nn=(-sum(cov[a][b]*j[a]*j[b] for a in range(4) for b in range(4))).sqrt()/sqrtg
            BB=sum(gi[a][b]*S[a][b] for a in range(4) for b in range(4))
            density=-sqrtg*(nn*M0+d*BB)
            extracted=-density.d*(2 if i==jidx else 1)
            assert extracted==T[i][jidx]
            expected_dn=n*(F(signs[i] if i==jidx else 0)+Ucov[i]*Ucov[jidx])
            if i==jidx: expected_dn/=2
            assert nn.d==expected_dn
            out.append({"component":[i,jidx],"T_lower":fs(extracted),
                        "metric_density_derivative":fs(nn.d)})
    # Same moving jet, exact current variation of n.
    for a in range(4):
        jd=[D(j[b],int(a==b)) for b in range(4)]
        nn=(-sum(signs[b]*jd[b]*jd[b] for b in range(4))).sqrt()
        assert nn.d==-Ucov[a]
        assert (-M0*nn).d==M0*Ucov[a]
    return {"normal_number_density":fs(n),"U_upper":list(map(fs,U)),
            "Bphi":fs(B),"all_ten_metric_components":out,
            "current_variation_chemical_mass":"M0, not M0+d*B/n"}

def exchange_checks():
    signs=[-1,1,1,1]
    U=[F(5,3),F(4,3),F(0),F(0)]
    n=F(2,3); m=F(11,10); B=F(7,5)
    grad=[F(2,7),F(-3,11),F(5,13),F(-7,17)]
    dotd=sum(U[i]*grad[i] for i in range(4))
    dotm=-B*dotd/n
    accel=[B/(n*m)*(signs[i]*grad[i]+U[i]*dotd) for i in range(4)]
    assert sum(signs[i]*U[i]*accel[i] for i in range(4))==0
    for i in range(4):
        lhs=n*(m*accel[i]+dotm*U[i])
        assert lhs==B*signs[i]*grad[i]
    return {"dot_M0":fs(dotm),"acceleration_upper":list(map(fs,accel)),
            "source_divergence_upper":[fs(B*signs[i]*grad[i]) for i in range(4)],
            "total_new_exchange":"0"}

def adm_checks():
    gamma=F(6,5); alpha=F(1,16); C0=F(3,10)
    P=F(-1,1000); g=F(5,7); gp=F(3,11)
    root_h=F(7,3); momentum_square=F(13,17); Bsp=F(19,23)
    out={}
    for name,dot_g,dot_P in [("P",0,1),("theta",gp,0)]:
        gd=D(g,dot_g); pd=D(P,dot_P)
        d=alpha*gd*pd/C0; Gamma=gamma+2*d
        h=momentum_square/(2*root_h*Gamma)+root_h*Gamma*Bsp/2
        Bfull=Bsp-momentum_square/(root_h**2*Gamma.v**2)
        assert h.d==root_h*d.d*Bfull
        out[name]={"Hamiltonian_derivative":fs(h.d),"full_Bphi":fs(Bfull),
                   "spatial_only_difference":fs(h.d-root_h*d.d*Bsp)}
        assert h.d!=root_h*d.d*Bsp
    return out


def hda_local_checks():
    # Local algebra diagnostic for the analytic arbitrary-test-function proof.
    # d is a scalar transported by X; its X derivative contains no lapse gradient.
    N,Nx,M,Mx=map(F,(2,3,5,7))
    Gamma=F(11,10); root_h=F(7,3); Pi=F(5,11); rod_grad=F(13,17)
    divergence_remainder=F(19,23) # includes nonzero Gamma/metric/field gradients
    dHN_q=-Nx*root_h*Gamma*rod_grad-N*divergence_remainder
    dHM_q=-Mx*root_h*Gamma*rod_grad-M*divergence_remainder
    dHN_p=N*Pi/(root_h*Gamma); dHM_p=M*Pi/(root_h*Gamma)
    rod_bracket=dHN_q*dHM_p-dHN_p*dHM_q
    assert rod_bracket==Pi*rod_grad*(N*Mx-M*Nx)
    W=F(3,5); Bfull=F(7,11); grad_d=[F(2,13),F(3,17),F(-5,19)]
    vel=[F(4,7),F(-2,9),F(1,11)]
    mechanical=sum((-N*W*Bfull*grad_d[i])*(M*vel[i])
                   -(-M*W*Bfull*grad_d[i])*(N*vel[i]) for i in range(3))
    assert mechanical==0
    dtheta=F(5,23); M0overE=F(3,5)
    internal=(N*W*Bfull*dtheta)*(M*M0overE)-(M*W*Bfull*dtheta)*(N*M0overE)
    assert internal==0
    assert Gamma*(1/Gamma)==1
    return {"field_lapse_gradient_bracket":fs(rod_bracket),
            "transported_scalar_mechanical_cross":"0",
            "clock_internal_cross":"0",
            "field_kinetic_gradient_coefficient_product":"1",
            "scope":"local check of the same no-caustic label-chart continuum identity"}

def p_add(a,b):
    out=[F(0)]*max(len(a),len(b))
    for i,x in enumerate(a): out[i]+=x
    for i,x in enumerate(b): out[i]+=x
    while len(out)>1 and not out[-1]: out.pop()
    return out
def p_mul(a,b):
    out=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): out[i+j]+=x*y
    while len(out)>1 and not out[-1]: out.pop()
    return out
def charpoly(A):
    n=len(A); out=[F(0)]
    for perm in permutations(range(n)):
        sign=(-1)**sum(perm[i]>perm[j] for i in range(n) for j in range(i+1,n))
        term=[F(sign)]
        for i,j in enumerate(perm):
            term=p_mul(term,[-A[i][j],F(1)] if i==j else [-A[i][j]])
        out=p_add(out,term)
    return out
def mmul(A,B):
    return [[sum(x*y for x,y in zip(row,col)) for col in zip(*B)] for row in A]
def matvec(A,v): return [sum(x*y for x,y in zip(row,v)) for row in A]

def principal_checks():
    cases=[]
    for n in (F(1),F(1,8)):
        for dp in (F(3,7),F(-2,5)):
            gamma=F(11,10); d=F(1,1000); G=gamma+2*d
            s=F(1,3); B=3*s*s; M=F(9,10); dt=F(2,11)
            b=B/(n*M)
            A=[[F(0) for _ in range(6)] for _ in range(6)]
            A[0][1]=A[1][0]=F(1)
            A[0][4]=2*s*dt/G; A[0][5]=2*s*dp/G
            A[2][4]=b*dt; A[2][5]=b*dp
            A[3][2]=-n
            poly=charpoly(A)
            assert poly==[F(0),F(0),F(0),F(0),F(-1),F(0),F(1)]
            # A zero-eigenvalue generalized vector: cancel the rod force using z.
            v=[F(0),-2*s*dt/G,F(0),F(0),F(1),F(0)]
            Av=matvec(A,v); A2v=matvec(A,Av); A3v=matvec(A,A2v)
            assert any(Av) and any(A2v) and not any(A3v)
            # Finite frozen-mode bounds, not actual device parameters.
            q=F(10); T=F(1,100); delta=F(1,1000000)
            vel=b*q*T*delta
            density=n*b*q*q*T*T*delta/2
            assert vel>0 and density>0
            cases.append({"n":fs(n),"d_P":fs(dp),
                          "characteristic_coefficients_low_to_high":list(map(fs,poly)),
                          "Jordan_chain_lengths_verified":[3],
                          "finite_example_q":fs(q),"finite_example_T":fs(T),
                          "finite_example_control_amplitude":fs(delta),
                          "velocity_bound":fs(vel),"density_bound":fs(density)})
    return cases

def homogeneous_checks():
    alpha=F(1,16); C0=F(3,10)
    rows=[]
    for a in (F(1),F(2),F(8,3)):
        for P in (F(1,1000),F(-1,1000)):
            n=a**-3; B=C0/a**2; g=F(7,11); gp=F(5,13)
            d=alpha*g*P/C0; dd=alpha*gp*P/C0
            H=F(3,10); M0=F(93,100)
            rho_f=n*M0; rho_i=d*B; pressure_i=-d*B/3
            dotM=-B*dd/n
            dotrho_f=-3*H*rho_f+n*dotM
            dotrho_i=dd*B-2*H*rho_i
            residual=dotrho_f+dotrho_i+3*H*(rho_f+rho_i+pressure_i)
            assert residual==0
            assert a**3*d*B==alpha*g*P*a
            assert dotM==-alpha*gp*P*a
            phase_rate=B*(alpha*g/C0)/n
            assert phase_rate==alpha*g*a
            # The n-derivative at fixed B,d,J is exactly M0.
            nd=D(n,1)
            per_particle=M0+d*B/nd
            assert (nd*per_particle).d==M0
            rows.append({"a":fs(a),"P_phys":fs(P),"h_int":fs(a**3*d*B),
                         "rho_int":fs(rho_i),"pressure_int":fs(pressure_i),
                         "J_new_rate":fs(dotM),"pointer_phase_rate":fs(phase_rate),
                         "continuity_residual":"0"})
    tau=F(1,10**8); gmax=F(35,16)/tau
    pulse_integral=140*(F(1,4)-3*F(1,5)+3*F(1,6)-F(1,7))
    assert pulse_integral==1 and gmax==218750000
    jcoeff=alpha*(3*gmax+4)
    gamma_coeff=2*alpha*gmax/C0
    assert jcoeff==F(164062501,4) and gamma_coeff==F(273437500,3)
    eps=F(1,2**100)
    assert jcoeff*eps<F(1,2**74)
    assert gamma_coeff*eps<F(1,2**73)
    # Endpoint controller remainder assumes the stated actual |a_dot|<=4 tube.
    assert 4*alpha==F(1,4)
    old=jread(STAGE/"1268/readout_results.json")
    assert old["task"]["N_decimal_digits"]>=3
    assert old["task"]["alpha"]=="1/16" and old["task"]["tau"]=="1/100000000"
    return {"samples":rows,"pulse_integral":fs(pulse_integral),
            "pulse_peak":fs(gmax),"conditional_clock_coefficient":fs(jcoeff),
            "rod_coefficient_perturbation_bound_coefficient":fs(gamma_coeff),
            "sufficient_epsilon_upper":"2^-100",
            "conditional_clock_strict_error":"2^-74",
            "rod_coefficient_strict_error":"2^-73",
            "endpoint_clock_coefficient":"1/4"}

def sha(path): return sha256(path.read_bytes()).hexdigest()
def jread(path): return json.loads(path.read_text(encoding="utf-8-sig"))
def result():
    sources={
        "1269/admission.md":HERE/"admission.md",
        "1269/material_proof.md":HERE/"material_proof.md",
        "1269/material_check.py":Path(__file__).resolve(),
        "1245/proof.md":STAGE/"1245/proof.md",
        "1247/proof.md":STAGE/"1247/proof.md",
        "1249/canonical.md":STAGE/"1249/canonical.md",
        "1268/material_proof.md":STAGE/"1268/material_proof.md",
        "1268/material_results.json":STAGE/"1268/material_results.json",
        "1268/readout_proof.md":STAGE/"1268/readout_proof.md",
        "1268/readout_results.json":STAGE/"1268/readout_results.json"}
    return {"round":1269,"part":"volumetric_internal_current_and_source",
        "status":"PASS_AUTHOR_WITH_SCOPE","author_groups":1,
        "arithmetic":"exact Fraction automatic differentiation and polynomial determinant",
        "metric_variation":source_metric_checks(),"exchange":exchange_checks(),
        "ADM_internal_derivatives":adm_checks(),"matter_HDA_local":hda_local_checks(),
        "principal":principal_checks(),
        "homogeneous_and_resource":homogeneous_checks(),
        "scope":[
            "same conserved label current with a new explicit density-dependent material extension",
            "P is physical action; original logical code uses P=epsilon*sigma_z",
            "exact new source exchange cancels; old source exchanges remain",
            "homogeneous dynamic principal jet only; no negative longitudinal branch in that block",
            "zero-speed Jordan chain retained; finite examples are not a full PDE or device certificate",
            "conditional resource bounds use actual homogeneous a<=3 and |a_dot|<=4",
            "old writing-pulse source and clock contributions are not deleted",
            "Q-only bounded completion agrees with local action on the declared core only",
            "new matter HDA proved in the stated no-caustic label chart; not arbitrary quantum completions",
            "no local quantum pointer/output implementation or general nonlinear PDE theorem"
        ],
        "source_bindings":{k:sha(p) for k,p in sources.items()}}
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--write-results",action="store_true")
    args=ap.parse_args()
    current=result()
    if args.write_results:
        with RESULT.open("x",encoding="utf-8",newline="\n") as f:
            json.dump(current,f,ensure_ascii=False,indent=2,sort_keys=True)
            f.write("\n")
        print("PASS_AUTHOR_WITH_SCOPE FIRST_WRITE_EXCLUSIVE")
    else:
        saved=jread(RESULT)
        if saved!=current:
            keys=sorted(k for k in set(saved)|set(current) if saved.get(k)!=current.get(k))
            raise SystemExit("MISMATCH "+", ".join(keys))
        print("PASS_AUTHOR_WITH_SCOPE READ_ONLY_MATCH")
    print("author_groups=1")
    print("10 metric stress components; current exchange, ADM and matter HDA cross terms: exact")
    print("principal lambda^4*(lambda^2-1); zero-speed Jordan chain retained")
    print("homogeneous source/clock identities and conditional finite bounds: exact")
if __name__=="__main__": main()
