"""1191: correlated Cartan parameters and exact finite-source linear response.
Default checks stored output; --write creates it exclusively.
No SciPy, finite-Fock dynamical truncation, or GR recovery claim.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse,json,math
import numpy as np

def rect_source(N,t=0):
    S=np.zeros((N+3,N+1),complex)
    for n in range(N+1):
        S[n,n]=-(n+.5)/2
        if n>=2:S[n-2,n]=.75*math.sqrt(n*(n-1))*np.exp(-2j*t)
        S[n+2,n]=.75*math.sqrt((n+1)*(n+2))*np.exp(2j*t)
    return S

def source_bound(N):
    return (N+.5)/2+.75*(math.sqrt(max(0,N*(N-1)))+math.sqrt((N+1)*(N+2)))

def kernel(mu,nu,t):
    if nu==0:return 2*math.sin(mu*t/2)**2/mu**2
    if abs(mu-nu)<1e-13:
        return t*math.sin(nu*t)/(2*nu)+1j*(t*math.cos(nu*t)/(2*nu)-math.sin(nu*t)/(2*nu**2))
    return (np.exp(-1j*nu*t)-math.cos(mu*t)+1j*nu*math.sin(mu*t)/mu)/(mu*mu-nu*nu)

def response(N,K,u,t):
    mu=math.sqrt(u/(2*K));I0=kernel(mu,0,t);I2=kernel(mu,2,t)
    B=np.zeros((N+3,N+1),complex)
    for n in range(N+1):
        B[n,n]=-(n+.5)*I0/2
        if n>=2:B[n-2,n]=.75*math.sqrt(n*(n-1))*I2
        B[n+2,n]=.75*math.sqrt((n+1)*(n+2))*I2.conjugate()
    return -B/(2*K)

def cold_root(r):
    lo,hi=1.,max(2.,(2*r/3)**(2/3)+2)
    fun=lambda x:1.5*(x**1.5-x**-2.5)
    for _ in range(100):
        mid=(lo+hi)/2
        if fun(mid)<r:lo=mid
        else:hi=mid
    return (lo+hi)/2

def calculate():
    rows=[]
    for b,v in [(F(1),F(1)),(F(2),F(3,2)),(F(1,3),F(2))]:
        kap=2*b*v;Lam=24/v**2;u=288*b/v**3
        for z in [F(1,10),F(1),F(10)]:
            K=3*v*z*z/(4*b);mass2=u/(2*K);alpha2=1/(8*v*v*K)
            ratio=alpha2/mass2
            assert ratio==1/(24*kap*Lam)
            # Zero-frequency contact coefficient stays fixed along the z path.
            assert alpha2/(2*mass2)==1/(48*kap*Lam)
            rows.append(dict(b=str(b),v=str(v),z=str(z),K=str(K),mass_squared=str(mass2),alpha_squared=str(alpha2),ratio=str(ratio)))
    exact_cold=[]
    for y in [F(1),F(11,10),F(3,2),F(2)]:
        x=y*y;r=F(3,2)*(y**3-y**-5)
        derivative=F(3,4)*(1-x**-4)-r/(2*y**3)
        convexity=3*x**-5+F(3,4)*r*y**-5
        assert derivative==0 and convexity>0
        energy=(x**-3+3*x)/4+r/y
        exact_cold.append(dict(x=str(x),density_ratio=str(r),stationarity=str(derivative),convexity=str(convexity),normalized_energy=str(energy)))
    small=[]
    for r in [1e-2,3e-3,1e-3]:
        x=cold_root(r);e=(x**-3+3*x)/4+r/math.sqrt(x)
        xerr=abs(x-(1+r/6+r*r/36));eerr=abs(e-(1+r-r*r/24))
        assert xerr<r**3 and eerr<r**3
        small.append(dict(r=r,x=x,x_expansion_error=xerr,energy_expansion_error=eerr,energy=e))
    # Independent numerical integral checks the closed retarded kernel,
    # including the resonance that would make an untreated quotient 0/0.
    nodes,weights=np.polynomial.legendre.leggauss(128)
    max_kernel=0.
    for mu in [.1,.7,2.,3.,7.]:
        for nu in [0.,2.]:
            for t in [.2,1.,2.3]:
                ss=(nodes+1)*t/2
                num=np.sum(weights*np.sin(mu*(t-ss))/mu*np.exp(-1j*nu*ss))*t/2
                err=abs(num-kernel(mu,nu,t));max_kernel=max(max_kernel,float(err));assert err<2e-12
    quantum=[]
    for N in [0,1,2,5,12]:
        C=source_bound(N);norm0=float(np.linalg.norm(rect_source(N),2));assert norm0<=C+1e-13
        for t in [.3,1.,2.1]:assert abs(np.linalg.norm(rect_source(N,t),2)-norm0)<2e-13
        quantum.append(dict(initial_max_occupation=N,output_max_occupation=N+2,exact_source_norm=norm0,analytic_upper=C))
    limit_rows=[]
    N=2;C=source_bound(N);u=288.;t=1.
    for z in [1.,10.,100.,1000.]:
        K=.75*z*z;mu=math.sqrt(u/(2*K));B=response(N,K,u,t)
        value=float(np.linalg.norm(B,2));upper=C*t*t/(4*K)
        assert value<=upper+1e-13
        # A joint derivative/energy vector integrates cos and sin, whose
        # stacked action is an exact isometry at each retarded time.
        ss=(nodes+1)*t/2;vel=np.zeros_like(B);pos=np.zeros_like(B)
        for w,tt in zip(weights,ss):
            source=rect_source(N,float(tt));fac=-w*t/(4*K)
            vel+=fac*math.cos(mu*(t-tt))*source
            pos+=fac*math.sin(mu*(t-tt))*source
        energy_matrix=K*(vel.conj().T@vel+pos.conj().T@pos)
        energy=float(np.linalg.eigvalsh(energy_matrix)[-1]);energy_upper=C*C*t*t/(4*K)
        assert np.linalg.norm(pos/mu-B)<2e-13
        assert energy<=energy_upper+1e-12
        limit_rows.append(dict(z=z,K=K,mass=mu,T=t,full_source_response_norm=value,analytic_upper=upper,induced_linear_energy_upper=energy_upper,computed_energy_norm=energy,zero_frequency_susceptibility=1/u))
    # Fixed nonzero k suppressed at large K; k=0 static response unchanged.
    static=[dict(z=z,k=k,susceptibility=1/(288+1.5*z*z*k*k)) for z in [.1,1.,10.,100.] for k in [0.,1.]]
    return dict(status='PASS_WITH_SCOPE',check_groups=4,parameter_identity=rows,cold_exact=exact_cold,cold_small_density=small,retarded_kernel_max_error=max_kernel,full_Fock_source=quantum,finite_time=limit_rows,static_response=static,scope={'full_source_maps_into':'entire Fock; quadratic source gives exact output N<=input N+2','reference':'operator norm tensor-stable for arbitrary R','dynamics':'free-source linear retarded response only; no coupled geometry/material error certified','initial_geometry':'zero induced component; homogeneous/free component must be retained','cold_branch':'separate diagonal adiabatic fixed-density task; not arbitrary unknown state','GR_or_dark_sector_generated':False})

def equal(a,b):
    if isinstance(a,float):return isinstance(b,(int,float)) and math.isclose(a,b,rel_tol=2e-9,abs_tol=2e-12)
    if isinstance(a,dict):return a.keys()==b.keys() and all(equal(a[k],b[k]) for k in a)
    if isinstance(a,list):return len(a)==len(b) and all(equal(x,y) for x,y in zip(a,b))
    return a==b

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    data=calculate();path=Path(__file__).with_name('results.json')
    if args.write:
        with path.open('x',encoding='utf8') as f:json.dump(data,f,ensure_ascii=False,indent=2);f.write('\n')
    else:assert equal(data,json.loads(path.read_text('utf8'))),'saved results differ'
    print(json.dumps(dict(status=data['status'],check_groups=4,N2_z1000_response=data['finite_time'][-1]['full_source_response_norm'],N2_z1000_bound=data['finite_time'][-1]['analytic_upper'],coupled_GR_certificate=False)))
if __name__=='__main__':main()
