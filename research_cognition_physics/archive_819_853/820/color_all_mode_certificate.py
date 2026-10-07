"""820: all integer Fourier modes of the fixed-geometry color response.
High modes have an analytic subspace-gap bound; the 18 remaining nonzero
modes have exact rational Hermitian LDL certificates. The diagnostic full
source solve retains the original 753 conformal metric and all configurations.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse,json,sys,itertools
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
import color_restricted_symbol as probe
previous=probe.previous;bg=probe.bg;geo=previous.geo
TARGET=HERE/'color_all_mode_certificate_results.json'
ZERO=(Q(0),Q(0))
def add(a,b):return (a[0]+b[0],a[1]+b[1])
def sub(a,b):return (a[0]-b[0],a[1]-b[1])
def mul(a,b):return (a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def conjugate(a):return (a[0],-a[1])
def summed(values):
    out=ZERO
    for v in values:out=add(out,v)
    return out
def scale(a,c):return (a[0]*c,a[1]*c)
def matrix_product(a,b):
    return [[summed(mul(a[i][k],b[k][j]) for k in range(3)) for j in range(3)] for i in range(3)]
def matrix_scale(a,c):return [[scale(z,c) for z in row] for row in a]
def cross(a,b):
    ab=matrix_product(a,b);ba=matrix_product(b,a)
    return [[mul((Q(0),Q(1)),sub(ab[i][j],ba[i][j])) for j in range(3)] for i in range(3)]
def inner(a,b):
    return scale(summed(mul(conjugate(a[i][j]),b[i][j]) for i in range(3) for j in range(3)),Q(2))
def exact_basis():
    pairs=[]
    for i,j in ((0,1),(0,2),(1,2)):
        x=[[ZERO]*3 for _ in range(3)];y=[[ZERO]*3 for _ in range(3)]
        x[i][j]=x[j][i]=(Q(1,2),Q(0))
        y[i][j]=(Q(0),Q(-1,2));y[j][i]=(Q(0),Q(1,2));pairs.extend((x,y))
    d3=[[ZERO]*3 for _ in range(3)];d8=[[ZERO]*3 for _ in range(3)]
    d3[0][0]=(Q(1,2),Q(0));d3[1][1]=(Q(-1,2),Q(0))
    for i,value in enumerate((Q(1,2),Q(1,2),Q(-1))):d8[i][i]=(value,Q(0))
    return [pairs[0],pairs[1],d3,*pairs[2:],d8]
def exact_symbol(k):
    basis=exact_basis();norm=[Q(1)]*7+[Q(3)]
    A=[matrix_scale(basis[j],Q(v)) for j,v in zip((0,1,3),('.31','.27','.21'))]
    E=[matrix_scale(basis[j],Q(v)) for j,v in zip((0,1,3),('.23','.19','.17'))]
    B=[[ZERO for j in range(24)] for i in range(12)]
    for i in range(3):
        for j in range(8):
            action=cross(A[i],basis[j])
            for a in range(8):
                value=scale(inner(basis[a],action),1/norm[a])
                B[a][8*i+j]=add(value,(Q(0),Q(k[i]) if a==j else Q(0)))
            for a in range(3):B[8+a][8*i+j]=inner(cross(A[a],A[i]),basis[j])
            B[11][8*i+j]=scale(inner(E[i],basis[j]),Q(1,10))
    return B

def exact_ldl(k):
    B=probe.symbol(k).copy();B[-1]*=.1
    # Use diag(1,...,1,sqrt(3)) to replace normalized T8 by a rational basis.
    S=np.ones(8);S[-1]=np.sqrt(3)
    B[:8]/=S[:,None];B*=np.tile(S,3)[None,:]
    b=exact_symbol(k)
    floating=np.array([[float(z[0])+1j*float(z[1]) for z in row] for row in b])
    assert np.max(abs(B-floating))<1e-14
    gram=[[summed(mul(b[i][a],conjugate(b[j][a])) for a in range(24))
           for j in range(12)] for i in range(12)]
    L=[[ZERO for j in range(12)] for i in range(12)];pivots=[]
    for i in range(12):
        d=sub(gram[i][i],summed(mul(mul(L[i][a],conjugate(L[i][a])),(pivots[a],Q(0))) for a in range(i)))
        assert d[1]==0 and d[0]>0,(k,i,d)
        pivots.append(d[0]);L[i][i]=(Q(1),Q(0))
        for j in range(i+1,12):
            val=sub(gram[j][i],summed(mul(mul(L[j][a],conjugate(L[i][a])),(pivots[a],Q(0))) for a in range(i)))
            L[j][i]=(val[0]/d[0],val[1]/d[0])
    return dict(k=list(k),positive_exact_pivots=[str(x) for x in pivots])

def analytic_bounds():
    # F12,F13,F23 have three orthogonal color generators.
    c=[Q('0.31')*Q('0.27'),Q('0.31')*Q('0.21')/2,Q('0.27')*Q('0.21')/2]
    e=[Q('.23'),Q('.19'),Q('.17')]
    delta_squared=min(min(v*v for v in c),(e[1]**2+e[2]**2)/100)
    R_squared=max(c[0]**2+c[1]**2,c[0]**2+c[2]**2,c[1]**2+c[2]**2,
                  sum(x*x for x in e)/100)
    connection_squared=Q('0.31')**2+Q('0.27')**2+Q('0.21')**2
    # ||ad A_i|| <= a_i in this orthonormal su(3) convention.
    assert delta_squared>R_squared*connection_squared/3
    gap=(np.sqrt(float(delta_squared))-np.sqrt(float(R_squared*connection_squared/3)))**2
    A,E,F,C=probe.data()
    assert np.linalg.norm(np.hstack(C),2)**2<=float(connection_squared)+1e-14
    assert abs(np.linalg.norm(np.vstack((F.reshape(3,24),.1*E.reshape(1,24))),2)**2-float(R_squared))<1e-14
    return dict(transverse_min_squared=str(delta_squared),response_norm_squared=str(R_squared),
                connection_norm_squared_upper=str(connection_squared),
                strict_rational_margin=str(delta_squared-R_squared*connection_squared/3),
                high_mode_radius_squared_at_least=3,
                reduced_Gram_uniform_lower=gap)

def solve_color(g,m,energy):
    """Solve G e=g, M e=m, E_bg dot e=energy, on the original periodic slice."""
    A,E,F,C=probe.data();N=g.shape[0];waves=geo.waves(N)
    target=np.concatenate((g,m,energy[...,None]),axis=-1)
    coeff=geo.fft(target);out=np.zeros(g.shape[:-1]+(24,),complex)
    for index in np.ndindex((N,N,N)):
        k=waves[index];B=probe.symbol(k)
        if not np.any(k):
            assert np.max(abs(coeff[index][8:11]-A@coeff[index][:8]))<1e-9
            reduced=np.vstack((B[:8],B[11:12]));rhs=np.r_[coeff[index][:8],coeff[index][11]]
            out[index]=reduced.conj().T@np.linalg.solve(reduced@reduced.conj().T,rhs)
        else:
            scaled=B.copy();scaled[11]*=.1;rhs=coeff[index].copy();rhs[11]*=.1
            out[index]=scaled.conj().T@np.linalg.solve(scaled@scaled.conj().T,rhs)
    result=np.fft.ifftn(out,axes=(0,1,2))
    assert np.max(abs(result.imag))<1e-10
    return result.real.reshape((N,N,N,3,8))

def fixed_geometry_diagnostic():
    d=previous.setup();q=d['q'];psi=d['psi']
    A,E,F,C=probe.data()
    # Reuse the old EW particular and its zero-mean total-momentum repair.
    # Original prescribed color source has total color in T3,T8, orthogonal
    # to A_i in T1,T2,T4, so the required zero-mode M=A G is zero here.
    partial_m=previous.old.momentum(q,d['dp'],d['dE'],d['dE0'])
    g=-d['sigma_c'];m=-d['J']-partial_m
    bc,bw,b0=geo.old.PAR['b']
    partial_energy=(psi**-12*np.sum(previous.old.kinverse(q,q['p'])*d['dp'],axis=-1)
                   +2*psi**-8*(bw*np.sum(q['f']['E']*d['dE'],axis=(-1,-2))
                                +b0*np.sum(q['f']['E0']*d['dE0'],axis=-1)))
    energy=-(d['rho']+partial_energy)*psi**8/(2*bc)
    ec=solve_color(g,m,energy)
    charge=sum(previous.D(ec[...,i,:],i,C) for i in range(3))
    mom=np.einsum('...ja,ija->...i',ec,F)
    colored=2*bc*psi**-8*np.sum(E*ec,axis=(-1,-2))
    residuals=dict(color_Gauss=float(np.max(abs(charge-g))),momentum=float(np.max(abs(mom-m))),
                   physical_energy=float(np.max(abs(colored+partial_energy+d['rho']))))
    assert max(residuals.values())<1e-10
    quadratic=(.5*psi**-12*np.sum(d['dp']*previous.old.kinverse(q,d['dp']),axis=-1)
               +psi**-8*(bw*np.sum(d['dE']**2,axis=(-1,-2))
                           +b0*np.sum(d['dE0']**2,axis=-1)+bc*np.sum(ec**2,axis=(-1,-2))))
    assert quadratic.min()>=0
    rows=[]
    for eps in (.02,.01,.005):
        # Fixed configuration makes Gauss and momentum exactly linear.
        # Uncorrected energy is quadratic; this is not an exact finite solve.
        error=float(np.max(abs(eps**2*quadratic)))
        rows.append(dict(epsilon=eps,uncorrected_energy_residual=error,
                         residual_over_epsilon_squared=float(np.max(quadratic))))
    return dict(N=q['N'],residuals=residuals,geometry_and_all_configurations_fixed=True,
                actual819_source_not_numerically_supplied=True,
                declared_original754_diagnostic_source=True,
                color_electric_correction_max=float(np.max(abs(ec))),rows=rows,
                exact_finite_strength_constraints_claimed=False)

def run():
    low=[exact_ldl(k) for k in itertools.product(range(-1,2),repeat=3)
         if sum(x*x for x in k) in (1,2)]
    assert len(low)==18
    return dict(round=820,all_checks_passed=True,analytic_high_mode_bounds=analytic_bounds(),
                exact_low_mode_certificates=low,fixed_geometry_diagnostic=fixed_geometry_diagnostic(),
                full_smooth_linear_fixed_geometry_compensation_proven=True,
                actual819_source_on_transition_slice_compensated=False,
                autonomous_preparation_proven=False)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args();r=run()
    if args.write:
        assert not TARGET.exists()
        TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps({k:v for k,v in r.items() if k!='exact_low_mode_certificates'},ensure_ascii=False,indent=2))
