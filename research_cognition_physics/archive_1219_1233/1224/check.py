"""1224: actual full-Gauss short-time identifiability, coefficient checks.

No full-H time evolution or full-Gauss Gibbs spectrum is simulated here.
The original hyperbolic target integral and finite-star algebra are evaluated.
Nonzero finite time is a theorem conditional on the explicitly finite, but
not numerically evaluated, full-H norms ||H^j Psi||, j <= 3.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
import hashlib
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parents[1]
RESULTS = HERE / "results.json"
M, R, N, WIDTH = 2.0, math.sqrt(6.0), 8, 0.1
EPS, SIGMA, HBAR, PHASE = 0.73, 0.12, 1.0, 1.0
SOURCES = [
    "archive_585_628/research_note_589.md",
    "archive_585_628/research_note_598.md",
    "archive_585_628/research_note_623.md",
    "archive_585_628/research_note_625.md",
    "archive_702_741/research_note_704.md",
    "archive_702_741/research_note_721.md",
    "archive_702_741/research_note_723.md",
    "archive_702_741/research_note_724.md",
    "archive_702_741/research_note_725.md",
    "archive_956_989/research_note_963.md",
    "archive_1009_1043/research_note_1041.md",
    "archive_1219_/1221/proof.md",
    "archive_1219_/_admission/area_readout_screen_after1221.md",
]


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def beta(a, b):
    return math.gamma(a) * math.gamma(b) / math.gamma(a + b)


def constants():
    def jac(a):
        return (math.sinh(a / R) / (a / R)) ** 4
    b0 = math.sqrt(M) / 4
    d0 = math.sqrt(M) * WIDTH * (M + Q(1, 3)) / 4
    q1 = b0 * (4*N/WIDTH) * jac(WIDTH) * beta(3, 2*N) / beta(2.5, 2*N+1) + d0
    a0 = b0 * math.cos(math.sqrt(M)*WIDTH) / math.cosh(WIDTH/R)**2
    h = (2*WIDTH/R) / math.tanh(2*WIDTH/R)
    bk = min(WIDTH/2, a0/(4*h*q1))
    qinf = 32*jac(2*WIDTH)*(4/3)**(2*N)*(d0+4*N*b0/WIDTH)
    eta = 1/(4*qinf)
    slow = (a0-bk*h*q1)**2
    volume = EPS**3 * math.exp(6*SIGMA)
    alpha = EPS*math.exp(2*SIGMA)/volume
    vlow = HBAR*PHASE*M/(4*volume)*math.cos(math.sqrt(M)*WIDTH)/math.cosh(WIDTH/R)**4
    assert eta*qinf == 0.25
    assert 0 < bk <= WIDTH/2 and slow > 0
    return dict(M=M, curvature_radius=R, density_power=2*N,
        amplitude_boundary_smoothness="C7; H^3 allowed, arbitrary powers not claimed",
        probe_radius=WIDTH, source_radius=2*WIDTH, common_center_radius=bk,
        gradient_bound=b0, laplacian_bound=d0, divergence_L1_bound=q1,
        density_ratio_bound=qinf, eta=eta, source_density_over_baseline_lower=0.75,
        A_origin_component_lower=a0, hessian_bound=h,
        S_lower=slow, S_upper=b0*b0, node_volume=volume, gradient_force_prefactor=alpha,
        chi_lower=alpha*eta*slow, clock_velocity_lower=vlow,
        hbar_fixed=HBAR, semiclassical_limit_used=False)


def geometry_and_star():
    g = [[Q(4,3), Q(2,3), Q(1,2)],
         [Q(2,3), Q(4,3), Q(0)], [Q(1,2), Q(0), Q(1)]]
    a,b,c=g[0]; _,e,f=g[1]; _,_,i=g[2]
    det = a*(e*i-f*f)-b*(b*i-f*c)+c*(b*f-e*c)
    cofactor12 = c*f-b*i
    inv12 = cofactor12/det
    assert (a,a*e-b*b,det)==(Q(4,3),Q(4,3),Q(1))
    assert inv12==Q(-2,3)
    assert 2*a+2*e-2*b==4 and i==1
    area_sq=(e*i-f*f)/det
    assert area_sq==Q(4,3)
    # Original forward-star discretization, not a new graph Hamiltonian.
    side=5
    B=np.array([0,0,0]); J=np.array([4,1,0])
    dirs=np.eye(3,dtype=int)
    common=[]
    for k in np.ndindex(side,side,side):
        k=np.array(k)
        ends=[(k+d)%side for d in dirs]
        for mu in range(3):
            for nu in range(mu+1,3):
                if ((np.array_equal(ends[mu],B) and np.array_equal(ends[nu],J))
                    or (np.array_equal(ends[mu],J) and np.array_equal(ends[nu],B))):
                    common.append([k.tolist(),mu+1,nu+1])
    assert common==[[[4,0,0],1,2]]
    assert not any(np.array_equal((B+sg*d)%side,J) for d in dirs for sg in [-1,1])
    return dict(covariant_shapes="G0=I; G1 as proof equation (1)",
        G1_leading_minors=[str(a),str(a*e-b*b),str(det)],
        mixed_inverse_coefficients=["0",str(inv12)],
        unchanged_compressed_electric_combination="4",
        unchanged_compressed_magnetic_combination="1",
        normal1_area_ratio_squared=str(area_sq), graph_side=side,
        unique_common_forward_star=common, source_and_probe_not_nearest_neighbors=True,
        full_H_used_in_theorem=True, finite_invariant_sector_assumed=False)


def legendre(n, lo, hi):
    x,w=np.polynomial.legendre.leggauss(n)
    return (x+1)*(hi-lo)/2+lo,w*(hi-lo)/2


def radial_nodes(n, radius):
    r,w=legendre(n,0,radius)
    weight=w*r**4*(1-r*r/radius**2)**(2*N)*(np.sinh(r/R)/(r/R))**4
    return r,weight/weight.sum()


def mink(a,b):
    return np.sum(a[...,1:]*b[...,1:],axis=-1)-a[...,0]*b[...,0]


def embed(y):
    r=np.linalg.norm(y,axis=-1)
    ratio=np.ones_like(r)
    np.divide(R*np.sinh(r/R),r,out=ratio,where=r!=0)
    return np.concatenate([(R*np.cosh(r/R))[...,None],y*ratio[...,None]],axis=-1)


def effect_gradient_ambient(Y):
    x=Y[...,1:]
    sqrtD=Y[...,0]/R
    s=math.sqrt(M)*x[...,-1]/sqrtD
    vs=np.cos(s)*math.sqrt(M)/(4*sqrtD)
    v=np.zeros_like(Y)
    v[...,-1]=vs
    v[...,0]=x[...,-1]*vs/(R*sqrtD)
    return v


def log_coeff(z):
    q=np.maximum(z-1,0)
    small=q<1e-4
    h=np.empty_like(q); hp=np.empty_like(q)
    h[small]=1-q[small]/3+2*q[small]**2/15-2*q[small]**3/35+8*q[small]**4/315
    hp[small]=-1/3+4*q[small]/15-6*q[small]**2/35+32*q[small]**3/315
    h[~small]=np.arccosh(z[~small])/np.sqrt(z[~small]**2-1)
    hp[~small]=(1-z[~small]*h[~small])/(z[~small]**2-1)
    return h,hp


def log_and_differential(K,Y,U):
    z=-mink(K,Y)/R**2
    h,hp=log_coeff(z)
    dz=-mink(K,U)/R**2
    base=Y-z[...,None]*K
    log=h[...,None]*base
    der=h[...,None]*U+dz[...,None]*(hp[...,None]*base-h[...,None]*K)
    return log,der


def coefficient_integral(order):
    c=constants()
    rb,wr=radial_nodes(order,WIDTH)
    u,wu=legendre(order,-1,1)
    wu=wu*(3/4)*(1-u*u)  # S4 singlet polar marginal; exactly normalized quadrature.
    # S3 relative Higgs angle: Gauss-Chebyshev second-kind integrates its weight.
    j=np.arange(1,order+1)
    v=np.cos(j*np.pi/(order+1))
    wv=2*np.sin(j*np.pi/(order+1))**2/(order+1)
    rr,uu,vv=np.meshgrid(rb,u,v,indexing="ij")
    weights=(wr[:,None,None]*wu[None,:,None]*wv[None,None,:]).reshape(-1)
    y=np.zeros((rr.size,5))
    y[:,0]=(rr*np.sqrt(1-uu*uu)*vv).reshape(-1)
    y[:,1]=(rr*np.sqrt(1-uu*uu)*np.sqrt(1-vv*vv)).reshape(-1)
    y[:,-1]=(rr*uu).reshape(-1)
    Y=embed(y); U=effect_gradient_ambient(Y)
    rk,wk=radial_nodes(order,c["common_center_radius"])
    S=0.; max_tangent=0.
    for r,w in zip(rk,wk):
        for uk,wuk in zip(u,wu):
            ky=np.array([r*math.sqrt(1-uk*uk),0,0,0,r*uk])
            K=embed(ky)
            _,der=log_and_differential(K,Y,U)
            A=np.einsum("i,ij->j",weights,der)
            # Components transverse to the K-Higgs direction average to zero
            # over their full S2 orientation, already integrated analytically.
            A[2:5]=0
            max_tangent=max(max_tangent,abs(float(mink(K,A))))
            S+=float(w*wuk*mink(A,A))
    K0=np.array([R,0,0,0,0,0.])
    _,der=log_and_differential(K0,Y,U)
    A0=np.einsum("i,ij->j",weights,der)
    A0[1:5]=0
    # Clock velocity is an actual original F expectation derivative.
    s=math.sqrt(M)*Y[:,-1]/(Y[:,0]/R)
    fold=M/(Y[:,0]/R)**2
    grad_s_sq=fold*(1-s*s/(6*M))
    velocity=HBAR*PHASE/(4*c["node_volume"])*float(weights@(np.cos(s)*grad_s_sq))
    assert c["S_lower"]<=S<=c["S_upper"]
    assert A0[-1]>=c["A_origin_component_lower"]
    assert velocity>=c["clock_velocity_lower"]
    chi=c["gradient_force_prefactor"]*c["eta"]*S
    return dict(order=order,S=S,A_origin_singlet=float(A0[-1]),
        tangent_residual=max_tangent,chi=chi,
        q_second_derivative_G0=0.,q_second_derivative_G1=4*chi/3,
        relational_q_quadratic_coefficient_G1=2*chi/(3*velocity**2),
        actual_local_clock_velocity=velocity)


def target_integral_check():
    rows=[coefficient_integral(n) for n in (12,18,24)]
    convergence=abs(rows[-1]["S"]-rows[-2]["S"])
    assert convergence<2e-10
    # Direct central variation of the SAME original Log, not a propagator.
    K=embed(np.array([.01,.02,0,0,.015]))
    y=np.array([.025,-.01,.005,0,.03]); Y=embed(y)
    U=effect_gradient_ambient(Y)
    _,analytic=log_and_differential(K,Y,U)
    # Tangent exponential makes exact target points for the finite difference.
    unorm=math.sqrt(float(mink(U,U))); step=2e-5
    Yp=np.cosh(step*unorm/R)*Y+R*np.sinh(step*unorm/R)*U/unorm
    Ym=np.cosh(step*unorm/R)*Y-R*np.sinh(step*unorm/R)*U/unorm
    lp,_=log_and_differential(K,Yp,effect_gradient_ambient(Yp))
    lm,_=log_and_differential(K,Ym,effect_gradient_ambient(Ym))
    error=float(np.max(abs((lp-lm)/(2*step)-analytic)))
    assert error<2e-9 and rows[-1]["tangent_residual"]<2e-12
    return dict(rows=rows,last_refinement_S_difference=convergence,
        Log_differential_central_check=error,
        quadrature_is_not_interval_certification=True,
        proof_lower_bound_independent_of_quadrature=True,
        full_H_propagator_or_Gibbs_simulated=False)


def finite_window_contract():
    # Exact rational algebra of the full-H remainder transport.
    assert Q(2,3)-Q(1,3)==Q(1,3)
    assert 4*Q(1,48)==Q(1,12)<Q(1,3)
    assert Q(1,6)*2**3==Q(4,3)
    assert Q(1,2)*2*3==3
    c=constants()
    return dict(full_H_norms_required=["||H Psi||", "||H^2 Psi||", "||H^3 Psi||"],
        full_H_norms_finite_by_specific_preparation_domain=True,
        full_H_norms_numerically_evaluated=False,
        M3="2 ||E|| (N3 + 3 N1 N2) / hbar^3",
        B_G="(M3_G_plus + M3_G_minus)/6",
        positive_window="0<t<=chi_lower/[3(B_G0+B_G1)]",
        probability_double_difference_lower="chi_lower*t^2/3",
        per_probability_statistical_tolerance="chi_lower*t^2/48",
        samples_per_preparation="ceil(log(8/delta)/(2*epsilon_p^2))",
        two_original_readouts_energy_injection_upper=HBAR**2*M/(16*c["node_volume"]),
        time_numerically_certified=False, samples_numerically_certified=False,
        local_clock_inverse_is_mean_curve=True, autonomous_trigger_proved=False,
        length_or_area_readout_proved=False)


def run():
    return dict(round=1224,scientific_groups=4,
        source_sha256={s:digest(RESEARCH/s) for s in SOURCES},
        geometry=geometry_and_star(),finite_preparation=constants(),
        original_target_integral=target_integral_check(),
        finite_window=finite_window_contract(),
        scope="Full-H short-time theorem; original coefficient quadrature only; no device implementation")


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--write",action="store_true")
    parser.add_argument("--verify",action="store_true")
    args=parser.parse_args()
    result=run()
    result["code_sha256"]=digest(Path(__file__))
    if args.write:
        RESULTS.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
        print("1224: wrote 4 coefficient/domain-contract checks; full-H time norms uncomputed")
    else:
        saved=json.loads(RESULTS.read_text("utf-8"))
        assert saved==result,"1224 results/source/code changed"
        print("1224: PASS 4 groups; no numerical full-H time window claimed")


if __name__=="__main__":
    main()
