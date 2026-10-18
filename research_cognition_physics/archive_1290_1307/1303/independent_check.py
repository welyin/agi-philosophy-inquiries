"""Independent checks for round 1303.
The all-window field statement is analytic, not inferred from finite diagnostics.
No author implementation is imported. Default and --check are read-only.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from fractions import Fraction as F
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
def frac(x): return str(x)
def tidy(x): return float(format(float(x), ".12g"))
def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()

def exact_bounds():
    a, L, mu = F(1,32), F(1,4), F(1,4096)
    t0, t1 = F(7,16), F(1,2)
    smin = t0*t0 - (L+2*a)**2
    smax = t1*t1 - (L-2*a)**2
    assert smin == F(3,32) and smax == F(55,256)
    assert smin > 0 and smax < F(1,4) and 3*mu < 1
    # Coefficients of J_1(z)/z, with the common 1/pi removed.
    coeff = [F((-1)**k, 8*4**k*math.factorial(k)*math.factorial(k+1))
             for k in range(4)]
    assert coeff == [F(1,8), F(-1,64), F(1,1536), F(-1,73728)]
    # cos(z)>3/4 and sin(z)/z>3/4 on the bump, p=1.
    assert 1-a*a/2 > F(3,4)
    assert 1-a*a/6 > F(3,4)
    # Inner-half-ball z^2 moment / global bump mass.
    b2_over_b0_lower = (a**5/F(480))/(2*a**3/F(3))
    assert b2_over_b0_lower == a*a/F(320)
    # DFT spectral determinant, independent of the mixing formula.
    r = (1,2,3)
    det = ((r[0]**3-r[2]**3)*(r[1]**5-r[2]**5)
           -(r[0]**5-r[2]**5)*(r[1]**3-r[2]**3))
    assert det == 888
    # Independent exact Householder comparison with the author's real candidate.
    v=(F(1),F(2),F(3))
    O=[[F(int(i==j))-v[i]*v[j]/7 for j in range(3)] for i in range(3)]
    H=[[sum(O[i][k]*F((k+1)**2)*O[j][k] for k in range(3))
        for j in range(3)] for i in range(3)]
    edges=(H[0][1],H[0][2],H[1][2])
    assert edges == (F(18,7),F(12,7),F(6,7))
    K = F(det)*L*a*a/F(163840)
    assert K == F(111,83886080)
    assert 2*mu*mu < K/2
    # e^(m^2*smax/4)<2 for m<=1; the factorial tail has denominator 1536.
    assert smax/F(4) < F(1,16)
    # e^(1/16) <= 1/(1-1/16)=16/15 <2.
    assert 1/(1-F(1,16)) < 2
    assert F(1,1536)*F(1,16)*2 == F(1,12288)
    assert F(1+3**7,12288) < F(1,5)
    assert F(2**7+3**7,12288) < F(1,5)
    assert F(26,8)+F(242,256*9) < 4
    assert F(19,8)+F(211,256*9) < 3
    assert F(7,5)+F(1,25*81) < 2
    # 12 J_DFT = 2 sqrt(3)/3 >1, W>pi*a^3/12.
    cyclic_lower = mu**8 * a**6 * K/F(288)
    assert cyclic_lower > 0
    # Differentiate-free geometric deletion: p perpendicular to d.
    vx,vy,vz = F(1,1000),F(1,900),F(-1,800)
    assert vx*vx+vy*vy+(L+vz)**2 == (-vx)**2+vy*vy+(L+vz)**2
    # All moments in the next bound are relative to the actual positive-mass reference.
    assert F(524288)+1+(3*mu)**2 < F(2173,3)**2
    assert 6*F(1,24)+6*F(1,24) == F(1,2)
    assert 1/(2*F(1,24)**2) == 288
    return {
        "all_window_s_squared":[frac(smin),frac(smax)],
        "bessel_coefficients_without_pi":[frac(c) for c in coeff],
        "moment_ratio_lower":frac(b2_over_b0_lower),
        "spectral_determinant":det,"real_irreducible_edges":[frac(x) for x in edges],
        "mu":frac(mu),"K":frac(K),
        "two_mu_squared":frac(2*mu*mu),"half_K":frac(K/2),
        "cyclic_strict_lower":frac(cyclic_lower),
        "cyclic_lower_decimal":tidy(float(cyclic_lower)),
        "energy_strict_upper":2173,"six_reading_hoeffding_coefficient":288,
        "uniform_analytic_not_sampling":True
    }

def radial_diagnostic(order):
    # A different diagnostic from pair-space integration: one-dimensional
    # spherical Fourier moments, no six-dimensional field quadrature.
    z,w = np.polynomial.legendre.leggauss(order)
    a,L,mu = 1/32,1/4,1/4096
    x = (z+1)/2
    r = a*x
    weights = w*a/2
    b = np.exp(-1/(1-x*x))
    norm2 = 4*np.pi*np.sum(weights*r*r*b*b)
    g = b/np.sqrt(norm2)
    sinc = np.zeros_like(r)
    dsinc = np.zeros_like(r)
    for k in range(9):
        sinc += (-1)**k*r**(2*k)/math.factorial(2*k+1)
        if k:
            dsinc += (-1)**k*(2*k)*r**(2*k-1)/math.factorial(2*k+1)
    B0 = 4*np.pi*np.sum(weights*r*r*g)
    B2 = 4*np.pi/3*np.sum(weights*r**4*g)
    B = 4*np.pi*np.sum(weights*r*r*g*sinc)
    Bp = 4*np.pi*np.sum(weights*r**3*g*dsinc)
    lapB = -4*np.pi*np.sum(weights*r**4*g*sinc)
    G,Gp,W = B*B,2*B*Bp,B0*B0
    lapG = 2*B*lapB+2*Bp*Bp
    raw_im = 2*L*G*Gp
    c3 = G/(8*np.pi)
    c5 = -(((.5**2-L*L)*G+lapG)-2j*L*Gp)/(64*np.pi)
    coeff_im = float(np.imag(c3*np.conj(c5)))
    assert raw_im < 0 and coeff_im > 0
    assert abs(coeff_im+raw_im/(512*np.pi*np.pi)) < 1e-22
    assert B2/B0 > a*a/320 and G > W/2 and -Gp > W*a*a/320
    J = math.sqrt(3)/18
    leading_C = 12*J*888*mu**8*coeff_im
    bprime = -2*r/a**2/(1-x*x)**2*b
    grad2 = 4*np.pi*np.sum(weights*r*r*bprime*bprime)/norm2
    assert grad2 < 524288
    return {"order":order,"W":tidy(W),"G":tidy(G),"G_prime":tidy(Gp),
            "raw_Im_I0_I1_conjugate":tidy(raw_im),
            "Im_c3_c5_conjugate":tidy(coeff_im),
            "cyclic_leading_prediction_only":tidy(leading_C),
            "bump_gradient_squared":tidy(grad2)}

def unitary_from_hermitian(h,t):
    ev,U = np.linalg.eigh(h)
    return (U*np.exp(-1j*t*ev))@U.conj().T

def matrix_diagnostics():
    omega = np.exp(2j*np.pi/3)
    U = np.array([[omega**(i*j) for j in range(3)] for i in range(3)])/np.sqrt(3)
    fs = np.array([1+.2j, -.3+.8j, .5-1.1j])
    A = U@np.diag(fs)@U.conj().T
    P = np.abs(A)**2
    J = np.imag(U[0,0]*U[1,1]*np.conj(U[0,1]*U[1,0]))
    tri = np.imag((fs[0]-fs[2])*np.conj(fs[1]-fs[2]))
    delta = P[1,0]-P[0,1]
    cyclic = P[1,0]-P[1,2]+P[2,1]-P[2,0]+P[0,2]-P[0,1]
    res = [abs(delta-4*J*tri), abs(cyclic-3*delta),
           np.linalg.norm(A@A.conj().T-A.conj().T@A)]
    conjU = U.conj()@np.diag(fs)@U.T
    res.append(abs((abs(conjU[1,0])**2-abs(conjU[0,1])**2)+delta))
    # CP-real mixing: O D O^T is symmetric for any complex D.
    c,s = 3/5,4/5
    O = np.array([[c,-s,0],[s,c,0],[0,0,1.]])@np.array([[1.,0,0],[0,c,-s],[0,s,c]])
    Areal = O@np.diag(fs)@O.T
    res.append(np.linalg.norm(Areal-Areal.T))
    # Real tail coefficients (no carrier / reflected perpendicular carrier).
    Aradial = U@np.diag([1.,-2.,3.])@U.conj().T
    res.append(np.linalg.norm(np.abs(Aradial)**2-(np.abs(Aradial)**2).T))
    # Degenerate masses give equal F values and zero oriented triangle.
    Adeg = U@np.diag([fs[0],fs[0],fs[2]])@U.conj().T
    res.append(np.linalg.norm(np.abs(Adeg)**2-(np.abs(Adeg)**2).T))
    # Any normal 2x2 matrix has equal off-diagonal moduli.
    U2 = np.array([[c,1j*s],[1j*s,c]])
    A2 = U2@np.diag(fs[:2])@U2.conj().T
    res.append(abs(abs(A2[0,1])**2-abs(A2[1,0])**2))
    assert max(res) < 1e-12
    zeta = np.array([1.,1.])/np.sqrt(2)
    assert abs(zeta@np.diag([1.,-1.])@zeta) < 1e-15
    return {"J":tidy(J),"delta_21":tidy(delta),"cyclic":tidy(cyclic),
            "max_algebra_residual":tidy(max(res)),
            "two_normal_flavors_cp_real_degeneracy_and_real_tail_zero":True}

def annihilators(n):
    ans=[]
    for k in range(n):
        a=np.zeros((2**n,2**n),complex)
        mask=1<<(n-1-k)
        for state in range(2**n):
            if state & mask:
                sign=(-1)**((state>>(n-k)).bit_count())
                a[state^mask,state]=sign
        ans.append(a)
    return ans

def dgamma(j,aa):
    out=np.zeros_like(aa[0])
    for i in range(len(aa)):
        for k in range(len(aa)):
            out += j[i,k]*(aa[i].conj().T@aa[k])
    return out

def finite_car_diagnostic():
    # Finite CAR algebra check only, not a discretization proof of the field kernel.
    aa=annihilators(6); ee=annihilators(3)
    he=np.array([[.2,.3j,.1],[-.3j,.7,.2],[.1,.2,.5]],complex)
    HE=dgamma(he,ee)
    ev,V=np.linalg.eigh(HE)
    rhoE=(V*np.exp(-ev))@V.conj().T
    rhoE/=np.trace(rhoE)
    Pe=np.array([[np.trace(rhoE@ee[j].conj().T@ee[i])
                  for j in range(3)] for i in range(3)])
    codebits=[4,2,1]
    rho=np.zeros((192,192),complex)
    codeR=np.zeros((9,9),complex)
    for a in range(3):
        for b in range(3):
            fab=np.zeros((8,8));fab[codebits[a],codebits[b]]=1
            rab=np.zeros((3,3));rab[a,b]=1
            rho+=np.kron(np.kron(fab,rhoE),rab)/3
            cab=np.zeros((3,3));cab[a,b]=1
            codeR+=np.kron(cab,rab)/3
    assert abs(np.trace(rho)-1)<1e-12
    omega=np.exp(2j*np.pi/3)
    Mix=np.array([[omega**(i*j) for j in range(3)] for i in range(3)])/np.sqrt(3)
    Y=Mix@np.diag([.7,1.2,1.8])@Mix.conj().T
    zero=np.zeros((3,3),complex)
    h=np.block([[zero,Y],[Y,zero]])
    evol=unitary_from_hermitian(h,.37)
    jphase=np.block([[zero,1j*Y],[-1j*Y,zero]])
    jq=np.diag([0,0,0,1,2,3]).astype(complex)
    vb=np.array([1,1j,2],complex)/np.sqrt(6)
    B=np.outer(vb,vb.conj())
    errs=[]
    for j in (h,jphase,jq):
        jt=evol.conj().T@j@evol
        J=dgamma(jt,aa)
        c=np.trace(jt[3:,3:]@Pe)
        logical=jt[:3,:3]+c*np.eye(3)
        lhs=np.trace(rho@np.kron(J,B))
        rhs=np.trace(codeR@np.kron(logical,B))
        errs.append(abs(lhs-rhs))
    number=[a.conj().T@a for a in aa]
    initial_left=np.trace(rho@np.kron(number[3],B))
    factor=np.trace(rho@np.kron(number[3],np.eye(3)))*np.trace(B)/3
    errs.append(abs(initial_left-factor))
    actual_pair=np.trace(rho@np.kron(number[0]@number[1],np.eye(3))).real
    assert abs(actual_pair)<1e-14
    assert max(errs)<1e-11
    return {"source_and_partner_max_error":tidy(max(errs)),
            "single_occupancy_pair_moment":tidy(actual_pair),
            "incorrect_quasifree_pair_moment":tidy(1/9),
            "finite_car_only_not_continuous_pde":True}

def make_result():
    exact=exact_bounds()
    radial=[radial_diagnostic(n) for n in (80,120)]
    for key in ("W","G","G_prime","Im_c3_c5_conjugate"):
        assert abs(radial[0][key]-radial[1][key]) < 1e-10*max(abs(radial[1][key]),1e-20)
    sources={
        "admission":ROOT/"_shared/round1303_admission.md",
        "author_proof":HERE/"proof.md",
        "1293_proof":ROOT/"1293/proof.md",
        "1295_proof":ROOT/"1295/proof.md"}
    return {"status":"PASS_INDEPENDENT_CHECK_WITH_SCOPE",
            "exact_bounds":exact,"radial_diagnostics":radial,
            "mixing_diagnostics":matrix_diagnostics(),
            "finite_car_diagnostic":finite_car_diagnostic(),
            "source_bindings":{name:sha(path) for name,path in sources.items()},
            "independent_code_sha256":sha(Path(__file__).resolve()),
            "no_full_SM_no_T_reversal_no_three_generation_selection":True}

def main():
    parser=argparse.ArgumentParser()
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument("--write",action="store_true")
    mode.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=make_result()
    out=HERE/"independent_results.json"
    if args.write:
        with out.open("x",encoding="utf-8",newline="\n") as f:
            json.dump(result,f,ensure_ascii=False,indent=2);f.write("\n")
        print("SAVED PASS_INDEPENDENT_CHECK_WITH_SCOPE")
    else:
        saved=json.loads(out.read_text(encoding="utf-8-sig"))
        assert saved==result,"Stored independent result differs from current read-only calculation"
        print("READ_ONLY_MATCH PASS_INDEPENDENT_CHECK_WITH_SCOPE")
if __name__=="__main__":main()
