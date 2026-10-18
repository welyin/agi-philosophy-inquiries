"""Independent checks for 1305. No author code is imported.
Default/--check read only. --write creates results exclusively.
Finite CAR and quadrature are diagnostics; analytic existence and PDE/source
bounds are reviewed in independent_review.md, not inferred from these tests.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import json
from math import factorial, pi, sqrt
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ARCHIVE = HERE.parent

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def annihilators(n):
    dim = 1 << n
    out = []
    for j in range(n):
        a = np.zeros((dim, dim), complex)
        for state in range(dim):
            if (state >> j) & 1:
                sign = (-1)**((state & ((1 << j)-1)).bit_count())
                a[state ^ (1 << j), state] = sign
        out.append(a)
    return out

def car_source():
    # Three coded modes plus three complement modes.  The lower three bits
    # are the code, so the matrix order is complement x code x partner.
    aa = annihilators(6)
    h = np.kron(np.array([[0, 1], [1, 0]], complex), np.eye(3))
    P = (np.eye(6)-h)/2
    rng = np.random.default_rng(130507)
    psi = rng.normal(size=9)+1j*rng.normal(size=9)
    psi /= np.linalg.norm(psi)
    rho_qr = np.outer(psi, psi.conj())
    W = np.zeros((8, 3), complex)
    for i in range(3):
        W[1 << i, i] = 1
    wr = np.kron(W, np.eye(3))
    state = np.kron(np.eye(8)/8, wr@rho_qr@wr.conj().T)
    rho_q = np.trace(rho_qr.reshape(3,3,3,3),axis1=1,axis2=3)
    rho_r = np.trace(rho_qr.reshape(3,3,3,3),axis1=0,axis2=2)
    rho_field = np.trace(state.reshape(64,3,64,3),axis1=1,axis2=3)
    state_r = np.trace(state.reshape(64,3,64,3),axis1=0,axis2=2)
    C = np.empty((6,6),complex)
    for i in range(6):
        for j in range(6):
            C[i,j] = np.trace(rho_field@aa[j].conj().T@aa[i])
    expected = np.zeros((6,6),complex)
    expected[:3,:3] = rho_q
    expected[3:,3:] = np.eye(3)/2
    trace_flavor = np.trace((C-P).reshape(2,3,2,3),axis1=1,axis2=3)
    expected_trace = np.array([[-.5,1.5],[1.5,0]],complex)
    errors = [
        np.max(np.abs(C-expected)),
        np.max(np.abs(state_r-rho_r)),
        np.max(np.abs(trace_flavor-expected_trace))
    ]
    # Full field source and arbitrary partner observable, without replacing
    # the non-Gaussian one-particle code by a quasi-free state.
    source_errors = []
    for unused in range(4):
        D = rng.normal(size=(6,6))+1j*rng.normal(size=(6,6))
        D = (D+D.conj().T)/2
        B = rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
        B = (B+B.conj().T)/2
        source = sum(D[i,j]*(aa[i].conj().T@aa[j])
                     for i in range(6) for j in range(6))
        left = np.trace(state@np.kron(source,B))
        left -= np.trace(D@P)*np.trace(rho_r@B)
        sq = D[:3,:3] + (np.trace(D[3:,3:])/2-np.trace(D@P))*np.eye(3)
        right = np.trace(rho_qr@np.kron(sq,B))
        source_errors.append(float(abs(left-right)))
    assert max(errors+source_errors)<3e-12
    return {
      "dimension_field":64,"dimension_partner":3,
      "trace_flavor_delta":[[-0.5,1.5],[1.5,0.0]],
      "max_covariance_partner_trace_error":float(max(errors)),
      "max_source_partner_error":max(source_errors),
      "scope":"Finite CAR diagnostic of the analytic three-mode source identity."
    }

def radial_moments(nodes, max_l):
    # Only one radial integration is used.  Convolution moments follow from
    # the Taylor coefficients of the squared radial characteristic function.
    a = 1/32
    x,w = np.polynomial.legendre.leggauss(nodes)
    r = a*(x+1)/2
    weights = a*w/2*4*pi*r*r
    b = np.exp(-1/(1-(r/a)**2))
    norm = sqrt(float(np.dot(weights,b*b)))
    R = [float(np.dot(weights,b*r**(2*j)))/norm for j in range(max_l+2)]
    M = []
    for k in range(max_l+2):
        value = sum(R[j]*R[k-j]/(factorial(2*j+1)*factorial(2*(k-j)+1))
                    for j in range(k+1))
        M.append(factorial(2*k+1)*value)
    return M, R[0]**2

def bessel_derivatives(t, moments, max_l):
    # D=t^2-L^2-r^2-2 L z.  Radial isotropy gives
    # integral r^(2u) z^(2v) rho = M[ u+v ]/(2v+1).
    L=1/4
    d0=t*t-L*L
    ap=app=bp=bpp=0.0
    for ell in range(max_l+1):
        I=J=0.0
        for u in range(ell+1):
            for v in range(ell-u+1):
                q=ell-u-v
                cf = factorial(ell)/(factorial(u)*factorial(v)*factorial(q))
                cf *= (-1)**u*(-2*L)**v*d0**q
                if v%2==0:
                    I += cf*moments[u+v//2]/(v+1)
                else:
                    J += cf*moments[u+(v+1)//2]/(v+2)
        cf = (-1)**ell/(2**(2*ell+3)*pi*factorial(ell)*factorial(ell+1))
        aval = cf*I
        bval = -cf*J  # derivative of exp(-i p z) at p=0
        power=2*ell+3
        ap += power*aval
        app += power*(power-1)*aval
        bp += power*bval
        bpp += power*(power-1)*bval
    return ap,app,bp,bpp,bp*app-ap*bpp

def small_s_diagnostic(nodes, terms):
    M,W0=radial_moments(nodes,terms)
    x,w=np.polynomial.legendre.leggauss(40)
    # One fixed smooth positive time weight; all normalization is independent
    # of s.  It vanishes to all orders at the ends of the admitted interval.
    times=(7/16+1/2)/2+(1/2-7/16)*x/2
    tw=w*np.exp(-1/(1-x*x))
    tw/=np.sum(tw)
    values=np.array([bessel_derivatives(float(t),M,terms) for t in times])
    avg=np.sum(values*tw[:,None],axis=0)
    mu=2**-12
    J=1/(6*sqrt(3))
    return {
      "radial_nodes":nodes,"bessel_last_index":terms,
      "W0":W0,
      "A_prime_at_1":float(avg[0]),"A_second_at_1":float(avg[1]),
      "B_prime_at_1":float(avg[2]),"B_second_at_1":float(avg[3]),
      "weighted_wronskian":float(avg[4]),
      "C_s4_coefficient":float(12*J*mu**3*avg[4]),
      "method":"Radial Fourier-square moments; no nearly degenerate amplitude subtraction.",
      "scope":"Numerical derivative diagnostic, not the nonzero/common-window proof."
    }

def exact_path():
    mu=F(1,4096)
    slopes=[j*mu-1 for j in (1,2,3)]
    x1=slopes[0]-slopes[2]
    y1=slopes[1]-slopes[2]
    x2=(slopes[0]**2-slopes[2]**2)/2
    y2=(slopes[1]**2-slopes[2]**2)/2
    determinant=x1*y2-x2*y1
    assert determinant==mu**3
    # Source before taking any center or null contraction.
    n=3
    assert F(1)-F(n,2)==F(-1,2)
    K=F(111,83886080)
    assert 2*mu**2<K/2  # inherited endpoint-tail margin
    cstar=K*F(1,32)**6*mu**8/288
    assert cstar>0
    return {
      "mass_second_order_determinant":str(determinant),
      "expected_mu_cubed":str(mu**3),
      "three_mode_ff_coefficient":"-1/2",
      "three_mode_cross_coefficient":"3/2",
      "endpoint_cstar_exact":str(cstar),
      "endpoint_cstar_decimal":float(cstar),
      "parameter_order":[
        "Fix the three-mode s=0 source, H0, Lambda0 and thermal/macroscopic positive margins.",
        "Choose positive s inside source/thermal continuity margins with C0(s) nonzero.",
        "Choose positive kappa below the phase transport bound and all fixed geometric bounds.",
        "Set a separate positive finite phase readout and sampling error budget."
      ],
      "scope":"Order and algebra checks; no numerical value is assigned to the PDE constants."
    }

def cyclic_identity():
    # Generic normal three-flavor amplitudes test the six-probability identity,
    # independently of small-s continuum amplitudes.
    z=np.exp(2j*pi/3)
    U=np.array([[z**(a*j)/sqrt(3) for j in range(3)] for a in range(3)])
    Fm=np.array([.17+.08j,-.11+.09j,.06-.05j])
    A=U@np.diag(Fm)@U.conj().T
    base=np.array([.2,.3,.4])
    P=base[:,None]+abs(A)**2
    C=P[1,0]-P[1,2]+P[2,1]-P[2,0]+P[0,2]-P[0,1]
    J=np.imag(U[0,0]*U[1,1]*U[0,1].conjugate()*U[1,0].conjugate())
    expected=12*J*np.imag((Fm[0]-Fm[2])*(Fm[1]-Fm[2]).conjugate())
    assert abs(C-expected)<2e-15
    return {"cyclic_value":float(C),"formula_value":float(expected),
            "baseline":"Synthetic input-independent row offsets only."}

def compute():
    d1=small_s_diagnostic(64,12)
    d2=small_s_diagnostic(96,16)
    rel=abs(d1["C_s4_coefficient"]/d2["C_s4_coefficient"]-1)
    assert rel<2e-8
    assert d2["weighted_wronskian"]>0
    paths=[
      "_shared/round1305_admission.md",
      "1298/proof.md","1300/proof.md","1300/acceptance.json",
      "1303/proof.md","1303/acceptance.json",
      "_shared/nondegenerate_joint_corollary_after1300.md",
      "1305/proof.md","1305/check.py","1305/results.json","1305/NEXT.md",
      "research_note_1305.md"
    ]
    return {
      "status":"PASS",
      "independent_code_sha256":sha(Path(__file__)),
      "bindings":{p:sha(ARCHIVE/p) for p in paths},
      "exact":exact_path(),
      "finite_car":car_source(),
      "cyclic":cyclic_identity(),
      "small_s_diagnostics":[d1,d2],
      "diagnostic_relative_difference":rel,
      "limits":[
        "Finite matrices do not prove the continuous null source or PDE estimates.",
        "The common-window proof uses real analyticity and the inherited positive endpoint.",
        "No all-state norm continuity of different-mass vacua is claimed.",
        "Thermal/macroscopic and phase-record errors are separate.",
        "No autonomous device, full CC parent, nonzero four-fermion or full SM is certified."
      ]
    }

def main():
    ap=argparse.ArgumentParser()
    group=ap.add_mutually_exclusive_group()
    group.add_argument("--check",action="store_true")
    group.add_argument("--write",action="store_true")
    args=ap.parse_args()
    result=compute()
    target=HERE/"independent_results.json"
    if args.write:
        with target.open("x",encoding="utf-8",newline="\n") as f:
            json.dump(result,f,ensure_ascii=False,indent=2,sort_keys=True)
            f.write("\n")
    else:
        saved=json.loads(target.read_text(encoding="utf-8"))
        assert result==saved,"Saved diagnostic or bindings changed"
    print(json.dumps({"status":"PASS","mode":"write" if args.write else "read-only",
       "C_s4_diagnostic":result["small_s_diagnostics"][-1]["C_s4_coefficient"],
       "source_error":result["finite_car"]["max_source_partner_error"]},
       ensure_ascii=False))
if __name__=="__main__":
    main()
