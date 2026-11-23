"""Independent coherent neutral-kaon scalar mixing at mu=163 GeV."""
import argparse
import cmath
import itertools
import json
import math
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE/path).read_bytes())


def couplings(conjugate=False):
    fit=read("../1489/results.json")["reconstruction"]
    a,b,c,e,f=fit["normalized_real_edges_A_B_C_E_F"]
    r,d=fit["two_real_coefficients"]
    d*=cmath.exp(1j*math.radians(fit["eta_example_deg"]))
    if conjugate: d=d.conjugate()
    Y=np.array([[0,a,0],[b,c,r+d],[0,e,f]])
    m2,V=np.linalg.eigh(Y@Y.conj().T)
    R=Y.conj().T@V/np.sqrt(m2)
    np.testing.assert_allclose(V.conj().T@Y@R,np.diag(np.sqrt(m2)),atol=1e-11)
    source=V[0,1].conjugate()*V[0,0]
    unit=source/abs(source)
    # Coefficients of the barred-strange down operators; tree reference is K_us^* K_ud.
    def shape(tr,tl):
        tmp=np.outer(tr[1,:].conj(),tl[2,:])
        return tmp[1,0]/unit,tmp[0,1].conj()/unit
    gR,gL=shape(V,R)
    yb=read("../1386/inputs.json")["yukawa"]["down"][2][0]
    alpha0=.1181
    alpha=alpha0/(1+(23/3)*alpha0*math.log(163/91.1876)/(2*math.pi))
    norm=(yb/math.sqrt(2))*(alpha/alpha0)**(12/23)
    # Rephase physical mass columns and recompute the unchanged normalized shapes.
    re=np.diag(np.exp(1j*np.array([.31,.8,-.74])))
    V2,R2=V@re,R@re
    unit=(V2[0,1].conjugate()*V2[0,0]);unit/=abs(unit)
    np.testing.assert_allclose(shape(V2,R2),[gR,gL],atol=1e-12)
    return r,d,norm*gR,norm*gL,alpha


def match(alpha):
    q=alpha/(4*math.pi)
    L=2*math.log(125/163)
    # Table 2, Buras et al. 1303.3723. Matrix elements already divided by 2 m_K.
    return np.array([-.064*(1+q*(-3*L+4.5)), -.107*q*(-L/12+1/8),
                     -.11*(-1.5*q), .18*(1-q)])


def contributions(w,gR,gL,hadronic):
    R,L=w*gR,w.conjugate()*gL
    # Conjugate the (sbar d)^2 amplitude to report M12 in the paper's eqs. (51),(52).
    factors=np.array([.5*(R**2+L**2)]*2+[R*L]*2)
    return -(hadronic*factors).conjugate()/125**2


def cart(z):
    return [float(z.real),float(z.imag)]


def solve():
    r,d,R,L,alpha=couplings()
    H=match(alpha)
    yb=read("../1386/inputs.json")["yukawa"]["down"][2]
    normal=[(.95*(1-2*yb[1]/yb[0]))**2,(1.05*(1+2*yb[1]/yb[0]))**2]
    masses=[3.484e-15, .006e-15]
    tauS=[8.954e-11, .004e-11]
    tauL=[5.116e-8, .021e-8]
    hbar=6.582119569e-25
    # CPT two-state trace and positive width matrix imply 2 |M12| <= sqrt(dm²+GammaS GammaL).
    twice_radius=math.hypot(masses[0]+2*masses[1],hbar/math.sqrt(
        (tauS[0]-2*tauS[1])*(tauL[0]-2*tauL[1])))
    # Same physical hadronic responses, common mass and transfer in both candidates.
    terms=[contributions(2*w,R,L,H) for w in (r,d)]
    _,dc,Rc,Lc,ac=couplings(True)
    for w0,z in zip((r,dc),terms):
        np.testing.assert_allclose(contributions(2*w0,Rc,Lc,match(ac)),z.conj(),rtol=1e-10,atol=1e-27)
    weights=np.array(list(itertools.product((.8,1.2),repeat=4)))
    corners=normal[0]*((terms[1]-terms[0])@weights.T)
    # The linear image is a convex polygon. Project on the central difference as a certified bound.
    axis=sum(terms[1]-terms[0]);axis/=abs(axis)
    lower=min((corners/axis).real)
    assert lower>twice_radius
    def row(z):
        zs=np.r_[normal[0]*(weights@z),normal[1]*(weights@z)]
        return {"M12_GeV":cart(sum(z)),"Re_range":[min(zs.real),max(zs.real)],
                "Im_range":[min(zs.imag),max(zs.imag)],"terms":list(map(cart,z))}
    # For all real slopes n_R+n_D=2, w=n_R r+n_D d sweeps one straight line.
    wmin=2*abs((r*d.conjugate()).imag)/abs(d-r)
    lr=.8*abs(sum(H[2:]))*abs(R*L)
    same=.5*1.2*abs(sum(H[:2]))*(abs(R)**2+abs(L)**2)
    lower_all=normal[0]*wmin**2*(lr-same)/125**2
    assert lower_all>twice_radius/2
    assert row(terms[0])["Im_range"][0]>0>row(terms[1])["Im_range"][1]
    return {"alpha163":alpha,"gR_for_unit_w":cart(R),"gL_for_unit_w_conjugate":cart(L),
            "matrix_element_products":H.tolist(),"real_path":row(terms[0]),"phase_path":row(terms[1]),
            "difference_center_GeV":cart(sum(terms[1]-terms[0])),
            "shared_box_projection_lower_GeV":lower,"two_radii_upper_GeV":twice_radius,
            "minimum_model_specific_baseline_difference_GeV":lower-twice_radius,
            "min_abs_w_all_real_affine_slopes":wmin,"LR_minus_SLL_over_w2_GeV":(lr-same)/125**2,
            "all_affine_slopes_abs_M12_lower_GeV":lower_all}


def stable(value):
    if isinstance(value,dict): return {k:stable(v) for k,v in value.items()}
    if isinstance(value,list): return [stable(x) for x in value]
    return float(f"{value:.12g}")


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    cli.add_argument("--check",action="store_true")
    args=cli.parse_args()
    result=stable(solve())
    path=HERE/"independent_results.json"
    if args.check:
        assert result==json.loads(path.read_bytes())
        print("1649 independent PASS: matched kaon Wilson coefficients, rephased scalar amplitudes, common-baseline distance")
    else:
        with path.open("x",encoding="utf8") as file:
            file.write(json.dumps(result,indent=2)+"\n")

