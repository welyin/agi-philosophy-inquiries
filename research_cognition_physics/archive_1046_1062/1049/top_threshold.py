"""1049: one fixed top threshold subtraction, physical scattering-bin interference.

The benchmark is a fixed-order, flux-averaged cross-section coefficient, not an
all-order finite-time instrument or a fit to measured Standard Model parameters.
The rest of the SM and its CKM matrix are common; the neutral slice does not set
CKM to the identity or compute the complete weak/Wess-Zumino threshold.
"""
from pathlib import Path
from fractions import Fraction
import argparse
import hashlib
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULT = HERE / "top_threshold_results.json"
HISTORY = [
    "archive_585_628/research_note_627.md",
    "archive_585_628/research_note_628.md",
    "archive_629_652/research_note_629.md",
    "archive_702_741/research_note_735.md",
    "archive_764_796/research_note_789.md",
    "archive_956_989/981/drafts/common_parent_contract_v1.md",
    "archive_990_1008/993/common_candidate_v1.md",
]

ALPHA = 1 / 137
E2 = 4 * math.pi * ALPHA
SW2 = .23
VEV = math.sqrt(2.)
GZ2 = E2 / (4 * SW2 * (1-SW2))
MZ2 = GZ2 * VEV**2
MH2 = .72**2
ME, MMU = 3e-6, 6e-4
MT0 = 1.
KAPPA = 2 * ALPHA * 3 * (2/3)**2 / math.pi
SLO, SHI = .004, .009
CLO, CHI = -.6, .4
RLO, RHI = .9, 1.1
ORDER = 2

I2 = np.eye(2, dtype=complex)
Z2 = np.zeros((2,2), dtype=complex)
PAULI = [np.array([[0,1],[1,0]],complex),
         np.array([[0,-1j],[1j,0]],complex), np.diag([1.,-1.]).astype(complex)]
GAMMA = [np.block([[I2,Z2],[Z2,-I2]])]
GAMMA += [np.block([[Z2,s],[-s,Z2]]) for s in PAULI]
G5 = 1j * GAMMA[0] @ GAMMA[1] @ GAMMA[2] @ GAMMA[3]
G0 = GAMMA[0]
METRIC = np.diag([1.,-1.,-1.,-1.])
VF, AF = -.5 + 2*SW2, -.5
ZVERT = [g @ (VF*np.eye(4)-AF*G5) for g in GAMMA]


def dot(a,b):
    return a[0]*b[0]-np.dot(a[1:],b[1:])


def p2(s):
    return ((s-ME**2-MMU**2)**2-4*ME**2*MMU**2)/(4*s)


def spinors(momentum,mass):
    energy=math.sqrt(mass**2+float(np.dot(momentum,momentum)))
    sigma=sum((p*s for p,s in zip(momentum,PAULI)),start=np.zeros((2,2),complex))
    return np.vstack([math.sqrt(energy+mass)*I2, sigma/math.sqrt(energy+mass)])


def bilinear(out, incoming, vertices):
    return np.array([out.conj().T@G0@g@incoming for g in vertices])


def tree_data(s,costheta):
    p=math.sqrt(p2(s))
    sin=math.sqrt(1-costheta**2)
    p1=np.array([0.,0.,p]); p2v=-p1
    p3=np.array([p*sin,0.,p*costheta]); p4=-p3
    us=[spinors(p1,ME),spinors(p2v,MMU),spinors(p3,ME),spinors(p4,MMU)]
    je=bilinear(us[2],us[0],GAMMA); jm=bilinear(us[3],us[1],GAMMA)
    ze=bilinear(us[2],us[0],ZVERT); zm=bilinear(us[3],us[1],ZVERT)
    he=us[2].conj().T@G0@us[0]; hm=us[3].conj().T@G0@us[1]
    q=np.concatenate([[0.],p3-p1]); t=dot(q,q)
    photon=[]; zboson=[]; higgs=[]
    for a in range(2):
        for b in range(2):
            for f in range(2):
                for g in range(2):
                    photon.append(E2/t*dot(je[:,f,a],jm[:,g,b]))
                    zboson.append(GZ2/(t-MZ2)*(dot(ze[:,f,a],zm[:,g,b])
                                   -dot(q,ze[:,f,a])*dot(q,zm[:,g,b])/MZ2))
                    higgs.append(-(ME*MMU/VEV**2)/(t-MH2)*he[f,a]*hm[g,b])
    photon=np.array(photon); extra=np.array(zboson)+np.array(higgs)
    baseline=photon+extra
    gamma_sq=float(np.vdot(photon,photon).real/4)
    sigma=ME**2+MMU**2
    u=2*sigma-s-t
    analytic=2*E2**2/t**2*((s-sigma)**2+(u-sigma)**2+2*t*sigma)
    assert math.isclose(gamma_sq,analytic,rel_tol=2e-12,abs_tol=1e-13)
    denom=32*math.pi*s
    W=float(2*np.vdot(baseline,photon).real/4/denom)
    return {
        "t":float(t),"W":W,"qed_weight":2*gamma_sq/denom,
        "born_sm":float(np.vdot(baseline,baseline).real/4/denom),
        "extra_amplitude_ratio":float(np.linalg.norm(extra)/np.linalg.norm(photon)),
    }


def beta_integral(power):
    return Fraction(math.factorial(power)**2,math.factorial(2*power+1))


def kernels(t,r,degree=ORDER,nquad=80):
    nodes,weights=np.polynomial.legendre.leggauss(nquad)
    x=(nodes+1)/2; weights=weights/2; u=x*(1-x)
    q2=-np.asarray(t,dtype=float)
    w=q2[...,None]*u/(MT0*r)**2
    P=KAPPA/3*math.log(r)+KAPPA*np.sum(weights*u*np.log1p(w),axis=-1)
    S=2*KAPPA*np.sum(weights*u/(1+w),axis=-1)
    Pn=np.full_like(q2,KAPPA/3*math.log(r))
    Sn=np.full_like(q2,KAPPA/3)
    for n in range(1,degree+1):
        power=(q2/(MT0*r)**2)**n*float(beta_integral(n+1))
        Pn += KAPPA*(-1)**(n+1)*power/n
        Sn += 2*KAPPA*(-1)**n*power
    return P,S,Pn,Sn


def analytic_box():
    sigma=ME**2+MMU**2
    qmin=2*p2(SLO)*(1-CHI)
    qmax=SHI*(1-CLO)/2 # conservative bound p^2 <= s/4
    numlow=(SLO-sigma)**2-2*qmax*sigma
    assert numlow>0
    gmin=math.sqrt(2*E2**2*numlow/qmax**2)
    gmaxsq=4*E2**2*SHI**2/qmin**2
    b=abs(VF)+abs(AF)
    z_each=GZ2/MZ2*(16+12*qmax/MZ2)*SHI/4*b*b
    h_each=ME*MMU/VEV**2*SHI/MH2
    extra_norm=2*(z_each+h_each)
    eta=extra_norm/gmin
    assert eta<.16
    Wlo=2*(1-eta)*gmin**2/(32*math.pi*SHI)*(CHI-CLO)
    Whi=2*(1+eta)*gmaxsq/(32*math.pi*SLO)*(CHI-CLO)
    zstar=qmax/(4*(MT0*RLO)**2)
    remainder=KAPPA*(4*zstar)**(ORDER+1)*float(beta_integral(ORDER+2))/(ORDER+1)
    source_remainder=2*(ORDER+1)*remainder
    Plo=KAPPA*qmin/(30*MT0**2*(1+qmax/(4*MT0**2)))
    return {"q_squared_min":qmin,"q_squared_max_upper":qmax,
            "weak_higgs_amplitude_ratio_upper":eta,
            "bin_W_integral_lower":Wlo,"bin_W_integral_upper":Whi,
            "z_star":zstar,"P_uniform_error_upper":remainder,
            "S_uniform_error_upper":source_remainder,
            "bin_cross_section_error_upper":Whi*remainder,
            "bin_source_error_upper":Whi*source_remainder,
            "nonzero_top_bin_lower_at_r1":Wlo*Plo,
            "lost_source_bin_lower_at_r1":KAPPA/3*Wlo,
            "lost_source_kernel_exact":KAPPA/3}


def box_grid(ns,nc):
    xs,ws=np.polynomial.legendre.leggauss(ns)
    xc,wc=np.polynomial.legendre.leggauss(nc)
    ss=SLO+(xs+1)*(SHI-SLO)/2
    cc=CLO+(xc+1)*(CHI-CLO)/2
    rows=[]; quad=[]
    for s,a in zip(ss,ws/2): # normalized uniform luminosity profile in s
        for c,b in zip(cc,wc*(CHI-CLO)/2):
            rows.append(tree_data(float(s),float(c)));quad.append(float(a*b))
    return rows,np.array(quad)


def integrate(ns,nc):
    rows,qw=box_grid(ns,nc)
    t=np.array([r["t"] for r in rows]); W=np.array([r["W"] for r in rows])
    weight=qw*W
    assert W.min()>0
    output={"W_integral":float(weight.sum()),
            "born_sm_bin":float(sum(w*r["born_sm"] for w,r in zip(qw,rows))),
            "qed_interference_weight_integral":float(sum(w*r["qed_weight"] for w,r in zip(qw,rows))),
            "max_extra_amplitude_ratio_sample":max(r["extra_amplitude_ratio"] for r in rows)}
    values=[]
    for r in (RLO,1.,RHI):
        P,S,Pn,Sn=kernels(t,r)
        lost=KAPPA/3
        values.append({"r":r,"top_bin_coefficient":float(weight@P),
                       "matched_top_bin_coefficient":float(weight@Pn),
                       "radial_top_kernel_bin":float(weight@S),
                       "matched_radial_top_kernel_bin":float(weight@Sn),
                       "per_background_subtracted_bin":float(weight@(P-lost*math.log(r))),
                       "per_background_subtracted_source_bin":float(weight@(S-lost)),
                       "P_error_max":float(np.max(abs(P-Pn))),
                       "S_error_max":float(np.max(abs(S-Sn)))})
    output["kernel_bins_on_fixed_reference_external_basis"]=values
    # The r=1 derivative difference is independent of all smooth common
    # variations of external states/flux/geometry: DeltaP(1,t) identically zero.
    output["total_source_difference_at_r1"]=float(weight.sum()*KAPPA/3)
    return output


def run():
    for mu in range(4):
        for nu in range(4):
            assert np.linalg.norm(GAMMA[mu]@GAMMA[nu]+GAMMA[nu]@GAMMA[mu]
                                  -2*METRIC[mu,nu]*np.eye(4))<1e-14
    bounds=analytic_box()
    fine=integrate(24,32);coarse=integrate(16,24)
    quad_difference=abs(fine["W_integral"]-coarse["W_integral"])
    assert quad_difference<2e-11
    assert bounds["bin_W_integral_lower"]<fine["W_integral"]<bounds["bin_W_integral_upper"]
    assert fine["max_extra_amplitude_ratio_sample"]<bounds["weak_higgs_amplitude_ratio_upper"]
    for row in fine["kernel_bins_on_fixed_reference_external_basis"]:
        assert row["P_error_max"]<=bounds["P_uniform_error_upper"]+1e-16
        assert row["S_error_max"]<=bounds["S_uniform_error_upper"]+1e-16
        assert abs(row["top_bin_coefficient"]-row["matched_top_bin_coefficient"])<=bounds["bin_cross_section_error_upper"]+1e-16
        assert abs(row["radial_top_kernel_bin"]-row["matched_radial_top_kernel_bin"])<=bounds["bin_source_error_upper"]+1e-16
    r1=fine["kernel_bins_on_fixed_reference_external_basis"][1]
    assert r1["top_bin_coefficient"]>bounds["nonzero_top_bin_lower_at_r1"]>20*bounds["bin_cross_section_error_upper"]
    assert r1["top_bin_coefficient"]==r1["per_background_subtracted_bin"]
    assert fine["total_source_difference_at_r1"]>bounds["lost_source_bin_lower_at_r1"]>1e4*bounds["bin_source_error_upper"]
    # Check mass differentiation of the SAME fixed-subtraction scalar kernel;
    # finite differences are diagnostics; the uniform derivative is analytic.
    t=-.003; step=1e-4
    f=lambda log_r:float(kernels(t,math.exp(log_r))[0])
    derivative=(f(-2*step)-8*f(-step)+8*f(step)-f(2*step))/(12*step)
    source=float(kernels(t,1.)[1])
    assert abs(derivative-source)<2e-14
    return {"round":1049,"new_scientific_groups":1,"new_cognitive_axioms":0,
            "all_scientific_checks_passed":True,
            "parameters":{"alpha":ALPHA,"m_top_reference":MT0,"vev_reference":VEV,
                          "m_e":ME,"m_mu":MMU,"sin_theta_w_squared":SW2,
                          "m_Z_squared":MZ2,"m_H_squared":MH2,"kappa":KAPPA,
                          "s_interval":[SLO,SHI],"cos_theta_interval":[CLO,CHI],
                          "rho_ratio_interval":[RLO,RHI],"matching_order":ORDER,
                          "empirical_parameter_fit":False},
            "analytic_uniform_bounds":bounds,"finite_flux_bin":fine,
            "quadrature_W_difference":quad_difference,
            "source_derivative_check":{"t":t,"five_point_derivative":derivative,
                                       "same_integral_source":source,
                                       "absolute_difference":abs(derivative-source)},
            "scope":{"baseline":"Complete gamma,Z,H e-mu tree baseline for the displayed one-loop insertion; other one-loop terms common",
                     "observable":"Finite-resolution luminosity-averaged cross-section coefficient, not an arbitrary coherent wavepacket probability",
                     "source":"Top neutral two-point radial zero-Higgs-momentum insertion; at r=1 counterterm defect survives arbitrary common smooth external variations",
                     "full_nonzero_momentum_higgs_vertex":False,
                     "complete_weak_WZ_threshold_computed":False,
                     "CKM_set_to_identity":False,
                     "finite_time_all_order_CP_instrument":False,
                     "muon_ir_higher_order_certified":False,
                     "full_SM_prediction_or_roadmap_completed":False},
            "historical_sha256":{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in HISTORY}}


def compare(a,b,path="root"):
    if isinstance(a,dict):
        assert a.keys()==b.keys(),path
        for k in a:compare(a[k],b[k],path+"."+k)
    elif isinstance(a,list):
        assert len(a)==len(b),path
        for i,(x,y) in enumerate(zip(a,b)):compare(x,y,path+"."+str(i))
    elif isinstance(a,float):
        assert math.isclose(a,b,rel_tol=3e-10,abs_tol=2e-13),(path,a,b)
    else:assert a==b,(path,a,b)


def main():
    parser=argparse.ArgumentParser();parser.add_argument("--write",action="store_true")
    args=parser.parse_args();result=run()
    if args.write:
        with RESULT.open("x",encoding="utf-8") as f:
            json.dump(result,f,ensure_ascii=False,indent=2);f.write("\n")
    else:compare(result,json.loads(RESULT.read_text(encoding="utf-8")))
    print(json.dumps({"round":1049,"all_scientific_checks_passed":True,
                      "mode":"exclusive_write" if args.write else "read_only_compare",
                      "analytic_uniform_bounds":result["analytic_uniform_bounds"],
                      "r1":result["finite_flux_bin"]["kernel_bins_on_fixed_reference_external_basis"][1],
                      "full_baseline_W_integral":result["finite_flux_bin"]["W_integral"]},ensure_ascii=False))


if __name__=="__main__":main()
