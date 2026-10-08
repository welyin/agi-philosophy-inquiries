"""1016: scalar three/four-point vertices from a common local gauge identity.

Calibrates the Noether Euler coefficient E-div(D)+box(Q), not merely a
pointwise density variation. S below is the canonical S_eff=S_raw-T_raw
after quotienting the allowed derivative improvements by total divergences.
This is a specified bilinear scalar vertex class, not a full EFT classification.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent
BASE=HERE.parents[1]
OUT=HERE/"scalar_gauge_vertex_selection_results.json"
ETA=np.array([-1.,1.,1.,1.])


def norm(a):
    return float(np.linalg.norm(np.asarray(a).reshape(-1)))


def real(z):
    assert abs(complex(z).imag)<1e-9*max(1.,abs(z)),z
    return float(complex(z).real)


def su2():
    sigma=np.array([[[0,1],[1,0]],[[0,-1j],[1j,0]],[[1,0],[0,-1]]],complex)
    f=np.zeros((3,3,3))
    for a,b,c in [(0,1,2),(1,2,0),(2,0,1)]:
        f[a,b,c]=1
        f[a,c,b]=-1
    return sigma/2,f


def bracket(f,a,b):
    return np.einsum("abc,...b,...c->...a",f,a,b)


def anticommutator_completion(c):
    return np.array([[(a@b+b@a)/2 for b in c] for a in c])


def hermitian(rng,n):
    a=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n))
    return (a+a.conj().T)/2


def model(c,f,g,r=None,k=None,s=None,v=None,quartic=.2):
    d,n,_=c.shape
    return dict(C=c,R=c.copy() if r is None else r,
                K=anticommutator_completion(c) if k is None else k,
                S=np.zeros_like(c) if s is None else s,
                V=np.eye(n) if v is None else v,lam=quartic,f=f,g=g)


def empty_jet(n,d):
    return dict(phi=np.zeros(n,complex),p=np.zeros((4,n),complex),
                pp=np.zeros((4,4,n),complex),a=np.zeros((4,d)),da=np.zeros((4,4,d)))


def tangent(m,z,e,de,dde):
    r=m["R"]
    generator=np.einsum("a,aij->ij",e,r)
    dp=1j*np.einsum("ma,aij,j->mi",de,r,z["phi"])+1j*np.einsum("ij,mj->mi",generator,z["p"])
    return dict(phi=1j*generator@z["phi"],p=dp,
        a=de+m["g"]*bracket(m["f"],z["a"],e),
        da=dde+m["g"]*(bracket(m["f"],z["da"],e)+bracket(m["f"],z["a"][None,:,:],de[:,None,:])))


def density(m,z):
    phi,p,a=z["phi"],z["p"],z["a"]
    c,k,s=m["C"],m["K"],m["S"]
    answer=-sum(ETA[mu]*np.vdot(p[mu],p[mu]) for mu in range(4))
    for mu in range(4):
        for b in range(len(c)):
            answer-=1j*ETA[mu]*a[mu,b]*(np.vdot(phi,c[b]@p[mu])-np.vdot(p[mu],c[b]@phi))
            answer+=ETA[mu]*z["da"][mu,mu,b]*np.vdot(phi,s[b]@phi)
            for d in range(len(c)):
                answer-=ETA[mu]*a[mu,b]*a[mu,d]*np.vdot(phi,k[b,d]@phi)
    rho=np.vdot(phi,phi).real
    answer-=np.vdot(phi,m["V"]@phi)+m["lam"]*rho*rho
    return real(answer)


def variation(m,z,e,de,dde):
    """Independent directional derivative of every retained density term."""
    dz=tangent(m,z,e,de,dde)
    phi,p,a=z["phi"],z["p"],z["a"]
    vp,vpder,va=dz["phi"],dz["p"],dz["a"]
    c,k,s=m["C"],m["K"],m["S"]
    answer=-2*sum(ETA[mu]*np.vdot(p[mu],vpder[mu]).real for mu in range(4))
    for mu in range(4):
        for b in range(len(c)):
            bil=np.vdot(phi,c[b]@p[mu])-np.vdot(p[mu],c[b]@phi)
            dbil=(np.vdot(vp,c[b]@p[mu])+np.vdot(phi,c[b]@vpder[mu])
                  -np.vdot(vpder[mu],c[b]@phi)-np.vdot(p[mu],c[b]@vp))
            answer-=1j*ETA[mu]*(va[mu,b]*bil+a[mu,b]*dbil)
            answer+=ETA[mu]*(dz["da"][mu,mu,b]*np.vdot(phi,s[b]@phi)
                +z["da"][mu,mu,b]*(np.vdot(vp,s[b]@phi)+np.vdot(phi,s[b]@vp)))
            for d in range(len(c)):
                answer-=ETA[mu]*((va[mu,b]*a[mu,d]+a[mu,b]*va[mu,d])*np.vdot(phi,k[b,d]@phi)
                    +a[mu,b]*a[mu,d]*(np.vdot(vp,k[b,d]@phi)+np.vdot(phi,k[b,d]@vp)))
    answer-=2*np.vdot(phi,m["V"]@vp).real
    answer-=4*m["lam"]*np.vdot(phi,phi).real*np.vdot(phi,vp).real
    return real(answer)


def parameter_coefficients(m,z):
    """delta L = epsilon E + (partial_mu epsilon) D^mu + box(epsilon) Q."""
    c,r,k,s,f,g=m["C"],m["R"],m["K"],m["S"],m["f"],m["g"]
    phi,p,a=z["phi"],z["p"],z["a"]
    d=len(c)
    e_coeff=np.array([variation(m,z,np.eye(d)[b],np.zeros((4,d)),np.zeros((4,4,d))) for b in range(d)])
    d_coeff=np.zeros((4,d))
    q_coeff=np.array([real(np.vdot(phi,s[b]@phi)) for b in range(d)])
    for mu in range(4):
        for b in range(d):
            delta=r[b]-c[b]
            value=1j*(np.vdot(phi,delta@p[mu])-np.vdot(p[mu],delta@phi))
            for aa in range(d):
                h=c[aa]@r[b]+r[b]@c[aa]-2*k[aa,b]
                value+=a[mu,aa]*np.vdot(phi,h@phi)
                for cc in range(d):
                    value+=g*f[cc,aa,b]*a[mu,aa]*np.vdot(phi,s[cc]@phi)
            d_coeff[mu,b]=real(ETA[mu]*value)
    return e_coeff,d_coeff,q_coeff


def noether(m,z):
    """Analytic Euler coefficient after integration by parts in the parameter."""
    c,r,k,s,f,g=m["C"],m["R"],m["K"],m["S"],m["f"],m["g"]
    phi,p,a,da=z["phi"],z["p"],z["a"],z["da"]
    lap=np.einsum("m,mmi->i",ETA,z["pp"])
    diva=np.einsum("m,mma->a",ETA,da)
    e_coeff,_,_=parameter_coefficients(m,z)
    answer=e_coeff.astype(complex)
    for b in range(len(c)):
        delta=r[b]-c[b]
        divergence=1j*(np.vdot(phi,delta@lap)-np.vdot(lap,delta@phi))
        for aa in range(len(c)):
            h=c[aa]@r[b]+r[b]@c[aa]-2*k[aa,b]
            divergence+=diva[aa]*np.vdot(phi,h@phi)
            for mu in range(4):
                divergence+=ETA[mu]*a[mu,aa]*(np.vdot(p[mu],h@phi)+np.vdot(phi,h@p[mu]))
            for cc in range(len(c)):
                temp=diva[aa]*np.vdot(phi,s[cc]@phi)
                temp+=sum(ETA[mu]*a[mu,aa]*(np.vdot(p[mu],s[cc]@phi)+np.vdot(phi,s[cc]@p[mu])) for mu in range(4))
                divergence+=g*f[cc,aa,b]*temp
        boxq=np.vdot(lap,s[b]@phi)+np.vdot(phi,s[b]@lap)
        boxq+=2*sum(ETA[mu]*np.vdot(p[mu],s[b]@p[mu]) for mu in range(4))
        answer[b]+=-divergence+boxq
    return np.array([real(x) for x in answer])


def covariant_rewrite(m,z):
    """Used only after testing independent coefficients, to check sufficiency."""
    dp=z["p"]-1j*np.einsum("ma,aij,j->mi",z["a"],m["C"],z["phi"])
    rho=np.vdot(z["phi"],z["phi"]).real
    return real(-sum(ETA[mu]*np.vdot(dp[mu],dp[mu]) for mu in range(4))
                -np.vdot(z["phi"],m["V"]@z["phi"])-m["lam"]*rho*rho)


def closure_defect(m):
    c,f,g=m["C"],m["f"],m["g"]
    return np.array([[a@b-b@a-1j*g*np.einsum("c,cij->ij",f[:,aa,bb],c)
                      for bb,b in enumerate(c)] for aa,a in enumerate(c)])


def random_jet(rng,n,d):
    z=empty_jet(n,d)
    for name in ["phi","p","pp"]:
        z[name]=.25*(rng.normal(size=z[name].shape)+1j*rng.normal(size=z[name].shape))
    z["pp"]=(z["pp"]+np.swapaxes(z["pp"],0,1))/2
    z["a"]=.25*rng.normal(size=(4,d))
    z["da"]=.25*rng.normal(size=(4,4,d))
    return z


def periodic_fields(x):
    z=empty_jet(2,3)
    z["phi"]=np.array([.4+.2*np.sin(x)+.13j*np.cos(2*x),-.3+.17*np.cos(x)+.11j*np.sin(3*x)])
    z["p"][1]=[.2*np.cos(x)-.26j*np.sin(2*x),-.17*np.sin(x)+.33j*np.cos(3*x)]
    z["pp"][1,1]=[-.2*np.sin(x)-.52j*np.cos(2*x),-.17*np.cos(x)-.99j*np.sin(3*x)]
    z["a"][1]=[.22*np.cos(x),.18*np.sin(2*x),.16*np.cos(3*x)]
    z["da"][1,1]=[-.22*np.sin(x),.36*np.cos(2*x),-.48*np.sin(3*x)]
    return z


def spectral_derivative(values,order=1):
    freq=np.fft.fftfreq(len(values),d=1/len(values))
    return np.fft.ifft(np.fft.fft(values,axis=0)*(1j*freq[:,None])**order,axis=0).real


def compare(fresh,saved):
    if isinstance(fresh,dict):
        assert fresh.keys()==saved.keys()
        for key in fresh:
            compare(fresh[key],saved[key])
    elif isinstance(fresh,list):
        assert len(fresh)==len(saved)
        for a,b in zip(fresh,saved):
            compare(a,b)
    elif isinstance(fresh,float):
        assert math.isclose(fresh,saved,rel_tol=3e-8,abs_tol=3e-10),(fresh,saved)
    else:
        assert fresh==saved,(fresh,saved)


def run():
    t,f=su2()
    g=.5
    c=g*t
    standard=model(c,f,g)
    wrong_r=model(c,f,g,r=.75*t)
    wrong_k=model(c,f,g,k=1.5*anticommutator_completion(c))
    wrong_g=model(t,f,g)
    s=np.zeros_like(c);s[0]=.2*np.eye(2)
    wrong_s=model(c,f,g,s=s)
    wrong_v=model(c,f,g,v=np.diag([1.,2.]))
    invalid=[("R_minus_C",wrong_r),("independent_four_point",wrong_k),
             ("matter_vector_coupling_mismatch",wrong_g),
             ("nontrivial_effective_derivative_improvement",wrong_s),
             ("noninvariant_potential",wrong_v)]

    # Independent jet obstructions to the action identity, after integrating
    # derivatives off an arbitrary compactly supported gauge parameter.
    probes=[]
    z=empty_jet(2,3);z["phi"][0]=1;z["pp"][1,1,0]=1j
    probes.append(("R_minus_C",wrong_r,z,2,.25))
    z=empty_jet(2,3);z["phi"][0]=1;z["da"][1,1,2]=1
    probes.append(("independent_four_point",wrong_k,z,2,.0625))
    z=empty_jet(2,3);z["phi"][0]=1;z["p"][1,0]=1j;z["a"][1,0]=1
    probes.append(("matter_vector_coupling_mismatch",wrong_g,z,1,-.5))
    z=empty_jet(2,3);z["p"][1,0]=1
    probes.append(("nontrivial_effective_derivative_improvement",wrong_s,z,0,.4))
    z=empty_jet(2,3);z["phi"]=np.array([1,1j])/math.sqrt(2)
    probes.append(("noninvariant_potential",wrong_v,z,0,-.25))
    witnesses=[]
    for label,m,z,b,expected in probes:
        value=float(noether(m,z)[b])
        assert abs(value-expected)<1e-13,(label,value,expected)
        assert abs(value)>0
        witnesses.append(dict(obstruction=label,parameter_component=b,
                              noether_euler_coefficient=value,expected=expected,
                              proof_object="Euler coefficient after integration by parts, not delta L at one point"))

    # A second species may be neutral or couple with g; giving it a different
    # nonzero multiple of the same simple-factor generators fails the identity.
    two_c=np.array([np.kron(np.diag([g,1.]),a) for a in t])
    species_mismatch=model(two_c,f,g)
    z=empty_jet(4,3);z["phi"][2]=1;z["p"][1,2]=1j;z["a"][1,0]=1
    two_species_obstruction=float(noether(species_mismatch,z)[1])
    assert abs(two_species_obstruction+.5)<1e-13

    abelian=model(np.array([np.diag([.3,-.9])]),np.zeros((1,1,1)),.8,v=np.diag([2.,3.]))
    neutral=model(np.zeros((3,2,2),complex),f,g,v=np.diag([2.,3.]))
    multiplicity_mass=np.array([[2,.3j],[-.3j,3]],complex)
    reducible=model(np.array([np.kron(np.eye(2),a) for a in c]),f,g,
                    v=np.kron(multiplicity_mass,np.eye(2)))
    valid=[("nonAbelian_doublet",standard),("Abelian_unequal_charges",abelian),
           ("neutral_nonAbelian_spectator",neutral),("equivalent_copies_nontrivial_mass_commutant",reducible)]
    rng=np.random.default_rng(1016)
    maxima=dict(noether=0.,pointwise=0.,parameter_decomposition=0.,directional_difference=0.,completed_density=0.)
    samples=[]
    for name,m in valid:
        d,n,_=m["C"].shape
        assert norm(closure_defect(m))<1e-13
        assert norm([m["V"]@a-a@m["V"] for a in m["R"]])<1e-13
        for index in range(6):
            z=random_jet(rng,n,d)
            e=.3*rng.normal(size=d);de=.3*rng.normal(size=(4,d))
            dde=.3*rng.normal(size=(4,4,d));dde=(dde+np.swapaxes(dde,0,1))/2
            value=variation(m,z,e,de,dde)
            nc=noether(m,z)
            ec,dc,qc=parameter_coefficients(m,z)
            decomposition=float(e@ec+np.sum(de*dc)+np.einsum("m,mma,a->",ETA,dde,qc))
            dz=tangent(m,z,e,de,dde)
            step=1e-5
            plus={key:z[key]+step*dz[key] if key in dz else z[key] for key in z}
            minus={key:z[key]-step*dz[key] if key in dz else z[key] for key in z}
            difference=(density(m,plus)-density(m,minus))/(2*step)
            completion=abs(density(m,z)-covariant_rewrite(m,z))
            assert norm(nc)<2e-12 and abs(value)<2e-12
            assert abs(value-decomposition)<2e-12
            assert abs(difference-value)<2e-8
            assert completion<2e-12
            for key,val in [("noether",norm(nc)),("pointwise",abs(value)),
                            ("parameter_decomposition",abs(value-decomposition)),
                            ("directional_difference",abs(difference-value)),("completed_density",completion)]:
                maxima[key]=max(maxima[key],val)
            samples.append(dict(family=name,sample=index,noether_residual=norm(nc),
                                density_variation=value,directional_difference_residual=abs(difference-value)))

    # Check general coefficient identities too, so zeros from valid families
    # do not hide an incorrectly omitted term.
    for _ in range(6):
        random_c=np.array([hermitian(rng,2) for _ in range(3)])*.2
        random_r=np.array([hermitian(rng,2) for _ in range(3)])*.2
        random_k=np.array([[hermitian(rng,2) for _ in range(3)] for _ in range(3)])*.1
        random_k=(random_k+np.swapaxes(random_k,0,1))/2
        random_s=np.array([hermitian(rng,2) for _ in range(3)])*.1
        arbitrary=model(random_c,f,g,r=random_r,k=random_k,s=random_s,v=hermitian(rng,2))
        z=random_jet(rng,2,3)
        e=rng.normal(size=3);de=rng.normal(size=(4,3))
        dde=rng.normal(size=(4,4,3));dde=(dde+np.swapaxes(dde,0,1))/2
        value=variation(arbitrary,z,e,de,dde)
        ec,dc,qc=parameter_coefficients(arbitrary,z)
        decomposition=float(e@ec+np.sum(de*dc)+np.einsum("m,mma,a->",ETA,dde,qc))
        assert abs(value-decomposition)<2e-12
        maxima["parameter_decomposition"]=max(maxima["parameter_decomposition"],abs(value-decomposition))

    # Independently differentiate coefficient functions over one periodic
    # spatial coordinate. This verifies integration by parts and nonzero
    # action witnesses. It is an algebraic calibration, not a new spacetime.
    xs=2*np.pi*np.arange(96)/96
    fields=[periodic_fields(x) for x in xs]
    periodic=[]
    for label,m in [("valid",standard)]+invalid:
        coefficients=[parameter_coefficients(m,z) for z in fields]
        es=np.array([a[0] for a in coefficients])
        ds=np.array([a[1][1] for a in coefficients])
        qs=np.array([a[2] for a in coefficients])
        direct=np.array([noether(m,z) for z in fields])
        spectral=es-spectral_derivative(ds)+spectral_derivative(qs,2)
        error=float(np.max(np.abs(direct-spectral)))
        assert error<2e-10
        # Choose epsilon=N only to supply a finite integrated witness of a
        # nonzero Euler coefficient; this is not a new physical control rule.
        epsilon=direct
        depsilon=spectral_derivative(epsilon)
        ddepsilon=spectral_derivative(epsilon,2)
        variations=[]
        for idx,z in enumerate(fields):
            de=np.zeros((4,3));de[1]=depsilon[idx]
            dde=np.zeros((4,4,3));dde[1,1]=ddepsilon[idx]
            variations.append(variation(m,z,epsilon[idx],de,dde))
        action_change=float(2*np.pi*np.mean(variations))
        euler_square=float(2*np.pi*np.mean(np.sum(direct*direct,axis=1)))
        assert abs(action_change-euler_square)<2e-11
        if label!="valid":
            assert euler_square>1e-6
        periodic.append(dict(model=label,noether_differentiation_residual=error,
                             integrated_density_variation=action_change,
                             integrated_euler_square=euler_square,
                             integration_by_parts_residual=abs(action_change-euler_square)))

    # Both raw improvements may be nonzero and equal: together they are a
    # total divergence. Only S_eff is constrained in the chosen representative.
    t_improve=np.array([.2*np.eye(2),.1*t[0],.3*t[2]])
    boundary_current=[];raw_density=[]
    for z in fields:
        phi,p,a=z["phi"],z["p"][1],z["a"][1]
        bval=sum(a[b]*np.vdot(phi,t_improve[b]@phi) for b in range(3))
        dval=sum(z["da"][1,1,b]*np.vdot(phi,t_improve[b]@phi)
                 +a[b]*(np.vdot(p,t_improve[b]@phi)+np.vdot(phi,t_improve[b]@p)) for b in range(3))
        boundary_current.append(real(bval));raw_density.append(real(dval))
    boundary_derivative=spectral_derivative(np.array(boundary_current)[:,None])[:,0]
    boundary_error=float(np.max(np.abs(boundary_derivative-raw_density)))
    boundary_integral=float(2*np.pi*np.mean(raw_density))
    assert boundary_error<2e-12 and abs(boundary_integral)<2e-13

    sources=[BASE/"archive_1009_/research_note_1015.md",
             BASE/"archive_1009_/1015/vector_consistency_selection_results.json",
             BASE/"archive_1009_/1015/next_selection_audit.md",
             BASE/"archive_585_628/research_note_613.md",
             BASE/"archive_531_553/research_note_531.md",
             BASE/"archive_585_628/research_note_598.md",
             BASE/"archive_956_989/981/drafts/common_parent_contract_v1.md",
             BASE/"archive_1009_/1009/input_dependency_ledger_v0_1.md"]
    return dict(round=1016,new_calibration_groups=1,cumulative_test_groups=3794,
        new_cognitive_axioms=0,all_scientific_calibrations_passed=True,
        scope="canonical complex scalar; complex-linear Hermitian field generators; constant bilinear leading dimension-four vertices; exact local gauge action modulo total derivatives",
        complete_arbitrary_EFT_classification_claimed=False,
        group_or_matter_content_generated=False,physical_data_fit_claimed=False,
        pointwise_density_used_as_general_action_failure_proof=False,
        conditions_selected=["canonical derivative improvement S_eff=0","R=C",
            "K_ab={C_a,C_b}/2","[C_a,C_b]=i g f_ab^c C_c","potential invariant under R"],
        independent_noether_jet_witnesses=witnesses,
        two_species_mismatched_coupling=dict(vector_g=g,matter_scalings=[g,1.],
            second_species_noether_obstruction=two_species_obstruction,
            neutral_alternative_permitted=True),
        valid_families=[name for name,_ in valid],valid_finite_jet_samples=samples,
        arbitrary_independent_coefficient_samples=6,maximum_residuals=maxima,
        periodic_action_checks=periodic,
        derivative_improvement_quotient=dict(canonical_coefficient="S_eff=S_raw-T_raw",
            nonzero_S_equals_T_permitted=True,total_derivative_identity_residual=boundary_error,
            periodic_total_derivative_integral=boundary_integral),
        residual_freedom=dict(Abelian_charges=[.3,-.9],Abelian_charges_not_required_equal=True,
            neutral_matter_permitted=True,
            equivalent_copy_mass_eigenvalues=np.linalg.eigvalsh(reducible["V"]).tolist(),
            representation_and_field_count_not_selected=True,
            different_simple_ideal_couplings_not_related=True,
            global_charge_lattice_and_higher_operators_not_fixed=True),
        analytic_obligations="independent jets and local Noether Euler identity prove necessity for all fields; completed covariant derivative plus invariant potential proves sufficiency",
        historical_source_sha256={str(p.relative_to(BASE)).replace("\\","/"):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--write",action="store_true")
    args=parser.parse_args();result=run()
    if args.write:
        with OUT.open("x",encoding="utf8") as stream:
            json.dump(result,stream,ensure_ascii=False,indent=2);stream.write("\n")
    else:
        compare(result,json.loads(OUT.read_text(encoding="utf8")))
    print(json.dumps(dict(round=1016,passed=result["all_scientific_calibrations_passed"],
        noether_witnesses=result["independent_noether_jet_witnesses"],
        maximum_residuals=result["maximum_residuals"],
        periodic_action_checks=result["periodic_action_checks"],
        derivative_improvement_quotient=result["derivative_improvement_quotient"]),indent=2))
