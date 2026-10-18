#!/usr/bin/env python3
"""Independent 1310 audit: explicit Dirac null-polarization and angular rates.
Default/--check recompute without writing. --write-new creates results once.
No author calculation module is imported.
"""
import argparse, hashlib, json, math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
RESULT = HERE / "independent_results.json"

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

I2 = np.eye(2, dtype=complex)
Z2 = np.zeros((2, 2), dtype=complex)
PAULI = [
    np.array([[0, 1], [1, 0]], complex),
    np.array([[0, -1j], [1j, 0]], complex),
    np.array([[1, 0], [0, -1]], complex),
]
GAMMA = [np.block([[I2, Z2], [Z2, -I2]])]
GAMMA += [np.block([[Z2, x], [-x, Z2]]) for x in PAULI]
G5 = 1j * GAMMA[0] @ GAMMA[1] @ GAMMA[2] @ GAMMA[3]
PL, PR = (np.eye(4)-G5)/2, (np.eye(4)+G5)/2
METRIC = np.array([1., -1., -1., -1.])
UP, DOWN = np.array([1., 0.]), np.array([0., 1.])

def mdot(a, b):
    return np.sum(METRIC * a * b)

def slash(p):
    return sum(METRIC[j] * p[j] * GAMMA[j] for j in range(4))

def spinor(p, mass, chi):
    t = math.sqrt(p[0]+mass)
    sigmap = sum(p[j+1]*PAULI[j] for j in range(3))
    return np.concatenate([t*chi, sigmap@chi/t])

def current(out, inc):
    bar = out.conj() @ GAMMA[0]
    return np.array([bar @ gamma @ PL @ inc for gamma in GAMMA])

def transverse_spinor(pout, mass_out, inc, q):
    # Solve q_mu J^mu=0 directly in the outgoing 2-spinor,
    # rather than using the author's inverse T formula.
    basis = [spinor(pout, mass_out, UP), spinor(pout, mass_out, DOWN)]
    coeff = np.array([mdot(q, current(u, inc)) for u in basis])
    chi = np.array([-coeff[1].conjugate(), coeff[0].conjugate()])
    norm = np.linalg.norm(chi)
    assert norm > 1e-18
    chi /= norm
    return spinor(pout, mass_out, chi), chi

def sample(m, n, E, cosine):
    r, s = math.sqrt(E*E-m*m), math.sqrt(E*E-n*n)
    sine = math.sqrt(1-cosine*cosine)
    momenta = [np.array([E,0.,0.,r]), np.array([E,0.,0.,-r]),
               np.array([E,s*sine,0.,s*cosine]),
               np.array([E,-s*sine,0.,-s*cosine])]
    p1,p2,p3,p4 = momenta
    q = p1-p3
    u1, u2 = spinor(p1,m,UP), spinor(p2,m,DOWN)
    u3, chi3 = transverse_spinor(p3,n,u1,q)
    u4, chi4 = transverse_spinor(p4,n,u2,-q)
    currents = [current(u3,u1),current(u4,u2)]
    z = float(-mdot(q,q))
    contraction = mdot(*currents)
    dirac = max(float(np.linalg.norm((slash(p)-mass*np.eye(4))@u))
                for p,mass,u in zip(momenta,[m,m,n,n],[u1,u2,u3,u4]))
    ward = max(abs(mdot(q,j)) for j in currents)
    energy = max(abs(np.vdot(u,u).real-2*E) for u in [u1,u2,u3,u4])
    goldstone = []
    for out,inc in [(u3,u1),(u4,u2)]:
        bar = out.conj() @ GAMMA[0]
        goldstone.append(m*(bar@PR@inc)-n*(bar@PL@inc))
    # Separate analytic comparison; phases of the direct null vectors
    # are arbitrary, so compare the magnitude of J1.J2.
    A, B = s*sine, r-s*cosine
    b = math.sqrt(E+n)
    T = b*I2-(p3[1]*PAULI[0]+p3[3]*PAULI[2])/b
    N = np.linalg.norm(np.linalg.solve(T,np.array([A,B])))
    f = (math.sqrt(E+m)-r/math.sqrt(E+m))/(2*N)
    analytic = 2*f*f*z
    return {
        "z":z, "jj":contraction, "analytic_abs_jj":analytic,
        "dirac":dirac, "ward":float(ward), "norm":float(energy),
        "goldstone":float(max(map(abs,goldstone))),
        "current":currents, "momenta":momenta,
    }

def coefficient(g,z,v):
    return g*g/(2*(z+g*g*v*v/4))

def angular_rate(m,n,E,v,g,c2,order=64):
    # Direct dOmega integration, not the author's z quadrature.
    nodes,weights=np.polynomial.legendre.leggauss(order)
    cosines=nodes/2
    weights=weights/2
    r,s=math.sqrt(E*E-m*m), math.sqrt(E*E-n*n)
    values=[]
    kernels=[]
    for x in cosines:
        point=sample(m,n,E,float(x))
        jj2=abs(point["jj"])**2
        base=c2*jj2*(s/r)/(128*math.pi*E*E)
        kernels.append(base)
        values.append(base*coefficient(g,point["z"],v)**2)
    return float(np.dot(weights,values)), float(np.dot(weights,kernels))

def bounds(m,n,E,c2):
    r,s=math.sqrt(E*E-m*m),math.sqrt(E*E-n*n)
    z0,z1=r*r+s*s-r*s,r*r+s*s+r*s
    u=math.sqrt(E+m)-r/math.sqrt(E+m)
    b=math.sqrt(E+n)
    hm,hp=u*(b-s/b)/2,u*(b+s/b)/2
    lam=16*E*E*r*r
    km=(z1-z0)*c2*(2*hm*hm)**2/(16*math.pi*lam)
    kp=(z1-z0)*c2*(2*hp*hp)**2/(16*math.pi*lam)
    return z0,z1,km,kp

def solve_monotone(fn,target,hi):
    lo=0.
    for _ in range(70):
        mid=(lo+hi)/2
        if fn(mid)<target: lo=mid
        else: hi=mid
    return (lo+hi)/2

def run():
    sources={
        "admission":STAGE/"_admission"/"round1310_admission.md",
        "author_proof":HERE/"proof.md",
        "author_check":HERE/"check.py",
        "author_results":HERE/"results.json",
        "note":STAGE/"research_note_1310.md",
        "next":HERE/"NEXT.md",
        "1309_proof":STAGE/"1309"/"proof.md",
        "1303_proof":STAGE.parent/"archive_1290_1307"/"1303"/"proof.md",
        "independent_check":Path(__file__).resolve(),
    }
    bindings={key:{"path":str(path.relative_to(STAGE.parent)),
                   "sha256":sha(path)} for key,path in sources.items()}
    # Gamma algebra and explicit covariant normalizations.
    clifford=max(float(np.max(np.abs(GAMMA[a]@GAMMA[b]+GAMMA[b]@GAMMA[a]
                     - (2*METRIC[a] if a==b else 0)*np.eye(4))))
                 for a in range(4) for b in range(4))
    assert clifford < 1e-13
    diagnostics=[]
    for m,n,E in [(3.,4.,5.),(.3,.4,1.),(2**-12,2**-11,1.)]:
        for cosine in [-.5,-.25,0.,.25,.5]:
            x=sample(m,n,E,cosine)
            ratio=abs(abs(x["jj"])-x["analytic_abs_jj"])/x["analytic_abs_jj"]
            assert x["dirac"]<2e-12 and x["ward"]<2e-12
            assert x["norm"]<2e-12 and x["goldstone"]<2e-12
            assert ratio<2e-8 and abs(x["jj"])>0
            diagnostics.append({"m":m,"n":n,"E":E,"cosine":cosine,
                "q2_spacelike":x["z"],"dirac_residual":x["dirac"],
                "ward_residual":x["ward"],"goldstone_residual":x["goldstone"],
                "relative_current_formula_error":float(ratio)})
    calibration=sample(3.,4.,5.,0.)
    assert abs(abs(calibration["jj"])-800/197)<1e-12

    # Full-propagator contrast for generic, non-transverse polarizations.
    point=sample(.3,.4,1.,0.)
    p1,p2,p3,p4=point["momenta"]
    q=p1-p3; z=point["z"]; v=1.7
    j1=current(spinor(p3,.4,UP),spinor(p1,.3,UP))
    j2=current(spinor(p4,.4,DOWN),spinor(p2,.3,DOWN))
    divergences=mdot(q,j1)*mdot(q,j2)
    limit=-2*divergences/(v*v*z)
    assert abs(limit)>1e-6
    generic=[]
    for g in [.1,.01,.001,.0001]:
        amp=coefficient(g,z,v)*(mdot(j1,j2)-divergences/(g*g*v*v/4))
        generic.append({"g":g,"abs_full_amplitude":float(abs(amp)),
                        "distance_to_longitudinal_limit":float(abs(amp-limit))})
    assert generic[-1]["distance_to_longitudinal_limit"]<1e-7

    m,n,E,v,c2=.3,.4,1.,1.,1/9
    z0,z1,km,kp=bounds(m,n,E,c2)
    rate32,kint32=angular_rate(m,n,E,v,.6,c2,32)
    rate64,kint=angular_rate(m,n,E,v,.6,c2,64)
    assert abs(rate32-rate64)<1e-12*rate64
    assert km<kint<kp
    assert km*coefficient(.6,z1,v)**2<rate64<kp*coefficient(.6,z0,v)**2
    rate=lambda g: angular_rate(m,n,E,v,g,c2,32)[0]
    samples=[rate(g) for g in [.1,.2,.4,.6,.8]]
    assert all(a<b for a,b in zip(samples,samples[1:]))
    ceiling=4*kint/v**4
    assert samples[-1]<ceiling

    # Fixed actual-kernel threshold and independent pointwise necessary bound.
    g_task=.45
    target=rate(g_task)
    pstar=.6; ell=-math.log1p(-pstar)
    luminosity=ell/target
    recovered=solve_monotone(rate,target,.8)
    astar=math.sqrt(target/kp)
    gstar=math.sqrt(2*astar*z0/(1-astar*v*v/2))
    assert abs(recovered-g_task)<1e-10 and 0<gstar<g_task
    # Removing a uniform exposure cap admits all positive couplings.
    exposures=[]
    for g in [.4,.2,.1,.05]:
        sigma=rate(g); exposure=ell/sigma
        prob=-math.expm1(-exposure*sigma)
        assert abs(prob-pstar)<2e-15
        exposures.append({"g":g,"sigma":sigma,"required_exposure":exposure})
    assert all(a["required_exposure"]<b["required_exposure"]
               for a,b in zip(exposures,exposures[1:]))

    # Actual-model discrepancy is a declared bound, not a computed error.
    epsilon=.1
    target_necessary=-math.log1p(-(pstar-epsilon))/luminosity
    assert 0<target_necessary<target
    aerr=math.sqrt(target_necessary/kp)
    gerr=math.sqrt(2*aerr*z0/(1-aerr*v*v/2))
    assert 0<gerr<gstar
    # Exact zero-endpoint integrated bound for a bounded kernel (not a loophole).
    zero_bin=[]
    for g in [.2,.1,.05]:
        zmax=2.
        direct=g*g/v**2-g**4/(4*(zmax+g*g*v*v/4))
        assert 0<direct<g*g/v**2
        zero_bin.append({"g":g,"K0_one_integral":direct,"upper_g2":g*g/v**2})

    # Coupled fixed-kernel task/RG conditions, including a counterexample to
    # using the loose pointwise floor as a sufficient task certificate.
    gmax=.8; L=50.
    Bactual=8*math.pi**2/L*(g_task**-2-gmax**-2)
    Bloose=8*math.pi**2/L*(gstar**-2-gmax**-2)
    assert Bloose>Bactual
    structure=[]
    loose_only=[]
    for Nc in range(3,80,4):
        for ng in range(3,41):
            beta=(-43+2*ng*(Nc+1))/6
            exact=beta<=Bactual
            loose=beta<=Bloose
            # For the exact smallest initial coupling, analytic x(t) remains
            # above gmax^{-2} iff the endpoint does (or beta<=0).
            end=g_task**-2-beta*L/(8*math.pi**2)
            actual_flow=(end>=gmax**-2)
            assert exact==actual_flow
            if exact: structure.append([Nc,ng])
            if loose and not exact: loose_only.append([Nc,ng])
    assert structure and loose_only

    # Common positive-mass/phase data; no rerun or new certification of the
    # frozen continuous-field 1303 phase experiment.
    mu=2**-12
    omega=np.exp(2j*math.pi/3)
    U=np.array([[omega**(a*b) for b in range(3)] for a in range(3)])/math.sqrt(3)
    masses=np.array([mu,2*mu,3*mu])
    charged=np.array([mu,2*mu,4*mu])
    M=U@np.diag(masses)@U.conj().T
    J=float(np.imag(U[0,0]*U[1,1]*U[0,1].conj()*U[1,0].conj()))
    c=U[1,0]*U[0,1].conj()
    assert np.min(np.linalg.eigvalsh(M))>0
    assert abs(J)>0.09 and abs(abs(c)-1/3)<1e-14
    assert masses[0]==charged[0] and masses[1]==charged[1]
    assert np.all(np.diff(charged)>0)
    # The same phase-compatible positive-mass point, with its OWN rate and
    # finite exposure. The earlier .3/.4 illustrative budget is not reused.
    p_v=math.sqrt(2.)
    p_rate32,p_k32=angular_rate(mu,2*mu,1.,p_v,.2,1/9,32)
    p_rate64,p_ktot=angular_rate(mu,2*mu,1.,p_v,.2,1/9,64)
    p_z0,p_z1,p_km,p_kp=bounds(mu,2*mu,1.,1/9)
    assert p_rate64>0 and p_km<=p_ktot*(1+1e-9) and p_ktot<p_kp
    # At this tiny-mass point the analytic lower bound is almost saturated;
    # roundoff in 4-spinor subtraction is allowed, not used as a proof.
    assert abs(p_rate32-p_rate64)<1e-8*p_rate64
    p_exposure=math.log(2)/p_rate64
    p_sufficient=math.log(2)/(p_km*coefficient(.2,p_z1,p_v)**2)
    assert math.isfinite(p_sufficient) and p_exposure<=p_sufficient
    more=[]
    for count in (4,5):
        Un=np.eye(count,dtype=complex);Un[:3,:3]=U
        R=np.eye(count)
        for j in range(3,count):
            th=.0002/(j+1);cs,sn=math.cos(th),math.sin(th)
            V=np.eye(count);V[2,2]=V[j,j]=cs;V[2,j]=sn;V[j,2]=-sn
            R=R@V
        Un=Un@R
        assert np.max(abs(Un[:,:2]-np.vstack([U[:,:2],np.zeros((count-3,2))])))<1e-14
        Hn=Un@np.diag((mu*np.arange(1,count+1))**2)@Un.conj().T
        seen={0}
        while True:
            expanded=seen|{j for i in seen for j in range(count) if abs(Hn[i,j])>1e-22}
            if expanded==seen:break
            seen=expanded
        assert len(seen)==count
        more.append({"ng":count,"first_two_mass_columns_fixed":True,
                     "full_flavor_graph_connected":True})

    return {"status":"PASS_WITH_SCOPE","bindings":bindings,
        "independent_method":"4-spinor Ward-null solution; angular dOmega quadrature; no author module imported",
        "clifford_residual":clifford,"transverse_samples":diagnostics,
        "calibration_abs_contraction":float(abs(calibration["jj"])),
        "generic_longitudinal_limit_abs":float(abs(limit)),"generic_samples":generic,
        "fixed_kernel":{"K_lower":km,"K_integral":kint,"K_upper":kp,
            "rate32":rate32,"rate64":rate64,"strict_monotone_samples":samples,
            "g_task_constructed":g_task,"g_task_recovered":recovered,
            "loose_necessary_gstar":gstar,"exposure_limit":luminosity,
            "probability_error_input":epsilon,"weakened_necessary_gstar":gerr,
            "RG_B_exact":Bactual,"RG_B_loose":Bloose,
            "exact_finite_scan":structure,"loose_only_examples":loose_only[:8]},
        "unbounded_exposure_counterfamily":exposures,
        "zero_endpoint_kernel_bound":zero_bin,
        "joint_phase_data":{"positive_neutrino_masses":masses.tolist(),
             "positive_charged_masses":charged.tolist(),"J":J,"abs_cc_mixing":float(abs(c)),
             "same_P_point":{"m":mu,"n":2*mu,"E":1.,"v":p_v,"g":.2,
                 "z0":p_z0,"z1":p_z1,"K_integral":p_ktot,"K_lower":p_km,"K_upper":p_kp,
                 "sigma32":p_rate32,"sigma64":p_rate64,"target_probability":.5,
                 "own_exposure":p_exposure,"own_analytic_sufficient_exposure":p_sufficient},
             "more_generation_transport":more},
        "scope":["Finite samples check the analytic proof, not a continuum proof.",
          "Poisson, exposure, polarization, matching and total errors are declared effective inputs.",
          "No certification of actual instrument, full interacting CP experiment, or higher-loop remainder."]}

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    parser.add_argument("--write-new",action="store_true")
    args=parser.parse_args()
    result=run()
    if args.write_new:
        with RESULT.open("x",encoding="utf-8",newline="\n") as f:
            json.dump(result,f,ensure_ascii=False,indent=2,allow_nan=False);f.write("\n")
    elif RESULT.exists():
        stored=json.loads(RESULT.read_text(encoding="utf-8"))
        if stored!=result:
            raise AssertionError("Frozen independent results or source bindings changed")
    print(json.dumps({"status":result["status"],"transverse_samples":len(result["transverse_samples"]),
         "bound_sources":len(result["bindings"]),"write_new":args.write_new,
         "K_integral":result["fixed_kernel"]["K_integral"],
         "g_task":result["fixed_kernel"]["g_task_recovered"]},ensure_ascii=False))
if __name__=="__main__":
    main()
