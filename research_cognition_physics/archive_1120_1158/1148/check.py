"""1148 finite algebra and rational-budget checks; not a field simulation."""
from pathlib import Path
from fractions import Fraction as Q
import argparse, hashlib, json, math, importlib.util
import numpy as np
HERE=Path(__file__).resolve().parent
def pack(q): return {"fraction":str(q),"float":float(q)}
def budget():
    p=Q(1,39999)
    old=4*p+Q(22,7)*(Q(1,10000)+(1+Q(1,10000))*Q(1,37500))+Q(2,1000)
    # Uniform upper bounds from parameter selection; h=1, Delta_mu<1.
    aw=Q(41,10**8)
    added=Q(4,10**6)+2*aw
    # sqrt(2)<1.415 suffices; no sampled parameters.
    position=Q(1415,1000)*Q(405,10**9)/Q(1,10)
    full=old+added+position
    clock=Q(121,10**6)+added+position+Q(2,1000)
    # sin(5pi/8)=sqrt(2+sqrt2)/2<.924, hence p_wrong<.962.
    assert Q(1415,1000)**2>2
    assert (Q(924,1000)*2)**2 > 2+Q(1415,1000)
    sqrt2_hi=Q(1414214,10**6); sin_hi=Q(9239,10000)
    assert sqrt2_hi**2>2 and (2*sin_hi)**2>2+sqrt2_hi
    wrong_actual=(1+sin_hi)/2+5*Q(22,7)/1200000
    wrong=Q(962,1000)
    assert wrong_actual<wrong
    gap=1-2*full
    clock_gap=1-clock-wrong
    assert full<Q(252,100000) and gap>Q(9949,10000)
    assert clock_gap>Q(35,1000)
    assert added<Q(483,10**8) and position<Q(5731,10**9)
    return {"name":"exact_uniform_budget","pass":True,
      "prior_record_instrument_bound":pack(old),"new_motion_bound":pack(added),
      "position_instrument_bound":pack(position),"full_instrument_bound":pack(full),
      "record_probability_gap_lower":pack(gap),"clock_total_bound":pack(clock),
      "wrong_clock_gap_lower":pack(clock_gap),"actual_frequency_wrong_probability_upper":pack(wrong_actual),
      "no_concrete_heat_trace_gap_or_gate_time_substituted":True,
      "finite_parameter_recipe":"proof (11), Gp/Go independent of I/M0"}
def spectral_window():
    # A three-point mass spectral measure with a genuine small far tail.
    # This is a 1+1 spin-zero numerical diagnostic of (8), not the material.
    eta=math.log(2); ch=math.cosh(eta); sh=math.sinh(eta)
    mu=50.; masses=np.array([49.99,50.01,65.])
    tail=1e-4; coeff=np.array([math.sqrt((1-tail**2)/2)]*2+[tail])
    sigma=.15; h=.5; lower=49.9
    rho=float(np.linalg.norm((masses-mu)*coeff))
    mask=np.abs(masses-mu)>h
    tail_norm=float(np.linalg.norm(coeff[mask]))
    assert tail_norm<=rho/h
    # Analytic normalized momentum Gaussian; Y norm = sigma in one dimension.
    norm=(2*sigma*sigma/math.pi)**.25
    def boosted(x,m):
        E=np.sqrt(m*m+x*x)
        k=ch*x-sh*E
        Ek=np.sqrt(m*m+k*k)
        return np.sqrt(Ek/E)*norm*np.exp(-sigma*sigma*k*k)
    samples=[]
    for count in (30001,60001):
        x=np.linspace(-30.,110.,count); step=x[1]-x[0]
        ref=boosted(x,mu)
        v=np.array([boosted(x,m) for m in masses])
        trap=getattr(np,"trapezoid",np.trapz)
        norms=trap(v*v,x,axis=1)
        actual=float(np.sqrt(np.sum(coeff**2*trap((v-ref)**2,x,axis=1))))
        near=float(np.sqrt(np.sum(coeff[~mask]**2*trap((v[~mask]-ref)**2,x,axis=1))))
        bound=2*rho/h+h*(2*sigma+1/lower)
        assert np.max(np.abs(norms-1))<1e-10 and actual<bound
        assert near<=h*(2*sigma+1/lower)
        # Exact Lorentz endpoint phase identity, without dropping spatial phase.
        T=.37
        residual=0.
        for m in masses:
            E=np.sqrt(m*m+x*x); k=ch*x-sh*E; Ek=np.sqrt(m*m+k*k)
            lhs=np.exp(-1j*ch*T*E)
            rhs=np.exp(-1j*sh*T*x)*np.exp(-1j*T*Ek)
            residual=max(residual,float(np.max(np.abs(lhs-rhs))))
        assert residual<1e-12
        samples.append({"count":count,"norm_max_error":float(np.max(np.abs(norms-1))),
            "boost_difference":actual,"near_window_difference":near,
            "proven_bound":bound,"endpoint_phase_residual":residual})
    assert abs(samples[0]["boost_difference"]-samples[1]["boost_difference"])<1e-10
    return {"name":"normal_mass_spectral_tail_and_continuous_boost_diagnostic",
        "pass":True,"mu":mu,"masses":masses.tolist(),"rho":rho,"window_h":h,
        "actual_tail_norm":tail_norm,"tail_certificate":rho/h,"samples":samples,
        "full_3D_field_and_spin_claim":False}
def car():
    # Four actual positive-mode CAR algebra diagnostic: two flavors x two spin.
    dim=16
    def annih(j):
        a=np.zeros((dim,dim),complex)
        for n in range(dim):
            if n>>j&1:
                a[n^(1<<j),n]=(-1)**((n&((1<<j)-1)).bit_count())
        return a
    b=[annih(j) for j in range(4)]; eye=np.eye(dim)
    phi=.73
    modes=[(b[s]+np.exp(1j*phi)*b[2+s])/math.sqrt(2) for s in range(2)]
    ns=[m.conj().T@m for m in modes]
    E=eye-(eye-ns[0])@(eye-ns[1])
    paulis=[np.array([[0,1],[1,0]],complex),np.array([[0,-1j],[1j,0]],complex),np.diag([1,-1])]
    spins=[]
    for pauli in paulis:
        S=np.zeros((dim,dim),complex)
        for a in range(2):
            for s in range(2):
                for t in range(2):
                    S+=.5*pauli[s,t]*b[2*a+s].conj().T@b[2*a+t]
        spins.append(S)
    W=np.zeros((16,4)); W[[1,2,4,8],np.arange(4)]=1
    v=np.array([1,np.exp(-1j*phi)])/math.sqrt(2)
    target=np.kron(np.outer(v,v.conj()),np.eye(2))
    errs={"projection":float(np.linalg.norm(E@E-E,2)),
          "rotation":max(float(np.linalg.norm(E@S-S@E,2)) for S in spins),
          "code":float(np.linalg.norm(W.conj().T@E@W-target,2)),
          "leakage":float(np.linalg.norm(E@W-W@target,2))}
    assert max(errs.values())<1e-13
    single=float(np.linalg.norm(ns[0]@spins[0]-spins[0]@ns[0],2))
    count_norm=float(np.linalg.norm(ns[0]+ns[1],2))
    assert abs(single-.5)<1e-13 and abs(count_norm-2)<1e-13
    # Full two-outcome isometry and spectator-complete compression.
    R=np.vstack([E,eye-E]); R0=np.vstack([target,np.eye(4)-target])
    assert np.linalg.norm(R.conj().T@R-eye)<1e-13
    WW=np.kron(np.eye(2),W)
    full=float(np.linalg.norm(R@W-WW@R0,2))
    assert full<1e-13
    phase_rows=[]
    for dp in (.001,.1,.7):
        other=[(b[s]+np.exp(1j*(phi+dp))*b[2+s])/math.sqrt(2) for s in range(2)]
        nn=[m.conj().T@m for m in other]
        Ep=eye-(eye-nn[0])@(eye-nn[1])
        error=float(np.linalg.norm(np.vstack([Ep,eye-Ep])-R,2))
        exact=math.sqrt(2*(1-math.cos(dp/2)**4))
        assert abs(error-exact)<1e-12 and error<=abs(dp)+1e-14
        phase_rows.append({"phase_error":dp,"full_instrument_norm":error,"exact":exact})
    return {"name":"spin_blind_full_CAR_effect_and_complete_instrument","pass":True,
      "full_Fock_phase_tolerance_checks":phase_rows,
      "errors":errs,"instrument_compression_residual":full,
      "single_spin_bad_commutator":single,"summed_occupancy_norm_not_effect":count_norm,
      "not_a_sea_or_interacting_field_simulation":True}
def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    manifest=HERE/'source_hashes.json'
    sources=json.loads(manifest.read_text(encoding='utf-8'))
    for name,digest in sources.items():
        assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest,name
    result={"round":1148,"all_passed":True,"check_groups":3,
            "scientific_check_increment":0,"source_sha256":sources,
            "groups":[budget(),spectral_window(),car()]}
    if args.write_results:
        (HERE/'results.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:
        saved=json.loads((HERE/'results.json').read_text(encoding='utf-8'))
        assert saved==result,"saved result differs"
    print(json.dumps({"round":1148,"all_passed":True,"check_groups":3,
        "readonly":not args.write_results,"source_count":len(sources),
        "full_error_upper":result["groups"][0]["full_instrument_bound"]["float"]}))
if __name__=="__main__": main()
