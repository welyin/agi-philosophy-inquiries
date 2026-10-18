"""1296 independent Nambu/CAR checks; no author or historical code imports.
Exact rational certificates and small independent finite-matrix diagnostics.
Default/--check read only; --write creates results exclusively.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse, hashlib, json, math
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=next(p for p in HERE.parents if (p/"research_cognition_physics").is_dir())
OUT=HERE/"independent_results.json"
ALPHA=F(31,2**27)
EPS=ALPHA/128

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def op(a): return float(np.linalg.norm(a,2))
def err(a,b): return float(np.max(np.abs(np.asarray(a)-np.asarray(b))))
def spectral(h,fn):
    w,v=np.linalg.eigh(h)
    return (v*fn(w))@v.conj().T
def evolve(h,t): return spectral(h,lambda x:np.exp(-1j*t*x))
def block(a,b):
    z=np.zeros((a.shape[0],b.shape[0]),complex)
    return np.block([[a,z],[z.conj().T,b]])
def bdg(h,d):
    assert err(d.T,-d)<1e-13
    return np.block([[h,d],[-d.conj(),-h.conj()]])
def swap(m):
    z=np.zeros((m,m));i=np.eye(m)
    return np.block([[z,i],[i,z]])
def encode(x):
    if isinstance(x,F):return {"numerator":x.numerator,"denominator":x.denominator}
    if isinstance(x,dict):return {k:encode(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [encode(v) for v in x]
    if isinstance(x,np.integer):return int(x)
    if isinstance(x,np.floating):return float(x)
    return x

def exact_budget():
    e=F(7,256)*ALPHA;b=ALPHA/256
    prob=ALPHA**2-2*ALPHA*e-b*b
    hs=(ALPHA-e)**2-b*b
    assert prob==F(61951,65536)*ALPHA**2
    assert hs==F(3875,4096)*ALPHA**2
    # All finite tests below are diagnostics; this bound uses the inherited continuum bump estimate.
    mass_upper=1+6*EPS+EPS
    energy_squared=4*(F(524288)+mass_upper**2)
    assert mass_upper<2 and energy_squared<1450**2
    return dict(alpha=ALPHA,epsilon=EPS,delta_interval=[F(0),EPS],
                normal_error=e,anomalous_error=b,
                probability_factor=prob/ALPHA**2,probability_lower=prob,
                HS_injectivity_factor=hs/ALPHA**2,HS_injectivity_lower=hs,
                relative_energy_squared_upper=energy_squared,
                inherited_relative_energy_strict_upper=1450,
                whole_window_by_analysis=True)

def fourier_diagnostics():
    sx=np.array([[0,1],[1,0]],complex)
    sy=np.array([[0,-1j],[1j,0]],complex)
    sz=np.diag([1.,-1.])
    mom=np.array([.23,-.41,.17])
    spin=mom[0]*sx+mom[1]*sy+mom[2]*sz
    n=3; e=float(EPS)
    d=np.diag(np.arange(1,4))
    path=np.diag(np.ones(2),1)+np.diag(np.ones(2),-1)
    mass=np.eye(n)+e*(d+path)
    mm=np.kron(np.eye(2),mass)
    ss=np.kron(spin,np.eye(n))
    hplus=np.block([[-ss,mm],[mm,ss]])
    hminus=np.block([[ss,mm],[mm,-ss]])
    h=block(hplus,hminus)
    hbase=block(np.block([[-ss,np.eye(6)],[np.eye(6),ss]]),
                np.block([[ss,np.eye(6)],[np.eye(6),-ss]]))
    # Translation-invariant pairing couples +p with -p, including the spin antisymmetry.
    delta0=np.zeros((12,12),complex)
    delta0[6:,6:]=np.kron(1j*sy,np.eye(n))
    pair= np.block([[np.zeros((12,12)),delta0],[delta0,np.zeros((12,12))]])
    assert err(pair.T,-pair)==0 and abs(op(pair)-1)<1e-14
    zeta=np.array([1,1j])/math.sqrt(2)
    vin=np.zeros((24,2),complex);vout=np.zeros((24,2),complex)
    for j in range(2):
        v=np.kron(zeta,np.eye(n)[:,j])/math.sqrt(2)
        vin[6:12,j]=v;vin[18:24,j]=v
        vout[0:6,j]=v;vout[12:18,j]=v
    ff=vin@vin.conj().T;fh=block(ff,ff.conj())
    cc=swap(24);charge=np.diag([1.]*24+[-1.]*24)
    h0=bdg(h,np.zeros_like(h));hb=bdg(hbase,np.zeros_like(h))
    data=[];mx=0.
    for fraction in (F(0),F(1,2),F(1)):
        delta=float(EPS*fraction)
        hh=bdg(h,delta*pair)
        pp=spectral(hh,lambda x:(x<0).astype(float))
        absh=spectral(hh,np.abs)
        assert np.min(np.abs(np.linalg.eigvalsh(hh)))>.9
        residual=max(err(hh,hh.conj().T),err(cc@hh.conj()@cc,-hh),
            err(cc@pp.conj()@cc,np.eye(48)-pp),err(pp@pp,pp),
            err(cc@fh.conj()@cc,fh),err(fh@hh@fh,0),
            err(fh@fh,fh))
        u=evolve(hh,.5);u0=evolve(h0,.5);ub=evolve(hb,.5)
        aa=vout.conj().T@u[:24,:24]@vin
        bb=vout.conj().T@u[:24,24:]@vin.conj()
        a0=vout.conj().T@ub[:24,:24]@vin
        normal=op(aa-a0);anomalous=op(bb)
        assert normal<=float(7*ALPHA/256)+2e-14
        assert anomalous<=delta/2+2e-14
        pairing_distance=op(u-u0)
        assert pairing_distance<=delta/2+2e-14
        comm=op(hh@charge-charge@hh)
        assert abs(comm-2*delta)<2e-14
        rho=np.array([[.43,.12+.09j],[.12-.09j,.57]])
        code=block(vin@rho@vin.conj().T,
                   vin.conj()@(np.eye(2)-rho.T)@vin.T)
        qq=np.eye(48)-fh
        reset=qq@pp@qq+code
        assert err(cc@reset.conj()@cc,np.eye(48)-reset)<2e-13
        energy=float((np.trace(hh@(reset-pp))/2).real)
        half=float((np.trace(fh@absh)/2).real)
        residual=max(residual,abs(energy-half))
        assert energy>0
        # Local algebraic S-phase Ward identity: a compensating background phase balances charge.
        sim=bdg(np.zeros_like(h),1j*delta*pair)
        ward=err(1j*(charge@hh-hh@charge),2*sim)
        residual=max(residual,ward);mx=max(mx,residual)
        data.append(dict(delta_fraction=str(fraction),normal_leg_distance=normal,
            anomalous_leg_norm=anomalous,pairing_Duhamel_distance=pairing_distance,
            fixed_fermion_charge_commutator=comm,half_weight_energy=half,
            reset_relative_energy=energy,S_phase_compensation_residual=ward,
            max_algebraic_residual=residual))
    assert mx<2e-12
    return dict(real_particle_Cauchy_dimension=24,Nambu_dimension=48,
        physical_code_modes=2,Nambu_code_projection_rank=4,
        symmetric_momentum_pair=mom.tolist(),candidate_samples=data,
        max_residual=mx,not_a_continuum_kernel_computation=True,
        no_extra_species_from_Nambu=True)

def car(m):
    ans=[]
    for j in range(m):
        a=np.zeros((2**m,2**m),complex)
        mask=1<<(m-1-j)
        for col in range(2**m):
            if col&mask:
                # Fermion order uses earlier (higher binary) slots.
                parity=sum((col>>(m-1-k))&1 for k in range(j))
                a[col^mask,col]=(-1)**parity
        ans.append(a)
    return ans

def many_quadratic(j,ops):
    z=ops+[a.conj().T for a in ops]
    out=np.zeros_like(ops[0])
    for i in range(len(z)):
        for k in range(len(z)):
            out+=.5*j[i,k]*(z[i].conj().T@z[k])
    return out

def covariance(state,ops):
    z=ops+[a.conj().T for a in ops]
    return np.array([[np.trace(state@z[j].conj().T@z[i])
                      for j in range(len(z))] for i in range(len(z))])

def trace_code(state):
    # Physical ordering: two code modes, two environment modes.
    return np.einsum("aiaj->ij",state.reshape(4,4,4,4))

def finite_car_diagnostic():
    # Deliberately larger coefficients test the identities without tiny-pair cancellation.
    # This is an algebra diagnostic, not another candidate or replacement physical model.
    ops=car(4);z=ops+[a.conj().T for a in ops]
    h=np.array([[0,0,1,.21],[0,0,-.16,1.13],
                [1,-.16,.34,.07j],[.21,1.13,-.07j,-.27]],complex)
    pair=np.array([[0,0,.17,.03],[0,0,-.04,.13],
                   [-.17,.04,0,.11],[-.03,-.13,-.11,0]],complex)
    hh=bdg(h,pair);ham=many_quadratic(hh,ops)
    vals,vecs=np.linalg.eigh(ham);ground=vecs[:,0:1]
    reference=ground@ground.conj().T
    pp=spectral(hh,lambda x:(x<0).astype(float))
    ref_res=err(covariance(reference,ops),pp)
    assert ref_res<1e-12 and vals[1]-vals[0]>.2
    env=trace_code(reference)
    assert abs(np.trace(env)-1)<1e-13
    pp_env=covariance(env,car(2))
    anomalous_environment=op(pp_env[:2,2:])
    assert anomalous_environment>1e-3
    f=np.diag([1,1,0,0]);fh=block(f,f)
    qq=np.eye(8)-fh;cc=swap(4)
    assert err(fh@hh@fh,0)<1e-14
    # A true pure entangled Q-R state, with nonzero code coherence.
    angle=.37
    uq=np.array([[math.cos(angle),-math.sin(angle)],[math.sin(angle),math.cos(angle)]])
    ur=np.array([[1,1j],[1j,1]])/math.sqrt(2)
    vqr=np.kron(uq,ur)@np.array([math.sqrt(.38),0,0,math.sqrt(.62)],complex)
    qr=np.outer(vqr,vqr.conj())
    rq=np.einsum("aibi->ab",qr.reshape(2,2,2,2))
    rr=np.einsum("aiaj->ij",qr.reshape(2,2,2,2))
    w=np.zeros((4,2),complex);w[2,0]=1;w[1,1]=1
    def code_block(i,j):
        return np.kron(np.outer(w[:,i],w[:,j].conj()),env)
    joint=np.zeros((32,32),complex)
    for i in range(2):
        for j in range(2):
            joint+=np.kron(code_block(i,j),qr.reshape(2,2,2,2)[i,:,j,:])
    state=np.einsum("aibi->ab",joint.reshape(16,2,16,2))
    remoteR=np.einsum("aijakl->ijkl",joint.reshape(4,4,2,4,4,2)).reshape(8,8)
    blind=err(remoteR,np.kron(env,rr))
    ccode=np.zeros((8,8),complex)
    ccode[:2,:2]=rq;ccode[4:6,4:6]=np.eye(2)-rq.T
    cstate=qq@pp@qq+ccode
    mx=max(ref_res,blind,err(covariance(state,ops),cstate),
           err(cc@cstate.conj()@cc,np.eye(8)-cstate))
    double=float(np.trace(state@ops[0].conj().T@ops[0]@ops[1].conj().T@ops[1]).real)
    wrong=float(np.linalg.det(rq).real)
    assert abs(double)<1e-14 and wrong>.2
    # Direct Fock evolution independently checks the normal/anomalous CAR legs.
    t=.31;uf=evolve(ham,t);u=evolve(hh,t)
    for j in range(8):
        evolved=uf.conj().T@z[j]@uf
        rhs=sum((u[j,k]*z[k] for k in range(8)),np.zeros((16,16),complex))
        mx=max(mx,err(evolved,rhs))
    aa=u[2:4,0:2];bb=u[2:4,4:6]
    out_effect=ops[2].conj().T@ops[2]
    probs=[]
    for j in range(2):
        branch=code_block(j,j)
        probs.append(float(np.trace(uf@branch@uf.conj().T@out_effect).real))
    difference=probs[0]-probs[1]
    formula=abs(aa[0,0])**2-abs(aa[0,1])**2-abs(bb[0,0])**2+abs(bb[0,1])**2
    wrong_formula=abs(aa[0,0])**2-abs(aa[0,1])**2
    mx=max(mx,abs(difference-formula))
    assert abs(difference-wrong_formula)>1e-5
    # Full matrix response, using actual four code matrix units, including off-diagonal ones.
    baseline_cov=qq@pp@qq+block(np.zeros((4,4)),f)
    for i in range(2):
        for j in range(2):
            ex=np.zeros((2,2),complex);ex[i,j]=1
            actual=covariance(uf@code_block(i,j)@uf.conj().T,ops)[2:4,2:4]
            pred=(u@baseline_cov@u.conj().T)[2:4,2:4]*np.trace(ex)
            pred+=aa@ex@aa.conj().T-bb@ex.T@bb.conj().T
            mx=max(mx,err(actual,pred))
    # Fixed sources include two S directions and a separate normal bilinear.
    normal=np.array([[.2,.1j,.3,0],[-.1j,-.4,.2,.07],
                     [.3,.2,.6,-.13j],[0,.07,.13j,-.2]],complex)
    sources=[("S_real",bdg(np.zeros((4,4)),pair)),
             ("S_imaginary",bdg(np.zeros((4,4)),1j*pair)),
             ("normal_source",bdg(normal,np.zeros((4,4)))),
             ("energy",hh)]
    source_data=[];cabs=3/7
    br=np.array([[.3,.2j],[-.2j,-.7]])
    for name,js in sources:
        jt=u.conj().T@js@u
        obs=uf.conj().T@many_quadratic(js,ops)@uf+cabs*np.eye(16)
        # Direct pullback from the true code/environment state.
        pulled=np.array([[np.trace(code_block(j,i)@obs) for j in range(2)] for i in range(2)])
        baseline=.5*np.trace(jt@baseline_cov)+cabs
        kernel=.5*(jt[:2,:2]-jt[4:6,4:6].T)
        predicted=baseline*np.eye(2)+kernel
        residual=max(err(pulled,predicted),
            abs(np.trace(joint@np.kron(obs,br))-np.trace(qr@np.kron(predicted,br))))
        mx=max(mx,residual)
        source_data.append(dict(name=name,max_pullback_residual=float(residual),
              traceless_logical_norm=op(pulled-np.trace(pulled)*np.eye(2)/2),
              baseline_real=float(baseline.real)))
    resetenergy=float(np.trace((state-reference)@ham).real)
    halfenergy=float((np.trace(fh@spectral(hh,np.abs))/2).real)
    mx=max(mx,abs(resetenergy-halfenergy))
    assert resetenergy>0
    charge=np.diag([1.]*4+[-1.]*4)
    charge_defect=op(hh@charge-charge@hh)
    assert charge_defect>0
    ward=err(1j*(charge@hh-hh@charge),2*bdg(np.zeros((4,4)),1j*pair))
    mx=max(mx,ward)
    assert mx<3e-12
    return dict(real_CAR_modes=4,Fock_dimension=16,partner_dimension=2,
         generic_coefficients_only_algebra_diagnostic=True,
         nonquasifree_actual_double_occupation=double,
         incorrect_Wick_double_occupation=wrong,
         environment_anomalous_norm=anomalous_environment,
         full_initial_remote_R_product_residual=blind,
         normal_leg_norm=op(aa),anomalous_leg_norm=op(bb),
         actual_receiver_probability_difference=difference,
         both_leg_formula_difference=float(formula),
         omitted_pairing_error=float(abs(difference-wrong_formula)),
         sources=source_data,reset_relative_energy=resetenergy,
         half_weight_energy=halfenergy,
         incorrect_no_half_energy=2*halfenergy,
         fixed_fermion_charge_defect=charge_defect,
         compensating_background_phase_Ward_residual=ward,
         max_residual=mx,
         no_source_variance_or_apparatus_budget_claim=True)

def bindings():
    files={
        "author_proof":HERE/"proof.md","author_check":HERE/"check.py",
        "author_results":HERE/"results.json","independent_check":Path(__file__),
        "admission":HERE.parent/"_shared/round1296_admission.md",
        "1295_proof":HERE.parent/"1295/proof.md",
        "730_note":ROOT/"research_cognition_physics/archive_702_741/research_note_730.md"}
    return {k:{"path":str(p.relative_to(ROOT)).replace("\\","/"),"sha256":sha(p)}
            for k,p in files.items()}
def compare(a,b):
    if isinstance(a,dict):
        assert isinstance(b,dict) and a.keys()==b.keys()
        for k in a:compare(a[k],b[k])
    elif isinstance(a,list):
        assert isinstance(b,list) and len(a)==len(b)
        for x,y in zip(a,b):compare(x,y)
    elif isinstance(a,float):
        assert math.isfinite(a) and math.isfinite(b)
        assert math.isclose(a,b,rel_tol=2e-8,abs_tol=5e-14),(a,b)
    else:assert a==b,(a,b)
def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    parser.add_argument("--write",action="store_true")
    parser.add_argument("--diagnose",action="store_true")
    args=parser.parse_args()
    if sum((args.check,args.write,args.diagnose))>1:parser.error("choose one mode")
    payload=encode(dict(round=1296,independent_check_groups=1,
       verdict="PASS_INDEPENDENT_SCOPED",
       exact_budget=exact_budget(),candidate_fourier=fourier_diagnostics(),
       independent_full_CAR=finite_car_diagnostic(),
       scope=dict(continuum_quantifiers_by_analytical_proof=True,
           own_selfdual_reference_and_full_environment=True,
           logical_code_dimension=2,Nambu_code_dimension=4,
           response_difference_not_a_standalone_CP_channel=True,
           fixed_background_quadratic_class=True,
           full_interacting_SM_certified=False,actual_apparatus_certified=False,
           uniform_unbounded_source_error_certified=False)))
    if args.diagnose:
        print("DIAGNOSTIC_ONLY_NO_WRITE")
        print(json.dumps({"exact_budget":payload["exact_budget"],
           "Fourier_max_residual":payload["candidate_fourier"]["max_residual"],
           "CAR_max_residual":payload["independent_full_CAR"]["max_residual"],
           "pairing_omission_error":payload["independent_full_CAR"]["omitted_pairing_error"]},
           ensure_ascii=False,indent=2))
        return
    payload["source_bindings"]=bindings()
    if args.write:
        with OUT.open("x",encoding="utf-8",newline="\n") as f:
            json.dump(payload,f,ensure_ascii=False,indent=2,sort_keys=True);f.write("\n")
        print("PASS_INDEPENDENT_SCOPED EXCLUSIVE_WRITE")
    else:
        compare(payload,json.loads(OUT.read_text(encoding="utf-8")))
        print("PASS_INDEPENDENT_SCOPED READ_ONLY_MATCH")
    print("probability factor 61951/65536; HS factor 3875/4096; real code 2/Nambu code 4")
    print("algebra max residual",max(payload["candidate_fourier"]["max_residual"],
        payload["independent_full_CAR"]["max_residual"]))
if __name__=="__main__":main()
