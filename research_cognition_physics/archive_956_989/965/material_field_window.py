"""965: unchanged 958 electronic material coupled to one quantum EM mode.
The infinite-occupation Hamiltonian is primary; finite Fock computations are
compared to it by an explicit residual/Duhamel bound, not just cutoff agreement.
"""
from pathlib import Path
import argparse, hashlib, importlib.util, json, math
from decimal import Decimal, localcontext
from fractions import Fraction
import numpy as np
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
TARGET=HERE/"material_field_window_results.json"
def read(p):return json.loads(p.read_text("utf-8-sig"))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def norm(a):return float(np.linalg.norm(a,2))
def fro(a):return float(np.linalg.norm(a))
def gamma(n):
    u=np.finfo(float).eps/2
    return n*u/(1-n*u)
def phases(values,t):
    # Independent high precision range reduction, including the large T.
    pi="3.141592653589793238462643383279502884197169399375105820974944592307816406286"
    out=[]
    with localcontext() as ctx:
        ctx.prec=70;twopi=2*Decimal(pi)
        for val in values:
            x=(-Decimal.from_float(float(val))*Decimal.from_float(float(t)))%twopi
            if x>twopi/2:x-=twopi
            if x<-twopi/2:x+=twopi
            ss=term=x;cc=cterm=Decimal(1)
            for j in range(1,65):
                term*=-x*x/Decimal((2*j)*(2*j+1));ss+=term
                cterm*=-x*x/Decimal((2*j-1)*(2*j));cc+=cterm
            out.append(complex(float(cc),float(ss)))
    return np.array(out)

def material():
    spec=importlib.util.spec_from_file_location("old956",STAGE/"956/native_material_interface.py")
    old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
    m=old.material(1.,.1);h=m["H"];W=m["W"]
    cs=[old.annihilate(i) for i in range(4)];ns=[c.T@c for c in cs]
    sec=np.eye(16)[:,[s for s in range(16) if s.bit_count()==2]]
    Q=sec.T@(ns[0]+ns[1]-ns[2]-ns[3])@sec/2
    Hm=np.kron(h,np.eye(6))+np.kron(np.eye(6),h)+.2*np.kron(Q,Q)
    S=np.kron(Q,np.eye(6))+np.kron(np.eye(6),Q)
    e=np.eye(16)
    code=np.stack([(e[5]-e[6]-e[9]+e[10])/2,
         (2*e[3]+2*e[12]-e[5]-e[6]-e[9]-e[10])/np.sqrt(12)],axis=1)
    local=np.kron(W,np.eye(4))@code
    return m,Q,Hm,S,np.kron(W,W),local

def construction_error(H,Hm,S,n,g,omega):
    # Compare actual stored matrix entries with the exact rational/radical
    # defining Hamiltonian. A row-sum bound controls the symmetric operator error.
    d=H.shape[0];errs=np.zeros(d);s=np.diag(S).astype(int)
    with localcontext() as ctx:
        ctx.prec=70
        gd=Decimal(3)/1000;wd=Decimal(1)/2
        for i in range(36):
            for j in range(36):
                rat=Fraction(float(Hm[i,j])).limit_denominator(100)
                hd=Decimal(rat.numerator)/Decimal(rat.denominator)
                assert abs(float(hd)-Hm[i,j])<1e-15
                for p in range(n):
                    target=hd+(wd*p+gd*gd/wd*int(s[i])**2 if i==j else 0)
                    err=abs(Decimal.from_float(float(H[i*n+p,j*n+p]))-target)
                    errs[i*n+p]+=float(err)+1e-60
            for p in range(n-1):
                target=gd*int(s[i])*Decimal(p+1).sqrt()
                for a,b in ((p,p+1),(p+1,p)):
                    err=abs(Decimal.from_float(float(H[i*n+a,i*n+b]))-target)
                    errs[i*n+a]+=float(err)+1e-60
    return float(max(errs))*(1+1e-10)+1e-55

def calculate(cutoff,m,Q,Hm,S,W0,local,T):
    n=cutoff+1;g=.003;omega=.5
    ann=np.diag(np.sqrt(np.arange(1,n)),1)
    number=np.diag(np.arange(n,dtype=float));X=ann+ann.T
    A=np.kron(Hm,np.eye(n));F=np.kron(np.eye(36),omega*number)
    C=g*np.kron(S,X);D=g*g/omega*np.kron(S@S,np.eye(n));H=A+F+C+D
    val,V=np.linalg.eigh(H)
    # Whole 16-dimensional old active low sector and BOTH photon inputs 0,1.
    Vin=np.kron(W0,np.eye(n)[:,:2]);co=V.T@Vin
    u=np.finfo(float).eps/2;dim=len(val)
    coeff_round=gamma(dim)*(abs(V.T)@abs(Vin))
    coeff_upper=np.linalg.norm(co,axis=1)+np.linalg.norm(coeff_round,axis=1)+1e-13
    errorH=construction_error(H,Hm,S,n,g,omega)
    r=H@V-V*val
    rround=gamma(dim)*(abs(H)@abs(V))+3*u*(abs(H@V)+abs(V*val))
    finite_res=np.linalg.norm(r,axis=0)+np.linalg.norm(rround,axis=0)+errorH*np.linalg.norm(V,axis=0)
    top=V.reshape(36,n,dim)[:,-1,:]
    boundary=g*math.sqrt(n)*(S@top)
    boundary_norm=np.linalg.norm(boundary,axis=0)*(1+gamma(100))+1e-15
    initial_defect=fro(V@co-Vin)
    initial_round=fro(gamma(dim)*(abs(V)@abs(co)))+fro(abs(V)@coeff_round)+1e-11
    # All exact definitions and endpoint evaluation have extra conservative guards.
    eps_initial=initial_defect+initial_round
    eps_boundary=T*float(boundary_norm@coeff_upper)
    eps_finite=T*float(finite_res@coeff_upper)
    eps=eps_initial+eps_boundary+eps_finite+1e-8
    assert errorH<1e-12
    phase=phases(val,T);U=(V*phase)@V.T
    assert float(np.max(abs(phase-np.exp(-1j*T*val))))<2e-10
    phi=np.zeros(n,complex);phi[:2]=1/np.sqrt(2)
    targetB=(np.exp(1j*m["J"]*T)*local[:,0]+local[:,1])/np.sqrt(2)
    targetF=np.zeros(n,complex);targetF[0]=1/np.sqrt(2)
    targetF[1]=1j*phases(np.array([omega]),T)[0]/np.sqrt(2)
    pf=[];pb=[];field_rhos=[];energy_rows=[];means=[]
    for label in (0,1):
        psi=np.kron(np.kron(local[:,label],local@np.array([1.,1.])/np.sqrt(2)),phi)
        block=psi.reshape(6,4,6,4,n).transpose(0,2,4,1,3).reshape(36*n,16)
        evolved=U@block
        tensor=evolved.reshape(6,6,n,4,4).transpose(0,3,1,4,2).reshape(24,24,n)
        rb=np.einsum("abn,adn->bd",tensor,tensor.conj())
        x=tensor.reshape(576,n);rf=x.T@x.conj()
        pb.append(float(np.vdot(targetB,rb@targetB).real))
        pf.append(float(np.vdot(targetF,rf@targetF).real))
        field_rhos.append(rf)
        def expectation(op,state):return float(np.sum(state.conj()*(op@state)).real)
        before=[expectation(op,block) for op in (A,F,C,D)]
        after=[expectation(op,evolved) for op in (A,F,C,D)]
        means.append(dict(label=label,
            dipole=expectation(np.kron(S,np.eye(n)),evolved),
            dipole_square=expectation(np.kron(S@S,np.eye(n)),evolved)))
        assert abs(sum(after)-sum(before))<1e-11
        energy_rows.append(dict(label=label,before=before,after=after,
            changes=[b-a for a,b in zip(before,after)],total_error=abs(sum(after)-sum(before))))
    dfield=float(np.sum(abs(np.linalg.eigvalsh(field_rhos[0]-field_rhos[1])))/2)
    optical=abs(pf[0]-pf[1]);record=pb[1]-pb[0]
    currents=[1j*(H@op-op@H) for op in (A,F,C,D)]
    balance=fro(sum(currents))
    # One gauge-specific parameter source from the SAME Hamiltonian.
    dhdg=np.kron(S,X)+(2*g/omega)*np.kron(S@S,np.eye(n))
    dg=1e-5
    def hg(z):return A+F+z*np.kron(S,X)+z*z/omega*np.kron(S@S,np.eye(n))
    source_error=fro((hg(g+dg)-hg(g-dg))/(2*dg)-dhdg)
    return dict(cutoff=cutoff,dimension=dim,parameters=dict(g=g,omega=omega),
        residual_bound=dict(total=eps,initial=eps_initial,
            infinite_Fock_boundary=eps_boundary,finite_spectral_residual=eps_finite,
            defining_matrix_operator_error=errorH,
            floating_endpoint_guard=1e-8,input_isometry_defect=fro(Vin.T@Vin-np.eye(32))),
        probabilities=dict(receiver=pb,field_phase=pf,receiver_contrast=record,
            fixed_field_effect_contrast=optical,field_trace_distance=dfield,
            infinite_receiver_contrast_lower=record-2*eps,
            infinite_field_contrast_lower=optical-2*eps),
        sources=dict(energy_order=["original_material","photon","linear_coupling","self_polarization"],
            energy_rows=energy_rows,means=means,current_identity_residual=balance,
            same_parameter_source_difference=source_error),
        algebra=dict(charge_internal_commutator=norm(Hm@S-S@Hm),
            low_charge_projection=norm(W0.T@S@W0),
            projected_square_norm=norm(W0.T@S@S@W0),
            self_polarization_operator_norm=norm(D)),
        scope=dict(single_mode_effective_model=True,
            infinite_occupation_bound_uniform_on_declared_input=True,
            endpoint_apparatus_and_photon_phase_reference_are_inputs=True,
            local_propagation_or_full_Maxwell_recovered=False,
            full_microscopic_gauge_matching_proven=False,
            unbounded_energy_observable_error_transported=False,
            full_goal_completed=False))

def run():
    m,Q,Hm,S,W0,local=material()
    old=read(STAGE/"958/capacitive_material_write_results.json")
    T=old["exact_interface"]["controlled_phase_time"]
    rows=[calculate(n,m,Q,Hm,S,W0,local,T) for n in (7,10)]
    best=rows[-1]
    assert best["residual_bound"]["total"]<.002
    assert best["probabilities"]["infinite_receiver_contrast_lower"]>.94
    assert best["probabilities"]["infinite_field_contrast_lower"]>.02
    for row in rows:
        assert row["sources"]["current_identity_residual"]<1e-12
        assert row["sources"]["same_parameter_source_difference"]<1e-8
    delta=max(abs(rows[0]["probabilities"][key]-best["probabilities"][key])
              for key in ("receiver_contrast","fixed_field_effect_contrast","field_trace_distance"))
    assert delta<rows[0]["residual_bound"]["total"]+best["residual_bound"]["total"]
    files=[Path(__file__),STAGE/"956/native_material_interface.py",
        STAGE/"958/capacitive_material_write_results.json",
        STAGE/"960/material_dispersion_bridge_results.json",
        STAGE/"964/material_cosmology_results.json"]
    return dict(round=965,date="2026-10-07",all_scientific_checks_passed=True,
        T=T,old_material_parameters=dict(U=1.,v=.1,kappa=.2),
        old_receiver_effect_and_preparation_retained=True,
        photon_initial_state="(vacuum + one photon)/sqrt(2), with declared phase reference",
        fixed_menu_rows=rows,cutoff_comparison_max_difference=delta,
        adoption="finite same-charge field readout and original record coexist; stop mode optimization",
        scope=dict(new_fundamental_mediator_introduced=False,
            whole_common_model_completed=False,cosmological_photon_matching_completed=False,
            physical_cavity_source_or_optical_apparatus_completed=False),
        references=["https://arxiv.org/html/2211.04241","https://arxiv.org/abs/2006.03191"],
        source_hashes={str(p.relative_to(ROOT)):sha(p) for p in files})

if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--write",action="store_true")
    args=parser.parse_args()
    if args.write:assert not TARGET.exists()
    out=run()
    if args.write:
        with TARGET.open("x",encoding="utf-8") as f:
            json.dump(out,f,ensure_ascii=False,indent=2);f.write("\n")
    else:
        before=read(TARGET);assert before["source_hashes"]==out["source_hashes"]
        assert abs(before["fixed_menu_rows"][-1]["probabilities"]["receiver_contrast"]
                   -out["fixed_menu_rows"][-1]["probabilities"]["receiver_contrast"])<1e-9
    print(json.dumps({k:v for k,v in out.items() if k!="source_hashes"},ensure_ascii=False,indent=2))

