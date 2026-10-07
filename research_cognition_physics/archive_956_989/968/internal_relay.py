"""968: fixed internal spin exchange changes a native relay's participation.
No external switch, graph register, new charge interaction or mean-state H.
"""
from pathlib import Path
from decimal import Decimal,localcontext
import argparse,hashlib,importlib.util,json,math
import numpy as np
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
OUT=HERE/"internal_relay_results.json"
def load(name,path):
    s=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def norm(a):return float(np.linalg.norm(a,2))
def kron(*args):
    r=np.ones((1,1))
    for a in args:r=np.kron(r,a)
    return r
def comm(a,b):return a@b-b@a
def gamma(n):
    u=np.finfo(float).eps/2;return n*u/(1-n*u)

def material():
    m956=load("m956",STAGE/"956/native_material_interface.py")
    m=m956.material(1.,.1)
    cs=[m956.annihilate(i) for i in range(4)]
    ns=[c.T@c for c in cs]
    sector=np.eye(16)[:,[k for k in range(16) if k.bit_count()==2]]
    q=sector.T@(ns[0]+ns[1]-ns[2]-ns[3])@sector/2
    sigma=(np.array([[0,1],[1,0]],complex),
           np.array([[0,-1j],[1j,0]]),np.diag([1.,-1.]))
    s1=[sector.T@sum(.5*p[a,b]*cs[a].T@cs[b] for a in range(2) for b in range(2))@sector
        for p in sigma]
    F=sum(kron(s,p/2,np.eye(2)) for s,p in zip(s1,sigma)).real
    b=np.eye(16);ss=np.array([0.,1.,-1.,0.])/math.sqrt(2)
    code=np.stack([(b[5]-b[6]-b[9]+b[10])/2,
         (2*b[3]+2*b[12]-b[5]-b[6]-b[9]-b[10])/math.sqrt(12)],axis=1)
    trip=kron(m["single"],np.eye(4))@code[:,1]
    B=np.column_stack([np.kron(m["s"],ss),np.kron(m["d"],ss),
                       np.kron(q@m["d"],ss),trip])
    h=np.array([[0.,-.2,0.,0.],[-.2,1.,0.,0.],[0.,0.,1.,0.],[0.,0.,0.,0.]])
    Q=np.zeros((4,4));Q[1,2]=Q[2,1]=1
    F4=np.zeros((4,4));F4[0,3]=F4[3,0]=-math.sqrt(3)/4;F4[3,3]=-.5
    residuals=dict(isometry=norm(B.T@B-np.eye(4)),
        internal=norm(kron(m["H"],np.eye(4))@B-B@h),
        charge=norm(kron(q,np.eye(4))@B-B@Q),
        spin_exchange=norm(F@B-B@F4))
    assert max(residuals.values())<1e-14
    W=np.array([[math.sqrt(1-m["occupancy"]),0],
                [math.sqrt(m["occupancy"]),0],[0.,0.],[0.,1.]])
    actual=kron(m["W"],np.eye(4))@code
    assert norm(B@W-actual)<1e-14
    return m,h,Q,F4,W,residuals

def certificate(H,values,V,time,eta):
    # Independently reconstruct the defining rational/radical H entries.
    # H=Hrat + eta*(-sqrt(3)/4) Fmask. All rational terms are multiples 1/40.
    mask=np.zeros_like(H)
    for a in range(4):
        for c in range(4):
            i=16*a+c;j=16*a+12+c;mask[i,j]=mask[j,i]=1.
    rational=H+eta*math.sqrt(3)/4*mask
    ints=np.rint(40*rational).astype(int)
    assert np.max(abs(rational-ints/40))<2e-15
    with localcontext() as ctx:
        ctx.prec=75;sqrt3=Decimal(3).sqrt();ed=Decimal.from_float(float(eta))
        # eta is exactly 0 or 1/20 in the mathematical model.
        ed=Decimal(0) if eta==0 else Decimal(1)/20
        row=[]
        for i in range(64):
            row.append(sum(float(abs(Decimal.from_float(float(H[i,j]))-
                (Decimal(int(ints[i,j]))/40-ed*sqrt3/4*int(mask[i,j]))))
                +1e-65 for j in range(64)))
    deltaH=max(row)*(1+1e-10)+1e-60
    gram=V.T@V-np.eye(64)
    gram_bound=float(np.linalg.norm(gram))+gamma(64)*float(np.linalg.norm(abs(V).T@abs(V)))+1e-13
    assert gram_bound<1e-10
    Vbound=math.sqrt(1+gram_bound)
    residual=H@V-V*values
    residual_bound=(float(np.linalg.norm(residual))+
        gamma(64)*float(np.linalg.norm(abs(H)@abs(V)))+
        gamma(2)*float(np.linalg.norm(abs(V)*abs(values)))+
        gamma(2)*float(np.linalg.norm(abs(H@V)+abs(V*values)))+
        deltaH*float(np.linalg.norm(V))+1e-14)
    # ||I-VV^t|| has the same singular/eigen deviations as ||I-V^t V||.
    total=gram_bound+abs(time)*residual_bound*Vbound+1e-7
    return dict(defining_H_operator_error=deltaH,gram_defect_bound=gram_bound,
                eigen_residual_bound=residual_bound,endpoint_guard=1e-7,
                total_operator_error_bound=total)

def run():
    m,h,q,f,w,identity=material()
    phase=load("phase965",STAGE/"965/material_field_window.py").phases
    old=json.loads((STAGE/"958/capacitive_material_write_results.json").read_text("utf-8"))
    chi=old["exact_interface"]["conditional_energy_chi"]
    # One predeclared numerical time; its exact stored double is the model's time.
    time=float(m["J"]/chi**2);eta=.05;k=.2
    I=np.eye(4)
    locals_=[kron(h,I,I),kron(I,h,I),kron(I,I,h)]
    edges=[k*kron(q,q,I),k*kron(I,q,q)]
    FB=kron(I,f,I)
    H0=sum(locals_)+sum(edges)
    H=H0+eta*FB
    Ps=np.diag([1.,1.,1.,0.])
    PA=kron(Ps,I,I);PB=kron(I,Ps,I);PC=kron(I,I,Ps)
    assert max(norm(comm(H,PA)),norm(comm(H,PC)))<1e-14
    assert norm(comm(H,PB))>.02
    bare=kron(w,w,w)
    plus=np.array([1.,1.])/math.sqrt(2)
    # Fixed free-material receiving effect (same for all source labels and eta).
    target=(phase(np.array([-m["J"]]),time)[0]*w[:,0]+w[:,1])/math.sqrt(2)
    E=kron(I,I,np.outer(target,target.conj()))
    rows=[]
    for e in (0.,eta):
        He=H0+e*FB;vals,V=np.linalg.eigh(He)
        cert=certificate(He,vals,V,time,e)
        assert cert["total_operator_error_bound"]<1e-4
        U=(V*phase(vals,time))@V.T
        probs=[];detail=[]
        for src in (0,1):
            initial=bare@np.kron(np.kron(np.eye(2)[:,src],np.array([0.,1.])),plus)
            state=U@initial
            prob=float(np.vdot(state,E@state).real)
            active=float(np.vdot(state,PB@state).real)
            norms=norm(U.conj().T@U-np.eye(64))
            probs.append(prob)
            parts=locals_+edges+[e*FB]
            before=[float(np.vdot(initial,A@initial).real) for A in parts]
            after=[float(np.vdot(state,A@state).real) for A in parts]
            assert abs(sum(before)-sum(after))<1e-11
            detail.append(dict(source=src,receiver_probability=prob,
                relay_spin_singlet_probability=active,before_energies=before,
                after_energies=after,unitarity_defect=norms,
                endpoint_records=[float(np.vdot(state,P@state).real) for P in (PA,PC)]))
        contrast=abs(probs[0]-probs[1])
        # Conservative fixed per-state effect budget, including preparation/effect rounding.
        eps=1e-4
        halfwidth=.001
        # Analytic triangle bound: ||h||<1.05, two .2 charge terms,
        # and eta*||s_1.s_3|| <= .05*(3/4), hence ||H||<3.6.
        Hnorm_bound=3.6
        assert norm(He)<Hnorm_bound
        window_lower=contrast-4*eps-4*Hnorm_bound*halfwidth
        currents=sum(1j*comm(He,A) for A in locals_+edges+[e*FB])
        assert norm(currents)<1e-12
        rows.append(dict(eta=e,receiver_probabilities=probs,contrast=contrast,
            certificate=cert,reported_single_state_vector_error=eps,
            reported_fixed_effect_contrast_lower=contrast-4*eps,
            half_window=halfwidth,window_contrast_lower=window_lower,
            analytic_H_norm_bound=Hnorm_bound,
            accepted_nonzero_window=window_lower>.01,
            details=detail,energy_current_balance=norm(currents),
            relay_participation_current_norm=norm(1j*comm(He,PB))))
    # Exact dark relay: Q|1_L>=h|1_L>=0, and eta=0 never leaves it.
    dark=np.eye(4)[:,3];assert np.array_equal(q@dark,np.zeros(4))
    assert np.array_equal(h@dark,np.zeros(4))
    assert rows[0]["contrast"]<1e-8 and max(d["relay_spin_singlet_probability"]
           for d in rows[0]["details"])<1e-12
    # Counterfactual cuts prove the A->C effect depends on BOTH physical edges.
    cuts=[]
    for remove in (0,1):
        Hcut=sum(locals_)+edges[1-remove]+eta*FB
        val,vec=np.linalg.eigh(Hcut);U=(vec*phase(val,time))@vec.T
        probs=[]
        for src in (0,1):
            initial=bare@np.kron(np.kron(np.eye(2)[:,src],np.array([0.,1.])),plus)
            out=U@initial;probs.append(float(np.vdot(out,E@out).real))
        # Analytic tensor factorization proves zero, numerical long-time evaluation
        # is only a check; allow residual phase amplification.
        assert abs(probs[0]-probs[1])<2e-6
        cuts.append(dict(removed_edge=["AB","BC"][remove],probabilities=probs,
                         exact_contrast=0.))
    # Analytical initial acceleration of relay participation.
    initial_relay=w[:,1]
    acceleration=float(np.vdot(initial_relay,-comm(h+eta*f,comm(h+eta*f,Ps))
                              @initial_relay).real)
    assert abs(acceleration-3*eta**2/8)<1e-16
    files=[STAGE/"956/native_material_interface.py",
        STAGE/"958/capacitive_material_write_results.json",
        STAGE/"965/material_field_window.py",
        STAGE/"research_note_967.md",STAGE/"research_note_966.md",
        STAGE.parent/"archive_429_466/research_note_458.md",
        STAGE.parent/"archive_429_466/research_note_459.md",
        HERE/"drafts/internal_relation_decision.md"]
    return dict(round=968,all_scientific_checks_passed=True,
        parameters=dict(U=1.,v=.1,kappa_AB=.2,kappa_BC=.2,kappa_AC=0.,
                        eta=eta,time=time,old_chi=chi,old_J=m["J"],
                        physical_local_dimension=24,exact_invariant_local_dimension=4,
                        physical_total_dimension=13824,exact_invariant_total_dimension=64),
        exact_restriction_checks=identity,restricted_exchange=f.tolist(),
        relay_initial_participation_acceleration=acceleration,
        endpoint_record_projectors_conserved=True,fixed_menu_rows=rows,edge_cut_controls=cuts,
        scope=dict(autonomous_material_participation_rule=True,
            external_time_dependent_coupling_used=False,new_graph_register_added=False,
            geometry_change_derived=False,local_lightcone_derived=False,
            full_QED_matching=False,full_goal_completed=False,
            real_material_lifetime_certified=False,
            signal_window_adopted=rows[1]["accepted_nonzero_window"]),
        source_hashes={str(p.relative_to(ROOT)):sha(p) for p in files},
        references=["https://arxiv.org/abs/quant-ph/0005116"])
def compare(a,b,path=""):
    if isinstance(a,dict):
        assert a.keys()==b.keys(),path
        for k in a:compare(a[k],b[k],path+"/"+k)
    elif isinstance(a,list):
        assert len(a)==len(b),path
        for i,(x,y) in enumerate(zip(a,b)):compare(x,y,path+f"/{i}")
    elif isinstance(a,float):assert math.isclose(a,b,rel_tol=5e-5,abs_tol=2e-7),(path,a,b)
    else:assert a==b,(path,a,b)
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--write",action="store_true");a=p.parse_args()
    result=run()
    if a.write:
        with OUT.open("x",encoding="utf-8") as f:json.dump(result,f,ensure_ascii=False,indent=2);f.write("\n")
    else:compare(result,json.loads(OUT.read_text("utf-8")))
    print(json.dumps(dict(round=968,passed=True,time=result["parameters"]["time"],
      rows=[dict(eta=r["eta"],probabilities=r["receiver_probabilities"],
             contrast=r["contrast"],window_lower=r["window_contrast_lower"],
             adopted=r["accepted_nonzero_window"],bound=r["certificate"]["total_operator_error_bound"])
             for r in result["fixed_menu_rows"]]),ensure_ascii=False))

