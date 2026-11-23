"""Transport the same charged and neutral matter to the gauge crossing or the first new pole."""
import argparse
import hashlib
import json
import math
import runpy
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
CFG = json.loads((HERE/"inputs.json").read_bytes())
LOOP = 16*math.pi**2


def read(key):
    return json.loads((HERE/CFG[key]).read_bytes())


def sv(a):
    u,s,_=np.linalg.svd(a,full_matrices=False)
    return u[:,::-1],s[::-1]


def packed(y,k,flavour_integral=0):
    return np.r_[y.ravel(),k.ravel(),flavour_integral]


def unpack(z):
    return z[:48].reshape(4,3,4),z[48:57].reshape(3,3)


class Running:
    def __init__(self):
        low=runpy.run_path(str(HERE/CFG["matrix_flow"]))
        self.flow,self.rep=low["charged_background"](.1)
        self.start=CFG["start_GeV"]
        self.g0sq=self.flow["segments"](.1)[-1][3][-1]
        self.b=np.array([33/5,1.,-3.])+np.sum([bb for _,bb in self.flow["physics"]["thresholds"]],axis=0)
        np.testing.assert_allclose(self.b,[43/5,3,-1])
        self.stop=read("gauge_spectrum")["shared_witness"]["U_GeV"]
        dispersion=np.ptp(4*math.pi/self.g2(math.log(self.stop/self.start)))
        assert dispersion<2e-8
        dec=low["unpack"]
        ueh=dec(self.rep["charged_left_10PeV"])
        self.ports=ueh.conj().T@dec(self.rep["neutral_columns_10PeV"])
        w,_=low["congruence"](self.flow,self.rep,.1)
        q=ueh.conj().T@w
        self.lower=q*np.conj(np.diag(q)/abs(np.diag(q)))[None,:]
        self.ckm=self.flow["low_ckm"](self.rep["gamma_10PeV_deg"],self.rep["CKM_angles_at_10PeV_rad"])
        self.up=np.array(self.rep["up_texture_diag_10PeV"])
        self.down=np.array(self.rep["Y_down_at_10PeV"])
        self.lepton=np.array(self.rep["Y_charged_at_10PeV"])

    def g2(self,log_ratio):
        return 1/(1/self.g0sq-2*self.b*log_ratio/LOOP)

    def initial(self):
        y=np.zeros((4,3,4),complex)
        y[2,:,:3]=np.diag(self.lepton)
        y[3]=self.ports
        return y

    def pure_seed(self,witness):
        y=self.initial()
        d=self.lepton[-1]*(np.array(witness["symmetric_down_in_tau_units_real"])
                            +1j*np.array(witness["symmetric_down_in_tau_units_imag"]))
        # The R1798 witness D is holomorphic, while the matrix-flow left rows use D*.
        y[1,:,:3]=d.conj()
        ud,sd=sv(y[1])
        uu=ud@self.ckm.conj().T
        y[0,:,:3]=(uu*self.up)@uu.T
        k=y[2,:,:3]-y[1,:,:3]
        np.testing.assert_allclose(sd,self.down,rtol=1e-8)
        return packed(y,k)

    def mixed_seed(self,x):
        mixed=runpy.run_path(str(HERE/CFG["mixed_matching"]))
        data=(self.up,self.down,self.lepton,self.ckm,self.lower)
        u,d,e,f,_,ue,_,_=mixed["matrices"](x,data)
        a=(np.exp(1j*x[6])*self.lepton[-1]
           /np.linalg.norm(d/self.down[-1]+complex(x[4],x[5])*u/self.up[-1],2)/self.down[-1])
        b=a*complex(x[4],x[5])*self.down[-1]/self.up[-1]
        np.testing.assert_allclose(e,a*d+b*u,atol=1e-14)
        rF=(-3-a)/b
        yn=u-4*rF*f
        ki=yn@np.linalg.solve(f,yn.T)
        zeta=complex(x[7],x[8])
        k=f/np.linalg.norm(f)+zeta*ki/np.linalg.norm(ki)
        t=ue@self.lower
        unit=mixed["neutral"](t.conj().T@k@t.conj())[1]["common_light_scale_eV"]
        tan=self.flow["TAN"]
        vu_over_root2=246.22/math.sqrt(2)*tan/math.hypot(1,tan)
        poles=sv(f)[1]*vu_over_root2**2*np.linalg.norm(ki)/(unit*1e-9*abs(zeta))
        y=self.initial()
        for i,m in enumerate((u,d,e)):
            y[i,:,:3]=ue.conj().T@m@ue.conj()
        k=ue.conj().T@k@ue.conj()
        majorana=-vu_over_root2**2*np.linalg.norm(ki)/(unit*1e-9*zeta)*f
        v,m=sv(majorana)
        # Phase-correct the SVD to get M=V diag(m) V^T. The right columns then give the actual seesaw.
        phases=np.diag(v.conj().T@majorana@v.conj())
        v=v*np.exp(.5j*np.angle(phases))
        np.testing.assert_allclose((v*m)@v.T,majorana,rtol=1e-12)
        columns=ue.conj().T@yn@v.conj()
        np.testing.assert_allclose(-vu_over_root2**2/(unit*1e-9)*(columns/m)@columns.T,
                                   k-ue.conj().T@f@ue.conj()/np.linalg.norm(f),atol=2e-13)
        return packed(y,k),poles,columns,ue.conj().T@f@ue.conj()/np.linalg.norm(f)

    def rhs(self,z,x):
        y,k=unpack(z)
        h=y@y.conj().transpose(0,2,1)
        u,d,e,n=h
        # beta above the original thresholds; k denotes kappa*, with its common scale removed.
        p=e+n
        wb=sv(y[1])[0][:,-1]
        wt=sv(y[2])[0][:,-1]
        yuk=(wb.conj()@(3*d+u)@wb-wt.conj()@(3*e+n)@wt).real
        return packed(self.flow["beta"](y,self.g2(x),True),p@k+k@p.T,yuk)/LOOP

    def step(self,z,x,h):
        a=self.rhs(z,x)
        b=self.rhs(z+h*a/2,x+h/2)
        c=self.rhs(z+h*b/2,x+h/2)
        d=self.rhs(z+h*c,x+h)
        return z+h*(a+2*b+2*c+d)/6


def observables(z):
    y,k=unpack(z)
    _,sd=sv(y[1])
    ue,se=sv(y[2])
    c=ue.conj().T@k@ue.conj()
    nu,ms=sv(c)
    p=abs(nu)**2
    tau=(se[-1]-se[-2])*math.sqrt(max(0,1-np.dot(p[-1],p[-1])))
    full=math.sqrt(max(0,se@se-np.linalg.norm(se@p)**2))
    # Min ||E+a K|| over a and free charged-row phases, for the actual Majorana phases.
    fixed=math.sqrt(max(0,se@se-(se@abs(np.diag(c)))**2/np.linalg.norm(c)**2))
    phase_free=math.sqrt(max(0,se@se-(se@p@ms)**2/(ms@ms)))
    dn=np.linalg.norm(sd)
    assert fixed+1e-12>=phase_free and full+1e-12>=tau
    return {"down_singulars":sd,"charged_singulars":se,"tau_probabilities":p[-1],
            "full_PMNS_moduli":abs(nu),"neutral_singular_ratios":ms/ms[-1],
            "bounds_over_actual_down_norm":{"tau":tau/dn,"three_flavours":full/dn,
                       "mass_direction_phase_free":phase_free/dn,"fixed_Majorana":fixed/dn},
            "fixed_direction_missing_Frobenius":max(0,fixed-dn),
            "fixed_direction_delta_b_lower":max(0,math.sqrt(max(0,fixed**2-sd[0]**2-sd[1]**2))/sd[-1]-1)}


def run_seed(run,seed,end,h,locate=True):
    hi=math.log(end/run.start)
    n=max(1,math.ceil(hi/h))
    dt=hi/n
    z=seed.copy()
    first=observables(seed)
    crossings={}
    last=first
    for j in range(n):
        new=run.step(z,j*dt,dt)
        this=observables(new)
        if locate:
            for key,ratio in this["bounds_over_actual_down_norm"].items():
                prev=last["bounds_over_actual_down_norm"][key]
                if key not in crossings and prev<=1<ratio:
                    lo,up=0.,dt
                    for _ in range(28):
                        mid=(lo+up)/2
                        probe=run.step(z,j*dt,mid)
                        if observables(probe)["bounds_over_actual_down_norm"][key]<1:
                            lo=mid
                        else:
                            up=mid
                    crossings[key]=run.start*math.exp(j*dt+(lo+up)/2)
        z,last=new,this
    gi=np.log(run.g2(hi)/run.g0sq)/(2*run.b)
    gauge=4/3*gi[0]-16/3*gi[2]
    actual=math.log((last["down_singulars"][-1]/last["charged_singulars"][-1])/
                    (first["down_singulars"][-1]/first["charged_singulars"][-1]))
    assert abs(actual-z[-1].real-gauge)<8e-10
    return {"start":first,"stop":last,"end_GeV":end,"sign_changing_bracket_roots_GeV":crossings,
            "log_bottom_tau_change":{"gauge":gauge,"yukawa":z[-1].real,"observed":actual},
            "max_abs_PMNS_modulus_drift":np.max(abs(last["full_PMNS_moduli"]-first["full_PMNS_moduli"]))}


def mixed_uv(run,seed,poles,columns,ii,h):
    """A holomorphic light transport specifies each new column at its own tree mass threshold."""
    y0,total=unpack(seed)
    y=np.zeros((4,3,7),complex)
    y[:,:,:4]=y0
    v=np.r_[y.ravel(),np.eye(3).ravel()]
    comparisons=[]
    jump=0.

    def state(z):
        return z[:84].reshape(4,3,7),z[84:].reshape(3,3)

    def rhs(z,x):
        y,gain=state(z)
        g=run.g2(x)
        u,e,n=(y[i]@y[i].conj().T for i in (0,2,3))
        external=e+n+(3*np.trace(u).real+np.trace(n).real-.6*g[0]-3*g[1])*np.eye(3)
        return np.r_[run.flow["beta"](y,g,True).ravel(),(external@gain).ravel()]/LOOP

    def data(z):
        y,gain=state(z)
        current_ii=gain@ii@gain.T
        current_total=gain@total@gain.T
        d,e,u=(y[j,:,:3] for j in [1,2,0])
        def view(k):
            # Only three charged matrices and k enter observables; no active neutral columns are dropped in the RG.
            return observables(packed(y[:,:,:4],k))
        def relative_mismatch(a,b,k):
            # Minimize over arbitrary complex coefficient of k and an arbitrary relative global phase of a,b.
            khat=k.ravel()/np.linalg.norm(k)
            a=a.ravel()-khat*np.vdot(khat,a.ravel())
            b=b.ravel()-khat*np.vdot(khat,b.ravel())
            phase=np.vdot(a,b)
            phase=phase/abs(phase) if abs(phase)>0 else 1.
            return np.linalg.norm(phase*a-b)/np.linalg.norm(d)
        pencil=np.column_stack((d.ravel(),u.ravel()))
        fit=np.linalg.lstsq(pencil,e.ravel(),rcond=None)[0]
        return {"pure_II_component":view(current_ii),"full_mass":view(current_total),
                "fixed_tag_D_minus_E_versus_II_relative_Frobenius":relative_mismatch(d,e,current_ii),
                "fixed_tag_E_outside_complex_span_D_and_U_fraction":np.linalg.norm(pencil@fit-e.ravel())/np.linalg.norm(e),
                "cancellation_Frobenius":(np.linalg.norm(current_ii)+np.linalg.norm(current_total-current_ii))/np.linalg.norm(current_total),
                "light_transport_condition_number":np.linalg.cond(gain)}

    initial=data(v)
    Q=run.start
    for j,end in enumerate([*poles,run.stop]):
        if not end>Q:
            raise ValueError("Registered poles must lie in the running hierarchy")
        steps=math.ceil(math.log(end/Q)/h)
        dt=math.log(end/Q)/steps
        for i in range(steps):
            x=math.log(Q/run.start)+i*dt
            a=rhs(v,x)
            b=rhs(v+dt*a/2,x+dt/2)
            c=rhs(v+dt*b/2,x+dt/2)
            d=rhs(v+dt*c,x+dt)
            v+=dt*(a+2*b+2*c+d)/6
        if j<len(poles):
            y,gain=state(v)
            new_column=gain@columns[:,j]
            # Transfer this rank-one mass from the contact operator to explicit exchange; all quantities dimensionless.
            source0=np.outer(columns[:,j],columns[:,j])/poles[j]
            now=np.outer(new_column,new_column)/poles[j]
            jump=max(jump,np.linalg.norm(now-gain@source0@gain.T)/np.linalg.norm(now))
            y[3,:,4+j]=new_column
            comparisons.append({"scale_GeV":end,"active_columns_after_threshold":5+j,
                                "new_column_norm":np.linalg.norm(new_column)})
        Q=end
    return {"initial":initial,"at_gauge_crossing":data(v),"steps":comparisons,
            "maximum_relative_tree_source_jump":jump}


def low_to_high_certificate(run,seed,h):
    """Prove a flavour- and absolute-mass-uniform high-scale bound using actual low-scale oscillations."""
    start,_=unpack(seed)
    v=np.r_[start.ravel(),np.eye(3).ravel()]
    steps=math.ceil(math.log(run.stop/run.start)/h)
    dt=math.log(run.stop/run.start)/steps
    def rhs(z,x):
        y=z[:48].reshape(4,3,4)
        gain=z[48:].reshape(3,3)
        e,n=(y[i]@y[i].conj().T for i in (2,3))
        return np.r_[run.flow["beta"](y,run.g2(x),True).ravel(),((e+n)@gain).ravel()]/LOOP
    for j in range(steps):
        x=j*dt
        a=rhs(v,x)
        b=rhs(v+dt*a/2,x+dt/2)
        c=rhs(v+dt*b/2,x+dt/2)
        d=rhs(v+dt*c,x+dt)
        v+=dt*(a+2*b+2*c+d)/6
    y=v[:48].reshape(4,3,4)
    ue,e=sv(y[2])
    down_norm=np.linalg.norm(y[1])
    S=run.lower.conj().T@np.linalg.solve(v[48:].reshape(3,3),ue)
    gram=S.conj().T@S
    e_min_sq=np.sum((e*gram.diagonal().real)**2)-sum(
        2*e[i]*e[j]*abs(gram[i,j])**2 for i in range(3) for j in range(i))
    W=(S*e)@S.conj().T
    w=W.diagonal().real
    off=np.linalg.norm(W-np.diag(w))
    H2=np.linalg.norm(S,2)**2
    hist=runpy.run_path(str(HERE/"../1453/check.py"))
    cfg=hist["read"]("inputs.json")
    osc=hist["read"](cfg["oscillation"])
    juno=hist["read"](cfg["juno"])["juno_207_days"]
    box=hist["box"]
    def tiles(bounds,n):
        grid=np.linspace(*bounds,n+1)
        return np.stack((grid[:-1],grid[1:]))
    reactor=[(1-math.sqrt(1-f))/2 for f in box(osc["DayaBay_amplitude"])]
    out={}
    for order in ("NO","IO"):
        solar=box(cfg["NO_solar_sin2"] if order=="NO" else osc[order]["solar_sine_squared"])
        atm=box(osc[order]["atmospheric_sine_squared"])
        a,A=tiles(solar,16)[:,:,None,None]
        b,B=tiles(reactor,10)[:,None,:,None]
        c,C=tiles(atm,32)[:,None,None,:]
        # Since w_e < w_mu, a lower electron probability gives an upper charged expectation.
        pe=[(1-A)*(1-B),a*(1-B),b]
        u1,U1=np.sqrt(a*c),np.sqrt(A*C)
        v1,V1=np.sqrt((1-A)*(1-C)*b),np.sqrt((1-a)*(1-c)*B)
        u2,U2=np.sqrt((1-A)*c),np.sqrt((1-a)*C)
        v2,V2=np.sqrt(a*(1-C)*b),np.sqrt(A*(1-c)*B)
        plus=lambda u,U,v,V:(U+V)**2
        minus=lambda u,U,v,V:np.maximum(abs(u-V),abs(U-v))**2
        qmax_u=[]
        for op1,op2 in ((plus,minus),(minus,plus)):
            pt=[op1(u1,U1,v1,V1),op2(u2,U2,v2,V2),(1-c)*(1-b)]
            q=[w[1]+(w[0]-w[1])*pe[j]+(w[2]-w[1])*pt[j] for j in range(3)]
            qmax_u.append(np.stack(np.broadcast_arrays(*q),axis=-1).reshape(-1,3))
        coefficients=np.concatenate(qmax_u)
        ss=box(juno["solar_21"])
        aa=box(juno["normal_atmospheric_31" if order=="NO" else "inverted_abs_31"])
        r,R=ss[0]/aa[1],ss[1]/aa[0]
        t,T=tiles((0.,1.),512)
        ones=np.ones_like(t)
        if order=="NO":
            low=np.array([t,np.sqrt(r+(1-r)*t*t),ones])
            high=np.array([T,np.sqrt(R+(1-R)*T*T),ones])
        else:
            low=np.array([ones,np.sqrt(1+r*(1-T*T)),t])
            high=np.array([ones,np.sqrt(1+R*(1-t*t)),T])
        upper_q=high/np.linalg.norm(low,axis=0)
        overlap=max(np.max(coefficients@slab) for slab in np.array_split(upper_q,16,axis=1))+off+1e-13
        required=math.sqrt(max(0,e_min_sq-overlap**2))/H2
        out[order]={"minimum_full_down_norm":required,"available_down_norm":down_norm,
                    "minimum_fractional_correction":required/down_norm-1,"threshold_error_Frobenius":max(0,required-down_norm)}
    return {"high_to_MZ_congruence_singular_values_up_to_scalar":sv(S)[1],
            "charged_offdiagonal_over_tau_budget":off/e[-1],"joint_angular_boxes":16*10*32,
            "absolute_mass_ratio_boxes":512,"ordering_results":out}


def calculate(step):
    run=Running()
    pure=read("two_type_II_examples")
    out={}
    for name in ("standard_z_symmetric_oscillation_calibrated_witness","another_standard_z_symmetric_oscillation_calibrated_witness"):
        seed=run.pure_seed(pure[name])
        out[name]=run_seed(run,seed,run.stop,step)
    global_bound=low_to_high_certificate(run,seed,step)
    mixed={}
    for label,x in read("mixed_sources")["mixed_type_I_II"].items():
        seed,poles,columns,ii=run.mixed_seed(x)
        # Leading estimates at the original effective matching, not self-consistently evolved pole masses.
        mixed[label]={"tree_seesaw_scales_GeV":poles,"up_to_first_estimated_threshold":run_seed(run,seed,min(poles[0],run.stop),step,False),
                      "sequential_tree_matching_model":mixed_uv(run,seed,poles,columns,ii,step)}
    return {"actual_threshold_gauge_crossing_GeV":run.stop,"gauge_slopes_above_10PeV":run.b,
            "inverse_alphas_at_crossing":4*math.pi/run.g2(math.log(run.stop/run.start)),
            "whole_low_energy_oscillation_box":global_bound,
            "two_pure_II_histories":out,"mixed_controls":mixed}


def clean(a):
    if isinstance(a,dict):
        return {k:clean(x) for k,x in a.items()}
    if isinstance(a,(list,tuple,np.ndarray)):
        return [clean(x) for x in a]
    if isinstance(a,(float,np.floating)):
        return float(format(a,".10g"))
    return a


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    main=calculate(CFG["step_log_mu"])
    other=calculate(CFG["step_log_mu"]/2)
    errors=[]
    for key,item in main["two_pure_II_histories"].items():
        fine=other["two_pure_II_histories"][key]
        errors.append(abs(item["stop"]["fixed_direction_delta_b_lower"]-fine["stop"]["fixed_direction_delta_b_lower"]))
        for k,c in item["sign_changing_bracket_roots_GeV"].items():
            errors.append(abs(c/fine["sign_changing_bracket_roots_GeV"][k]-1))
    for key,item in main["mixed_controls"].items():
        a=item["sequential_tree_matching_model"]["at_gauge_crossing"]
        b=other["mixed_controls"][key]["sequential_tree_matching_model"]["at_gauge_crossing"]
        errors.append(abs(a["fixed_tag_D_minus_E_versus_II_relative_Frobenius"]-
                          b["fixed_tag_D_minus_E_versus_II_relative_Frobenius"]))
    main["max_absolute_delta_b_or_relative_crossing_step_error"]=max(errors)
    assert max(errors)<2e-7
    data=json.dumps(clean(main),indent=2,allow_nan=False)+"\n"
    target=HERE/"results.json"
    if args.check:
        assert json.loads(target.read_bytes())==json.loads(data)
    else:
        target.write_text(data,encoding="utf-8")
    print("1800",hashlib.sha256(data.encode()).hexdigest()[:12])
