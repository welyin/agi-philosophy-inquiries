"""1240: three continuous-clock photon exchange cycles, conditional on g.

Exact Fraction certificates carry the error theorem. Composite Simpson free
segments and event roots diagnose two same-parameter histories; no SciPy.
Default is read-only comparison, --write creates results.json once.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parents[1]
Q=1e-8
END=.82
DISTANCE=.03
EMISSIONS=np.array([0.,.331,.728])


def need(ok,label):
    if not bool(ok): raise AssertionError(label)


def rat(x): return {"fraction":str(x),"decimal":float(x)}


def determinant(e,r):
    return (e[1]-e[0])*(r[2]-r[0])-(e[2]-e[0])*(r[1]-r[0])


def certificates():
    q=F(1,100000000); T=F(41,50); A=F(3,2); d=F(3,100)
    e=[F(0),F(331,1000),F(728,1000)]
    returns=[F(7651,125000),F(6327,15625),F(101981,125000)]
    V=9*q/(1-6*q)
    motion=4*A*T*V+(A+1)*T*V*V
    need(V<F(1,10000000) and 2*q<1-6*q,"finite mass and velocity domain")
    need((1+T)**2<A**3,"a=(1+ht)^(2/3) is below 1.5")
    need(motion<F(4429,10000000000),"analytic motion bound below 4.429e-7")
    trigger_margin=T-T*V*V-e[-1]
    out=(T-e[-1]-T*V*V)/A-d-2*T*V
    back=(T-returns[-1])/A-T*V*V-4*T*V
    order=[e[1]-returns[0]-motion,e[2]-returns[1]-motion,T-returns[2]-motion]
    need(min([trigger_margin,out,back]+order)>0,"first-exit and chronological margins")
    old=json.loads((HERE.parent/"1236/results.json").read_text(encoding="utf-8"))
    previous=json.loads((HERE.parent/"1239/results.json").read_text(encoding="utf-8"))
    geometry=F(previous["group_1"]["per_return_geometry_budget"]["fraction"])
    D=F(old["group_2"]["D"]["fraction"])
    coeff=F(old["group_2"]["linear_error_coefficient"]["fraction"])
    gap=F(old["group_3"]["analytic_Peano_gap_lower"]["fraction"])
    eta=F(3,10000000)
    total=geometry+motion+eta
    flat_error=coeff*(2*total)+8*(2*total)**2
    source_error=10*total+16*total*total
    need(total<F(1,1000000),"full declared six-record budget below old threshold")
    need(flat_error<-D and source_error<gap,"old finite comparisons retain positive margins")
    return dict(status="PASS",certificate_type="analytic, exact rational constants",
                photon_rest_energy=rat(q),mass_lower=rat(1-6*q),velocity_bound=rat(V),
                proper_lag_bound=rat(T*V*V),coordinate_drift_bound=rat(T*V),
                motion_record_bound=rat(motion),trigger_existence_margin=rat(trigger_margin),
                outgoing_first_exit_margin=rat(out),returning_first_exit_margin=rat(back),
                chronological_margins=[rat(x) for x in order],
                previous_homogeneous_geometry_bound=rat(geometry),
                remaining_unimplemented_contract=rat(eta),total_six_record_budget=rat(total),
                arbitrary_flat_two_box_error=rat(flat_error),flat_margin=rat(-D-flat_error),
                both_source_branches_error=rat(source_error),source_margin=rat(gap-source_error),
                same_metric_family_triangle_only=True)


class History:
    def __init__(self,h,panels):
        self.h=float(h); self.panels=panels
        self.t=0.; self.x=np.array([0.,DISTANCE]); self.lag=np.zeros(2)
        self.dm=np.zeros(2); self.Pi=np.zeros(2)
        self.photon=0. # signed comoving photon momentum; energy=abs(photon)/a
        self.work=0.; self.vertices=[]; self.records=[]; self.events=[]
        self.max_energy_vertex=0.; self.max_momentum_vertex=0.; self.max_work_residual=0.
        self.max_velocity=0.; self.max_momentum=0.; self.max_clock_lag=0.

    def scale(self,t): return (1+self.h*np.asarray(t))**(2/3)
    def hubble(self,t): return (2*self.h/3)/(1+self.h*np.asarray(t))

    def conformal_difference(self,t0,t1):
        z0=(1+self.h*t0)**(1/3); z1=(1+self.h*t1)**(1/3)
        return 3*(t1-t0)/(z1*z1+z1*z0+z0*z0)

    def future(self,end):
        """Time-only geodesic integrands, independent composite Simpson rule."""
        if end==self.t: return self.x.copy(),self.lag.copy(),0.
        need(end>=self.t,"forward-only event evolution")
        ts=np.linspace(self.t,end,self.panels+1)
        a=self.scale(ts)
        masses=1+self.dm
        p=self.Pi[:,None]/a[None,:]
        energy=np.sqrt(masses[:,None]**2+p*p)
        dx=p/(a[None,:]*energy)
        # Stable 1 - M/E = p^2 / [E(E+M)], including tiny proper-time lags.
        dl=p*p/(energy*(energy+masses[:,None]))
        power=self.hubble(ts)*(np.sum(p*p/energy,axis=0)+abs(self.photon)/a)
        weights=np.ones(self.panels+1);weights[1:-1:2]=4;weights[2:-1:2]=2
        factor=(end-self.t)/(3*self.panels)
        return (self.x+factor*(dx@weights),self.lag+factor*(dl@weights),
                float(factor*(power@weights)))

    def kinetic(self):
        p=self.Pi/self.scale(self.t); M=1+self.dm
        E=np.sqrt(M*M+p*p)
        return p*p/(E+M)

    def excess_energy(self):
        return float(np.sum(self.dm)+np.sum(self.kinetic())+abs(self.photon)/self.scale(self.t))

    def total_momentum(self): return float((np.sum(self.Pi)+self.photon)/self.scale(self.t))

    def move(self,end):
        self.x,self.lag,dwork=self.future(end);self.work+=dwork;self.t=float(end)
        masses=1+self.dm;p=self.Pi/self.scale(self.t)
        velocity=p/np.sqrt(masses*masses+p*p)
        self.max_velocity=max(self.max_velocity,float(np.max(np.abs(velocity))))
        self.max_momentum=max(self.max_momentum,float(np.max(np.abs(self.Pi))))
        self.max_clock_lag=max(self.max_clock_lag,float(np.max(self.lag)))
        self.max_work_residual=max(self.max_work_residual,abs(self.excess_energy()+self.work))
        need(velocity[0]<=1e-22 and velocity[1]>=-1e-22,"persistent endpoint momentum signs")

    def exp_rapidity(self,which):
        M=1+self.dm[which];p=self.Pi[which]/self.scale(self.t)
        kinetic=p*p/(np.sqrt(M*M+p*p)+M)
        return 1+(kinetic+p)/M

    def root(self,function):
        lo=self.t;hi=END
        fl=function(lo);fh=function(hi)
        need(fl<=1e-15 and fh>0,"event bracket with positive finite upper margin")
        if abs(fl)<1e-17: return lo
        for _ in range(54):
            mid=(lo+hi)/2
            if function(mid)<0:lo=mid
            else:hi=mid
        return (lo+hi)/2

    def vertex(self,kind,cycle):
        before_energy=self.excess_energy();before_momentum=self.total_momentum()
        a=float(self.scale(self.t));mass_before=(1+self.dm).copy()
        incoming=abs(self.photon)/a
        loss=0.;gain=0.;outgoing=0.
        if kind=="emit":
            M=1+self.dm[0]
            need(2*Q<M,"positive emitted mass shell")
            loss=2*Q/(1+np.sqrt(1-2*Q/M))
            outgoing=Q*self.exp_rapidity(0)
            self.dm[0]-=loss;self.Pi[0]-=a*outgoing;self.photon=a*outgoing
        elif kind=="reflect":
            M=1+self.dm[1];em=self.exp_rapidity(1)
            rest_in=incoming/em
            outgoing=(rest_in/(1+2*rest_in/M))/em
            self.Pi[1]+=a*(incoming+outgoing);self.photon=-a*outgoing
        elif kind=="absorb":
            M=1+self.dm[0];rest_in=incoming*self.exp_rapidity(0)
            gain=2*rest_in/(1+np.sqrt(1+2*rest_in/M))
            self.dm[0]+=gain;self.Pi[0]-=a*incoming;self.photon=0.
        else:raise ValueError(kind)
        en=self.excess_energy()-before_energy
        mom=self.total_momentum()-before_momentum
        self.max_energy_vertex=max(self.max_energy_vertex,abs(en))
        self.max_momentum_vertex=max(self.max_momentum_vertex,abs(mom))
        self.max_work_residual=max(self.max_work_residual,abs(self.excess_energy()+self.work))
        masses=1+self.dm
        p=self.Pi/a
        velocity=p/np.sqrt(masses*masses+p*p)
        self.max_velocity=max(self.max_velocity,float(np.max(np.abs(velocity))))
        self.max_momentum=max(self.max_momentum,float(np.max(np.abs(self.Pi))))
        need(velocity[0]<=1e-22 and velocity[1]>=-1e-22,"post-vertex momentum signs")
        need(min(masses)>=1-6*Q,"finite cumulative mass lower bound")
        need(max(incoming,outgoing)<=Q*(1+1e-14),"no photon replenishment; redshift and recoil retained")
        need(loss<=2*Q and gain>=0,"fuel loss and absorption storage correctly signed")
        self.vertices.append(dict(cycle=cycle,kind=kind,time=self.t,
                                  A_proper_time=float(self.t-self.lag[0]),
                                  masses_before=mass_before.tolist(),masses_after=masses.tolist(),
                                  stable_mass_loss=float(loss),stable_mass_gain=float(gain),
                                  normal_photon_in=float(incoming),normal_photon_out=float(outgoing),
                                  energy_residual=float(en),momentum_residual=float(mom),
                                  A_position=float(self.x[0]),B_position=float(self.x[1]),
                                  comoving_momenta=self.Pi.tolist(),clock_lags=self.lag.tolist()))
        self.events.append(self.t)

    def run(self):
        for cycle,e in enumerate(EMISSIONS):
            need(self.t-self.lag[0]<=e+1e-15,"previous return precedes next proper trigger")
            trigger=self.root(lambda t:t-self.future(t)[1][0]-e)
            self.move(trigger);self.vertex("emit",cycle)
            emit_t=self.t;emit_x=float(self.x[0])
            reflection=self.root(lambda t:self.conformal_difference(emit_t,t)
                                  -(self.future(t)[0][1]-emit_x))
            self.move(reflection);self.vertex("reflect",cycle)
            reflect_t=self.t;reflect_x=float(self.x[1])
            receive=self.root(lambda t:self.conformal_difference(reflect_t,t)
                               -(reflect_x-self.future(t)[0][0]))
            self.move(receive);self.vertex("absorb",cycle)
            self.records.append(float(self.t-self.lag[0]))
        V=9*Q/(1-6*Q)
        need(self.max_velocity<=V and self.max_momentum<=9*Q,"same-history motion domain")
        need(self.max_clock_lag<=END*V*V,"same-history proper lag domain")
        need(self.max_energy_vertex<2e-22 and self.max_momentum_vertex<2e-22,
             "stable local complete four-momentum identities")
        need(self.max_work_residual<2e-21,"whole branch energy plus redshift work")
        need(np.all(np.diff(self.events)>0),"all nine interactions are chronologically distinct")
        return self


def original_returns(h):
    z=1+h*EMISSIONS;d=DISTANCE
    return EMISSIONS+2*d*z**(2/3)+(4/3)*h*d*d*z**(1/3)+(8/27)*h*h*d**3


def three_histories():
    rows=[]
    for h in (1.,np.sqrt(.8)):
        fine=History(h,256).run();coarse=History(h,128).run()
        reference=original_returns(h)
        records=np.array(fine.records);shifts=records-reference
        need(min(shifts)>0 and max(shifts)<4.429e-7,"actual finite recoil shifts inside certified bound")
        difference=max(abs(np.array(fine.records)-np.array(coarse.records)))
        need(difference<3e-14,"independent Simpson resolution diagnostic")
        need(fine.records[0]<EMISSIONS[1] and fine.records[1]<EMISSIONS[2]
             and fine.records[2]<END,"same-history menu completed without reset")
        rows.append(dict(h=h,panels=[128,256],vertices=fine.vertices,
                         proper_returns=fine.records,original_returns=reference.tolist(),
                         actual_motion_shifts=shifts.tolist(),return_grid_difference=float(difference),
                         final_A_mass=float(1+fine.dm[0]),final_B_mass=float(1+fine.dm[1]),
                         final_A_mass_increment=float(fine.dm[0]),
                         final_comoving_momenta=fine.Pi.tolist(),final_clock_lags=fine.lag.tolist(),
                         final_positions=fine.x.tolist(),maximum_velocity=float(fine.max_velocity),
                         maximum_comoving_momentum=float(fine.max_momentum),
                         maximum_vertex_energy_residual=float(fine.max_energy_vertex),
                         maximum_vertex_momentum_residual=float(fine.max_momentum_vertex),
                         integrated_redshift_work=float(fine.work),
                         final_branch_energy_excess=float(fine.excess_energy()),
                         maximum_full_branch_work_residual=float(fine.max_work_residual),
                         actual_D=float(determinant(EMISSIONS,records)),
                         reference_D=float(determinant(EMISSIONS,reference))))
    need(rows[1]["actual_D"]-rows[0]["actual_D"]>3.3e-5,"source contrast survives real cumulative history")
    return dict(status="PASS",interpretation="given-g classical event diagnostics, no self-consistent geometry claim",cases=rows)


def deletions(histories):
    rows=[]
    for case in histories["cases"]:
        reflection=next(v for v in case["vertices"] if v["kind"]=="reflect")
        frozen_p_residual=reflection["normal_photon_in"]+reflection["normal_photon_out"]
        need(frozen_p_residual>1.9e-8,"freezing mirror violates local spatial momentum")
        rr=np.array(case["proper_returns"])
        D=determinant(EMISSIONS,rr)
        kick=-D/EMISSIONS[1]
        kicked=rr.copy();kicked[-1]+=kick
        erased=determinant(EMISSIONS,kicked)
        need(abs(erased)<1e-16 and kick>1e-3,"last absorption phase freedom can erase negative witness")
        need(abs(D-case["reference_D"])>1e-10,"freezing all endpoints changes the actual three-event record")
        # Vertex boundary terms for any test covector cancel by the same local P.
        # Finite covectors only diagnose this algebra; distribution proof is analytic.
        ward=max(abs((j+1)*v["energy_residual"]+(2-j)*v["momentum_residual"])
                 for j,v in enumerate(case["vertices"]))
        need(ward<2e-21,"test-covector boundary residuals")
        rows.append(dict(h=case["h"],frozen_mirror_momentum_residual=frozen_p_residual,
                         actual_minus_frozen_D=float(D-case["reference_D"]),
                         last_absorption_phase_kick=float(kick),D_after_kick=float(erased),
                         earlier_emissions_and_returns_unchanged_by_final_kick=True,
                         vertex_test_covector_residual=float(ward)))
    return dict(status="PASS",cases=rows,
                phase_counterexample_reuses_1231_not_a_new_phase_theorem=True,
                no_claim_kick_control_is_built=True)


def calculate():
    first=certificates();second=three_histories();third=deletions(second)
    result=dict(status="PASS",science_groups=3,group_1=first,group_2=second,group_3=third,
                scope=dict(CT="candidate conjecture, not adopted",
                           instantaneous_classical_local_vertex_is_input=True,
                           continuous_theta_and_fixed_I_are_extra_contract=True,
                           theorem_is_on_declared_metric_family=True,
                           branch_distributional_stress_conserved=True,
                           branch_includes_unbuilt_control_and_readout=False,
                           baseline_localized_mass_geometry_matched=False,
                           directional_source_geometry_feedback_obtained=False,
                           actual_finite_contact_or_quantum_packet_obtained=False,
                           arbitrary_unknown_I_and_quantum_reference_obtained=False,
                           remaining_eta_is_unimplemented=True,
                           full_GR_proved=False))
    files=[HERE/name for name in ("admission.md","proof.md","check.py","NEXT.md")]
    files += [HERE.parent/"1236"/name for name in ("proof.md","results.json")]
    files += [HERE.parent/"1239"/name for name in ("proof.md","results.json")]
    files += [RESEARCH/"archive_956_989/research_note_971.md",
              RESEARCH/"archive_1161_1200/research_note_1162.md",
              RESEARCH/"archive_1161_1200/research_note_1169.md",
              RESEARCH/"archive_1219_1233/1231/proof.md"]
    result["source_sha256"]={str(p.relative_to(RESEARCH)).replace("\\","/"):
                             hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    return result


def main():
    parser=argparse.ArgumentParser();parser.add_argument("--write",action="store_true")
    args=parser.parse_args();result=calculate();path=HERE/"results.json"
    if args.write:
        if path.exists():raise FileExistsError("results.json already exists; default is read-only")
        path.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
        mode="CREATED"
    else:
        need(result==json.loads(path.read_text(encoding="utf-8")),"read-only result and source binding")
        mode="READ_ONLY_MATCH"
    print(json.dumps(dict(status="PASS",mode=mode,science_groups=3,
                          motion_bound=result["group_1"]["motion_record_bound"]["decimal"],
                          total_budget=result["group_1"]["total_six_record_budget"]["decimal"],
                          source_margin=result["group_1"]["source_margin"]["decimal"])))


if __name__=="__main__":main()
