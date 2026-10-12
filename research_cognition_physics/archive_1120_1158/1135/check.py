"""1135: bare support / original-clock finite boost endpoint.

Default checks are read-only; --write-results writes this results.json only.
Finite noncommuting masses and continuous one-dimensional wave packets are
mechanism diagnostics, not simulations of the actual normal Fock source.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FROZEN = {
    "archive_1120_/1134/proof.md": "4b75693a9efeb7c2747152763da6fb3bed2f1e9b25418e9827cdd615f49c5666",
    "archive_1120_/1134/results.json": "e5fb6abfac5ade09cbdb02f3a3719b2114464ebb6f97c4289e537e9175403894",
    "archive_1120_/1132/proof.md": "f05daf49972a0b7fe13f5d593a032998bf1c9008182c051f2e12bb6ea958df63",
    "archive_1120_/1132/results.json": "95dcbd80a2a9e5d3d7ba03e52f1d1aa7e931606c9735f0359623f7b4f969430c",
    "archive_1106_1119/1106/proof.md": "3cf66b629a49e00a6503b68e516f67bdd5b9fecd258447e90156b15cf87fb8cb",
    "archive_1106_1119/1106/results.json": "1fa96fcf58ce22ccf529aebf0fd8d0eb8dd4edb1bdb449998e0ccaec2bbea7e2",
    "archive_702_741/research_note_720.md": "d5a25eddbaa9ea85ac305ba7b6bea68aa0d36717a53e68227ddfea6dda65d75b",
}


def opnorm(a):
    return float(np.linalg.norm(a, 2))


def sqrtpos(a):
    e, v = np.linalg.eigh(a)
    assert np.min(e) > 0
    return (v*np.sqrt(e))@v.conj().T


def exp_h(a,t):
    e,v=np.linalg.eigh(a)
    return (v*np.exp(-1j*t*e))@v.conj().T


def sources():
    hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in FROZEN}
    assert hashes==FROZEN, hashes
    old=json.loads((ROOT/"archive_1106_1119/1106/results.json").read_text(encoding="utf-8"))
    material=old["group3_independently_defined_local_material"]
    assert material["local_dimension"]==288 and material["data_dimension"]==72
    assert material["launch_spectrum"]==[0.,1.,2.,3.]
    assert material["full_72_data_endpoint_residual"]<1e-10
    assert material["launch_probabilities"]==[.125,.375,.375,.125]
    clockrow=old["group5_clock_record_and_real_receiver"]["clock_record_rows"][1]
    assert abs(clockrow["M_probability"]-.125)<1e-12
    assert clockrow["all_input_effect_residual"]<1e-12
    return dict(source_sha256=hashes,old_dimension=288,old_mass_block_multiplicity=72,
        old_gamma_residual=material["full_72_data_endpoint_residual"],
        original_Mr_all_input_identity_residual=clockrow["all_input_effect_residual"],
        old_large_device_reexecuted=False,H_BT_changed=False,
        new_contract="M1135: one fixed normal boosted input isometry and one terminal smooth Y / original Mr instrument; no live acceleration or autonomous detector claim.",
        primary_text_read=dict(url="https://arxiv.org/html/1208.5840",section="8, equations 228-234; canonical generators only, not physical-position permission"))


def sine_interval(x):
    # Alternating series: 12 terms end negative, 13 terms end positive.
    terms=[(-1)**j*x**(2*j+1)/math.factorial(2*j+1) for j in range(13)]
    assert all(abs(terms[j+1])<abs(terms[j]) for j in range(12))
    return sum(terms[:12],F(0)),sum(terms,F(0))


def budgets():
    bt=json.loads((ROOT/"archive_1120_/1134/results.json").read_text(encoding="utf-8"))["exact_budgets"]
    old=bt["fractions"]
    row=bt["preparation_rows"]["code"]["fractions"]
    M,A,sigma,t=F(10**80),F(10**34),F(1,10**30),F(5,4)
    m=F(old["internal_mass_lower"])
    A1=F(row["A1"])
    p1=268+94*A
    n_input=464+2*sigma*p1
    N=F(2*10**6)
    y=2*sigma+t*A1/m
    ry=2*y+(N+F(1,2))/m
    eb=6*y+(6*N+3)/m
    ell=F(1,1000)
    leak=(ry/ell)**2
    d0=F(old["old_clock_effect_error"])
    clock=d0+2*eb
    joint=clock+leak
    # sqrt(2)<2 is used only for the complete instrument estimate.
    instrument=d0+eb+2*ry/ell
    assert n_input<N and m>M/2
    assert ry<F(5,10**30) and 2*eb<F(3,10**29)
    assert leak<F(25,10**54)
    assert joint<F(3126437,10**14)
    assert instrument<F(3126437,10**14)
    slo,shi=sine_interval(t)
    wlo,whi=sine_interval(F(25,16))
    plo,phi=((1+slo)/2)**3,((1+shi)/2)**3
    qlo,qhi=((1+wlo)/2)**3,((1+whi)/2)**3
    gap=qlo-phi-joint
    assert gap>F(7453608,10**8)
    vals=dict(mass=M,amplitude=A,sigma=sigma,proper_time=t,lab_time=F(25,16),
        endpoint_z=F(15,16),gamma=F(5,4),velocity=F(3,5),coarse_radius=ell,
        internal_mass_lower=m,total_first_momentum_norm=A1,
        LLP_internal_spin_input_bound=n_input,internal_spin_upper=N,
        rest_position_norm_upper=y,boost_position_norm_upper=ry,
        clock_blind_boost_norm_error=eb,bare_Mr_boost_error=2*eb,
        coarse_position_out_probability=leak,rest_full_comparison_error=d0,
        bare_clock_probability_error=clock,joint_in_clock_probability_error=joint,
        complete_terminal_instrument_error=instrument,
        correct_probability_lower=plo,correct_probability_upper=phi,
        wrong_lab_clock_probability_lower=qlo,wrong_lab_clock_probability_upper=qhi,
        finite_wrong_dictionary_gap_lower=gap)
    return dict(scope="Exact Fractions on the original LLP-correlated normal preparation; no reapplication of rest BT budget to boosted momentum.",
        fractions={k:str(v) for k,v in vals.items()},decimal_views={k:float(v) for k,v in vals.items()},
        time_window_legal=True,old_W_reset_used=False,
        source_pulse_replacement_error_added=False,
        position_probability_uses_squared_tail=True,instrument_uses_sqrt_tail=True,
        terminal_instrument_energy_in_autonomous_H=False)


def mass_lemma():
    # Three noncommuting positive-mass pairs. Sylvester is independently solved
    # with a Kronecker linear system and checked by finite interval quadrature.
    x=np.array([[0.,1.,0.],[1.,0.,.3],[0.,.3,0.]])
    y=np.array([[.2,.1j,.2],[-.1j,-.3,.4],[.2,.4,.1]],complex)
    rows=[]
    for angle in [.2,.7,1.1]:
        rot=exp_h(x,angle)
        aa=rot@np.diag([2.,4.,7.])@rot.conj().T
        bb=aa+.35*y
        mm=min(np.linalg.eigvalsh(aa).min(),np.linalg.eigvalsh(bb).min())
        cc=aa-bb
        delta=opnorm(cc)
        comm=opnorm(aa@bb-bb@aa)
        assert comm>.01
        for kval in [0.,.1,3.,30.]:
            ha=sqrtpos(aa@aa+kval*kval*np.eye(3))
            hb=sqrtpos(bb@bb+kval*kval*np.eye(3))
            ga,gb=ha-aa,hb-bb
            dd=ha-hb
            qq=ga@cc+cc@gb
            # vec(HA Z + Z HB) = (I kron HA + HB.T kron I) vec(Z).
            matrix=np.kron(np.eye(3),ha)+np.kron(hb.T,np.eye(3))
            zz=np.linalg.solve(matrix,-qq.reshape(-1,order="F")).reshape((3,3),order="F")
            residual=opnorm(dd-cc-zz)
            h0=np.sqrt(mm*mm+kval*kval)
            q=h0-mm
            assert opnorm(ga)<=q+2e-13 and opnorm(gb)<=q+2e-13
            assert opnorm(dd)<=2*delta+1e-12 and residual<1e-12
            deriv=-kval*np.linalg.solve(ha,dd)@np.linalg.inv(hb)
            step=2e-4
            def diff(k):
                return sqrtpos(aa@aa+k*k*np.eye(3))-sqrtpos(bb@bb+k*k*np.eye(3))
            numerical=(-diff(kval+2*step)+8*diff(kval+step)-8*diff(kval-step)+diff(kval-2*step))/(12*step)
            derivative_residual=opnorm(numerical-deriv)
            assert derivative_residual<2e-9 and opnorm(deriv)<=delta/mm+1e-12
            nodes,weights=np.polynomial.legendre.leggauss(160)
            end=24/h0
            ts=(nodes+1)*end/2
            weights=weights*end/2
            ea,va=np.linalg.eigh(ha); eb,vb=np.linalg.eigh(hb)
            integ=np.zeros((3,3),complex)
            for t,w in zip(ts,weights):
                left=(va*np.exp(-t*ea))@va.conj().T
                right=(vb*np.exp(-t*eb))@vb.conj().T
                integ-=w*left@qq@right
            quadrature_residual=opnorm(integ-zz)
            assert quadrature_residual<1e-12
            rows.append(dict(angle=angle,k=kval,mass_gap=mm,mass_difference=delta,
                noncommutativity=comm,sqrt_difference=opnorm(dd),
                sylvester_residual=residual,integral_residual=quadrature_residual,
                k_derivative_norm=opnorm(deriv),finite_difference_residual=derivative_residual))
    return dict(scope="Finite matrix diagnostics of the noncommuting lemma; the actual unbounded theorem uses bounded g and a weak Sylvester equation.",
        rows=rows,unbounded_A_semigroup_norm_integrated=False)


def old_clock():
    hc=1.5*np.eye(4)
    for j in range(3):
        hc[j,j+1]=hc[j+1,j]=.5*np.sqrt((j+1)*(3-j))
    e,v=np.linalg.eigh(hc)
    assert np.max(abs(e-np.arange(4)))<1e-12
    clock=v.conj().T@exp_h(hc,np.pi/2)@np.eye(4)[:,0]
    record=v.conj().T@np.diag([0.,0.,0.,1.])@v
    return e+15+2*np.sqrt(19),clock,record


def terminal_instrument_diagnostic(record):
    # Sample the smooth readout at inside / middle-shell / outside positions.
    # All four clock input matrix units are kept; this is an instrument algebra
    # check, not a discretization offered as a continuum position detector.
    angles=np.array([0.,np.pi/4,np.pi/2])
    lin,lout=np.diag(np.cos(angles)),np.diag(np.sin(angles))
    effects=[np.eye(4)-record,record]
    kraus=[np.kron(loc,e) for loc in [lin,lout] for e in effects]
    v=np.vstack(kraus)
    v0=np.vstack([np.kron(np.eye(3),e) for e in effects]+[np.zeros((12,12))]*2)
    completeness=opnorm(v.conj().T@v-np.eye(12))
    assert completeness<1e-12
    p=1e-8
    loc=np.sqrt(np.array([1-p-p*p,p*p,p]))
    code=np.kron(loc[:,None],np.eye(4))
    leak=float(np.vdot(loc,lout@lout@loc).real)
    postnorm=opnorm((v-v0)@code)
    fullnorm=opnorm(code.conj().T@v.conj().T@v@code-np.eye(4))
    assert fullnorm<1e-12
    assert postnorm<=np.sqrt(2*leak)+1e-12
    assert postnorm>1000*leak
    bell=code.reshape(12,4)/2
    assert abs(np.linalg.norm(v@bell)**2-1)<1e-12
    return dict(scope="Four Kraus completeness and complete clock input / reference diagnostic at three samples of the smooth effect; not detector construction.",
        completeness_residual=completeness,full_clock_code_Gram_residual=fullnorm,
        out_probability=leak,complete_instrument_difference=postnorm,
        sqrt_tail_bound=float(np.sqrt(2*leak)),probability_tail_cannot_bound_poststate=True)


def continuous_boost():
    # Normal 1D wave packet, full four clock masses, a two-dimensional field
    # register carrying unknown Q and its passive reference. The K dependence
    # contains packet(K-Pfield), so it is never replaced by an independent packet.
    gamma,sh,t=1.25,.75,1.25
    tp,xp=gamma*t,sh*t
    sigma=.005
    pfield=np.array([-.7,1.3])
    hs,clock,record=old_clock()
    field=np.array([[1.2,.31+.17j],[.31-.17j,2.1]])
    rows=[]
    for mass in [1000.,4000.]:
        am=[field+np.diag(np.sqrt((mass+h)**2+pfield**2)) for h in hs]
        eig=[np.linalg.eigh(a) for a in am]
        assert opnorm(am[0]@am[-1]-am[-1]@am[0])>1e-10
        def base(k):
            k=np.asarray(k)
            out=np.zeros((len(k),2,4,2),complex)
            normal=(2*sigma*sigma/np.pi)**.25
            for f in range(2):
                out[:,f,:,f]=normal*np.exp(-sigma*sigma*(k-pfield[f])**2)[:,None]*clock
            return out
        def evolve(arr,k,tt):
            out=np.empty_like(arr)
            for n,(mu,vv) in enumerate(eig):
                energy=np.sqrt(k[:,None]**2+mu[None,:]**2)
                vvstate=np.einsum("jf,njq->nfq",vv.conj(),arr[:,:,n,:])
                vvstate*=np.exp(-1j*tt*energy)[:,:,None]
                out[:,:,n,:]=np.einsum("fj,njq->nfq",vv,vvstate)
            return out
        def rest(k):
            return evolve(base(k),k,t)
        def apply_record(arr):
            return np.einsum("ab,nfbq->nfaq",record,arr)
        def boost(fn,k,reference=False):
            out=np.zeros((len(k),2,4,2),complex)
            for n in range(4):
                mu,vv=eig[0 if reference else n]
                for j in range(2):
                    energy=np.sqrt(k*k+mu[j]*mu[j])
                    kp=gamma*k-sh*energy
                    ep=gamma*energy-sh*k
                    arr=fn(kp)[:,:,n,:]
                    component=np.einsum("f,nfq->nq",vv[:,j].conj(),arr)
                    out[:,:,n,:]+=vv[:,j][None,:,None]*np.sqrt(ep/energy)[:,None,None]*component[:,None,:]
            return out
        def actual(k):
            return evolve(boost(base,k),k,tp)
        def phase_removed(k):
            return np.exp(1j*xp*k)[:,None,None,None]*actual(k)
        def gram(arr,weights):
            return np.einsum("nfaq,nfar,n->qr",arr.conj(),arr,weights)
        center=sh*(mass+hs.mean()+np.trace(field).real/2)
        span=2100.
        previous=None
        for order in [320,480]:
            xx,ww=np.polynomial.legendre.leggauss(order)
            kk=center+span*xx; ww=span*ww
            lab=actual(kk)
            transported=np.exp(-1j*xp*kk)[:,None,None,None]*boost(rest,kk)
            end_residual=np.sqrt(max(0.,np.linalg.eigvalsh(gram(lab-transported,ww)).max()))
            assert end_residual<2e-11
            normal_res=opnorm(gram(lab,ww)-np.eye(2))
            assert normal_res<1e-10
            rest_boost=boost(rest,kk)
            ref=boost(rest,kk,True)
            full_diff=np.sqrt(max(0.,np.linalg.eigvalsh(gram(rest_boost-ref,ww)).max()))
            comm=apply_record(rest_boost)-boost(lambda k:apply_record(rest(k)),kk)
            commnorm=np.sqrt(max(0.,np.linalg.eigvalsh(gram(comm,ww)).max()))
            probs=gram(apply_record(lab),ww)
            # Independent position derivative on the phase-removed end state.
            step=.002
            der=(-phase_removed(kk+2*step)+8*phase_removed(kk+step)-8*phase_removed(kk-step)+phase_removed(kk-2*step))/(12*step)
            ypos=np.sqrt(max(0.,np.linalg.eigvalsh(gram(der,ww)).max()))
            a1=np.sqrt(1/(4*sigma*sigma)+max(pfield**2))
            mmin=min(v[0].min() for v in eig)
            yr=sigma+t*a1/mmin
            ybound=2*yr+.5/mmin
            eb=6*yr+3/mmin
            assert ypos<=ybound and full_diff<=eb and commnorm<=2*eb
            # Bell Q/R input probability equals half the full input effect trace.
            bell_p=float(np.trace(probs).real/2)
            value=dict(mass=mass,quadrature_order=order,normalization_residual=normal_res,
                full_unknown_input_endpoint_residual=end_residual,
                complete_clock_blind_boost_difference=full_diff,bare_clock_commutator_norm=commnorm,
                analytic_clock_blind_bound=eb,actual_position_rms=ypos,analytic_position_bound=ybound,
                original_record_on_Bell_QR=bell_p,
                full_input_record_effect_eigenvalues=np.linalg.eigvalsh(probs).tolist())
            if previous is not None:
                convergence=max(abs(bell_p-previous["original_record_on_Bell_QR"]),abs(ypos-previous["actual_position_rms"]))
                assert convergence<2e-8
                value["quadrature_convergence"]=convergence
            previous=value
            rows.append(value)
    # Clock formula itself uses the original 4D Hc compression, bound by Gamma.
    hcvals=hs-(15+2*np.sqrt(19))
    proper=clock*np.exp(-1j*t*hcvals)
    wrong=clock*np.exp(-1j*tp*hcvals)
    pc=float(np.vdot(proper,record@proper).real)
    pw=float(np.vdot(wrong,record@wrong).real)
    assert abs(pc-((1+np.sin(t))/2)**3)<1e-12
    assert abs(pw-((1+np.sin(tp))/2)**3)<1e-12
    return dict(scope="Continuous 1D spin-zero external diagnostic, noncommuting finite internal masses, complete original four clock spectrum, all unknown-Q matrix units and passive Bell reference. Not the 3D Fock proof; spin is bounded analytically there.",
        rows=rows,old_clock_correct_probability=pc,old_clock_wrong_lab_probability=pw,
        finite_dictionary_gap=pw-pc,mass_momentum_correlations_kept=True,
        full_input_Gram_checked=True,original_288_device_reexecuted=False,
        terminal_instrument=terminal_instrument_diagnostic(record))


def compare(a,b,path=""):
    if isinstance(a,dict):
        assert a.keys()==b.keys(),path
        for k in a: compare(a[k],b[k],path+"/"+k)
    elif isinstance(a,list):
        assert len(a)==len(b),path
        for i,(x,y) in enumerate(zip(a,b)): compare(x,y,path+f"/{i}")
    elif isinstance(a,float):
        assert math.isclose(a,b,rel_tol=2e-7,abs_tol=2e-10),(path,a,b)
    else:
        assert a==b,(path,a,b)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    result=dict(round=1135,status="PASS",groups=4,new_axioms=0,science_experiments=0,stage="active",
        source_identity=sources(),exact_budgets=budgets(),noncommuting_mass_lemma=mass_lemma(),
        continuous_boost_and_old_clock=continuous_boost())
    target=HERE/"results.json"
    if args.write_results:
        target.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    else:
        compare(result,json.loads(target.read_text(encoding="utf-8")))
    print("1135 PASS: 4 mathematical groups; frozen sources; exact budgets; noncommuting mass lemma; continuous boost / original clock.")


if __name__=="__main__":
    main()
