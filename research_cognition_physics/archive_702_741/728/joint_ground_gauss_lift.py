"""728: isolated original-matter bands and their Gauss-compatible local lift.

Analytic graph-size-independent mass bound for the declared flat-link branch.
This does not replace574's non-flat Einstein source or prove its fermionic gap.
"""
import argparse
from fractions import Fraction as Q
from math import isqrt
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_matter_ground_source as old
HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_ground_gauss_lift_results.json'


def sqrt_bounds(x):
    scale=10**10
    k=isqrt(x.numerator*scale*scale//x.denominator)
    return Q(k,scale),Q(k+1,scale)


def bounds_check():
    M=Q(2);hlo=Q('.95')*Q('.65');hhi=Q('1.05')*Q('.67')
    slo=Q('.94')*Q('.54');shi=Q('1.06')*Q('.55')
    h0lo=Q('1.05')*Q('.65');h0hi=hhi;s0lo=Q('.54');s0hi=Q('.55')
    def F(h,s):return M-(h*h+s*s)/6
    def coef(lo,hi,hlo,slo,hhi,shi):
        return lo/sqrt_bounds(F(hlo,slo))[1],hi/sqrt_bounds(F(hhi,shi))[0]
    a=coef(hlo,hhi,hlo,slo,hhi,shi);b=coef(slo,shi,hlo,slo,hhi,shi)
    a0=coef(h0lo,h0hi,h0lo,s0lo,h0hi,s0hi)
    b0=coef(s0lo,s0hi,h0lo,s0lo,h0hi,s0hi)
    da=max(a[1]-a0[0],a0[1]-a[0]);db=max(b[1]-b0[0],b0[1]-b[0])
    Ys={k:sqrt_bounds(x) for k,x in dict(u=Q('1.53'),d=Q('.26'),e=Q('.0449'),
                                        nu=Q('.1744'),s=Q('.1042')).items()}
    lower_e=Ys['e'][0]*a0[0];Dmin=Ys['nu'][0]*a0[0];MRmax=Ys['s'][1]*b0[1]
    # sqrt(M^2+4D^2)-M is increasing in D and decreasing in M.
    neutral_lower=(sqrt_bounds(MRmax**2+4*Dmin**2)[0]-MRmax)/2
    uniform_lower=min(lower_e,neutral_lower)
    perturbation=Ys['u'][1]*da+Ys['s'][1]*db
    certified=uniform_lower-perturbation
    assert Ys['u'][0]>Ys['d'][1]>Ys['e'][1]
    assert certified>Q('.013')
    src=old.old.wall.old.fields(32);phi=src['q']['phi']
    h=np.linalg.norm(phi[...,:4],axis=-1);s=phi[...,4];f=old.matter.original.F(phi)
    f0=old.matter.original.F(np.array([0,old.CENTER[0],0,0,old.CENTER[1]]))
    da_actual=float(np.max(abs(h/np.sqrt(f)-old.CENTER[0]/np.sqrt(f0))))
    db_actual=float(np.max(abs(s/np.sqrt(f)-old.CENTER[1]/np.sqrt(f0))))
    assert da_actual<float(da) and db_actual<float(db)
    values=dict(uniform_mass_gap_lower=uniform_lower,mass_perturbation_upper=perturbation,
                graph_independent_gap_lower=certified,neutral_mass_lower=neutral_lower,
                radial_coefficient_deviation_upper=da,singlet_coefficient_deviation_upper=db)
    return dict(rational_bounds={k:dict(exact=str(v),decimal=float(v)) for k,v in values.items()},
                sampled32_profile_deviations=[da_actual,db_actual],
                profile_sites=32**3,proof_not_based_on_sampled_minimum=True,
                bound_only_for_flat_links_with_original_three_direction_kinetic=True)


def charges():
    q=np.zeros(32)
    q[0:12]=np.tile([Q(2,3),Q(2,3),Q(-1,3),Q(-1,3)],3)
    q[12:18]=2/3;q[18:24]=-1/3;q[24:28]=[0,0,-1,-1];q[28:30]=-1
    quark=q.copy();quark[24:]=0
    lepton=q-quark
    return q,quark,lepton


def moments(P,G,normal_trace=0.):
    I=np.eye(len(P))
    mean=float(np.trace(P@G).real/2+normal_trace/2)
    variance=float(np.trace(P@G@(I-P)@G).real/2)
    return mean,variance


def band_and_charge_check():
    B1,_=old.star();dim=len(B1);n=dim//2
    # Separate the original star mass and edge matrices, preserving Nambu order.
    _,G=old.star();K=-G;Bm=B1-K
    q,quark,lepton=charges()
    qs=[np.diag(np.tile(x,4)) for x in (q,quark,lepton)]
    Qs=[old.bdg(x,np.zeros_like(x)) for x in qs]
    rows=[];errs=[]
    for t in (0.,.1,.3,.6,1.):
        B=Bm+t*K;d=old.data(B,Qs[0]);P=d['P']
        full,noise=moments(P,Qs[0],float(np.trace(qs[0])))
        qm,qn=moments(P,Qs[1],float(np.trace(qs[1])))
        lm,ln=moments(P,Qs[2],float(np.trace(qs[2])))
        errs.extend([abs(full),abs(noise),abs(qm-8),abs(lm+8),
                     float(np.max(abs(B@Qs[0]-Qs[0]@B)))])
        assert d['gap']>.013
        rows.append(dict(hopping_scale=t,quasiparticle_gap=d['gap'],total_EM_charge=full,
                         EM_variance=noise,quark_EM_charge=qm,lepton_EM_charge=lm))
    assert max(errs)<1e-10
    # Direct mass/kinetic anticommutation in the constant-mass comparison.
    center=old.bdg(*old.radial(*old.CENTER)[0])
    hs,ds=old.radial(*old.CENTER)[0]
    comparison=old.bdg(old.block([hs]*4),old.block([ds]*4))
    anti=float(np.max(abs(comparison@K+K@comparison)))
    assert anti<1e-12
    return dict(rows=rows,max_charge_and_commutator_error=max(errs),
                homogeneous_mass_edge_anticommutator_error=anti,
                physical_modes=n,comparison_one_node_gap=float(min(abs(np.linalg.eigvalsh(center)))),
                gap_path_is_analytically_controlled_not_just_five_samples=True,
                quark_and_lepton_character_cancellation_uses_original_complete_generation=True)


def representation_and_record_check():
    B,G=old.star();d=old.data(B,G);P=d['P'];q,_,_=charges()
    rng=np.random.default_rng(728)
    C=old.matter.gauge.group_exp(rng.normal(size=8),3)
    beta=.173
    W=np.diag([np.exp(3j*beta),np.exp(-3j*beta)])
    R=old.matter.representation(C,W,np.exp(1j*beta))
    RR=old.block([R]*4);RN=old.block([RR,RR.conj()])
    # H-invariant projector and zero EM character; SU3 has no continuous characters.
    errors=[float(np.max(abs(RN@B@RN.conj().T-B))),
            float(np.max(abs(RN@P@RN.conj().T-P)))]
    zeta=old.matter.representation(np.eye(3)*np.exp(2j*np.pi/3),-np.eye(2),np.exp(1j*np.pi/3))
    errors.append(float(np.max(abs(zeta-np.eye(32)))))
    # Electroweak residual center is the same action as a color center.
    ew=old.matter.representation(np.eye(3),-np.eye(2),np.exp(1j*np.pi/3))
    color=old.matter.representation(np.eye(3)*np.exp(-2j*np.pi/3),np.eye(2),1.)
    errors.append(float(np.max(abs(ew-color))))
    # Sterile occupation record is gauge invariant for independent endpoints too.
    Rlocal=[]
    for _ in range(4):
        Rlocal.append(old.matter.representation(old.matter.gauge.group_exp(rng.normal(size=8),3),
            old.matter.gauge.group_exp(rng.normal(size=3),2),np.exp(1j*rng.normal())))
    L=old.block(Rlocal);LN=old.block([L,L.conj()])
    reflection=np.eye(256);reflection[30,30]=-1;reflection[158,158]=-1
    errors.append(float(np.max(abs(reflection@LN-LN@reflection))))
    neutral=old.neutral(*old.CENTER)
    parity=np.array([(-1)**i.bit_count() for i in range(16)])
    parity_error=float(np.linalg.norm(parity*neutral['chi']-neutral['chi']))
    errors.append(parity_error)
    assert max(errors)<2e-12
    Qnormal=np.diag(np.tile(q,4));Qn=old.bdg(Qnormal,np.zeros_like(Qnormal))
    Pr=reflection@P@reflection
    charge,noise=moments(Pr,Qn,float(np.trace(Qnormal)))
    assert max(abs(charge),abs(noise))<1e-10
    return dict(max_representation_error=max(errors),neutral_ground_even_parity_error=parity_error,
                reflected_ground_total_EM_charge=charge,reflected_ground_EM_variance=noise,
                actual_sterile_record_preserves_Gauss_and_parity=True,
                full_Fock_line_character_proof_is_analytic=True,
                old574_nonflat_stabilizer_reused_not_reproved=True)


def run():
    results=dict(rational_global_band_bound=bounds_check(),
                 complete_generation_character=band_and_charge_check(),
                 quotient_and_actual_record=representation_and_record_check())
    deps=('research_note_574.md','research_note_591.md','research_note_598.md',
          'research_note_604.md','research_note_634.md','research_note_727.md',
          'joint_gravity_material_coordinates.py','joint_matter_ground_source.py',
          'round728_drafts/macroscopic_source_entry.py')
    return dict(round=728,tests_run=3,failures=0,errors=0,results=results,
                dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
                scope='Original matter conditional ground-line Gauss lift:574 nonflat stabilizer removes character obstruction conditional on isolated ground. Original flat-link Weyl branch has a rational graph-independent mass gap bound and complete-generation EM cancellation, supplying a genuine local Gauss band and admissible record states. The flat branch is not the original nonflat Einstein initial state. No nonflat full-graph gap, invariant slow band, autonomous reset, continuum chirality or Einstein backreaction claimed.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');a=p.parse_args();r=run()
    if a.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(r['results'],ensure_ascii=False,indent=2))
