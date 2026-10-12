"""1146 joint certificate. Default readonly; --write records new results.
Finite CAR diagnostics check algebra, never simulate the continuous field.
"""
from __future__ import annotations
import argparse, hashlib, importlib.util, json, math
from fractions import Fraction as Q
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent

def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    mod=importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    return mod

def fq(x):
    return {'fraction':str(x),'float':float(x)}

def provenance():
    sources=json.loads((HERE/'source_hashes.json').read_text(encoding='utf-8'))
    for name,sha in sources.items():
        assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==sha,name
    return len(sources)

def joint_budget(packet):
    lam=Q(11,7); delta=Q(1,1000); m=10**10
    eps=Q(1,10**16); A=10**100; tc=Q(1,4); r=Q(4,5)
    # Independently recompute the common two-spectrum bound using total variation.
    pkt=(delta+16*lam+Q(3,2)*(delta**2+2*delta*lam+16*lam**2)
         +(1+delta+lam)*20+500)/m
    kernel=2*lam*eps*(Q(3*m,4)+20)
    dc=(delta+Q(3,2)*tc*delta**2+(1+tc*delta)*20+tc*500)/m
    d=pkt+kernel
    C=33536/eps**3; q=8/eps**3
    source=lam*C*512*3/(eps*A)
    back=lam*512*C*q/(2*A*A)
    endpoint=lam*C*512*(3/(eps*A)+back)
    full=3*d+source; early=3*dc+tc*source
    x=delta*tc
    clock_gap=(1-r)*x/2*(1-x*x/2)-early
    sqrt5=Q(223607,100000); assert sqrt5*sqrt5>5
    source_gap=(3-sqrt5)/8-2*delta-2*full
    causal=4*lam**10/Q(math.factorial(10))/(1-lam/11)+2*source
    sea=(Q(400)+delta**2)/(m+Q(5*m,4))**2
    Zlow=(1-sea)**2
    invZ=1+Q(2,10**18); assert invZ**2*Zlow>1
    energy_col=(Q(5*m,4)+21)*(1+4*invZ)
    sqrt2=Q(99,70); assert sqrt2**2>2
    material_energy=sqrt2*energy_col
    source_energy=1073*A*A; vnorm=16*lam*C
    total_free=material_energy+source_energy+2*vnorm
    assert d<Q(246,10**8) and dc<Q(1451,10**11)
    assert full<Q(738,10**8) and early<Q(4353,10**11)
    assert source<Q(81,10**30) and back<Q(11,10**97)
    assert endpoint<Q(81,10**30) and vnorm<9*10**53
    assert sea<Q(8,10**19) and material_energy<9*10**10
    assert clock_gap>Q(24,10**6) and source_gap>Q(934,10000)
    assert causal<Q(12,100000)
    return {'pass':True,'kind':'joint_exact_rational_full_CAR_source_budget',
      'packet_error_upper':fq(pkt),'finite_kernel_error_upper':fq(kernel),
      'one_particle_error_upper':fq(d),'early_one_particle_upper':fq(dc),
      'quantum_source_error_upper':fq(source),'full_isometry_error_upper':fq(full),
      'early_full_isometry_upper':fq(early),'actual_clock_gap_lower':fq(clock_gap),
      'actual_scalar_vector_probability_gap_lower':fq(source_gap),
      'retarded_normalized_source_backreaction_upper':fq(back),
      'endpoint_interaction_mean_upper':fq(endpoint),'finite_menu_causal_upper':fq(causal),
      'negative_sea_occupancy_upper':fq(sea),'material_input_graph_upper':fq(material_energy),
      'source_initial_energy_upper':str(source_energy),'actual_total_free_energy_upper':fq(total_free),
      'unbounded_energy_independently_controlled':True,'full_apparatus_four_momentum_closed':False}

def car_group(old):
    n=4; e=np.eye(n,dtype=complex); ann=old.annihilators(n)
    ident=np.eye(2**n,dtype=complex); vac=np.zeros(2**n,complex); vac[0]=1
    create=lambda v:old.creation(v,ann)
    s0,s1=.17,.43
    w0=math.sqrt(s0)*e[:,0]+math.sqrt(1-s0)*e[:,2]
    w1=math.sqrt(s1)*e[:,1]+math.sqrt(1-s1)*e[:,3]
    sea=create(w0)@create(w1)@vac
    ns=[ann[i].conj().T@ann[i] for i in range(n)]
    empty=(ident-ns[0])@(ident-ns[1])@sea
    prob=float(np.vdot(empty,empty).real)
    assert abs(prob-(1-s0)*(1-s1))<1e-14
    omega=empty/math.sqrt(prob)
    J=np.column_stack([create(e[:,i])@omega for i in (0,1)])
    assert old.opnorm(J.conj().T@J-np.eye(2))<1e-14
    internal=old.evolution(np.array([[.02,.6],[.6,.13]],complex))
    target=np.eye(n,dtype=complex); target[:2,:2]=internal
    hh=np.array([[.11,.02j,.4,.09],[-.02j,-.2,.12j,.31],
                 [.4,-.12j,.1,.07j],[.09,.31,-.07j,-.01]],complex)
    U=old.evolution(hh,.013)@target
    V=U[:,:2]; W=target[:,:2]
    d0=float(np.linalg.norm(V[:,0]-W[:,0])); d1=float(np.linalg.norm(V[:,1]-W[:,1]))
    def correction(v,w):
        a=np.vdot(v,w); b2=max(0.,1-abs(a)**2)
        assert b2>1e-12  # this diagnostic uses the generic branch
        b=math.sqrt(b2); companion=(w-a*v)/b
        basis=np.column_stack([v,companion])
        rot=np.array([[a,-b],[b,a.conjugate()]],complex)
        out=np.eye(n)+basis@(rot-np.eye(2))@basis.conj().T
        assert old.opnorm(out.conj().T@out-np.eye(n))<2e-11
        assert np.linalg.norm(out@v-w)<2e-11
        assert abs(np.linalg.det(out)-1)<2e-11
        return out
    R1=correction(V[:,0],W[:,0]); v1=R1@V[:,1]
    R2=correction(v1,W[:,1]); R=R2@R1
    G1,G2,GR,GU=(old.lift(x) for x in (R1,R2,R,U))
    assert np.linalg.norm(R2@W[:,0]-W[:,0])<3e-11
    assert old.opnorm(R@V-W)<3e-11
    assert old.opnorm(G1-ident)<=d0+3e-11
    second=float(np.linalg.norm(v1-W[:,1]))
    assert second<=d0+d1+3e-11
    assert old.opnorm(G2-ident)<=second+3e-11
    assert old.opnorm(GR-ident)<=2*d0+d1+4e-11
    Xi=GR@GU@omega
    assert max(np.linalg.norm(ann[i]@Xi) for i in (0,1))<3e-11
    ideal=np.column_stack([create(W[:,i])@Xi for i in (0,1)])
    actual=GU@J
    correction_residual=old.opnorm(GR@actual-ideal)
    assert correction_residual<4e-11
    gap=old.opnorm(actual-ideal); bound=2*d0+d1
    assert gap<=bound+4e-11
    # Keep arbitrary reference explicitly; induced isometry norm tensor I_R is identical.
    ref=3; ent=np.array([1,1j,.3,-.2j,.7,.2],complex); ent/=np.linalg.norm(ent)
    actual_R=np.kron(actual,np.eye(ref))@ent
    ideal_R=np.kron(ideal,np.eye(ref))@ent
    assert np.linalg.norm(actual_R-ideal_R)<=bound+4e-11
    E=ns[1]@(ident-ns[0]); ER=np.kron(E,np.eye(ref))
    pa=float(np.vdot(actual_R,ER@actual_R).real)
    pi=float(np.vdot(ideal_R,ER@ideal_R).real)
    assert abs(pa-pi)<=bound+4e-11
    # Pure phase branch requires a companion; it must not become a scalar vacuum phase.
    theta=.17; RP=np.eye(n,dtype=complex)
    RP[0,0]=np.exp(1j*theta); RP[2,2]=np.exp(-1j*theta)
    phase_error=abs(np.exp(1j*theta)-1)
    assert abs(old.opnorm(old.lift(RP)-ident)-phase_error)<2e-14
    return {'pass':True,'kind':'unequal_sea_full_CAR_two_mixed_columns_and_reference',
      'finite_modes':n,'fock_dimension':2**n,'reference_dimension':ref,
      'sea_occupancies':[s0,s1],'empty_probability':prob,
      'column_errors':[d0,d1],'full_correction_bound':bound,
      'complete_isometry_difference':gap,'corrected_column_residual':correction_residual,
      'bare_readout_probability_difference':abs(pa-pi),'pure_phase_companion_checked':True,
      'continuous_field_simulation':False}

def run():
    frozen=provenance()
    pack=module('packet1146',HERE/'packet_check.py').run()
    saved=json.loads((HERE/'packet_results.json').read_text(encoding='utf-8'))
    assert pack==saved,'packet saved results mismatch'
    old=module('car1128',HERE.parent/'1128/check.py')
    joint=joint_budget(pack); car=car_group(old)
    return {'round':1146,'status':'PASS','check_groups':pack['check_groups']+2,
      'frozen_sources_checked':frozen,'packet_results_match':True,
      'packet_groups':pack['groups'],'joint_budget':joint,'full_CAR_diagnostic':car,
      'new_adopted_axioms':0,'scientific_calibration_increment':0,
      'goal_complete':False,'scope':'SAME_NORMAL_MOVING_MATERIAL_CLOCK_SCALAR_OPERATION_AND_QUANTUM_SOURCE'}

if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('--write',action='store_true')
    args=parser.parse_args(); data=run(); destination=HERE/'results.json'
    if destination.exists():
        assert json.loads(destination.read_text(encoding='utf-8'))==data,'saved result mismatch'
    elif args.write:
        destination.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':data['status'],'check_groups':data['check_groups'],
      'saved_results_match':destination.exists(),'joint_budget':data['joint_budget'],
      'full_CAR_diagnostic':data['full_CAR_diagnostic']},ensure_ascii=False,indent=2))
