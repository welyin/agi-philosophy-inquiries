"""761: loop-scaled finite-time matrix transport and first source response.

The scaling H_e=Hb,e+e*B is a NEW asymptotic input, identical to the original
operator only at e=1. Numerical propagation below uses a stated local
quadratic calibration with original neutral CAR coefficients, not the full graph.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_record_mass_feedback as old
import joint_reference_constraint_strata as background
import joint_background_source_closure as alg

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_loop_scale_transport_results.json'

def fock(h,d):
    dim=1<<len(h);B=np.zeros((dim,dim),complex)
    for j in range(dim):
        for i,c in old.car.quadratic({j:1.+0j},h,d).items():B[i,j]=c
    assert np.max(abs(B-B.conj().T))<1e-13
    return B

def original_data():
    src,psi,*_=background.completed(12,1.)
    phi=src['phi'][0,3,2]
    F=float(old.geom.F(phi));x=phi/np.sqrt(F)
    ids=[24,25,30,31];ix=np.ix_(ids,ids)
    h,d=old.mass_x(x);mh,md=old.mass_x(np.array([0.,0.,0.,0.,1.]))
    B=fock(h[ix],d[ix]);M=fock(mh[ix],md[ix])
    phase=M[12,0]/abs(M[12,0])
    vp=np.zeros(16,complex);vm=vp.copy()
    vp[0]=vm[0]=1/np.sqrt(2)
    vp[12]=phase/np.sqrt(2);vm[12]=-phase/np.sqrt(2)
    return phi,x,F,B,M,vp,vm,float(psi[0,3,2])

def scaled_operator_check():
    I=alg.I
    M=np.array([[Q(1),Q(2)],[Q(2),Q(-1)]],dtype=object)
    N=np.array([[Q(0),Q(1)],[Q(1),Q(0)]],dtype=object)
    B={(0,1):M,(0,0):N};D={(1,0):I};Bp={(0,0):M}
    rows=[]
    for e in (Q(1,7),Q(2,5),Q(1)):
        Hb={(2,0):-e*e*I/2,(2,2):-e*e*I/12,
            (1,1):-5*e*e*I/12,(0,4):I*Q(3,8),(0,2):I*Q(1,9)}
        H=alg.add(Hb,alg.scaled(B,e))
        rhs=alg.add(alg.comm(Hb,D),alg.scaled(Bp,-e))
        assert not alg.add(alg.comm(H,D),rhs,Q(-1))
        oldH=alg.add(Hb,B)
        oldrhs=alg.add(alg.comm(Hb,D),Bp,Q(-1))
        assert not alg.add(alg.comm(oldH,D),oldrhs,Q(-1))
        mismatch=alg.add(H,oldH,Q(-1))
        assert bool(mismatch)==(e!=1)
        # General exact telescoping commutator: extra generator is i[B,A],
        # while force on p=-i*e*D carries e.
        A=alg.sym(D,B)
        extra=alg.add(alg.comm(H,A),alg.comm(Hb,A),Q(-1))
        assert not alg.add(extra,alg.scaled(alg.comm(B,A),e),Q(-1))
        rows.append(dict(epsilon_exact=str(e),force_order_epsilon=True,
                         old_and_new_family_equal=(e==1),noncommuting_matrix_retained=True))
    return dict(rows=rows,exact_rational_curved_differential_identity=True,
                new_scaling_not_inferred_from_cognition=True)

def original_source_check(data):
    phi,x,F,B,M,vp,vm,psi=data
    n=32;vac={0:1/np.sqrt(2)}
    pair=(1<<30)|(1<<31)
    h,d=old.mass_x(np.array([0.,0.,0.,0.,1.]))
    amp=old.car.quadratic({0:1.+0j},h,d)[pair]
    vac[pair]=amp/abs(amp)/np.sqrt(2)
    means=[]
    for axis in range(5):
        e=np.eye(5)[axis];h,d=old.mass_x(e)
        means.append(float(old.dot(vac,old.car.quadratic(vac,h,d)).real))
    G=np.eye(5)+np.outer(x,x)/6
    u=1+float(x@x)/6
    ds=np.sqrt(2/u)*(np.eye(5)[4]-x[4]*x/(6*u))
    contracted=float(ds@G@np.array(means))
    expected=np.sqrt(F)*abs(old.matter.Y['s'])
    assert max(abs(np.array(means[:4])))<1e-14
    assert abs(means[4]-abs(old.matter.Y['s']))<1e-14
    assert abs(contracted-expected)<1e-14
    coeff=-np.cos(phi[4])*contracted/8
    assert coeff<0
    assert abs(np.vdot(vp,M@vp).real-means[4])<1e-14
    assert abs(np.vdot(vm,M@vm).real+means[4])<1e-14
    return dict(original_full_32_mode_mass_gradient_means=means,
                original_H5_force_contraction=contracted,
                analytic_contraction=expected,
                single_plus_record_t_squared_coefficient_for_unit_node_weight=coeff,
                opposite_input_gap_t_squared_coefficient_for_unit_node_weight=2*coeff,
                no_other_species_removed_from_analytic_graph_result=True,
                no_full_graph_evolution_simulated=True)

def exact_calibration(eps,nosc,data,t=.37,omega=1.):
    phi,x,F,B,M,vp,vm,psi=data
    annih=np.zeros((nosc,nosc),complex)
    for n in range(1,nosc):annih[n-1,n]=np.sqrt(n)
    metric=1+x[4]*x[4]/6
    X=np.sqrt(metric/(2*omega))*(annih+annih.conj().T)
    q=np.sqrt(eps)*X
    gen=np.kron(np.diag(omega*(np.arange(nosc)+.5)),np.eye(16))
    gen=gen+np.kron(np.eye(nosc),B)+np.sqrt(eps)*np.kron(X,M)
    eig,V=np.linalg.eigh(gen)
    phase=np.exp(-1j*t*eig)
    initial=np.zeros((16*nosc,2),complex);initial[:16,0]=vp;initial[:16,1]=vm
    state=V@(phase[:,None]*(V.conj().T@initial))
    be,bv=np.linalg.eigh(B)
    ferm=bv@(np.exp(-1j*t*be)[:,None]*(bv.conj().T@np.column_stack((vp,vm))))
    approx=np.zeros_like(state);approx[:16]=np.exp(-1j*t*omega/2)*ferm
    qev,qv=np.linalg.eigh(q)
    x5=x[4]+qev
    s=np.sqrt(2)*x5/np.sqrt(1+(float(x[:4]@x[:4])+x5*x5)/6)
    E=(qv*(.5+np.sin(s)/4))@qv.conj().T
    base=float(E[0,0].real)
    rows=[]
    for j in range(2):
        a=state[:,j].reshape(nosc,16)
        qmean=float(np.vdot(a,q@a).real)
        prob=float(np.vdot(a,E@a).real)
        source=float(np.vdot(a,a@M.T).real)
        predsource=float(np.vdot(ferm[:,j],M@ferm[:,j]).real)
        rows.append(dict(q_mean_over_epsilon=qmean/eps,record_delta_over_epsilon=(prob-base)/eps,
                         vector_error=float(np.linalg.norm(state[:,j]-approx[:,j])),
                         unscaled_source_error=abs(source-predsource)))
    return dict(epsilon=eps,oscillator_modes=nosc,rows=rows,baseline_probability=base)

def finite_time_check(data):
    phi,x,F,B,M,vp,vm,psi=data;t=.37;omega=1.
    be,bv=np.linalg.eigh(B)
    # Independent quadrature of the retarded classical oscillator response
    # driven by the full original neutral CAR matrix evolution.
    nodes,weights=np.polynomial.legendre.leggauss(100)
    times=(nodes+1)*t/2;weights=weights*t/2
    metric=1+x[4]*x[4]/6
    coeff=[]
    for v in (vp,vm):
        fv=bv@(np.exp(-1j*be[:,None]*times)*(bv.conj().T@v)[:,None])
        force=np.einsum('it,ij,jt->t',fv.conj(),M,fv).real
        coeff.append(-metric*float(np.sum(weights*np.sin(omega*(t-times))/omega*force)))
    u=1+float(x@x)/6
    ds=np.sqrt(2/u)*(1-x[4]*x[4]/(6*u))
    derivative=np.cos(phi[4])*ds/4
    recordcoeff=[derivative*c for c in coeff]
    values=[exact_calibration(e,18,data,t,omega) for e in (.16,.08,.04,.02)]
    errs=[];ratio=[]
    for val in values:
        err=max(abs(row['record_delta_over_epsilon']-recordcoeff[j]) for j,row in enumerate(val['rows']))
        errs.append(err)
        ratio.append(max(row['vector_error']/np.sqrt(val['epsilon']) for row in val['rows']))
    finer=exact_calibration(.04,24,data,t,omega)
    trunc=max(abs(finer['rows'][j][key]-values[2]['rows'][j][key])
              for j in range(2) for key in ('q_mean_over_epsilon','record_delta_over_epsilon','vector_error'))
    assert trunc<1e-10,trunc
    assert all(errs[j+1]<errs[j] for j in range(3)),errs
    assert errs[-1]<.001 and max(ratio)<.5,(errs,ratio)
    assert abs(recordcoeff[0]-recordcoeff[1])>.001
    return dict(time=t,retarded_q_coefficients=coeff,retarded_record_coefficients=recordcoeff,
                source_mass_commutator_norm=float(np.linalg.norm(B@M-M@B)),
                values=values,record_coefficient_errors=errs,vector_errors_over_sqrt_epsilon=ratio,
                oscillator_truncation_comparison=trunc,
                scope='Only a frozen local quadratic bosonic calibration with original neutral four-mode CAR masses. The full original H5 graph, all gauge links and its canonical trajectory are not simulated here; the original full-graph statement is analytic. This is not a replacement theory.')

def run():
    data=original_data()
    a=scaled_operator_check();b=original_source_check(data);c=finite_time_check(data)
    deps=('research_note_741.md','research_note_753.md','research_note_758.md',
          'research_note_759.md','research_note_760.md','joint_record_mass_feedback.py',
          'joint_reference_constraint_strata.py','joint_background_source_closure.py')
    return dict(round=761,tests_run=3,failures=0,errors=0,scaling_identity=a,
                actual_original_source=b,finite_time_local_calibration=c,
                dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
                scope='Conditional fixed-graph finite-time transport for the explicitly new family Hb,epsilon+epsilon B, with the original scalar background, full finite Gauss fiber and unchanged B coefficients. A first-order matter-induced response is compared with the same quantum bosonic baseline; original scalar-record probabilities have the same first correction. This does not prove the scaling follows from H1-H3, is accurate at epsilon=1, is the continuum limit, or yields the full Einstein equations.')
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(result,ensure_ascii=False,indent=2))

