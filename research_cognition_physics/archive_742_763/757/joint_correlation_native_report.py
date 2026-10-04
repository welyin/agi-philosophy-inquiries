"""757: the same neutral correlation, source variance and original scalar report.

Analytic full-graph fourth-commutator theorem; finite exact differential algebra,
full96-mode CAR contractions and normal-packet integrals check its ingredients.
No full graph time simulation, new detector coupling or continuum limit.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from math import comb
from pathlib import Path
import numpy as np
import joint_record_mass_feedback as old

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_correlation_native_report_results.json'

def clean(p):return {k:v for k,v in p.items() if v}
def plus(a,b,scale=F(1)):
    out=dict(a)
    for k,v in b.items():out[k]=out.get(k,F(0))+scale*v
    return clean(out)
def multiply(a,b):
    out={}
    for i,x in a.items():
        for j,y in b.items():out[i+j]=out.get(i+j,F(0))+x*y
    return clean(out)
def derivative(a,n=1):
    for _ in range(n):a={i-1:i*v for i,v in a.items() if i}
    return a
def scale(a,s):return clean({k:s*v for k,v in a.items()})
def op_add(A,B,s=F(1)):
    out={k:dict(v) for k,v in A.items()}
    for k,p in B.items():out[k]=plus(out.get(k,{}),p,s)
    return {k:p for k,p in out.items() if p}
def compose(A,B):
    out={}
    for i,a in A.items():
        for j,b in B.items():
            for k in range(i+1):
                p=scale(multiply(a,derivative(b,k)),F(comb(i,k)))
                order=i-k+j;out[order]=plus(out.get(order,{}),p)
    return {k:p for k,p in out.items() if p}
def comm(A,B):return op_add(compose(A,B),compose(B,A),-1)

def exact_differential_check():
    # Original H5 Laplacian restricted to functions of x5 alone. The optional
    # drift=1 comparison checks that the quartic force term is drift-independent.
    metric={0:F(1),2:F(1,6)};metric_prime=derivative(metric)
    potential={0:F(2,7),1:F(-1,5),2:F(3,2),4:F(7,9)}
    cases=0
    for dim in (1,5):
        for w in (F(1),F(3,2)):
            hbar=F(7,10)
            for degree in range(1,7):
                f={degree:F(1)};answers=[]
                for coupling in (-1,0,1):
                    H={2:scale(metric,-hbar*hbar/(2*w)),
                       1:{1:-hbar*hbar*dim/(12*w)},
                       0:plus(potential,{1:F(coupling)})}
                    A={0:f};jets=[A]
                    for _ in range(4):A=comm(H,A);jets.append(A)
                    answers.append(jets)
                for j in range(4):
                    assert not op_add(op_add(answers[0][j],answers[2][j]),answers[1][j],-2)
                actual=op_add(op_add(answers[0][4],answers[2][4]),answers[1][4],-2)
                R=plus(scale(multiply(multiply(metric,metric),derivative(f,2)),3),
                       scale(multiply(multiply(metric,metric_prime),derivative(f)),2))
                expected={0:scale(R,2*hbar**4/w**2)}
                assert actual=={k:v for k,v in expected.items() if v}
                cases+=1
    return dict(exact_fraction_cases=cases,lower_orders_zero=True,
                fourth_commutator_multiplication_identity=True,
                no_boson_packet_or_time_discretization_in_identity=True)

def full_car_check():
    rng=np.random.default_rng(757);count=3;n=32*count
    keys=(0,1<<30,1<<31,(1<<30)|(1<<31));sign=(1,-1,-1,1)
    Ys=old.matter.Y['s'];max_error=0.;rows=[];gauge_error=0.
    for trial in range(4):
        x=rng.normal(size=(count,5))*.4
        h=np.zeros((n,n),complex);d=np.zeros_like(h);operators=[]
        for v in range(count):
            sl=slice(32*v,32*v+32);h[sl,sl],d[sl,sl]=old.mass_x(x[v])
        hop=np.zeros_like(h)
        for v,w in ((0,1),(1,2),(2,0)):
            R=old.matter.representation(old.matter.gauge.group_exp(rng.normal(size=8),3),
                                       old.matter.gauge.group_exp(rng.normal(size=3),2),np.exp(.31j))
            gauge_error=max(gauge_error,float(np.linalg.norm(R[:,30:32]-np.eye(32)[:,30:32])))
            sv=slice(v*32,(v+1)*32);sw=slice(w*32,(w+1)*32)
            spin=np.array([[.21,.13+.07j],[.13-.07j,-.21]])
            edge=(.31+.17j)*R@np.kron(np.eye(16),spin)
            hop[sv,sw]=edge;hop[sw,sv]=edge.conj().T
        for axis,(hi,di) in enumerate(old.LINEAR_MASS):
            a=np.zeros_like(h);b=np.zeros_like(d);a[:32,:32]=hi;b[:32,:32]=di
            operators.append((a,b,float(axis==4)))
        operators.extend(((h+hop,d,float(x[0,4])),(hop,np.zeros_like(d),0.)))
        images=[[old.car.quadratic({key:1.+0j},a,b) for key in keys] for a,b,_ in operators]
        error=0.
        for i,(a,b,ci) in enumerate(operators):
            means=sum(s*images[i][r].get(key,0.) for r,(s,key) in enumerate(zip(sign,keys)))
            error=max(error,float(abs(means)))
            for j,(_,_,cj) in enumerate(operators):
                value=sum(sign[r]*old.dot(images[i][r],images[j][r]) for r in range(4))
                error=max(error,float(abs(value-2*abs(Ys)**2*ci*cj)))
        max_error=max(max_error,error)
        rows.append(dict(full_CAR_modes=n,edges=3,quartic_contraction_error=error,
                         all_original_species_and_nonzero_hopping_retained=True))
    assert max_error<2e-12 and gauge_error<2e-14
    return dict(rows=rows,max_error=max_error,neutral_gauge_columns_error=gauge_error,
                two_bilinear_contraction='Tr(Pi A B) = 2 |Ys|^2 a5 b5 in the original selected-mode sector',
                full_H_not_time_evolved=True)

def packet(order):
    z,w=np.polynomial.legendre.leggauss(order)
    r=.7*(z+1)/2;wr=.7*w/2;y=.25*z;wy=.25*w
    r,y=np.meshgrid(r,y,indexing='ij');u=1+(r*r+y*y)/6
    weights=wr[:,None]*wy[None,:]*r**3/np.sqrt(u)
    weights*=np.exp(-2/(1-(r/.7)**2)-2/(1-(y/.25)**2))
    Z=float(weights.sum());weights/=Z
    F0=2/u;s=np.sqrt(F0)*y
    R=3*F0/4*np.cos(2*s)-5*s/48*np.sin(2*s)
    Ys=old.matter.Y['s'];kappa=.25
    moment=lambda a:float(np.sum(weights*a))
    fourth=2*kappa*abs(Ys)**2*moment(R)
    return dict(normalization=Z,mean_R=moment(R),mean_x5_squared=moment(y*y),
                fourth_probability_derivative_w1=fourth,
                leading_t4_probability_coefficient_w1=fourth/24,
                initial_total_energy_variance_difference=2*kappa*abs(Ys)**2*moment(y*y),
                pair_force_variance_difference=2*kappa*abs(Ys)**2)

def source_report_check():
    # Strict lower bound uses support alone: |x_H|<7/10, |x5|<1/4.
    Fmin=F(4800,2621);lower=F(9,16)*Fmin-F(5,192)
    assert lower>0
    small=packet(48);large=packet(80)
    err=max(abs(small[k]-large[k]) for k in large)
    assert err<2e-8 and large['mean_R']>float(lower)
    # Coordinate/metric connection at arbitrary points: grad_K x5=sqrt(F) d/ds.
    rng=np.random.default_rng(7571);error=0.
    for _ in range(20):
        x=rng.normal(size=5)*.4;u=1+x@x/6;phi=np.sqrt(2/u)*x
        ff=old.geom.F(phi)
        dx5=np.eye(5)[4]/np.sqrt(ff)+phi[4]*phi/(6*ff**1.5)
        error=max(error,float(np.linalg.norm(old.geom.inverse(phi)@dx5-np.sqrt(ff)*np.eye(5)[4])))
    assert error<2e-14
    return dict(normal_Gauss_packet=large,quadrature_48_to_80_max_change=err,
                exact_support_lower_bound_R=str(lower),strict_lower_bound_float=float(lower),
                curved_target_coordinate_identity_error=error,
                p=.5,q=.5,comparison_q=.25,w=1.,
                original_two_read_injection_bound='0 <= D_two <= hbar^2/(8 w)',
                strict_sign_from_support_not_quadrature=True,
                finite_time_window_not_certified=True)

def run():
    a=exact_differential_check();b=full_car_check();c=source_report_check()
    deps=('research_note_598.md','research_note_623.md','research_note_643.md','research_note_704.md',
          'research_note_706.md','research_note_717.md','research_note_743.md','research_note_746.md',
          'research_note_756.md','joint_record_mass_feedback.py','joint_fermion_gauss_completion.py')
    return dict(round=757,tests_run=3,failures=0,errors=0,
        exact_differential_identity=a,original_full_graph_CAR=b,common_source_and_report=c,
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope='On the original full finite graph, two physical neutral states with identical two-point functions and initial mean sources have different original scalar double-read probabilities. A strict local fourth derivative is fixed by the same pair-source variance and curved target, without ideal occupation measurement or new interaction. Proof keeps full H and fixed edges; calculations check algebra and normal preparations, not full time evolution. No autonomous terminal apparatus, continuous full-Gauss process map, finite-time resource window or quantum gravity is proved.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');a=p.parse_args();r=run()
    if a.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(round=757,tests_run=3,failures=0,coefficient=r['common_source_and_report']['normal_Gauss_packet']['leading_t4_probability_coefficient_w1'])))
