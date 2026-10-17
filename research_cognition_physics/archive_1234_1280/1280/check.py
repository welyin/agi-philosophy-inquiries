"""1280: fixed natural conformal reference; exact bounds plus non-certifying diagnostics.
Default mode is read-only. --write creates results.json exactly once.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse, hashlib, json, math
import numpy as np

HERE=Path(__file__).resolve().parent

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def bump(u):
    x=u+2.0
    if not 0.0<x<1.0:
        return 0.0,0.0
    h=4*x*(1-x)
    exponent=1-1/h
    if exponent < -740:
        return 0.0,0.0
    f=math.exp(exponent)
    fp=f*(1-2*x)/(4*x*x*(1-x)*(1-x))
    return f,fp

def derivative(u,y):
    r,rp,eta=y[:3]
    f,fp=bump(u)
    q=fp*fp/36000
    z=-u
    xi=-eta*r*r
    psip=rp*f+r*fp
    return np.array([rp,-q*r,1/(r*r),fp*fp,z*r*r*fp*fp,
                     z*psip*psip,xi*psip*psip,
                     (r*rp+z*q*r*r)*f*f,
                     (2*u+2)*r*r*(2*q)-(u+2)*2*rp*rp],dtype=float)

def integrate(n):
    y=np.zeros(9);y[0]=1.0;y[2]=-1.0
    h=-1.0/n
    min_ratio=1.0;max_rate=1.0
    for i in range(n):
        u=-1+i*h
        k1=derivative(u,y)
        k2=derivative(u+h/2,y+h*k1/2)
        k3=derivative(u+h/2,y+h*k2/2)
        k4=derivative(u+h,y+h*k3)
        y+=h*(k1+2*k2+2*k3+k4)/6
        un=u+h
        ratio=-y[2]*y[0]*y[0]/(-un)
        min_ratio=min(min_ratio,float(ratio))
        max_rate=max(max_rate,float(1/ratio))
    r,rp,eta=y[:3]
    E,B,J,M,C,Eg=-y[3:]
    area=1-r*r
    return dict(steps=n,E=float(E),B=float(B),J=float(J),M_eta=float(M),C=float(C),
        r_at_minus_2=float(r),rprime_at_minus_2=float(rp),eta_at_minus_2=float(eta),
        past_caustic_coordinate=float(-2-r/rp),area_density_change=float(area),
        geometry_correction=float(Eg),rescaling_identity_residual=float(J-B-C),
        geometry_identity_residual=float(area-B/18000-Eg),
        relative_modular_source_difference=float((M-B)/B),
        minimum_weight_ratio_diagnostic=min_ratio,maximum_clock_rate_diagnostic=max_rate)

def null_tensor_checks():
    # Finite exact sanity checks of the analytic null-tensor argument, not a proof by sampling.
    rows=[]
    for r,q in ((F(1),F(1,36000)),(F(499,500),F(1,1000)),(F(3,4),F(2,7))):
        g=[[F(0) for _ in range(4)] for _ in range(4)]
        g[0][1]=g[1][0]=F(-1);g[2][2]=g[3][3]=r*r
        partner=(1,0,2,3);weight=(F(-1),F(-1),1/(r*r),1/(r*r))
        ric=[[F(0) for _ in range(4)] for _ in range(4)];ric[0][0]=2*q
        inv=lambda a,b:weight[a] if b==partner[a] else F(0)
        scalar=sum((inv(a,b)*ric[a][b] for a in range(4) for b in range(4)),F(0))
        ric_square=[[sum((ric[a][c]*inv(c,d)*ric[d][b]
                      for c in range(4) for d in range(4)),F(0))
                      for b in range(4)] for a in range(4)]
        def riem(a,b,c,d):
            return (g[a][c]*ric[b][d]+g[b][d]*ric[a][c]
                    -g[a][d]*ric[b][c]-g[b][c]*ric[a][d])/2
        contraction=sum((riem(a,b,c,d)*weight[a]*weight[b]*weight[c]*weight[d]
                    *riem(partner[a],partner[b],partner[c],partner[d])
                    for a in range(4) for b in range(4)
                    for c in range(4) for d in range(4)),F(0))
        reconstructed=[[sum((inv(a,c)*riem(a,b,c,d) for a in range(4) for c in range(4)),F(0))
                        for d in range(4)] for b in range(4)]
        assert scalar==0 and contraction==0 and reconstructed==ric
        assert all(v==0 for row in ric_square for v in row)
        rows.append(dict(r=str(r),q=str(q),scalar_curvature='0',ricci_square_tensor='0',
                         riemann_square='0',ricci_reconstructed=True))
    return rows

def compute():
    A=F(1,1000);rmin=1-2*A
    dc=3*A/(4*rmin*rmin)
    de=1-rmin*rmin
    total=dc+de*(1+dc)
    target=F(1,200)
    assert dc==F(375,498002) and de==F(999,250000)
    assert total==F(591628623,124500500000)<target
    margin=target-total
    assert margin==F(30873877,124500500000)
    clock_excess=de/(1-de)
    assert clock_excess==F(999,249001)
    tensor_rows=null_tensor_checks()
    gauss=[]
    for n in (128,256,512):
        z,w=np.polynomial.legendre.leggauss(n)
        val=sum(float(wi)*bump(float(zi/2-1.5))[1]**2/2 for zi,wi in zip(z,w))
        gauss.append(dict(nodes=n,E=val))
    rows=[integrate(n) for n in (4096,8192)]
    for row in rows:
        assert 4<row['E']<36
        assert row['r_at_minus_2']>float(rmin)
        assert row['past_caustic_coordinate'] < -500
        assert abs(row['relative_modular_source_difference']) < float(total)
        assert abs(row['rescaling_identity_residual'])<1e-9
        assert abs(row['geometry_identity_residual'])<1e-10
        assert row['minimum_weight_ratio_diagnostic']>float(rmin*rmin)
        assert row['maximum_clock_rate_diagnostic']<float(1+clock_excess)
    files=['admission.md','fixed_spec.md','proof.md','check.py','NEXT.md',
           '../research_note_1280.md','../1279/proof.md','../1279/results.json',
           '../_shared/ct_hg_adoption_after1278.md',
           '../_shared/jacobson_reconstruction_after1278.md']
    return dict(schema='round1280_natural_conformal_reference_v1',round=1280,status='PASS_WITH_SCOPE',
        author_check_groups=1,new_adopted_axioms=0,
        strict_certificate=dict(method='Analytic continuum proof plus exact rational bounds; diagnostic samples do not prove the QFT identities',
            A_q_upper=str(A),r_pulse_lower=str(rmin),delta_C=str(dc),delta_eta=str(de),
            relative_total_bound=str(total),target=str(target),strict_margin=str(margin),
            clock_rate_excess_bound=str(clock_excess),null_tensor_checks=tensor_rows),
        numerical_diagnostics=dict(scope='Fixed RK4/Gauss only; grid differences are not certified numerical errors',
            gauss=gauss,rk4=rows,RK_B_difference=abs(rows[0]['B']-rows[1]['B']),
            RK_M_difference=abs(rows[0]['M_eta']-rows[1]['M_eta'])),
        physical_scope=dict(reference='Natural eta conformal Minkowski vacuum on noncompact transverse covering',
            bulk_Hadamard=True,absolute_minimal_reference_stress_zero_in_fixed_compatible_prescription=True,
            reference_is_1279_artanh_vacuum=False,torus_source_zero_certified=False,
            working_state_locally_normal=True,global_Fock_finite_energy_certified=False,
            response_unit='per transverse coordinate area; finite lightfront subalgebra only as scoped in proof',
            mean_source_semiclassical_identity=True,fluctuation_or_stability_certificate=False,
            exact_proper_time_KMS_certified=False,actual_clock_instrument_certified=False,
            complete_global_finite_resources_certified=False,full_goal_completed=False),
        primary_sources=['https://arxiv.org/html/0806.0803',
            'https://doi.org/10.1103/PhysRevD.16.1712',
            'https://arxiv.org/pdf/hep-th/0008211','https://arxiv.org/html/1502.04106'],
        bindings={name:sha(HERE/name) for name in files})

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args()
    result=compute();p=HERE/'results.json'
    if args.write:
        with p.open('x',encoding='utf-8',newline='\n') as f:
            json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
        mode='CREATED_ONCE'
    else:
        saved=json.loads(p.read_text('utf-8'))
        assert saved==result,'Saved results differ; default mode never writes'
        mode='READ_ONLY_MATCH'
    fine=result['numerical_diagnostics']['rk4'][-1]
    print(json.dumps(dict(status=result['status'],mode=mode,
        strict_relative_bound=result['strict_certificate']['relative_total_bound'],
        diagnostic_relative_difference=fine['relative_modular_source_difference'],
        diagnostic_E=fine['E'],full_goal_completed=False),ensure_ascii=False))
