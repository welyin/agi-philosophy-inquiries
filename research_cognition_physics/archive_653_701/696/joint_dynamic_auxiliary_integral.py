"""696: full original S9 integral for a specified dynamic spatial-loop background.
Exact moments and analytic Clifford reduction certify the result. Numerical
original matrices and independent cubature verify the object mapping, not RP.
"""
import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_dynamic_auxiliary_integral_results.json'
spec=importlib.util.spec_from_file_location('center696',HERE/'round696_drafts/dynamic_center_probe.py')
probe=importlib.util.module_from_spec(spec);spec.loader.exec_module(probe)
b=probe.b

def rising(x,n):
    out=F(1)
    for j in range(n):out*=x+j
    return out

def multiply(p,q):
    out={}
    for a,x in p.items():
        for c,y in q.items():
            key=tuple(i+j for i,j in zip(a,c));out[key]=out.get(key,F(0))+x*y
    return out

def exact_integral(m=1):
    # x=sqrt(rs)u, Y=sqrt(tau*(1-r)(1-s))v, z=rs.
    tau=F(1,1+m*m)
    d={(0,0,0):F(1,4),(1,0,0):F(1,2),(0,1,0):F(1,2),
       (2,0,0):tau/4,(1,1,0):F(1,2),(0,2,0):F(1,4),(0,0,1):(1-tau)/4}
    p={(0,0,0):F(1)}
    for _ in range(8):p=multiply(p,d)
    total=F(0);kept=0
    for (a,c,e),coef in p.items():
        if a%2 or c%2:continue
        k,l=a//2,c//2;radial_power=k+e
        radius=rising(F(3),radial_power)*rising(F(2),l)/rising(F(5),radial_power+l)
        angular=rising(F(1,2),k)/rising(F(3),k)*rising(F(1,2),l)/rising(F(2),l)
        total+=coef*tau**l*radius**2*angular;kept+=1
    mean=F(17,50)-F(13,200)*tau
    assert total>mean**8>0 and total<1
    # Exact Q(sqrt(1+m^2)) physical prefactor ((1+m/sqrt(1+m^2))/4)^16.
    radical=1+m*m;factor=(F(1),F(0));one=(F(1,4),F(m,4*radical))
    for _ in range(16):
        a,c=factor;x,y=one;factor=(a*x+radical*c*y,a*y+c*x)
    result=(factor[0]*total,factor[1]*total)
    if m==0:assert total==rising(F(9,2),16)/rising(F(9),16)
    return dict(average_spatial_Wilson_value=m,radicand=radical,
        auxiliary_mean=str(total),auxiliary_mean_decimal=float(total),
        mean_quartic=str(mean),strict_Jensen_lower=str(mean**8),
        physical_factor_algebraic=[str(x) for x in factor],
        full_scalar_mean_algebraic=[str(x) for x in result],
        full_scalar_mean_decimal=float(result[0])+np.sqrt(radical)*float(result[1]),
        polynomial_monomials=len(p),nonzero_integrated_monomials=kept)

def determinant_formula(e,pc,pw,m=1):
    r=float(e[0,pc]@e[0,pc]);s=float(e[1,pc]@e[1,pc])
    x=float(e[0,pc]@e[1,pc]);y=float(e[0,pw]@e[1,pw])
    d=((1+x)/2+y/(2*np.sqrt(1+m*m)))**2+(r*s-x*x)*m*m/(4*(1+m*m))
    return float(d)**8

def exact_clifford_check(pc,pw):
    # Original matrices have small Gaussian integer components. Products have
    # length16 and magnitude bounds far below2^53; conversion must be exact.
    def ip(a):
        r=np.rint(a.real).astype(np.int64);i=np.rint(a.imag).astype(np.int64)
        assert np.array_equal(a,r+1j*i);return r,i
    def mm(a,c):return a[0]@c[0]-a[1]@c[1],a[0]@c[1]+a[1]@c[0]
    def adj(a):return a[0].T,-a[1].T
    def eq(a,c):return all(np.array_equal(x,y) for x,y in zip(a,c))
    ident=(np.eye(16,dtype=np.int64),np.zeros((16,16),dtype=np.int64));zero=(ident[1],ident[1])
    ts=[ip(t) for t in b.internal.T];r=ip(b.prior.rep(np.eye(3),-np.eye(2),1))
    t0=ts[pc[0]];js=[mm(adj(t0),ts[a]) for a in (pc[1],pw[0],pw[1])]
    assert eq(mm(r,r),ident)
    for j in js:assert eq(mm(j,j),(-ident[0],ident[1]))
    assert eq(mm(r,js[0]),mm(js[0],r))
    for j in js[1:]:
        a,c=mm(r,j),mm(j,r);assert eq((a[0]+c[0],a[1]+c[1]),zero)
    for i in range(3):
        for k in range(i):
            a,c=mm(js[i],js[k]),mm(js[k],js[i]);assert eq((a[0]+c[0],a[1]+c[1]),zero)
    q=mm(r,js[0]);assert eq(mm(q,q),(-ident[0],ident[1]))
    assert np.trace(q[0])==np.trace(q[1])==0
    for j in (r,*js):assert eq(mm(q,j),mm(j,q))
    # All connected rotations within each6/4 subspace commute with the actual
    # centre link: these generate the accidental Spin6 x Spin4 stabilizer.
    pairs=0
    for block in (pc,pw):
        for ii,a in enumerate(block):
            for c in block[:ii]:
                generator=mm(adj(ts[a]),ts[c])
                assert eq(mm(r,generator),mm(generator,r));pairs+=1
    return dict(exact_Clifford_relations=True,central_Q_eigenmultiplicities=[8,8],
        stabilizer_generators_checked=pairs,color_axes=pc.tolist(),weak_axes=pw.tolist(),
        accidental_background_stabilizer_not_added_gauge_group=True)

def negative_kernel_certificate(integrals):
    endpoints={1:(F(1),F(1)),2:(F(1414213,10**6),F(1414214,10**6)),
               5:(F(2236067,10**6),F(2236068,10**6))}
    bounds=[]
    for row in integrals:
        d=row['radicand'];lo,hi=endpoints[d];assert lo*lo<=d<=hi*hi
        a,c=map(F,row['full_scalar_mean_algebraic']);assert a>0 and c>=0
        bounds.append((a+c*lo,a+c*hi))
    k00,k01,k11=bounds
    det_upper=k00[1]*k11[1]-k01[0]**2
    test_upper=k00[1]-F(2,5)*k01[0]+F(1,25)*k11[1]
    assert det_upper<0 and test_upper<F(-5,10**10)
    assert k01[0]**2>28**2*k00[1]*k11[1]
    return dict(exact_algebraic_entry_bounds=[[str(v) for v in row] for row in bounds],
        determinant_strict_upper=str(det_upper),rational_test=[1,'-1/5'],
        test_strict_upper=str(test_upper),test_upper_decimal=float(test_upper),
        normalized_offdiagonal_strictly_greater_than=28,
        complete_S9_averaged_fixed_seam_kernel_not_positive=True,
        NOT_full_Haar_Hb_physical_RP_certificate=True)

def run():
    cases=[];integrals=[];algebra=None
    for m,(first,second) in enumerate(((False,False),(True,False),(True,True))):
        exact=exact_integral(m);integrals.append(exact)
        coeff,pc,pw,factor,background=probe.background(first,second)
        if algebra is None:algebra=exact_clifford_check(pc,pw)
        reference=np.zeros((2,10));reference[:,pc[0]]=1
        frame_phase=np.linalg.det(np.einsum('xa,xaij->ij',reference,coeff))
        assert abs(abs(frame_phase)-1)<2e-12
        factor*=frame_phase
        expected_factor=float(F(exact['physical_factor_algebraic'][0]))+np.sqrt(1+m*m)*float(F(exact['physical_factor_algebraic'][1]))
        assert abs(factor/expected_factor-1)<3e-12
        rng=np.random.default_rng(69612);errors=[];fields=[]
        for _ in range(12):
            e=rng.normal(size=(2,10));e/=np.linalg.norm(e,axis=1)[:,None];fields.append(e)
        for a in (pc[0],pw[0]):
            for sign in (1,-1):
                e=np.zeros((2,10));e[0,a]=1;e[1,a]=sign;fields.append(e)
        for e in fields:
            compressed=np.linalg.det(np.einsum('xa,xaij->ij',e,coeff))/frame_phase;formula=determinant_formula(e,pc,pw,m)
            error=float(abs(compressed-formula)/max(abs(formula),1e-12));errors.append(error)
            assert error<2e-11
        rows=[]
        for nr,na in ((9,17),(10,18)):
            value,lo,hi=probe.integral(coeff,pc,pw,nr,na)
            value/=frame_phase
            error=float(abs(value/exact['auxiliary_mean_decimal']-1));assert error<2e-11
            rows.append(dict(radial_nodes=nr,angular_nodes=na,original_determinants=nr*nr*na*na,
                original_compressed_mean=[float(value.real),float(value.imag)],
                exact_rational_relative_error=error,diagnostic_node_real_range=[lo,hi]))
        cases.append(dict(first_weak_center=first,second_weak_center=second,
            original_matrix=background,formula_checks=len(fields),
            solver_frame_phase=[float(frame_phase.real),float(frame_phase.imag)],
            max_compressed_formula_relative_error=max(errors),independent_original_matrix_cubature=rows))
    certificate=negative_kernel_certificate(integrals)
    dependencies=('research_note_653.md','research_note_657.md','research_note_671.md',
        'research_note_675.md','research_note_695.md','joint_gauss_boundary_functional.py',
        'round696_drafts/dynamic_center_probe.py','round696_drafts/exterior_average_entry.py',
        'round696_drafts/exterior_average_entry_results.json')
    return dict(date='2026-10-02',round=696,tests_run=2,failures=0,errors=0,
        original_background=dict(spatial_sites=1,AP_time_sites=2,internal_channels=16,
            original_full_N_dimension=256,half_configurations=['identity_spatial_loop','weak_center_spatial_loop'],
            temporal_links=['identity','identity'],offdiagonal_plaquette_weak_center=-1,
            actual_spatial_hopping_between_distinct_sites=False),
        original_tensor_reduction=dict(**algebra,determinant_auxiliary='D_m^8',
            D='((1+x)/2+y/(2*sqrt(1+m^2)))^2+m^2*(r*s-x^2)/(4*(1+m^2))'),
        complete_S9_integrals=integrals,original_matrix_cases=cases,
        fixed_seam_negative_kernel=certificate,
        dependency_hashes={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in dependencies},
        scope=dict(complete_auxiliary_integral_on_specified_dynamic_background=True,
            no_group_or_Hb_integral_claim=True,no_general_RP_or_HF_identity=True,
            not_same_volume_as695_counterexample=True,not_general_space_propagation=True,
            no_new_spacetime_dimension_or_GR_theorem=True,all_original16_channels_retained=True),
        all_checks_passed=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args();r=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(round=696,all_checks_passed=True,
        auxiliary_means=[x['auxiliary_mean'] for x in r['complete_S9_integrals']],
        negative_test_upper=r['fixed_seam_negative_kernel']['test_upper_decimal'])))
