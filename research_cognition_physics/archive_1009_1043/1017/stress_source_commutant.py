"""1017: exact polynomial criterion for constant weighted scalar stress sources.

Fraction arithmetic finds the symmetric common commutant of EVERY coefficient
matrix of the polynomial Hessian. Floating point samples only cross-check the
on-shell divergence and an exact field-redefinition counterexample.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent
BASE=HERE.parents[1]
OUT=HERE/"stress_source_commutant_results.json"
SIGNS=np.array([-1.,1.,1.,1.])
ETA=np.diag(SIGNS)


def clean(p):
    return {power:F(value) for power,value in p.items() if value}


def add(*polys):
    result={}
    for p in polys:
        for power,value in p.items():
            result[power]=result.get(power,F(0))+value
    return clean(result)


def scale(p,a):
    return clean({power:value*F(a) for power,value in p.items()})


def multiply(p,q):
    result={}
    for power,a in p.items():
        for other,b in q.items():
            key=tuple(x+y for x,y in zip(power,other))
            result[key]=result.get(key,F(0))+a*b
    return clean(result)


def power(p,n):
    dim=len(next(iter(p))) if p else 2
    out={(0,)*dim:F(1)}
    for _ in range(n):
        out=multiply(out,p)
    return out


def derivative(p,axis):
    result={}
    for powers,value in p.items():
        if powers[axis]:
            target=list(powers);target[axis]-=1
            result[tuple(target)]=value*powers[axis]
    return clean(result)


def hessian(p,n=2):
    return [[derivative(derivative(p,i),j) for j in range(n)] for i in range(n)]


def evaluate(p,x):
    return sum(float(a)*math.prod(float(xx)**n for xx,n in zip(x,powers)) for powers,a in p.items())


def eval_exact(p,x):
    return sum((a*math.prod(xx**n for xx,n in zip(x,powers)) for powers,a in p.items()),F(0))


def serial_poly(p):
    return {",".join(map(str,key)):str(value) for key,value in sorted(p.items())}


def serial_matrix(q):
    return [[str(v) for v in row] for row in q]


def matmul(a,b):
    return [[sum((x*y for x,y in zip(row,col)),F(0)) for col in zip(*b)] for row in a]


def transpose(a):
    return [list(row) for row in zip(*a)]


def matrix_add(a,b,other_scale=F(1)):
    return [[x+other_scale*y for x,y in zip(row,other)] for row,other in zip(a,b)]


def symmetric_basis(n):
    result=[]
    for i in range(n):
        for j in range(i,n):
            q=[[F(0) for _ in range(n)] for _ in range(n)]
            q[i][j]=q[j][i]=F(1)
            result.append(q)
    return result


def commutator_polys(q,h):
    n=len(q)
    return [[add(*[add(scale(h[k][j],q[i][k]),scale(h[i][k],-q[k][j]))
                   for k in range(n)]) for j in range(n)] for i in range(n)]


def rref_nullspace(rows,ncols):
    a=[list(map(F,row)) for row in rows if any(row)]
    pivots=[];row=0
    for col in range(ncols):
        pivot=next((r for r in range(row,len(a)) if a[r][col]),None)
        if pivot is None:
            continue
        a[row],a[pivot]=a[pivot],a[row]
        value=a[row][col]
        a[row]=[x/value for x in a[row]]
        for r in range(len(a)):
            if r!=row and a[r][col]:
                value=a[r][col]
                a[r]=[x-value*y for x,y in zip(a[r],a[row])]
        pivots.append(col);row+=1
        if row==len(a):
            break
    basis=[]
    for free in range(ncols):
        if free not in pivots:
            v=[F(0)]*ncols;v[free]=F(1)
            for rr,pivot in enumerate(pivots):
                v[pivot]=-a[rr][free]
            basis.append(v)
    return len(pivots),basis,a[:len(pivots)]


def common_commutant(p,n=2):
    h=hessian(p,n)
    qb=symmetric_basis(n)
    comms=[commutator_polys(q,h) for q in qb]
    rows=[];labels=[]
    for i in range(n):
        for j in range(i+1,n):
            powers=sorted(set().union(*(set(c[i][j]) for c in comms)))
            for monomial in powers:
                rows.append([c[i][j].get(monomial,F(0)) for c in comms])
                labels.append(dict(matrix_entry=[i,j],field_power=list(monomial)))
    rank,vectors,reduced=rref_nullspace(rows,len(qb))
    matrices=[]
    for vec in vectors:
        q=[[sum((a*base[i][j] for a,base in zip(vec,qb)),F(0)) for j in range(n)] for i in range(n)]
        assert all(not entry for row in commutator_polys(q,h) for entry in row)
        matrices.append(q)
    return dict(rank=rank,dimension=len(vectors),matrices=matrices,rows=rows,row_labels=labels,rref=reduced)


def q_gradient(q,p):
    grad=[derivative(p,i) for i in range(len(q))]
    return [add(*[scale(grad[j],q[i][j]) for j in range(len(q))]) for i in range(len(q))]


def radial_primitive(field):
    # On R^n: integral_0^1 alpha(t phi).phi dt. It is a primitive only when
    # the coefficient one-form is closed; gradient equality is checked below.
    result={}
    for axis,component in enumerate(field):
        for monomial,value in component.items():
            target=list(monomial);target[axis]+=1
            result=add(result,{tuple(target):value/F(sum(monomial)+1)})
    return result


def primitive(q,p):
    target=q_gradient(q,p)
    w=radial_primitive(target)
    assert all(derivative(w,i)==target[i] for i in range(len(q)))
    return w


def rotate_poly(p,o):
    n=len(o)
    linear=[]
    for row in o:
        item={}
        for j,value in enumerate(row):
            monomial=[0]*n;monomial[j]=1
            item[tuple(monomial)]=value
        linear.append(clean(item))
    out={}
    for monomial,value in p.items():
        term={(0,)*n:value}
        for j,exponent in enumerate(monomial):
            term=multiply(term,power(linear[j],exponent))
        out=add(out,term)
    return out


def potential(m1,m2,a,b,c):
    return clean({(2,0):F(m1)/2,(0,2):F(m2)/2,(4,0):F(a),(0,4):F(b),(2,2):F(c)})


def symbolic_classification():
    # Build the coefficient relations from five independent parameter basis
    # polynomials, then group by field monomial and unknown A,B,D. No samples
    # of x,y, or numeric ranks at selected field points, enter this step.
    parameters={"m1_squared":{(2,0):F(1,2)},"m2_squared":{(0,2):F(1,2)},
                "a":{(4,0):F(1)},"b":{(0,4):F(1)},"c":{(2,2):F(1)}}
    equations={}
    for name,p in parameters.items():
        h=hessian(p)
        for label,q in zip(["A","B","D"],symmetric_basis(2)):
            for monomial,value in commutator_polys(q,h)[0][1].items():
                key=",".join(map(str,monomial))
                equations.setdefault(key,{})[f"{label}*{name}"]=str(value)
    expected={"0,0":{"B*m1_squared":"-1","B*m2_squared":"1"},
              "2,0":{"B*a":"-12","B*c":"2"},
              "0,2":{"B*b":"12","B*c":"-2"},
              "1,1":{"A*c":"4","D*c":"-4"}}
    assert equations==expected
    return equations


def on_shell_divergence(q,potential_poly,w,phi,p,pp):
    # Directly differentiate the unreduced T^{mu nu} tensor. pp[mu,nu,i]
    # is symmetric in spacetime derivatives and chosen to obey box(phi)=dV.
    q=np.asarray(q,float)
    grad_v=np.array([evaluate(derivative(potential_poly,i),phi) for i in range(len(phi))])
    grad_w=np.array([evaluate(derivative(w,i),phi) for i in range(len(phi))])
    lap=np.einsum("m,mmi->i",SIGNS,pp)
    assert np.max(np.abs(lap-grad_v))<2e-11
    direct=np.zeros(4)
    for nu in range(4):
        for mu in range(4):
            direct[nu]+=SIGNS[mu]*SIGNS[nu]*(pp[mu,mu]@q@p[nu]+p[mu]@q@pp[mu,nu])
            if mu==nu:
                direct[nu]-=SIGNS[nu]*(sum(SIGNS[rho]*(pp[mu,rho]@q@p[rho]) for rho in range(4))
                                      +grad_w@p[mu])
    predicted=(q@grad_v-grad_w)@p.T*SIGNS
    return direct,predicted


def compare(fresh,saved):
    if isinstance(fresh,dict):
        assert fresh.keys()==saved.keys()
        for key in fresh:
            compare(fresh[key],saved[key])
    elif isinstance(fresh,list):
        assert len(fresh)==len(saved)
        for a,b in zip(fresh,saved):compare(a,b)
    elif isinstance(fresh,float):
        assert math.isclose(fresh,saved,rel_tol=3e-8,abs_tol=3e-10),(fresh,saved)
    else:
        assert fresh==saved,(fresh,saved)


def run():
    equations=symbolic_classification()
    families=[
        ("separated_quartic",(1,2,1,2,0),2),
        ("genuine_mixed_quartic",(1,2,1,2,1),1),
        ("pure_cross_quartic",(1,2,0,0,1),1),
        ("rotated_decoupled_quartic",(1,1,F(1,2),F(1,2),3),2),
        ("rotated_quartic_mass_split",(1,2,F(1,2),F(1,2),3),1),
        ("equal_mass_free",(2,2,0,0,0),3),
        ("unequal_mass_free",(1,2,0,0,0),2),
        ("equal_mass_separated_quartic",(1,1,1,1,0),2),
        ("O2_invariant_radial_quartic",(1,1,1,1,2),1),
    ]
    o=[[F(3,5),F(-4,5)],[F(4,5),F(3,5)]]
    identity=[[F(1),F(0)],[F(0),F(1)]]
    assert matmul(transpose(o),o)==identity
    rng=np.random.default_rng(1017)
    rows=[];max_divergence=max_formula=0.;jet_count=0
    for name,params,expected_dim in families:
        v=potential(*params)
        space=common_commutant(v)
        assert space["dimension"]==expected_dim,(name,space["dimension"])
        vrot=rotate_poly(v,o)
        rotated=common_commutant(vrot)
        assert rotated["dimension"]==expected_dim
        primitives=[]
        for q in space["matrices"]:
            w=primitive(q,v)
            qrot=matmul(transpose(o),matmul(q,o))
            assert all(not entry for row in commutator_polys(qrot,hessian(vrot)) for entry in row)
            assert primitive(qrot,vrot)==rotate_poly(w,o)
            primitives.append(serial_poly(w))
            for _ in range(2):
                phi=.3*rng.normal(size=2);p=.3*rng.normal(size=(4,2))
                pp=.3*rng.normal(size=(4,4,2));pp=(pp+np.swapaxes(pp,0,1))/2
                grad_v=np.array([evaluate(derivative(v,i),phi) for i in range(2)])
                pp[0,0]=sum(pp[i,i] for i in range(1,4))-grad_v
                direct,predicted=on_shell_divergence(q,v,w,phi,p,pp)
                max_divergence=max(max_divergence,float(np.max(np.abs(direct))))
                max_formula=max(max_formula,float(np.max(np.abs(direct-predicted))))
                assert np.max(np.abs(direct))<2e-11
                assert np.max(np.abs(direct-predicted))<2e-11
                jet_count+=1
        rows.append(dict(name=name,parameters={key:str(F(value)) for key,value in zip(["m1_squared","m2_squared","a","b","c"],params)},
            potential=serial_poly(v),exact_constraint_matrix=[[str(x) for x in row] for row in space["rows"]],
            coefficient_labels=space["row_labels"],constraint_rank=space["rank"],
            symmetric_commutant_dimension=space["dimension"],
            symmetric_commutant_basis=[serial_matrix(q) for q in space["matrices"]],
            corresponding_W=primitives,rotated_commutant_dimension=rotated["dimension"],
            all_polynomial_coefficients_checked_exactly=True))

    # A visibly mixed polynomial which is precisely two decoupled modes.
    pplus=[[F(1,2),F(1,2)],[F(1,2),F(1,2)]]
    pminus=[[F(1,2),F(-1,2)],[F(-1,2),F(1,2)]]
    zero=[[F(0),F(0)],[F(0),F(0)]]
    assert matmul(pplus,pplus)==pplus and matmul(pminus,pminus)==pminus
    assert matmul(pplus,pminus)==zero and matrix_add(pplus,pminus)==identity
    rank_sources,_,_=rref_nullspace([[pplus[0][0],pplus[0][1],pplus[1][1]],
                                  [pminus[0][0],pminus[0][1],pminus[1][1]]],3)
    assert rank_sources==2
    vhidden=potential(1,1,F(1,2),F(1,2),3)
    wp=primitive(pplus,vhidden);wm=primitive(pminus,vhidden)
    plus={(1,0):F(1),(0,1):F(1)}
    minus={(1,0):F(1),(0,1):F(-1)}
    expected_wp=add(scale(power(plus,2),F(1,4)),scale(power(plus,4),F(1,4)))
    expected_wm=add(scale(power(minus,2),F(1,4)),scale(power(minus,4),F(1,4)))
    assert wp==expected_wp and wm==expected_wm and add(wp,wm)==vhidden

    # Independent nonidentical background metrics remain in their original
    # sectors. Only scalar field coordinates change; no Einstein-Hilbert
    # identity or gravitational generation is claimed by this calculation.
    density_residuals=[]
    for _ in range(8):
        phi=.4*rng.normal(size=2);p=.4*rng.normal(size=(4,2))
        u=(phi[0]+phi[1])/math.sqrt(2);v=(phi[0]-phi[1])/math.sqrt(2)
        du=(p[:,0]+p[:,1])/math.sqrt(2);dv=(p[:,0]-p[:,1])/math.sqrt(2)
        lengths_p=np.exp(.2*rng.normal(size=4));lengths_m=np.exp(.2*rng.normal(size=4))
        sqrt_p=float(np.prod(lengths_p));sqrt_m=float(np.prod(lengths_m))
        inverse_p=SIGNS/(lengths_p**2);inverse_m=SIGNS/(lengths_m**2)
        original=-sqrt_p*(.5*np.sum(inverse_p*du**2)+.5*u*u+u**4)
        original-=sqrt_m*(.5*np.sum(inverse_m*dv**2)+.5*v*v+v**4)
        qp=np.asarray(pplus,float);qm=np.asarray(pminus,float)
        rotated=-sqrt_p*(.5*sum(inverse_p[mu]*(p[mu]@qp@p[mu]) for mu in range(4))+evaluate(wp,phi))
        rotated-=sqrt_m*(.5*sum(inverse_m[mu]*(p[mu]@qm@p[mu]) for mu in range(4))+evaluate(wm,phi))
        error=abs(original-rotated);assert error<2e-12
        density_residuals.append(error)

    # Genuine interaction: diag(1,0) has no allowed scalar W. Show a closed
    # loop / path dependence exactly, then a nonconserved candidate tensor.
    vcross={(2,2):F(1)}
    badq=[[F(1),F(0)],[F(0),F(0)]]
    field=q_gradient(badq,vcross)
    curl=add(derivative(field[1],0),scale(derivative(field[0],1),-1))
    assert curl=={(1,1):F(-4)}
    path_x_then_y=F(0);path_y_then_x=F(4)  # endpoint (1,2)
    # Integral along the final x-segment at y=2 from 0 to 1.
    independent_path=sum((value*F(2)**j/F(i+1) for (i,j),value in field[0].items()),F(0))
    assert independent_path==path_y_then_x
    badw=radial_primitive(field)
    defect=[add(field[i],scale(derivative(badw,i),-1)) for i in range(2)]
    defect_at=[eval_exact(poly,[F(1),F(2)]) for poly in defect]
    assert defect_at==[F(4),F(-2)]
    phi=np.array([1.,2.]);p=np.zeros((4,2));p[1]=[1.,0.]
    pp=np.zeros((4,4,2));pp[0,0]=-np.array([evaluate(derivative(vcross,i),phi) for i in range(2)])
    direct,predicted=on_shell_divergence(badq,vcross,badw,phi,p,pp)
    assert np.max(np.abs(direct-predicted))<1e-13 and direct[1]==4

    # A symmetric commutant is a vector space, not generally a commutative
    # multiplication algebra: the equal-mass free example makes this visible.
    b0,b1,_=symmetric_basis(2)
    free_commutator=matrix_add(matmul(b0,b1),matmul(b1,b0),F(-1))
    assert free_commutator!=zero
    sources=[BASE/"archive_301_341/research_note_326.md",
             BASE/"archive_342_369/research_note_358.md",
             BASE/"archive_1009_/research_note_1016.md",
             BASE/"archive_1009_/1016/NEXT.md",
             BASE/"archive_1009_/1009/input_dependency_ledger_v0_1.md"]
    return dict(round=1017,new_calibration_groups=1,cumulative_test_groups=3795,
        new_cognitive_axioms=0,all_scientific_calibrations_passed=True,
        scope="canonical real scalar fields on flat background; target R^n; constant real symmetric Q; specified weighted stress tensor and scalar W",
        full_multigraviton_no_go_claimed=False,gravity_generated=False,
        target_dimension_or_matter_content_selected=False,
        symbolic_family_commutator_equations=equations,
        exact_polynomial_families=rows,
        rational_rotation=serial_matrix(o),
        rotation_preserves_all_exact_commutant_dimensions=True,
        on_shell_jet_calibration=dict(samples=jet_count,maximum_divergence=max_divergence,
                                      maximum_direct_formula_difference=max_formula),
        mixed_vertex_counterexample=dict(potential=serial_poly(vhidden),
            Q_plus=serial_matrix(pplus),Q_minus=serial_matrix(pminus),
            W_plus=serial_poly(wp),W_minus=serial_poly(wm),
            projector_identities_exact=True,source_span_rank=rank_sources,
            two_sectors_decoupled_in_u_v=True,
            visible_x2_y2_coefficient="3",
            different_metric_scalar_density_checks=len(density_residuals),
            maximum_scalar_density_redefinition_residual=max(density_residuals),
            both_gravity_actions_unchanged_by_scalar_coordinate_rotation=True,
            gravity_action_numerically_rederived=False,
            ordinary_vertex_graph_is_basis_dependent=True),
        impossible_weighted_source=dict(potential=serial_poly(vcross),Q=serial_matrix(badq),
            exact_curl=serial_poly(curl),endpoint=["1","2"],
            path_x_then_y=str(path_x_then_y),path_y_then_x=str(path_y_then_x),
            candidate_radial_W=serial_poly(badw),gradient_defect_at_endpoint=[str(x) for x in defect_at],
            direct_on_shell_divergence=direct.tolist(),predicted_on_shell_divergence=predicted.tolist()),
        equal_mass_free_boundary=dict(symmetric_commutant_dimension=3,
            commutator_of_two_symmetric_solutions=serial_matrix(free_commutator),
            unique_sector_decomposition_claimed=False),
        retained_freedom=["additive constant in each W","coefficient of scalar total stress",
            "number of independent spin2 fields not selected by linear source analysis",
            "nonlinear gravitational consistency beyond the stated counterexample",
            "other stress improvements, derivative interactions, quantum effects and matter representations"],
        analytic_obligations="on-shell divergence and Poincare lemma give gradient/commutant equivalence; spectral projectors give global orthogonal separation on R^n; polynomial coefficient systems are exact certificates for stated families",
        historical_source_sha256={str(p.relative_to(BASE)).replace("\\","/"):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--write",action="store_true")
    args=parser.parse_args();result=run()
    if args.write:
        with OUT.open("x",encoding="utf8") as stream:
            json.dump(result,stream,ensure_ascii=False,indent=2);stream.write("\n")
    else:compare(result,json.loads(OUT.read_text(encoding="utf8")))
    print(json.dumps(dict(round=1017,passed=result["all_scientific_calibrations_passed"],
        dimensions={row["name"]:row["symmetric_commutant_dimension"] for row in result["exact_polynomial_families"]},
        on_shell=result["on_shell_jet_calibration"],
        independent_rotated_sources=result["mixed_vertex_counterexample"]["source_span_rank"],
        density_redefinition_residual=result["mixed_vertex_counterexample"]["maximum_scalar_density_redefinition_residual"],
        path_difference=result["impossible_weighted_source"]["path_y_then_x"]),indent=2))
