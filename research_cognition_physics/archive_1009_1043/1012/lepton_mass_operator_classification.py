"""1012: direct two-Weyl lepton mass operators at dimension <= 5.

Unknown hypercharges are retained throughout.  This is a finite algebraic
calibration at a constant nonzero Higgs background, not a complete EFT basis,
an experimental fit, or a proof of a high-energy completion.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations, combinations_with_replacement
import argparse
import hashlib
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
OUT = HERE / "lepton_mass_operator_classification_results.json"
LABELS = [f"{kind}{i+1}" for kind in "NDE" for i in range(3)]
SYMMETRIC_INDICES = list(combinations_with_replacement(range(9), 2))


def opnorm(a):
    return float(np.linalg.norm(a, ord=2))


def fraction_rank(rows):
    if not rows:
        return 0
    a = [[F(x) for x in row] for row in rows]
    rank = 0
    for col in range(len(a[0])):
        pivot = next((j for j in range(rank, len(a)) if a[j][col]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        value = a[rank][col]
        a[rank] = [x/value for x in a[rank]]
        for j in range(len(a)):
            if j != rank:
                value = a[j][col]
                a[j] = [x-value*y for x,y in zip(a[j], a[rank])]
        rank += 1
        if rank == len(a):
            break
    return rank


def operators():
    out = []
    specifications = [
        ("LL_bare", 3, "0", list(combinations(range(3), 2))),
        ("EE_bare", 3, "0", list(combinations_with_replacement(range(3), 2))),
        ("AE", 4, "H", [(i,j) for i in range(3) for j in range(3)]),
        ("BE", 4, "Hbar", [(i,j) for i in range(3) for j in range(3)]),
        ("AA", 5, "HH", list(combinations_with_replacement(range(3), 2))),
        ("BB", 5, "HbarHbar", list(combinations_with_replacement(range(3), 2))),
        ("AB", 5, "HHbar", [(i,j) for i in range(3) for j in range(3)]),
        ("EE_Hnorm", 5, "HHbar", list(combinations_with_replacement(range(3), 2))),
    ]
    for kind, dimension, scalar_type, indices in specifications:
        for i,j in indices:
            out.append(dict(kind=kind,dimension=dimension,scalar_type=scalar_type,
                            i=i,j=j,key=f"{kind}_{i+1}{j+1}"))
    assert len(out) == 54
    return out


def charge(item, h, t):
    x = [t,-t,F(0)]
    i,j = item["i"],item["j"]
    kind = item["kind"]
    if kind in ("LL_bare", "AB"):
        return x[i]+x[j]-2*h
    if kind in ("EE_bare", "EE_Hnorm"):
        return 4*h-x[i]-x[j]
    if kind == "AE":
        return x[i]+2*h-x[j]
    if kind == "BE":
        return x[i]-x[j]
    if kind == "AA":
        return x[i]+x[j]
    if kind == "BB":
        return x[i]+x[j]-4*h
    raise ValueError(kind)


def elementary_spinors():
    # Each Weyl component is represented by a linear form in 18 independent
    # Grassmann generators.  Arrays retain spin and weak indices separately.
    f = np.zeros((9,2,18),dtype=complex)
    for i in range(9):
        for s in range(2):
            f[i,s,2*i+s] = 1
    ls = np.stack([f[:3],f[3:6]],axis=1)
    es = f[6:9]
    return f,ls,es


def wedge(u,v):
    # Upper triangular entries are polynomial coefficients of theta_a theta_b.
    return np.outer(u,v)-np.outer(v,u)


def spinor_bilinear(u,v):
    return wedge(u[0],v[1])-wedge(u[1],v[0])


def invariant_spinors(ls,es,higgs):
    a = higgs[1]*ls[:,0]-higgs[0]*ls[:,1]
    b = higgs[0].conjugate()*ls[:,0]+higgs[1].conjugate()*ls[:,1]
    return a,b,es


def polynomial(item,ls,es,higgs):
    a,b,e = invariant_spinors(ls,es,higgs)
    i,j,kind = item["i"],item["j"],item["kind"]
    if kind == "LL_bare":
        return spinor_bilinear(ls[i,0],ls[j,1])-spinor_bilinear(ls[i,1],ls[j,0])
    if kind in ("EE_bare","EE_Hnorm"):
        result = spinor_bilinear(e[i],e[j])
        if kind == "EE_Hnorm":
            result = result*float(np.vdot(higgs,higgs).real)
    elif kind == "AA":
        result = spinor_bilinear(a[i],a[j])
    elif kind == "BB":
        result = spinor_bilinear(b[i],b[j])
    elif kind == "AB":
        return spinor_bilinear(a[i],b[j])
    elif kind == "AE":
        return spinor_bilinear(a[i],e[j])
    elif kind == "BE":
        return spinor_bilinear(b[i],e[j])
    else:
        raise ValueError(kind)
    # Convention: a diagonal symmetric-matrix coefficient multiplies 1/2 psi psi.
    return result*(.5 if i==j else 1.)


def matrix_from_polynomial(p):
    m = np.zeros((9,9),dtype=complex)
    reconstructed = np.zeros((18,18),dtype=complex)
    for a,b in SYMMETRIC_INDICES:
        value = p[2*a,2*b+1]
        m[a,b] = m[b,a] = value
        reconstructed[2*a,2*b+1] += value
        reconstructed[2*b+1,2*a] -= value
        if a != b:
            reconstructed[2*a+1,2*b] -= value
            reconstructed[2*b,2*a+1] += value
    assert np.linalg.norm(p-reconstructed) < 1e-11
    return m


def mass_basis(item,v=1.):
    _,ls,es = elementary_spinors()
    return matrix_from_polynomial(polynomial(item,ls,es,np.array([0.,v],dtype=complex)))


def exact_mass_vector(item):
    mat = mass_basis(item)
    result = []
    for i,j in SYMMETRIC_INDICES:
        assert abs(mat[i,j].imag)<1e-14 and abs(mat[i,j].real-round(mat[i,j].real))<1e-14
        result.append(F(int(round(mat[i,j].real))))
    return result


def symmetric_representation(generator):
    n = len(generator)
    basis = list(combinations_with_replacement(range(n),2))
    locations = {pair:i for i,pair in enumerate(basis)}
    action = np.zeros((len(basis),len(basis)),dtype=complex)
    # Differentiate commuting products. Lorentz-contracted Weyl bilinears are
    # symmetric in their combined species labels, owing to Grassmann statistics.
    for k,(a,b) in enumerate(basis):
        for c in range(n):
            action[locations[tuple(sorted((c,b)))],k] += generator[a,c]
            action[locations[tuple(sorted((a,c)))],k] += generator[b,c]
    return action


def su2_invariant_dimensions():
    sigmas = [np.array([[0,1],[1,0]],dtype=complex),
              np.array([[0,-1j],[1j,0]],dtype=complex),
              np.diag([1.,-1.])]
    generators = [s/2 for s in sigmas]
    types = {"0":1,"H":2,"Hbar":2,"HH":3,"HbarHbar":3,"HHbar":4}
    results = {}
    expected = {"0":9,"H":9,"Hbar":9,"HH":6,"HbarHbar":6,"HHbar":15}
    for scalar_type,s in types.items():
        rows=[]
        for t in generators:
            g = np.zeros((9,9),dtype=complex)
            for i in range(3):
                g[np.ix_([i,3+i],[i,3+i])] = t
            gf = symmetric_representation(g)
            if scalar_type == "0":
                gs = np.zeros((1,1),dtype=complex)
            elif scalar_type == "H":
                gs = t.T
            elif scalar_type == "Hbar":
                gs = -t.conj().T
            elif scalar_type == "HH":
                gs = symmetric_representation(t)
            elif scalar_type == "HbarHbar":
                gs = symmetric_representation(-t.conj())
            else:
                gs = np.kron(t.T,np.eye(2))+np.kron(np.eye(2),-t.conj().T)
            rows.append(np.kron(gf,np.eye(s))+np.kron(np.eye(45),gs))
        sv = np.linalg.svd(np.vstack(rows),compute_uv=False)
        dim = 45*s-int(np.count_nonzero(sv>1e-10))
        assert dim == expected[scalar_type],(scalar_type,dim)
        results[scalar_type] = dict(symmetric_fermion_dimension=45,scalar_dimension=s,
            independent_su2_invariant_dimension=dim)
    return results


def anomaly_certificate(h,t):
    q = h/3
    x = [t,-t,F(0)]
    l = [xi-h for xi in x]
    e = [2*h-xi for xi in x]
    uc,dc = -q-h,-q+h
    # Every fermion below is left-handed, including the conjugate E=e_R^c.
    fields=[]
    for i in range(3):
        fields.extend([dict(name=f"Q{i+1}",multiplicity=6,Y=q),
                       dict(name=f"uc{i+1}",multiplicity=3,Y=uc),
                       dict(name=f"dc{i+1}",multiplicity=3,Y=dc),
                       dict(name=f"L{i+1}",multiplicity=2,Y=l[i]),
                       dict(name=f"E{i+1}",multiplicity=1,Y=e[i])])
    traces = {
        "SU3_cubed":F(3*(2-1-1)),
        "SU3_squared_U1":3*(2*q+uc+dc),
        "SU2_squared_U1":9*q+sum(l),
        "gravity_squared_U1":sum(row["multiplicity"]*row["Y"] for row in fields),
        "U1_cubed":sum(row["multiplicity"]*row["Y"]**3 for row in fields),
    }
    assert all(value==0 for value in traces.values())
    assert sum(row["multiplicity"] for row in fields)==45
    return dict(traces={key:str(value) for key,value in traces.items()},
                weyl_components=45,weak_doublets=12,ordinary_su2_parity_anomaly=False,
                complete_global_anomaly_claimed=False)


def serialize_matrix(a):
    if opnorm(a.imag) < 1e-13:
        return a.real.tolist()
    return [[[float(z.real),float(z.imag)] for z in row] for row in a]


def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for key in a:
            compare(a[key],b[key])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b):
            compare(x,y)
    elif isinstance(a,float):
        assert math.isclose(a,b,rel_tol=3e-9,abs_tol=3e-11),(a,b)
    else:
        assert a==b,(a,b)


def run():
    ops = operators()
    by_key = {item["key"]:item for item in ops}
    f,ls,es = elementary_spinors()
    rng = np.random.default_rng(1012)
    max_su2_error = max_mixed_identity = max_symmetric_bilinear = 0.
    eps = np.array([[0.,1.],[-1.,0.]])
    for _ in range(8):
        z = rng.normal(size=4)
        z /= np.linalg.norm(z)
        u = np.array([[z[0]+1j*z[3],z[2]+1j*z[1]],
                      [-z[2]+1j*z[1],z[0]-1j*z[3]]])
        assert opnorm(u.conj().T@u-np.eye(2)) < 1e-14
        assert abs(np.linalg.det(u)-1) < 1e-14
        hv = rng.normal(size=2)+1j*rng.normal(size=2)
        l1 = np.einsum("ab,ibsk->iask",u,ls)
        for item in ops:
            p0 = polynomial(item,ls,es,hv)
            p1 = polynomial(item,l1,es,u@hv)
            err = float(np.linalg.norm(p0-p1))
            assert err < 2e-12
            max_su2_error = max(max_su2_error,err)
        a,b,_ = invariant_spinors(ls,es,hv)
        for i,j in combinations(range(3),2):
            lhs = spinor_bilinear(a[i],b[j])-spinor_bilinear(a[j],b[i])
            rhs = polynomial(by_key[f"LL_bare_{i+1}{j+1}"],ls,es,hv)*np.vdot(hv,hv)
            err = float(np.linalg.norm(lhs-rhs))
            assert err < 2e-12
            max_mixed_identity = max(max_mixed_identity,err)
            symmetry = float(np.linalg.norm(spinor_bilinear(a[i],b[j])-spinor_bilinear(b[j],a[i])))
            assert symmetry < 1e-14
            max_symmetric_bilinear=max(max_symmetric_bilinear,symmetry)
        assert abs(hv@eps@hv)<1e-12

    invariant_dimensions=su2_invariant_dimensions()
    # Exact polynomial charge equations a*r+b=0 determine all exceptional r.
    exceptional=set()
    charges={}
    for item in ops:
        b=charge(item,F(1),F(0))
        a=charge(item,F(1),F(1))-b
        charges[item["key"]]=dict(r_coefficient=str(a),constant=str(b))
        if a:
            exceptional.add(-b/a)
    assert exceptional=={F(x) for x in [-4,-2,-1,0,1,2,4]}
    rows=[]
    for r in sorted(exceptional | {F(-3),F(1,3),F(3),F(5)}):
        h,t=F(1,2),r/2
        x=[t,-t,F(0)]
        qtilde=x+[xi-2*h for xi in x]+[2*h-xi for xi in x]
        allowed=[item for item in ops if charge(item,h,t)==0]
        vectors=[exact_mass_vector(item) for item in allowed]
        rank=fraction_rank(vectors)
        ward_pairs=[(i,j) for i,j in SYMMETRIC_INDICES if qtilde[i]+qtilde[j]==0]
        assert rank==len(ward_pairs)
        # Every vector is in the exact Ward kernel; equal dimension proves
        # equality for this calibration point, not a continuum-by-sampling proof.
        for vector in vectors:
            assert all((qtilde[i]+qtilde[j])*entry==0
                       for (i,j),entry in zip(SYMMETRIC_INDICES,vector))
        counts={kind:sum(item["kind"]==kind for item in allowed)
                for kind in ["LL_bare","EE_bare","AE","BE","AA","BB","AB","EE_Hnorm"]}
        expected = {F(0):15,F(1):7,F(2):11,F(4):7}.get(abs(r),5)
        assert rank==expected,(r,rank)
        zero_dimension=sum(value==0 for value in qtilde)
        assert zero_dimension==(3 if r in (F(0),F(2),F(-2)) else 1)
        # NN restriction can be checked without assuming closure of that block.
        nn_support = {(i,j) for item in allowed for i in range(3) for j in range(i,3)
                      if abs(mass_basis(item)[i,j])>1e-12}
        aa_support = {(item["i"],item["j"]) for item in allowed if item["kind"]=="AA"}
        assert nn_support==aa_support
        assert all(not np.any(np.abs(mass_basis(item)[:3,:3])>1e-12)
                   for item in allowed if item["kind"]!="AA")
        rows.append(dict(r=str(r),exceptional=r in exceptional,
            allowed_operator_count=len(allowed),operator_counts=counts,
            independent_mass_matrix_dimension=rank,ward_kernel_dimension=len(ward_pairs),
            zero_charge_dimension=zero_dimension,
            compatible_with_three_simple_positive_masses_necessary=zero_dimension>=3,
            upper_component_pairs=[f"N{i+1}N{j+1}" for i,j in sorted(nn_support)],
            allowed_operator_keys=[item["key"] for item in allowed],
            exact_anomalies=anomaly_certificate(h,t)))

    # Strong counterexample requested by the shared proof audit. Both theories
    # have all nine singular values identical; the upper weak components differ.
    rho=F(1,100)
    coefficients={
        "AA_12":12*rho,"AA_33":4*rho,
        "AB_31":rho,"AB_13":rho,"BB_11":4*rho,
        "EE_bare_11":4*rho,"AE_23":rho,
        "BE_11":rho,"BE_22":17*rho,"BE_33":12*rho,
    }
    nonstandard=np.zeros((9,9),dtype=complex)
    exact_nonstandard=[[F(0) for _ in range(9)] for _ in range(9)]
    for key,value in coefficients.items():
        item=by_key[key]
        assert charge(item,F(1,2),F(1))==0,key
        basis=mass_basis(item)
        nonstandard += float(value)*basis
        for i in range(9):
            for j in range(9):
                exact_nonstandard[i][j]+=value*F(int(round(basis[i,j].real)))
    neutral=np.array([[4.,1.,0.],[1.,4.,1.],[0.,1.,4.]])/100
    cross=np.array([[12.,1.],[1.,12.]])/100
    target=np.zeros((9,9),dtype=complex)
    target[np.ix_([2,3,6],[2,3,6])]=neutral
    target[np.ix_([0,8],[1,5])]=cross
    target[np.ix_([1,5],[0,8])]=cross.T
    target[4,7]=target[7,4]=.17
    assert opnorm(nonstandard-target)<1e-14

    standard_coefficients={
        "AA_11":4*rho,"AA_12":rho,"AA_22":4*rho,
        "AA_23":rho,"AA_33":4*rho,
        "BE_11":11*rho,"BE_22":13*rho,"BE_33":17*rho,
    }
    standard=np.zeros((9,9),dtype=complex)
    for key,value in standard_coefficients.items():
        assert charge(by_key[key],F(1,2),F(0))==0
        standard+=float(value)*mass_basis(by_key[key])
    expected_masses=np.sort(np.array([4-math.sqrt(2),4.,4+math.sqrt(2),11.,11.,13.,13.,17.,17.])/100)
    standard_masses=np.sort(np.linalg.svd(standard,compute_uv=False))
    nonstandard_masses=np.sort(np.linalg.svd(nonstandard,compute_uv=False))
    assert np.max(np.abs(standard_masses-expected_masses))<1e-13
    assert np.max(np.abs(nonstandard_masses-expected_masses))<1e-13
    multiplicities=[]
    for mass in expected_masses:
        if not multiplicities or abs(mass-multiplicities[-1][0])>1e-12:
            multiplicities.append([float(mass),1])
        else:
            multiplicities[-1][1]+=1
    simple_positive=sum(mass>0 and multiplicity==1 for mass,multiplicity in multiplicities)
    assert simple_positive==3
    qns=[F(1),F(-1),F(0),F(0),F(-2),F(-1),F(0),F(2),F(1)]
    qs=[F(0)]*3+[F(-1)]*3+[F(1)]*3
    exact_ward=[[ (qns[i]+qns[j])*exact_nonstandard[i][j] for j in range(9)] for i in range(9)]
    assert all(x==0 for row in exact_ward for x in row)
    assert opnorm(np.diag([float(q) for q in qs])@standard+standard@np.diag([float(q) for q in qs]))<1e-13
    pn=np.diag([1.,1.,1.,0.,0.,0.,0.,0.,0.])
    gn=nonstandard.conj().T@nonstandard
    gs=standard.conj().T@standard
    leakage_mass=opnorm(pn@nonstandard@(np.eye(9)-pn))
    leakage_gram=opnorm(pn@gn-gn@pn)
    assert leakage_mass>0 and leakage_gram>0
    assert opnorm(pn@gs-gs@pn)<1e-14
    nonstandard_upper=np.sort(np.linalg.svd(nonstandard[:3,:3],compute_uv=False))
    assert np.max(np.abs(nonstandard_upper-np.array([.04,.12,.12])))<1e-14
    assert opnorm(np.diag([.01,.17,.12]))>0
    exact_gram= [[sum(exact_nonstandard[k][i]*exact_nonstandard[k][j] for k in range(9))
                 for j in range(9)] for i in range(9)]
    gram_cross_entries={f"{LABELS[i]}:{LABELS[j]}":str(exact_gram[i][j])
                        for i in range(3) for j in range(3,9) if exact_gram[i][j]}
    assert gram_cross_entries

    sources=[BASE/"archive_1009_/research_note_1010.md",
             BASE/"archive_1009_/research_note_1011.md",
             BASE/"archive_1009_/1011/next_selection_audit.md",
             BASE/"archive_1009_/1010/hypercharge_neutrino_selection_results.json",
             BASE/"archive_1009_/1011/flavor_charge_selection_results.json"]
    return dict(round=1012,new_calibration_groups=1,cumulative_test_groups=3790,
        scope="two left-handed Weyl fields, one Higgs doublet, no derivatives, scalar insertions <=2, direct masses at fixed nonzero Higgs vev",
        all_scientific_calibrations_passed=True,
        complete_general_EFT_basis_claimed=False,physical_spectrum_fit_claimed=False,
        physical_UV_completion_claimed=False,new_cognitive_axioms=0,
        lepton_field_order=LABELS,
        operator_count_before_U1=54,independent_su2_invariant_dimensions=invariant_dimensions,
        exact_charge_polynomials=charges,exceptional_r=[str(x) for x in sorted(exceptional)],
        parameter_branches=rows,
        su2_grassmann_calibration=dict(rotations=8,operators_per_rotation=54,
            maximum_invariance_residual=max_su2_error,
            maximum_mixed_contraction_identity_residual=max_mixed_identity,
            maximum_lorentz_bilinear_symmetry_residual=max_symmetric_bilinear),
        mass_operator_normalization="1/2 diagonal Weyl bilinear; offdiagonal coefficients equal symmetric mass matrix entries",
        isospectral_counterexample=dict(h="1/2",t_nonstandard="1",r_nonstandard="2",v=1.,
            nonstandard_operator_coefficients={k:str(v) for k,v in coefficients.items()},
            standard_operator_coefficients={k:str(v) for k,v in standard_coefficients.items()},
            nonstandard_electric_charges=[str(q) for q in qns],
            standard_electric_charges=[str(q) for q in qs],
            nonstandard_mass_matrix=serialize_matrix(nonstandard),
            standard_mass_matrix=serialize_matrix(standard),
            nonstandard_masses=nonstandard_masses.tolist(),standard_masses=standard_masses.tolist(),
            expected_nine_singular_values=expected_masses.tolist(),
            distinct_masses_and_multiplicities=multiplicities,
            simple_positive_mass_count=simple_positive,
            maximum_isospectral_residual=float(np.max(np.abs(nonstandard_masses-standard_masses))),
            nonstandard_upper_block_singular_values=nonstandard_upper.tolist(),
            standard_upper_block_singular_values=np.sort(np.linalg.svd(standard[:3,:3],compute_uv=False)).tolist(),
            nonstandard_upper_mass_mixing_norm=leakage_mass,
            nonstandard_upper_gram_commutator_norm=leakage_gram,
            exact_nonstandard_upper_gram_cross_entries=gram_cross_entries,
            original_BE_yukawa_diagonal=[.01,.17,.12],original_BE_yukawa_full_rank=True,
            nonstandard_anomalies=anomaly_certificate(F(1,2),F(1)),
            standard_anomalies=anomaly_certificate(F(1,2),F(0)),
            mass_ward_identity_exact=True,
            same_weak_currents_or_gauge_charges_claimed=False,
            upper_block_closure_derived_from_full_spectrum=False),
        analytic_obligations="representation decomposition proves contraction completeness; nonzero vev polynomial lifts prove Ward-kernel surjectivity; analytic charge equations cover all real r",
        historical_source_sha256={str(path.relative_to(BASE)).replace("\\","/"):hashlib.sha256(path.read_bytes()).hexdigest() for path in sources})


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--write",action="store_true")
    args=parser.parse_args()
    result=run()
    if args.write:
        with OUT.open("x",encoding="utf8") as stream:
            json.dump(result,stream,ensure_ascii=False,indent=2)
            stream.write("\n")
    else:
        compare(result,json.loads(OUT.read_text("utf8")))
    print(json.dumps({"round":result["round"],
        "passed":result["all_scientific_calibrations_passed"],
        "operator_count_before_U1":54,"exceptional_r":result["exceptional_r"],
        "branch_dimensions":[[row["r"],row["allowed_operator_count"],row["independent_mass_matrix_dimension"]] for row in result["parameter_branches"]],
        "su2":result["su2_grassmann_calibration"],
        "isospectral_residual":result["isospectral_counterexample"]["maximum_isospectral_residual"],
        "nonstandard_upper_mass_mixing_norm":result["isospectral_counterexample"]["nonstandard_upper_mass_mixing_norm"],
        "nonstandard_upper_gram_commutator_norm":result["isospectral_counterexample"]["nonstandard_upper_gram_commutator_norm"]},ensure_ascii=False,indent=2))
