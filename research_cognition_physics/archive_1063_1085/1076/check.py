"""Candidate 1076: exact-integer Taylor certificate for the original fixed reader.
No image generation, no filesystem writes by this program.
"""
from pathlib import Path
from fractions import Fraction as F
import itertools,math,json
import numpy as np

LEAVES=(0,3,4,5)
def system():
    choices=list(itertools.combinations(LEAVES,2))
    edges=[[(1,2)]+[(1 if j in s else 2,j) for j in LEAVES] for s in choices]
    rows=[]
    for data in range(64):
        for g,s in enumerate(choices):
            adj=[6*data+h for h,ss in enumerate(choices) if len(set(s)&set(ss))==1]
            for a,b in edges[g]:
                mask=(1<<(5-a))|(1<<(5-b))
                d=data^mask if ((data>>(5-a))&1)!=((data>>(5-b))&1) else data
                adj.append(6*d+g)
            assert len(adj)==9
            rows.append(adj)
    rows=np.array(rows,dtype=int)
    h=np.zeros((384,384),dtype=np.int64)
    for i,rr in enumerate(rows):
        for j in rr:h[i,j]+=1
    assert np.array_equal(h,h.T) and np.all(h.sum(axis=1)==9)
    return choices,edges,rows,h

def input_columns():
    # For each mixed input (data0,data2), V=C/2 is an isometry QxG -> DG.
    # The two |+X>,|+Y> vectors are unnormalized here; all mixed weights are 1/4.
    cr=np.zeros((384,48),dtype=object);ci=cr.copy()
    for m in range(4):
        d0,d2=m//2,m%2
        for q,g,b3,b4 in itertools.product(range(2),range(6),range(2),range(2)):
            bits=[d0,q,d2,b3,b4,0]
            data=sum(v<<(5-k) for k,v in enumerate(bits))
            if b4:ci[6*data+g,12*m+6*q+g]=1
            else:cr[6*data+g,12*m+6*q+g]=1
    c=cr.astype(float)+1j*ci.astype(float)
    assert np.array_equal(c.conj().T@c,4*np.eye(48))
    return cr,ci

def exact_reader(t=F(7,10),degree=40):
    choices,edges,rows,h=system()
    cr,ci=input_columns()
    den=t.denominator**degree*math.factorial(degree)
    out_r=np.zeros(cr.shape,dtype=object);out_i=out_r.copy()
    vr,vi=cr.copy(),ci.copy()
    for k in range(degree+1):
        coeff=t.numerator**k*t.denominator**(degree-k)*(math.factorial(degree)//math.factorial(k))
        ar,ai=((vr,vi),(vi,-vr),(-vr,-vi),(-vi,vr))[k%4]
        out_r+=coeff*ar;out_i+=coeff*ai
        if k<degree:
            vr=np.sum(vr[rows],axis=1);vi=np.sum(vi[rows],axis=1)
    mask=np.array([((data>>4)&1)==0 for data in range(64) for g in range(6)])
    nr=np.zeros((12,12),dtype=object);ni=nr.copy()
    for m in range(4):
        rr=out_r[mask,12*m:12*(m+1)];ii=out_i[mask,12*m:12*(m+1)]
        nr+=rr.T@rr+ii.T@ii
        ni+=rr.T@ii-ii.T@rr
    common=16*den*den
    assert np.array_equal(nr,nr.T) and np.array_equal(ni,-ni.T)
    # Pauli coefficients B_mu in N=sum sigma_mu tensor B_mu; coefficient means no extra /2.
    br=[nr[:6,:6]+nr[6:,6:], nr[:6,6:]+nr[6:,:6],
        ni[6:,:6]-ni[:6,6:], nr[:6,:6]-nr[6:,6:]]
    bi=[ni[:6,:6]+ni[6:,6:], ni[:6,6:]+ni[6:,:6],
        nr[:6,6:]-nr[6:,:6], ni[:6,:6]-ni[6:,6:]]
    bd=2*common
    for rr,ii in zip(br,bi):
        assert np.array_equal(rr,rr.T) and np.array_equal(ii,-ii.T)
    # Coordinates x_k=p_k-p_(5-k); fixed each complementary pair total 1/3.
    response=[[F(int(rr[k,k]-rr[5-k,5-k]),2*bd) for k in range(3)] for rr in br]
    # W=(|0><1|+|1><0|)/2 has trace norm 1 and zero population.
    coh=[F(int(rr[0,1]+rr[1,0]),2*bd) for rr in br]
    a=9*abs(t)
    tail=a**(degree+1)/math.factorial(degree+1)/(1-a/F(degree+2))
    # V_m is isometric and averaged, so ||N-Npoly|| <= 2tail+tail².
    effect_error=2*tail+tail*tail
    return choices,rows,h,np.array(response,dtype=object),coh,tail,effect_error,(nr,ni,common)

def det(a):
    a=np.asarray(a,dtype=object);n=len(a)
    if n==1:return a[0,0]
    return sum(((-1)**j)*a[0,j]*det(np.delete(np.delete(a,0,axis=0),j,axis=1)) for j in range(n))

def inverse(a):
    n=len(a);d=det(a);assert d
    cof=np.array([[(-1)**(i+j)*det(np.delete(np.delete(a,i,axis=0),j,axis=1)) for j in range(n)] for i in range(n)],dtype=object)
    inv=cof.T/d
    assert np.array_equal(inv@a,np.eye(n,dtype=object))
    return inv

def maxrow(a):return max(sum(abs(x) for x in row) for row in a)
def serialized_matrix(a):return [[str(x) for x in row] for row in a]
def run():
    choices,rows,h,A,coh,tail,err,raw=exact_reader()
    bloch=A[1:,:]
    inv=inverse(bloch);invbound=maxrow(inv);pert=3*err*invbound
    assert pert<F(1,10**8)
    inv_true_bound=invbound/(1-pert)
    full=np.column_stack((A,np.array(coh,dtype=object)))
    invfull=inverse(full);fullbound=maxrow(invfull);fullpert=4*err*fullbound
    assert fullpert<F(1,10**8)
    # Source I/6 +/- sW is strictly positive for 0<s<1/3.
    s=F(1,10)
    # The Z Bloch difference yields an actual probability gap on at least one +/-Z qubit.
    gap=2*s*(abs(coh[3])-err)
    assert gap>F(48,10000) # > .0048, provided this chosen pair normalization is retained.
    assert inv_true_bound<479
    assert fullbound/(1-fullpert)<71
    assert err<F(1,20000000000000000)
    plus_gap=2*s*(coh[0]+coh[3]-2*err)
    assert plus_gap>F(37,5000)
    # Numerical representation only for independent diagnostics.
    result={
      "candidate_round":1076,"status":"strict_certificate_passed",
      "t":"7/10","Taylor_degree":40,"Hilbert_dimension":384,
      "graph_leaf_choices":[list(x) for x in choices],
      "coordinate_definition":"x_k=p_k-p_(5-k), k=0,1,2; each pair sums to 1/3",
      "effect_convention":"E=b0 I+bx X+by Y+bz Z",
      "Taylor_operator_tail_bound":str(tail),"effect_operator_error_upper":"1/20000000000000000",
      "effect_error_float":float(err),
      "response_matrix_4x3_float":np.asarray(A,dtype=float).tolist(),
      "bloch_singular_values_diagnostic":np.linalg.svd(np.asarray(bloch,dtype=float),compute_uv=False).tolist(),
      "bloch_determinant_sign":1 if det(bloch)>0 else -1,
      "bloch_det_float_diagnostic":float(det(bloch)),
      "bloch_inverse_infinity_certified_upper":479,
      "bloch_inverse_infinity_bound_float":float(inv_true_bound),
      "bloch_invertibility_perturbation":float(pert),
      "response_4x4_with_same_population_coherence_float":np.asarray(full,dtype=float).tolist(),
      "four_coefficient_inverse_infinity_certified_upper":71,
      "four_coefficient_inverse_infinity_bound_float":float(fullbound/(1-fullpert)),
      "four_coefficient_invertibility_perturbation":float(fullpert),
      "same_population_graph_pair_s":"1/10",
      "same_population_at_least_one_of_plus_minus_Z_gap_certified_lower":"48/10000",
      "same_population_at_least_one_of_plus_minus_Z_gap_lower_float":float(gap),
      "same_population_plus_Z_gap_certified_lower":"37/5000",
      "same_population_plus_Z_gap_lower_float":float(plus_gap),
      "rational_response_sha256":__import__("hashlib").sha256(json.dumps(serialized_matrix(full)).encode()).hexdigest(),
      "scope":{"fixed_preparation_section_bloch_rank":3,
        "declared_four_source_directions_effect_rank":4,
        "entire_population_quotient_reader_well_defined":False,
        "spatial_dimension_derived":False,"empirical_results":0,
        "new_coordinate_to_Pauli_feedback":False,"all_ranks_certified_by_exact_arithmetic":True}
    }
    return result

if __name__=="__main__":
    print(json.dumps(run(),ensure_ascii=False,indent=2))
