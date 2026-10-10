"""Round 1077: exact source-reader bridge; stdout JSON only.
Original 480 preparation, ideal pulse as certificate intermediary, then
a Duhamel error returns the statement to the finite H-on pulse.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib,importlib.util,itertools,json,math
import numpy as np

BASE=Path(__file__).resolve().parent
def load_reader():
    path=BASE.parent/"1076/check.py"
    spec=importlib.util.spec_from_file_location("frozen1076",path)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    return mod

def initial_columns():
    real=np.zeros((384,28),dtype=object);imag=real.copy();col=0
    for source,graphs in ((0,[None]),(1,list(range(6)))):
        for d0,d2 in itertools.product(range(2),repeat=2):
            for gmode in graphs:
                for q,d3,d4 in itertools.product(range(2),repeat=3):
                    data=sum(v<<(5-j) for j,v in enumerate((d0,q,d2,d3,d4,0)))
                    phase=(d4+source*q)%4
                    ar,ai=((1,0),(0,1),(-1,0),(0,-1))[phase]
                    for g in (range(6) if gmode is None else (gmode,)):
                        real[6*data+g,col]=ar;imag[6*data+g,col]=ai
                col+=1
    assert col==28 and np.sum(real*real+imag*imag)==384
    return real,imag

def rational_source(reader,degree=40):
    choices,edges,rows,h=reader.system()
    den=5**degree*math.factorial(degree)
    swap=[]
    for data in range(64):
        new=data^24 if ((data>>4)&1)!=((data>>3)&1) else data
        for g in range(6):swap.append(6*new+g)
    swap=np.array(swap)
    def ham(v):return np.sum(v[0][rows],axis=1),np.sum(v[1][rows],axis=1)
    def mi(v):return v[1],-v[0]
    def sk(v):return v[0][swap],v[1][swap]
    def prop(v):
        r,i=v
        out_r=np.zeros(r.shape,dtype=object);out_i=out_r.copy()
        for k in range(degree+1):
            coeff=5**(degree-k)*(math.factorial(degree)//math.factorial(k))
            ar,ai=((r,i),(i,-r),(-r,-i),(-i,r))[k%4]
            out_r+=coeff*ar;out_i+=coeff*ai
            if k<degree:r,i=ham((r,i))
        return out_r,out_i
    start=prop(initial_columns())
    y=prop(prop(sk(start)))
    dt=mi(prop(prop(sk(ham(start)))))
    du=mi(ham(y))
    da=mi(prop(prop(start)))
    yr,yi=(v.reshape(64,6,28) for v in y)
    gr=np.zeros((6,6),dtype=object);gi=gr.copy()
    for a,b in zip(yr,yi):
        gr+=a@a.T+b@b.T;gi+=b@a.T-a@b.T
    deriv=[]
    for dr,di in (dt,du,da):
        rr=np.zeros((6,6),dtype=object);ii=rr.copy()
        for a,b,c,d in zip(dr.reshape(64,6,28),di.reshape(64,6,28),yr,yi):
            rr+=a@c.T+b@d.T+c@a.T+d@b.T
            ii+=b@c.T-a@d.T+d@a.T-c@b.T
        assert np.array_equal(rr,rr.T) and np.array_equal(ii,-ii.T)
        deriv.append((rr,ii))
    scale=384*den**6
    a=F(9,5)
    tail=a**(degree+1)/math.factorial(degree+1)/(1-a/F(degree+2))
    werr=(1+tail)**3-1
    return choices,(gr,gi),deriv,scale,tail,werr

def run():
    reader=load_reader()
    choices,g,dg,gs,tail,werr=rational_source(reader)
    _,_,_,_,_,_,nerr,raw=reader.exact_reader()
    nr,ni,nscale=raw
    br=[nr[:6,:6]+nr[6:,6:],nr[:6,6:]+nr[6:,:6],
        ni[6:,:6]-ni[:6,6:],nr[:6,:6]-nr[6:,6:]]
    bi=[ni[:6,:6]+ni[6:,6:],ni[:6,6:]+ni[6:,:6],
        nr[:6,6:]-nr[6:,:6],ni[:6,:6]-ni[6:,6:]]
    bs=2*nscale
    def contract(r,i,a,b):
        return F(int(np.sum(r*a.T-i*b.T)),gs*bs)
    jac=np.array([[contract(r,i,a,b) for a,b in dg] for r,i in zip(br,bi)],dtype=object)
    effect=[contract(r,i,g[0],g[1]) for r,i in zip(br,bi)]
    inverse=reader.inverse(jac[1:]);invnorm=reader.maxrow(inverse)
    tau=F(1,10**10)
    # Bound all four rows, including B0. No missing cross-error term.
    errors=[2*m*nerr+(1+nerr)*2*m*werr*(2+werr)+36*m*tau for m in (9,9,1)]
    perturb=sum(errors)*invnorm
    assert perturb<F(1,1000)
    upper=invnorm/(1-perturb)
    assert upper<2000
    radius=F(1,246240000)
    lower=F(1,2000)-2052*radius
    assert 2000*2052*radius==F(1,60)
    assert lower==F(59,120000) and lower/38==F(59,4560000)
    q=np.array([[int(0 in s)-int(leaf in s) for s in choices] for leaf in (3,4,5)])
    qjac=np.array([[F(int(sum(int(qj[k])*a[k,k] for k in range(6))),gs) for a,b in dg] for qj in q],dtype=object)
    qvalue=[F(int(sum(int(qj[k])*g[0][k,k] for k in range(6))),gs) for qj in q]
    # Existing 480 theorem already certifies the finite-pulse q chart.
    # This comparison is diagnostic; do not count its theorem anew.
    change=jac[1:]@reader.inverse(qjac)
    graph=np.array([[complex(float(F(int(g[0][i,j]),gs)),float(F(int(g[1][i,j]),gs))) for j in range(6)] for i in range(6)])
    assert abs(np.trace(graph)-1)<1e-13 and np.min(np.linalg.eigvalsh(graph))>-1e-12
    return {
      "round":1077,"status":"strict_certificate_passed",
      "source":"unchanged 480/481 correlated seed and fixed local control context; full graph coherence retained",
      "point":["1/5","2/5","pi/2"],"finite_pulse_tau":"1/10000000000",
      "reader_time":"7/10","source_Taylor_degree":40,
      "source_unitary_tail_bound":str(tail),"source_W_error_bound":str(werr),
      "reader_effect_error_bound":str(nerr),
      "ideal_pulse_effect_coefficients_diagnostic":[float(v) for v in effect],
      "ideal_pulse_effect_jacobian_4x3_diagnostic":np.asarray(jac,float).tolist(),
      "finite_pulse_jacobian_entry_error_by_column":[str(e) for e in errors],
      "finite_pulse_jacobian_entry_error_by_column_float":[float(e) for e in errors],
      "bloch_inverse_infinity_certified_upper":2000,
      "reused_480_control_radius":str(radius),
      "uniform_bloch_per_parameter_lower":str(lower),
      "uniform_bloch_per_original_Q_lower":str(lower/38),
      "bloch_inverse_infinity_bound_diagnostic":float(upper),
      "neumann_perturbation_bound":str(perturb),
      "neumann_perturbation_float":float(perturb),
      "ideal_pulse_bloch_determinant_diagnostic":float(reader.det(jac[1:])),
      "ideal_pulse_bloch_singular_values_diagnostic":np.linalg.svd(np.asarray(jac[1:],float),compute_uv=False).tolist(),
      "ideal_pulse_original_Q_diagnostic":[float(v) for v in qvalue],
      "ideal_pulse_original_Q_jacobian_diagnostic":np.asarray(qjac,float).tolist(),
      "ideal_pulse_bloch_per_original_Q_diagnostic":np.asarray(change,float).tolist(),
      "retained_graph_coherence_frobenius_diagnostic":float(np.linalg.norm(graph-np.diag(np.diag(graph)))),
      "rational_jacobian_sha256":hashlib.sha256(json.dumps([[str(v) for v in row] for row in jac]).encode()).hexdigest(),
      "frozen_reader_code_sha256":hashlib.sha256((BASE.parent/"1076/check.py").read_bytes()).hexdigest(),
      "scope":{"actual_source_bridge":True,"finite_H_on_pulse_certified":True,
        "graph_dephasing_used":False,"source_coordinate_to_Pauli_feedback":False,
        "all_equal_Q_histories_same_effect_claimed":False,
        "actual_spatial_dimension_derived":False,"full_reorientation_generated":False,
        "source_preparation_and_control_are_inputs":True,
        "empirical_results":0,"new_accepted_cognitive_axioms":0}
    }
if __name__=="__main__":
    print(json.dumps(run(),ensure_ascii=False,indent=2))
