"""Independent round-1216 replay. Does not import author/earlier programs.
Rebuilds the fixed-point matrices from original h and the old exact cold source,
then checks supplied trial inverses and Choi direction by integer arithmetic.
NumPy diagonalization is a separate diagnostic, not the strict certificate.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse, hashlib, json, math
import numpy as np

HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=STAGE.parent
D=2**80

def ensure(c,label):
    if not bool(c): raise AssertionError(label)

def close(a,b,label,tol=3e-9):
    e=float(np.max(np.abs(np.asarray(a)-np.asarray(b))))
    ensure(e<tol,label+": "+str(e))
    return e

def fstr(x): return str(F(x))
def scalar(x): return round(float(x),13)
def zero(n): return np.zeros((n,n),dtype=object)
def ident(n): return np.eye(n,dtype=object)

def nearest_integer(n,d):
    q,r=divmod(int(n),int(d))
    return q+int(2*r>=d)

def quantize_fraction(x,scale=D):
    q=F(x)
    return nearest_integer(q.numerator*scale,q.denominator)

def quantize_array(a,den):
    return np.array([nearest_integer(x,den) for x in a.flat],dtype=object).reshape(a.shape)

def rawmul(a,b):
    return (a[0]@b[0]-a[1]@b[1],a[0]@b[1]+a[1]@b[0])

def fixedmul(a,b,divisor=D):
    real,imag=rawmul(a,b)
    return quantize_array(real,divisor),quantize_array(imag,divisor)

def adj(a): return a[0].T,-a[1].T
def decode(a): return np.array(a["real"],dtype=object),np.array(a["imag"],dtype=object)
def numeric(a,scale): return np.asarray(a[0],float)/scale+1j*np.asarray(a[1],float)/scale
def sqnorm(a,den):
    return F(sum(int(x)**2 for x in a[0].flat)+sum(int(x)**2 for x in a[1].flat),den*den)
def csubtract(a,b): return a[0]-b[0],a[1]-b[1]

def rational_reference(cold):
    pops=[F(x) for x in cold["one_marginal_rational"]]
    ensure(sum(pops)==1 and min(pops)>0,"original exact cold source")
    rootden=int(cold["preparation"]["root_denominator"])
    boxes=[(F(a,rootden),F(b,rootden)) for a,b in cold["preparation"]["root_boxes"]]
    roots=[(lo+hi)/2 for lo,hi in boxes]
    n=math.isqrt(3*D*D)
    ensure(n*n<=3*D*D<(n+1)*(n+1),"sqrt3 interval")
    h=np.array([[F(0),F(-1,5),F(0),-F(n,80*D)],
                [F(-1,5),F(1),F(0),F(0)],
                [F(0),F(0),F(1),F(0)],
                [-F(n,80*D),F(0),F(0),F(-1,40)]],dtype=object)
    b=zero(4)
    for i,p in enumerate(pops):
        projector=ident(4)
        for j in range(4):
            if j!=i:
                projector=projector@((h-roots[j]*ident(4))/(roots[i]-roots[j]))
        b+=p*projector
    bfixed=np.array([quantize_fraction(x) for x in b.flat],dtype=object).reshape(4,4)
    gap=F(38,1000);width=F(1094,1000);eroot=F(1,10**18);eh=F(1,80*D)
    ensure(min(boxes[j][0]-boxes[i][1] for i in range(4) for j in range(i+1,4))>gap,"root gap")
    ensure(boxes[3][1]-boxes[0][0]<width,"root width")
    factor=(eh+eroot)/(gap-2*eroot)+width*2*eroot/(gap*(gap-2*eroot))
    ensure(factor<F(2,10**15),"factor perturbation")
    ensure((width+eh+eroot)/(gap-2*eroot)<30,"factor norms")
    eb=4*3*30**2*F(2,10**15)+F(16,D)
    ensure(eb<F(3,10**11),"reference trace error")
    Q=zero(4);Q[1,2]=Q[2,1]=1
    hfull=np.kron(h,ident(4))+np.kron(ident(4),h)+F(1,10)*np.kron(Q,Q)
    return pops,bfixed,hfull,eb,factor

def integer_exponential(h,t):
    # A uses nearest rounding; the /k is merged with the single matrix rounding.
    ai=np.array([quantize_fraction(-F(t,4096)*x) for x in h.flat],dtype=object).reshape(16,16)
    ensure(np.array_equal(ai,ai.T),"exact anti-Hermiticity")
    anorm=F(max(sum(abs(int(x)) for x in row) for row in ai),D)
    ensure(anorm<=F(1,2),"scaled input")
    A=zero(16),ai
    term=ident(16)*D,zero(16)
    total=term[0].copy(),term[1].copy()
    for k in range(1,31):
        term=fixedmul(term,A,D*k)
        total=total[0]+term[0],total[1]+term[1]
    tail=F(2,2**31*math.factorial(31))
    error=F(960,D)+tail
    bounds=[error]
    for j in range(12):
        ensure(error<F(1,1000),"squaring premise")
        total=fixedmul(total,total)
        error=F(2001,1000)*error+F(16,D)
        bounds.append(error)
    phase=F(2*t,80*D)+F(4096*16,D)
    ensure(error+phase<F(1,10**12),"full exponential error")
    return total,{"t":t,"scaled_norm":str(anorm),"tail":str(tail),
                  "squaring_bounds":[str(x) for x in bounds],
                  "phase_error":str(phase),"total_error":str(error+phase)}

def integer_channel(u,b):
    sr,si=zero(16),zero(16)
    for i in range(4):
        for j in range(4):
            initial=zero(16)
            initial[4*i:4*i+4,4*j:4*j+4]=b
            out=fixedmul(fixedmul(u,(initial,zero(16))),adj(u))
            for a in range(4):
                for c in range(4):
                    sr[a+4*c,i+4*j]=sum(out[0][4*a+k,4*c+k] for k in range(4))
                    si[a+4*c,i+4*j]=sum(out[1][4*a+k,4*c+k] for k in range(4))
    return sr,si

def choi(a):
    n=16
    r,i=zero(n),zero(n)
    for p in range(4):
        for q in range(4):
            for x in range(4):
                for y in range(4):
                    r[p*4+x,q*4+y]=a[0][x+4*y,p+4*q]
                    i[p*4+x,q*4+y]=a[1][x+4*y,p+4*q]
    return r,i

def independent_float_maps(pops):
    h=np.array([[0,-.2,0,-np.sqrt(3)/80],[-.2,1,0,0],
                [0,0,1,0],[-np.sqrt(3)/80,0,0,-.025]])
    e,basis=np.linalg.eigh(h)
    rhoB=(basis*np.array([float(x) for x in pops]))@basis.T
    Q=np.zeros((4,4));Q[1,2]=Q[2,1]=1
    H=np.kron(h,np.eye(4))+np.kron(np.eye(4),h)+.1*np.kron(Q,Q)
    eig,U=np.linalg.eigh(H)
    maps=[]
    for t in (300,600):
        prop=(U*np.exp(-1j*t*eig))@U.conj().T
        mat=np.zeros((16,16),complex)
        for i in range(4):
            for j in range(4):
                E=np.zeros((4,4));E[i,j]=1
                state=prop@np.kron(E,rhoB)@prop.conj().T
                output=np.trace(state.reshape(4,4,4,4),axis1=1,axis2=3)
                mat[:,i+4*j]=output.reshape(-1,order="F")
        maps.append(mat)
    S,T=maps
    inv=np.linalg.inv(S);mid=T@inv
    C=np.zeros((16,16),complex)
    for i in range(4):
        for j in range(4):
            C[4*i:4*i+4,4*j:4*j+4]=mid[:,i+4*j].reshape(4,4,order="F")/4
    # Both physical endpoint maps have the same input free evolution.
    ef,Vf=np.linalg.eigh(h)
    free=(Vf*np.exp(-1j*ef))@Vf.conj().T
    Afree=np.kron(free.conj(),free)
    cancellation=close((T@Afree)@np.linalg.inv(S@Afree),mid,"common free input phase cancels")
    return S,T,{"smin_S":scalar(np.linalg.svd(S,compute_uv=False)[-1]),
               "inverse_F":scalar(np.linalg.norm(inv)),
               "intermediate_F":scalar(np.linalg.norm(mid)),
               "minimum_normalized_Choi_eigenvalue":scalar(np.linalg.eigvalsh((C+C.conj().T)/2)[0]),
               "input_free_phase_cancellation_error":scalar(cancellation)}

def compute():
    author=json.loads((HERE/"results.json").read_text(encoding="utf-8"))
    cold=json.loads((STAGE/"1202/results.json").read_text(encoding="utf-8-sig"))
    clock=json.loads((STAGE/"1203/results.json").read_text(encoding="utf-8-sig"))
    cert=author["exact_matrix_certificate"]; ag=author["groups"]
    ensure(int(cert["S_That_denominator"])==D,"channel dyadic scale")
    dx=int(cert["X_Y_denominator"]);ensure(dx==2**40,"trial scale")
    S,T,X,Y=[decode(cert[k]) for k in ("Shat","That","X","Y")]

    # Group 1: independently reconstruct the entire fixed-point certificate.
    for den in (2,3,4,7,D):
        for n in range(-25,26):
            q=nearest_integer(n,den)
            ensure(abs(F(q)-F(n,den))<=F(1,2),"nearest round proof checks")
    pops,b,H,eb,factor=rational_reference(cold)
    ensure(np.array_equal(b,np.array(ag["native_reference_polynomial_certificate"]["reference_fixed_real"],dtype=object)),
           "exact original reference reconstruction")
    expdata=[]
    for t,expected in ((300,S),(600,T)):
        U,bounds=integer_exponential(H,t)
        rebuilt=integer_channel(U,b)
        ensure(np.array_equal(rebuilt[0],expected[0]) and np.array_equal(rebuilt[1],expected[1]),
               "complete integer channel reconstruction "+str(t))
        expdata.append(bounds)
    g1={"name":"independent_full_integer_reconstruction","passed":True,
        "rounding":"nearest; Taylor product/division rounded together",
        "exact_ShThat_match":True,"reference_trace_error":str(eb),
        "factor_error":str(factor),"exponentials":expdata}

    # Group 2: separately diagonalize the physical matrices, in the original site basis.
    Sf,Tf,diag=independent_float_maps(pops)
    diag["S_integer_vs_independent_float"]=scalar(close(numeric(S,D),Sf,"S independent float"))
    diag["T_integer_vs_independent_float"]=scalar(close(numeric(T,D),Tf,"T independent float"))
    g2={"name":"independent_original_site_quantum_channel_diagnostic","passed":True,**diag,
        "floating_eigenvalues_used_as_strict_proof":False}

    # Group 3: exact residuals for both inverse orders and the actual needed YS-T.
    common=dx*D
    I=ident(16)*common,zero(16)
    rx= sqnorm(csubtract(rawmul(X,S),I),common)
    rxright=sqnorm(csubtract(rawmul(S,X),I),common)
    ry=sqnorm(csubtract(rawmul(Y,S),(T[0]*dx,T[1]*dx)),common)
    nx,ny=sqnorm(X,dx),sqnorm(Y,dx)
    ensure(nx<36 and ny<16,"trial norm certificates")
    ensure(max(rx,rxright,ry)<F(1,10**16),"strict 1e-8 residuals")
    g=ag["exact_inverse_and_intermediate_residuals"]
    for val,key in ((nx,"X_Frobenius_squared"),(ny,"Y_Frobenius_squared"),
                    (rx,"X_Shat_minus_I_Frobenius_squared"),
                    (rxright,"Shat_X_minus_I_Frobenius_squared"),
                    (ry,"Y_Shat_minus_That_Frobenius_squared")):
        ensure(val==F(g[key]),"author exact residual crosscheck "+key)
    g3={"name":"exact_trial_norms_and_two_sided_inverse_residuals","passed":True,
        "X_norm_squared":str(nx),"Y_norm_squared":str(ny),
        "X_S_minus_I_squared":str(rx),"S_X_minus_I_squared":str(rxright),
        "Y_S_minus_T_squared":str(ry)}

    # Group 4: exact normalized Choi Rayleigh quotient, retaining complex phases.
    wg=ag["rational_normalized_Choi_witness"]
    zr=np.array(wg["vector_real"],dtype=object)
    zi=np.array(wg["vector_imag"],dtype=object)
    dz=int(wg["vector_denominator"])
    Cr,Ci=choi(Y)
    wr,wi=Cr@zr-Ci@zi,Cr@zi+Ci@zr
    numerator=sum(int(a)*int(b) for a,b in zip(zr,wr))+sum(int(a)*int(b) for a,b in zip(zi,wi))
    zn=sum(int(a)**2+int(b)**2 for a,b in zip(zr,zi))
    ensure(zn>0,"nonzero witness")
    ray=F(numerator,4*dx*zn)
    ensure(ray==F(wg["exact_real_Rayleigh"]) and ray<F(-3,50),"strict rational negative Choi")
    ensure(F(zn,dz*dz)==F(wg["vector_norm_squared"]),"vector normalization")
    g4={"name":"independent_integer_normalized_Choi_witness","passed":True,
        "ordering":"vec column-major; Choi input-major/output; divide by 4",
        "norm_squared":str(F(zn,dz*dz)),"exact_Rayleigh":str(ray)}

    # Group 5: recalculate the analytic numerical enclosure from rational bounds.
    q=F(1,10**8);r=F(1,10**8);eu=F(1,10**12);b_safe=F(3,10**11)
    column=2*eu+eu*eu+b_safe*(1+eu)**2+F(256,D)*(2+eu)
    ensure(4*column<q,"Liouville rounding and state error")
    saved=ag["two_finite_integer_exponential_certificates"]["certificates"]
    for own,old in zip(expdata,saved):
        ensure(F(own["total_error"])==F(old["total_operator_error_bound"]),"exponential bound equality")
        ensure(own["squaring_bounds"]==old["squaring_error_bounds"],"full rational recurrence")
    ensure(4*column==F(ag["reduced_channel_error_certificate"]["Liouville_operator_error_bound"]),
           "channel analytic error equality")
    g5={"name":"analytic_rounding_projector_exponential_and_channel_bounds","passed":True,
        "channel_operator_error":str(4*column),"q_used":str(q),
        "partial_trace_and_nonhermitian_matrix_units_counted":True}

    # Group 6: same v=1 father H at 301 and 601, including all actual J tails.
    epsilon=F(clock["bounds_by_program"]["1"]["isometry_error_upper"])
    ensure(epsilon==F("0.000109379990"),"unchanged physical clock budget")
    mass=F(10**24);sigma=F(1,10**10);L=F(150)+F(3,2)*F(22,7)
    endpoint=[]
    for t in (301,601):
        local=L*t*sigma*(1+F(t*t)/(8*mass*mass*sigma**4))+F(1,10000)
        ensure(local<epsilon and t<603,"same endpoint full gate bound")
        endpoint.append({"t":t,"coarse_local_bound":str(local)})
    xi=4*epsilon+q
    inverse_residual=r+6*xi
    ensure(inverse_residual<F(1,4),"actual invertibility")
    inverse_bound=F(6)/(1-inverse_residual)
    ensure(inverse_bound<8,"actual inverse norm")
    perturbation=8*(r+5*xi)
    ensure(perturbation<F(18,1000),"actual intermediate perturbation")
    ray_actual=ray+perturbation
    ensure(ray_actual<F(-21,500),"actual non-CP margin")
    ensure(ray_actual==F(ag["same_clock_transport_and_actual_non_CP"]["rational_actual_Rayleigh_upper"]),
           "actual signed certificate")
    g6={"name":"same_clock_actual_invertibility_and_non_CP_transport","passed":True,
        "physical_times":[301,601],"program_v":1,"endpoint_bounds":endpoint,
        "epsilon_each_half_diamond":str(epsilon),"Liouville_bound_xi":str(xi),
        "actual_inverse_norm_upper_from_exact_formula":str(inverse_bound),
        "normalized_Choi_perturbation":str(perturbation),
        "actual_Rayleigh_upper":str(ray_actual),
        "strict_common_conclusion":"actual normalized Choi Rayleigh < -0.042"}

    paths=["archive_370_428/research_note_390.md","archive_1201_/1202/results.json",
           "archive_1201_/1203/proof.md","archive_1201_/1203/results.json"]
    for name in ("proof.md","check.py","results.json","NEXT.md"):
        paths.append("archive_1201_/1216/"+name)
    paths.append("archive_1201_/research_note_1216.md")
    hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}
    return {"round":1216,"status":"INDEPENDENT_PASS_SCOPED_NOT_ACCEPTANCE",
            "all_checks_passed":True,"check_groups":6,"groups":[g1,g2,g3,g4,g5,g6],
            "source_path_base":"research_cognition_physics","source_sha256":hashes,
            "scope":{
                "same_v1_fixed_initial_environment":True,"arbitrary_A_reference_independent_of_B_factory":True,
                "actual_copy_Gaussian_tails_retained":True,"initial_free_phase_retained_then_cancelled":True,
                "A_only_two_time_non_CP_divisibility":True,"pair_trace_distance_backflow_claimed":False,
                "P_nondivisibility_claimed":False,"full_AB_non_CP_divisibility_claimed":False,
                "all_time_local_generators_impossible_claimed":False,
                "continuous_clock_grid_simulation":False,"new_intermediate_measurement":False}}

def compare(a,b,path="root"):
    ensure(type(a)==type(b),"saved type "+path)
    if isinstance(a,dict):
        ensure(a.keys()==b.keys(),"saved keys "+path)
        for k in a:compare(a[k],b[k],path+"."+k)
    elif isinstance(a,list):
        ensure(len(a)==len(b),"saved length "+path)
        for i,(x,y) in enumerate(zip(a,b)):compare(x,y,path+"."+str(i))
    elif isinstance(a,float):
        ensure(path.startswith("root.groups.1."),"unexpected floating field "+path)
        ensure(math.isfinite(a) and math.isfinite(b) and abs(a-b)<=1e-12,
               "saved diagnostic "+path)
    else:ensure(a==b,"saved value "+path)

def main():
    p=argparse.ArgumentParser();p.add_argument("--write",action="store_true");args=p.parse_args()
    target=HERE/"independent_results.json"
    if args.write and target.exists():raise SystemExit("Refuse to overwrite independent_results.json")
    out=compute()
    if args.write:
        with target.open("x",encoding="utf-8") as f:
            json.dump(out,f,ensure_ascii=False,indent=2,allow_nan=False);f.write("\n")
        print(json.dumps({"status":"INDEPENDENT_RESULTS_WRITTEN","all_checks_passed":True,"check_groups":6}))
    else:
        if not target.exists():raise SystemExit("Missing independent_results.json; first creation requires --write")
        compare(out,json.loads(target.read_text(encoding="utf-8")))
        print(json.dumps({"status":"PASS","all_checks_passed":True,"results_match":True,
                          "check_groups":6,"mode":"read_only_replay"}))
if __name__=="__main__":main()
