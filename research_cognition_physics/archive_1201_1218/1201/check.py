"""1201: finite native-contact thermometry and an unused-time energy prediction.

NumPy is diagnostic only. Authoritative K and energy coefficients use 80-digit
Decimal scaled Taylor evaluation, with explicit rational error allowances.
No Gibbs reference is assumed by the unknown-state reconstruction.
"""
from decimal import Decimal as D, localcontext
from fractions import Fraction as F
from pathlib import Path
import json, math
import numpy as np

HERE = Path(__file__).resolve().parent
ZERO = (D(0), D(0))
def z(x=0, y=0): return D(x), D(y)
def add(a,b): return a[0]+b[0], a[1]+b[1]
def mul(a,b): return a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0]
def scale(a,s): return a[0]*s,a[1]*s
def conj(a): return a[0],-a[1]
def zeros(n,m): return [[ZERO for _ in range(m)] for _ in range(n)]
def eye(n): return [[z(int(i==j)) for j in range(n)] for i in range(n)]
def dag(a): return [[conj(a[j][i]) for j in range(len(a))] for i in range(len(a[0]))]
def mm(a,b):
    # This guard also bounds each unrounded product (<100), difference
    # (<200) and partial dot-product sum (<3201); checking only the final
    # answer would not exclude large intermediate cancellations.
    assert all(abs(v)<10 for mat in (a,b) for row in mat for x in row for v in x)
    bt=list(zip(*b)); out=[]
    for row in a:
        line=[]
        for col in bt:
            re=im=D(0)
            for x,y in zip(row,col):
                re += x[0]*y[0]-x[1]*y[1]
                im += x[0]*y[1]+x[1]*y[0]
            assert abs(re)<3201 and abs(im)<3201
            line.append((re,im))
        out.append(line)
    return out
def plus(a,b): return [[add(x,y) for x,y in zip(ar,br)] for ar,br in zip(a,b)]
def sm(a,s): return [[scale(x,s) for x in row] for row in a]
def kron(a,b):
    return [[mul(a[i][j],b[k][l]) for j in range(len(a[0])) for l in range(len(b[0]))]
            for i in range(len(a)) for k in range(len(b))]
def dec(f): return D(f.numerator)/D(f.denominator)
def poly(x): return x**3-F('.975')*x*x-F('.06546875')*x-F('.00053125')
def roots():
    boxes=[]
    for a,b in [('-.054072','-.054071'),('-.009461','-.009460'),('1.038532','1.038533')]:
        a,b=F(a),F(b)
        assert poly(a)*poly(b)<0
        for _ in range(240):
            c=(a+b)/2
            if poly(a)*poly(c)<0:b=c
            else:a=c
        boxes.append((a,b))
    return boxes[:2]+[(F(1),F(1))]+boxes[2:]
def native():
    rt=D(3).sqrt()
    q=zeros(4,4);q[1][2]=q[2][1]=z(1)
    h=[[z(0),z('-0.2'),z(0),z(-rt/80)],
       [z('-0.2'),z(1),z(0),z(0)],
       [z(0),z(0),z(1),z(0)],
       [z(-rt/80),z(0),z(0),z('-.025')]]
    boxes=roots(); e=[dec((a+b)/2) for a,b in boxes]
    R=zeros(4,4)
    for k,x in enumerate(e):
        v=[D(0),D(0),D(1),D(0)] if k==2 else [D(1),D('.2')/(1-x),D(0),-rt/(80*(x+D('.025')))]
        norm=sum(y*y for y in v).sqrt()
        for i in range(4):R[i][k]=z(v[i]/norm)
    # Exact rational residual checks on the stored decimals. Bound sqrt(3)
    # by a rational interval, then bound h*R-R*diag(e) by entry sums.
    s=F(rt); ds=F('1e-78'); assert (s-ds)**2<3<(s+ds)**2
    rf=[[F(x[0]) for x in row] for row in R]
    hf=[[F(x[0]) for x in row] for row in h]
    assert abs(hf[0][3]+s/80)<F('1e-78')
    h_entry_error=F('1e-77')
    res=F(0); gram=F(0)
    for i in range(4):
        for j in range(4):
            rij=sum(hf[i][k]*rf[k][j] for k in range(4))-rf[i][j]*F(e[j])
            res+=abs(rij)+h_entry_error*sum(abs(rf[k][j]) for k in range(4))
            gram+=abs(sum(rf[k][i]*rf[k][j] for k in range(4))-int(i==j))
    assert res<F('1e-65') and gram<F('1e-65')
    # Every distinct spectral separation exceeds .038, by the old rational
    # isolating brackets. Residual/gap + normalization gives ||R-Rexact||<1e-60.
    assert 4*res/F('.038')+2*gram<F('1e-60')
    return h,q,e,R,{'residual_l1_upper':float(res),'gram_l1':float(gram),
                    'basis_operator_error_upper':'1e-60','root_bracket_width_max':float(max(b-a for a,b in boxes))}
def exp_native(H,t):
    n=16;s=14;N=24; unit=F('1e-70')
    # ||H|| <= ||H||_infty <3, so scaled argument norm <1/8 at t<=600.
    assert max(sum(abs(a[0])+abs(a[1]) for a in row) for row in H)<3
    assert 3*t/(2**s)<.125
    X=[[scale(mul(z(0,-1),v),D(t)/D(2**s)) for v in row] for row in H]
    term=eye(n);P=eye(n)
    for k in range(1,N+1):
        term=sm(mm(term,X),D(1)/k);P=plus(P,term)
    # exp tail <=2*(1/8)^25/25!. 1e8*unit covers input construction,
    # 24 polynomial matrix products, additions and divisions.
    err=2*F(1,8)**25/math.factorial(25)+10**8*unit
    for _ in range(s):
        P=mm(P,P)
        assert err<1
        # 2*err+err**2 <= 3*err; this avoids exponentially large rational
        # denominators while keeping a strict conservative certificate.
        err=3*err+4096*unit
    assert err<F('1e-38')
    assert max(abs(a[0])+abs(a[1]) for row in P for a in row)<2
    return P,{'scaling_squarings':s,'taylor_degree':N,'unit_round_absolute_allowance':'1e-70',
              'computed_operator_error_upper':float(err),'declared_operator_error_upper':'1e-38'}
def inv_fraction(a):
    n=len(a);m=[[F(x) for x in row]+[F(int(i==j)) for j in range(n)] for i,row in enumerate(a)]
    for k in range(n):
        pivot=next(j for j in range(k,n) if m[j][k]);m[k],m[pivot]=m[pivot],m[k]
        v=m[k][k];m[k]=[x/v for x in m[k]]
        for j in range(n):
            if j!=k:
                v=m[j][k];m[j]=[x-v*y for x,y in zip(m[j],m[k])]
    return [r[n:] for r in m]
def rounded(x):return str(x.quantize(D('1e-12')))
def contact(t,h,q,e,R):
    H0=plus(kron(h,eye(4)),kron(eye(4),h));V=kron(q,q);H=plus(H0,sm(V,D('.1')))
    U,cert=exp_native(H,t);T=kron(R,R);Ue=mm(mm(dag(T),U),T)
    K=[[D(0)]*4 for _ in range(4)];EA=[];EB=[];W=[]
    for b in range(4):
        probs=[]
        for j in range(4):
            probs.append([sum(x*x for x in Ue[4*j+d][8+b]) for d in range(4)])
            K[j][b]=sum(probs[j])
        assert abs(sum(sum(row) for row in probs)-1)<D('1e-35')
        a=sum(e[j]*probs[j][d] for j in range(4) for d in range(4))-1
        bdelta=sum(e[d]*probs[j][d] for j in range(4) for d in range(4))-e[b]
        vout=mm(mm(dag(U),V),U)
        cb=[T[i][8+b] for i in range(16)]
        vv=ZERO
        for i in range(16):
            for j in range(16):vv=add(vv,mul(mul(conj(cb[i]),vout[i][j]),cb[j]))
        w=-D('.1')*vv[0]
        assert abs(a+bdelta-w)<D('1e-35')
        EA.append(a);EB.append(bdelta);W.append(w)
    # All stored K/energy numbers are rounded to 12 decimals, hence the
    # deliberately loose common interval radius 1e-10 dominates all errors.
    return {'t':t,'K':[[rounded(x) for x in row] for row in K],
            'delta_EA_coefficients':[rounded(x) for x in EA],
            'delta_EB_coefficients':[rounded(x) for x in EB],
            'switch_work_coefficients':[rounded(x) for x in W],
            'entry_error_upper':'1e-10','exponential_certificate':cert}
def entropy(r):
    w=np.linalg.eigvalsh((r+r.conj().T)/2);w=w[w>1e-14]
    return float(-np.sum(w*np.log(w)))
def diagnostics(h,q,e,R,contacts):
    h=np.array([[float(x[0]) for x in row] for row in h]);q=np.array([[float(x[0]) for x in row] for row in q]);R=np.array([[float(x[0]) for x in row] for row in R]);e=np.array([float(x) for x in e])
    H0=np.kron(h,np.eye(4))+np.kron(np.eye(4),h);V=np.kron(q,q);H=H0+.1*V
    ev,Z=np.linalg.eigh(H);T=np.kron(R,R)
    control=np.exp(-1.3*e);control/=sum(control)
    rows={}
    for t in [600,300]:
        U=(Z*np.exp(-1j*t*ev))@Z.T;Ue=T.T@U@T
        K=np.array([[sum(abs(Ue[4*j+d,8+b])**2 for d in range(4)) for b in range(4)] for j in range(4)])
        difference=float(np.max(abs(K-np.array(contacts[str(t)]['K'],dtype=float))));assert difference<1e-10
        vals={}
        for name,r in [('nonthermal',np.array([3,1,1,3])/8),('gibbs_calibration_only',control)]:
            ini=np.kron(np.diag([0,0,1,0]),np.diag(r));out=Ue@ini@Ue.conj().T
            a=np.einsum('abcb->ac',out.reshape(4,4,4,4));b=np.einsum('abad->bd',out.reshape(4,4,4,4))
            p=np.real(np.diag(a)); vals[name]={'reference_population':r.tolist(),'probabilities':p.tolist(),
                'SA_before_pointer':entropy(a),'SB':entropy(b),'SAB':entropy(out),
                'I_AB':entropy(a)+entropy(b)-entropy(out),'record_entropy':float(-sum(p*np.log(p))),
                'delta_EA':float(np.dot(e,np.diag(a).real)-1),'delta_EB':float(np.dot(e,np.diag(b).real-r))}
            vals[name]['switch_work']=vals[name]['delta_EA']+vals[name]['delta_EB']
        rows[str(t)]={'independent_float_K_deviation':difference,'examples':vals}
    # Parity deletes c/outer coherences from this entire probability menu.
    S=np.diag([1,1,-1,1]);assert np.linalg.norm(H@np.kron(S,S)-np.kron(S,S)@H)<1e-14
    return rows
def main():
    with localcontext() as ctx:
        ctx.prec=80
        h,q,e,R,basis_cert=native(); contacts={str(t):contact(t,h,q,e,R) for t in [600,300]}
        K=[[F(x) for x in row] for row in contacts['600']['K']]; L=inv_fraction(K)
        for i in range(4):
            for j in range(4):assert sum(L[i][k]*K[k][j] for k in range(4))==int(i==j)
        ln=max(sum(abs(x) for x in row) for row in L);assert ln<F('3.14')
        matrix_error=F('1e-10'); neumann=ln*4*matrix_error
        assert ln/(1-neumann)<F('3.14')
        prob_error=F('.005');rerror=F('.015701');assert ln*(prob_error+matrix_error)<rerror
        prediction={}
        for key in ['delta_EA_coefficients','delta_EB_coefficients','switch_work_coefficients']:
            coeff=[F(x) for x in contacts['300'][key]]
            prediction[key]={'coefficient_l1':float(sum(abs(x) for x in coeff)),
                'error_upper_from_reconstructed_r':float(sum(abs(x) for x in coeff)*rerror+matrix_error)}
        # The two inequalities r3>r2 and r0>r1 jointly reject every real beta.
        nonthermal_margin=F(1,8)-2*rerror;assert nonthermal_margin>F('.0935')
        diag=diagnostics(h,q,e,R,contacts)
        result={'schema':'finite_native_thermometry_1201_v1','status':'PASS_FINITE_CONDITIONAL_CERTIFICATE',
            'g':'1/10','hbar':1,'native_dimensions':[4,4],'stationarity_is_input':True,
            'basis_certificate':basis_cert,'eigenvalues':[str(x) for x in e], 'contacts':contacts,
            'inverse_inf_norm':float(ln),'certified_actual_inverse_inf_upper':'3.14',
            'inverse_rational':[[str(x) for x in row] for row in L],
            'finite_error_contract':{'probability_error':str(prob_error),'population_error_upper':str(rerror),
                'systematic_probability_budget':'.001','statistical_probability_budget':'.004',
                'iid_fresh_reference_copies':250000,'failure_probability_upper':8*math.exp(-8),
                'nonthermal_distance_from_all_real_beta_gibbs_lower':str(nonthermal_margin),
                'holdout_t300_energy_prediction_errors':prediction},
            'floating_diagnostics_not_interval_certificates':diag,
            'remaining_inputs':['independent stationary B preparation; fresh iid copies for sampling',
                'known c preparation, finite control/timing and energy readout with stated systematic error',
                'controller/battery and blank-record preparation costs not autonomously constructed',
                'actual clock/three-dimensional region identity not established by tensor labels'],
            'not_claimed':['unknown reference is Gibbs','full quantum-state tomography','new entropy law',
                'autonomous internal control','geometry, area entropy, Unruh or GR']}
    (HERE/'results.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':result['status'],'inverse_inf':float(ln),'nonthermal_margin':float(nonthermal_margin),
                      'holdout_errors':prediction},ensure_ascii=False))
if __name__=='__main__':main()
