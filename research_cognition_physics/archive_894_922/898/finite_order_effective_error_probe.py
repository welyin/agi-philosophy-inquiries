"""898 working: actual854 Q_eff versus its retained E coefficients.
This is NOT an original-Q/E bridge and not a remainder theorem for E(t).
"""
from pathlib import Path
import json,sys,math,itertools
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'854'))
import finite_task_jet_bridge as old
TARGET=HERE/'finite_order_effective_error_probe_results.json'
D=old.D;uv=((0,0),(1,0),(0,1),(2,0),(1,1),(0,2))

def mul(a,b):
    out={}
    for i,x in a.items():
        for j,y in b.items():
            k=tuple(v+w for v,w in zip(i,j))
            if sum(k)<=2:out[k]=out.get(k,0)+x@y
    return out

def exponential_at(h,t):
    d=next(iter(h.values())).shape[0];gen={}
    for k,a in h.items():
        if k[1]+k[2]<=2:gen[k[1:]]=gen.get(k[1:],0)+1j*t**k[0]*a
    out={(0,0):np.eye(d,dtype=complex)};term=out
    for n in range(1,51):
        term={k:a/n for k,a in mul(term,gen).items()}
        for k,a in term.items():out[k]=out.get(k,0)+a
    return out

def models():
    # Identical coefficients and word-closed projection as frozen854.
    ne=12;nr=D+1;q=np.zeros((ne,ne),complex)
    for n in range(ne-1):q[n,n+1]=q[n+1,n]=n+1
    num=np.diag(np.arange(ne,dtype=float));eye=np.eye(ne)
    h1={(1,0,0):np.kron(old.X,q)/3,(2,0,0):np.kron(old.Z,num)/5,
        (2,1,0):np.kron(old.Z,q)/7,(6,2,0):np.kron(old.Y,eye)/11}
    h2={(1,0,0):np.kron(old.Z,q)/4,(2,0,0):np.kron(old.Y,q@q)/6,
        (2,0,1):np.kron(old.X,q)/9,(6,1,1):np.kron(old.Z,eye)/13}
    embedding=np.kron(old.I,np.eye(ne)[:,:nr]);w=np.kron(old.I,np.eye(nr)[:,[0]])
    hs=[{k:embedding.conj().T@a@embedding for k,a in h.items()} for h in (h1,h2)]
    return nr,w,hs

def histories(series,nr,w,time_graded=False):
    out={}
    for a,b in itertools.product((-1,1),repeat=2):
        k1=np.kron((old.I+a*old.Z)/2,np.eye(nr));k2=np.kron((old.I+b*old.X)/2,np.eye(nr))
        a1={k:k1@v for k,v in series[0].items()}
        product=(old.pmul if time_graded else mul)(series[1],a1)
        out[a,b]={k:k2@v@w for k,v in product.items()}
    return out

def run():
    nr,w,hs=models();poly=histories([old.exp_series(h) for h in hs],nr,w,True)
    # Bilateral source coefficients of each record effect on arbitrary input.
    jets={}
    for outcome,p in poly.items():
        jets[outcome]={}
        for ka,a in p.items():
            for kb,b in p.items():
                n=ka[0]+kb[0]
                if n<=D:
                    key=(n,*ka[1:],*kb[1:]);jets[outcome][key]=jets[outcome].get(key,0)+a.conj().T@b
    radius=.15;source_radius=.2
    M=sum(radius**k[0]*source_radius**sum(k[1:])*float(np.linalg.norm(a,2)) for h in hs for k,a in h.items())
    majorant=float(np.exp(2*M))
    alphas=((0,0,0,0),(1,0,0,0),(0,0,1,0),(1,0,0,1),(0,0,2,0))
    rows=[]
    for t in (.005,.01,.02):
        actual=histories([exponential_at(h,t) for h in hs],nr,w)
        sourcefree=histories([{(0,0):old.exact_u(h,t,0.,0.)} for h in hs],nr,w)
        assert max(np.max(abs(actual[r][0,0]-sourcefree[r][0,0])) for r in actual)<2e-14
        for alpha in alphas:
            k=sum(alpha);factor=math.prod(math.factorial(i) for i in alpha)
            bound=factor*majorant*(t/radius)**(D+1)/(1-t/radius)/source_radius**k/t**k
            errors=[]
            for outcome,p in actual.items():
                l=p.get(alpha[:2],np.zeros_like(w));rr=p.get(alpha[2:],np.zeros_like(w))
                exact=l.conj().T@rr
                truncated=sum(t**key[0]*a for key,a in jets[outcome].items() if key[1:]==alpha)
                err=factor*float(np.linalg.norm(exact-truncated,2))/t**k
                # Floating evaluation error is separate from analytic truncation.
                assert err<=bound+3e-13/t**k
                errors.append(err)
            rows.append(dict(parameter=t,bilateral_source_derivative=list(alpha),physical_insertions=k,
                maximum_record_effect_error=max(errors),analytic_each_record_bound=bound,
                sum_record_effect_error_bound=4*bound))
    return dict(working_round=898,status='working_not_formal',formal_round=897,cumulative_numbered_groups=3682,
        fresh_numbered_groups=0,model='Frozen854 two-stage ladder Q_eff; diagnostic, not original Q or full continuum E',
        retained_time_degree=D,finite_dimension=2*nr,complex_time_radius=radius,complex_each_source_radius=source_radius,
        sum_generator_norm_majorant=M,common_bilateral_majorant=majorant,rows=rows,
        comparison_is_actual_Q_eff_to_retained_E_polynomial=True,
        original_Q_mapping_proved=False,exact_uncut_E_probability_assumed=False,
        full_goal_completed=False,floating_point_not_interval_arithmetic=True)

if __name__=='__main__':
    result=run();assert not TARGET.exists()
    TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
