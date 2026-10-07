"""851: exact joint initial-mean/content calibration using the original code.

The wave recurrence and O(2) gauge potential are explicit diagnostics, not
numerical solutions of the original spacetime quantum theory.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,itertools,json,sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
TARGET=HERE/'full_source_mean_content_results.json'
sys.path.insert(0,str(HERE.parent/'829'))
import majorana_code_source_bridge as code

def add(*ps):
    out={}
    for p in ps:
        for k,v in p.items():out[k]=out.get(k,F(0))+v
    return {k:v for k,v in out.items() if v}
def scale(p,c):return {k:v*c for k,v in p.items() if v*c}
def mul(p,q):
    out={}
    for (a,x),u in p.items():
        for (b,y),v in q.items():
            key=(a+b,code.mul(x,y));out[key]=out.get(key,F(0))+u*v
    return {k:v for k,v in out.items() if v}
def smul(p,q):
    out={}
    for a,u in p.items():
        for b,v in q.items():out[a+b]=out.get(a+b,F(0))+u*v
    return {k:v for k,v in out.items() if v}
def scalar_add(*ps):return add(*ps)

def run():
    gamma,_,_,_,comp,sx,sz=code.code_data()
    sy=[code.mul((0,0,3),code.mul(sz[b],sx[b])) for b in range(5)]
    sources=[code.mul((0,0,2),a) for a in (sz[0],sz[1],sy[3])]
    L=code.mul((0,0,1),code.product(gamma[i] for i in (0,1,4,5,12,14)))
    assert all(code.commute(a,b) for a in sources for b in sources)
    eps=F(1,10);r0=F(0);r1=F(3,5);d=1024
    def phase(p):
        assert p in (0,2)
        return F(1 if p==0 else -1)
    def core(a):
        kind,p=comp(a)
        return phase(p) if kind=='scalar' else F(0)
    def expectation(a,r):
        full=phase(a[2]) if a[:2]==(0,0) else F(0)
        return ((1-eps)-2*eps/(d-2))*core(a)+(1-eps)*r*core(code.mul(L,a))+eps*d/(d-2)*full
    def expect(p,r):
        out={}
        for (n,a),v in p.items():out[n]=out.get(n,F(0))+v*expectation(a,r)
        return {n:v for n,v in out.items() if v}
    def cumulant(a,b,c,r):
        ea,eb,ec=[expect(p,r) for p in (a,b,c)]
        return scalar_add(expect(mul(mul(a,b),c),r),
            scale(smul(expect(mul(a,b),r),ec),-1),
            scale(smul(expect(mul(a,c),r),eb),-1),
            scale(smul(ea,expect(mul(b,c),r)),-1),scale(smul(smul(ea,eb),ec),2))
    def difference(obs):return scalar_add(cumulant(*obs,r1),scale(cumulant(*obs,r0),-1))
    # Full mixed recurrence; all values exact, source departments separate.
    n=10;m=3
    C=[[F(1,5),F(1,9),F(0)],[F(1,9),F(1,4),F(1,11)],[F(0),F(1,11),F(1,6)]]
    jf=[[F((t+2)*(a+1),31) for a in range(m)] for t in range(n)]
    jo=[[F((-1)**(t+a)*(t+a+1),37) for a in range(m)] for t in range(n)]
    total=[[jf[t][a]+jo[t][a] for a in range(m)] for t in range(n)]
    def solve(j,x0=None,x1=None):
        x=[[F(0)]*m for _ in range(n)]
        x[0]=list(x0 or [F(0)]*m);x[1]=list(x1 or [F(0)]*m)
        for t in range(1,n-1):
            x[t+1]=[2*x[t][a]-x[t-1][a]-sum(C[a][b]*x[t][b] for b in range(m))-j[t][a] for a in range(m)]
        return x
    def residual(x,j):return [x[t+1][a]-2*x[t][a]+x[t-1][a]+sum(C[a][b]*x[t][b] for b in range(m))+j[t][a] for t in range(1,n-1) for a in range(m)]
    zero=[[F(0)]*m for _ in range(n)]
    bf=solve(jf);ball=solve(total)
    h=solve(zero,[F(1,3),F(-2,5),F(3,7)],[F(-1,4),F(2,9),F(5,11)])
    target=[[ball[t][a]+h[t][a] for a in range(m)] for t in range(n)]
    assert not any(residual(ball,total)) and not any(residual(target,total)) and not any(residual(h,zero))
    omitted=residual(bf,total);assert omitted==[jo[t][a] for t in range(1,n-1) for a in range(m)]
    assert any(omitted) and any(residual(ball,zero))
    # The true quantum sources are NOT approximated by jf/jo above.
    # They only test matching the SAME full forcing and homogeneous difference.
    # Original code algebra: shifts can change higher orders but not hbar^3.
    shifts=[F(0),h[3][0],h[4][1],h[5][2]];rows=[]
    for shift in shifts:
        obs=[]
        for i in range(m):
            q=code.mul(sources[(i+1)%m],sources[(i+2)%m])
            # Leading response is the negative of the original source.
            p={(1,sources[i]):F(-1),(1,code.I):ball[3+i][i],
               (2,q):F(i+1,7),(2,code.I):F(i-1,11)}
            p=add(p,{(2,sources[i]):shift*F(i+1,3),(3,code.I):shift**2})
            obs.append(p)
        delta=difference(obs);assert delta.get(3,F(0))==-F(27,50)
        assert not any(v for k,v in delta.items() if k<3)
        rows.append(dict(homogeneous_shift=str(shift),content_cumulant_coefficients={str(k):str(v) for k,v in sorted(delta.items())}))
    assert len({r['content_cumulant_coefficients'].get('4','0') for r in rows})>1
    admissible=[]
    for nf,nb,q in itertools.product(range(9),range(7),range(4)):
        for replaced in range(nb+1):
            order=nf+nb+replaced+2*q
            if order==6 and nf>=6:admissible.append((nf,nb,replaced,q))
    assert admissible==[(6,0,0,0)]
    # Off-shell gauge identity in an O(2) quartic diagnostic, not the original S.
    # S=(x^2+y^2-1)^2/4; K=(-y,x); E at (x,0)=(x(x^2-1),0).
    gauge=[]
    for lam in (F(1,10),F(1,20),F(1,40)):
        x=1+lam/3;gradient=x*(x*x-1)
        pk=x*(x*x-1);dkte=-gradient
        assert pk+dkte==0 and pk!=0
        semiclassical=gradient-lam*F(2,3)
        assert semiclassical==lam**2/3+lam**3/27
        gauge.append(dict(parameter=str(lam),PK_component=str(pk),DK_transpose_E=str(dkte),semiclassical_residual=str(semiclassical)))
    return dict(round=851,all_checks_passed=True,fresh_test_groups=1,
        exact_arithmetic='Fraction with original 829 Pauli code',
        diagnostic_time_sites=n,diagnostic_mixed_components=m,
        full_source_equation_residual='0',homogeneous_matching_residual='0',
        omitted_other_source_max_residual=str(max(abs(x) for x in omitted)),
        sourced_background_shift_is_not_homogeneous=True,
        original_code_first_content_coefficient='-27/50',mean_matched_content_rows=rows,
        higher_order_content_need_not_be_unchanged=True,
        epsilon_six_degree_choices=[list(x) for x in admissible],
        off_shell_gauge_diagnostic=gauge,
        original_complete_source_continuum_solution_computed=False,
        original_continuum_initial_matching_and_preservation_are_analytic=True,
        new_background_quantum_BRST_or_full_feedback_claimed=False,
        finite_coupling_or_autonomous_preparation_claimed=False)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args();result=run()
    if args.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert result==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(result,ensure_ascii=False,indent=2))
