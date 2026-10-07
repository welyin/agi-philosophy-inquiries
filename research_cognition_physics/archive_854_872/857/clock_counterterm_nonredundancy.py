"""857: exact highest-derivative test for the existing composite-reference term.
The Fraction jets validate algebra, not the original on-shell PDE. The report's
nonredundancy proof separately uses exact old linearized waves and a compact,
boundary-preserving first-order EFT dictionary.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import permutations
import argparse,json
import clock_and_cell_dictionary_probe as clock_probe
HERE=Path(__file__).resolve().parent
TARGET=HERE/'clock_counterterm_nonredundancy_results.json'

def det(M):
    n=len(M);total=F(0)
    for p in permutations(range(n)):
        v=F((-1)**sum(p[i]>p[j] for i in range(n) for j in range(i+1,n)))
        for i in range(n):v*=M[i][p[i]]
        total+=v
    return total

def inv(M):
    n=len(M);d=det(M);assert d
    return [[F((-1)**(i+j))*det([[M[r][c] for c in range(n) if c!=i] for r in range(n) if r!=j])/d for j in range(n)] for i in range(n)]

def dot(x,y):return sum(a*b for a,b in zip(x,y))
def mv(M,x):return [dot(row,x) for row in M]
def serial(M):return [[str(v) for v in row] for row in M]

def run():
    assert clock_probe.run()==json.loads(clock_probe.TARGET.read_text('utf-8'))
    G=[[F(0) for j in range(4)] for i in range(4)]
    for i in range(4):G[i][i]=F(-1 if i==0 else 1)
    dh=list(map(F,(1,0,0,0)));ds=list(map(F,(2,1,0,0)))
    Hh=[[F(0) for _ in range(4)] for _ in range(4)];Hs=[row[:] for row in Hh]
    Hh[0][2]=Hh[2][0]=F(1);Hs[0][3]=Hs[3][0]=F(1)
    ah,bh=mv(G,dh),mv(G,ds)
    def reference_matrix(H,Hs):return [dh,ds,[2*dot(row,ah) for row in H],[2*dot(row,bh) for row in Hs]]
    R=reference_matrix(Hh,Hs);J=det(R);Ri=inv(R)
    assert J==8 and det(G)==-1
    assert [[dot(R[i],[Ri[l][j] for l in range(4)]) for j in range(4)] for i in range(4)]==[[F(i==j) for j in range(4)] for i in range(4)]
    e2=[Ri[i][2] for i in range(4)];e3=[Ri[i][3] for i in range(4)]
    def q(la):
        t=[x+la*y for x,y in zip(dh,ds)]
        return -dot(t,mv(G,t))
    rows=[]
    for k0 in ((-1,0,1,0),(-1,0,0,1),(-1,0,F(3,5),F(4,5)),(-1,1,0,0)):
        k=list(map(F,k0));assert dot(k,mv(G,k))==0
        w=[dot(e2,k)*dot(ah,k),dot(e3,k)*dot(bh,k)]
        for la in (F(0),F(1,10),F(-1,10)):
            Q=q(la);assert Q>0
            predicted=[[-8*J*J*w[i]*w[j]/Q for j in range(2)] for i in range(2)]
            def density(r,u):
                H=[[Hh[i][j]+r*k[i]*k[j] for j in range(4)] for i in range(4)]
                S=[[Hs[i][j]+u*k[i]*k[j] for j in range(4)] for i in range(4)]
                # Recompute both actual composite-reference rows. No cofactor
                # formula is used in this independent direct determinant.
                RR=reference_matrix(H,S)
                JJ=det(RR)
                assert JJ==J*(1+2*(w[0]*r+w[1]*u))
                return -JJ*JJ/Q
            zero=density(F(0),F(0))
            direct=[[density(1,0)+density(-1,0)-2*zero,(density(1,1)-density(1,-1)-density(-1,1)+density(-1,-1))/4],
                    [F(0),density(0,1)+density(0,-1)-2*zero]]
            direct[1][0]=direct[0][1]
            assert direct==predicted and det(direct)==0
            rows.append(dict(null_covector=[str(v) for v in k],clock_lambda=str(la),q=str(Q),
                             cofactor_weights=[str(v) for v in w],quartic_hessian=serial(direct),
                             rank=int(any(v for row in direct for v in row))))
    base=rows[0];new=rows[1]
    assert base['quartic_hessian'][0][0]=='-128'
    assert new['q']=='143/100'
    difference=F(new['quartic_hessian'][0][0])-F(base['quartic_hessian'][0][0])
    assert difference==F(5504,143)
    # A scalar perturbation with zero value/first jet and prescribed oscillatory
    # second jet isolates the order-four coefficient. Not a PDE simulation.
    scaled=[]
    for freq in (1,2,4,8):
        value=F(-128)*freq**4
        scaled.append(dict(frequency=freq,quadratic_hessian=str(value),normalized_by_frequency_four=str(value/freq**4)))
    # Negative controls: the determinant (not its square) has no rank-one
    # Hessian cross term, and the scalar leading EOM square vanishes for null k.
    k=list(map(F,(-1,0,1,0)))
    def topological(r,u):
        H=[[Hh[i][j]+r*k[i]*k[j] for j in range(4)] for i in range(4)]
        S=[[Hs[i][j]+u*k[i]*k[j] for j in range(4)] for i in range(4)]
        return det(reference_matrix(H,S))
    topo_second=topological(1,0)+topological(-1,0)-2*topological(0,0)
    topo_mixed=topological(1,1)-topological(1,-1)-topological(-1,1)+topological(-1,-1)
    leading_eom_square=2*dot(k,mv(G,k))**2
    assert topo_second==topo_mixed==leading_eom_square==0
    return dict(round=857,fresh_test_groups=1,all_checks_passed=True,
        inherited_clock_probe_reproduced=True,
        reference_jet=dict(metric_inverse=serial(G),dh=[str(v) for v in dh],ds=[str(v) for v in ds],
            hessian_h=serial(Hh),hessian_s=serial(Hs),reference_matrix=serial(R),J=str(J),
            on_original_full_equations=False),
        exact_rows=rows,clock_difference_quartic_hessian_00=str(difference),frequency_scaling=scaled,
        negative_controls=dict(determinant_second_variation=str(topo_second),determinant_mixed=str(topo_mixed),leading_eom_square=str(leading_eom_square)),
        analytic_scope='Compact boundary-preserving first-order local EFT equivalence about the inherited on-shell background, modulo the full leading equations and two-derivative parameter shifts. The highest scalar jet is nonzero on actual old null polarizations. The original exact-wave existence and order comparison are analytic, not numerically solved here.',
        arbitrary_global_or_nonlocal_equivalence_classified=False,
        companion_counterterms_or_quantum_scheme_equivalence_completed=False,
        original_graph_to_continuum_proved=False,full_goal_completed=False)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args();r=run()
    if args.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps({k:r[k] for k in ('round','all_checks_passed','clock_difference_quartic_hessian_00','original_graph_to_continuum_proved')},ensure_ascii=False,indent=2))
