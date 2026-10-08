"""Round 1036: finite geometry and exact mass-kernel elimination.

Python standard library + NumPy. Rational identities are exact; numerical
resolvents only calibrate the proved finite-window bounds. Default is read-only.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / 'spectral_mass_transport_results.json'
V = [(i, j) for i in range(1, 5) for j in range(1, 5) if i != j]
G = np.array([[0,-1,-1,1],[1,0,-1,-1],[1,1,0,-1],[-1,1,1,0]], dtype=int)
LIGHT = [V.index((1,3)), V.index((3,1))]
HEAVY = [i for i in range(12) if i not in LIGHT]


def fm(a):
    return [[F(int(x)) for x in row] for row in a]


def transpose(a):
    return [list(row) for row in zip(*a)]


def mm(a, b):
    return [[sum((x*y for x,y in zip(row,col)), F(0))
             for col in zip(*b)] for row in a]


def inverse_det(a):
    n=len(a)
    m=[list(row)+[F(int(i==j)) for j in range(n)] for i,row in enumerate(a)]
    det=F(1)
    for j in range(n):
        p=next((i for i in range(j,n) if m[i][j]), None)
        if p is None:
            return None, F(0)
        if p!=j:
            m[j],m[p]=m[p],m[j];det=-det
        pivot=m[j][j];det*=pivot
        m[j]=[x/pivot for x in m[j]]
        for i in range(n):
            if i!=j:
                c=m[i][j]
                m[i]=[x-c*y for x,y in zip(m[i],m[j])]
    return [row[n:] for row in m],det


def comm(a,b):
    return a@b-b@a


def model():
    gamma=np.diag([G[i-1,j-1] for i,j in V])
    J=np.array([[int(v==w[::-1]) for w in V] for v in V])
    left=[np.diag([int(v[0]==i) for v in V]) for i in range(1,5)]
    right=[J@a@J for a in left]
    D=np.array([[int(gamma[a,a]!=gamma[b,b] and
                        (v[0]==w[0] or v[1]==w[1]))
                  for b,w in enumerate(V)] for a,v in enumerate(V)],dtype=int)
    for v,w in [((1,2),(1,4)),((2,1),(4,1))]:
        a,b=V.index(v),V.index(w);D[a,b]=D[b,a]=2
    return gamma,J,left,right,D


def components(D):
    unseen=set(range(len(D)));answer=[]
    while unseen:
        todo=[min(unseen)];part=set()
        while todo:
            a=todo.pop()
            if a in part:continue
            part.add(a)
            todo.extend(int(x) for x in np.flatnonzero(D[a]) if x not in part)
        answer.append(sorted(part));unseen-=part
    return answer


def run():
    gamma,J,left,right,D=model();eye=np.eye(12,dtype=int)
    identities={'selfadjoint':D-D.T,'J_square':J@J-eye,
                'gamma_square':gamma@gamma-eye,'KO6_J_gamma':J@gamma+gamma@J,
                'J_D':J@D-D@J,'odd_D':gamma@D+D@gamma}
    for i,a in enumerate(left):
        identities[f'gamma_left_{i}']=comm(gamma,a)
        for j,b in enumerate(right):
            identities[f'zero_order_{i}_{j}']=comm(a,b)
            identities[f'first_order_{i}_{j}']=comm(comm(D,a),b)
    orientation=sum(G[i,j]*left[i]@right[j] for i in range(4) for j in range(4))
    identities['orientation']=orientation-gamma
    assert all(not np.any(x) for x in identities.values())
    intersection=np.array([[np.trace(gamma@a@b) for b in right] for a in left])
    assert np.array_equal(intersection,G)
    _,det_intersection=inverse_det(fm(intersection))
    _,det_D=inverse_det(fm(D))
    B=D[np.ix_(HEAVY,HEAVY)];C=D[np.ix_(LIGHT,HEAVY)]
    invB,detB=inverse_det(fm(B))
    zero_kernel=[[-x for x in row] for row in mm(mm(fm(C),invB),transpose(fm(C)))]
    invB_frob_sq=sum(x*x for row in invB for x in row)
    assert det_intersection==1 and det_D==16 and detB==-4
    assert zero_kernel==[[F(0),F(2)],[F(2),F(0)]]
    assert invB_frob_sq==F(41,2)
    assert np.array_equal(C@C.T,2*np.eye(2,dtype=int))
    graph=components(D)
    assert len(graph)==1
    path=[(1,3),(1,4),(3,4),(3,1)]
    assert all(D[V.index(a),V.index(b)] for a,b in zip(path,path[1:]))
    # A global epsilon must be a real diagonal matrix constant on this graph.
    # The J-odd condition then makes its only linear solution zero.
    eps_equations=[]
    for a in range(12):
        for b in range(a+1,12):
            if D[a,b]:
                row=[F(0)]*12;row[a]=F(1);row[b]=F(-1);eps_equations.append(row)
        row=[F(0)]*12;row[a]=F(1);row[V.index(V[a][::-1])]+=F(1)
        eps_equations.append(row)
    gram=mm(transpose(eps_equations),eps_equations)
    _,eps_gram_det=inverse_det(gram)
    assert eps_gram_det!=0
    numerical=[]
    for t,M in [(1.,1.),(.1,5.),(.25,3.)]:
        full=np.zeros((12,12));full[np.ix_(HEAVY,HEAVY)]=M*B
        full[np.ix_(LIGHT,HEAVY)]=t*C;full[np.ix_(HEAVY,LIGHT)]=t*C.T
        for name,res in [('J_D',J@full-full@J),('odd_D',gamma@full+full@gamma)]:
            assert np.linalg.norm(res)<1e-14
        s0=(t*t/M)*np.array([[0.,2.],[2.,0.]])
        for ratio in [-.01,-.005,0.,.005,.01]:
            z=ratio*M
            sigma=-t*t*C@np.linalg.solve(M*B-z*np.eye(10),C.T)
            actual=float(np.linalg.norm(sigma-s0,2))
            bound=2*t*t*abs(z)/((M/5)*(M/5-abs(z)))
            assert actual<=bound+1e-12
            assert abs(sigma[0,1])+1e-12 >= (28/19)*t*t/M
            numerical.append({'t':t,'M':M,'z':z,'difference_norm':actual,
                              'analytic_bound':bound,'offdiagonal':float(sigma[0,1])})
        for ratio in [.005j,.05j,.1j]:
            z=ratio*M
            sigma=-t*t*C@np.linalg.solve(M*B-z*np.eye(10),C.T)
            exact=np.linalg.inv(full-z*np.eye(12))[np.ix_(LIGHT,LIGHT)]
            residual=float(np.linalg.norm(exact-np.linalg.inv(sigma-z*np.eye(2)),2))
            assert residual<1e-9
            numerical.append({'t':t,'M':M,'complex_z':[z.real,z.imag],
                              'compressed_resolvent_residual':residual})
    reduced_left=[a[np.ix_(LIGHT,LIGHT)] for a in left]
    reduced_right=[a[np.ix_(LIGHT,LIGHT)] for a in right]
    s0=np.array([[0,2],[2,0]])
    first_order_defect=comm(comm(s0,reduced_left[0]),reduced_right[0])
    assert np.array_equal(first_order_defect,-s0)
    return {'round':1036,'status':'local analytical and computational result',
            'scientific_group':'finite_geometry_mass_transport',
            'global_cumulative_count':None,'new_cognitive_axioms':0,
            'vertices':V,'gamma':np.diag(gamma).tolist(),'D_base':D.tolist(),
            'exact_identity_count':len(identities),'all_exact_identities_zero':True,
            'intersection_form':intersection.tolist(),'intersection_determinant':str(det_intersection),
            'D_determinant':str(det_D),'heavy_determinant':str(detB),
            'heavy_inverse_frobenius_squared':str(invB_frob_sq),
            'heavy_inverse':[[str(x) for x in row] for row in invB],
            'light_vertices':[V[i] for i in LIGHT],
            'zero_frequency_kernel':[[str(x) for x in row] for row in zero_kernel],
            'components':graph,'conjugate_path':path,'epsilon_gram_determinant':str(eps_gram_det),
            'first_order_defect_after_elimination':first_order_defect.tolist(),
            'heavy_gap_diagnostic':float(min(abs(np.linalg.eigvalsh(B)))),
            'rigorous_heavy_gap_lower_bound':'M/5',
            'window':'|z| <= M/100', 'offdiagonal_lower_bound':'28 t^2/(19 M)',
            'numerical_cases':numerical,
            'limits':['not a physical neutrino model','not a pole-mass identification',
                      'not the full classification with minimality and physical filters',
                      'not the full A1_D-A7 operational contract']}


def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:compare(a[k],b[k])
    elif isinstance(a,(list,tuple)):
        assert len(a)==len(b)
        for x,y in zip(a,b):compare(x,y)
    elif isinstance(a,float):assert abs(a-b)<1e-9*(1+abs(a))
    else:assert a==b,(a,b)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    result=run()
    if args.write:
        with OUT.open('x',encoding='utf8') as f:json.dump(result,f,ensure_ascii=False,indent=2)
    elif OUT.exists():compare(result,json.loads(OUT.read_text('utf8')))
    print(json.dumps({'round':1036,'exact_identities':result['exact_identity_count'],
                      'intersection_determinant':result['intersection_determinant'],
                      'heavy_determinant':result['heavy_determinant'],
                      'schur_kernel':result['zero_frequency_kernel'],
                      'numerical_cases':len(result['numerical_cases']),'passed':True},ensure_ascii=False))
