"""Local independent-algorithm checks; not independent-agent peer review."""
from pathlib import Path
from fractions import Fraction
from urllib.parse import unquote
import argparse
import ast
import hashlib
import json
import re
import numpy as np
import spectral_mass_transport as science

HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
BASE=STAGE.parent
OUT=HERE/'research_round_1036_checks.json'
NOTE=STAGE/'research_note_1036.md'
OWN=['spectral_mass_transport.py','spectral_mass_transport_results.json',
     'drafts/spectral_mass_transport_derivation.md','selection_audit.md',
     'input_dependency_update_v0_25.md','review.md','NEXT.md','verify_round1036.py']
INPUTS=[BASE/'archive_531_553/research_note_531.md',
        BASE/'archive_956_989/957/drafts/unified_operation_hypotheses_v0_2.md',
        STAGE/'1009/input_dependency_ledger_v0_1.md',
        STAGE/'research_note_1012.md',STAGE/'research_note_1014.md',
        STAGE/'research_note_1025.md',STAGE/'1025/NEXT.md',
        STAGE/'research_note_1035.md',STAGE/'1035/drafts/joint_selection_derivation.md',
        STAGE/'1035/input_dependency_update_v0_24.md']


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def independent_matrix_checks():
    # Explicit edge table, not the construction rule in the science script.
    vertices=['12','13','14','21','23','24','31','32','34','41','42','43']
    edges=[('12','14',2),('12','32',1),('12','42',1),('13','14',1),
           ('13','43',1),('14','24',1),('14','34',1),('21','23',1),
           ('21','24',1),('21','41',2),('23','43',1),('31','34',1),
           ('31','41',1),('32','34',1),('41','42',1),('41','43',1)]
    D=np.zeros((12,12),complex)
    for v,w,c in edges:
        a,b=vertices.index(v),vertices.index(w);D[a,b]=D[b,a]=c
    J=np.zeros((12,12))
    for a,v in enumerate(vertices):J[a,vertices.index(v[::-1])]=1
    gamma=np.diag([-1,-1,1,1,-1,-1,1,1,-1,-1,1,1])
    A=[np.diag([int(v[0]==str(i)) for v in vertices]) for i in range(1,5)]
    saved=json.loads(science.OUT.read_text('utf8'))
    assert np.array_equal(D,np.array(saved['D_base']))
    assert np.linalg.norm(J@D.conj()-D@J)==0
    assert np.linalg.norm(gamma@D+D@gamma)==0
    eps_basis=[]
    for i in range(12):
        e=np.zeros((12,12),complex);e[i,i]=1;eps_basis.append(e)
    for i in range(12):
        for j in range(i+1,12):
            e=np.zeros((12,12),complex);e[i,j]=e[j,i]=1;eps_basis.append(e)
            e=np.zeros((12,12),complex);e[i,j]=1j;e[j,i]=-1j;eps_basis.append(e)
    def linear_residual(e):
        terms=[e@a-a@e for a in A]
        terms += [e@D-D@e,e@J+J@e.conj()]
        val=np.concatenate([r.ravel() for r in terms])
        return np.concatenate([val.real,val.imag])
    constraint=np.column_stack([linear_residual(e) for e in eps_basis])
    singular=np.linalg.svd(constraint,compute_uv=False)
    assert len(eps_basis)==144 and min(singular)>0.1
    # Direct two-column solve, independently verifies the printed CB^-1.
    l=[1,6];h=[i for i in range(12) if i not in l]
    B=D[np.ix_(h,h)];C=D[np.ix_(l,h)]
    printed=np.array([[1,0,0,2,0,0,-1,-1,0,0],
                      [0,-1,1,0,0,2,0,0,0,-1]])
    assert np.array_equal(printed@B,C)
    assert np.array_equal(-printed@C.T,np.array([[0,2],[2,0]]))
    invB=np.array([[float(Fraction(x)) for x in row] for row in saved['heavy_inverse']])
    inv_res=float(np.linalg.norm(B@invB-np.eye(10),2))
    assert inv_res==0
    defect=(np.array([[0,2],[2,0]])@np.diag([1,0])-np.diag([1,0])@np.array([[0,2],[2,0]]))
    defect=defect@np.diag([0,1])-np.diag([0,1])@defect
    assert np.array_equal(defect,-np.array([[0,2],[2,0]]))
    # Deleting the weighted J-edge pair leaves every slot incident: nonminimal.
    trimmed=D.copy()
    for a,b in [(0,2),(3,9)]:trimmed[a,b]=trimmed[b,a]=0
    assert np.all(np.sum(abs(trimmed),axis=1)>0)
    return {'explicit_edges':len(edges),'full_hermitian_unknowns':144,
            'epsilon_constraint_rank':int(np.sum(singular>1e-9)),
            'epsilon_constraint_min_singular':float(min(singular)),
            'heavy_inverse_residual':inv_res,'printed_CB_inverse_exact':True,
            'nonminimal_graph_control':True}


def link_check(path):
    text=path.read_text('utf8')
    text=re.sub(r'```.*?```|\$\$.*?\$\$','',text,flags=re.S)
    count=0
    for target in re.findall(r'(?<!!)\[[^\]\n]*\]\(([^)\n]+)\)',text):
        target=unquote(target.split('#',1)[0])
        if not target or re.match(r'^[A-Za-z]+:',target):continue
        p=(path.parent/target).resolve()
        assert p.exists(),(str(path),target)
        count+=1
    return count


def verify():
    fresh=science.run();science.compare(fresh,json.loads(science.OUT.read_text('utf8')))
    second=independent_matrix_checks()
    for p in [HERE/'spectral_mass_transport.py',HERE/'verify_round1036.py']:
        ast.parse(p.read_text('utf8'))
    docs=[NOTE]+[HERE/p for p in OWN if p.endswith('.md')]
    links=sum(link_check(p) for p in docs)
    for p in docs:
        s=p.read_text('utf8')
        assert '\ufffd' not in s
        assert not any(ord(c)<32 and c not in '\n\r\t' for c in s)
    artifacts=[NOTE]+[HERE/p for p in OWN]
    return {'round':1036,'passed':True,'independent_agent_review':'pending main-line assignment',
            'author_second_algorithm_is_not_independent_review':True,
            'scientific_groups_proposed':1,'global_cumulative_count':None,
            'exact_identity_count':fresh['exact_identity_count'],
            'numerical_cases':len(fresh['numerical_cases']),
            'independent_algorithm':second,'local_links_checked':links,
            'artifact_sha256':{str(p.relative_to(STAGE)):sha(p) for p in artifacts},
            'input_sha256':{str(p.relative_to(BASE)):sha(p) for p in INPUTS}}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    # Let local links refer to this receipt during its first creation.
    if args.write:
        with OUT.open('x',encoding='utf8') as f:f.write('{}\n')
    result=verify()
    if args.write:OUT.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n','utf8')
    elif OUT.exists():science.compare(result,json.loads(OUT.read_text('utf8')))
    print(json.dumps({'round':1036,'passed':True,'links':result['local_links_checked'],
                      'hermitian_epsilon_rank':result['independent_algorithm']['epsilon_constraint_rank'],
                      'artifacts':len(result['artifact_sha256'])},ensure_ascii=False))
