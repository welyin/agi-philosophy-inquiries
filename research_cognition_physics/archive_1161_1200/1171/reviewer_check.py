"""1171 independent review, without importing the author's implementation.

Physical H is applied to sparse integer Fourier labels with no momentum grid.
Its invariant 20-column sector is assembled from explicit prefix formulae,
then cross-checked against two independent five-dimensional head propagators.
Default is read-only; --write-once creates review_checks.json exclusively.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math
import re
import subprocess
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
EXPECTED = {
 'archive_1161_/research_note_1171.md': '47840ba4fcd53a456ddd68584070b049eca2c84e6d404d00ee4447222b66f16a',
 'archive_1161_/1171/admission.md': '4b1b025cd58c81aef92c2ca8659477a8716a8a6a7c18bc6b18020be3f8790ad4',
 'archive_1161_/1171/proof.md': 'ecf48b81b9b42a21bcf7007b386bae16132cb9849e3fc939f02c398860e8d4ea',
 'archive_1161_/1171/check.py': '7a88684bad79f61f7f436bcfb678f0d6a5f28620074a3d1b512463ff4a2e0892',
 'archive_1161_/1171/results.json': 'fe43f6ec9c43352ea989605538fb90b2c7c4325ce2c148a3be09020a6d100d3c',
 'archive_1161_/1171/NEXT.md': 'ee88f0bd493ee621fa9d7e8a9335b02a8970a91663d8c8671c936f39b8a52794',
}
ID = np.eye(2, dtype=complex)
Y = np.array([[0,-1j],[1j,0]],complex)
X = np.array([[0,1],[1,0]],complex)
CY = np.kron((ID+Y)/2,ID)+np.kron((ID-Y)/2,X)
INERTIA = 1000.
TIME = math.pi


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def opnorm(a):
    return float(np.linalg.norm(a,2))


def must_zero(a, tol=3e-12):
    error=opnorm(a)
    assert error<tol,error
    return error


def add(out,key,value):
    out[key]=out.get(key,0j)+value


def gate_on_label(p,a,l,edge,inverse=False,n=0):
    """Actual edge link on one Fourier basis vector, on the infinite rotor."""
    if edge==0 or edge==2:
        sign=1 if edge==0 else -1
        if inverse: sign=-sign
        return {(p+sign*(n-a),a,l):1.+0j}
    if edge==1:
        sign=1 if inverse else -1
        return {(p,a,l):np.exp(sign*.5j*math.pi*p*p)}
    assert edge==3
    return {(p,b,r):CY[2*b+r,2*a+l] for b in range(2) for r in range(2)
            if CY[2*b+r,2*a+l]!=0}


def physical_h(state,n):
    out={}
    for (j,p,a,l),value in state.items():
        add(out,(j,p,a,l),(2+p*p/(2*INERTIA))*value)
        for direction in [-1,1]:
            dest=j+direction
            if not 0<=dest<5: continue
            edge=min(j,dest)
            hopping=-math.sqrt((edge+1)*(4-edge))/2
            for (pp,aa,ll),phase in gate_on_label(p,a,l,edge,direction<0,n).items():
                add(out,(dest,pp,aa,ll),hopping*phase*value)
    return out


def sector_columns(n):
    """Closed prefix formulae, not a repeated application of author's gates."""
    columns=[]
    for j in range(5):
        for a in range(2):
            for l in range(2):
                b=a-n
                p=-b if j in [1,2] else 0
                phase=np.exp(-.5j*math.pi*b*b) if j>=2 else 1.+0j
                if j<4: col={(j,p,a,l):phase}
                else:
                    col={(j,0,aa,ll):phase*CY[2*aa+ll,2*a+l]
                         for aa in range(2) for ll in range(2)
                         if CY[2*aa+ll,2*a+l]!=0}
                columns.append(col)
    return columns


def matrix_on_rows(columns,rows):
    return np.array([[c.get(row,0j) for c in columns] for row in rows],complex)


def exp_h(h,t):
    values,vectors=np.linalg.eigh(h)
    return (vectors*np.exp(-1j*t*values))@vectors.conj().T


def full_physical_check(saved):
    xhead=np.zeros((5,5),complex)
    for j in range(4):
        xhead[j,j+1]=xhead[j+1,j]=math.sqrt((j+1)*(4-j))/2
    base=2*np.eye(5)-xhead
    start5=np.eye(5)[:,0]
    target5=np.eye(5)[:,4]
    vec0=exp_h(base,TIME)@start5
    vec1=exp_h(base+np.diag([0,1,1,0,0])/(2*INERTIA),TIME)@start5
    z=vec1[4]
    scalar_probabilities=[float(abs(1+z)**2/4),float(abs(1-z)**2/4)]
    scalar_isometry=float(np.linalg.norm(vec1-target5))
    eps=math.pi/(2*INERTIA)
    assert scalar_isometry<eps
    cases=[]
    for n in [0,1]:
        columns=sector_columns(n)
        images=[physical_h(c,n) for c in columns]
        # Include every label produced by actual H; never discard outside labels.
        rows=sorted(set().union(*(set(c) for c in columns+images)))
        embed=matrix_on_rows(columns,rows)
        hembed=matrix_on_rows(images,rows)
        iso=must_zero(embed.conj().T@embed-np.eye(20))
        h=embed.conj().T@hembed
        hermitian=must_zero(h-h.conj().T)
        leakage=must_zero(hembed-embed@h)
        q=np.zeros(20)
        for j in [1,2]:
            for a in range(2): q[4*j+2*a:4*j+2*a+2]=(a-n)**2
        target=np.kron(base,np.eye(4))+np.diag(q)/(2*INERTIA)
        generator_error=must_zero(h-target)
        ein=np.zeros((20,2),complex); ein[0,0]=ein[2,1]=1
        ideal=np.zeros_like(ein); ideal[16,0]=ideal[18,1]=1
        evolution=exp_h(h,TIME)
        actual=evolution@ein
        operator_error=opnorm(actual-ideal)
        assert abs(operator_error-scalar_isometry)<3e-12
        output=embed@actual
        plus=np.ones(2)/math.sqrt(2)
        vec=output@plus
        record1=np.array([r[3]==1 for r in rows])
        probability=float(np.vdot(vec,record1*vec).real)
        probability_error=abs(probability-scalar_probabilities[n])
        assert probability_error<3e-12
        recorded=saved['autonomous_positive_completion']['cases'][n]
        assert abs(probability-recorded['record1_probability_all_heads'])<3e-12
        assert abs(operator_error-recorded['prepared_isometry_error'])<3e-12
        mean=ein.conj().T@h@ein
        variance=ein.conj().T@h@h@ein-mean@mean
        variance_error=must_zero(variance-np.eye(2))
        mean_error=must_zero(mean-2*np.eye(2))
        final_mean_error=must_zero(actual.conj().T@h@actual-mean)
        final_h2_error=must_zero(actual.conj().T@h@h@actual-ein.conj().T@h@h@ein)
        # Compare complete, untraced output Gram maps, uniform on all input/reference.
        full_isometry_error=must_zero(output.conj().T@output-np.eye(2))
        actual_kinetic=output.conj().T@(np.array([r[1]**2/(2*INERTIA) for r in rows])[:,None]*output)
        cases.append({'n':n,'exact_sector_dimension':20,'physical_labels_in_h_action':len(rows),
            'physical_modes_in_h_action':sorted(set(r[1] for r in rows)),
            'embedding_error':iso,'physical_H_sector_leakage':leakage,
            'physical_H_hermiticity_error':hermitian,'compressed_H_formula_error':generator_error,
            'prepared_full_isometry_error':operator_error,
            'record1_probability_no_postselection':probability,
            'independent_5x5_record_probability_error':probability_error,
            'complete_output_isometry_error':full_isometry_error,
            'initial_energy_mean':2.,'initial_energy_variance':1.,
            'energy_mean_error':mean_error,'energy_variance_error':variance_error,
            'complete_mean_conservation_error':final_mean_error,
            'complete_second_moment_conservation_error':final_h2_error,
            'final_rotor_energy_norm':opnorm(actual_kinetic)})
    gap=cases[0]['record1_probability_no_postselection']-cases[1]['record1_probability_no_postselection']
    assert gap>=1-math.pi/INERTIA
    return {'status':'PASS','cases':cases,'independent_5x5_PST_error':float(np.linalg.norm(vec0-target5)),
        'independent_5x5_transfer':[float(z.real),float(z.imag)],
        'uniform_prepared_error_bound':eps,'actual_gap':gap,'guaranteed_gap':1-math.pi/INERTIA,
        'rotor_truncation':False,'postselection':False,
        'general_reference_reason':'Operator norm of the whole prepared isometry tensors unchanged with I_R.'}


def source_check():
    a=np.diag([0.,1.]); b=(ID+X)/2
    q=.37
    # Exact projector exponential and analytic derivatives, no author routine.
    u=ID+(np.exp(-1j*q)-1)*a
    v=ID+(np.exp(-1j*q)-1)*b
    du=-1j*np.exp(-1j*q)*a; dv=-1j*np.exp(-1j*q)*b
    direct=1j*(u@v).conj().T@(du@v+u@dv)
    composed=v.conj().T@a@v+b
    composition=must_zero(direct-composed)
    omission=opnorm(direct-(a+b))
    assert omission>.1
    matrix_super=-1j*(np.kron(ID,composed)-np.kron(composed.T,ID))
    shifted=composed+7.3*ID
    shifted_super=-1j*(np.kron(ID,shifted)-np.kron(shifted.T,ID))
    scalar_invisible=must_zero(matrix_super-shifted_super)
    vals,vs=np.linalg.eigh(composed)
    psi=(vs[:,0]+vs[:,-1])/math.sqrt(2)
    rho=np.outer(psi,psi.conj())
    tangent=-1j*(composed@rho-rho@composed)
    trace_norm=float(np.linalg.svd(tangent,compute_uv=False).sum())
    diameter=float(vals[-1]-vals[0])
    assert abs(trace_norm-diameter)<3e-12
    bounds=[]
    for n in [16,64,256]:
        theta=1/math.sqrt(n)
        u=ID+(np.exp(-1j*theta)-1)*a
        plus=np.ones(2)/math.sqrt(2)
        rho=np.outer(plus,plus)
        difference=u@rho@u.conj().T-rho
        attained=float(np.linalg.svd(difference,compute_uv=False).sum())
        exact=2*math.sin(theta/2)
        assert abs(attained-exact)<3e-12
        bounds.append({'N':n,'finite_channel_difference':attained,'source_diameter_at_zero':math.sqrt(n)})
    return {'status':'PASS','composition_error':composition,'error_if_pullback_omitted':omission,
        'scalar_source_superoperator_error':scalar_invisible,'source_diameter':diameter,
        'extremal_state_tangent_trace_norm':trace_norm,'C0_counterexample':bounds,
        'diamond_norm_scope':'Self-contained complete upper bound plus attained pure-state lower bound; no sampled SDP claim.'}


def text_checks():
    local_links=0; docs=[]
    for rel in EXPECTED:
        path=BASE/rel
        if path.suffix!='.md': continue
        text=path.read_text(encoding='utf-8-sig')
        assert all(ord(c)>=32 or c in '\n\r\t' for c in text),rel
        assert text.count('$$')%2==0,rel
        for raw in re.findall(r'\]\(([^)]+)\)',text):
            target=raw.strip().strip('<>').split('#')[0]
            if not target or re.match(r'\w+://',target): continue
            assert (path.parent/target).resolve().exists(),(rel,target)
            local_links+=1
        assert '\\!left' not in text,rel
        docs.append(rel)
    return {'status':'PASS','markdown_files':docs,'valid_local_link_occurrences':local_links,
            'display_math_delimiters_paired':True,'control_characters_absent':True,'image_inspection':False}


def calculate():
    before={rel:sha(BASE/rel) for rel in EXPECTED}
    assert before==EXPECTED,before
    saved=json.loads((HERE/'results.json').read_text(encoding='utf-8-sig'))
    sources=saved['source_sha256']
    assert len(sources)==9
    for rel,want in sources.items(): assert sha(BASE/rel)==want,rel
    completed=subprocess.run([sys.executable,'-B','-X','utf8',str(HERE/'check.py')],
        capture_output=True,text=True,encoding='utf-8',check=False)
    assert completed.returncode==0,(completed.stdout,completed.stderr)
    default=json.loads(completed.stdout.strip())
    assert default['mode']=='READ_ONLY_MATCH' and default['groups']==3
    physical=full_physical_check(saved)
    source=source_check()
    checks=text_checks()
    assert {rel:sha(BASE/rel) for rel in EXPECTED}==before
    for rel,want in sources.items(): assert sha(BASE/rel)==want,rel
    return {'schema':'round1171_independent_review_v1','round':1171,'date':'2026-10-13',
        'verdict':'PASS_WITH_SCOPE','blockers':0,
        'reviewer':'/root/local_dynamics_goal_scope_audit','blind_review':False,
        'disclosure':'Pre-admission source/phase scope and loop review; suggested energy-expectation wording. No author core edits.',
        'core_sha256':EXPECTED,'source_sha256':sources,'source_count':len(sources),
        'author_default':{'exit_code':completed.returncode,'stdout':default,'stderr':completed.stderr},
        'independent_physical_Fourier_and_5x5':physical,'independent_source_check':source,
        'text_and_links':checks,'frozen_inputs_unchanged':True,
        'scope':{'new_adopted_axioms':0,'empirical_calibrations':0,'GR_hard_deliverables':0,
          'hardware_accuracy_certified':False,'unknown_black_box_control_derived':False,
          'complete_channel_equality_domain':'Original D at fixed q and known n only.',
          'coherent_implementation_permission_is_input':True},
        'reviewer_script_sha256':sha(Path(__file__))}


def compare(a,b,path='root'):
    if isinstance(a,dict):
        assert isinstance(b,dict) and a.keys()==b.keys(),path
        for key in a: compare(a[key],b[key],path+'.'+key)
    elif isinstance(a,list):
        assert isinstance(b,list) and len(a)==len(b),path
        for i,(x,y) in enumerate(zip(a,b)): compare(x,y,path+'.'+str(i))
    elif isinstance(a,float):
        assert math.isclose(a,b,rel_tol=3e-10,abs_tol=3e-11),(path,a,b)
    else: assert a==b,(path,a,b)


if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('--write-once',action='store_true')
    args=parser.parse_args(); result=calculate(); path=HERE/'review_checks.json'
    if args.write_once:
        with path.open('x',encoding='utf-8',newline='\n') as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
        mode='WRITE_ONCE'
    else:
        compare(result,json.loads(path.read_text(encoding='utf-8-sig')))
        mode='READ_ONLY_MATCH'
    print(json.dumps({'verdict':result['verdict'],'mode':mode,'blockers':0,
          'actual_gap':result['independent_physical_Fourier_and_5x5']['actual_gap'],
          'source_count':result['source_count']},ensure_ascii=False))
