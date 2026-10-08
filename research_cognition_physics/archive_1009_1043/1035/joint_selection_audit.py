"""1035: typed bridge audit and exact joint-condition comparison family.

Only the stated local anomalies, pure Spin-Z4 obstruction, representation
completion and mass identities are certified. No full global QFT is constructed.
"""
from pathlib import Path
from fractions import Fraction as F
from math import comb
import argparse
import copy
import hashlib
import json
import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
OUT = HERE/'joint_selection_audit_results.json'
LEDGER = HERE/'bridge_ledger_v0_2.json'
OLD = 'archive_1009_/1024/bridge_ledger.json'
HISTORY = [OLD, 'archive_531_553/research_note_531.md',
           'archive_956_989/981/drafts/common_parent_contract_v1.md',
           'archive_990_1008/993/common_candidate_v1.md']
HISTORY += [f'archive_1009_/research_note_{n}.md' for n in range(1010, 1015)]
HISTORY += [f'archive_1009_/research_note_{n}.md' for n in range(1025, 1035)]
HISTORY += ['archive_1009_/1032/drafts/electric_completion_derivation.md',
            'archive_1009_/1034/drafts/discrete_spin_selection_derivation.md']
# d, integer q0=6Y_SM, b=3(B-L), twice color and weak indices, triality, 2j
FIELDS = {'q':(6,1,1,2,3,1,1), 'uc':(3,-4,-1,1,0,2,0),
          'dc':(3,2,-1,1,0,2,0), 'l':(2,-3,-3,0,1,0,1),
          'ec':(1,6,3,0,0,0,0), 'N':(1,0,3,0,0,0,0)}


def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()


def polynomial_anomalies():
    rows = {key:[F(0)]*4 for key in ('gravity_Y','Y3','SU3_2Y','SU2_2Y')}
    for d, q0, b, t3, t2, _, _ in FIELDS.values():
        a, slope = F(q0,6), F(b)  # Y = Y_SM + 3 k (B-L)
        for j in range(4):
            rows['Y3'][j] += d*comb(3,j)*a**(3-j)*slope**j
        for key, weight in (('gravity_Y',d), ('SU3_2Y',t3), ('SU2_2Y',t2)):
            rows[key][0] += weight*a
            rows[key][1] += weight*slope
    assert all(all(v==0 for v in r) for r in rows.values())
    return {key:[str(v) for v in r] for key,r in rows.items()}


def joint_family():
    cases=[]
    for k in range(-6,7):
        charges = {name:q0+6*k*b for name,(_,q0,b,*_) in FIELDS.items()}
        rows=[]; s1=s3=0
        for name,(d,q0,b,t3,t2,triality,n) in FIELDS.items():
            q=charges[name]; x=F(5*b-2*q,3)
            assert x.denominator==1 and x%4==1
            assert (2*triality+3*n+q)%6==0
            s1+=3*d*x; s3+=3*d*x**3
            rows.append(dict(field=name,q_6Y=q,X=int(x),Z4=int(x%4)))
        assert F(11*s3-5*s1,96)%1==0 and F(2*s3+s1,4)%1==0
        yukawa=[charges['q']+3+charges['uc'],charges['q']-3+charges['dc'],
                charges['l']-3+charges['ec'],charges['l']+3+charges['N']]
        assert yukawa==[0]*4
        assert -charges['uc']+charges['dc']==6
        assert 2*charges['N']==36*k  # real neutral S cannot compensate this
        # All central characters, at the six possible H-trivial U(1) phases.
        kernel=[]
        for c in range(3):
            for w in range(2):
                for u in range(6):
                    if (F(w,2)+F(3*u,6))%1: continue
                    if all((F(triality*c,3)+F(n*w,2)+F(charges[name]*u,6))%1==0
                           for name,(_,_,_,_,_,triality,n) in FIELDS.items()):
                        kernel.append([c,w,u])
        assert kernel==sorted([[m%3,m%2,m] for m in range(6)])
        cases.append(dict(k=k,charges=rows,central_kernel=kernel,
                          pure_twisted_obstruction=0,Dirac_Yukawa_charges=yukawa,
                          real_S_Majorana_charge_6Y=36*k,
                          neutrino_electric_charge=-3*k))
    return cases


def representation_words():
    count=0; max_length=0
    for k in range(-3,4):
        for a in range(4):
            for b in range(4):
                for n in range(4):
                    q0=(-2+6*k)*a+(2-6*k)*b+3*n
                    for q in range(-18,19):
                        if (2*(a+2*b)+3*n+q)%6: continue
                        copies=F(q-q0,6)
                        assert copies.denominator==1
                        assert q0+6*copies==q
                        length=a+b+n+2*abs(int(copies))
                        count+=1; max_length=max(max_length,length)
    return dict(cases=count,maximum_tested_word_length=max_length,
                all_irrep_quantifier_proved_analytically=True,
                completion_only_for_internal_G6=True,physical_endpoint_preparation=False)


def mass_blocks():
    cases=[]
    masses=np.array([.1,.2,.3]); charged=np.array([1.,2.,3.])
    mat=np.zeros((12,12))
    for j in range(3):
        mat[j,9+j]=mat[9+j,j]=masses[j]
        mat[3+j,6+j]=mat[6+j,3+j]=charged[j]
    singular=np.linalg.svd(mat,compute_uv=False)
    expected=np.sort(np.repeat(np.r_[masses,charged],2))[::-1]
    assert np.max(np.abs(singular-expected))<1e-13
    for k in (-3,-1,0,1,2,4):
        Q=np.diag([-3*k]*3+[-1-3*k]*3+[1+3*k]*3+[3*k]*3)
        ward=float(np.max(np.abs(Q.T@mat+mat@Q)))
        assert ward<1e-13
        with_majorana=mat.copy()
        with_majorana[9:,9:]=np.diag([10.,13.,17.])
        defect=float(np.max(np.abs(Q.T@with_majorana+with_majorana@Q)))
        assert abs(defect-102*abs(k))<1e-12
        # r maps H to -H and every Weyl field to i times itself.
        U=1j*np.eye(12)
        assert np.max(np.abs(U.T@(-mat)@U-mat))<1e-13
        cases.append(dict(k=k,Dirac_Ward_residual=ward,Majorana_Ward_defect=defect,
                          singular_values=[float(x) for x in singular],
                          positive_neutrino_masses_are_doubled=True))
    return cases


# These classifications are explicit main-agent judgments, NOT inferred from
# file hashes or case counts. The derivation explains the semantic obligations.
NEW_ROWS = [
 (1025,'finite_internal_geometry','conditional_tool','有限母代数比较类','四元数实形式及最小性仍独立；非SM有限见证不直接是P981'),
 (1026,'speed_attraction','conditional_tool','两物种一圈速度流','到全父物种的RG、阈值及误差未映射；不提供1027的精确共同主部'),
 (1027,'charged_source','mapped_interface','P981/B993领先带荷源','明确源类内同源；新增N/S后的味连接不能自动继承'),
 (1028,'gravity_parent_bridge','mapped_interface','同父径向Higgs解及一阶作用','二阶必要条件已接；自由自旋2和共同主部仍输入'),
 (1029,'lift_equivalence','mapped_interface','固定S0/S1的正规局部类','提升模冗余已依赖于作用；S1本身未生成'),
 (1030,'finite_scattering','conditional_tool','稳定归一通道和HEFT切片','实际W/Z/h到压缩S的桥未认证'),
 (1031,'parent_coupling_freedom','mapped_interface','同父领先经典作用解族','有限规范不变量不同；未构造全量子共同过程族'),
 (1032,'electric_completion','conditional_tool','内部群Gp及明确端点表示菜单','E_alg不是实际制备；菜单开放改变商的选择'),
 (1033,'classical_channel','conditional_tool','二态Markov测量反馈族','真实双端引力及同源噪声匹配未认证'),
 (1034,'twisted_spin_mass','conditional_tool','额外Spin-Z4背景及新增质量菜单','原P981固定C5不相容；N/S属于扩展而非原父对象'),
]


def ledger():
    old=json.loads((BASE/OLD).read_text('utf8'))
    rows=copy.deepcopy(old['rows'])
    revisions={
        'scalar_source_to_parent': ('常权领先完整物质源', '1027已补完整规范—Higgs—Weyl来源桥；非零门户纳入B993中性标量；范围外源类仍未分类'),
        'gravity_ideal_to_parent': ('固定一阶作用的共同父对象接口', '1028合法径向见证及1029提升等价补必要条件桥；共同主部、自由自旋2与S1仍采用'),
    }
    for row in rows:
        if row['id'] in revisions:
            row['previous_status']=row['status']
            row['status']='mapped_interface'
            row['object'],row['claim']=revisions[row['id']]
            row['requires']='1035解析稿明确保留1027—1029的领先阶、共同背景、正规性及作用假设'
            row['not_proved']='从认知选择该作用、全部高阶或完整量子过程尚未证明'
            row['source_rounds'] += [1027] if row['id']=='scalar_source_to_parent' else [1028,1029]
            row['source_paths'] += [f'archive_1009_/research_note_{n}.md' for n in row['source_rounds'][1:]]
    for n,ident,status,obj,gap in NEW_ROWS:
        rows.append(dict(id=ident,source_rounds=[n],status=status,object=obj,
                         claim=gap,requires='见本轮共同条件矩阵及原轮完整前提',
                         not_proved='不自动给认知选择或完整父过程',
                         source_paths=[f'archive_1009_/research_note_{n}.md'],
                         cognitive_origin_closed=False,full_parent_process_certified=False))
    return dict(round=1035,kind='main_agent_audited_typed_dependency_ledger',
                base_ledger=OLD,base_sha256=sha(BASE/OLD),rows=rows,
                locally_mapped_old_gaps=list(revisions),
                new_adopted_cognitive_axioms=0,full_parent_process_certified=False)


def run():
    return dict(round=1035,date='2026-10-08',new_calibration_groups=1,
                cumulative_research_groups=3812,new_adopted_cognitive_axioms=0,
                code_sha256=sha(Path(__file__)),
                historical_source_sha256={p:sha(BASE/p) for p in HISTORY},
                anomaly_polynomial_coefficients=polynomial_anomalies(),
                family=joint_family(),representation_words=representation_words(),
                mass_cases=mass_blocks(),bridge_rows=len(ledger()['rows']),
                original_P981_replaced=False,full_mixed_anomalies_certified=False,
                actual_physical_predictions_matched=False,
                external_independent_agent_review_completed=False,goal_completed=False)


def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:compare(a[k],b[k])
    elif isinstance(a,(list,tuple)):
        assert len(a)==len(b)
        for x,y in zip(a,b):compare(x,y)
    elif isinstance(a,float):assert np.isclose(a,b,atol=1e-12,rtol=1e-10),(a,b)
    else:assert a==b,(a,b)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true')
    args=p.parse_args();r=run();l=ledger()
    if args.write:
        assert not OUT.exists() and not LEDGER.exists()
        for path,data in ((OUT,r),(LEDGER,l)):
            with path.open('x',encoding='utf8') as f:
                json.dump(data,f,ensure_ascii=False,indent=2,allow_nan=False);f.write('\n')
    else:
        compare(r,json.loads(OUT.read_text('utf8')))
        compare(l,json.loads(LEDGER.read_text('utf8')))
    print(json.dumps(dict(round=1035,family_cases=len(r['family']),
                         representation_words=r['representation_words'],
                         mass_cases=len(r['mass_cases']),bridge_rows=len(l['rows'])),ensure_ascii=False))
