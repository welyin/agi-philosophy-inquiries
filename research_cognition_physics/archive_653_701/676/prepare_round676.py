"""Prepare676; exact finite seed/support audit with historical preservation."""
import hashlib
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent

def write(name,text):
    p=HERE/name;p.parent.mkdir(exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)

def replace(text,mapping):
    pats=[]
    for key in sorted(mapping,key=len,reverse=True):
        p=re.escape(key)
        if key.isdecimal():p=r'(?<!\d)'+p+r'(?!\d)'
        pats.append(p)
    return re.sub('|'.join(pats),lambda m:mapping[m.group()],text)

ledger=(HERE/'unified_physics_condition_ledger_675.md').read_text('utf8')
ledger=ledger.replace('# 联合条件总账：675实际标量边缘、辅助平均与共同正性','# 联合条件总账：676非平坦边界与原泛函零支撑')
ledger=ledger.replace('接[674全账](unified_physics_condition_ledger_674.md)，回填[675报告](research_note_675.md)。[结果](joint_physical_boundary_average_results.json)、[核验](research_round_675_checks.json)。',
    '接[675全账](unified_physics_condition_ledger_675.md)，回填[676报告](research_note_676.md)。[结果](joint_nonflat_seed_support_results.json)、[核验](research_round_676_checks.json)。')
ledger=ledger.replace('675实际标量边缘、精确辅助多项式及大τ必要条件|',
    '675实际标量边缘、精确辅助多项式及大τ必要条件；676指定非平坦边界表示的秩障碍与原零支撑|')
ledger=ledger.replace('实际无限动力学存在|','676第五方向正转移不等于物理时间；实际无限动力学存在|')
ledger=ledger.replace('任意费米来源反射、动态正性、指定归一与连续局域性仍缺',
    '676保留两方向时的spin分区不平衡使全部纯物理来源为零；任意费米来源反射、动态正性、指定归一与连续局域性仍缺')
ledger=ledger.replace('实际空间缩放与原物质嵌入连接，完整动力学收敛及接触重整化',
    '676固定Jminus边界像不能覆盖部分非平坦目标谱，须保留原泛函零支撑；实际空间缩放与原物质嵌入连接，完整动力学收敛及接触重整化')
at=ledger.index('## 本轮合并与下一项')
ledger=ledger[:at]+'''- 676在当前保留γ1、γ4传播的有限盒中，旧固定Jminus初始空间与目标负谱的S分区秩可能不同；整数通量证书给72/56对64/64，边界投影误差精确为1。Wplus/minus与第五方向T仍正，正性不能替代初始重叠合同。
- 同一分区不平衡使辅助配对留下不可由物理来源饱和的零模。原候选在此完整扇区的所有纯物理来源系数为零，故不能把投影表示失败升级成物理边缘反例。
- 该S来自两条未用spin方向，不推广到任意四维规范背景。非空开集只在当前保留的两方向、全部原内部群的配置空间中声明。

## 本轮合并与下一项

676连接C01、C04、C14、C16、C17、C20：指定非平坦边界表示的秩问题与原物理来源的支撑同步审计。旧强场近零权重和来源数值有明确的未饱和辅助零模机制，不能用其商构造归一物理响应。原矩形逆传播子无关字典和675有限边缘继续有效。

整数矩阵有理合同给精确存在证书；原浮点强场样本只是同机制诊断。既没有推翻整个domain-wall方法，也未证明现有局部体减除已经正确保留原零支撑。总Gauss积分的正性、严格归一、原H_F身份和物理时间仍缺。

接[677](round677_drafts/STATUS.md)：检验不依赖固定初始空间的全谱有理符号近似，保持原固定尺寸观测、相位和零支撑，结合673配置支配审查完整平均与极限。第五方向不改名物理时间，量子GR及四分支共同极限保持开放。目标及旧空间合同不改。
'''
write('unified_physics_condition_ledger_676.md',ledger)
mapping={'joint_physical_boundary_average':'joint_nonflat_seed_support',
         '674':'675','675':'676','676':'677','3255':'3257','3257':'3259',
         '1334':'1337','1337':'1340','2433':'2444','2444':'2460'}
verify=replace((HERE/'verify_round675.py').read_text('utf8'),mapping)
verify=verify.replace('Verify physical marginal, original nonflat auxiliary contraction and ground-kernel scope.',
    'Verify exact nonflat seed obstruction and original physical null support.')
verify=verify.replace('import joint_gauss_reflection_reality as prior_model','import joint_physical_boundary_average as prior_model')
a=verify.index('    import importlib.util');b=verify.index('    result=core.read(model.TARGET)',a)
verify=verify[:a]+verify[b:]
verify=verify.replace('    assert result==model.run()',"""    assert result==model.run()
    section=result['positive_fifth_transfer_and_fixed_seed_obstruction']
    assert section['transfer']==core.read(HERE/'round676_drafts/nonflat_boundary_probe_results.json')
    entry=core.read(HERE/'round676_drafts/exact_flux_seed_results.json')
    assert section['certificate']['exact_sectors']==entry['exact_sectors']
    assert section['certificate']['negative_counts']==[72,56]
""")
a=verify.index('    names=');b=verify.index('    preserved=',a)
verify=verify[:a]+"""    names=('unified_physics_condition_ledger_676.md','round676_drafts/research_note_676_draft.md',
           'round676_drafts/final_review.txt','round677_drafts/STATUS.md',
           'round676_drafts/nonflat_boundary_probe.py','round676_drafts/nonflat_boundary_probe_results.json',
           'round676_drafts/seed_sector_probe.py','round676_drafts/seed_sector_probe_results.json',
           'round676_drafts/exact_flux_seed.py','round676_drafts/exact_flux_seed_results.json',
           'round676_drafts/joint_nonflat_seed_support_before_json_key_fix.py',
           'round676_drafts/final_review_before_json_key_fix.txt','round676_drafts/reproduction_failure.json')
"""+verify[b:]
verify=verify.replace("assert text['display_formulas']==14","assert text['display_formulas']==12")
write('verify_round676.py',verify)
publish=replace((HERE/'publish_round675.py').read_text('utf8'),mapping|{'## 321.':'## 322.','## 226.':'## 227.'})
a=publish.index('summary=');b=publish.index('planned={}',a)
publish=publish[:a]+'''summary=('**第676轮完成：** [非平坦边界初始空间与原物理泛函的零权重支撑]({p}research_note_676.md)'
         '原整数通量精确给72／56负谱分区，659固定64／64边界像误差始终为1；'
         '同一分区同时使辅助权重及全部纯物理来源为零。'
         '因此是指定表示障碍，不是完整物理边缘反例；第五方向正性保留。'
         '两组、十二式通过，最新676／3259，1340份编号科学文件、2460份保护证据。'
         '[核验]({p}research_round_676_checks.json)、[全条件账]({p}unified_physics_condition_ledger_676.md)。'
         '限当前两传播方向盒，完整物理时间、连续及量子GR仍开放。')
order=('**当前执行顺序（676后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[677全谱有理符号与共同正则表示]({p}round677_drafts/STATUS.md)，'
       '保留原观测、零支撑及完整平均，不把第五方向当物理时间；目标不改。')
'''+publish[b:]
publish=publish.replace('物理标量边缘、完整辅助平均与共同正性判据','非平坦边界初始空间与原物理泛函的零权重支撑')
publish=publish.replace('实际标量边缘与旧基态的共同正性限制','非平坦边界表示与原物理来源支撑')
publish=publish.replace('非平坦辅助收缩与实际共同时间的边界','固定初始空间失效与物理边缘反例的区别')
write('publish_round676.py',publish)
write('postcheck_round676.py',replace((HERE/'postcheck_round675.py').read_text('utf8'),mapping))
write('round676_drafts/research_note_676_draft.md',(HERE/'research_note_676.md').read_text('utf8'))
review='''676 primary-agent proof/code/scope review; no independent agent.
Previous675 completed and published,676 executed immediately. Progress.
All prior evidence preserved; no task/automation/install/image checks.
Old space382-386,425,522-523 inherited, goal unchanged.
Wilson shifts unitary imply Re X>=-I; chiral W>=1-a, so fifth transfer positive.
S=i gamma2 gamma3 is real symmetric and commutes with retained gamma1,gamma4
propagation and gamma5. This is a two-direction box, NOT generic4d symmetry.
M and inverse(2+aX) preserve S. Fixed seed has64/64 sector dimensions.
Exact flux phases in{1,-1,i}, full original16 charge multiplicities retained.
Integer basis Gram2I gives raw8x8 complex blocks; realification doubles inertia.
Exact Fraction congruence uses1x1 or offdiagonal2x2 pivots with correct signs;
all determinants nonzero, row-norm determinant bound certifies finite gap.
Exact negative ranks72/56 persist for0<a<1 because Ha nonzero continuous.
Actual boundary and target sector ranks differ, projection distance1 for allL.
Open Haar neighborhood means all original-group links in retained two directions.
Original arbitrary extra spacetime links not included; early label clarified.
S^T M+M S=0 is an exact local Clifford identity, independent of original E.
Auxiliary skew pairing only crosses S sectors; nullity>=abs(nminus-nplus).
673 pure physical source independent of auxiliary Grassmann b, so ALL pure
physical coefficients vanish on imbalanced sector; no division by small weights.
Original spin-epsilon mass obeys same skew pairing; physical NW grading
diag(Sv,-S_Jplus) implies additional nullity, full N>=2abs imbalance.
Random original examples are diagnostics; exact flux is existence certificate.
Boundary seed failure does not refute physical marginal which is zero there,
does not refute domain-wall methods, actual HF, cognition principles or GR.
No physical time positivity inferred from fifth-direction positivity.
2 groups,12 equations.677 must check full-spectrum bounded rational alternative,
complete sources/integrals and support; no mere layer-count optimization.
'''
for name in ('research_note_676.md','joint_nonflat_seed_support.py','joint_nonflat_seed_support_results.json','unified_physics_condition_ledger_676.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round676_drafts/final_review.txt',review)
print('676 preparation completed')

