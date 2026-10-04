"""Prepare751 source-retention audit and cognitive joint-hypothesis ledger."""
import hashlib,json,re
from pathlib import Path
HERE=Path(__file__).resolve().parent

def write(name,text):
    p=HERE/name;p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)

def remap(text,mapping):
    return re.sub('|'.join(re.escape(k) for k in sorted(mapping,key=len,reverse=True)),
                  lambda m:mapping[m[0]],text)

ledger=(HERE/'unified_physics_condition_ledger_750.md').read_text('utf8')
ledger=ledger.replace('# 联合条件总账：750实际标记来源与共同概率',
                      '# 联合条件总账：751共同报告、内部差异与来源资源',1)
ledger=ledger.replace(ledger.split('\n')[2],
    '2026-10-04。接[750全账](unified_physics_condition_ledger_750.md)，回填[751报告](research_note_751.md)。[结果](joint_record_source_retention_results.json)、[核验](research_round_751_checks.json)。联合目标保持。',1)
updates={
'C03':'751原占据报告不充分决定完整来源；两未知输入的全cq后态相同而注能相反，理想instrument与真实内部实现区分',
'C05':'751索引更正：747原跳跃对远端读口的六阶分量已严格认证；748合成原标量边得总系数数值同号，严格总符号仍缺',
'C10':'751报告一致可与不同实际来源共存，不能把共同报告直接视为完整几何源；未建动态量子几何',
'C21':'751在原非零Majorana源块上，精确Luders、固定独立装置、端点可加能源守恒的组合被排除；有限精度或保端点相互作用等分支保留',
'C22':'751同一128维原来源的24符号保留下来的未知输入差与记录删除的成本相干共同核验；共同末态不足回推资源成本'}
rows=ledger.splitlines();seen=[]
for i,row in enumerate(rows):
    if row.startswith('|C15 '):
        parts=row.split('|');cut=parts[2].find('；744')
        assert cut>=0 and '748' in parts[2][cut:]
        parts[2]=parts[2][:cut]+'；751索引更正：744—748读出／传播移回C03／C05，不计作手征或三代问题完成'
        rows[i]='|'.join(parts)
    for key,value in updates.items():
        if row.startswith('|'+key+' '):
            parts=rows[i].split('|');parts[2]+='；'+value;rows[i]='|'.join(parts);seen.append(key)
assert set(seen)==set(updates)
ledger='\n'.join(rows)+'\n'
ledger=ledger.rsplit('## 本轮合并与下一项',1)[0]+'''## 751共同报告与完整过程的对象限制

- 对原占据读口，未知输入扩展必须另验；750已知参考的来源表不能直接当任意未知输入的普遍反馈律。
- 原规范中性、偶宇称子空间里，报告恒为0仍保留非标量能源与singlet来源；此校准的仅报告能源误差至少0.129229449417。
- 两原输入的整个cq后态相同，实际注能却为正负0.091138464916；共同来源／资源账不能仅从该后态追回。
- 原Majorana非零矩阵元使成熟WAY精确重复测量限制适用到限定实现类；增加平均电池账不构成充分完成。
- 以上只是有限原符号的对象见证与限定解析实施障碍，不构造全量子GR约束态，不否定近似记录或含端点相互作用的原候选。
- 工作假说允许共同报告R、未共享量子部门Q、内部资源及关联A共同存在；认知观察提供独立动机，数学与物理继续检验。不是宣称假说已由FUCP推出或已统一所有部门。
- C15仍为物种手征；读出归C03，传播归C05。[更正记录](round751_drafts/ledger_mapping_erratum.md)保留冻结旧账。
- [障碍／假说报告](round751_drafts/joint_hypothesis_obstacle_review.md)区分具体失败和未完成证明，保留Q／E两类候选，不把所有补丁相加。

## 本轮合并与下一项

C03／C10／C21／C22的输入量词、来源信息与内部补偿必须一起审计。净限制是排除“报告／物质后态已经包含全部未来来源与成本”的指定接法；维数、群、三代、共同尺度和引力约束未因此解决。

接[752](round752_drafts/STATUS.md)：比较同一量子整体与有效经典几何方案能否接入同一过程、来源、约束和尺度；优先减少独立分支和输入，不再扩充读口精度反例。旧空间、604、649／699边界保持。
'''
write('unified_physics_condition_ledger_751.md',ledger)
write('round751_drafts/research_note_751_draft.md',(HERE/'research_note_751.md').read_text('utf8'))
write('round752_drafts/STATUS.md','''# 第752轮入口：先检验联合假说组合的共同对象

接[751](../research_note_751.md)、[条件账](../unified_physics_condition_ledger_751.md)、[工作假说报告](../round751_drafts/joint_hypothesis_obstacle_review.md)。

1. 用户方法建议保持：障碍分类、提炼共同原因、相容假说组合、联合模型、检验充分性。认知观察可以独立提出工作假说，不要求先由FUCP推出。
2. 751证明原记录／完整cq后态不含所有未知输入来源与成本；不要重做无广播、WAY或更精细的读口扫描。
3. 比较Q方案（同一量子整体、共同报告子系统）与E方案（受控有效经典几何）。需共同准备、未知输入、全部来源、内部资源及约束；不是只分别写一份方程。
4. 先复用590／594量子几何控制与内部相互作用、649约束商失败、702—706同一有限过程、724—725实际区域、730—741条件连续源；599/699等反常和正性限制不能遗漏。
5. 一个真正有用的新增映射应同步连接C01/C03/C10/C11/C20/C22，或严格排除一套明确组合。590正控制Hamiltonian不是Einstein约束；741一阶响应不是精确共同量子几何。
6. 空间382—386、425、522—523直接复用，只补实际共同对象所缺接口，不重加已消去条件。C15专指物种手征，通信为C05。
7. 若暂时只有整理，保存工作报告而不虚增完成轮次。目标不改，不创建任务、定时或图像。
''')
write('round751_drafts/literature_scope_audit.json',json.dumps(dict(
    sources=[dict(title='Ozawa, Conservation laws, uncertainty relations, and quantum limits of measurements',
                  url='https://arxiv.org/html/quant-ph/0112154',
                  checked='WAY review and equations6-7: additive conservation, repeatability versus accurate probabilities, independent apparatus.',
                  use='Mature theorem only; original finite-source matrix-element obstruction independently derived. No general quantitative-error or unbounded GR theorem imported.')],
    inherited='252 common record/internal phase;350 mean feedback;593/634 compensation;594 interacting endpoint;649/699 restricted failures;730/732 original sources;750 fixed reference.',
    new='Actual source-retention and erased-cost witnesses in original gauge-even sterile sector; restrictive unknown-input completion test.',
    excluded='All physical GR instruments; exact continuum detector; generation of SM groups or space dimension; cognitive assumptions as already proven universals.'
),ensure_ascii=False,indent=2)+'\n')
write('round751_drafts/scope_and_dedup_review.md','''# 751范围、去重与方法衔接审查

本轮重读三导航、750笔记和结果、751入口、634／593及关键原代码；检索全部编号笔记相关关键词。Get-CimInstance过程命令权限不足，改用Get-Process未发现Python进程；研究代码复算无执行失败。没有使用子代理。

- 252曾有同标签相位影响图输出；751不把该一般现象当新发现，新增原物质完整来源与资源补偿的明确矩阵证据。
- 350的非线性平均反馈、593／634供能与跳换源、590／594实际相互作用、649物理商、699负测度均复用，不扩大失败范围。
- 新量词为同一固定原仪器的未知输入扩展；与750已知纯参考的校准不同，不推翻该旧公式。
- 128维是完整原局部符号；3个sterile模式的8维和2模的4维计算是矩阵元压缩，不是另造物质模型，也不是封闭动力学或全GR物理子空间。
- 所有源差以完整矩阵迹核对；码外原耦合存在，瞬时矩阵元之外不声称压缩H代替完整H。
- 第二见证的全cq相同而成本相反，不只报告概率相同；负注能允许未知非真空输入，不改变原真空注能正结论。
- WAY为成熟结论；当前矩阵元证明限定固定独立装置、精确Luders、可加端点总能源守恒及有限期望域，不将普通QM的总能源移成任意GR全局能源。
- 最小扩张保全局相位但不守恒指定可加能源；其结果是诊断而非成功装置。
- 认知观察仅动机／候选假说；报告之外保留量子与资源不自动补完手征连续、引力约束或跨尺度。
- 冻结751入口不修改。旧C15通信归类在新账明确更正，旧史不覆盖。整体目标和定时设置不改。
''')
main=('research_note_751.md','joint_record_source_retention.py','joint_record_source_retention_results.json','unified_physics_condition_ledger_751.md')
write('round751_drafts/final_review.txt',
      'Primary-agent review only. Original source sufficiency and erased-cost witnesses; known WAY theorem explicitly attributed; fixed-reference and unknown-input tasks separate. Cognitive hypotheses not declared proven physics. Frozen histories unchanged.\n'+
      '\n'.join(n+' '+hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in main)+'\n')
mapping={'751':'752','750':'751','749':'750','3439':'3441','3436':'3439','1562':'1565',
         '3506':'3521','3493':'3506','396':'397','301':'302',
         'joint_marked_source_completion':'joint_record_source_retention'}
pub=remap((HERE/'publish_round750.py').read_text('utf8'),mapping)
summary='**第751轮完成：** [共同报告与内部来源账]({p}research_note_751.md)原未知输入可有相同报告而不同来源，甚至同一完整cq后态却注能相反；限定精确端点可加实现受原Majorana的守恒障碍。保留有限精度／联合关联方案，不扩大到全GR。两组、十六式通过，最新751／3441，1565份编号科学文件、3521份保护证据。[核验]({p}research_round_751_checks.json)、[条件账]({p}unified_physics_condition_ledger_751.md)。'
order='**当前执行顺序（751后，优先于下方历史安排）：** 接[752联合假说的共同对象]({p}round752_drafts/STATUS.md)，沿障碍分类、共同原因、假说组合审计统一候选；复用旧源／约束边界，不再延伸读口精度反例。[工作报告]({p}round751_drafts/joint_hypothesis_obstacle_review.md)。旧空间、604、649／699与目标保持。'
pub=re.sub(r'^summary=.*$',lambda m:'summary='+repr(summary),pub,flags=re.M)
pub=re.sub(r'^order=.*$',lambda m:'order='+repr(order),pub,flags=re.M)
pub=pub.replace('实际标记来源与共同概率','共同报告与内部来源账').replace('旧空间合同保持，原条件源与共同概率相容','旧空间合同保持，联合假说与完整过程审计')
write('publish_round751.py',pub)
write('postcheck_round751.py',remap((HERE/'postcheck_round750.py').read_text('utf8'),mapping))
ver=remap((HERE/'verify_round750.py').read_text('utf8'),mapping)
start=ver.index('    names=(');end=ver.index('    new=',start)
names=('unified_physics_condition_ledger_751.md','round751_drafts/research_note_751_draft.md',
       'round751_drafts/final_review.txt','round751_drafts/literature_scope_audit.json',
       'round751_drafts/scope_and_dedup_review.md','round752_drafts/STATUS.md',
       'round751_drafts/record_source_sufficiency.py','round751_drafts/record_source_sufficiency_results.json',
       'round751_drafts/erased_source_balance.py','round751_drafts/erased_source_balance_results.json',
       'round751_drafts/joint_hypothesis_obstacle_review.md','round751_drafts/ledger_mapping_erratum.md')
ver=ver[:start]+'    names='+repr(names)+'\n'+ver[end:]
ver=ver.replace('==(3,0,0)','==(2,0,0)').replace('fresh_tests=dict(run=3','fresh_tests=dict(run=2')
ver=ver.replace('actual_marked_source_and_joint_normalization_checked=True,nonlinear_common_geometry_realized=False',
                'original_source_retention_and_conditional_completion_checked=True,nonlinear_common_geometry_realized=False')
write('verify_round751.py',ver)
print('Prepared751 source-retention and joint-hypothesis audit.')

