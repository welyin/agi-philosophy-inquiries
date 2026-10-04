"""Prepare752 conditional relational-interaction/common-cone compatibility record."""
import hashlib,json,re
from pathlib import Path
HERE=Path(__file__).resolve().parent

def write(name,text):
    p=HERE/name;p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)

def remap(text,mapping):
    return re.sub('|'.join(re.escape(k) for k in sorted(mapping,key=len,reverse=True)),
                  lambda m:mapping[m[0]],text)

ledger=(HERE/'unified_physics_condition_ledger_751.md').read_text('utf8')
ledger=ledger.replace(ledger.split('\n')[0],'# 联合条件总账：752关系定位、实际交互与共同光锥',1)
ledger=ledger.replace(ledger.split('\n')[2],
    '2026-10-04。接[751全账](unified_physics_condition_ledger_751.md)，回填[752报告](research_note_752.md)。[结果](joint_relational_interaction_cones_results.json)、[核验](research_round_752_checks.json)。联合目标保持。',1)
updates={
'C03':'752原关系窗口作为观测与作为实际作用区分；新增耦合只是候选，没有构造自治量子instrument',
'C05':'752指定导数窗口耦合下，两个径向极化精确保原度规锥要求窗口对A/B仿射；并非因果性本身要求精确共锥',
'C09':'752复用573/647满秩X=(h,s,A,B)，其导数部分参与作用会改变主符号；参考存在不等于可直接控制原动力学',
'C10':'752仿射窗口且有效目标度规正为此受限类的共锥相容方案；非零四标签紧支撑与强共锥合同不相容',
'C22':'752新增局部作用的度规来源与径向主符号由同一泛函共同变分；数值限原初始jet，未求新全约束解'}
rows=ledger.splitlines();seen=[]
for i,row in enumerate(rows):
    for key,value in updates.items():
        if row.startswith('|'+key+' '):
            parts=row.split('|');parts[2]+='；'+value;rows[i]='|'.join(parts);seen.append(key)
assert set(seen)==set(updates)
ledger='\n'.join(rows)+'\n'
ledger=ledger.rsplit('## 本轮合并与下一项',1)[0]+'''## 752相容类与限定排除

- 候选假说：内部参考参与指定实际交互；新增作用ΔL=-εW(s)f(h,s,A,B)。它不是647旧诊断泛函已经具备的动力学，也不是从FUCP推出。
- 强合同：每个g-null方向有两个独立径向物理极化。该合同强于共同有效度规近似，不能由微分同胚协变性替代。
- 原参考满秩、εW非零时，窗口的A/B Hessian必须为零。仿射类改变目标度规和势，正目标度规时保共同主锥。
- 非零四标签紧支撑窗口不在此类内。限制只作用于声明的局部作用接法，不否定原关系观测、其他操作或所有内部参考。
- 原573/730同一点的独立梯度Hessian、度规来源与特征根已核；局部动能和空间能量为正且径向速度不同于原度规速度。没有求新增模型的约束初值或全发展。
- 旧352一般多锥限制、573主部、647关系来源直接复用；此次新增为原参考窗口的允许类及实际数值映射。
- 路径积分590/623—625、全量子物理态和BFR纯引力文献范围已审计，不能拿既有Gram正性或形式BV直接签收完整量子引力。
- C15仍专指手征物种。空间、604、649/699、量子连续映射、共同资源及自由物理参数缺口不因此改变。

## 本轮合并与下一项

C03/C05/C09/C10/C22现在共享一项可核的交叉限制：内部定位一旦进入实际作用，其来源和传播须共同改算。减少的是独立指定定位耦合与共同光锥的自由，不是新增一条普遍认知公理。

接[753](round753_drafts/STATUS.md)：回到原完整动力学，比较关系窗口只作读出描述的候选和修改作用的候选，检验共同过程与来源的对象身份；不继续优化窗口或再造仪器小模型。统一目标保持。
'''
write('unified_physics_condition_ledger_752.md',ledger)
write('round752_drafts/research_note_752_draft.md',(HERE/'research_note_752.md').read_text('utf8'))
write('round753_drafts/STATUS.md','''# 第753轮入口：原过程、关系报告与实际局域交互

接[752](../research_note_752.md)、[条件账](../unified_physics_condition_ledger_752.md)、[去重审计](../round752_drafts/scope_and_dedup_review.md)。

1. 752完成的是指定接法的相容分类：原导数四参考进入局部窗口作用，又要求两径向精确共锥，只允许A/B仿射。不是所有认知交互失败，未求新全约束解。
2. 保留一条不修改原作用的候选：实际相互作用由原完整动力学承担，关系窗口用于事件描述与报告。不能把描述函数当作已实现的局部instrument。
3. 先回查590、624—625、647、652、723—725及730—751：双历史、区域读口、来源变分已有，直接复用；尤其751报告不是完整资源账。
4. 核“同一量子过程、共同来源与物理约束、内部参考与报告”的共同对象。BFR纯引力微扰路线有额外前提，不能直接声称已量子化原全物质和记录。
5. 认知观察可独立提出相容假说组合；继续区分限定失败、未完成证明和真正公理反证。C15是物种手征，传播为C05。
6. 不继续窗口剖面、锥差精度或读口资源优化。仅组织／成熟定理复述保存工作稿，不虚增完成轮次。
7. 旧空间382—386、425、522—523及604、649/699保持。不改目标，不设定时、不建任务、不做图像。
''')
write('round752_drafts/literature_scope_audit.json',json.dumps(dict(
    sources=[
      dict(title='Quantum gravity from the point of view of locally covariant quantum field theory',
           url='https://arxiv.org/html/1306.1058',
           checked='Sections3.2,5,6: on-shell expansion; positive physical free BRST representation premise; pure gravity treatment and proposed matter extensions.',
           use='Audit only. Not claimed as completed full original gravity-matter-record construction.'),
      dict(title='k-Essence, superluminal propagation, causality and emergent geometry',
           url='https://arxiv.org/html/0708.0561',
           checked='Section2: derivative Lagrangian principal tensor and effective causal cone; common spacelike initial surface condition.',
           use='Known method, independently differentiated for original two-radial target and new relational window. General multiple-cone fact already in352.')
    ],
    inherited='352 multiple cones;573 original classical constrained reference;590/623-625 double histories;647 relational sources;649 restricted quotient obstruction;730/731 actual background and sources;751 resource retention.',
    new='Conditional classification of one original-reference interaction candidate jointly with exact two-polarization cone contract; same action stress and Hessian audit on original jet.',
    excluded='No autonomous instrument, no new on-shell Einstein solution, no quantum-gravity physical state, no FUCP-wide no-go, no all-support stability certificate.'
),ensure_ascii=False,indent=2)+'\n')
main=('research_note_752.md','joint_relational_interaction_cones.py','joint_relational_interaction_cones_results.json','unified_physics_condition_ledger_752.md')
write('round752_drafts/final_review.txt',
      'Primary-agent review only. Explicit new action versus old diagnostic observable separated. Strong complete-polarization cone contract stated. Exact Hessian necessity and compact-fiber obstruction proved; positive affine class retained. Numerical original jet is not a new constrained solution. Historical results and whole goal unchanged.\n'+
      '\n'.join(n+' '+hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in main)+'\n')
mapping={'752':'753','751':'752','750':'751','3441':'3443','3439':'3441',
         '1565':'1568','3521':'3530','3506':'3521','397':'398','302':'303',
         'joint_record_source_retention':'joint_relational_interaction_cones'}
pub=remap((HERE/'publish_round751.py').read_text('utf8'),mapping)
summary='**第752轮完成：** [关系定位、实际交互与共同光锥]({p}research_note_752.md)原导数参考窗口若进入指定局部作用，两径向精确共锥只允许A/B仿射；非零四标签紧支撑接法被排除，正目标度规的仿射类保留。数值限原初始jet，未求新全约束解。两组、十六式通过，最新752／3443，1568份编号科学文件、3530份保护证据。[核验]({p}research_round_752_checks.json)、[条件账]({p}unified_physics_condition_ledger_752.md)。'
order='**当前执行顺序（752后，优先于下方历史安排）：** 接[753原过程与关系报告]({p}round753_drafts/STATUS.md)，回到原动力学检验描述性窗口、实际交互和共同来源；旧双历史及关系读口直接复用，不继续窗口优化。[范围审计]({p}round752_drafts/scope_and_dedup_review.md)。旧空间、604、649／699与目标保持。'
pub=re.sub(r'^summary=.*$',lambda m:'summary='+repr(summary),pub,flags=re.M)
pub=re.sub(r'^order=.*$',lambda m:'order='+repr(order),pub,flags=re.M)
pub=pub.replace('共同报告与内部来源账','关系定位、实际交互与共同光锥').replace('旧空间合同保持，联合假说与完整过程审计','旧空间合同保持，原参考交互与传播相容审计')
write('publish_round752.py',pub)
write('postcheck_round752.py',remap((HERE/'postcheck_round751.py').read_text('utf8'),mapping))
ver=remap((HERE/'verify_round751.py').read_text('utf8'),mapping)
start=ver.index('    names=(');end=ver.index('    new=',start)
names=('unified_physics_condition_ledger_752.md','round752_drafts/research_note_752_draft.md',
       'round752_drafts/final_review.txt','round752_drafts/literature_scope_audit.json',
       'round752_drafts/scope_and_dedup_review.md','round753_drafts/STATUS.md')
ver=ver[:start]+'    names='+repr(names)+'\n'+ver[end:]
ver=ver.replace('original_source_retention_and_conditional_completion_checked=True',
                'original_reference_interaction_and_common_cone_contract_checked=True')
write('verify_round752.py',ver)
print('Prepared752 compatibility classification and next joint-process entry.')

