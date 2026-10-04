"""Prepare737 preservation, ledger and publication without changing old files."""
import hashlib
import json
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent


def write(name,text):
    p=HERE/name;p.parent.mkdir(exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)


def remap(text,values):
    return re.sub('|'.join(re.escape(k) for k in sorted(values,key=len,reverse=True)),
                  lambda m:values[m.group()],text)


ledger=(HERE/'unified_physics_condition_ledger_736.md').read_text('utf8')
ledger='# 联合条件总账：737响应零方向的局部补足与实际径向闭合\n'+ledger.split('\n',1)[1]
ledger=ledger.replace('接[735全账](unified_physics_condition_ledger_735.md)，回填[736报告](research_note_736.md)。[结果](joint_causal_log_inverse_results.json)、[核验](research_round_736_checks.json)。',
                      '接[736全账](unified_physics_condition_ledger_736.md)，回填[737报告](research_note_737.md)。[结果](joint_response_null_completion_results.json)、[核验](research_round_737_checks.json)。')
marker='## 当前共同对象及仍存在的分支';assert marker in ledger
ledger=ledger.replace(marker,'**737当前增量：** 原谱的质量／共形零方向具有非零领先局部主部；R²有限项要求四阶组织。632匹配真空的共同缩放射线不受原标量动力学保持，下一项须保两个实际中性径向的混合谱。局部互补命题不替代约束消元。\n\n'+marker)
updates={'C04':'737限定局部互补命题；零／非零总四阶系数分别验收，非线性约束仍开放',
         'C19':'737不把630诊断射线替换为632匹配真空；同一背景与真实准备保持',
         'C22':'737排除单一共同质量缩放作为匹配真空的闭合动力学；需完整中性混合响应'}
rows=ledger.splitlines();seen=[]
for i,row in enumerate(rows):
    for key,value in updates.items():
        if row.startswith('|'+key+' '):
            parts=row.split('|');parts[2]+='；'+value;rows[i]='|'.join(parts);seen.append(key)
assert set(seen)==set(updates)
ledger='\n'.join(rows)+'\n'
ledger=ledger.rsplit('## 本轮合并与下一项',1)[0]+'''## 737局部补足与真实方程块

- 583联合二次型、623质量坐标及630—631谱零方向直接复用，新增三者的对象映射。
- 原F>0领先作用在非局部零方向有A=6F/M>0的二阶幅值；这是受限主部，未作全约束消元。
- 原632匹配真空的K^-1H不保持共同质量缩放射线；真空驻定不等于射线动力学闭合。
- R²总局部有限系数非零时，零方向是四阶；不能藏进按二阶未知量定义的有界余项。
- 明确二阶／四阶两种局部互补接法及连续因果合同；不保证系数趋零的统一极限或长期稳定。
- 静态均匀背景的一阶lapse约束退化不等于完整初值相容；570／601及731原约束成果继续保留。
- 不用有限质量Gram替代连续混合谱，也不把所有原标量坐标当作独立物理粒子。

## 本轮合并与下一项

C04／C19／C22共同限定可继续求解的真实方程块，统一目标仍开放。

接[738](round738_drafts/STATUS.md)：在632同一背景上保两个实际中性径向，核非交换质量顶点的混合阈值、对数主部和局部余项；之后接全约束及730变化参考。旧空间、604、649／699保持。
'''
write('unified_physics_condition_ledger_737.md',ledger)
write('round737_drafts/research_note_737_draft.md',(HERE/'research_note_737.md').read_text('utf8'))
write('round738_drafts/STATUS.md','''# 第738轮入口：同一匹配真空的完整中性质量响应

接[737](../research_note_737.md)、[条件账](../unified_physics_condition_ledger_737.md)。领先局部作用补足非局部零方向，但原匹配真空的共同缩放射线不闭合。

1. 直接复用583、623、630—632、730及735—737，停止重复局部Schur、共同缩放谱或高频参数扫描。
2. 保632的实际h、s真空、整代复杂Yukawa和两个Majorana质量；它与630的q★诊断背景分开。
3. 给出两个实际中性质量顶点的完整谱矩阵，包含不对易／混合质量跃迁；不能以两份同质量缩放谱代替。
4. 核正定对数系数、完整阈值余核及736逆核的矩阵推广条件。有限质量Gram只能认证代数注入，不能单独证明连续响应主部。
5. 同一来源接632局部接触与737二阶／四阶互补条件，明确有限常数、有效阶次和初始历史。
6. 不把中性受限部门当全部规范约束或实际730动态背景；后续再核这些连接。旧空间、604及649／699范围、统一目标均保持。
''')
write('round737_drafts/literature_scope_audit.json',json.dumps(dict(
    new_external_theorem_imported=False,
    reused_notes=[570,583,601,623,630,631,632,730,735,736],
    method='Direct algebra from the saved original action; the conditional inverse proof uses the already audited736 logarithmic causal inverse.',
    distinction='The derivative order and finite coefficients are not new cognitive principles. No fresh external nonlinear existence theorem is invoked.',
    excluded='Full mixed-mass continuum spectrum, physical mode count, complete constraints, general nonlinear semiclassical existence or global stability.'),ensure_ascii=False,indent=2)+'\n')
write('round737_drafts/scope_and_dedup_review.md','''# 737证明、去重与范围复核

主代理审查，无独立代理，无图像检查。先核736报告、入口、结果和目标，未修改目标或设置定时任务。

- 583已经给原共形／标量混合及Schur工具；本轮只新增其与当前质量响应零方向的映射，不称首次发现。
- x坐标及矩阵质量线性结构来自623／735；G由原K和Jacobian独立校准。
- A严格正只针对原领先作用；有限局部动能修正未默认为小，实际互补系数须非零。
- 匹配真空的两个中性径向仍保原L、u与632同一W Hessian。驻定条件不蕴含原动力学保持共同质量缩放射线。
- 解析闭合条件和原参数数值失败分开；没有将诊断不闭合推广为所有认知模型无解。
- R²二次变分由完整曲率作用独立求出；有限频率表不承担不有界证明。
- 零／非零总四阶系数对应不同最高导数，短时常数不保证在系数趋零时一致；601微扰降阶不被替代。
- 局部互补Volterra逆需a非零、L1余核及规定的活跃log核。全混合谱尚未证明，不用有限Gram冒充。
- Euclidean负方向、lapse约束、真空线性化及非线性可积性分清；不宣称已解全引力约束。
- 全部历史与冻结资料保持。两组代码核验是局部对象与阶数检查，不是时空数值求解。
''')
main=('research_note_737.md','joint_response_null_completion.py','joint_response_null_completion_results.json','unified_physics_condition_ledger_737.md')
write('round737_drafts/final_review.txt','Primary-agent review only. Original local complement and actual radial nonclosure; conditional causal criterion, no full constrained solution.\n'+
      '\n'.join(n+' '+hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in main)+'\n')

# Derive only NEW scripts. Simultaneous mapping prevents cascading round IDs.
mapping={'737':'738','736':'737','735':'736','3411':'3413','3409':'3411',
         '1520':'1523','3318':'3327','3304':'3318','382':'383','287':'288',
         'joint_causal_log_inverse':'joint_response_null_completion'}
publication=remap((HERE/'publish_round736.py').read_text('utf8'),mapping)
publication=re.sub(r"^summary=.*$", "summary='**第737轮完成：** [响应零方向的局部补足与实际径向闭合]({p}research_note_737.md)原领先局部作用补足谱零方向；R²有限项须按四阶组织。632匹配真空的共同缩放射线不闭合，须保实际中性混合响应。两组、十四式通过，最新737／3413，1523份编号科学文件、3327份保护证据。[核验]({p}research_round_737_checks.json)、[条件账]({p}unified_physics_condition_ledger_737.md)。完整约束与非线性自洽仍开放。'",publication,flags=re.M)
publication=re.sub(r"^order=.*$", "order='**当前执行顺序（737后，优先于下方历史安排）：** 接[738实际中性混合质量响应]({p}round738_drafts/STATUS.md)，保632同一匹配背景及两个真实径向，核混合阈值、矩阵对数主部和共同局部项；不能把单一缩放射线当闭合动力学。旧空间、604、649／699及统一目标保持。'",publication,flags=re.M)
publication=publication.replace('原完整谱的短时因果逆与反馈适用条件','响应零方向的局部补足与实际径向闭合')
publication=publication.replace('旧空间合同保持，短时逆不替代完整约束','旧空间合同保持，局部补足不替代完整约束')
write('publish_round737.py',publication)
write('postcheck_round737.py',remap((HERE/'postcheck_round736.py').read_text('utf8'),mapping))
print('Prepared737: ledger, frozen draft, review and next-round scope; old files untouched.')
