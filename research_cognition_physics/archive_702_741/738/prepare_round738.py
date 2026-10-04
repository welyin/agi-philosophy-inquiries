"""Prepare738 complete mixed response and immutable publication artifacts."""
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


ledger=(HERE/'unified_physics_condition_ledger_737.md').read_text('utf8')
ledger='# 联合条件总账：738同一匹配真空的完整中性混合响应\n'+ledger.split('\n',1)[1]
ledger=ledger.replace('接[736全账](unified_physics_condition_ledger_736.md)，回填[737报告](research_note_737.md)。[结果](joint_response_null_completion_results.json)、[核验](research_round_737_checks.json)。',
                      '接[737全账](unified_physics_condition_ledger_737.md)，回填[738报告](research_note_738.md)。[结果](joint_mixed_neutral_response_results.json)、[核验](research_round_738_checks.json)。')
marker='## 当前共同对象及仍存在的分支';assert marker in ledger
ledger=ledger.replace(marker,'**738当前增量：** 632同一匹配真空的完整两中性来源谱已接通，保原所有带荷物种及两个Majorana的混合阈值。矩阵对数主部正定，完整记忆有短时因果逆；737的局部互补合同因此已有真实活跃核，全部约束与实际动态反馈仍开放。\n\n'+marker)
updates={'C04':'738原两中性活跃矩阵响应有短时因果逆；全局部项及全约束另接',
         'C19':'738使用632同一匹配真空、原复耦合与整代半权，不重置730实际参考',
         'C22':'738真实非对角Majorana通道补足共同缩放遗漏；矩阵对数系数与599热核同源'}
rows=ledger.splitlines();seen=[]
for i,row in enumerate(rows):
    for key,value in updates.items():
        if row.startswith('|'+key+' '):
            parts=row.split('|');parts[2]+='；'+value;rows[i]='|'.join(parts);seen.append(key)
assert set(seen)==set(updates)
ledger='\n'.join(rows)+'\n'
ledger=ledger.rsplit('## 本轮合并与下一项',1)[0]+'''## 738完整中性矩阵谱与因果逆

- 原632背景的两个独立中性质量顶点经同一Takagi合同处理，全部复耦合、颜色重数、Majorana半权保留。
- 实际混合成对阈值为m1+m2；共同缩放方向精确抵消该通道，不能用其标量谱代替两来源矩阵。
- 完整谱按正半定矩阵项相加；矩阵对数系数严格正定，并与599原导数Gram计数一致。
- 每个同／异质量阈值均保留，得到矩阵log主部、常数及局部可积记忆的精确拆分。
- 736逆对数合同推广为非对易矩阵范数合同，原T=.01给q≤.267598的充分证书。
- 此证书不把全部引力局部常数设零，不解决完整规范约束、实际变化参考、非线性自洽或长期稳定。

## 本轮合并与下一项

C04／C19／C22获得同一背景与同一整代物质的实际联合响应，统一目标仍开放。

接[739](round739_drafts/STATUS.md)：保同一背景合并两径向、共形、剪切及必要规范约束，核共同局部接触、有限系数和真实初值；停止逆核精度优化。旧空间、604、649／699保持。
'''
write('unified_physics_condition_ledger_738.md',ledger)
write('round738_drafts/research_note_738_draft.md',(HERE/'research_note_738.md').read_text('utf8'))
write('round739_drafts/STATUS.md','''# 第739轮入口：同一退迟方程的全局部项与约束

接[738](../research_note_738.md)、[条件账](../unified_physics_condition_ledger_738.md)。原632匹配真空的完整中性矩阵谱及活跃短时因果逆已核，737提供局部互补阶数合同。

1. 复用570、583、601、630—632、731—735及737—738；停止重复缩放谱、Takagi或混合阈值积分。
2. 在632同一常真空上保两径向、共形及五剪切；剪切质量必须取该背景，不能沿用630 q★诊断数值。
3. 把632全部接触、735同源守恒处方与总有限局部系数放回同一因果方程。区分低能降阶与完整高阶分支，不以额外导数作为有界小项。
4. 从共同协变方程核lapse／shift及相关规范条件；静态均匀一阶约束退化不等于非线性可积。
5. 明确有限空间模式、非零动量与完整连续空间估计的差别，不用逐模式逆替代统一PDE界。
6. 之后再拓展至730实际动态参考、记录来源和非线性发展。旧空间、604、649／699范围及统一目标保持。
''')
write('round738_drafts/literature_scope_audit.json',json.dumps(dict(sources=[dict(
    url='https://arxiv.org/html/1007.3664',authors='de Gouvea, Huang, Shalgar',
    locations='Section III; Section IV.1, equations IV.11-IV.12',
    checked='Mass-basis Majorana Higgs vertices need not be diagonal; unequal-mass pair production depends on phases.',
    scope='Only the mechanism is reused. No triplet sector, external diagnostic parameter, decay normalization or nonlinear spacetime theorem is imported.')],
    original_derivation='Canonical unequal-mass projectors and630 phase-space convention fix the complete original one-generation spectrum. The logarithm-plus-memory matrix inverse extends736.',
    excluded='General interacting gauge response, full constraints, dynamic-background or nonlinear semiclassical existence.'),ensure_ascii=False,indent=2)+'\n')
write('round738_drafts/scope_and_dedup_review.md','''# 738谱、因果逆与范围审查

主代理审查，无独立代理，无图像检验。承接已正式完成737，不另开任务或改目标。

- 非对角Majorana顶点是成熟现象；原论文III与IV.1已核。代码从原复Y和632同一背景固定本项目归一。
- 保全部原带荷物种以及两个中性质量。中性有序i,j求和与总体半权配套，不再乘spin或Nambu重数。
- 式(5)由实际Dirac投影矩阵独立检查；解析谱不是用自身外积生成唯一验证。
- 共同缩放投影必须复现旧形式；混合通道抵消来自H T_H+S T_S=diag m，不能推广为任意中性变化。
- 对数矩阵正定来自原非零Y与精确Gram；有限采样不承担全频正性或右半平面无零点证明。
- 完整阈值gS、gP均在[0,1]且1-g≤2x²；余核绝对收敛和T²log界已逐项积分核对。
- 全部Cj不必对易；逆证明使用C正平方根的合同、标量K0卷积和矩阵算子范数。
- 数值q是纯非局部核证书，物理局部常数、互补变量及约束没有默设为零。
- 632常真空与730实际变化参考分开，所有旧空间、图连续及649／699缺口保持。
- 已保存两组结果；报告不宣称全部物理、广义相对论或完整非线性自洽已完成。
''')
main=('research_note_738.md','joint_mixed_neutral_response.py','joint_mixed_neutral_response_results.json','unified_physics_condition_ledger_738.md')
write('round738_drafts/final_review.txt','Primary-agent review only. Complete original neutral-source stationary fermion spectrum and exact-memory matrix inverse, not full constrained dynamics.\n'+
      '\n'.join(n+' '+hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in main)+'\n')
mapping={'738':'739','737':'738','736':'737','3413':'3415','3411':'3413',
         '1523':'1526','3327':'3336','3318':'3327','383':'384','288':'289',
         'joint_response_null_completion':'joint_mixed_neutral_response'}
publication=remap((HERE/'publish_round737.py').read_text('utf8'),mapping)
publication=re.sub(r'^summary=.*$',"summary='**第738轮完成：** [实际中性混合质量谱与完整矩阵因果逆]({p}research_note_738.md)原匹配真空的两中性来源补齐混合Majorana阈值；完整矩阵对数主部正定，保全部记忆的短时因果逆成立。两组、十六式通过，最新738／3415，1526份编号科学文件、3336份保护证据。[核验]({p}research_round_738_checks.json)、[条件账]({p}unified_physics_condition_ledger_738.md)。全局部项、约束及非线性自洽仍开放。'",publication,flags=re.M)
publication=re.sub(r'^order=.*$',"order='**当前执行顺序（738后，优先于下方历史安排）：** 接[739同一退迟方程与完整约束]({p}round739_drafts/STATUS.md)，在同一匹配真空合并两径向、共形、剪切及共同局部系数，核约束传播和实际初值；不再优化谱积分。旧空间、604、649／699及统一目标保持。'",publication,flags=re.M)
publication=publication.replace('响应零方向的局部补足与实际径向闭合','实际中性混合质量谱与完整矩阵因果逆')
publication=publication.replace('旧空间合同保持，局部补足不替代完整约束','旧空间合同保持，矩阵响应不替代完整约束')
write('publish_round738.py',publication)
write('postcheck_round738.py',remap((HERE/'postcheck_round737.py').read_text('utf8'),mapping))
verification=remap((HERE/'verify_round737.py').read_text('utf8'),mapping)
verification=verification.replace("checks['display_formulas']==14","checks['display_formulas']==16")
verification=verification.replace('inherited_original_metric_and_mass_maps_checked=True','original_unequal_mass_projectors_and_common_scaling_projection_checked=True')
write('verify_round738.py',verification)
print('Prepared738 original mixed-response report, ledger, review and next scope.')
