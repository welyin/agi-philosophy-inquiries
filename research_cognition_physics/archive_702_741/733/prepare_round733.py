"""Prepare733 evidence without modifying frozen preceding artifacts."""
import hashlib
import json
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent

def write(name,text):
    p=HERE/name;p.parent.mkdir(exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)

def remap(text,mapping):
    return re.sub('|'.join(re.escape(k) for k in sorted(mapping,key=len,reverse=True)),lambda m:mapping[m.group()],text)

ledger=(HERE/'unified_physics_condition_ledger_732.md').read_text('utf8')
ledger='# 联合条件总账：733有限记录模式的非线性边界与完整参考响应\n'+ledger.split('\n',1)[1]
ledger=ledger.replace('接[731全账](unified_physics_condition_ledger_731.md)，回填[732报告](research_note_732.md)。[结果](joint_relative_source_development_results.json)、[核验](research_round_732_checks.json)。',
 '接[732全账](unified_physics_condition_ledger_732.md)，回填[733报告](research_note_733.md)。[结果](joint_reference_polarization_boundary_results.json)、[核验](research_round_733_checks.json)。')
marker='## 当前共同对象及仍存在的分支';assert marker in ledger
ledger=ledger.replace(marker,'**733当前增量：** 合法Hadamard参考可精确保原有限记录块；固定局部操作字典下，明确非共形剪切的参考主符号变化非紧，有限光滑补丁不能吸收。完整分支反馈须保参考及共同减除响应；732领先相对解仍有效。此结论不否定有效约化，不从非紧性单独推出指定重整化应力非零。\n\n'+marker)
updates={'C04':'733有限记录模式不能一般替代绝对参考反馈，完整因果来源须续接',
 'C19':'733保记录的混合Hadamard嵌入可行，但参考随背景的改变另有来源与准备问题',
 'C22':'733完整分支差含同背景记录差及跨背景参考响应；共同减除只在前项自动相消'}
rows=ledger.splitlines();seen=[]
for i,row in enumerate(rows):
    for key,value in updates.items():
        if row.startswith('|'+key+' '):
            parts=row.split('|');parts[2]+='；'+value;rows[i]='|'.join(parts);seen.append(key)
assert set(seen)==set(updates)
ledger='\n'.join(rows)+'\n'
ledger=ledger.rsplit('## 本轮合并与下一项',1)[0]+'''## 733有限记录与完整参考的反馈边界

- 纯原参考给秩至多4的不变记录块；新Hadamard参考的块压缩可保正性、自对偶、Hadamard及原记录差。一般为新混态，制备未完成。
- 保局部标架和半密度字典，原完整Clifford主部在明确体积保持剪切下给非零零阶协方差符号差。高频弱趋零波包证明其非紧，任何有限光滑补丁不能消除。
- 此剪切是允许的背景比较，不宣称731共形初始修正已生成它。纯共形的领先方向不变，不纳入该零阶障碍。
- 完整来源分解为同背景记录差与跨背景参考差。仅第一项自动消去同一局部反项；后一项仍需背景、态、减除核和局部有限项共同响应。
- 非紧性不单独证明某个指定应力分量非零，也不意味着无穷自由参数或不存在受控有效理论。
- 以量子源ε记账且几何差O(ε)时，参考响应从二阶进入；732一阶结果保持。若只缩小测量而不缩小绝对源，不可沿用此阶次。
- 三组原完整矩阵校准分别核主符号、有限补丁边界及原动态参考的非零二阶来源；不是连续重整化应力数值预测。

## 本轮合并与下一项

C04／C19／C22明确了有限记录初始化与绝对反馈之间不可省略的参考接口，统一目标仍开放。

接[734](round734_drafts/STATUS.md)：在同一过去准备下核完整因果参考响应及共同减除，区别动力响应和保记录的新准备；停止有限补丁、秩及精度优化。旧空间、604、649／699范围保持。
'''
write('unified_physics_condition_ledger_733.md',ledger)
write('round733_drafts/research_note_733_draft.md',(HERE/'research_note_733.md').read_text('utf8'))
write('round734_drafts/STATUS.md','''# 第734轮入口：同一过去准备的因果参考响应与联合来源

接[733](../research_note_733.md)、[全账](../unified_physics_condition_ledger_733.md)。原有限记录可合法嵌入新参考，但不能替代完整背景参考变化；732领先相对发展保持。

1. 回查580／599／601、623—625、630—635、703—704和730—733。一般Duhamel／Kubo及反项相消不重报新轮。
2. 先核同一过去准备下，原完整Dirac／Majorana参考对时空背景变化的因果响应。明确变形的时间支撑、初态项、旋量识别及实际记录算符。
3. 733块嵌入为保记录重选参考；其准备响应不能直接叫作共同过去的退迟响应。核这两个允许处方对实际源的差异。
4. 来源须同时含态、局部源算符、Hadamard减除核及有限局部项的变化。不要对奇异连续协方差取裸有限矩阵迹。
5. 原共同Ward及线性约束传播直接复用；只在补出实际全源响应、共同正则性或明确反例时再编号，不继续符号幅度或补丁秩扫描。
6. 区分同背景记录差、跨背景完整参考差与两份绝对解存在。保旧空间、604、649／699、全量子Gauss和物理仪器边界，目标不改。
''')
write('round733_drafts/literature_scope_audit.json',json.dumps(dict(sources=[
 dict(url='https://arxiv.org/html/2108.11630',locations='Equation6.6; Proposition6.8; Theorem6.9',
 use='Cauchy Hadamard projections have the sign projection of the Dirac principal symbol. Lower-order corrections and smooth same-background state differences retain this symbol.',
 caution='Original730 finite matrix potentials and self-dual identification are retained. Local symbol numerical projectors are not substituted for the dynamic Hadamard covariance.'),
 dict(url='https://arxiv.org/html/1504.01034',locations='Section3.3; Lemma3.20; Theorem3.23',
 use='Classical Einstein-Dirac-Maxwell Cauchy/variational comparison only.',
 caution='Not a theorem for complete project chiral non-Abelian/Majorana/curved-target quantum renormalized feedback.')],
 new='Nonconformal reference-symbol obstruction mapped to the original legitimate record embedding, with exact source decomposition and the original dynamic second-order reference response.',
 excluded='Full nonlinear quantum-gravity closure; automatic nonzero renormalized stress from covariance noncompactness; every metric deformation; failure of all effective theories.'
),ensure_ascii=False,indent=2)+'\n')
write('round733_drafts/scope_and_dedup_review.md','''# 733范围与重复性审查

主代理审查；无独立代理、无图像检验。

- 732与733入口冻结不改。阅读最新报告及580／599、601、630—635、663／667、730—732，并检索旧编号稿的参考主符号／非紧记录论证。
- 新结果是当前原记录嵌入与完整参考变化的共同边界，不是重新发现Hadamard理论或参考有代价。
- 纯性对有限K不变性必要；新参考允许混合。补丁光滑与新背景Hadamard相容。
- 正负谱与原Clifford符号一致，全部原质量和连接保留为低阶项。
- 明确剪切给允许的反例，不称它是731共形初值生成的实际未来。纯共形不具有该零阶符号差。
- 连续非紧证明靠局部高频波包。有限矩阵校准不能代替连续证明。
- 非紧协方差不单独推出任何指定重整化应力非零，不能推出全统一计划失败。有限参数有效描述仍可研究。
- 二阶仅在源强ε且背景差O(ε)及连续来源可微时成立。矩阵族光滑不证明后者；732领先相对结果保持。
- 733保记录块的准备族不是天然退迟动力响应；下一项必须区分初态准备与同一过去演化。
- 原数值所有来源均为有限矩阵量，不作为重整化连续物理预测。经典Einstein-Dirac存在结果不覆盖本量子来源接口。
- 本轮没有新增认知公理、修改原目标或重新要求旧空间已消去的Lipschitz假设。
''')
main=('research_note_733.md','joint_reference_polarization_boundary.py','joint_reference_polarization_boundary_results.json',
 'unified_physics_condition_ledger_733.md')
write('round733_drafts/final_review.txt','733 primary review; no independent agent review.\n'+
 '\n'.join(n+' '+hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in main)+'\n'+
 'Three checks pass. Nonconformal reference-symbol boundary, legitimate finite record patch, full branch decomposition and original dynamic second-order reference source. Finite-symbol data are not renormalized stress predictions. Goal remains open.\n')
v=(HERE/'verify_round732.py').read_text('utf8')
v=remap(v,{'joint_relative_source_development':'joint_reference_polarization_boundary',
 'joint_source_constraint_response':'joint_relative_source_development',
 'round733':'round734','round732':'round733','round731':'round732',
 '_732':'_733','_731':'_732','range(584,732)':'range(584,733)',
 '3242':'3256','3256':'3270','round=732':'round=733','round=731':'round=732',
 '3401':'3404','1508':'1511','source_feedback_entry':'reference_embedding_entry'})
write('verify_round733.py',v.replace('Reproduce732','Reproduce733'))
post=(HERE/'postcheck_round732.py').read_text('utf8')
post=remap(post,{'range(584,733)':'range(584,734)','round732':'round733','_732':'_733',
 '3256':'3270','第732':'第733','round=732':'round=733'}).replace('Check732','Check733')
write('postcheck_round733.py',post)
pub=(HERE/'publish_round732.py').read_text('utf8')
pub=remap(pub,{'733':'734','732':'733','731':'732','3401':'3404','3398':'3401',
 '1508':'1511','3256':'3270','## 378.':'## 379.','## 283.':'## 284.',
 'joint_relative_source_development':'joint_reference_polarization_boundary',
 '同一记录的相对量子来源与短时联合发展':'有限记录模式的非线性边界与完整参考响应',
 '旧空间合同保持，相对线性发展不替代完整非线性反馈':'旧空间合同保持，完整参考响应不由有限记录替代'})
lines=pub.splitlines()
for i,line in enumerate(lines):
    if line.startswith('summary='):
        lines[i]="summary='**第733轮完成：** [有限记录模式的非线性边界与完整参考响应]({p}research_note_733.md)有限记录可合法嵌入新参考，但明确非共形变化的完整参考差非紧，不能用有限光滑模式吸收；跨背景反馈须计参考及减除响应。三组、十八式通过，最新733／3404，1511份编号科学文件、3270份保护证据。[核验]({p}research_round_733_checks.json)、[全条件账]({p}unified_physics_condition_ledger_733.md)。732一阶结果保持，绝对反馈开放。'"
    if line.startswith('order='):
        lines[i]="order='**当前执行顺序（733后，优先于下方历史安排）：** 接[734同一过去准备的因果参考响应]({p}round734_drafts/STATUS.md)，核完整态、源算符、共同减除及初态项；区别保记录的新准备与真实退迟响应。停止有限补丁、秩及精度优化；旧空间、604、649／699及统一目标保持。'"
write('publish_round733.py','\n'.join(lines)+'\n')
