"""Prepare740 without changing any prior evidence or goal."""
import hashlib
import json
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent
def write(name,text):
    p=HERE/name;p.parent.mkdir(exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)
def remap(text,values):
    return re.sub('|'.join(re.escape(k) for k in sorted(values,key=len,reverse=True)),lambda m:values[m.group()],text)

ledger=(HERE/'unified_physics_condition_ledger_739.md').read_text('utf8')
ledger='# 联合条件总账：740原增长极点与一致有限阶因果响应\n'+ledger.split('\n',1)[1]
ledger=ledger.replace('接[738全账](unified_physics_condition_ledger_738.md)，回填[739报告](research_note_739.md)。[结果](joint_covariant_response_closure_results.json)、[核验](research_round_739_checks.json)。',
    '接[739全账](unified_physics_condition_ledger_739.md)，回填[740报告](research_note_740.md)。[结果](joint_causal_eft_branch_results.json)、[核验](research_round_740_checks.json)。')
marker='## 当前共同对象及仍存在的分支';assert marker in ledger
ledger=ledger.replace(marker,'**740当前增量：** 原匹配真空的全频一费米圈重求和有规范不变增长极点。势与完整记忆同阶展开可保因果和线性引力约束，并有明确导数预算的二阶方程残差界；该界不是不稳定精确解误差。原校准脉冲的一圈／领先窗口比约51%，实际参数有效窗口尚未关闭。\n\n'+marker)
updates={'C04':'740一致有限阶保记忆、因果及约束；全频重求和原分支有增长极点',
         'C19':'740数学校准不等于实际记录来源；原实际参数共同小误差仍开放',
         'C22':'740完整共同圈阶有二阶方程残差界，不等于精确非线性自洽误差'}
rows=ledger.splitlines();seen=[]
for i,row in enumerate(rows):
    for key,value in updates.items():
        if row.startswith('|'+key+' '):
            parts=row.split('|');parts[2]+='；'+value;rows[i]='|'.join(parts);seen.append(key)
assert set(seen)==set(updates)
ledger='\n'.join(rows)+'\n'
ledger=ledger.rsplit('## 本轮合并与下一项',1)[0]+'''## 740原总极点及有限阶合同

- 全频精确重求和的原分支至少两次零模交叉；第一根约2241，第二根以log q约757保存，不作物理UV预测。
- 树级势与匹配一圈势分开计数，保原中性混合、共形交叉及剪切记忆。
- 仅逆经典退迟算符构造有限阶系数，来源相容时逐阶保持线性引力约束。
- 固定γ与有限K、六个来源导数预算下，二阶方程残差受控；完整重求和解及未知高圈的误差另计。
- 实际λ=1校准的一圈修正约51%；缩小来源幅值不能改变相对修正。未认证共同物理窗口。
- 不重复601一般高阶警告、732相对源发展或739协变投影。

## 本轮合并与下一项

C04／C19／C22的有效阶次及约束连接进一步收紧，统一目标未完成。

接[741](round741_drafts/STATUS.md)：回到730—735实际历史和绝对源，分清首阶求源与下一阶反馈，核实际初值、守恒与有效窗口。旧空间、604、649／699保持。
'''
write('unified_physics_condition_ledger_740.md',ledger)
write('round740_drafts/research_note_740_draft.md',(HERE/'research_note_740.md').read_text('utf8'))
write('round741_drafts/STATUS.md','''# 第741轮入口：实际历史与绝对来源的一致阶次

接[740](../research_note_740.md)、[条件账](../unified_physics_condition_ledger_740.md)。原全频重求和失败不等于一致有限阶失败；实际物理小误差仍未认证。

1. 回查570—574、601、630—635及730—735，复用相对源发展、实际参考、初始约束补偿和共同守恒处方；不重新证明普通线性波方程。
2. 在原实际背景、同一准备和完整物质上区分绝对源Q(B0)、完整响应DQ[b]及它们真正进入的圈阶，避免对首阶问题多做全频重求和。
3. 核绝对源与原初始Gauss／引力约束是否共同可解。绝对有限常数保持输入；不能以相对反项抵消冒充绝对分支。
4. 审核沿实际背景族及过去准备的正则性；固定光滑路径可微与未知非线性方程适定是不同命题。
5. 实际记录装置、参考能量、有效窗口及全量子／连续边界保留。不要扫脉冲幅值、积分阶数或有限常数凑稳定样本。
6. 当前目标不变，旧空间382—386、425、522—523，604及649／699全部保持。
''')
write('round740_drafts/literature_scope_audit.json',json.dumps(dict(sources=[dict(
    url='https://arxiv.org/html/gr-qc/9211002',authors='Parker and Simon',
    use='Consistent perturbative order versus nonperturbative branches; method only.',
    boundary='Their explicit local/state/conformal examples do not prove a massive mixed nonlocal-memory result.')],
    new='Original full-symbol growing poles; matched-potential and full-memory finite-order response; causal Ward-compatible derivative-budgeted equation residual; actual pulse calibration.',
    inherited='601 order-reduction distinction;632 matching;735 conserving prescription;738 full spectra;739 covariant symbols and constraints.',
    excluded='Uniform error relative to unstable exact inverse, physical lambda=1 accuracy, actual record apparatus, nonlinear unified existence.'),ensure_ascii=False,indent=2)+'\n')
write('round740_drafts/scope_and_dedup_review.md','''# 740范围、证明与去重审查

前一目标轮次完成739并执行740入口，归类为进展。当前复用冻结入口，不覆盖其5份文件。主代理审查，无代理或图像检验。

- 入口的惯性证明承担极点存在性，浮点根不称区间认证。全频一圈近似分支失败不扩大为微观物理或认知原则失败。
- λ同时乘一圈势与响应；没有把总势无限重求和伪称严格圈阶。
- P1保质量—度规交叉、原完整混合阈值及剪切。经典正质量给P0在右半平面退迟逆。
- 固定γ、K的粗五阶乘子界留导数余量，未恢复736排除的恰好四阶无损失界。
- 因果性由退迟算符支集证明；空间带限情形不滥称来源空间局域。
- 二阶方程残差不能推成不稳定精确解误差，未知高圈与UV另计。
- 数字51%是指定窗口范数比，不是所有来源算子范数；小幅度不改善相对圈修正。
- 初次160／256积分分辨率未通过原要求，增至384／512后通过，未放宽门槛；失败与修复保存在结果。
- 真实记录制备未由校准源构造；下一轮复用732、735，须有真实新增连接才计轮次。
- 旧空间、604、649／699与所有冻结文件保留，目标未改。
''')
main=('research_note_740.md','joint_causal_eft_branch.py','joint_causal_eft_branch_results.json','unified_physics_condition_ledger_740.md')
write('round740_drafts/final_review.txt','Primary-agent review only. Original poles and consistent causal finite-order response; equation residual is not solution error.\n'+
    '\n'.join(n+' '+hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in main)+'\n')
mapping={'740':'741','739':'740','738':'739','3417':'3419','3415':'3417','1529':'1532',
         '3345':'3359','3336':'3345','385':'386','290':'291',
         'joint_covariant_response_closure':'joint_causal_eft_branch'}
publication=remap((HERE/'publish_round739.py').read_text('utf8'),mapping)
publication=re.sub(r'^summary=.*$',"summary='**第740轮完成：** [原增长极点与一致有限阶因果响应]({p}research_note_740.md)原全频一圈重求和有增长极点；一致有限阶保完整记忆、因果与线性约束，并有二阶方程残差界。原校准脉冲修正约51%，实际物理窗口尚未认证。两组、十八式通过，最新740／3419，1532份编号科学文件、3359份保护证据。[核验]({p}research_round_740_checks.json)、[条件账]({p}unified_physics_condition_ledger_740.md)。'",publication,flags=re.M)
publication=re.sub(r'^order=.*$',"order='**当前执行顺序（740后，优先于下方历史安排）：** 接[741实际历史与绝对来源的共同阶次]({p}round741_drafts/STATUS.md)，复用730—735，核实际初值、守恒及有效窗口。方程残差不等于不稳定精确解误差。旧空间、604、649／699及统一目标保持。'",publication,flags=re.M)
publication=publication.replace('同一中性—几何退迟方程与有限空间频段的约束闭合','原总反馈的增长极点与一致有限阶因果响应')
publication=publication.replace('旧空间合同保持，线性约束不替代非线性自洽','旧空间合同保持，有限阶残差不替代完整解误差')
write('publish_round740.py',publication)
write('postcheck_round740.py',remap((HERE/'postcheck_round739.py').read_text('utf8'),mapping))
verification=remap((HERE/'verify_round739.py').read_text('utf8'),mapping)
verification=verification.replace("'round741_drafts/STATUS.md')","'round741_drafts/STATUS.md',\n           'round740_drafts/total_pole_entry.py','round740_drafts/total_pole_entry_results.json',\n           'round740_drafts/research_note_740_working.md','round740_drafts/check_and_publish_entry.py',\n           'round740_drafts/entry_checks.json')",1)
verification=verification.replace("    result=model.run();", "    entry=core.read(HERE/'round740_drafts/entry_checks.json')\n    for name,digest in entry['artifact_hashes'].items():assert core.digest(HERE/'round740_drafts'/name)==digest,name\n    result=model.run();")
verification=verification.replace('original_local_curvature_and_spatial_spectral_shift_checked=True','original_total_poles_and_causal_pulse_residual_checked=True')
write('verify_round740.py',verification)
print('Prepared740: finite-order proof, ledger, review, publication and next scope.')
