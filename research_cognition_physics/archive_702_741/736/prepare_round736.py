"""Prepare736 exact-memory inverse scope, ledger and publication scripts."""
import hashlib
import json
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent


def write(name,text):
    p=HERE/name;p.parent.mkdir(exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as handle:handle.write(text)


def remap(text,mapping):
    return re.sub('|'.join(re.escape(k) for k in sorted(mapping,key=len,reverse=True)),lambda m:mapping[m.group()],text)


ledger=(HERE/'unified_physics_condition_ledger_735.md').read_text('utf8')
ledger='# 联合条件总账：736原完整谱的短时因果逆\n'+ledger.split('\n',1)[1]
ledger=ledger.replace('接[734全账](unified_physics_condition_ledger_734.md)，回填[735报告](research_note_735.md)。[结果](joint_local_source_normalization_results.json)、[核验](research_round_735_checks.json)。',
    '接[735全账](unified_physics_condition_ledger_735.md)，回填[736报告](research_note_736.md)。[结果](joint_causal_log_inverse_results.json)、[核验](research_round_736_checks.json)。')
marker='## 当前共同对象及仍存在的分支';assert marker in ledger
ledger=ledger.replace(marker,'**736当前增量：** 原630—631完整阈值核精确分解为对数主部和局部可积记忆；活跃质量／剪切来源可给短时唯一因果逆。附加有限局部项有明确有界／压缩条件；原7来源块的非局部零方向、约束消元和730动态背景仍未签收。\n\n'+marker)
updates={'C04':'736定常活跃谱块有保完整记忆的短时因果逆；不等于一般动态自洽发展',
         'C19':'736反演限零过去来源差；既往历史必须保留，不重置实际准备',
         'C22':'736谱正性、精确记忆与局部项接入判据合并；非局部秩六不能直接逆全部七来源'}
rows=ledger.splitlines();seen=[]
for index,row in enumerate(rows):
    for key,value in updates.items():
        if row.startswith('|'+key+' '):
            parts=row.split('|');parts[2]+='；'+value;rows[index]='|'.join(parts);seen.append(key)
assert set(seen)==set(updates)
ledger='\n'.join(rows)+'\n'
ledger=ledger.rsplit('## 本轮合并与下一项',1)[0]+'''## 736精确记忆与因果逆

- 630—631原全部非零质量与半权保持；谱及运行系数不重算为新发现。
- 两个精确阈值核分别拆成c log s、固定阈值常数和L1因果余核，得到O(T² log(1/T))显式余核界。
- 归一原谱核是s²变量的正测度变换，右半平面无零点；辅助逆对数的极点不能当作完整模型极点。
- 局部可积逆对数与完整余核给短时唯一弱因果响应解；有限共同局部项需要满足明确阶数及有界条件。
- 数值核原复频分解、正测度身份、逆核Laplace身份和一个固定短时证书；没有模拟非线性时空。
- 原质量／共形块秩一及五个剪切意味着只有六个活跃非局部通道；剩余方向必须接经典方程、632共同接触和引力约束。
- 原比较背景未被宣称为自洽真空；不将局部逆替代长期稳定、ħ统一极限或实际动态背景定理。

## 本轮合并与下一项

C04／C19／C22增加一份可实际用于求解的限定响应合同，统一目标仍开放。

接[737](round737_drafts/STATUS.md)：把原共同局部项、非局部零方向和约束放回同一方程块，检验真实阶数和可逆性；再推广到730变化背景。旧空间、604、649／699保持。
'''
write('unified_physics_condition_ledger_736.md',ledger)
write('round736_drafts/research_note_736_draft.md',(HERE/'research_note_736.md').read_text('utf8'))
write('round737_drafts/STATUS.md','''# 第737轮入口：非局部零方向、共同局部项与实际约束

接[736](../research_note_736.md)、[条件账](../unified_physics_condition_ledger_736.md)。原定常完整谱的六个活跃来源通道已有保完整质量记忆的短时因果逆；不能直接反演七来源全块。

1. 回查601、630—632、731—735；631非局部质量／共形秩一和632完整接触不同，不能删去零方向或把它当规范冗余。
2. 明确原q／sigma局部主部、混合项、lapse及约束。若使用632匹配真空，必须使用那个背景的真实质量和完整标量方向，不把630诊断射线原封不动换名。
3. 检验约束消元后补足非局部零方向的局部块是否可逆，所得剩余项是否符合736连续函数空间的有界条件。存在额外导数就如实指出，不能假设该条件已自动成立。
4. 旧601共同ħ阶次作为对照保留；736短时逆不自动给ħ统一极限、物理分支或长期稳定。
5. 在同一原模型中完成这项再拓展至730实际变化背景、装置来源和非线性反馈；避免只增加各自可解但输入不一致的子模型。
6. 不重复质量谱、阈值积分或对数高频渐近，不以扫参数替代原约束审计。旧空间和统一目标保持。
''')
write('round736_drafts/literature_scope_audit.json',json.dumps(dict(sources=[
    dict(url='https://arxiv.org/html/2007.14665',locations='Proposition5.3, equations61-64; Proposition5.8 scope',
         use='Causal local inverse of a logarithmic Volterra operator; original massive fermion remainder is derived independently in note736.',
         caution='Cosmological scalar-field nonlinear existence theorem does not apply automatically to the original complete gauge/fermion/gravity model.'),
    dict(url='https://arxiv.org/abs/gr-qc/0209075',locations='Abstract and original linear-response framework',
         use='Distinguish common-action response and scale-qualified stability from merely positive spectra.',
         caution='Its explicit matter example is scalar. Fermion spectra and coefficients are reused from630/631.')],
    new='Exact complete-mass memory decomposition and local causal inversion on the active response block, with a bounded-local-term criterion and explicit null-direction boundary.',
    inherited='630/631 spectra and coefficients,632 contact conditions,601 order-reduction distinction,735 chosen local Ward prescription.',
    excluded='Nonlinear semiclassical spacetime, actual dynamic-background inverse, full constraints, global stability, or uniform classical limit.'),ensure_ascii=False,indent=2)+'\n')
write('round736_drafts/scope_and_dedup_review.md','''# 736范围与证明审查

主代理审查，无独立代理，无图像检查。前一目标轮次完成735并执行736入口，属于进展。当前未见运行中的Python；无重复启动。

- 先读入口、735和相关已存结果，回查601、630—632及734；不重复已有谱或把旧秩事实计为新发现。
- 精确分解由原全谱积分恒等得到，矩阵码只校准积分；质量阈值和全部半权不删。
- 式(3)的零点排除来自正测度符号，四个复频数值点不是无零点证明。
- 余核时间积分绝对收敛；界使用1-a≤6m²/omega²及min(omega t,1)，系数与积分原函数已核。
- 对数逆包含留数和切割；辅助极点不得当完整核增长模。完整I在右半平面无零点与辅助极点相容。
- 短时压缩证明在零过去连续u及因果分布意义成立；不暗加任意经典正则性、初态自由度或记忆重置。
- 有限局部作用必须进入同一V_T；数值纯核证书不把Einstein项默设为零。额外导数、空间无界算子或非线性态依赖尚须证明。
- ħ趋零不保证时间窗口一致；没有用短时弱解替代有效理论物理分支选择。
- 原七来源块只在六个活跃方向具非局部逆，经典互补块与引力约束仍开放；零方向不是凭空视为规范。
- 保存735全部冻结证据和736入口；旧空间与649／699范围保持。没有修改目标或定时任务。
''')
main=('research_note_736.md','joint_causal_log_inverse.py','joint_causal_log_inverse_results.json','unified_physics_condition_ledger_736.md')
write('round736_drafts/final_review.txt','Primary-agent review only. Exact old-spectrum causal inverse on active channels; no full nonlinear gravity conclusion.\n'+
      '\n'.join(n+' '+hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in main)+'\n')

mapping={'735':'737','734':'736','733':'735','3270':'3304','3284':'3318','3407':'3411','3404':'3409','1514':'1520','380':'382','285':'287',
         'joint_retarded_reference_response':'joint_causal_log_inverse',
         'joint_reference_polarization_boundary':'joint_local_source_normalization',
         'causal_reference_entry':'response_regularity_entry'}
verify=remap((HERE/'verify_round734.py').read_text('utf8'),mapping)
verify=verify.replace('(3,0,0)','(2,0,0)').replace('run=3,','run=2,').replace("==20","==18")
verify=verify.replace('    import joint_local_source_normalization as previous\n    assert previous.run()==core.read(previous.TARGET)',
    "    import sys\n    sys.path.insert(0,str(HERE/'round736_drafts'))\n    import response_regularity_entry as previous\n    assert previous.run()==core.read(previous.TARGET)")
write('verify_round736.py',verify)
publish=remap((HERE/'publish_round734.py').read_text('utf8'),mapping)
start=publish.index("summary='");end=publish.index('planned={}',start)
publish=publish[:start]+'''summary='**第736轮完成：** [原完整谱的短时因果逆]({p}research_note_736.md)原定常完整质量核精确分解为对数主部和可积记忆，活跃来源有短时唯一因果逆；有限局部项须满足明确有界合同。两组、十八式通过，最新736／3411，1520份编号科学文件、3318份保护证据。[核验]({p}research_round_736_checks.json)、[条件账]({p}unified_physics_condition_ledger_736.md)。完整约束、变化背景及非线性自洽仍开放。'
order='**当前执行顺序（736后，优先于下方历史安排）：** 接[737零方向、局部项与实际约束]({p}round737_drafts/STATUS.md)，把632共同接触和经典互补块放回七来源方程；不能仅逆六个通道就签收完整反馈。旧空间、604、649／699及统一目标保持。'
'''+publish[end:]
publish=publish.replace('同一过去准备的退迟来源与完整接触条件','原完整谱的短时因果逆与反馈适用条件')
publish=publish.replace('旧空间合同保持，退迟来源与接触须共同变分','旧空间合同保持，短时逆不替代完整约束')
write('publish_round736.py',publish)
post=remap((HERE/'postcheck_round734.py').read_text('utf8'),mapping)
write('postcheck_round736.py',post)
print('Prepared736 ledger, scope and reproducibility scripts.')
