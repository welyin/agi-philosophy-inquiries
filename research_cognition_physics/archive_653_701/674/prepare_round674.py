"""Prepare674, preserving all published evidence and the unchanged goal."""
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

ledger=(HERE/'unified_physics_condition_ledger_673.md').read_text('utf8')
ledger=ledger.replace('# 联合条件总账：673共同边界、变秩字典与有限动态泛函',
    '# 联合条件总账：674互补配对与完整标量泛函实性')
ledger=ledger.replace('接[672全账](unified_physics_condition_ledger_672.md)，回填[673报告](research_note_673.md)。[结果](joint_gauss_boundary_functional_results.json)、[核验](research_round_673_checks.json)。',
    '接[673全账](unified_physics_condition_ledger_673.md)，回填[674报告](research_note_674.md)。[结果](joint_gauss_reflection_reality_results.json)、[核验](research_round_674_checks.json)。')
ledger=ledger.replace('673原H_b与全边界的有限未归一泛函|',
    '673原H_b与全边界的有限未归一泛函；674完整标量核Hermitian与总权重实性|')
ledger=ledger.replace('673物理来源矩形字典保同一观测；原物理态',
    '673物理来源矩形字典保同一观测；674原标量／质量来源反射已证；原物理态')
ledger=ledger.replace('673全双边界Haar对象及全部字段输送已接；群选择',
    '673全双边界Haar对象及全部字段输送已接，674标量核Hermitian已证；群选择')
ledger=ledger.replace('一般反射身份、动态正性、指定归一与连续局域性仍缺',
    '674标量／质量来源反射身份已证；任意费米来源反射、动态正性、指定归一与连续局域性仍缺')
ledger=ledger.replace('673全配置质量／物理来源绝对可积；动态过程',
    '673全配置质量／物理来源绝对可积；674原质量反射含奇异质量的多项式延伸；动态过程')
ledger=ledger.replace('动态归一多项式的指定零点仍待排除；原H_F',
    '674总权重及反射配对实质量来源为实；总权重正性和指定零点仍待排除；原H_F')
ledger=ledger.replace('670非平坦来源及零模处正则求导|',
    '670非平坦来源及零模处正则求导；674有限实质量来源的共轭合同|')
at=ledger.index('## 本轮合并与下一项')
ledger=ledger[:at]+'''- 674采用成熟互补Pfaffian恒等式，证明原标量与局部实质量来源的反射身份，包含不同手征秩、零权重及奇异质量。临时质量逆仅在证明的稠密集使用；最终固定尺寸定义不含逆质量或逆传播子。一般Grassmann来源合同仍须另外证明。

## 本轮合并与下一项

674连接C01、C14、C16、C17、C19、C22：原辅助配对、原物理质量、规范时间反射与完整双Haar边界共同给标量核Hermitian；与673绝对可积性结合，完整总权重及反射配对的实质量来源为实。这是解析补足673撤回的一部分表述，未把旧失败相位诊断解释成反例。

辅助Pfaffian与物理质量的互补子式共同保留原方向、相位和变秩结构。近零样本不再用不可靠相对误差作判据；失败及首次通过版本均保留。没有取权重绝对值、删负方向、人工正化、降到低维子群或数值替代全群积分。

实性不等于非负、非零、全费米反射正性或原H_F身份。完整正半区代数、原有序CAR影响／真实记录、共同时间参数及连续／量子GR仍开放。辅助E是否为物理观测仍须遵守原合同，不能用固定E上的失败直接否定E积分后的物理泛函。

接[675](round675_drafts/STATUS.md)：检验完整平均后的正性，先核原H_b基态投影与双Haar的受控数学约化及固定λ范围，再接原有序过程；不把该候选参数极限当作已识别的原H_F低温过程。旧空间成果沿667 §1.1直接继承，认知设计后置，统一目标不改。
'''
write('unified_physics_condition_ledger_674.md',ledger)
mapping={'joint_gauss_boundary_functional':'joint_gauss_reflection_reality',
         '672':'673','673':'674','674':'675','3251':'3253','3253':'3255',
         '1328':'1331','1331':'1334','2403':'2419','2419':'2433'}
verify=replace((HERE/'verify_round673.py').read_text('utf8'),mapping)
verify=verify.replace('Verify rectangular physical sources and full double-boundary finite functional.',
    'Verify scalar reflection via complementary Pfaffians and full Gauss reality.')
verify=verify.replace('import joint_gauge_history_kernel as prior_model','import joint_gauss_boundary_functional as prior_model')
verify=verify.replace('boundary_holonomy_probe','unnormalized_source_probe')
a=verify.index('    names=');b=verify.index('    preserved=',a)
verify=verify[:a]+"""    names=('unified_physics_condition_ledger_674.md','round674_drafts/research_note_674_draft.md',
           'round674_drafts/final_review.txt','round675_drafts/STATUS.md',
           'round674_drafts/unnormalized_source_entry.md','round674_drafts/unnormalized_source_probe.py',
           'round674_drafts/unnormalized_source_probe_results.json',
           'round674_drafts/first_reflection_attempt.py','round674_drafts/first_reflection_diagnostic.json',
           'round674_drafts/first_pass_reflection.py','round674_drafts/first_pass_reflection_results.json')
"""+verify[b:]
verify=verify.replace("assert text['display_formulas']==18","assert text['display_formulas']==16")
write('verify_round674.py',verify)
publish=replace((HERE/'publish_round673.py').read_text('utf8'),mapping|{'## 319.':'## 320.','## 224.':'## 225.'})
a=publish.index('summary=');b=publish.index('planned={}',a)
publish=publish[:a]+'''summary=('**第674轮完成：** [互补手征配对、完整规范边界与标量泛函的实性]({p}research_note_674.md)'
         '互补Pfaffian证明原标量／实质量来源的反射身份，含变秩和奇异质量；'
         '接673可积性，完整双Gauss标量核Hermitian、总权重为实。'
         '正性、非零归一、任意费米来源反射与原过程身份仍待证明。'
         '两组、十六式通过，最新674／3255，1334份编号科学文件、2433份保护证据。'
         '[核验]({p}research_round_674_checks.json)、[全条件账]({p}unified_physics_condition_ledger_674.md)。'
         '连续及量子GR仍开放。')
order=('**当前执行顺序（674后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[675完整标量边界核的正性与共同时间]({p}round675_drafts/STATUS.md)，'
       '核全群平均、原H_b与实际时间参数；旧空间接口直接复用，目标不改。')
'''+publish[b:]
publish=publish.replace('共同规范边界、变秩手征字典与有限动态泛函','互补手征配对、完整规范边界与标量泛函的实性')
publish=publish.replace('原物理来源与完整边界积分的共同存在条件','原配对、完整边界与标量反射的共同合同')
publish=publish.replace('有限动态泛函与正物理过程的验收区别','完整标量泛函的实性与正性的区别')
write('publish_round674.py',publish)
write('postcheck_round674.py',replace((HERE/'postcheck_round673.py').read_text('utf8'),mapping))
write('round674_drafts/research_note_674_draft.md',(HERE/'research_note_674.md').read_text('utf8'))
review='''674 primary-agent proof/code/scope review; no independent agent.
Navigation, latest reports, saved entry, code and frozen673 receipt reviewed.
Goal unchanged; no new task, agent, automation, image check or installation.
Old space interfaces382-386,425,522-523 directly inherited through667.
No old existence result re-opened; only same-object embedding remains.
Proof: T swaps Jplus/Jminus with blocks +I,-I; verified original fixed matrices.
For even ru,rv and r multiple32, swapping frames contributes no sign.
Auxiliary skew unitary M has Pf1 and complementary Au=d*conj(Av).
Pfaffian Schur on physical q is temporary: no inverse Wilson, Kl or weight.
Complement identity holds at zero compressed Pfaffians, by polynomial minors.
q=-conj(p) and r/2 even fix Pf mass prefactor; original frame determinants1.
Odd auxiliary dimensions give zero scalar and pure physical-source coefficients.
Invertible skew p are dense; fixed original scalar identity extends to all p,
including original singular Higgs and zero masses. Local real mass sources
must be reflected together; arbitrary Grassmann-source reflection not proved.
Each temporal seam inverted, not swapped; actual links checked at four
full-group samples with all original mass channels, phi and E retained.
Full Haar and heat symmetry plus673 integrability imply Hermitian scalar
orbit kernel and real fully integrated weight. They do not imply RP or Z>0.
Fixed auxiliary E is not silently declared a physical observable.
First relative-Pfaffian attempt failed at near-zero compressions; preserved.
Balanced absolute residuals are diagnostics, not certified high accuracy;
resolved samples only get relative scalar comparison. Algebra carries proof.
Unequal-rank/empty-frame fixtures not claimed physical Wilson topologies.
Mass/Higgs test singular rank4 uses inverse-free original finite N.
2 check groups,16 equations. No Haar integral evaluated this round.
Next675 must inspect actual full average and time contract, not thermal tuning.
'''
for name in ('research_note_674.md','joint_gauss_reflection_reality.py','joint_gauss_reflection_reality_results.json','unified_physics_condition_ledger_674.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round674_drafts/final_review.txt',review)
print('674 preparation completed; scope audited')

