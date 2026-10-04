"""Prepare697 publication without altering any frozen evidence."""
import hashlib
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent

def write(name,text):
    p=HERE/name;p.parent.mkdir(exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)

def replace(text,mapping):
    pats=[r'(?<!\d)'+re.escape(k)+r'(?!\d)' if k.isdecimal() else re.escape(k) for k in sorted(mapping,key=len,reverse=True)]
    return re.sub('|'.join(pats),lambda m:mapping[m.group()],text)

names=('unified_physics_condition_ledger_697.md','round697_drafts/research_note_697_draft.md',
    'round697_drafts/final_review.txt','round698_drafts/STATUS.md',
    'round697_drafts/center_holonomy_entry.py','round697_drafts/center_holonomy_entry_results.json',
    'round697_drafts/center_holonomy_entry.md','round697_drafts/check_and_publish_entry.py',
    'round697_drafts/primary_source_audit.json','round697_drafts/entry_checks.json',
    'round697_drafts/prepare_entry_publication.py','round697_drafts/full_average_probe.py',
    'round697_drafts/importance_full_average.py','round697_drafts/importance_full_average_results.json',
    'round697_drafts/joint_full_holonomy_reduction_first.py')
protected=2720+3+len(names)
assert protected==2738
ledger=(HERE/'unified_physics_condition_ledger_696.md').read_text('utf8')
ledger=ledger.replace('# 联合条件总账：696完整辅助平均后的固定接缝边界核',
                      '# 联合条件总账：697完整规范平均的统一谱隙与数值候选')
ledger=ledger.replace('接[695全账](unified_physics_condition_ledger_695.md)，回填[696报告](research_note_696.md)。[结果](joint_dynamic_auxiliary_integral_results.json)、[核验](research_round_696_checks.json)。',
    '接[696全账](unified_physics_condition_ledger_696.md)，回填[697报告](research_note_697.md)。[结果](joint_full_holonomy_reduction_results.json)、[核验](research_round_697_checks.json)。')
at=ledger.index('## 本轮合并与下一项')
ledger=ledger[:at]+'''## 697全原群上的特殊族谱隙与原来源降维

- 仍为一空间自环、两AP时间片、全16通道、零物理质量。两个中心半区的完整规范闭合消去一份Haar，保留最终holonomy；只在完整辅助平均后使用类函数Weyl公式。原商群规范化保留。
- 实际2×2极分解块对全部群元有谱隙≥√2−1。此族的零谱集合为空，不再需要693的一般坏谱区域控制；不外推到694的一般背景。
- 原256维Pfaffian精确化为32维带相位行列式，显式谱帧的方向因子完整抵消。没有取绝对值或替换原物理来源。
- 完整群／两S9测度的数值估计对向量(1,−1/1000)约为−5.03×10⁻¹⁷，估计标准误约1.05×10⁻¹⁸。独立先导固定控制变量，旧653／688自由积分只作已知均值。
- 这不是确定性积分符号证书，也没有认证置信覆盖率；完整物理RP、有限τ、原H_F身份依旧开放。675仅在严格Bᴳ负号取得后适用。
- C07—C09直接继承697 §1.1：384额外Lipschitz已经消去；386和425替代；522—523共同实现已完成，实际h、s、CAR／Gauss映射仍待接合。E不是s记录。

## 本轮合并与下一项

C01、C16、C19的完整规范平均获得无近零谱的精确计算表示及全测度负候选。所补条件是当前特殊族的统一谱隙／原来源对象身份，不是空间、连续极限或参考仪器的新假设。

接[698](round698_drafts/STATUS.md)：利用该谱隙建立有限多项式与全平均误差，核实际可计算性和符号余量。不能把有限多项式的存在当成已经求积分，也不重复无认证抽样。原指定归一、物理RP、H_F同一性、实际空间传播、共同连续及量子GR仍开放；目标与四分支不改。
'''
write('unified_physics_condition_ledger_697.md',ledger)
mapping={'joint_dynamic_auxiliary_integral':'joint_full_holonomy_reduction',
    '695':'696','696':'697','697':'698','3297':'3299','3299':'3301',
    '1397':'1400','1400':'1403','2700':'2720','2720':str(protected)}
verify=replace((HERE/'verify_round696.py').read_text('utf8'),mapping)
verify=verify.replace('import joint_auxiliary_sign_certificate as prior_model','import joint_dynamic_auxiliary_integral as prior_model')
verify=verify.replace('full auxiliary integrals and strict fixed-seam negative kernel.','uniform-gap original source reduction and full-measure numerical candidate.')
a=verify.index('    names=');z=verify.index('    preserved=',a)
verify=verify[:a]+'    names='+repr(names)+'\n'+verify[z:]
write('verify_round697.py',verify)
pub=replace((HERE/'publish_round696.py').read_text('utf8'),mapping|{'## 342.':'## 343.','## 247.':'## 248.'})
a=pub.index('summary=');z=pub.index('planned={}',a)
pub=pub[:a]+'''summary=('**第697轮完成：** [完整规范平均的统一谱隙与数值候选]({p}research_note_697.md)'
         '原中心历史在全群上谱隙≥√2−1，256维带相位来源精确降为32维。'
         '完整测度出现负核数值候选，严格积分符号尚待认证。'
         '两组、十六式通过，最新697／3301，1403份编号科学文件、2738份保护证据。'
         '[核验]({p}research_round_697_checks.json)、[全条件账]({p}unified_physics_condition_ledger_697.md)。'
         '不构成原物理RP或H_F反例；旧空间接口直接复用。')
order=('**当前执行顺序（697后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[698全积分的有限多项式与严格余项]({p}round698_drafts/STATUS.md)，'
       '核实际全平均及认证误差，不以样本标准误代替证明；'
       '保留原H_b及旧空间接口，目标不改。')
'''+pub[z:]
pub=pub.replace('完整辅助平均后的固定接缝负核','完整规范平均的统一谱隙与数值候选')
pub=pub.replace('完整辅助平均与尚待完成的规范闭合','完整规范平均的解析约化与数值边界')
pub=pub.replace('旧空间合同复用，单空间自环不等于空间生成','旧空间合同复用，全群谱隙不等于连续空间')
write('publish_round697.py',pub)
write('postcheck_round697.py',replace((HERE/'postcheck_round696.py').read_text('utf8'),mapping))
with (HERE/'round697_drafts/research_note_697_draft.md').open('xb') as f:f.write((HERE/'research_note_697.md').read_bytes())
review='''697 primary-agent proof/code/scope review; no independent agent.
Read project navigation,696 proof/results,697 entry and user-designated382-386/425/522-523.
No concurrent Python research at start; one sequential replay session now completed.
All16 channels and original N256 retained. One spatial self-loop is not intersite propagation.
Weyl theorem1 read from Woit primary university notes. Only E-averaged source is a class function.
Original quotient Haar normalization12 checked; no spurious factor6.
Actual gamma4 blocks A and A†; singular values give whole-group gap sqrt2-1, no sampling proof.
General693/694 zero-spectrum issues remain outside this restricted family.
Canonical frames u=[I;-V†]/sqrt2,v=[I;V†]/sqrt2 account for full phase.
Pf auxiliary=det S; frame determinant=(det V†)^2 cancels physical determinant phase exactly.
Original256 Pf, canonical frame factor and32 reduction crosschecked for12 backgrounds.
First discriminant-rounding failure at exact degenerate free block saved; exact modulus1 identity used.
Full torus Weyl and both S9 measures sampled with strictly positive normalized proposals.
Ordinary importance estimates, no clipping, absolute determinant or self-normalization.
Known B00 inherited653/688; independent16384 pilot,524288 main, deterministic seeds and covariance.
47.7 empirical SE is NOT rigorous coverage or a deterministic sign certificate; flag remains false.
No full physical RP/normalization/HF identity claim.675 future implication remains conditional.
Space382/383/384/386/425/522/523 inherited explicitly; no restored Lipschitz, no doubled coordinate routes.
E auxiliary differs from actual s record; no new spatial assumption supplied by this computation.
Two new verification groups,16 numbered equations,no images; next698 controlled polynomial integral.
'''
for name in ('research_note_697.md','joint_full_holonomy_reduction.py','joint_full_holonomy_reduction_results.json','unified_physics_condition_ledger_697.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round697_drafts/final_review.txt',review)
print('697 prepared; protected evidence',protected)
