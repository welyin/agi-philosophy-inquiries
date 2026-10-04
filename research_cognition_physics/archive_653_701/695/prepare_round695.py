"""Prepare695 publication; retain all frozen historical and entry evidence."""
import hashlib
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent

def write(name,text):
    p=HERE/name;p.parent.mkdir(exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)

def replace(text,mapping):
    keys=sorted(mapping,key=len,reverse=True)
    pats=[r'(?<!\d)'+re.escape(k)+r'(?!\d)' if k.isdecimal() else re.escape(k) for k in keys]
    return re.sub('|'.join(pats),lambda m:mapping[m.group()],text)

ledger=(HERE/'unified_physics_condition_ledger_694.md').read_text('utf8')
ledger=ledger.replace('# 联合条件总账：694原物理来源临界跳变的严格证书',
                      '# 联合条件总账：695原辅助逐配置非负性的严格反例')
ledger=ledger.replace('接[693全账](unified_physics_condition_ledger_693.md)，回填[694报告](research_note_694.md)。[结果](joint_critical_source_certificate_results.json)、[核验](research_round_694_checks.json)。',
    '接[694全账](unified_physics_condition_ledger_694.md)，回填[695报告](research_note_695.md)。[结果](joint_auxiliary_sign_certificate_results.json)、[核验](research_round_695_checks.json)。')
at=ledger.index('## 本轮合并与下一项')
ledger=ledger[:at]+'''## 695原完整群背景中辅助权重确实有两种符号

- 同一原2空间点×2AP时间片、两传播方向模型，保留SU3／SU2／U1及全部16通道。保存的整数定义精确有理规范链路与两份单位S9配置，不依赖随机种子重新抽样。
- 自旋互补把原辅助Pfaffian压成同一64维行列式。六个实际物种块的有理谱残差认证真实共同谱帧；Gaussian整数Bareiss计算及严格相位余量认证两份行列式之比为负。
- 另认证原物理手征块可逆及原512维N的反酉实性，故两份实际标量权重均为非零实数且符号相反。至少一份严格为负；浮点小权重不承担符号证明。
- 原候选的普遍逐配置非负性被反例排除。原S9／双Haar／H_b完整平均及实际物理RP仍未判定，不能以固定E负值取代实际物理二次型。
- 该完整规范背景不同于694特定U1路径；694固定E跳变的辅助平均并未被决定。693完整平均误差及657自由Gram、671静态正性均按原范围保留。
- C07—C09及C20逐项复用382／383／384／386／425／522／523；不恢复已消去的额外Lipschitz，不把386和425叠加，不把辅助E认作实际s记录。本轮没有补齐新的空间合同。

## 本轮合并与下一项

C01、C16及C19现在有同一原来源的严格符号边界：共同实相位不足以把辅助处方变成逐配置正测度。停止随机符号搜寻与精度扫描；此普遍断言已经有反例。

接[696](round696_drafts/STATUS.md)：复用657／675原局部球面张量，核完整辅助积分的实际外代数表示及规范／时间拼接。任何正算符表示都须证明等于原来源；不以另一个正核或辅助过程替代原物理时间。动态RP、指定非零归一、原H_F、共同连续及量子GR继续开放。目标不改。
'''
write('unified_physics_condition_ledger_695.md',ledger)
mapping={'joint_critical_source_certificate':'joint_auxiliary_sign_certificate',
         '693':'694','694':'695','695':'696','3293':'3295','3295':'3297',
         '1391':'1394','1394':'1397','2668':'2681','2681':'2700'}
verify=replace((HERE/'verify_round694.py').read_text('utf8'),mapping)
verify=verify.replace('import joint_gauge_sublevel_control as prior_model',
                      'import joint_critical_source_certificate as prior_model')
verify=verify.replace('Verify695 exact original-source one-sided discontinuity certificate.',
                      'Verify695 exact opposite-sign original auxiliary configurations.')
a=verify.index('    names=');b=verify.index('    preserved=',a)
names=('unified_physics_condition_ledger_695.md','round695_drafts/research_note_695_draft.md',
       'round695_drafts/final_review.txt','round696_drafts/STATUS.md',
       'round695_drafts/self_complement_entry.py','round695_drafts/self_complement_entry_results.json',
       'round695_drafts/self_complement_entry.md','round695_drafts/check_and_publish_entry.py',
       'round695_drafts/entry_checks.json','round695_drafts/primary_source_audit.json',
       'round695_drafts/prepare_entry_publication.py','round695_drafts/compression_sign_probe.py',
       'round695_drafts/compression_sign_probe_results.json','round695_drafts/make_negative_fixture.py',
       'round695_drafts/negative_auxiliary_fixture.json','round695_drafts/sign_primary_source_audit.json')
verify=verify[:a]+'    names='+repr(names)+'\n'+verify[b:]
verify=verify.replace("text['display_formulas']==14","text['display_formulas']==16")
write('verify_round695.py',verify)
pub=replace((HERE/'publish_round694.py').read_text('utf8'),mapping|{'## 340.':'## 341.','## 245.':'## 246.'})
a=pub.index('summary=');b=pub.index('planned={}',a)
pub=pub[:a]+"""summary=('**第695轮完成：** [原辅助逐配置非负性的严格反例]({p}research_note_695.md)'
         '同一原完整规范背景中，两份合法辅助配置的原标量权重严格非零且符号相反，'
         '关闭普遍逐配置非负路线。两组、十六式通过，最新695／3297，'
         '1397份编号科学文件、2700份保护证据。'
         '[核验]({p}research_round_695_checks.json)、[全条件账]({p}unified_physics_condition_ledger_695.md)。'
         '完整辅助／规范平均、原物理RP、H_F与共同连续仍待证。')
order=('**当前执行顺序（695后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[696完整辅助积分与物理二次型]({p}round696_drafts/STATUS.md)，'
       '复用657／675局部张量，核实际算符身份与时间拼接；停止符号随机搜寻，'
       '旧空间接口逐项复用，目标不改。')
"""+pub[b:]
pub=pub.replace('原物理来源临界谱跳变的严格证书','原辅助逐配置非负性的严格反例')
pub=pub.replace('原物理来源临界支撑与精确证书','原辅助符号与完整物理平均的边界')
pub=pub.replace('旧空间合同复用，规范路径不是空间坐标','旧空间合同复用，辅助E不是实际方向仪器')
write('publish_round695.py',pub)
write('postcheck_round695.py',replace((HERE/'postcheck_round694.py').read_text('utf8'),mapping))
with (HERE/'round695_drafts/research_note_695_draft.md').open('xb') as f:
    f.write((HERE/'research_note_695.md').read_bytes())
review='''695 primary-agent proof/code/scope review; no independent agent.
Previous goal turn completed694 and verified/published695 entry: progress.
Read current navigation,694 note/results and695 saved entry. Get-Process found no Python.
Reviewed original673 rectangular physical factor and675 exact auxiliary factorization.
Hayata-Yamamoto sufficient Majorana positivity hypotheses checked, not automatically imported.
Random search only proposes input. Stored integers define exact SU3/SU2/U1 links and S9 fields.
No added species, abs Pfaffian, squared determinant, external reference or hard admissibility cutoff.
Exact rational group unitarity/determinants and Hermitian H sectors checked.
Floating eigenvectors/inverses propose dyadic witnesses only. Six full matter modules retained.
Polarization and projector intertwiner bound simultaneous true positive and negative frames.
Same actual64 compressed determinant for each E; scalar ratio cancels E-independent physical factor.
Integer Bareiss exact divisions and negative-real cone inequalities decide relative phase.
Relative determinant phase error<0.001788353, while both determinants are nonsingular.
Actual physical prefactor separately certified by the two square chiral projection blocks.
Original512-dimensional scalar Pfaffian real by actual antiunitary identity, not physical time reflection.
Opposite nonzero real weights imply a strict negative exists; individual float signs not certified separately.
Direct original512 matrix crosscheck agrees; tiny absolute weights not used as strict proof.
Universal pointwise positivity refuted only for this original finite two-direction candidate.
Full S9/Haar/Hb average, physical RP and694 specific flux-path auxiliary average remain undecided.
657 free Gram and671 static positivity retained. Auxiliary E not physical s/instrument.
382/383/384/386/425/522/523 inherited explicitly; no restored Lipschitz or stacked geometry routes.
Fixed malformed TeX rm escape before freezing. Two new groups/16 display equations, no images.
Next696 tests original full auxiliary tensor representation, not new random sign searches.
'''
for name in ('research_note_695.md','joint_auxiliary_sign_certificate.py','joint_auxiliary_sign_certificate_results.json','unified_physics_condition_ledger_695.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round695_drafts/final_review.txt',review)
print('695 publication preparation completed')
