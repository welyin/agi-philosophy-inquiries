"""Prepare694 publication without changing frozen693 or694-entry evidence."""
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

ledger=(HERE/'unified_physics_condition_ledger_693.md').read_text('utf8')
ledger=ledger.replace('# 联合条件总账：693完整规范平均的近零谱与来源误差',
                      '# 联合条件总账：694原物理来源临界跳变的严格证书')
ledger=ledger.replace('接[692全账](unified_physics_condition_ledger_692.md)，回填[693报告](research_note_693.md)。[结果](joint_gauge_sublevel_control_results.json)、[核验](research_round_693_checks.json)。',
    '接[693全账](unified_physics_condition_ledger_693.md)，回填[694报告](research_note_694.md)。[结果](joint_critical_source_certificate_results.json)、[核验](research_round_694_checks.json)。')
at=ledger.index('## 本轮合并与下一项')
ledger=ledger[:at]+'''## 694同一原来源的支撑变化不能普遍平滑消去

- 在原2空间点×2AP时间片、两传播方向盒内，原U1路径可精确有理化。十二个原荷／S块的整数惯性认证负谱分区从64／64变为65／63，全部16通道保持。
- 同一四维谱子空间可在临界区间解析延拓；这里的固定切口不是物理sign(H)的零切口，二者明确区分。候选浮点本征向量只生成dyadic证书，结论由精确残差核验。
- 固定合法E=e₀、零质量下，原512维N的八个精确64维块均有整数逆残差证书。延拓后最小奇异值在整个区间严格大于0.136441001。
- 结合676的原不平衡零支撑，得到原标量来源的非零左极限与右侧恒零；不需要临界根唯一或额外横截性。排除“每份物理来源必连续并自动消去近零谱贡献”的普遍接法。
- 固定E现象在辅助开邻域中存在，但尚未排除S9平均的复相位抵消。B来源跳变仅有入口数值，本轮不由N可逆自动签收。
- 693原全Haar误差及677／686全平均收敛仍成立；固定来源跳变不推出规范平均失败、RP负性或对认知原则的反例。C07—C09逐项继承上表，无新增空间缺口。

## 本轮合并与下一项

C01来源、C16测度支撑、C19未归一泛函获得同一原矩阵上的严格边界：不能以物理来源的普遍连续性替代693坏谱域处理。尚缺完整辅助／规范平均的符号、动态RP、H_F及共同尺度／量子几何。

接[695](round695_drafts/STATUS.md)：先回查原自旋配对、复共轭与固定背景球面多项式，核辅助平均是否保留已认证的支撑边界。692九点全平均规则不可移作固定规范背景恒等式；停止根精度扫描。目标不改。
'''
write('unified_physics_condition_ledger_694.md',ledger)
mapping={'joint_gauge_sublevel_control':'joint_critical_source_certificate',
         '692':'693','693':'694','694':'695','3291':'3293','3293':'3295',
         '1388':'1391','1391':'1394','2655':'2668','2668':'2681'}
verify=replace((HERE/'verify_round693.py').read_text('utf8'),mapping)
verify=verify.replace('import joint_auxiliary_orbit_quadrature as prior_model',
                      'import joint_gauge_sublevel_control as prior_model')
verify=verify.replace('Verify694 original Haar sublevel bound and full-source error control.',
                      'Verify694 exact original-source one-sided discontinuity certificate.')
a=verify.index('    names=');b=verify.index('    preserved=',a)
names=('unified_physics_condition_ledger_694.md','round694_drafts/research_note_694_draft.md',
       'round694_drafts/final_review.txt','round695_drafts/STATUS.md',
       'round694_drafts/critical_source_entry.py','round694_drafts/critical_source_entry_results.json',
       'round694_drafts/critical_source_entry.md','round694_drafts/check_and_publish_entry.py',
       'round694_drafts/entry_checks.json','round694_drafts/prepare_entry_publication.py')
verify=verify[:a]+'    names='+repr(names)+'\n'+verify[b:]
verify=verify.replace("text['display_formulas']==16","text['display_formulas']==14")
write('verify_round694.py',verify)
pub=replace((HERE/'publish_round693.py').read_text('utf8'),mapping|{'## 339.':'## 340.','## 244.':'## 245.'})
a=pub.index('summary=');b=pub.index('planned={}',a)
pub=pub[:a]+"""summary=('**第694轮完成：** [原物理来源临界谱跳变的严格证书]({p}research_note_694.md)'
         '整数惯性与原512维N的逆残差证明固定合法E的标量来源存在非零单侧跳变，'
         '排除来源自动平滑消去近零谱的普遍接法。两组、十四式通过，最新694／3295，'
         '1394份编号科学文件、2681份保护证据。'
         '[核验]({p}research_round_694_checks.json)、[全条件账]({p}unified_physics_condition_ledger_694.md)。'
         '尚非完整辅助／规范平均结论，不构成RP、H_F或连续时空证明。')
order=('**当前执行顺序（694后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[695完整辅助平均与临界支撑]({p}round695_drafts/STATUS.md)，'
       '核原自旋配对及相位，不能从固定E非零直接推积分非零；旧空间接口逐项复用，目标不改。')
"""+pub[b:]
pub=pub.replace('完整规范平均的近零谱与原来源误差','原物理来源临界谱跳变的严格证书')
pub=pub.replace('完整规范平均的近零谱与来源控制','原物理来源临界支撑与精确证书')
pub=pub.replace('旧空间合同复用，有限谱模量不是空间尺度','旧空间合同复用，规范路径不是空间坐标')
write('publish_round694.py',pub)
write('postcheck_round694.py',replace((HERE/'postcheck_round693.py').read_text('utf8'),mapping))
with (HERE/'round694_drafts/research_note_694_draft.md').open('xb') as f:
    f.write((HERE/'research_note_694.md').read_bytes())
review='''694 primary-agent proof/code/scope review; no independent agent.
Previous goal turn completed693 and verified/published694 entry: progress.
Read current navigation,693 note/results and694 saved entry. Get-Process found no Python.
Original676 path reparameterized rationally, not a new physical cutoff or group.
Twelve actual614 charge/S blocks use exact integer Laurent coefficients.
Exact realified rational symmetric-congruence inertia:64/64 to65/63.
Floating eigenvectors/inverses propose dyadic data only; all decisive norms checked with integer matrices/Fractions.
Polar-decomposition error accounts for nonunitary proposed frames; exact residual controls eigenvalue gaps.
Fixed rank-four cuts separated throughout rational interval. These cuts are not the sign(H) zero cut.
Resolvent identity bounds projector perturbation; original Wilson variation<=6|q| interval length.
Original signed J and normalized S frames used; direct-sum maximum error, not omitted multiplicities.
Original512-dimensional N uses all16 channels and E=e0,lambda0. Exact sparsity gives eight64 blocks.
Integer inverse-residual and norm certificate bounds minimum singular value; no inverse defines sources.
Strict continuation margin>0.136441001, implies scalar Pfaffian limit nonzero.
Analytic isolated roots finite on compact interval; endpoint count change guarantees a balance-to-imbalance boundary.
676 supplies exact right-open zero source; no transverse or unique-root premise added.
B source nonzero limit not proved. FixedE result and auxiliary-open-neighborhood extension not confused with S9 average.
693 full-Haar weightedL1 approximation retained. No RP/refutation/HF/continuum/GR claim.
382/383/384/386/425/522/523 inherited: removed Lipschitz not restored;386/425 alternative;E not s.
Two new groups/14 display equations; no images; next695 actual auxiliary phase/average, no precision-only scanning.
'''
for name in ('research_note_694.md','joint_critical_source_certificate.py','joint_critical_source_certificate_results.json','unified_physics_condition_ledger_694.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round694_drafts/final_review.txt',review)
print('694 publication preparation completed')
