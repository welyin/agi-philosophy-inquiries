"""Prepare702 while preserving all frozen science, entry evidence and navigation."""
import hashlib
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent


def write(name,text):
    p=HERE/name;p.parent.mkdir(exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)


def replace(text,mapping):
    pats=[r'(?<!\d)'+re.escape(k)+r'(?!\d)' if k.isdecimal() else re.escape(k)
          for k in sorted(mapping,key=len,reverse=True)]
    return re.sub('|'.join(pats),lambda m:mapping[m.group()],text)


names=('unified_physics_condition_ledger_702.md','round702_drafts/research_note_702_draft.md',
       'round702_drafts/final_review.txt','round703_drafts/STATUS.md',
       'round702_drafts/source_scale_entry.py','round702_drafts/source_scale_entry_results.json',
       'round702_drafts/source_scale_entry.md','round702_drafts/check_and_publish_entry.py',
       'round702_drafts/entry_checks.json','round702_drafts/prepare_entry_publication.py',
       'round702_drafts/literature_scope_audit.json')
protected=2809+3+len(names);assert protected==2823
ledger=(HERE/'unified_physics_condition_ledger_701.md').read_text('utf8')
ledger=ledger.replace('# 联合条件总账：701观测分辨率与共同正性极限',
                      '# 联合条件总账：702原完整过程的共同热态、演化与记录')
ledger=ledger.replace('接[700全账](unified_physics_condition_ledger_700.md)，回填[701报告](research_note_701.md)。[结果](joint_flow_resolution_limit_results.json)、[核验](research_round_701_checks.json)。',
    '接[701全账](unified_physics_condition_ledger_701.md)，回填[702报告](research_note_702.md)。[结果](joint_gibbs_preserving_transfer_results.json)、[核验](research_round_702_checks.json)。')
ledger=ledger.replace('## 当前共同对象及仍存在的分支',
    '**702当前增量：** 原完整有限图H_F的新束缚势分割具有统一热尾控制；其自身Gauss Gibbs态、真实时间、实际记录联合后态、全部热能量矩及热熵共同收敛。补齐655明确保留的变化热态缺口。有限步处方不同于655和673，699负证书不撤销；空间细化、一般几何来源及四分支同一映射仍开放。\n\n## 当前共同对象及仍存在的分支')
ledger=ledger.replace('原固定图完整量子过程|603／623—625／643／652／655／662；665接触扩展；666—667参考、空间传播、质量及来源同步|',
    '原固定图完整量子过程|603／623—625／643／652／655／662；665接触扩展；666—667参考、空间传播、质量及来源同步；702自身热态及同一真实记录、能源与热熵联合极限|')
updates={
    'C01 量子对象':'702完整Gauss自身Gibbs态在迹范数中回到原态',
    'C03 事件记录':'702变化热态与同族真实演化下的原记录联合后态共同收敛',
    'C04 内部演化':'702正转移对数保原全部相互作用，强预解及热迹共同收敛',
    'C09 参考系统':'702近似过程自己的热参考得到归一及迹范数控制；不另借固定原态',
    'C13 熵与面积':'702全部能量矩控制热熵极限，不领取面积律',
    'C20 尺度映射':'702固定图时间近似有统一热尾及共同Gibbs／记录／能源极限；空间一致界仍缺'}
rows=ledger.splitlines();seen=[]
for i,row in enumerate(rows):
    for key,value in updates.items():
        if row.startswith('|'+key+'|'):
            cells=row.split('|');cells[2]+='；'+value;rows[i]='|'.join(cells);seen.append(key)
assert set(seen)==set(updates),seen
ledger='\n'.join(rows)+'\n';at=ledger.index('## 本轮合并与下一项')
ledger=ledger[:at]+'''## 702原完整量子过程的统一热尾与同一记录极限

- 明确新分割Aθ=K+θW，Dθ=V−θW+J+C+C†，固定0<θ<1。所有原质量、跳跃、场和规范动能保留；θ只分配原势，总H不变。603的原势尾保证Aθ为Gibbs生成元。
- 点态Young界给Dθ≥−cθ，正分割Sa≤exp(acθ)exp(−aAθ)。log的算符单调性给Ha+cθ≥Aθ；没有使用错误的指数算符单调性，cθ不是几何真空扣除。
- 在同一Gauss空间，Tr exp(−βLa)≤ZA(β)，Tr Aθ exp(−βLa)≤2ZA(β/2)/(eβ)。固定Aθ低能有限秩块与统一尾界，将655强收敛升级为热算符迹范数收敛。
- Za趋向严格正原Z，近似自身Gibbs态ρa趋ρ。原624仪器与同一Ha的真实时间共同给记录联合后态的迹范数极限；零概率后选择不领取统一相对误差。
- 复热时间的迹类局部界和Vitali／Cauchy公式给每个固定阶能源矩；利用Gibbs身份给热熵极限。固定a的能源来源不等于对a=β/n的总导数，一般几何来源尚未证明。
- 有界偶Gauss Euclidean词的正性来自同一Hilbert过程；不能由此把原673费米来源、空间overlap物种控制或699候选全部搬入。不同固定θ回到同一个H。
- 数值只校准原64CAR条件纤维下界及256维径向—中性Majorana耦合，不声称算过全图Gauss热迹。解析定理不删原相互作用；诊断中的Dirac零耦合、有限盒和差分逐项声明。
- 702入口的699无量纲负裕量仍有效，未计新轮次；其来源放大没有实际资源证书。旧382—386、425、522—524直接复用，384额外Lipschitz不恢复，386／425不叠加，E不认作s。

## 本轮合并与下一项

C01态、C03记录、C04演化、C09热参考、C13热熵和C20时间近似在原完整有限图中共享一族自身归一的正过程。原正H_F、指定连续物质、给定作用经典几何与手征辅助候选四分支仍未自动统一。

接[703](round703_drafts/STATUS.md)：先核实际有界记录／规范来源的共同响应，再明确改变动能／域的一般几何来源需要什么。旧603已有来源有限性不重证；不能把态收敛直接当全部来源导数收敛。停止优化702有限格点及Young常数，目标不变。
'''
write('unified_physics_condition_ledger_702.md',ledger)
mapping={'joint_flow_resolution_limit':'joint_gibbs_preserving_transfer',
         '700':'701','701':'702','702':'703','3307':'3309','3309':'3311',
         '1412':'1415','1415':'1418','2795':'2809','2809':'2823'}
verify=replace((HERE/'verify_round701.py').read_text('utf8'),mapping)
verify=verify.replace('import joint_anisotropic_time_contract as prior_model',
                      'import joint_flow_resolution_limit as prior_model')
verify=verify.replace('joint magnetic-source continuum/RP bound and nonuniform negative limit.',
                      'original full-model Gibbs, real-time, records and energy/entropy joint limit.')
a=verify.index('    names=');b=verify.index('    preserved=',a)
verify=verify[:a]+'    names='+repr(names)+'\n'+verify[b:]
write('verify_round702.py',verify)
pub=replace((HERE/'publish_round701.py').read_text('utf8'),mapping|{'## 347.':'## 348.','## 252.':'## 253.'})
a=pub.index('summary=');b=pub.index('planned={}',a)
pub=pub[:a]+'''summary=('**第702轮完成：** [原完整过程的共同热态、演化与记录]({p}research_note_702.md)'
         '保留原全部相互作用的新正分割有统一热尾控制；自身Gauss Gibbs态、真实时间、'
         '记录联合后态、全部能量矩及热熵共同收敛，补齐655变化热态缺口。'
         '两组、十六式通过，最新702／3311，1418份编号科学文件、2823份保护证据。'
         '[核验]({p}research_round_702_checks.json)、[全条件账]({p}unified_physics_condition_ledger_702.md)。'
         '固定图，不修复699或证明空间连续；旧空间接口与目标不变。')
order=('**当前执行顺序（702后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[703共同热过程的实际来源响应]({p}round703_drafts/STATUS.md)，'
       '区分有界物理来源和改变动能的一般几何来源，核同一分割、态及来源导数；'
       '停止有限格点优化，保留四分支及699反例范围。')
'''+pub[b:]
pub=pub.replace('观测分辨率与共同正性极限','原完整过程的共同热态、演化与记录')
pub=pub.replace('物理来源分辨率与正性需要共同取极限','原完整过程的热参考、记录和演化共用一个极限')
pub=pub.replace('旧空间及局域探针继承，来源极限不替代维数前提','旧空间接口继承，热尾紧性不替代空间细化')
write('publish_round702.py',pub)
write('postcheck_round702.py',replace((HERE/'postcheck_round701.py').read_text('utf8'),mapping))
with (HERE/'round702_drafts/research_note_702_draft.md').open('xb') as f:
    f.write((HERE/'research_note_702.md').read_bytes())
review='''702 primary-agent proof/code/scope review; no independent agent.
Prior goal turn completed701 and702 entry: progress. No live Python process found.
603 supplies confining heat trace;655 explicitly did not prove varying own Gibbs convergence.
New split changes finite-step representation, preserves original full H, all masses/hopping and Gauss.
A=K+theta W is Gibbs for fixed0<theta<1; bare K is not falsely assumed Gibbs on noncompact target.
D bounded below by original full-configuration Young estimate; no scalar cutoff.
S positive injective, not bounded away from zero; log defined by spectral calculus.
Log monotonicity understood via epsilon resolvent integrals and closed-form limit.
Ha+c>=A and minmax yield trace bounds; no false exponential operator monotonicity.
Energy-weighted trace bound evaluates A<=La on La eigenvectors; noncommutation allowed.
Finite-rank A projection used only to prove tightness, not alter physical model.
Strong convergence plus tightness gives trace norm, then positive original Z gives normalization.
Real-time histories use original bounded instruments and same Ha/Gibbs; zero-probability branches scoped.
Passive reference needs joint-state convergence, not just convergent reduced states.
Banach-valued holomorphy, right-half-plane local trace bounds and real-axis convergence give Vitali.
Energy derivatives fixeda; Gibbs entropy follows energy plus logZ, not arbitrary trace-distance continuity.
Even Gauss RP from Hilbert trace does not sign original chiral odd-source dictionary.
Constants fixed graph; no spatial uniformity, thermalization, dimension or GR proof.
Numerical256 calibration declares radial box, original neutral Majorana and Diraczero choices.
Full64 fiber uses original particle-hole/BdG trace convention; general lower bound analytic.
699 remains fixed-candidate failure; four branches distinct, no claim new transfer repairs it.
384 Lipschitz removed;386/425 alternatives;523/524 realized probes inherited; E not s.
Two check groups and16 equations;703 will inspect actual physical-source response, not tune grids.
'''
for name in ('research_note_702.md','joint_gibbs_preserving_transfer.py','joint_gibbs_preserving_transfer_results.json',
             'unified_physics_condition_ledger_702.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round702_drafts/final_review.txt',review)
print('702 prepared; protected evidence',protected)
