"""Prepare705 publication; preserve all frozen history including entry attribution."""
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

names=('unified_physics_condition_ledger_705.md','round705_drafts/research_note_705_draft.md',
       'round705_drafts/final_review.txt','round706_drafts/STATUS.md',
       'round705_drafts/regional_compression_entry.py','round705_drafts/regional_compression_entry_results.json',
       'round705_drafts/regional_compression_entry.md','round705_drafts/check_and_publish_entry.py',
       'round705_drafts/entry_checks.json','round705_drafts/prepare_entry_publication.py',
       'round705_drafts/literature_scope_audit.json','round705_drafts/attribution_correction_363.md')
protected=2855+3+len(names);assert protected==2870
ledger=(HERE/'unified_physics_condition_ledger_704.md').read_text('utf8')
ledger=ledger.replace('# 联合条件总账：704准备几何与真实记录历史的共同响应',
                      '# 联合条件总账：705区域压缩、边界电荷与Gauss支持')
ledger=ledger.replace('接[703全账](unified_physics_condition_ledger_703.md)，回填[704报告](research_note_704.md)。[结果](joint_preparation_history_limit_results.json)、[核验](research_round_704_checks.json)。',
    '接[704全账](unified_physics_condition_ledger_704.md)，回填[705报告](research_note_705.md)。[结果](joint_regional_charge_compression_results.json)、[核验](research_round_705_checks.json)。')
ledger=ledger.replace('## 当前共同对象及仍存在的分支',
    '**705当前增量：** 原完整热态在非平凡切口有无限边界表示支持。严格单侧有限总输出与精确Gauss不相容；Gauss损失、远端改变及输出泄漏之和至少为遗漏扇区概率。保全部边界载体、压缩重数的局部CPTP族可保持Gauss和未知远端，且有原比较能源误差，但整体仍无限维。实际新动态／来源匹配和空间连续仍开放。入口压缩交换子直接归属363，补正不计轮次。\n\n## 当前共同对象及仍存在的分支')
updates={
    'C01 量子对象':'705固定原物理空间上的边界保留局部CPTP；每扇区有限，整体无限',
    'C02 区域组合':'705严格单侧通道保持完整Gauss所需的载体身份；两侧重数压缩相互交换',
    'C03 事件记录':'705未知远端及被动参考边缘保持；记录通道的时间拼接仍待核',
    'C13 熵与面积':'705同一局部保边界通道对态和参考给相对熵极限，直接复用638工具',
    'C14 规范群':'705原Z₆的无限Q=6n环表示给严格热支持；协变密度不等于Gauss支持',
    'C19 参考态':'705同一忠实全H Gibbs态使任何有限边界表示遗漏概率严格正，无须改成电热态',
    'C20 尺度映射':'705有限输出／Gauss／远端三误差权衡；保边界重数近似有参考一致能源界',
    'C21 主体存储':'705严格局部精确匹配必须保原无限边界支持，或改变有限容量／误差／访问合同',
    'C22 来源反作用':'705原K_A+K_B比较能源一二阶矩不增；完整H来源响应不能自动继承'}
rows=ledger.splitlines();seen=[]
for i,row in enumerate(rows):
    if row.startswith('|原固定图完整量子过程|'):
        cells=row.split('|');cells[2]+='；705区域局部压缩的边界支持合同及保载体替代';rows[i]='|'.join(cells)
    for key,value in updates.items():
        if row.startswith('|'+key+'|'):
            cells=row.split('|');cells[2]+='；'+value;rows[i]='|'.join(cells);seen.append(key)
assert set(seen)==set(updates)
ledger='\n'.join(rows)+'\n';at=ledger.index('## 本轮合并与下一项')
ledger=ledger[:at]+'''## 705区域、Gauss支持、热参考与容量的共同条件

- 363式(4)已有压缩交换子恒等式；705入口只作旧式热尾推论及原读口校准，未计完整轮次。新增补正见[归属记录](round705_drafts/attribution_correction_363.md)，冻结旧稿不覆盖。
- 617的同一切口表示给HA=直和(Mlambda张量Vlambda)，原物理态由两侧对偶载体匹配。保全部物理输入的严格单侧Kraus在配对部门必为直和(k_lambda张量I_V)；成熟超选结构直接复用。
- 636原字符环可取lambda_n=(0,0,0,6n)，正常且有限能源，各n都是原允许部门。603忠实全H Gibbs态在任何有限边界集合之外仍有严格正概率，未计算其数值。
- 对有限不变输出F_A，Gauss支持损失eta、远端迹距离delta及输出泄漏epsilon满足eta+delta+epsilon>=p_excluded。原单侧保迹通道远端不变，故严格有限输出无法同时完全保Gauss。全局625／704族与617精确等距不满足这组额外要求，不被否定。
- 原L/u/Q两切口载体维数为d²；只去极化A载体后物理存活为d^-4。标签未变、通道协变也不够，匹配载体必须一起保留。
- 新局部Phi_A,N按637真实K_A截断重数，每个lambda的尾回填到同扇区最低重数向量，保持整个Vlambda。单一主Kraus保留已保留扇区的相干，不先做全标签测量。
- 输出每扇区有限但整体无限；保全部Gauss、未知远端／参考、K_A一二阶矩，受能源预算的迹范数误差<=2sqrt(E_A/N)+E_A/N。两区域通道相互交换并共同保匹配。
- 同通道作用真实态和参考，638数据处理及下半连续性给局部相对熵极限。没有变成有限矩阵算法；K_A不是完整物理Hamiltonian，矩控制不冒充完整H下降或真实二阶来源。
- 数值三组保原允许标签及非Abel载体，重数能阶(0,1,7)为明确诊断输入；不求原全图热谱。无限范围、严格支持和误差量词由解析承担。
- 382—386、425、522—524均直接继承；不将内部边界表示维数认作空间维数，不恢复384已删Lipschitz，386／425仍替代。

## 本轮合并与下一项

C02／C14／C19／C20／C21的区域局部化和容量要求必须共同选择；可选保边界无限载体、跨边协调、近似支持或扩大系统边界。它们不是新的认知公理，严格有限区域容量也从未被自动加入原原则。

接[706](round706_drafts/STATUS.md)：核保边界局部族的原真实历史及来源图范数。592有界历史望远镜可作入口推论，不另计整轮；新压缩Hamiltonian、自身热准备及二阶来源须同对象验证。四分支和699边界不变，认知设计后置，目标不变。
'''
write('unified_physics_condition_ledger_705.md',ledger)
mapping={'joint_preparation_history_limit':'joint_regional_charge_compression',
         '703':'704','704':'705','705':'706','3313':'3315','3315':'3318',
         '1421':'1424','1424':'1427','2839':'2855','2855':'2870'}
verify=replace((HERE/'verify_round704.py').read_text('utf8'),mapping)
verify=verify.replace('import joint_geometry_thermal_limit as prior_model',
                      'import joint_preparation_history_limit as prior_model')
verify=verify.replace('joint preparation-Gibbs and actual dynamical record jets.',
                      'regional charge support, finite-output tradeoff and local multiplicity channels.')
a=verify.index('    names=');b=verify.index('    preserved=',a)
verify=verify[:a]+'    names='+repr(names)+'\n'+verify[b:]
verify=verify.replace("text['display_formulas']==18","text['display_formulas']==20")
verify=verify.replace("==(2,0,0)","==(3,0,0)").replace('fresh_tests=dict(run=2,','fresh_tests=dict(run=3,')
verify=verify.replace("'round687_drafts/attribution_correction_653.md')",
    "'round687_drafts/attribution_correction_653.md','round705_drafts/attribution_correction_363.md')")
write('verify_round705.py',verify)
pub=replace((HERE/'publish_round704.py').read_text('utf8'),mapping|{'## 350.':'## 351.','## 255.':'## 256.'})
a=pub.index('summary=');b=pub.index('planned={}',a)
pub=pub[:a]+'''summary=('**第705轮完成：** [区域压缩、边界电荷与Gauss支持]({p}research_note_705.md)'
         '原完整热态的无限边界支持使严格单侧有限输出与精确Gauss冲突，给出三误差权衡；'
         '保边界、压缩重数的局部通道可行，但整体仍无限维。'
         '三组、二十式通过，最新705／3318，1427份编号科学文件、2870份保护证据。'
         '[核验]({p}research_round_705_checks.json)、[全条件账]({p}unified_physics_condition_ledger_705.md)。'
         '新族动态来源及空间连续仍开放；入口交换子归属363已补正，旧空间与目标不变。')
order=('**当前执行顺序（705后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[706保边界局部近似与原真实过程]({p}round706_drafts/STATUS.md)，'
       '核比较能源到实际来源的共同域，复用592／623—625／637；'
       '保留四分支及699范围，停止有限容量反例与精度扫描。')
'''+pub[b:]
pub=pub.replace('准备几何与真实记录历史的共同响应','区域压缩、边界电荷与Gauss支持')
pub=pub.replace('准备、演化与记录来源共用一个有限过程','区域、边界支持及有限容量共同验收')
pub=pub.replace('旧空间接口继承，固定图混合响应不替代空间映射','旧空间接口继承，内部边界表示不当空间维数')
write('publish_round705.py',pub)
write('postcheck_round705.py',replace((HERE/'postcheck_round704.py').read_text('utf8'),mapping))
with (HERE/'round705_drafts/research_note_705_draft.md').open('xb') as f:
    f.write((HERE/'research_note_705.md').read_bytes())
review='''705 primary-agent mathematical/code/scope review; no independent agent.
Previous turn704 and705 entry completed and archived, constituting progress.
363 already proved compression commutator identity. Added explicit attribution correction, no new round.
360-365 shared center and finiteZ2 protocols;617 cut;636 loop labels;637 energy;638 entropy reused.
Original boundary group is full quotientG product at cuts, not an invented gauge group.
Strict-local means channel on extendedA tensor identityB under the same cut definition.
Positivity of Kraus outputs and Schur identify boundary-invariant Kraus on paired physical sectors.
Covariant CP channel alone need not preserve invariant-vector support; explicit carrier check.
Original pureQ=6n loop states give infinitely many nonzero paired sectors, full interactions retained.
Faithful603 Gibbs yields strictly positive excluded-sector probability without computing spectrum.
Finite output invariantF has finite irrep set; projection inequality proves eta+delta+epsilon>=tail.
Tradeoff applies fixed input and does not assert strict local finite outputs were cognition axioms.
Finite-tail approximation allowed; original637 Casimir energy yields upper tail bound.
Positive construction keeps all boundary labels/carriers, only truncates multiplicities.
Normal countable Kraus completeness and invariance; mainP preserves low-sector coherences.
Ground multiplicity vector exists from compact-resolvent originalKA restriction.
Output finite per sector but globally infinite; apparatus implementation not claimed.
ComparisonKA first/second moments decrease, not fullHF energy.
Gentle tail bound reference-consistent; two local maps commute and preserveGauss.
Nested local maps justify relative entropy monotonicity with same state/reference and trace limits.
No actual source derivative theorem for new family, no finite effective Hamiltonian claimed.
Actual original character/representation calibrations; multiplicity gaps0,1,7 openly declared.
Numerical middle cutoffs inaccurate but obey bounds; infinite quantifiers analytic, not sampled.
No continuum, chiral, quantumGR, group choice or physical dimension inference.
Old382-386,425,522-524 retained exactly;386/425 alternatives, removedLipschitz not restored.
Three scientific groups,20 formulas; next706 real-process/source interface.
'''
for name in ('research_note_705.md','joint_regional_charge_compression.py','joint_regional_charge_compression_results.json',
             'unified_physics_condition_ledger_705.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round705_drafts/final_review.txt',review)
print('705 prepared; protected evidence',protected)

