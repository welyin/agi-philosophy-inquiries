"""Prepare760 publication without modifying frozen evidence."""
import hashlib
import json
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent

def write(name,text):
    p=HERE/name;p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)

def remap(text,mapping):
    return re.sub('|'.join(re.escape(k) for k in sorted(mapping,key=len,reverse=True)),lambda m:mapping[m[0]],text)

ledger=(HERE/'unified_physics_condition_ledger_759.md').read_text('utf8')
ledger=ledger.replace(ledger.split('\n')[0],'# 联合条件总账：760混合Gauss过程与条件来源流',1)
ledger=ledger.replace(ledger.split('\n')[2],
 '2026-10-04。接[759全账](unified_physics_condition_ledger_759.md)，回填[760报告](research_note_760.md)。[结果](joint_conditional_background_results.json)、[核验](research_round_760_checks.json)。统一目标保持。',1)
updates={
 'C01':'760以有限列振幅表示原有限秩物理态；列不是免费外部系统，构形间相干保留',
 'C03':'760原sin s记录及任意有限等待历史由同一条件振幅精确保留，终端自主实现仍缺',
 'C04':'760精确条件动态无需能带隙；继续保量子慢变量，尚非压缩的宏观动力学',
 'C19':'760背景密度与条件振幅是原物质/规范配置的表示，外给空间度规仍固定',
 'C20':'760同密度/条件矩阵/总流的有限能源Gauss态可有不同来源流，有限能源不强制条件几何项消失',
 'C21':'760不允许将构形依赖右酉视为无代价的同态表示自由，它改变完整密度核',
 'C22':'760真实质量来源流含条件相干项eta，均值乘总流不足；原动能与实际记录注能仍同账'}
rows=ledger.splitlines();seen=[]
for i,row in enumerate(rows):
    for key,value in updates.items():
        if row.startswith('|'+key+' '):
            parts=row.split('|');parts[2]+='；'+value;rows[i]='|'.join(parts);seen.append(key)
assert set(seen)==set(updates)
ledger='\n'.join(rows)+'''
 
## 760共同动态对象与真正压缩的区别

原有限秩混合Gauss过程可用背景密度n和完整条件振幅Phi精确重构，曲目标动态、原字段记录与来源由同一对象给出，无需假定能源谱分离。此处接用成熟精确因子化方法，不把普遍恒等式报成新物理定理。

本项目新限制是：原753/758紧管内存在两份有限能源Gauss态，同n、逐点sigma和总流j，却有至少kappa/2的原质量来源—动量联合矩差。有限能源不能保证hbar平方乘条件几何项趋零。这排除指定摘要对扩大合法态族的统一来源合同，不是再次反证758特定初态的动态极限。

得到的是原过程的精确表示，复杂度未降低；空间度规、原群/维数/作用、连续与内部自主准备仍输入或开放。正替代是保量子慢变量的表示，不是完整Q到E宏观桥。

## 下一项：同一任务上的可删除资料

接[761入口](round761_drafts/STATUS.md)。先复用704/706的既有加权过程误差，核原记录、来源与物质反作用任务下真正可删除的方向。不能仅证明又一份恒等表示或一般范数界便声称宏观化；不能默认新增消相干、测量、热浴或谱隙。原空间和604、649/699范围保持。
'''
write('unified_physics_condition_ledger_760.md',ledger)
write('round760_drafts/research_note_760_draft.md',(HERE/'research_note_760.md').read_text('utf8'))
write('round761_drafts/STATUS.md','''# 第761轮入口：同一记录与来源任务下的实际压缩

接[760](../research_note_760.md)、[760条件账](../unified_physics_condition_ledger_760.md)及[共同候选](../round755_drafts/cognitive_joint_candidate_working_report.md)。

1. 760精确条件振幅保原完整动态、实际记录与来源，无需能带隙。它保全部量子慢变量，没有证明廉价宏观模型或Q到E。
2. 原紧管有限能源态可同n、sigma、j却异来源流。新增见证输入比758单相点准备宽；不得冒充其特定tau动态反例。
3. 下一项先声明原有限任务历史、来源菜单、输入及误差；复用704/706的来源加权收敛和625分支估计，不将一般范数不等式重算科学轮次。
4. 重点核一份有认知或物理动机的真实删除/组织方式，证明保任务或给受限反例。不能先加入免费能带测量、消相干、外浴，再声称原过程自动产生。
5. 若用hbar平方g趋零、粗时间、随机背景或记忆闭合，单列输入及验证对象。有限能源已不足以推出前者；不要重复扫相位频率或旧谱精度。
6. 保原同一作用、全部物种及边；避免只用新玩具替换原过程。所有新数值区分算术校验、局部模型和全图演化。
7. 原维数、群、作用、物种、外给空间gamma及连续/动力引力缺口继续在联合主账。旧382—386、425、522—523及604、649/699不改。
8. 不改目标，不建任务或调度，不做图像。若只能整理出条件而没有实质连接，写工作报告，不凑轮次。
''')
write('round760_drafts/scope_and_dedup_review.md','''# 760范围与去重审查

- 先读三导航、759笔记/结果、760入口及其数值；无在运行Python。旧编号科学文件和760冻结STATUS未覆盖。
- 回查精确因子化、条件矩阵、量子流体与纯化关键词；直接复用643、704/706、727、753、757—759。
- 文献精确因子化是成熟方法，本轮不报首次发现；将原有限秩混合、严格Gauss、曲目标和实际记录/来源放到同一振幅上。
- 实质摘要限制是原紧管两态同n、sigma、所有方向总流，却异质量来源流；这是更宽合法输入族，不是758单相点初态的又一个有限时间反证。
- 常数右酉不改密度核；构形依赖右酉会改变核。列标号不是免费物理辅助，准备仍未自主化。
- 用X流延伸等变框架使XT=0；Xr=1且原径向场与规范对易，质量迹在小管保持正。有限能源界不需要交换时间与半经典极限。
- 完整条件振幅保跨带相干但不降低求解复杂度。固定图精确改写不是连续、宏观、空间度规动力学或量子引力。
- 数值三组仅校验微分、CAR及原字段记录代数；未模拟整个周期图，未重复计数759迹公式。
''')
write('round760_drafts/literature_scope_audit.json',json.dumps(dict(
 sources=[dict(title='Correlated electron-nuclear dynamics: Exact factorization of the molecular wavefunction',
 url='https://arxiv.org/html/1208.4388',
 checked='SectionsII.1 andII.2: normalized conditional amplitude, marginal amplitude, exact coupled equations and phase freedom.',
 use='Mature factorization idea; current curved, mixed Gauss and source identities are derived in the note. No molecular time scale or automatic classical approximation imported.')],
 new='Original finite-energy compact Gauss states have identical configuration density, conditional matrix and total current but distinct original mass-source current; full conditional amplitudes preserve the original finite-history process.',
 excluded='No reduced-complexity macroscopic limit, no full graph time simulation, no new gravity equation or autonomous terminal apparatus.'
),ensure_ascii=False,indent=2)+'\n')
main=('research_note_760.md','joint_conditional_background.py','joint_conditional_background_results.json','unified_physics_condition_ledger_760.md')
write('round760_drafts/final_review.txt',
 'Primary-agent review. Exact mixed-state amplitudes retain the off-diagonal configuration kernel, unlike local conditional density matrices. Product-rule equations use the original curved measure and scalar kinetic principal part. Source current contains hbar Im<Phi,(C-c)X Phi>. The original mass has zero conditional trace and positive square trace; a configuration-dependent right phase gives two normal compact Gauss states with equal n,sigma,current and unequal source current, while the energy remains uniformly bounded. This phase is not a free purification gauge or a modification of the Hamiltonian. The flow-adapted equivariant frame satisfies XT=0; source continuity gives the integrated strict bound. Larger input family explicitly distinguished from758 packets. Exact full-amplitude dynamics is a representation, not a macroscopic simplification. Nodal sets use original W. Finite records and kinetic cost retained; no actual source-current measuring apparatus newly built.\n'+
 '\n'.join(n+' '+hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in main)+'\n')

mapping={'760':'761','759':'760','758':'759','3464':'3467','3461':'3464',
         '1589':'1592','3613':'3627','3600':'3613','405':'406','310':'311',
         'joint_background_source_closure':'joint_conditional_background'}
pub=remap((HERE/'publish_round759.py').read_text('utf8'),mapping)
summary='**第760轮完成：** [混合Gauss过程与条件来源流]({p}research_note_760.md)精确条件振幅保原动态/记录/来源，无需能带隙；原有限能源两态同密度、条件矩阵及总流却异来源流。三组、十六式通过，最新760／3467，1592份编号科学文件、3627份保护证据。[核验]({p}research_round_760_checks.json)、[条件账]({p}unified_physics_condition_ledger_760.md)。这是精确表示，尚非受控宏观或量子引力。'
order='**当前执行顺序（760后，优先于下方历史安排）：** 接[761同一任务与实际压缩]({p}round761_drafts/STATUS.md)，复用旧加权过程估计，检验原记录与来源任务下真正可删除的资料；不把完整相干重写算作宏观化。[范围审计]({p}round760_drafts/scope_and_dedup_review.md)。认知共同候选、旧空间、604、649／699及目标保持。'
pub=re.sub(r'^summary=.*$',lambda m:'summary='+repr(summary),pub,flags=re.M)
pub=re.sub(r'^order=.*$',lambda m:'order='+repr(order),pub,flags=re.M)
pub=pub.replace('**第760轮完成（限定反例）：**','**第760轮完成：**')
pub=pub.replace('确定背景动态闭合的联合来源检验','混合Gauss过程与条件来源流').replace('确定物质背景限制','条件相干与来源')
write('publish_round760.py',pub)
post=remap((HERE/'postcheck_round759.py').read_text('utf8'),mapping)
post=post.replace('**第760轮完成（限定反例）：**','**第760轮完成：**')
write('postcheck_round760.py',post)
ver=remap((HERE/'verify_round759.py').read_text('utf8'),mapping)
start=ver.index("    entry=core.read(");end=ver.index('    result=model.run();',start)
ver=ver[:start]+'''    entry=core.read(HERE/'round760_drafts/transport_entry_checks.json')
    for name in ('transport_source_entry.py','transport_source_entry_results.json','transport_applicability_working_report.md'):
        assert (HERE/'round760_drafts'/name).exists()
'''+ver[end:]
start=ver.index('    names=(');end=ver.index('    new=',start)
names=('unified_physics_condition_ledger_760.md','round760_drafts/research_note_760_draft.md',
       'round760_drafts/final_review.txt','round760_drafts/literature_scope_audit.json',
       'round760_drafts/scope_and_dedup_review.md','round761_drafts/STATUS.md',
       'round760_drafts/transport_applicability_working_report.md','round760_drafts/transport_source_entry.py',
       'round760_drafts/transport_source_entry_results.json','round760_drafts/publish_transport_entry.py',
       'round760_drafts/transport_entry_checks.json')
ver=ver[:start]+'    names='+repr(names)+'\n'+ver[end:]
ver=ver.replace('original_joint_integral_moment_identity_checked=True,stated_deterministic_material_background_contract_refuted=True,full_physical_theory_refuted=False,positive_dynamic_replacement_proven=False',
 'exact_original_conditional_amplitude_representation_checked=True,local_density_matrix_current_summary_insufficient_on_stated_input_family=True,reduced_complexity_macroscopic_limit_proven=False,full_continuum_process_equivalence_proven=False')
write('verify_round760.py',ver)
print('Prepared760 finite-graph representation and scoped summary witness.')

