"""Prepare757 native report/source bridge, preserving earlier working evidence."""
import hashlib
import json
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent

def write(name,text):
    p=HERE/name;p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)

def remap(text,mapping):
    return re.sub('|'.join(re.escape(k) for k in sorted(mapping,key=len,reverse=True)),
                  lambda m:mapping[m[0]],text)

# Preserve the initial unnumbered source before correcting two escaped controls.
note=HERE/'research_note_757.md';raw=note.read_text('utf8')
assert raw.count(chr(12))==2
write('round757_drafts/note_before_latex_escape_repair.txt',raw)
note.write_text(raw.replace(chr(12)+'rac',chr(92)+'frac'),encoding='utf8',newline='\n')

ledger=(HERE/'unified_physics_condition_ledger_756.md').read_text('utf8')
ledger=ledger.replace(ledger.split('\n')[0],'# 联合条件总账：757原相互作用、同一来源关联与实际报告',1)
ledger=ledger.replace(ledger.split('\n')[2],
 '2026-10-04。接[756全账](unified_physics_condition_ledger_756.md)，回填[757报告](research_note_757.md)。[结果](joint_correlation_native_report_results.json)、[核验](research_round_757_checks.json)。联合目标保持。',1)
updates={
 'C01':'757同二点的两份中性输入均是原正常全Gauss物理态；无需755非中性准自由替换',
 'C03':'757同一来源关联经原H进入原sin s联合报告，严格四阶差；保实际Kraus后态，不用理想占据测量',
 'C04':'757全H、固定边、完整物种保留，只有四阶差局部；后续过程与余项不被宣称独立于图',
 'C19':'757明确正常紧Gauss准备及关联资源；相同平均H仍有不同H方差，内部制备未构造',
 'C20':'757有限中性q连接756源字典与原报告，但不是两分支完整动态/尺度等价',
 'C21':'757同二点/均值摘要无法保所列原实际报告，须保关联或完整过程；不要求任意任务精确摘要',
 'C22':'757同一配对来源方差固定报告四阶系数及初始能源方差差；实际末读代价按原H计入'}
rows=ledger.splitlines();seen=[]
for i,row in enumerate(rows):
    for key,value in updates.items():
        if row.startswith('|'+key+' '):
            parts=row.split('|');parts[2]+='；'+value;rows[i]='|'.join(parts);seen.append(key)
assert set(seen)==set(updates)
ledger='\n'.join(rows)+'\n\n## 757原相互作用的共同来源—报告连接\n\n'
ledger+='''- 757入口证书的CP/全方差连接直接复用，不重新计轮次。本轮新核的是756同二点关联能否通过原H被原标量读口看到。
- 两份原正常Gauss输入共用全部二点、平均来源和玻色准备；同一q−p²固定配对来源方差、原报告四阶差与初始H方差差。
- 四阶系数保原完整图所有作用，仅凭CAR腿数和曲目标对易精确消去其他作用在该阶差中的贡献。禁止把这个局部系数理解为全动力学可删边。
- 原sin s两读的完整instrument有严格小正时间区分，正常准备支持直接给正号；无需新耦合或理想占据仪器。末读注能界保留，自治末读与制备尚未完成。
- q连接到756连续源的相同有限代数，但Q与E的完整状态/动力学仍未等同。维数、群、Y、作用与几何仍为输入。

## 本轮合并与下一项

本轮把同一内部来源关联和原实际报告接通，减少该有限任务的理想占据测量假设；不增加独立噪声或记录耦合。接[758](round758_drafts/STATUS.md)，返回原玻色—费米过程与有效连续来源的共同映射，不继续局部读口或时间窗口优化。应用目标、空间旧桥及604/649/699保持。
'''
write('unified_physics_condition_ledger_757.md',ledger)
write('round757_drafts/research_note_757_draft.md',note.read_text('utf8'))
write('round758_drafts/STATUS.md','''# 第758轮入口：原联合相互作用与有效来源的同一映射

接[757](../research_note_757.md)、[条件账](../unified_physics_condition_ledger_757.md)、[认知共同候选](../round755_drafts/cognitive_joint_candidate_working_report.md)。

1. 756同一中性q给背景Hadamard来源噪声，757原全图H把该q转成原sin s报告；后者不再需要理想占据测量。但两份玻色准备和完整动态过程尚未同一化。
2. 直接复用643原全Gauss有序历史、704/706固定图真实过程与来源极限、743原量子标量关联生成。固定图近似不是空间连续极限；不得抹掉604/649/699各自限定反例。
3. 下一项应核联合关联与完整来源在同一有效展开中的输送，找明确映射、受控近似或一个指定接法的失败。若需连续扰动工具，单列范围和正性，不把形式级数等同精确正态。
4. 不继续757读口寿命、微小时间窗或系数精度；这些不是当前主线。也不重复CP/全方差公式凑轮次。
5. 同一准备、内部交互、实际记录、后态、全部来源及资源要共同列明。末读与制备仍是输入，不能仅凭费用有限宣布封闭整体已实现。
6. 旧空间382—386、425、522—523保持，不重加Lipschitz；不改应用目标，不建新任务或调度，不做图像检查。
''')
write('round757_drafts/literature_scope_audit.json',json.dumps(dict(
 sources=[dict(title='Quantum fields and local measurements',url='https://arxiv.org/html/1810.06512v3',
  checked='Section3.3 and introductory scope: induced instruments; probe preparation/control/final readout are inputs; no quantum gravity or universal measurability theorem.',
  use='Distinguish actual Kraus process from effects and autonomous apparatus.'),
 dict(title='Stochastic Gravity: Theory and Applications',url='https://arxiv.org/html/0802.0658',
  use='Inherited source-noise/response framework; no imported gravitational action as cognition derivation.')],
 inherited='643 ordered Gauss history;704/706 fixed-graph process limit;717 force;743 fluctuation generation;746 commutator word decomposition;756 same-correlation noise.',
 new='Original full-graph fourth scalar-report difference for two strict-Gauss same-two-point neutral states, fixed by the same pair-source variance.',
 excluded='New probe coupling, full quantum time simulation, autonomous terminal implementation, full continuous state/process equivalence or gravity derivation.'
),ensure_ascii=False,indent=2)+'\n')
main=('research_note_757.md','joint_correlation_native_report.py','joint_correlation_native_report_results.json','unified_physics_condition_ledger_757.md')
write('round757_drafts/final_review.txt',
 'Primary-agent review only. Pi removes all constant/quadratic CAR terms; four selected legs require two local Majorana insertions in the two-B/two-kinetic words. Full graph H retained. Exact operator identity and rational support bound fix a nonzero fourth derivative, not a certified usable time window. Same initial means, differing H variance; correlation not a free preparation resource. Actual two-read Kraus process and original finite injection bound retained. Source variance to report is proved in finite Q; no full Q/E dynamical equality, autonomous detector or quantum gravity claim.\n'+
 '\n'.join(n+' '+hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in main)+'\n')
mapping={'757':'758','756':'757','755':'756','3455':'3458','3452':'3455',
 '1580':'1583','3574':'3589','3562':'3574','402':'403','307':'308',
 'joint_correlated_hadamard_noise':'joint_correlation_native_report'}
pub=remap((HERE/'publish_round756.py').read_text('utf8'),mapping)
summary='**第757轮完成：** [原相互作用、来源关联与实际报告]({p}research_note_757.md)两份同二点/平均来源的正常Gauss态，经原全H在原sin s联合报告上有严格四阶差，系数由同一配对来源方差固定；不增理想占据测量。三组、十八式通过，最新757／3458，1583份编号科学文件、3589份保护证据。[核验]({p}research_round_757_checks.json)、[条件账]({p}unified_physics_condition_ledger_757.md)。全连续过程与量子引力未完成。'
order='**当前执行顺序（757后，优先于下方历史安排）：** 接[758原相互作用与有效来源]({p}round758_drafts/STATUS.md)，复用643／704／706／743，核同一玻色—费米关联、实际记录与完整来源的有效映射；不继续读口时间窗或系数精度。[范围审计]({p}round757_drafts/scope_and_dedup_review.md)。认知共同候选、旧空间、604、649／699及应用目标保持。'
pub=re.sub(r'^summary=.*$',lambda m:'summary='+repr(summary),pub,flags=re.M)
pub=re.sub(r'^order=.*$',lambda m:'order='+repr(order),pub,flags=re.M)
pub=pub.replace('同一关联态、连续记录与来源涨落','原相互作用、来源关联与实际报告').replace('旧空间合同保持，同一记录与来源统计','旧空间合同保持，原关联与实际报告')
write('publish_round757.py',pub)
write('postcheck_round757.py',remap((HERE/'postcheck_round756.py').read_text('utf8'),mapping))
ver=remap((HERE/'verify_round756.py').read_text('utf8'),mapping)
start=ver.index('    names=(');end=ver.index('    new=',start)
names=('unified_physics_condition_ledger_757.md','round757_drafts/research_note_757_draft.md',
 'round757_drafts/final_review.txt','round757_drafts/literature_scope_audit.json',
 'round757_drafts/scope_and_dedup_review.md','round758_drafts/STATUS.md',
 'round757_drafts/common_instrument_certificate.md','round757_drafts/common_instrument_certificate_check.py',
 'round757_drafts/common_instrument_certificate_checks.json','round757_drafts/publish_common_instrument_entry.py',
 'round757_drafts/common_instrument_entry_navigation_checks.json','round757_drafts/note_before_latex_escape_repair.txt')
ver=ver[:start]+'    names='+repr(names)+'\n'+ver[end:]
ver=ver.replace('same_record_Hadamard_source_noise_connection_checked=True,conditional_joint_response_scope_checked=True,full_continuum_Gauss_state_constructed=False',
 'original_full_graph_source_to_report_identity_checked=True,strict_sign_and_actual_resource_scope_checked=True,full_continuum_process_equivalence_proven=False')
ver=ver.replace('    result=model.run();',"    import sys\n    sys.path.insert(0,str(HERE/'round757_drafts'))\n    import common_instrument_certificate_check as entry\n    assert entry.verify()==core.read(entry.TARGET)\n    result=model.run();")
write('verify_round757.py',ver)
print('Prepared757: original source correlation to native report, complete scopes.')
