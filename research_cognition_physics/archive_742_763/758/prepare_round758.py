"""Prepare the initial matrix bridge; never overwrite existing artifacts."""
import hashlib
import json
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent

def write(name,text):
    p=HERE/name;p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)

def remap(text,mapping):
    return re.sub('|'.join(re.escape(k) for k in sorted(mapping,key=len,reverse=True)),
                  lambda m:mapping[m[0]],text)

ledger=(HERE/'unified_physics_condition_ledger_757.md').read_text('utf8')
ledger=ledger.replace(ledger.split('\n')[0],'# 联合条件总账：758同一Gauss背景中的量子矩阵与初始任务',1)
ledger=ledger.replace(ledger.split('\n')[2],
 '2026-10-04。接[757全账](unified_physics_condition_ledger_757.md)，回填[758报告](research_note_758.md)。[结果](joint_gauss_matrix_background_results.json)、[核验](research_round_758_checks.json)。统一目标保持。',1)
updates={
 'C01':'758原严格Gauss等距嵌入保稳定子不变的整个有限Fock纤维，含原中性关联族；不要求单带或谱隙',
 'C03':'758同映射保初始光滑instrument、未归一后态及来源矩；含等待的真实过程尚未输送',
 'C04':'758原全H给初始矩阵残差与有界能源谱测试；物理时间需除以hbar，现有残差不足',
 'C09':'758原753非零颜色的完整周期环路稳定子只留中心，精确采样保Gauss；未证明量子Einstein约束',
 'C19':'758正常Gauss准备可置于753同背景中心，原757实际报告正号保留；未给内部制备',
 'C20':'758固定图半经典初始桥，常数依图；非固定物理hbar连续场论、动态Q/E同一映射',
 'C21':'758保完整有限量子矩阵及关联，避免只取均值或基态；成本及有效时间仍需核',
 'C22':'758同一初始菜单的均值和二阶来源矩共同输送，允许非零量子来源涨落'}
rows=ledger.splitlines();seen=[]
for i,row in enumerate(rows):
    for key,value in updates.items():
        if row.startswith('|'+key+' '):
            parts=row.split('|');parts[2]+='；'+value;rows[i]='|'.join(parts);seen.append(key)
assert set(seen)==set(updates)
ledger='\n'.join(rows)+'''

## 758同一Gauss矩阵初始桥

- H1/H3选择保量子内部资料的宏观背景分支，H2要求状态、记录及来源同对象。不是另添认知基本公理。
- 原非零颜色周期环路、574原电弱补全及其稳定子条件容许同一精确玻色Gauss相点；全有限Fock的稳定子不变子空间可严格等距提升。保原维数、群、作用、离散化及耦合输入。
- 原完整H、规范不变有限来源菜单及初始光滑instrument有同一矩阵初始极限；完整关联态无需有隙单带。去掉该初始桥的单带/谱隙要求，未消去物理结构包。
- 原q族的同二点和不同能源/来源方差均保留。背景经典不等于来源成为无涨落数字；这部分沿726/727，不重新计其一般反例。
- 有界能源谱参数不能当物理时间；O(sqrt(hbar))初始残差并未给固定正时间的动态误差界。未构造量子Einstein约束、连续全Gauss态、自治准备或完整Q/E等价。

## 下一项：动态合同而非局部精度

接[759](round759_drafts/STATUS.md)，先回查旧慢快、共同过程和来源结果，核原关联输入是否需要随时间变化的背景分支/记忆，能否保实际记录与完整来源。固定谱测试或初始态相容不重复扩成新轮。空间382—386、425、522—523及604/649/699的限定范围保持；应用目标未修改。
'''
write('unified_physics_condition_ledger_758.md',ledger)
write('round758_drafts/research_note_758_draft.md',(HERE/'research_note_758.md').read_text('utf8'))
write('round759_drafts/STATUS.md','''# 第759轮入口：同一量子矩阵与动态背景的合同

接[758初始桥](../research_note_758.md)、[758条件账](../unified_physics_condition_ledger_758.md)、[认知共同候选](../round755_drafts/cognitive_joint_candidate_working_report.md)。

1. 758已保同一Gauss初始准备、稳定子不变量子矩阵、来源菜单及初始光滑记录；无需有隙单带。不要重新证明初始态相容或重复谱参数扫描。
2. 物理演化有1/hbar，初始O(sqrt(hbar))残差不足以给固定时间误差。应明确实际动态合同：同一输入族、时间/分辨率、记录、完整来源与误差。
3. 先回查525、558、591、643、704/706、726—730、741、743及757，区分旧单带/均场限制和新联合过程问题。不能以此前来源有涨落再次凑轮次。
4. 优先核原持续相互作用是否需要多背景分支或记忆，保原量子关联及来源。若引入新的缩放、谱隙、平均或随机过程，分别登记；不把任意混合动力学默认成原Q的有效极限。
5. 有效初始矩阵不等于原完整过程；原有限图也不等于连续手征场或量子Einstein态。限定失败不得推广为全部统一计划失败。
6. 不优化局部读口/系数；不改目标，不建任务/调度，不做图像。空间382—386、425、522—523及604、649/699保持。
''')
write('round758_drafts/scope_and_dedup_review.md','''# 758范围与去重核查

- 初始认知选择：H1/H3允许共同背景与内部差异并存，H2要求来源同态；并未要求所有微观来源方差为零。
- 已回查574管状波包、726空CAR残差、727单带与带间来源、643全Gauss历史、704/706固定图过程、753背景、756关联噪声及757实际报告。
- 新增是稳定子不变整个有限Fock纤维的原Gauss等距桥及同一初始矩阵任务；不是重新证明一般Duhamel、Gauss投影或Born-Oppenheimer理论。
- 原753颜色场用完整轴向环路核稳定子；商群中心lift不产生额外非中心稳定子，离散中心仍保留。未把颜色Lie代数零核当全群平凡。
- 源矩由算符作用残差与内积控制，不能由记录迹范数直接推无界源。原二维径向数值只校准局部组件；全图结论由证明承担。
- 756/757原关联态通过同一桥保留；757正系数可在753字段管内成立。两极限不能交换，未给hbar一致的记录时间窗。
- u为谱测试变量，t/hbar才是物理时间参数；初始桥不足以证明动态闭合。连续Gauss/来源/几何过程仍开放。
- 后继进入动态合同和来源共同输送，不继续包宽、谱参数或单一读口优化。旧空间及限定失败保留。
''')
write('round758_drafts/literature_scope_audit.json',json.dumps(dict(
 sources=[dict(title='The time-dependent Born-Oppenheimer approximation',
 url='https://arxiv.org/html/0712.4369',checked='Sections1-3: projected Hamiltonians, gap hypotheses and actual time approximation are distinct.',
 use='Avoid promoting an initial isometric matrix/spectral bridge to physical-time accuracy; no molecular hierarchy imported.')],
 inherited='574 tubular Gauss packet;726 persistent fermion variance;727 one-band limitations;643/704/706 full fixed-graph objects;753 background;756/757 correlation dictionary.',
 new='Same stabilizer-invariant matrix initial isometry, original colored background applicability, initial source/record transport and bounded spectral tests.',
 excluded='Physical-time dynamics, graph-uniform constants, continuum quantum field limit, autonomous preparation, dynamic quantum Einstein constraints.'
),ensure_ascii=False,indent=2)+'\n')
main=('research_note_758.md','joint_gauss_matrix_background.py','joint_gauss_matrix_background_results.json','unified_physics_condition_ledger_758.md')
write('round758_drafts/final_review.txt',
 'Primary-agent review only. The whole stabilizer-fixed Fock fiber is reducing for each invariant principal-symbol matrix at the chosen phase orbit. Its equivariant tube frame is well-defined and smooth without a spectral gap. Fixed-graph Gaussian derivative bounds give operator-norm initial residuals, joint source second moments and initial instrument transport. Energy characteristic tests use bounded spectral u; t/hbar yields no fixed-time conclusion. Actual nonzero color is checked by full torus loops, retaining center and quotient restrictions. Original neutral correlations remain inside this same map. Local numerical fixtures do not stand for a full-graph simulation. Quantum Einstein constraints and continuum process still open.\n'+
 '\n'.join(n+' '+hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in main)+'\n')
mapping={'758':'759','757':'758','756':'757','3458':'3461','3455':'3458',
 '1583':'1586','3589':'3600','3574':'3589','403':'404','308':'309',
 'joint_correlation_native_report':'joint_gauss_matrix_background'}
pub=remap((HERE/'publish_round757.py').read_text('utf8'),mapping)
summary='**第758轮完成：** [同一Gauss背景中的量子矩阵与初始记录]({p}research_note_758.md)原严格Gauss等距桥同时保量子关联、能源/来源矩及初始光滑记录，无需有隙单带；有限物理时间仍待证。三组、十六式通过，最新758／3461，1586份编号科学文件、3600份保护证据。[核验]({p}research_round_758_checks.json)、[条件账]({p}unified_physics_condition_ledger_758.md)。固定图初始连接不是完整连续或量子引力。'
order='**当前执行顺序（758后，优先于下方历史安排）：** 接[759动态背景与同一量子矩阵]({p}round759_drafts/STATUS.md)，回查旧慢快及来源结果，核同一关联输入、实际过程及背景/记忆的动态合同；不继续初始谱参数或包宽扫描。[范围审计]({p}round758_drafts/scope_and_dedup_review.md)。认知共同候选、旧空间、604、649／699及应用目标保持。'
pub=re.sub(r'^summary=.*$',lambda m:'summary='+repr(summary),pub,flags=re.M)
pub=re.sub(r'^order=.*$',lambda m:'order='+repr(order),pub,flags=re.M)
pub=pub.replace('原相互作用、来源关联与实际报告','同一Gauss背景中的量子矩阵与初始记录').replace('旧空间合同保持，原关联与实际报告','旧空间合同保持，量子矩阵初始桥')
write('publish_round758.py',pub)
write('postcheck_round758.py',remap((HERE/'postcheck_round757.py').read_text('utf8'),mapping))
ver=remap((HERE/'verify_round757.py').read_text('utf8'),mapping)
start=ver.index('    names=(');end=ver.index('    new=',start)
names=('unified_physics_condition_ledger_758.md','round758_drafts/research_note_758_draft.md',
 'round758_drafts/final_review.txt','round758_drafts/literature_scope_audit.json',
 'round758_drafts/scope_and_dedup_review.md','round759_drafts/STATUS.md',
 'round758_drafts/matrix_background_entry.py','round758_drafts/matrix_background_entry_results.json')
ver=ver[:start]+'    names='+repr(names)+'\n'+ver[end:]
ver=ver.replace('import common_instrument_certificate_check as entry','import matrix_background_entry as entry').replace('assert entry.verify()==core.read(entry.TARGET)','assert entry.run()==core.read(entry.TARGET)')
ver=ver.replace("checks['display_formulas']==18","checks['display_formulas']==16")
ver=ver.replace('original_full_graph_source_to_report_identity_checked=True,strict_sign_and_actual_resource_scope_checked=True,full_continuum_process_equivalence_proven=False',
 'original_Gauss_matrix_initial_bridge_checked=True,initial_source_record_and_spectral_scope_checked=True,finite_physical_time_equivalence_proven=False,full_continuum_process_equivalence_proven=False')
write('verify_round758.py',ver)
print('Prepared758 initial Gauss matrix bridge and759 dynamic entry.')
