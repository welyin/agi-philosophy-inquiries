"""Prepare750 same-record sources and joint normalization audit."""
import hashlib,json,re
from pathlib import Path
HERE=Path(__file__).resolve().parent
def write(name,text):
    p=HERE/name;p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)
def remap(text,mapping):
    return re.sub('|'.join(re.escape(k) for k in sorted(mapping,key=len,reverse=True)),lambda m:mapping[m[0]],text)
ledger=(HERE/'unified_physics_condition_ledger_749.md').read_text('utf8')
ledger=ledger.replace('# 联合条件总账：749合并映射与非线性几何约束','# 联合条件总账：750实际标记来源与共同概率',1)
ledger=ledger.replace(ledger.split('\n')[2],'2026-10-04。接[749全账](unified_physics_condition_ledger_749.md)，回填[750报告](research_note_750.md)。[结果](joint_marked_source_completion_results.json)、[核验](research_round_750_checks.json)。联合目标保持。',1)
updates={'C03':'750原占据instrument的实际条件两点函数与完整标记来源明确；概率族须共同完成',
'C10':'750固定共同背景的逐记录首阶响应可加权；分别正的几何/参考分支不能自动拼成一次实验',
'C19':'750原背景依赖概率要求同一准备历史；逐结果选择不同过去会破坏归一',
'C22':'750实际24来源符号及条件源权重导数共同核验；原连续Ward差核与一阶响应复用'}
rows=ledger.splitlines();seen=[]
for i,row in enumerate(rows):
    for key,value in updates.items():
        if row.startswith('|'+key+' '):
            parts=row.split('|');parts[2]+='；'+value;rows[i]='|'.join(parts);seen.append(key)
assert set(seen)==set(updates)
ledger='\n'.join(rows)+'\n'
ledger=ledger.rsplit('## 本轮合并与下一项',1)[0]+'''## 750同一实际记录与共同概率完成

- 不再用749外给source标签充当实际量子记录。原730 occupancy仪器的P0/P1由同一纯准自由P明确给出，原全种128维对象保留。
- 所有双线性源从同一Pr−P得到；24个原局部源符号以独立4维Fock核验，不能将其当完整连续PDE数值解。
- 同一记录间来源核为p(1−p)d_alpha(x)d_beta(y)，不是全部量子噪声。634总体荷协方差与730/732局域光滑/Ward直接复用。
- 背景变化时须保概率权重导数；原能源中该项1.9264e−7，不能视为恒等于零。
- 741同一固定背景的一阶响应可逐条件态连接；但结果依赖背景历史时，共同归一新增必要式sum_r Dp_r[b_r]=0。
- 原历史±gamma诊断给概率总和不等于1；反馈规则只是声明的试接，不是已解出的Einstein分支。
- 固定共同准备先记录、后条件保迹演化保原权重，但该代数身份没有构造自主因果装置或约束相容未来。
- 禁止用依未知输入的事后除Z冒充确定CPTP完成；条件化不等于物理坍缩。

## 本轮合并与下一项

C03/C10/C19/C22的同一实际来源、记录权重与候选几何提升共同受限。减少的是独立选择分支源/准备后再拼接的自由，不是已消去维数、群、参数或连续输入。

接[751](round751_drafts/STATUS.md)：以一份共同过去和原权重审查原记录后的来源、初值补偿及未来响应；区分条件资料与可实施过程。旧空间、604、649/699和747/748局部总符号缺口均保持。
'''
write('unified_physics_condition_ledger_750.md',ledger)
write('round750_drafts/research_note_750_draft.md',(HERE/'research_note_750.md').read_text('utf8'))
write('round751_drafts/STATUS.md','''# 第751轮入口：共同过去、实际来源与未来响应

接[750](../research_note_750.md)、[全账](../unified_physics_condition_ledger_750.md)。

1. 750把原实际条件后态接到完整来源，并发现逐结果改变准备历史再拼概率不自动归一。共同过去/同一instrument完成条件必须保留。
2. 复用634真实记录的供能、动量与同记录条件下补偿关系；不重报一般测量扰动或归一恒等式。
3. 核730—741原条件源与731初值补偿：数学条件解族、主体掌握不同资料后的推断、可实施条件反馈、共同量子几何必须分别声明。
4. 不要求每个主体的量子态描述相同。要保的是同一个物理过程与可比较事件，不能以各自有利参考取代共同准备。
5. 优先检验保原未知输入时，几何/资源读数能否由同一线性正过程产生；不能将由未知条件态期望算出的经典背景当已实现的输出。
6. 若触及平均来源的非线性或经典-量子映射，先检索旧593、620、631、634、649、708、741及成熟理论；一般无克隆/非线性警告不另计轮次。
7. 每项新增须限制跨部门接口或排除具体原候选。不要无限延伸局部读口优化，旧空间及604、649/699边界保持。
''')
write('round750_drafts/literature_scope_audit.json',json.dumps(dict(sources=[],
inherited='634 actual record/source covariance and process boundaries;730 conditional smooth differences;732 joint Ward;735 shared normalization;741 first-order response;742 pure-reference partner.',
new='Explicit original conditional source bundle, probability-weight geometry derivative, and joint completion audit for branch-dependent histories.',
excluded='Causal detector realization; actual nonlinear conditional geometry; full quantum source or continuum limit.'),ensure_ascii=False,indent=2)+'\n')
write('round750_drafts/scope_and_dedup_review.md','''# 750共同对象与范围审查

上一目标轮747—749完成并冻结，属于有效进展。本轮先读四份导航、749笔记/结果与750入口，未发现活跃Python进程。无新代理、任务、定时或图像。

- CAR/Wick条件Pr公式以原纯参考伙伴的独立Fock投影核验，并逐一检验24源，非仅非选择两点匹配。
- Pr是两点函数，非选择完整态不因此Gaussian；742非高斯界保持。
- p在固定非零邻域内才有统一分母控制；未覆盖任意罕见记录。
- 实际源差的光滑性及Ward继承730/732，绝对源须沿735同一规范化。
- 原颜色不变参考和singlet记录的one-point颜色源零，使731当前颜色积分条件成立；不构造全量子Gauss几何态。
- 加权一阶响应仅同一固定线性右逆；未由此得到非线性共同几何。
- 背景变化时保Pr、J及p的导数。归一缺陷的两种规则为试接历史，不是宣称实际GR解失败。
- 逐分支CP/正性不推出总instrument归一；按输入重归一不是确定CPTP。
- 固定过去后反馈的公式只是概率完成结构，内部资源、因果与约束仍缺。
- 新结果限制共同模型拼接，不重算634来源噪声或旧空间，目标未改。
''')
main=('research_note_750.md','joint_marked_source_completion.py','joint_marked_source_completion_results.json','unified_physics_condition_ledger_750.md')
write('round750_drafts/final_review.txt','Primary-agent review only. Actual conditional sources, weight derivatives and joint normalization condition. Separate conditional backgrounds are not automatically one quantum experiment. No nonlinear geometry or detector completion.\n'+
'\n'.join(n+' '+hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in main)+'\n')
mapping={'750':'751','749':'750','748':'749','3436':'3439','3434':'3436','1559':'1562',
'3493':'3506','3480':'3493','395':'396','300':'301',
'joint_record_geometry_averaging':'joint_marked_source_completion'}
pub=remap((HERE/'publish_round749.py').read_text('utf8'),mapping)
summary='**第750轮完成：** [实际标记来源与共同概率]({p}research_note_750.md)原条件后态固定完整来源和概率权重响应；分别合法的背景/准备分支仍须满足同一instrument归一。原历史试接出现非零概率缺陷，未宣称实际GR分支失败。三组、十六式通过，最新750／3439，1562份编号科学文件、3506份保护证据。[核验]({p}research_round_750_checks.json)、[条件账]({p}unified_physics_condition_ledger_750.md)。'
order='**当前执行顺序（750后，优先于下方历史安排）：** 接[751共同过去与未来响应]({p}round751_drafts/STATUS.md)，保原未知输入及共同记录权重，核来源/初值/未来的同一过程；不把条件解族当真实实施。统一目标和旧空间、604、649／699保持。'
pub=re.sub(r'^summary=.*$',lambda m:'summary='+repr(summary),pub,flags=re.M)
pub=re.sub(r'^order=.*$',lambda m:'order='+repr(order),pub,flags=re.M)
pub=pub.replace('合并映射与非线性几何约束','实际标记来源与共同概率').replace('旧空间合同保持，合并资料与非线性约束共同检验','旧空间合同保持，原条件源与共同概率相容')
pub=pub.replace('**第750轮完成（限定反例）：**','**第750轮完成：**')
write('publish_round750.py',pub)
post=remap((HERE/'postcheck_round749.py').read_text('utf8'),mapping).replace('**第750轮完成（限定反例）：**','**第750轮完成：**')
write('postcheck_round750.py',post)
verification=remap((HERE/'verify_round749.py').read_text('utf8'),mapping)
start=verification.index('    names=(');end=verification.index('    new=',start)
names=('unified_physics_condition_ledger_750.md','round750_drafts/research_note_750_draft.md',
'round750_drafts/final_review.txt','round750_drafts/literature_scope_audit.json',
'round750_drafts/scope_and_dedup_review.md','round751_drafts/STATUS.md',
'round750_drafts/marked_source_entry.py','round750_drafts/marked_source_entry_results.json',
'round750_drafts/branch_probability_audit.py','round750_drafts/branch_probability_audit_results.json')
verification=verification[:start]+'    names='+repr(names)+'\n'+verification[end:]
verification=verification.replace("==(2,0,0)","==(3,0,0)").replace("fresh_tests=dict(run=2","fresh_tests=dict(run=3")
verification=verification.replace("checks['display_formulas']==12","checks['display_formulas']==16")
verification=verification.replace('original_nonlinear_mean_map_counterexample_checked=True,actual_quantum_source_realization_proven=False',
'actual_marked_source_and_joint_normalization_checked=True,nonlinear_common_geometry_realized=False')
write('verify_round750.py',verification)
print('Prepared750 actual source bundle and joint completion audit.')
