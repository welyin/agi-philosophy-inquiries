"""Prepare749 restricted counterexample and continuation of the joint goal."""
import hashlib,json,re
from pathlib import Path
HERE=Path(__file__).resolve().parent
def write(name,text):
    p=HERE/name;p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)
def remap(text,mapping):
    return re.sub('|'.join(re.escape(k) for k in sorted(mapping,key=len,reverse=True)),lambda m:mapping[m[0]],text)
ledger=(HERE/'unified_physics_condition_ledger_748.md').read_text('utf8')
ledger=ledger.replace('# 联合条件总账：748原标量边与异地总系数','# 联合条件总账：749合并映射与非线性几何约束',1)
ledger=ledger.replace(ledger.split('\n')[2],'2026-10-04。接[748全账](unified_physics_condition_ledger_748.md)，回填[749报告](research_note_749.md)。[结果](joint_record_geometry_averaging_results.json)、[核验](research_round_749_checks.json)。当前联合目标保持。',1)
updates={'C02':'749直接平均psi及canonical字段的指定合并映射不保原Hamiltonian约束；不排除全部几何粗化',
'C03':'749复用634实际记录分支噪声，未把新外给源标签当实际仪器',
'C10':'749原731正能量响应族给合并后严格正二阶约束缺陷，限既定Einstein初值模型',
'C19':'749外给source族与实际记录的同一性仍未证，不允许分别挑有利参考',
'C22':'749保平均来源不足保非线性几何，需共同涨落或修改合并／有效约束'}
rows=ledger.splitlines();seen=[]
for i,row in enumerate(rows):
    for key,value in updates.items():
        if row.startswith('|'+key+' '):
            parts=row.split('|');parts[2]+='；'+value;rows[i]='|'.join(parts);seen.append(key)
assert set(seen)==set(updates)
ledger='\n'.join(rows)+'\n'
ledger=ledger.rsplit('## 本轮合并与下一项',1)[0]+'''## 749主体资料合并与原非线性几何约束

- 634已给实际记录的来源分支协方差，706/724已保真实记录与拼接；本轮不重复这些结果。
- 检验明确映射：遗忘等权来源标签，直接平均psi、物质canonical动量和共形剪切，保持原Hamiltonian约束。
- 731原受约束非平坦背景可取能量相对源±epsilon rho_star，rho_star=2C+4A psi^-12+4Y psi^-8，其他来源及动量不变。
- 一阶响应u=psi。二阶S2=10C psi^5+68A psi^-7+52Y psi^-3严格正，平均后的物理约束缺陷为epsilon² psi^-5 S2+O(epsilon⁴)，故该映射不闭合。
- 原15³模型和另一已有多来源族复算与解析式相容；严格符号不依赖网格极值。
- 排除的是指定映射，不是全部粗化、认知原则、GR或统一候选。给定Einstein作用仍是输入。
- 新source族尚未由实际量子instrument实现；不把条件化解释为物理坍缩。
- 741一阶平均相容性保持，不能不经二阶检验提升到非线性几何。合并合同与来源、物质能源和几何必须同时满足。

## 本轮合并与下一项

C02/C03/C10/C19/C22得到一个有明确正缺陷的共同筛查：单保平均资料不够，需保联合相关或改变合并/有效约束。没有新设独立通信耦合或几何涨落参数。

接[750](round750_drafts/STATUS.md)：审计634实际标记来源与731/741的共同过程身份，核同一记录、补偿及响应是否可联立；不将外给source标签冒充物理记录。747/748总读口符号仍为局部开放，旧空间、604、649/699均保持。
'''
write('unified_physics_condition_ledger_749.md',ledger)
write('round749_drafts/research_note_749_draft.md',(HERE/'research_note_749.md').read_text('utf8'))
write('round750_drafts/STATUS.md','''# 第750轮入口：实际标记来源与同一几何响应

接[749](../research_note_749.md)和[全账](../unified_physics_condition_ledger_749.md)。

1. 749严格排除原731外给源族上的指定平均psi/canonical资料映射；不等于实际量子记录分支已生成这些几何。
2. 直接复用634实际能源—动量分支协方差及其不等于完整局域噪声的限制；不再重报条件化和物理坍缩的区别。
3. 回查730—741真实连续记录、Hadamard族、同一绝对源与约束补偿；哪些结果逐分支可用、哪些仅非选择或需额外统一域，要明确。
4. 寻找一项能共同约束记录、资源、来源与几何的实际合同。优先核同一marked source的光滑性、守恒及初值，不通过另造外给source回避原过程。
5. 若同一实际记录只确定全空间荷而不确定局域来源，保留真正缺口并查原波包/二点函数；不能由634两维协方差猜完整时空噪声核。
6. 只在确有新共同接口或限定反例时另计轮次，不反复验证Taylor/Jensen恒等式，不无限设计装置。
7. 旧空间382—386、425、522—523与604、649/699保持；747/748总信号严格符号仍开放。应用目标不变，不新任务或定时安排。
''')
write('round749_drafts/literature_scope_audit.json',json.dumps(dict(sources=[],
inherited='634 real record/source covariance;706 actual cq transport;708 common fluctuation data;724/725 material regions;731 nonlinear constraint existence and smooth dependence;741 first-order response.',
new='Explicit positive second-order mean-psi erasure defect on the original731 background.',
excluded='Actual quantum realization of prescribed source family; all coarse-grainings no-go; physical collapse; derivation of GR.'),ensure_ascii=False,indent=2)+'\n')
write('round749_drafts/scope_and_dedup_review.md','''# 749主线与推导审查

先查634，确认记录与来源协方差已有，未重复。706/708/724/725、731/741与空间旧结果直接复用。无新代理、应用任务、定时或图像检查。

- 原731 C>0、Y>0及A>=0给明确合法双向小来源族，固定配置及canonical物质动量，所以Gauss/动量保持。
- 来源rho_star依同一原解及系数定义；允许的只是外给相对源，不声称任意量子态可实现。
- L psi =2rho_star psi使u=psi；原正负幂展开所有项保留，S2系数为10、68、52且严格正。
- 取psi±平均只消奇数项，残差为epsilon²S2；物理残差还须乘平均psi^-5。
- 光滑参数依赖和紧初片保证统一余项，故正号是解析结论，不靠配点证明。
- 限制的是特定共形变量平均，不是平均度规、全部几何粗化或FUCP。实际条件记录族仍待接通。
- 634非零总体荷协方差不能被冒充局域应力核。下一轮只接同一对象，不重新求一般线性波方程。
- 741一阶结论不受影响；新共同限制在非线性二阶。749的新报告不重开局部读口高阶优化。
''')
main=('research_note_749.md','joint_record_geometry_averaging.py','joint_record_geometry_averaging_results.json','unified_physics_condition_ledger_749.md')
write('round749_drafts/final_review.txt','Primary-agent review only. Restricted mean-psi/canonical-field erasure counterexample, with analytic positive coefficient. Prescribed sources are not actual record branches; full joint goal remains open.\n'+
'\n'.join(n+' '+hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in main)+'\n')
mapping={'749':'750','748':'749','747':'748','3434':'3436','3432':'3434','1556':'1559',
'3480':'3493','3468':'3480','394':'395','299':'300',
'joint_remote_total_coefficient':'joint_record_geometry_averaging'}
pub=remap((HERE/'publish_round748.py').read_text('utf8'),mapping)
summary='**第749轮完成（限定反例）：** [合并映射与非线性几何约束]({p}research_note_749.md)原731两份合法外给来源解，遗忘标签后直接平均psi及canonical资料，会有严格正二阶Hamiltonian约束缺陷；只排除该合并映射，不把外给源当实际记录或否定全部粗化。两组、十二式通过，最新749／3436，1559份编号科学文件、3493份保护证据。[核验]({p}research_round_749_checks.json)、[条件账]({p}unified_physics_condition_ledger_749.md)。'
order='**当前执行顺序（749后，优先于下方历史安排）：** 接[750实际标记来源与共同几何响应]({p}round750_drafts/STATUS.md)，复用634、731及741，核同一记录、来源、补偿与几何的对象身份；不重新证明已知记录协方差。统一目标与旧空间、604、649／699保持。'
pub=re.sub(r'^summary=.*$',lambda m:'summary='+repr(summary),pub,flags=re.M)
pub=re.sub(r'^order=.*$',lambda m:'order='+repr(order),pub,flags=re.M)
pub=pub.replace('原标量边与异地总系数','合并映射与非线性几何约束').replace('旧空间合同保持，原总响应接回联合合同','旧空间合同保持，合并资料与非线性约束共同检验')
pub=pub.replace('**第749轮完成（合成公式）：**','**第749轮完成（限定反例）：**')
write('publish_round749.py',pub)
post=remap((HERE/'postcheck_round748.py').read_text('utf8'),mapping).replace('**第749轮完成（合成公式）：**','**第749轮完成（限定反例）：**')
write('postcheck_round749.py',post)
verification=remap((HERE/'verify_round748.py').read_text('utf8'),mapping)
start=verification.index('    names=(');end=verification.index('    new=',start)
names=('unified_physics_condition_ledger_749.md','round749_drafts/research_note_749_draft.md',
'round749_drafts/final_review.txt','round749_drafts/literature_scope_audit.json',
'round749_drafts/scope_and_dedup_review.md','round750_drafts/STATUS.md',
'round749_drafts/joint_coarse_source_entry.py','round749_drafts/joint_coarse_source_entry_results.json',
'round749_drafts/positive_constraint_defect.py','round749_drafts/positive_constraint_defect_results.json')
verification=verification[:start]+'    names='+repr(names)+'\n'+verification[end:]
verification=verification.replace('original_scalar_force_and_total_coefficient_formula_checked=True,total_remote_signal_proven=False',
'original_nonlinear_mean_map_counterexample_checked=True,actual_quantum_source_realization_proven=False')
write('verify_round749.py',verification)
print('Prepared749 restricted map counterexample and750 actual-record contract.')
