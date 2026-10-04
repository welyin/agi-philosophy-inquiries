"""Freeze747's coefficient result; preserve the unresolved full remote response."""
import hashlib,json,re
from pathlib import Path
HERE=Path(__file__).resolve().parent
def write(name,text):
    p=HERE/name;p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)
def remap(text,mapping):
    return re.sub('|'.join(re.escape(k) for k in sorted(mapping,key=len,reverse=True)),lambda m:mapping[m[0]],text)
ledger=(HERE/'unified_physics_condition_ledger_746.md').read_text('utf8')
ledger=ledger.replace('# 联合条件总账：746原连通图中的占据读出','# 联合条件总账：747原跳跃与异地读口六阶贡献',1)
ledger=ledger.replace(ledger.split('\n')[2],'2026-10-04。接[746全账](unified_physics_condition_ledger_746.md)，回填[747报告](research_note_747.md)。[结果](joint_remote_record_response_results.json)、[核验](research_round_747_checks.json)。目标及当前任务保持不变。',1)
updates={'C03':'747计算旧远端T菜单的跳跃六阶分量，完整总效果仍待合并',
'C15':'747原跳跃向远端实际读口的非零分量已严格认证；未签收总通信',
'C19':'747紧支撑正常Gauss准备的分量可继承Gaussian严格符号',
'C22':'747所有实际后态及旧注能界继续保留；原标量边不能删除'}
rows=ledger.splitlines();seen=[]
for i,row in enumerate(rows):
    for key,value in updates.items():
        if row.startswith('|'+key+' '):
            parts=row.split('|');parts[2]+='；'+value;rows[i]='|'.join(parts);seen.append(key)
assert set(seen)==set(updates)
ledger='\n'.join(rows)+'\n'
ledger=ledger.rsplit('## 本轮合并与下一项',1)[0]+'''## 747原跳跃的异地读口贡献

- 原空／对未知输入、完整有限图和共同几何保持；lambda只是算符计数标记，不是新的可调物理资源。
- 两动能、两质量、两跳跃给原远端T读口的六阶分量；A准备反演对称和CAR奇偶使其约化到B完整质量与原sterile边。
- 全五坐标、整数CAR和有理多项式无截断计算；原复spin边逐项系数为单位spin的177/1000。
- 新密度不是746局部P4的倍数，旧占据转移与本地读取概率不能直接相乘。
- 严格区间给J6约−0.50536308765372，紧支撑物理准备可继承此分量符号。
- 原标量边在相同六阶可能贡献a6；d6(1)=a6+c6的总符号尚未算出。不宣布完整异地记录。
- 原完整输出、资源、来源及末端仪器边界保持；未新增物种、免费测量、弱边或引力假设。

## 本轮合并与下一项

C03／C15获得原跳跃到实际远端菜单的严格分量，但整体通信合同尚未关闭。接[748](round748_drafts/STATUS.md)求原标量边及同一总响应。旧空间、604、649／699保持，不重开已经消去的空间假设。
'''
write('unified_physics_condition_ledger_747.md',ledger)
write('round747_drafts/research_note_747_draft.md',(HERE/'research_note_747.md').read_text('utf8'))
write('round748_drafts/STATUS.md','''# 第748轮入口：原标量边与实际异地总响应

接[747](../research_note_747.md)、[全账](../unified_physics_condition_ledger_747.md)。

1. 747已严格计算原跳跃的六阶贡献，尚非完整异地读口总和。保持原未知空／对输入、原全H、原固定几何与既有T菜单。
2. 回查577、746，不重复相位通信或局部四阶证明；复用局部P4密度。
3. 从原Duhamel或等价算符词求原标量测地边对异地读取的同阶贡献，证明时间次序系数，核其他词是否消失。
4. 不将原标量边设为零、不独立调弱k，不把形式lambda调节当可用资源；在同一w、k、tau下求总和。
5. 数值收敛与严格区间分开。若有抵消点，记录适用范围和下一非零阶；未有严格符号就不签收实际通信。
6. 保Gauss正常准备、紧支撑核、完整CP后态及共同来源代价。自治终端、连续和引力缺口原样保持。
7. 旧空间382—386、425、522—523和604、649／699直接复用，当前目标不变。
''')
write('round747_drafts/literature_scope_audit.json',json.dumps(dict(sources=[],
inherited='577 original scalar-edge phase transmission;746 original T readout, Gaussian density and interval method;745 exact CAR polynomial implementation.',
new='Exact remote hop-squared sixth coefficient, original complex-spin calibration, strict sign of that contribution.',
excluded='Nonzero total remote signal, absence of scalar-edge cancellation, autonomous detector, geometry generation, full gravity.'),ensure_ascii=False,indent=2)+'\n')
write('round747_drafts/scope_and_dedup_review.md','''# 747范围、去重与解析审查

已回查四份导航、746正式报告及冻结747入口。未发现活跃Python进程；无新代理、应用任务、自动化或图像检查。目标保持。

- 747入口完整64模占据转移不是实际T读数；正式报告保持区别。
- 两跳跃系数的动能／质量计数、CAR奇偶、Dirac模块切换、A反演奇项逐类审查。
- 径向缩约在完成全五坐标微分后执行，利用B规范不变性；不是删除横向动能。
- 原复spin矩阵精确计算，非仅靠Hilbert-Schmidt猜测；不声称任意新边的分类。
- 原标量边和共同几何始终保留；lambda仅分解系数，不能通过调lambda避开总信号问题。
- 证书旧字段含fourth的复用仅在内存中，公开结果已正确命名sixth；未修改旧代码或旧结果。
- 紧支撑逼近只继承本分量符号，不能推出总响应非零或有实用时间窗。
- 关掉Ys的草稿诊断明确是不同模型，只保留探索历史，不据此宣称本模型解析抵消。
- 原实际输出、平均来源和注能界直接复用；认知动机、模型输入、解析系数和数值／区间验证分开。
- 第748轮必须补原同阶标量边；当前不是统一目标结项或完整通信证明。
''')
main=('research_note_747.md','joint_remote_record_response.py','joint_remote_record_response_results.json','unified_physics_condition_ledger_747.md')
write('round747_drafts/final_review.txt','Primary-agent review only. Exact remote hop contribution certified; original scalar-edge total remains open. No independent agent review or full field time simulation claimed.\n'+
'\n'.join(n+' '+hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in main)+'\n')
mapping={'747':'748','746':'747','745':'746','3430':'3432','3428':'3430','1550':'1553',
'3449':'3468','3436':'3449','392':'393','297':'298',
'joint_connected_population_readout':'joint_remote_record_response'}
pub=remap((HERE/'publish_round746.py').read_text('utf8'),mapping)
summary='**第747轮完成（分量结论）：** [原跳跃与异地读口六阶贡献]({p}research_note_747.md)原费米跳跃对旧远端T读口的六阶系数经精确CAR与严格区间认证；原标量边同阶贡献未算，尚未签收总通信。保完整后态、共同几何与来源。两组、十六式通过，最新747／3432，1553份编号科学文件、3468份保护证据。[核验]({p}research_round_747_checks.json)、[条件账]({p}unified_physics_condition_ledger_747.md)。'
order='**当前执行顺序（747后，优先于下方历史安排）：** 接[748原标量边与异地总响应]({p}round748_drafts/STATUS.md)，合并两条原通道的完整时间次序，不改弱边、不插入中间理想测量。旧空间、604、649／699及统一目标保持。'
pub=re.sub(r'^summary=.*$',lambda m:'summary='+repr(summary),pub,flags=re.M)
pub=re.sub(r'^order=.*$',lambda m:'order='+repr(order),pub,flags=re.M)
pub=pub.replace('原连通图中的占据读出','原跳跃与异地读口六阶贡献').replace('旧空间合同保持，原固定图读出与资源合并','旧空间合同保持，异地分量与总响应区分')
# Existence guard and postcheck share the established short completion marker.
pub=pub.replace("assert '**第747轮完成：**' not in text","assert '**第747轮完成（分量结论）：**' not in text")
write('publish_round747.py',pub)
post=remap((HERE/'postcheck_round746.py').read_text('utf8'),mapping).replace("**第747轮完成：**","**第747轮完成（分量结论）：**")
write('postcheck_round747.py',post)
verification=remap((HERE/'verify_round746.py').read_text('utf8'),mapping)
start=verification.index('    names=(');end=verification.index('    new=',start)
names=('unified_physics_condition_ledger_747.md','round747_drafts/research_note_747_draft.md',
'round747_drafts/final_review.txt','round747_drafts/literature_scope_audit.json',
'round747_drafts/scope_and_dedup_review.md','round748_drafts/STATUS.md',
'round747_drafts/remote_population_entry.py','round747_drafts/remote_population_entry_results.json',
'round747_drafts/research_note_747_working.md','round747_drafts/check_and_publish_entry.py',
'round747_drafts/entry_checks.json','round747_drafts/remote_hop_coefficient.py',
'round747_drafts/remote_hop_coefficient_results.json','round747_drafts/remote_hop_sign_certificate.py',
'round747_drafts/remote_hop_sign_certificate_results.json','round747_drafts/majorana_coefficient_diagnostic.json')
verification=verification[:start]+'    names='+repr(names)+'\n'+verification[end:]
verification=verification.replace("checks['display_formulas']==18","checks['display_formulas']==16")
verification=verification.replace('original_fixed_graph_words_and_native_Higgs_readout_sign_checked=True',
'exact_remote_hop_contribution_checked=True,total_remote_signal_proven=False')
verification=verification.replace('    result=model.run();',"""    entry=core.read(HERE/'round747_drafts/entry_checks.json')
    for name,digest in entry['artifact_hashes'].items():
        assert core.digest(HERE/'round747_drafts'/name)==digest,name
    result=model.run();""")
write('verify_round747.py',verification)
print('Prepared747; total remote signal remains open.')
