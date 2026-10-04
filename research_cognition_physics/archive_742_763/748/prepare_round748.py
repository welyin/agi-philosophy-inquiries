"""Prepare748 exact channel-combination identity and return to joint contracts."""
import hashlib,json,re
from pathlib import Path
HERE=Path(__file__).resolve().parent
def write(name,text):
    p=HERE/name;p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)
def remap(text,mapping):
    return re.sub('|'.join(re.escape(k) for k in sorted(mapping,key=len,reverse=True)),lambda m:mapping[m[0]],text)
ledger=(HERE/'unified_physics_condition_ledger_747.md').read_text('utf8')
ledger=ledger.replace('# 联合条件总账：747原跳跃与异地读口六阶贡献','# 联合条件总账：748原标量边与异地总系数',1)
ledger=ledger.replace(ledger.split('\n')[2],'2026-10-04。接[747全账](unified_physics_condition_ledger_747.md)，回填[748报告](research_note_748.md)。[结果](joint_remote_total_coefficient_results.json)、[核验](research_round_748_checks.json)。目标保持，回接共同指认的联合合同。',1)
updates={'C03':'748原T读口两条原通道合成公式已得，不新设中间占据仪器',
'C15':'748总系数由原测地边力与旧局部密度固定；数值同号但严格总符号未认证',
'C19':'748继续使用同一正常Gauss准备，未把声明准备变成内部自主准备',
'C22':'748无独立通信耦合，原w/k/tau同时决定响应；动态来源与跨尺度匹配仍缺'}
rows=ledger.splitlines();seen=[]
for i,row in enumerate(rows):
    for key,value in updates.items():
        if row.startswith('|'+key+' '):
            parts=row.split('|');parts[2]+='；'+value;rows[i]='|'.join(parts);seen.append(key)
assert set(seen)==set(updates)
ledger='\n'.join(rows)+'\n'
ledger=ledger.rsplit('## 本轮合并与下一项',1)[0]+'''## 748原两通道的共同响应及研究优先级

- Duhamel展开与旧P4密度给原标量边六阶响应；完整时间因子为1/720。
- 原曲目标边力显式计算，12组五分量有限差分校准；不把导数检查当积分误差界。
- 与747严格跳跃分量合并。在epsilon=psi_A=psi_B=hbar=1、原复tau下，标量约−0.00690024874，跳跃约−0.11875444084，总和约−0.12565468958。
- 四档求积的收敛仅为数值证据；原总信号的严格符号仍缺。不能将此签成完整异地通信。
- 原同一物质、几何、准备与菜单固定两项，无需另加探针或通信耦合。原独立物理参数并未从认知原则消失。
- 同一后态、资源及来源继续沿已有接口处理，不把原固定几何读出当引力解。
- 依据用户主线提醒，停止继续优化特殊读口的小数位或任意高阶，转回联合总账选择跨部门净约束。当前应用目标保持活跃，不设置新任务或自动化。

## 本轮合并与下一项

本轮关闭的是“原标量边六阶系数是否需独立指定”的局部接口，答案是不需；不是关闭C15全合同。新增总系数公式与数值证据及时记录，未凑轮次重复局部信号。

接[749](round749_drafts/STATUS.md)：回查既有区域组合、共同记录、参考与尺度输送，检验共同指认能否同时约束多个部门。先去重并形成明确候选／失败判据；不把群重命名、有限速度、冗余或图关系直接等同已知物理。旧空间、604、649／699保留。
'''
write('unified_physics_condition_ledger_748.md',ledger)
write('round748_drafts/research_note_748_draft.md',(HERE/'research_note_748.md').read_text('utf8'))
write('round749_drafts/STATUS.md','''# 第749轮入口：共同指认的跨部门合同与去重审计

接[748](../research_note_748.md)、[联合账](../unified_physics_condition_ledger_748.md)及[用户方向提醒](../round748_drafts/joint_scope_steering.md)。

1. 先回查猜想开头与旧530以后联合研究、共同记录624/633/647、区域组合617/649/705/706、记忆与粗化707—709、722—725及来源730—741。以下编号是检索入口，不表示未经核对已完全覆盖问题。
2. 目标是同一量子过程同时承担主体、记录、物质、参考、几何、资源及反作用；共同指认不等于同一量子态描述。
3. 优先提出一个能同时限制若干部门的具体合同，写清事件／读口的对应、主体合并与尺度变换的操作对象、误差与资源及实际来源；标明哪些约束旧研究已证。
4. 若只是重命名旧CP、协变、能量保持或粗化交换图，不计新轮次。必须指出新合同排除何类原候选或减少哪个独立输入。
5. 747—748总六阶数值未获严格积分证书。保留该局部缺口；不无限延伸读口高阶／精度。仅在下一共同检验确需严格通信前提时补粗符号界。
6. 382—386、425、522—523直接复用；604、649、699仍约束共同连续物种、几何过程和特定候选正性。
7. 不改应用目标，不新任务或自动化，不做图像检查。每个真正实质单元写编号报告、可复算结果及账目净变化。
''')
write('round748_drafts/literature_scope_audit.json',json.dumps(dict(sources=[],
inherited='577 native geodesic force;745/746 exact local configuration density;747 exact original hopping coefficient; existing Duhamel expansion derived in note.',
new='Original scalar-edge sixth coefficient and full original channel combination, with declared numerical evidence.',
excluded='Strict total sign, autonomous records, universal readout, continuum or GR derivation.'),ensure_ascii=False,indent=2)+'\n')
write('round748_drafts/scope_and_dedup_review.md','''# 748解析范围与联合贡献审查

回查577、746、747、604、649和猜想开头；已读联合总账。没有重新证明旧空间或已有相位通信，无新代理或图像检查。

- eta及lambda仅是逐词提取标记，实际原H全部边保持。
- Duhamel一边插入分离A配置密度与B交换子。源端低0—3阶差为零，B端首项为t−s；时间因子精确1/720。
- 两势词在六阶没有足够动能让A端两个质量和B读口同时贡献；混合势与两跳跃同样失败。原现场势并未删出H。
- 角平均来自Gauss不变准备及Haar链路；S3角分布不是物理空间维数声明。
- 原曲目标q与force在q=1有光滑极限；五分量差分仅校准表达式。
- Laguerre/Hermite/相对角加密是经验收敛，没有报告为严格符号证书。
- 完整总和负号尚未严格认证，紧支撑实际信号只写条件性推论。
- 联合净贡献是移除该响应中的独立通信系数，不是减少所有原物理输入。完整C15、自治参考、动态几何和共同尺度仍开放。
- 用户方向提醒已保存；接回总账，避免继续特殊读口精度优化。未修改应用目标。
''')
main=('research_note_748.md','joint_remote_total_coefficient.py','joint_remote_total_coefficient_results.json','unified_physics_condition_ledger_748.md')
write('round748_drafts/final_review.txt','Primary-agent review only. Exact scalar-channel reduction and joint coefficient; numerical evidence is not a strict sign certificate. Joint-condition synthesis next.\n'+
'\n'.join(n+' '+hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in main)+'\n')
mapping={'748':'749','747':'748','746':'747','3432':'3434','3430':'3432','1553':'1556',
'3468':'3480','3449':'3468','393':'394','298':'299',
'joint_remote_record_response':'joint_remote_total_coefficient'}
pub=remap((HERE/'publish_round747.py').read_text('utf8'),mapping)
summary='**第748轮完成（合成公式）：** [原标量边与异地总系数]({p}research_note_748.md)旧局部密度与原测地边力固定标量通道，和原跳跃共同给总响应；校准实例数值同号约−0.12565，严格总符号尚缺。无需独立通信耦合。两组、十二式通过，最新748／3434，1556份编号科学文件、3480份保护证据。[核验]({p}research_round_748_checks.json)、[条件账]({p}unified_physics_condition_ledger_748.md)。'
order='**当前执行顺序（748后，优先于下方历史安排）：** 接[749共同指认与跨部门合同]({p}round749_drafts/STATUS.md)，先复用区域／尺度／来源旧结果，寻找联合净约束；保留局部总信号证书缺口，不继续特殊读口精度优化。统一目标及旧空间、604、649／699保持。'
pub=re.sub(r'^summary=.*$',lambda m:'summary='+repr(summary),pub,flags=re.M)
pub=re.sub(r'^order=.*$',lambda m:'order='+repr(order),pub,flags=re.M)
pub=pub.replace('原跳跃与异地读口六阶贡献','原标量边与异地总系数').replace('旧空间合同保持，异地分量与总响应区分','旧空间合同保持，原总响应接回联合合同')
pub=pub.replace('**第748轮完成（分量结论）：**','**第748轮完成（合成公式）：**')
write('publish_round748.py',pub)
post=remap((HERE/'postcheck_round747.py').read_text('utf8'),mapping).replace('**第748轮完成（分量结论）：**','**第748轮完成（合成公式）：**')
write('postcheck_round748.py',post)
verification=remap((HERE/'verify_round747.py').read_text('utf8'),mapping)
start=verification.index("    entry=core.read(");end=verification.index('    result=model.run();',start)
verification=verification[:start]+verification[end:]
start=verification.index('    names=(');end=verification.index('    new=',start)
names=('unified_physics_condition_ledger_748.md','round748_drafts/research_note_748_draft.md',
'round748_drafts/final_review.txt','round748_drafts/literature_scope_audit.json',
'round748_drafts/scope_and_dedup_review.md','round749_drafts/STATUS.md',
'round748_drafts/scalar_edge_response.py','round748_drafts/scalar_edge_response_results.json',
'round748_drafts/joint_scope_steering.md')
verification=verification[:start]+'    names='+repr(names)+'\n'+verification[end:]
verification=verification.replace("checks['display_formulas']==16","checks['display_formulas']==12")
verification=verification.replace('exact_remote_hop_contribution_checked=True','original_scalar_force_and_total_coefficient_formula_checked=True')
write('verify_round748.py',verification)
print('Prepared748; exact formula, numerical total, joint-ledger next.')
