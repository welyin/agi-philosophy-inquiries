"""Guarded addition of completed round 460 and the user's particle hypothesis."""
from pathlib import Path

HERE=Path(__file__).resolve().parent
BASE,ROOT=HERE.parent,HERE.parent.parent
paths=[ROOT/'README.md',BASE/'README.md',BASE/'research_direction.md',
       BASE/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md']
raw={p:p.read_bytes() for p in paths}
texts={p:v.decode('utf-8').replace('\r\n','\n') for p,v in raw.items()}
assert '**第460轮完成：**' not in texts[BASE/'RESEARCH_STATE.md']
for p in paths[:5]:
    prefix='research_cognition_physics/archive_231_/' if p==ROOT/'README.md' else ('' if p==HERE/'README.md' else 'archive_231_/')
    quote='> ' if p==ROOT/'README.md' else ''
    old=next(line for line in texts[p].splitlines() if '**第459轮完成：**' in line)
    new=quote+'**第460轮完成：** [固定接口中的记忆响应最低作用权限](%sresearch_note_460.md)完整分类四载体全部SU(2)不变至多三体项；精确保旧三体码时，私有L驱动外部响应的最低支撑为3，若另要求原始自旋时间反演偶性则为4。两种达到例均给未知G／参考一致的实际读数与有限窗；保持任务人口不等于不扰动完整L边缘。6项检查，累计2193项；691份编号科学文件、727份保护证据。高体项是额外输入，未由429自行生成。'%prefix
    particle=quote+'**用户新增方向：组织的稳定模式与粒子。** [猜想及验收记录](%srecursive_organization_particle_hypothesis.md)接续原物理生成纲领，检验同一底层动力学下的稳定、可传播模式是否具有粒子性质；内部SU(2)与组织层级尚不等于物理自旋、标准模型或空间维数。已核Levin–Wen现成路线作为后继参照，不新增科学轮次，不改变当前先落实三维空间的阶段目标。'%prefix
    texts[p]=texts[p].replace(old,old+'\n\n'+new+'\n\n'+particle,1).replace('231—459轮','231—460轮')
p=BASE/'research_direction.md'
texts[p]=texts[p].replace('最新科学轮次与检查数为459／2187','最新科学轮次与检查数为460／2193')
p=BASE/'RESEARCH_STATE.md'
texts[p]=texts[p].replace('最新459轮及累计2187项见本文开头','最新460轮及累计2193项见本文开头')
p=BASE/'README.md'
texts[p]=texts[p].replace('当前复算与冻结入口：','当前复算与冻结入口：[460科学核验](archive_231_/verify_protected_memory_actuation_round.py)、[460整合核验](archive_231_/verify_round460_integration.py)；历史入口：',1)
texts[p]=texts[p].replace('第三阶段累计2187项；旧科学证据保持原字节','第三阶段累计2193项；旧科学证据保持原字节')
p=HERE/'README.md'
texts[p]=texts[p].replace('当前完成459轮','当前完成460轮').replace('当前459不作为预定终点','当前460不作为预定终点')
row=next(line for line in texts[p].splitlines() if line.startswith('| [459：'))
texts[p]=texts[p].replace(row,row+'\n| [460：固定接口中的记忆响应](research_note_460.md) | 全部三局域作用分类；最低支撑3／T偶时4；实际响应与记忆扰动 | [代码](protected_memory_actuation_audit.py)、[结果](protected_memory_actuation_audit_results.json)、[核验](research_round_460_checks.json)；6项检查；科学基线459 |',1)
texts[p]=texts[p].replace('当前整合入口：','当前整合入口：[460整合核验](verify_round460_integration.py)；历史入口：',1)
texts[p]=texts[p].replace('当前2187项科学检查；724份保护证据，其中本阶段编号科学文件688份','当前2193项科学检查；727份保护证据，其中本阶段编号科学文件691份')
p=HERE/'spatial_premise_closure_audit.md'
assert '## 88.' not in texts[p]
texts[p]+='''
## 88. 第460轮：固定旧接口的最小解除权限

2026-09-24，基线459。三个原始载体A及一个对象B、无额外作用媒介、共同SU(2)不变、全部时间精确保住旧A自旋1／2代码。Pauli不变张量给全部至多三体作用的11维实Hermitian基：I、六SWAP、四手征χ＝i[Sij,Sjk]。对3P_A的跨边和跨手征整数Gram各秩2，核分别(1,1,1)及(1,−1,1)，实／纯虚部门无抵消；给完整保码充要条件。

均匀K仍只作用G。C＝χ013−χ023＋χ123则精确保码，压缩为sqrt(3) Y_L σ_G·σ_B。H3＝K＋C／(2sqrt(3))让Y−人口不作用目标、Y+人口驱动部分交换，初始对象及内部参考为singlet时，实际效果为I_G⊗[I−3sin²(2t)ΠY+／4]，覆盖全部未知GLR。中心π／4给差3／4，半宽1／8窗差≥45／64。

原始自旋反幺正Θ＝(iσy)^⊗4𝒦使全部χ为T奇。若额外要求T偶，全部至多三体不变作用退回两体，仍不能输出新L。H4＝(I−S01)(I＋S23)／2为T偶四体达到例，ZZZZ系数−1／8，精确保码并由L0驱动；中心π／2、半宽1／4给同一下界。因此最低原始支撑为3，另加T偶则4；不是摄动阶数或空间维数。

两例均只保持所用记录人口，非全部未知量子边缘。H3中心使指定纯L的输出本征值3／4、1／4，相干转入整体关联，完整等距仍保留未知信息。三／四体作用、指定系数、接触及singlet资源都是额外建模输入；434—436一般模拟不等于自然来源。此下界不排除媒介、近似或改变旧编码的两体路线。

6项检查、13公式通过，独立审查核完整性、Gram、交织、实际效果与时间反演范围；累计2193项、691份编号科学文件、727份保护证据。旧文件字节不变。

用户同期猜想大组织的稳定模式或对应标准模型基本粒子，已另存recursive_organization_particle_hypothesis.md，接续原生成纲领并核Levin–Wen来源。它是待检验方向：同一规则中的稳定传播、物理自旋、统计、荷和相互作用仍需建立，内部表示标签不能直接命名电子。未增加认知公理或科学轮次，保留先落实三维空间的阶段目标。
'''
for p in paths:
    assert p.read_bytes()==raw[p],f'Concurrent edit: {p}'
for p,s in texts.items():
    assert p.read_bytes()==raw[p],f'Concurrent edit: {p}'
    newline='\r\n' if b'\r\n' in raw[p] else '\n'
    p.write_bytes((s.rstrip()+'\n').replace('\n',newline).encode('utf-8'))
