"""Publish-ready746 materials: fixed original graph and existing Higgs readout."""
import hashlib,json,re
from pathlib import Path
HERE=Path(__file__).resolve().parent
def write(name,text):
    p=HERE/name;p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)
def remap(text,mapping):
    return re.sub('|'.join(re.escape(k) for k in sorted(mapping,key=len,reverse=True)),lambda m:mapping[m[0]],text)
ledger=(HERE/'unified_physics_condition_ledger_745.md').read_text('utf8')
ledger=ledger.replace('# 联合条件总账：745原联合历史的严格信号','# 联合条件总账：746原连通图中的占据读出',1)
first=ledger.split('\n')[2]
ledger=ledger.replace(first,'2026-10-04。接[745全账](unified_physics_condition_ledger_745.md)，回填[746报告](research_note_746.md)。[结果](joint_connected_population_readout_results.json)、[核验](research_round_746_checks.json)。目标及当前任务保持不变。',1)
updates={'C03':'746旧T菜单给原固定图未知偶编码的实际对角效果，保完整CP后态；不是理想占据投影',
         'C15':'746原固定边的本地占据读出已证，异地记录尚待接通',
         'C19':'746反射对称紧支撑Gauss准备同初始平均能源；未实现自治准备',
         'C22':'746同源注能不超过hbar²/(9w_v)，全部原跳跃与几何保留'}
rows=ledger.splitlines();seen=[]
for i,row in enumerate(rows):
    for key,value in updates.items():
        if row.startswith('|'+key+' '):
            parts=row.split('|');parts[2]+='；'+value;rows[i]='|'.join(parts);seen.append(key)
assert set(seen)==set(updates)
ledger='\n'.join(rows)+'\n'
ledger=ledger.rsplit('## 本轮合并与下一项',1)[0]+'''## 746原固定连通图的本地读取

- 原共形几何使k_vw=(w_v^(1/6)+w_w^(1/6))²/4，不能固定体积而免费缩弱所有边。未采用该比较代替共同几何。
- 改用624／723已经存在的T=|X|²/2平方根菜单；没有新增物种、读口或耦合。744／745的sin s结论保持其原范围。
- 原保持物种模块的规范协变跳跃、所有标量图势、完整CAR及正几何继续参与演化。
- 四阶占据差的算符词精确局部化：图势不进入，跨节点单跳远端奇偶杀掉，两跳加一个Dirac不能回到原物种对角。
- 原局部Gaussian比较给严格J_T约−0.15145884204633；正常反射偶紧支撑逼近把符号接入全图核准备，前3阶差为零。
- 原未知编码输入的T效果严格对角，充分小正时间可读占据；这是本地读取，不是跨节点通信或无扰测量。
- 同准备所有编码态初始平均能源相同；单读注能至多hbar²/(9w_v)，来源继续使用完整后态。
- 未给可执行时间／截止半径预算，未证明自治终端、连续或动态量子引力。旧空间、604、649／699保持。

## 本轮合并与下一项

C03／C15／C19／C22在原固定边图上合并本地可读与资源，不增加弱边假设。一般空间维数和连续缺口不重新研究。

接[747](round747_drafts/STATUS.md)：将未知占据信息沿原相互作用送入其他节点的实际记录，复用577而不把已有相位通信重复计轮次；控制、消息及来源须共用原过程。
'''
write('unified_physics_condition_ledger_746.md',ledger)
write('round746_drafts/research_note_746_draft.md',(HERE/'research_note_746.md').read_text('utf8'))
write('round747_drafts/STATUS.md','''# 第747轮入口：原未知占据输入与异地记录

接[746](../research_note_746.md)、[条件账](../unified_physics_condition_ledger_746.md)。

1. 746已在原固定边图上证明本地T读口可读未知偶编码占据；不能把它当成其他主体已经收到记录。
2. 先回查577（原完整图的A相位→B读口三阶通道）、598、623—625、718。不重复证明已有相位通信。
3. 输入继续为原sterile空／对未知编码。异地原T或s菜单、原完整量子作用和全部中间后态保留；不注入无成本条件控制或独立经典信道。
4. 区分原自治演化实际转导与先测后按结果另作准备的操作流程；若引入控制，说明来源及是否原模型已有。
5. 优先确定最早允许的跨节点响应、原边势／跳跃贡献及真正非零准备。不能只用一般连续性、互信息名称或局部系数作通信证明。
6. 共同资源、原几何和来源同时保留。若当前接口可直接复用旧结论而无新增，则回到统一条件总账寻找真实缺口，不虚增轮次。
7. 旧空间382—386、425、522—523及604、649／699保持，目标不变。
''')
write('round746_drafts/literature_scope_audit.json',json.dumps(dict(
    sources=[],inherited='574/589 shared geometry;598 original CAR, hopping and readout energy;623 common domain;624/723 existing T instrument;717 mass force;745 exact polynomial and certified-quadrature method.',
    new='Original fixed-graph fourth-order locality of occupation contrast for the existing Higgs readout; strict local sign connected to compact physical preparations and same energy/source.',
    not_repeated='577 phase transmission;723 incomplete single-reference marginals;744 reflection identity;generic CP or Cauchy theorem.',
    excluded='Remote communication, ideal number instrument, weak-edge model substitution, autonomous terminal, quantum continuum and gravity.'),ensure_ascii=False,indent=2)+'\n')
write('round746_drafts/scope_and_dedup_review.md','''# 746范围与解析审查

上一目标轮744／745完成并发布，属于有效进展。已读最新导航、笔记与结果，未发现活跃Python进程。无新代理、任务、自动化或图像检查。

- 对共同几何回查揭示固定w独立缩小k不合法；本轮未把数学弱耦合比较冒充原几何。
- 旧T菜单在624／723已经存在。更换当前使用的菜单项明确写出，没有说原sin s联合协议已推广全图。
- 原H写作标量动能、标量全部图势和矩阵质量／跳跃。图势可非局部于端点，但仍是Fock标量乘法。
- 四阶零B项CAR差为零；一B项始终无编码对角。二B项中的额外标量势插入逐词为零。
- 首个B插入只能在v；另一节点偶质量对易，单跳在邻端留奇CAR次数。对系数的玻色导数不改变该结构。
- 三B词只剩对局部Dirac力做二次矩阵交换子。奇Majorana改变粒子数；奇Dirac改变左右模块；一跳远端奇、两跳保模块，均不给占据对角。
- 低三阶也逐项检查。全局Θ使T效果对角，但不使实际仪器无扰。
- 精确P4从745原五维／32CAR结果复用，当前sinT角积分及误差重新计算，不把Gaussian宽度扫描当新定理。
- 全图存在性使用反射偶紧支撑核态，避免直接声称一节点Gaussian具有原图所有高阶矩。局部有限阶算符的多项式界使截止系数趋于严格负Gaussian值。
- 注能按原H计，V†HV为同一标量，条件后态也有限能源。准备和终端仍是输入。
- 输入与输出同节点；没有跨节点信号、图统一时间窗、动态几何或完整量子理论唯一性结论。
''')
main=('research_note_746.md','joint_connected_population_readout.py','joint_connected_population_readout_results.json','unified_physics_condition_ledger_746.md')
write('round746_drafts/final_review.txt','Primary-agent review only. Original fixed graph retains local population readability through the existing T menu; analytic word proof and strict interval. Remote signaling and autonomous detector remain open.\n'+
      '\n'.join(n+' '+hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in main)+'\n')
mapping={'746':'747','745':'746','744':'745','3428':'3430','3426':'3428','1547':'1550',
         '3436':'3449','3421':'3436','391':'392','296':'297',
         'joint_native_history_certificate':'joint_connected_population_readout'}
publication=remap((HERE/'publish_round745.py').read_text('utf8'),mapping)
summary='**第746轮完成：** [原连通图中的占据读出]({p}research_note_746.md)旧Higgs幅度读口在原固定边图上有严格四阶占据差，图势和跳跃仅在此阶差值中抵消；全演化、Gauss及同态资源保留。无需弱化边或新增探针，尚非异地通信、无扰测量或自治终端。两组、十八式通过，最新746／3430，1550份编号科学文件、3449份保护证据。[核验]({p}research_round_746_checks.json)、[条件账]({p}unified_physics_condition_ledger_746.md)。'
order='**当前执行顺序（746后，优先于下方历史安排）：** 接[747原占据输入与异地记录]({p}round747_drafts/STATUS.md)，复用577通道并保未知态、实际后态及共同代价，避免新增免费控制。旧空间、604、649／699及统一目标保持。'
publication=re.sub(r'^summary=.*$',lambda m:'summary='+repr(summary),publication,flags=re.M)
publication=re.sub(r'^order=.*$',lambda m:'order='+repr(order),publication,flags=re.M)
publication=publication.replace('原联合历史的严格信号','原连通图中的占据读出')
publication=publication.replace('旧空间合同保持，严格原历史信号进入图接口','旧空间合同保持，原固定图读出与资源合并')
write('publish_round746.py',publication)
write('postcheck_round746.py',remap((HERE/'postcheck_round745.py').read_text('utf8'),mapping))
verification=remap((HERE/'verify_round745.py').read_text('utf8'),mapping)
start=verification.index("    names=(");end=verification.index("    new=",start)
verification=verification[:start]+"""    names=('unified_physics_condition_ledger_746.md','round746_drafts/research_note_746_draft.md',
           'round746_drafts/final_review.txt','round746_drafts/literature_scope_audit.json',
           'round746_drafts/scope_and_dedup_review.md','round747_drafts/STATUS.md',
           'round746_drafts/local_population_jet_probe.py','round746_drafts/local_population_jet_probe_results.json',
           'round746_drafts/higgs_readout_certificate.py','round746_drafts/higgs_readout_certificate_results.json')
"""+verification[end:]
verification=verification.replace("checks['display_formulas']==16","checks['display_formulas']==18")
verification=verification.replace('exact_native_coefficients_and_strict_integral_sign_checked=True',
    'original_fixed_graph_words_and_native_Higgs_readout_sign_checked=True')
write('verify_round746.py',verification)
print('Prepared746 fixed-graph result and747 remote-record interface.')
