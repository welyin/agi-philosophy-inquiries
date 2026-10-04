"""Prepare744 publication without modifying frozen history or entry files."""
import hashlib,json,re
from pathlib import Path
HERE=Path(__file__).resolve().parent
def write(name,text):
    p=HERE/name;p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)
def remap(text,mapping):
    return re.sub('|'.join(re.escape(k) for k in sorted(mapping,key=len,reverse=True)),lambda m:mapping[m[0]],text)

ledger=(HERE/'unified_physics_condition_ledger_743.md').read_text('utf8')
ledger=ledger.replace('# 联合条件总账：743原量子标量与严格Gauss关联','# 联合条件总账：744原读口与联合历史',1)
first=ledger.split('\n')[2]
ledger=ledger.replace(first,'2026-10-04。接[743全账](unified_physics_condition_ledger_743.md)，回填[744报告](research_note_744.md)。[结果](joint_native_history_readout_results.json)、[核验](research_round_744_checks.json)。目标及当前任务保持不变。',1)
updates={'C03':'744原未知态仪器保持CP及完整物理输出；未实现理想占据仪器',
         'C15':'744单读占据盲性解析成立，偶联合历史有原单节点数值信号',
         'C19':'744全局反射对称准备给精确选择规则；新Gaussian准备与743及730不同',
         'C22':'744全部历史后态共用原来源；尚未认证联合信号和有限时反作用'}
rows=ledger.splitlines();seen=[]
for i,row in enumerate(rows):
    for key,value in updates.items():
        if row.startswith('|'+key+' '):
            parts=row.split('|');parts[2]+='；'+value;rows[i]='|'.join(parts);seen.append(key)
assert set(seen)==set(updates)
ledger='\n'.join(rows)+'\n'
ledger=ledger.rsplit('## 本轮合并与下一项',1)[0]+'''## 744原读口的联合历史

- 原完整有限图H保全局singlet反射与费米四分之一相位的联合对称；对称准备下单次sin s效果任意时刻不分编码空态与配对态。
- 保完整未知态CP输出；原准备和最终读口仍为输入，不以编码效果替代实际后态。
- 两次原平方根读取的结果乘积为偶报告，允许占据信息。概率对应sin²(s)/4，但后态不可替换成合并效果的平方根。
- 原一节点完整H5量子动能、势及32CAR质量，另用正常Gaussian准备，给六阶差约−0.55005077145的数值证据。
- 该信号不是严格符号证书、有限时结果或连通图结论；低阶消去也未全部解析证明。
- 来源、能量与反作用仍归属原完整历史；旧空间、604及649／699保持。

## 本轮合并与下一项

C03／C15／C19／C22共用原过程，单次限制不扩大为全历史不可读。

接[745](round745_drafts/STATUS.md)：严格认证原联合历史信号，或给出导致假信号的反例。优先误差来源及有限时间接口，不靠无限提高求积阶数替代证明；再接正间隔和连通图。
'''
write('unified_physics_condition_ledger_744.md',ledger)
write('round744_drafts/research_note_744_draft.md',(HERE/'research_note_744.md').read_text('utf8'))
write('round745_drafts/STATUS.md','''# 第745轮入口：原联合历史的严格信号与有限时接口

接[744](../research_note_744.md)、[全账](../unified_physics_condition_ledger_744.md)。

1. 保原读口、完整物理后态及同一实际H；单次选择规则已严格完成，不重复证明。
2. 先核单节点原Gaussian准备的六阶候选信号：稀疏多项式误差、实际积分权重、解析消去或有界近似。数值收敛不是严格符号证明。
3. 若能认证非零导数，再连接有限时间存在或明确余项；不误把六阶系数当已证的首项。若失败，记录真正障碍及可检验替代。
4. 随后处理两次读取之间正等待及连通图。原多物种、原势、五维内部目标、所有反作用继续保留。
5. 原准备和最终仪器仍是操作输入，来源域复用623—625、704、718。不是自治探测器、理想n_f读取或动态引力完成。
6. 不更改研究目标，不新增物种、几何或认知公理来迎合信号；旧空间382—386、425、522—523及604、649／699保持。
''')
write('round744_drafts/literature_scope_audit.json',json.dumps(dict(
    sources=[dict(url='https://arxiv.org/html/1810.06512v3',locations='Sections3.2,3.3',
        use='Induced effects and successive instruments as background; does not prove our native signal or autonomous detector.')],
    inherited='598 complete H;623 domain;624-625 histories;717 native force;718 full cq and source treatment;743 nonGaussian preparation.',
    new='Exact full-H single-read selection rule; actual parity-history instrument; sixth-jet evidence using native onsite H.',
    numerical_scope='One vertex, a newly declared Gaussian ready state, full five-coordinate differentiation and original32CAR; not interval certified.',
    excluded='A complete occupation instrument, finite-time discrimination certificate, connected-graph theorem, autonomous terminal, continuum or gravity.'),ensure_ascii=False,indent=2)+'\n')
write('round744_drafts/scope_and_dedup_review.md','''# 744范围与证明审查

上一轮743完成、744入口执行，属于进展。当前无活跃Python进程，无新代理、任务、图像检查或目标修改。

- 717读口力、718真实历史与623共同域复用；新的是原全局联合对称性对未知态实际效果的约束。
- 全局反射和CAR四分之一相位各反转Majorana，合起来保H；不能替换成逐节点反射或连续反常结论。
- ΘV=VZ的准备条件必须保留。单次只约束效果，未声称完整条件后态相同。
- 两次读取保L_q L_r；合并概率不等于更换为一个平方根仪器。
- 新Gaussian准备和一节点图是明确输入。全五方向微分后才能取Higgs径向代表。
- 欧氏共轭Laplacian是α²|x|²；H5部分为(Dtilde²+4Dtilde)/6，与代码一致。
- 六阶期望公式直接来自两个传播子的弱导数。所有H矩有限不代表时间Taylor级数收敛。
- 0至4阶数值抵消不替代解析消去，第五阶未纳入该检查，未将六阶称为首个非零阶。
- 浮点、系数舍弃和求积未区间化；三阶数收敛与恒等残差只支持数值证据。
- 全部原质量保留，没有CAR截断；无边是单节点图本身。未外推到连通图和连续Hadamard参考。
- 共同来源使用完整历史后态，原终端仪器仍需实现；不把CP存在式称自治测量或引力生成。
''')
main=('research_note_744.md','joint_native_history_readout.py','joint_native_history_readout_results.json','unified_physics_condition_ledger_744.md')
write('round744_drafts/final_review.txt','Primary-agent review only. Exact single-read symmetry; parity-history signal is native onsite numerical evidence, not interval certified.\n'+
      '\n'.join(n+' '+hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in main)+'\n')
mapping={'744':'745','743':'744','742':'743','3424':'3426','3422':'3424','1541':'1544',
         '3401':'3421','3387':'3401','389':'390','294':'295',
         'joint_native_quantum_record':'joint_native_history_readout','quantum_scalar_record_entry':'native_instrument_selection_entry'}
publication=remap((HERE/'publish_round743.py').read_text('utf8'),mapping)
summary='**第744轮完成：** [原读口与联合历史]({p}research_note_744.md)原全H的联合对称性严格限制单次读取；两次原读取的偶报告允许占据信息。原一节点完整量子H给六阶差约−0.55005的数值证据，严格符号、有限时及连通图结论仍待证。两组、十六式通过，最新744／3426，1544份编号科学文件、3421份保护证据。[核验]({p}research_round_744_checks.json)、[条件账]({p}unified_physics_condition_ledger_744.md)。'
order='**当前执行顺序（744后，优先于下方历史安排）：** 接[745联合历史的严格信号]({p}round745_drafts/STATUS.md)，核原数值信号的误差来源和有限时接口；随后接正等待及连通图，保完整后态与来源。旧空间、604、649／699及统一目标保持。'
publication=re.sub(r'^summary=.*$',lambda m:'summary='+repr(summary),publication,flags=re.M)
publication=re.sub(r'^order=.*$',lambda m:'order='+repr(order),publication,flags=re.M)
publication=publication.replace('原量子标量与严格Gauss记录关联','原读口与联合历史')
publication=publication.replace('旧空间合同保持，原相互作用与实际记录共同核验','旧空间合同保持，单次与联合读取范围分开')
write('publish_round744.py',publication)
write('postcheck_round744.py',remap((HERE/'postcheck_round743.py').read_text('utf8'),mapping))
verification=remap((HERE/'verify_round743.py').read_text('utf8'),mapping)
verification=verification.replace("checks['display_formulas']==18","checks['display_formulas']==16")
verification=verification.replace("           'round744_drafts/entry_checks.json')","           'round744_drafts/entry_checks.json',\n"+
    "           'round744_drafts/history_jet_probe.py','round744_drafts/history_jet_probe_results.json',\n"+
    "           'round744_drafts/history_sixth_jet_probe.py','round744_drafts/history_sixth_jet_probe_results.json',\n"+
    "           'round744_drafts/history_density_probe.py','round744_drafts/history_density_probe_results.json')")
verification=verification.replace('original_full_H_Gauss_preparation_four_point_and_scalar_correlation_checked=True',
    'original_full_H_selection_and_native_onsite_history_jet_checked=True')
write('verify_round744.py',verification)
print('Prepared744 report and745 strict-signal entry.')
