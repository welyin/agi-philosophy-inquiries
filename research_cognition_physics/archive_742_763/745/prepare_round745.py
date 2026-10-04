"""Prepare745 certified native history report and next connected-graph interface."""
import hashlib,json,re
from pathlib import Path
HERE=Path(__file__).resolve().parent
def write(name,text):
    p=HERE/name;p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)
def remap(text,mapping):
    return re.sub('|'.join(re.escape(k) for k in sorted(mapping,key=len,reverse=True)),lambda m:mapping[m[0]],text)
ledger=(HERE/'unified_physics_condition_ledger_744.md').read_text('utf8')
ledger=ledger.replace('# 联合条件总账：744原读口与联合历史','# 联合条件总账：745原联合历史的严格信号',1)
first=ledger.split('\n')[2]
ledger=ledger.replace(first,'2026-10-04。接[744全账](unified_physics_condition_ledger_744.md)，回填[745报告](research_note_745.md)。[结果](joint_native_history_certificate_results.json)、[核验](research_round_745_checks.json)。目标及当前任务保持不变。',1)
updates={'C03':'745保原实际两次CP输出；正读取间隔的可区分存在性已证，理想仪器仍开放',
         'C15':'745同一单节点原准备的联合占据信号经精确系数及区间积分认证',
         'C19':'745冻结十进制参数精确解释明确，未认证RG无限精度或730参考',
         'C22':'745同一原H及后态保持；具体有限时资源与连通图来源仍须核'}
rows=ledger.splitlines();seen=[]
for i,row in enumerate(rows):
    for key,value in updates.items():
        if row.startswith('|'+key+' '):
            parts=row.split('|');parts[2]+='；'+value;rows[i]='|'.join(parts);seen.append(key)
assert set(seen)==set(updates)
ledger='\n'.join(rows)+'\n'
ledger=ledger.rsplit('## 本轮合并与下一项',1)[0]+'''## 745联合历史信号的严格认证

- 原冻结表的十进制参数及Yukawa数作为有理实数，精确求u；未把RG数值表当成已认证的物理常数。
- 独立整数CAR和全五维微分给精确六阶差密度；不舍弃小系数，不冻结玻色场。
- 一维角积分化简、向外Decimal区间及解析求积／sin级数／尾部界共同给未归一化差J在约−0.41431411至−0.41431407之间，严格负。
- 六阶差严格非零推出任意短邻域内存在可区分等待；强连续给足够小正读取间隔保信号。未指定可用时间窗口。
- 这是原单节点正常准备的存在性，未扩为给定连通图、理想占据仪器、末读自治或连续／引力完成。
- 744单次全时间盲性继续成立；旧空间与604、649／699保持。

## 本轮合并与下一项

C03／C15／C19／C22在同一原菜单上接通严格可读存在性。解除单次限制无需立即新增物种。

接[746](round746_drafts/STATUS.md)：原连通图、完整跳跃、共同准备与有限资源；先区分弱边耦合存在性和原固定系数实例，停止一节点求积精度优化。
'''
write('unified_physics_condition_ledger_745.md',ledger)
write('round745_drafts/research_note_745_draft.md',(HERE/'research_note_745.md').read_text('utf8'))
write('round746_drafts/STATUS.md','''# 第746轮入口：严格联合信号与原连通图

接[745](../research_note_745.md)、[条件账](../unified_physics_condition_ledger_745.md)。

1. 原一节点联合占据信号已严格认证，不再做积分精度优化。单次盲性仍保。
2. 回查598、623—625、704和718：完整有限图、规范链路、Gauss、所有Dirac／Majorana及跳跃必须同时保留。
3. 从去耦参考接连通图时，明确哪些边系数可作为既有模型参数变化；弱耦合存在开集不是原给定边强度的证明。
4. 同一正常准备、实际正读取间隔、未知输入及完整cq后态必须共同运输。保有限能源与同源代价，不能重新注入无成本准备。
5. 若能证明小边耦合的信号稳定，再判断原实例是否可用同一证书；只取得存在结果时按该范围报告。
6. 终端读口仍是操作输入；连续、动态几何及旧649／699边界保持。涉及空间时复用382—386、425、522—523。
''')
write('round745_drafts/literature_scope_audit.json',json.dumps(dict(
    sources=[dict(url='https://docs.python.org/3/library/decimal.html',use='Correctly rounded arithmetic/exp/sqrt and neighboring decimal numbers; no general noninteger power used.')],
    inherited='623 domain and744 native selection/history objects; standard Cauchy estimate and finite differences proved in their needed form.',
    new='Computer-assisted strict sign in the original onsite instance with exact-decimal inputs; positive waits and positive read-gap existence.',
    excluded='Explicit useful time window, connected fixed graph, autonomous detector, continuous quantum SM, gravity completion.'),ensure_ascii=False,indent=2)+'\n')
write('round745_drafts/scope_and_dedup_review.md','''# 745范围与证书审查

744已正式发布，属于进展。未修改其冻结文件，无新代理或图像检查。

- 回查723：它处理四参考边缘相同但不同顺序记录，不是本轮未知sterile空／对输入及原sin s的时间响应；不重复其一般历史结论。
- 精确参数来自原冻结表十进制及原Yukawa，不声称表背后的RG误差已证明。u由有理数解出，与旧二进制计算的差仅作校准。
- 原质量块以整数重新构建并使用精确CAR符号；全部五维动能在微分之后才用径向内积。
- H幂分母D^k，六阶内积共同D^6。所有小系数保留；26项精确差密度已保存。
- 0、1、2、3、5阶密度为零；4阶密度非零，不误宣布六阶是首项。
- 角Jacobian为R^4(1-z²)，原度量权重未漏。符号判断不需要归一化值，因为Z严格正。
- cos实Taylor余项、全半轴Gamma矩、R≥10真函数尾部共同有显式界。
- 复圆盘不跨6+z²的零点，Cauchy界与插值20次矩身份实际核对；负权重绝对值保留。
- Decimal每项向外扩展，输入分数也包围；只用基本运算、sqrt、exp，不使用一般power。
- 六阶负号只推出任意短邻域中存在非零时刻；正读间隔由强连续给存在，不提供可执行时间值。
- 未证明连通图或完整自治。原势乘0、2的辅助诊断明确不是原模型且不参与严格结论。
''')
main=('research_note_745.md','joint_native_history_certificate.py','joint_native_history_certificate_results.json','unified_physics_condition_ledger_745.md')
write('round745_drafts/final_review.txt','Primary-agent review only. Exact-decimal native onsite history has a certified negative sixth contrast derivative. Positive-time existence only; no useful window or connected-graph completion.\n'+
      '\n'.join(n+' '+hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in main)+'\n')
mapping={'745':'746','744':'745','743':'744','3426':'3428','3424':'3426','1544':'1547',
         '3421':'3436','3401':'3421','390':'391','295':'296',
         'joint_native_history_readout':'joint_native_history_certificate'}
publication=remap((HERE/'publish_round744.py').read_text('utf8'),mapping)
summary='**第745轮完成：** [原联合历史的严格信号]({p}research_note_745.md)原一节点完整H的联合占据差经精确系数和带解析余项的区间积分认证；任意短邻域内存在可区分等待，足够小正读取间隔保差异。未给可用时间窗口或连通图结果。两组、十六式通过，最新745／3428，1547份编号科学文件、3436份保护证据。[核验]({p}research_round_745_checks.json)、[条件账]({p}unified_physics_condition_ledger_745.md)。'
order='**当前执行顺序（745后，优先于下方历史安排）：** 接[746原连通图与共同过程]({p}round746_drafts/STATUS.md)，保完整跳跃、Gauss、未知输入及资源，区分弱耦合存在性和原固定边系数。旧空间、604、649／699及统一目标保持。'
publication=re.sub(r'^summary=.*$',lambda m:'summary='+repr(summary),publication,flags=re.M)
publication=re.sub(r'^order=.*$',lambda m:'order='+repr(order),publication,flags=re.M)
publication=publication.replace('原读口与联合历史','原联合历史的严格信号')
publication=publication.replace('旧空间合同保持，单次与联合读取范围分开','旧空间合同保持，严格原历史信号进入图接口')
write('publish_round745.py',publication)
write('postcheck_round745.py',remap((HERE/'postcheck_round744.py').read_text('utf8'),mapping))
verification=remap((HERE/'verify_round744.py').read_text('utf8'),mapping)
start=verification.index("    entry=core.read(");end=verification.index("    result=model.run()",start)
verification=verification[:start]+verification[end:]
start=verification.index("    names=(");end=verification.index("    new=",start)
verification=verification[:start]+"""    names=('unified_physics_condition_ledger_745.md','round745_drafts/research_note_745_draft.md',
           'round745_drafts/final_review.txt','round745_drafts/literature_scope_audit.json',
           'round745_drafts/scope_and_dedup_review.md','round746_drafts/STATUS.md',
           'round745_drafts/history_mechanism_probe.py','round745_drafts/history_mechanism_probe_results.json',
           'round745_drafts/exact_history_density.py','round745_drafts/exact_history_density_results.json',
           'round745_drafts/certify_history_integral.py','round745_drafts/certify_history_integral_results.json')
"""+verification[end:]
verification=verification.replace('original_full_H_selection_and_native_onsite_history_jet_checked=True',
    'exact_native_coefficients_and_strict_integral_sign_checked=True')
write('verify_round745.py',verification)
print('Prepared745 certification and746 connected-graph interface.')
