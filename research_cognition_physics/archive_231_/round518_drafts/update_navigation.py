"""Preserve seven navigation files and atomically check their expected versions."""
from pathlib import Path
import hashlib
import json

HERE=Path(__file__).resolve().parents[1]
RESEARCH=HERE.parent
ROOT=RESEARCH.parent
checks=json.loads((HERE/'research_round_518_checks.json').read_text('utf8'))
assert checks['all_reported_checks_passed']
assert checks['fresh_tests']==dict(run=6,failures=0,errors=0)
assert len(checks['new_file_hashes'])==3 and len(checks['preserved_draft_hashes'])==2
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',
       HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
baseline={p:p.read_bytes() for p in paths}
snapshot=HERE/'round518_drafts/navigation_before_round518'
snapshot.mkdir(exist_ok=False)
planned={}
for p,raw in baseline.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    nl='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(enc).replace('\r\n','\n')
    prefix=('research_cognition_physics/archive_231_/' if p.parent==ROOT else
            'archive_231_/' if p.parent==RESEARCH else '')
    banner=(f'**第518轮完成：** [同一结束时刻的局部邻域来源]({prefix}research_note_518.md)'
            '在原持续换边H下，以四次实际端口记录联合认证叶根的当前半径2包络；'
            '全规模／未知图参考的成功与坏率界不含N倍因子。'
            '单位示例成功率至少2⁻¹⁶⁰、条件坏率小于1／5000，低来源率完整保留。'
            '读取、计时、交付及局部量子摘要仍需输入，未生成宏观空间。'
            '6项、16式及独立终审通过；累计2542项、865份编号科学文件、992份保护证据。'
            '三维和GR尚未完成。')
    first,rest=body.split('\n',1)
    body=first+'\n\n'+banner+'\n'+rest
    if p==RESEARCH/'research_direction.md':
        old='最新科学轮次与检查数为517／2536'
        assert body.count(old)==1
        body=body.replace(old,'最新科学轮次与检查数为518／2542')
        body=body.replace('完成231—517轮。','完成231—518轮。',1)
        body+='''

**518后的当前接续：** 实际四跳提供同刻随机P₂包络，接498／499的按记录预测与温和投影是直接推论，不另开519。下一项优先核定位任务真正需要哪些局部量子数据、边相干及读口，及这些如何在同一实际过程被使用；经典成员表不能替代量子摘要。更大包络、有效时间窗和跨尺度定位仍缺来源，不为半径或常数扩展自动增加轮次。

**518冻结核验：** [科学检查](archive_231_/research_round_518_checks.json)、[复算入口](archive_231_/verify_same_time_neighborhood_round.py)。旧987份证据保持，新增3份科学文件与2份稿件；延迟事件已明确改查图与记录，不要求包停在原末端口。
'''
    if p==RESEARCH/'RESEARCH_STATE.md':
        body=body.replace('已完成第231—517轮。','已完成第231—518轮。',1)
    if p==HERE/'spatial_premise_closure_audit.md':
        body+='''

## 156. 第518轮：真实记录生成同刻半径2包络

上一目标回合完成517的全规模来源及两项成熟接口去重，留下987份已核保护证据，判为有效进展。本回合先读导航、517及结果，检查518／519占用；没有凭进程查询的权限错误判定旧研究终止，也没有重启或中止不明进程。

接155.2候选，独立回读496、498—500、502。四条异时边的证书不能直接相乘；本轮将四段都比较到F相互作用绘景的同一最终时刻4t。完整端口记录的非停留映射范数≤3|J|t，四段替换误差≤|J|⁴(10368|κ|＋1944|J|)t⁵，不按N⁴条记录累加。变回实际绘景后领先项为B exp(−4itκF)，输出满足当前三边，且B†B＝2I。

[518](research_note_518.md)由此给覆盖任意未知G/R的成功下界及联合当前坏率。单位t＝2⁻²⁰时，成功至少2⁻¹⁶⁰、条件坏率<1／5000。该罕见来源不能免费条件化成高成功率；每条指定记录也不保证相同置信度。六图中一条不可能完整星形的记录严格发生，指定记录坏率为1；完整随机输出仍满足保证。

图数1、6、90、2520的领先效果与当前B₂逐整数核对；旧500原H同构、四段复整数Kraus及全部未知图Gershgorin证书通过；1296条完整历史保存失败及参考，较大t＝1／8的诊断坏率约0.4926，表明时间范围不能省略。解析证明承担任意规模，不由这些数值外推。

独立终审发现延迟时必须移去源末端口D=w的投影。正式稿明确使用Q_graph＝Σ_r I_D⊗π_r⊗|r><r|，包继续运动本身不算图失配；旧稿保留，代码／结果无需更改。最终稿与终审稿逐字节相同。6项、16式核验通过；累计2542项、865份编号科学文件、992份保护证据，旧987份不改。没有图像检查或新应用任务。

此轮只减少499的一份随机近似P₂来源输入，未给任意大R、量子摘要获取设施、预测器的自主实现、位置或维数。叶源、三度角色、单包、读口与时序仍是明确模型条件，不升级为新的认知原则。下一项已立即去重如下，不以局部技术补充代替三维目标。

### 156.1 与已有预测器的直接连接，不另开519

回读498原完整原子截断和499联合温和投影：在成功的经典记录报告上，按r控制同一个规则定义的Φ_{A_r}与局部预测器即可。联合坏质量q给完整根／参考报告误差≤E₂(s)＋2√q；有延迟则先使用518式(14)，有实施误差再逐项记账。这里只比较声明的报告，不要求位置预测全部私人未来。

Φ_{A_r}保留的是A内实际数据与全部潜在边的量子态，包括可读相干；四个经典标签不提供这些态。预测器仍必须用498完整原子截断，不能直接截取单包Laplacian后冒充同一局部通道。R=2时旧截断界只有floor(2/3)＋1＝1阶，不因得到当前列表而自动获得大的有效窗口。

上述是旧定理的直接合成，没有新代码或科学编号。真正后继须说明：所选实际定位菜单要使用哪些量子变量、怎样在原过程内访问／比较它们，或怎样由同一机制取得跨尺度有效窗口。若只扩大四跳的半径、复制列表或优化罕见成功率，不作为新的三维进展。当前已经完成一项明确来源增量，但三维与GR目标继续开放，不触发结项或阻塞状态。
'''
    if p==HERE/'three_dimensional_four_conditions_review.md':
        body+='''

## 61. 518：同刻邻域来源仍不等于空间坐标

[518](research_note_518.md)首次在本具体候选中用实际路径记录提供叶根的随机近似半径2包络，并把三边放在同一最终时刻检验；没有把不同停止时刻的记录拼图。来源稀有，当前列表有有限有效期，完整读口和交付仍输入。

这一来源可接498／499的旧局部预测工具，但列表不能替代A内量子态及边相干，也不生成位置群、一致半幅、保组合重定向或完整方向接口。扩大局部来源能力没有选择维数3。编号518／2542，阶段与GR继续开放；只记直接合成，不重复开519。
'''
    planned[p]=body.replace('\n',nl).encode(enc)
    label=('project_README.md' if p==ROOT/'README.md' else
           'research_'+p.name if p.parent==RESEARCH else 'archive_'+p.name)
    with (snapshot/label).open('xb') as f:f.write(raw)
for p,raw in baseline.items():
    assert p.read_bytes()==raw,str(p)
for p,data in planned.items():
    assert p.read_bytes()==baseline[p],str(p)
    p.write_bytes(data)
with (snapshot/'manifest.json').open('x',encoding='utf8') as f:
    json.dump({str(p.relative_to(ROOT)):{'before':hashlib.sha256(baseline[p]).hexdigest(),
        'after':hashlib.sha256(data).hexdigest()} for p,data in planned.items()},f,indent=2,ensure_ascii=False)
print(json.dumps({'navigation_files_updated':len(planned),'snapshots_preserved':len(baseline)}))
