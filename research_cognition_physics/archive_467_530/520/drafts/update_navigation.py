"""Publish round 520 navigation only after checked science; preserve snapshots."""
from pathlib import Path
import hashlib
import json

HERE=Path(__file__).resolve().parents[1]
RESEARCH=HERE.parent
ROOT=RESEARCH.parent
checks=json.loads((HERE/'research_round_520_checks.json').read_text('utf8'))
assert checks['all_reported_checks_passed']
assert checks['fresh_tests']==dict(run=6,failures=0,errors=0)
assert checks['previous_protected_evidence_hashes_verified']==1004
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',
       HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
baseline={p:p.read_bytes() for p in paths}
snapshot=HERE/'round520_drafts/navigation_before_round520'
snapshot.mkdir(exist_ok=False)
planned={}
for p,raw in baseline.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    nl='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(enc).replace('\r\n','\n')
    prefix=('research_cognition_physics/archive_231_/' if p.parent==ROOT else
            'archive_231_/' if p.parent==RESEARCH else '')
    banner=(f'**第520轮完成：** [实际记录与当前关系态学习]({prefix}research_note_520.md)'
        '在同一原H与完整端口仪器下，对每个有限规模给共同步长并排除全部暗子空间，'
        '接成熟轨迹定理；从数学满秩猜测更新的估计器，对任意真实初态及旧参考有统一平均报告误差趋零。'
        '无需重置真实图或预给初态说明；读者、汇集、计算及资源仍输入，未得到固定经典图或空间坐标。'
        '6项、11式及独立终审通过；累计2552项、871份编号科学文件、1008份保护证据。三维与GR仍开放。')
    first,rest=body.split('\n',1)
    body=first+'\n\n'+banner+'\n'+rest
    if p==RESEARCH/'research_direction.md':
        old='最新科学轮次与检查数为519／2546'
        assert body.count(old)==1
        body=body.replace(old,'最新科学轮次与检查数为520／2552')
        body=body.replace('完成231—519轮。','完成231—520轮。',1)
        body+='''

**520后的当前接续：** 全部实际端口记录可提供当前活动态的数学估计来源，不能把它自动压缩为位置或三坐标。下一项只接受具有明确实际定位菜单的稳定摘要，或能减少其来源／尺度输入的证明；一般预测态闭合、继续计算纯度和放大遍历字均复用520，不新开轮。原精确操作语义、宏观有效空间方向与三维优先级不变。

**520冻结核验：** [代码](archive_231_/selective_record_state_learning.py)、[结果](archive_231_/selective_record_state_learning_results.json)、[科学核验](archive_231_/research_round_520_checks.json)、[只读复算](archive_231_/verify_selective_record_learning_round.py)。继承519后CRT审查的4份证据，原1004份全保留，再增3份科学文件和1份终审稿，共1008份。无规模一致等待率，也不把固定每步误差下的无限时滤波当作已证。
'''
    if p==RESEARCH/'RESEARCH_STATE.md':
        body=body.replace('已完成第231—519轮。','已完成第231—520轮。',1)
        body+='''

**520核验入口：** [笔记](archive_231_/research_note_520.md)、[保存结果](archive_231_/selective_record_state_learning_results.json)、[科学检查](archive_231_/research_round_520_checks.json)、[复算](archive_231_/verify_selective_record_learning_round.py)、[整合入口](archive_231_/verify_round520_integration.py)。编号检查2546＋6＝2552，科学文件868＋3＝871，保护证据1004＋4＝1008。得到条件状态学习能力，未生成稳定位置、三维或GR，阶段不结项。
'''
    if p==HERE/'README.md':
        body+='''

**520：** [完整端口记录的关系态学习](research_note_520.md)；[代码](selective_record_state_learning.py)、[结果](selective_record_state_learning_results.json)、[核验](research_round_520_checks.json)。共同采样步长与暗投影排除覆盖全部有限规模；经典记录用于更新当前量子态，不等于记录一个固定经典图。后继回到有任务来源的宏观位置摘要。
'''
    if p==HERE/'spatial_premise_closure_audit.md':
        body+='''

## 158. 第520轮：选择性历史与当前关系态的学习来源

本回合首先复核全部导航、519及保存结果。前一目标回合仅作解释和状态核对，按目标审计记为无新增进展；没有确认正在运行的旧计算，也未凭旧状态文件重启任务。本回合并行只读去重确认：472的替换式读取、492／517的非选择混合、491／493来源链及518局部包络均未证明完整选择性轨迹的当前态学习。

[520](research_note_520.md)以每棵树、每个根的Euler返回字构造完整旋转基效果。全部字使用同一个δ＝|J|／(16dLb²)，归一化效果误差小于1／d；任何秩r≥2暗投影在某射线的压缩离标量至少r／(2d)，故不存在暗子空间。接Maassen–Kümmerer Corollary 5的既有纯化结论，不另发明一般轨迹定理。

满秩数学猜测不是真实来源重置。对全部未知输入及旧参考，真实记录权重下的完整活动／参考报告误差≤min(1,2√(d m_n)＋d m_n)，m_n是可由同一模型有限字计算的估计纯度缺陷，并趋于零。保留真实条件参考与全部记录；不声称恢复单份初态、学会未知参考或控制新测量环境。固定有限目标精度可计完整仪器误差，但固定非零每步误差的无限时稳定性未证。

6项复算和11式独立终审通过，含全部有限字身份、全有理共同步长、实际Kraus字、1296份完整短历史、未知参考及范围反例。短历史总体误差界尚大于1，罕见字小误差诊断保留约10⁻¹¹⁸至10⁻¹¹⁶概率；不以这些数值宣称整体收敛速率。累计2552项，871份编号科学文件，1008份保护证据，正式稿与终审稿同字节。

这项进展减少了“观察者须预给未知当前图态说明”的来源要求，但仍输入精确模型、完整记录、汇集和滤波能力。宏观位置如何选择仍未解；不把普通状态估计、后验纯化或成熟谱工具自动称为三维，也不以获取完整状态作为所有定位方案的必要条件。后继只核任务相对的可持续定位摘要与来源，不追加更长遍历、更多纯度或一般预测态闭合轮次。目标继续活动。
'''
    if p==HERE/'three_dimensional_four_conditions_review.md':
        body+='''

## 63. 520：可学习当前关系态，仍需选择空间接口

[520](research_note_520.md)在保留完整实际记录的合同下，提供当前活动态的统一估计误差。条件纯化不等于某张经典树固定；满秩滤波初值也不是重新制造来源。局部紧的位置群、一致半幅、保组合重定向和完整方向合同均未因估计能力而自动成立。不得把状态矩阵参数个数当空间维数，或重新要求位置必须预测所有私人未来。编号520／2552，总保护证据1008；三维及GR仍开放。
'''
    planned[p]=body.replace('\n',nl).encode(enc)
    label=('project_README.md' if p==ROOT/'README.md' else
           'research_'+p.name if p.parent==RESEARCH else 'archive_'+p.name)
    with (snapshot/label).open('xb') as f:f.write(raw)
for p,raw in baseline.items():assert p.read_bytes()==raw,str(p)
for p,data in planned.items():
    assert p.read_bytes()==baseline[p],str(p)
    p.write_bytes(data)
with (snapshot/'manifest.json').open('x',encoding='utf8') as f:
    json.dump({str(p.relative_to(ROOT)):{'before':hashlib.sha256(baseline[p]).hexdigest(),
        'after':hashlib.sha256(data).hexdigest()} for p,data in planned.items()},f,indent=2,ensure_ascii=False)
print(json.dumps({'navigation_files_updated':len(planned),'snapshots_preserved':len(baseline)}))
