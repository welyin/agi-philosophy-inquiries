"""One-shot navigation update with byte preservation and compare-before-write."""
from pathlib import Path
import hashlib
import json

HERE=Path(__file__).resolve().parents[1]
RESEARCH=HERE.parent
ROOT=RESEARCH.parent
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',
       HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
baseline={p:p.read_bytes() for p in paths}
assert (HERE/'research_round_517_checks.json').exists()
assert (HERE/'qca_finite_region_scope_checks.json').exists()
snapshot=HERE/'round517_drafts/navigation_before_round517'
snapshot.mkdir(exist_ok=False)
planned={}
for p,raw in baseline.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    nl='\r\n' if b'\r\n' in raw else '\n'
    s=raw.decode(enc).replace('\r\n','\n')
    prefix=('research_cognition_physics/archive_231_/' if p.parent==ROOT else
            'archive_231_/' if p.parent==RESEARCH else '')
    banner=(f'**第517轮完成：** [任意有限规模的实际图来源]({prefix}research_note_517.md)'
            '将492六主体的端口监测混合推广到任意有限的不同连通图族；'
            '并证明全标号二叉树换边连通。任意初始活动态／被动参考可得到最大混合图来源，'
            '无需逐步重置数据或图；读取、计时及隔离仍输入，未给规模一致等待预算。'
            '5项检查、12式及独立终审通过；累计2536项、862份编号科学文件、987份保护证据。'
            '三维和GR尚未完成。\n\n'
            f'**既有文献接续校正：** [完整有限区域的近似操作]({prefix}qca_finite_region_scope_addendum.md)'
            '已有Schaeffer普通物理普适性构造；404限制固定程序强版本，不否定逐精度换程序。'
            '不再把修补418相干程序当作必经步骤；离散操作能力不签收连续Time或宏观维数。'
            '此来源审计不增加轮次，2份文稿单列继承。')
    first,rest=s.split('\n',1)
    s=first+'\n\n'+banner+'\n'+rest
    if p==RESEARCH/'research_direction.md':
        old='最新科学轮次与检查数为516／2531'
        assert s.count(old)==1
        s=s.replace(old,'最新科学轮次与检查数为517／2536')
        s=s.replace('完成231—516轮。','完成231—517轮。',1)
        s+='''

**517后的当前接续：** 各有限规模的均匀图来源已获条件性实现，不再外给每步独立新图。先核它与491／493的实际端口统计及504均值替换障碍的衔接；直接迁移、I＋1标签慢模或重复随机树矩不新开518。下一项必须减少实际定位菜单、尺度或端点选择的输入。测量自主来源、规模预算及三维选择仍开放。

**517冻结核验：** [科学检查](archive_231_/research_round_517_checks.json)、[可复算入口](archive_231_/verify_all_scale_graph_source_round.py)。完整继承516后5份HQCA证据及2份QCA来源审计；旧稿保留，正式稿与终审稿相同。索引没有把新增来源能力当作三维结项。
'''
    if p==RESEARCH/'RESEARCH_STATE.md':
        s=s.replace('已完成第231—516轮。','已完成第231—517轮。',1)
    if p==HERE/'spatial_premise_closure_audit.md':
        s+='''

## 154. 516后成熟来源接续：普通物理普适性已经覆盖完整胞元

回读404—405及Schaeffer原文定义6—7、定理5，确认普通版本允许程序随目标和精度改变，作用于完整有限原生区域。输入线性性统一环境末态，未知参考及半diamond误差可以直接接入。有限仪器和预定报告沿用230的酉扩张，不重新写门模拟；[范围补充](qca_finite_region_scope_addendum.md)及[核验](qca_finite_region_scope_checks.json)已独立终审。

这修正153末尾的后继选择：不必先修418相干程序才能取得完整有限区域的近似操作。418的字母守恒障碍与404正常固定程序障碍都保留。另一规则的逐任务程序、无限空白背景、离散更新、实际末程序及后续过程仍明示；不签收精确FUCP／连续Time，也不把原生一维充作宏观三维反证。此次只是成熟文献的重新接续，不登记科学轮次。两份原字节文稿继承旧980份证据，总计982。

## 155. 第517轮：把六主体的实际图来源推广到任意有限规模

[517](research_note_517.md)接正向宏观空间主线，回用440／492原H和已声明的端口Lüders仪器。没有把516持续记忆H的结果移植到另一H；没有重置D/G或输入逐步独立新图。

新增证明固定任意不同连通图对：外围算符的端口块系数在共有边相等、在对称差边两端均零；并图连通迫全部图间相干系数为零。各图内剩余标量再被连通图转换F统一。配合492的HS等号及无混叠工具，得到每个有限图族的完整primitive分类和未知活动态／参考的diamond收敛。全标号二叉树用毛虫化、内部标签两次分支互换及侧叶排序证明F连通；不是从有限枚举外推。

有限验证覆盖1、6、90、2520张图，90图全部4005对整数约束、原720维单包H及未知参考；精确边界例保留图标签可区分、物理连通和采样前提。等待预算只逐有限规模存在，未给规模一致速度。式(11)—(12)是493来源选择工具的直接推广，不重复包装成新重试定理。

独立终审指出能量常数说明和可计算参数范围须明确：现保留ΣSWAP内的J|E|I，只省另一μ|E|I；有理搜索限定可计算系数。旧稿保留，科学代码／结果无需修改；最终稿与终审稿原字节一致。5项、12式通过，累计2536项、862份编号科学文件。完整继承旧982份，新增代码／结果／笔记与2份草稿，总计987份保护证据。没有图像检查或新应用任务。

这是可重复统计来源的规模缺口得到条件性解决，不是自主读取装置、位置、维数或引力的完成。立刻接续491／493与504做去重：若来源桥只让已知I＋1慢模和固定短时误差地板具备另一实施途径，则记推论，不新增518；真正后继须对同一实际定位菜单产生跨尺度新预测，或减少端点／尺度选择输入。
'''
    if p==HERE/'three_dimensional_four_conditions_review.md':
        s+='''

### 59.3 完整胞元近似能力已有另一文献来源

[QCA补充](qca_finite_region_scope_addendum.md)核清普通物理普适性覆盖完整有限区域及未知参考；404不排除这个版本。此能力不能代替实际定位、连续Time或宏观空间来源，原四项条件没有因此被证明；不增加编号检查。

## 60. 517：全有限规模来源仍不等于维数选择

[517](research_note_517.md)将实际监测产生的最大混合图来源推广到每个有限连通图族，移除六主体限制及预给均匀图边缘的一项输入。完整图族、读取、计时、隔离仍为合同；没有统一规模等待界，也没有让全部新读者与对象共同去关联。

旧491一般I族的I＋1慢方向未变，来源能力不会把它改成3。这里既不宣布最大混合意味着空间消失，也不把图族大小当维数。实际位置、一致半幅、保组合重定向和完整空间方向合同继续等待来源。编号517／2536，三维阶段和GR仍开放。
'''
    planned[p]=s.replace('\n',nl).encode(enc)
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
