"""One-time guarded navigation update for the local-halving bridge."""
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
ROOT = BASE.parent
PATHS = [ROOT/'README.md', BASE/'README.md', BASE/'research_direction.md',
         BASE/'RESEARCH_STATE.md', HERE/'README.md', HERE/'spatial_premise_closure_audit.md']
original = {p: p.read_bytes() for p in PATHS}
texts = {p: raw.decode('utf-8').replace('\r\n', '\n') for p, raw in original.items()}
assert '**第425轮完成' not in texts[BASE/'RESEARCH_STATE.md'], 'Already integrated'

summary = ('**第425轮完成：** [一致半幅与光滑局部坐标]({prefix}research_note_425.md)证明：'
           '局部有界成本与双侧重放一致性给无小子群，局部紧性再接成熟Lie定理；无需半幅保乘法或全局收缩。'
           'SU(2)正例保留非交换图册，三维须另给同一小壳的完整qubit方向合同。'
           '8项检查，累计1977项；586份编号科学文件、621份保护证据。'
           '实际端点、半幅及资源账的来源仍开放，未无条件生成三维。')
for path in PATHS[:5]:
    prefix = ('research_cognition_physics/archive_231_/' if path == ROOT/'README.md'
              else ('' if path == HERE/'README.md' else 'archive_231_/'))
    marker = next(line for line in texts[path].splitlines() if '**第424轮完成：**' in line)
    insertion = ('> ' if path == ROOT/'README.md' else '') + summary.format(prefix=prefix)
    texts[path] = texts[path].replace(marker, marker+'\n\n'+insertion, 1)
    texts[path] = texts[path].replace('231—424轮', '231—425轮')
    texts[path] = texts[path].replace('不能自动替代端点重放、相容收缩和完整方向合同；424已削去全比例缩放的要求。',
        '不能自动替代实际端点重放和完整方向合同；424削去全比例缩放，425另给一致局部半幅通向Lie图册的桥梁，其实际来源仍须核查。')

next_step = ('**下一项：实际关系上的重放与半幅一致性。** 425取得局部半幅→NSS→Lie图册的替代桥，'
             '不要求全局收缩或交换性。继续检验哪些认知关系改变具有同一端点身份、可跨基点重放，'
             '以及同一内部资源账能否支持一致半幅。局部紧性和该实际小壳上的完整qubit读取仍是分列输入；'
             '不把任选数学范数当成已给物理成本，不重复状态压缩、控制幅度调参及旧端点商反例。'
             '像与透镜暂缓；三维、阶段结项及GR仍未完成。')
direction = BASE/'research_direction.md'
texts[direction] = texts[direction].replace('最新科学轮次与检查数为424／1969', '最新科学轮次与检查数为425／1977')
old_heading = '**下一项：一个端点收缩的实际来源。**'
assert old_heading in texts[direction]
texts[direction] = texts[direction].replace(old_heading,
    next_step+'\n\n**424后的历史问题（原收缩分支保留）：**', 1)
state = BASE/'RESEARCH_STATE.md'
texts[state] = texts[state].replace('最新424轮及累计1969项见本文开头', '最新425轮及累计1977项见本文开头')

readme = BASE/'README.md'
assert '当前复算与冻结入口：' in texts[readme]
texts[readme] = texts[readme].replace('当前复算与冻结入口：',
    '当前复算与冻结入口：[425科学核验](archive_231_/verify_local_halving_round.py)、[425整合核验](archive_231_/verify_round425_integration.py)；历史入口：', 1)
texts[readme] = texts[readme].replace('第三阶段累计1969项；旧科学证据保持原字节', '第三阶段累计1977项；旧科学证据保持原字节')

archive = HERE/'README.md'
texts[archive] = texts[archive].replace('当前完成424轮', '当前完成425轮').replace('当前424不作为预定终点', '当前425不作为预定终点')
row = next(line for line in texts[archive].splitlines() if line.startswith('| [424：'))
texts[archive] = texts[archive].replace(row, row+'\n| [425：一致半幅与光滑局部坐标](research_note_425.md) | 双侧重放一致性给NSS；局部紧接Lie图册；SU(2)非交换例及局部紧／误差边界 | [代码](local_halving_coordinate_audit.py)、[结果](local_halving_coordinate_audit_results.json)、[核验](research_round_425_checks.json)；8项检查；三维方向合同仍独立 |', 1)
texts[archive] = texts[archive].replace('当前整合入口：',
    '当前整合入口：[425整合核验](verify_round425_integration.py)、[整合记录](round425_integration_checks.json)；历史入口：', 1)
texts[archive] = texts[archive].replace('当前1969项科学检查；618份保护证据，其中本阶段编号科学文件583份',
    '当前1977项科学检查；621份保护证据，其中本阶段编号科学文件586份')
assert old_heading in texts[archive]
texts[archive] = texts[archive].replace(old_heading,
    next_step+'\n\n**424后的历史问题（原收缩分支保留）：**', 1)
texts[archive] += '\n425以424为科学基线，只读复算：`python -B -X utf8 research_cognition_physics/archive_231_/verify_local_halving_round.py`。8项检查、16个编号公式；一般Lie图册不能沿用413的加法交会。\n'

audit = HERE/'spatial_premise_closure_audit.md'
assert '## 52.' not in texts[audit]
texts[audit] += '''
## 52. 第425轮：把幅度调节落实为可检验的一致半幅

2026-09-24。用户再提出连续变化、资源、大小幅度及方向可选；接续本审计第41、44、50节，不重做原猜想、385周期例或415控制拓扑。先读导航、424笔记与结果；两个独立只读审计分别检查收缩来源和新的局部半幅桥。前者主要落入旧范围，后者有实质新连接。

### 52.1 新桥及不能混淆的比较范围

[425](research_note_425.md)保留实际端点群与局部紧性，改查一个小邻域上的双侧半幅一致性和有界成本比例降低。若一个子群完全留在邻域，连续倍增结合半幅一致性迫使每个元素成本为零，所以该子群平凡。这直接给NSS；已核成熟Hilbert第五问题定理给有限维Lie结构及真实局部坐标图。

不要求半幅保持不同位移的组合，也不要求全局收缩自同构。SU(2)小主对数球满足新合同，根映射却非同态，唯一二阶元素阻止全局收缩；它具体表明局部坐标不必是全局加法坐标。这里是替代桥，不能宣称整个新合同在同一给定预算下严格弱于424；从Lie结构任选范数也不等于真实资源账已经有该成本。

三维仍须在同一完整指数小壳另给383—384的qubit反向分离、满对比、协变和极小重定向。一般预算水平集未被证明为球面。原413加法交会不沿用；SU(2)示例的数值坐标读取另计成本和概率访问。没有将内部姿态群自动认作现实空间。

### 52.2 被取消的噪声来源候选

普通去极化、扩大环境后恢复及摘要商下降主要由223、235、390、221、415、417覆盖，不开新轮。必须避免两项过强推断：通道本身没有CP逆，不足以否定它在某个操作群上诱导自同构；整体可逆也不禁止某个位移部门收缩。真正需要核对的是同一端点关系和实际权限。相关无限维平移正例预先输入了平移部门，并未减少来源假设，因此也不另计成果。

### 52.3 边界、核验与后继

无限圆环积有连续正定预算及单射精确收缩，却有任意小尾部子群，不能把它当一致半幅。无限Hilbert加法群在单位球上有一致半幅及有界成本，却无局部紧性，故有限维不是由幅度调节单独给出。正误差的递推只给分辨率下限，不能据有限容差认证精确NSS。

8项检查、16个编号公式通过，累计1977项；586份编号科学文件、621份保护证据。独立代理核对成熟定理与383—386接口，并只读复算代码及正式结果；终审补明指数图域、壳中心及Hilbert例的局部预算范围。历史科学文件和旧核验保持原字节，无图像检查。

后继集中于实际关系的重放和半幅来源，不延长收缩率、误差或预算扫描。一般认知状态的连续变化与实际定位关系仍须区分，但不引入认知之外的空间，也不要求宇宙有全局加法坐标。当前没有充分依据触发空间阶段收尾或完成GR目标；原342—343形式合同不足性保留，强化内部实施范围仍独立开放。
'''

prepared = {}
for path, content in texts.items():
    newline = '\r\n' if b'\r\n' in original[path] else '\n'
    prepared[path] = (content.rstrip()+'\n').replace('\n', newline).encode('utf-8')
for path in PATHS:
    assert path.read_bytes() == original[path], f'Concurrent edit: {path}'
for path in PATHS:
    assert path.read_bytes() == original[path], f'Concurrent edit: {path}'
    path.write_bytes(prepared[path])
    print(path.relative_to(ROOT))
