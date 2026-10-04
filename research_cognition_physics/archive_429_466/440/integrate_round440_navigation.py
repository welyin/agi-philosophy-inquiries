"""Guarded navigation update for simultaneous edge exchange and relational frames."""
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASE, ROOT = HERE.parent, HERE.parent.parent
PATHS = [ROOT/'README.md', BASE/'README.md', BASE/'research_direction.md',
         BASE/'RESEARCH_STATE.md', HERE/'README.md', HERE/'spatial_premise_closure_audit.md']
original = {p: p.read_bytes() for p in PATHS}
texts = {p: raw.decode('utf-8').replace('\r\n', '\n') for p, raw in original.items()}
assert '**第440轮完成：**' not in texts[BASE/'RESEARCH_STATE.md']
summary = ('**第440轮完成：** [同时换边与动态关系表示](%sresearch_note_440.md)'
           '复用Quantum Graphity保度数换边项，无需临时空槽即恢复长链和长环重组；'
           '给与439合并后的C＝2完整配置分类、规模一致端点范数及未知输入概率证书。'
           '沿链表示保留标签记录，数据作用和主体读数同时变换，192维整数核验通过。'
           '7项检查，累计2071项；631份编号科学文件、666份保护证据。'
           '初始图与换边项仍为输入；后继核真实传播接口，三维与GR未完成。')
for path in PATHS[:5]:
    prefix = ('research_cognition_physics/archive_231_/' if path == ROOT/'README.md'
              else ('' if path == HERE/'README.md' else 'archive_231_/'))
    old = next(line for line in texts[path].splitlines() if '**第439轮完成：**' in line)
    texts[path] = texts[path].replace(old, old+'\n\n'+('> ' if path == ROOT/'README.md' else '')+summary % prefix, 1)
    texts[path] = texts[path].replace('231—439轮', '231—440轮')
for path in (BASE/'research_direction.md', HERE/'README.md'):
    old = next(line for line in texts[path].splitlines() if line.startswith('**下一项：没有临时空槽时能否自主交换关系。**'))
    texts[path] = texts[path].replace(old,
        '**下一项：动态标签下的真实传播距离。** 440已给无需临时空槽的完整重组接口、未知数据证书和精确沿链表示，'
        '主体标签不能删去，真实读数须同步变换。先复用成熟局域传播定理及已有量子传播轮次，'
        '核动态标签、任意未知数据及实际主体读数下，哪些信号界沿当前关系成立、哪些只依赖给定路径部门。'
        '不重做单粒子链扩散、配置枚举或经典混合。如果只是旧模型重述，转查式(14)允许的虚过程／内部载体以削减换边项输入，'
        '不继续无来源地添加图规则。初始分区、维数选择和新主体接入仍开放。'
        '\n\n**439后同时换边任务（已完成模型与完整数据接口）：** '+old.split('** ', 1)[1], 1)
texts[BASE/'research_direction.md'] = texts[BASE/'research_direction.md'].replace('最新科学轮次与检查数为439／2064', '最新科学轮次与检查数为440／2071')
texts[BASE/'RESEARCH_STATE.md'] = texts[BASE/'RESEARCH_STATE.md'].replace('最新439轮及累计2064项见本文开头', '最新440轮及累计2071项见本文开头')
readme = BASE/'README.md'
texts[readme] = texts[readme].replace('当前复算与冻结入口：', '当前复算与冻结入口：[440科学核验](archive_231_/verify_simultaneous_edge_exchange_round.py)、[440整合核验](archive_231_/verify_round440_integration.py)；历史入口：', 1)
texts[readme] = texts[readme].replace('第三阶段累计2064项；旧科学证据保持原字节', '第三阶段累计2071项；旧科学证据保持原字节')
archive = HERE/'README.md'
texts[archive] = texts[archive].replace('当前完成439轮', '当前完成440轮').replace('当前439不作为预定终点', '当前440不作为预定终点')
row = next(line for line in texts[archive].splitlines() if line.startswith('| [439：'))
texts[archive] = texts[archive].replace(row, row+'\n| [440：同时换边与动态关系表示](research_note_440.md) | 无临时空槽的重组；C＝2组合分类；完整数据与标签读数的沿链变换 | [代码](simultaneous_edge_exchange_audit.py)、[结果](simultaneous_edge_exchange_audit_results.json)、[核验](research_round_440_checks.json)；7项检查；初图与换边生成元仍为输入 |', 1)
texts[archive] = texts[archive].replace('当前整合入口：', '当前整合入口：[440整合核验](verify_round440_integration.py)、[整合记录](round440_integration_checks.json)；历史入口：', 1)
texts[archive] = texts[archive].replace('当前2064项科学检查；663份保护证据，其中本阶段编号科学文件628份', '当前2071项科学检查；666份保护证据，其中本阶段编号科学文件631份')
audit = HERE/'spatial_premise_closure_audit.md'
assert '## 68.' not in texts[audit]
texts[audit] += '''
## 68. 第440轮：同时换边恢复长链重组，沿链读数必须保留主体标签

2026-09-24。接439，直接复用Quantum Graphity 0801.0861式(27)的保度数换边模式，以边qubit及1/4有序和明确归一化；每条无向三步路径给单位振幅，不同中间边可叠加。无临时资源槽，纯换边分支逐D_i、f_i、Q_i和分区守恒；与439合并只保Q与分区，不误称全部逐D仍守恒。端点翻转范数≤2|kappa|C(C−1)^2，固定C时规模一致。

路径换边生成全部内部相邻排列，s点固定端点路径可达(s−2)!图；s点环可达(s−1)!/2图。一般连通正则图可达性明确归于flip文献，经典混合未移植。C＝2与439合并后的全部配置数为4^r乘各长路径与长环因子，1858个六点合法图完整核验；长链和长环的原冻结障碍已解除，但主体所属分量、长路径端点标签及N仍不改变。

四点路径部门给所有未知数据及参考的严格短时转移下界18769/345960000；对称数据时p=sin²t，明确与随机链趋于1/2不同。未排除大系统局部弛豫。由完整联合酉性保留信息，不声称所有单主体边缘不变。

新增精确沿链框架：π将位置映到主体标签，受控重排V把数据边作用变成固定路径交换，同时把图换边变成数据SWAP乘标签邻位换序。标签i的读数则成为随π选择位置π^−1(i)的算符。192维全矩阵及各主体Z读数有整数证书；覆盖未知数据、关系叠加和参考。V是数学基底变换，不声称自主物理执行，不能删掉图标签而宣称已有裸固定空间。

两个逐顶点同度的不同简单图最少改变四条边，因此原边载体上≤3因子局域且精确逐度守恒的H不能改图。辅助编码、暂时违度的虚过程和436模拟不被排除；新增F并非429已自然选定。下一项核动态标签下真正的传播／测距接口，复用已有传播理论，若重复则转查关系项的内部来源，不继续任意叠加规则。

7项检查、14个公式；累计2071项、631份编号科学文件、666份保护证据。独立只读终审与本地证书完成，旧科学字节不变，实际三维与GR目标继续。
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
