"""Guarded navigation update for a moving receiver's propagation bound."""
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASE, ROOT = HERE.parent, HERE.parent.parent
PATHS = [ROOT/'README.md', BASE/'README.md', BASE/'research_direction.md',
         BASE/'RESEARCH_STATE.md', HERE/'README.md', HERE/'spatial_premise_closure_audit.md']
original = {p: p.read_bytes() for p in PATHS}
texts = {p: raw.decode('utf-8').replace('\r\n', '\n') for p, raw in original.items()}
assert '**第441轮完成：**' not in texts[BASE/'RESEARCH_STATE.md']
summary = ('**第441轮完成：** [动态标签下的主体接收与完整传播界](%sresearch_note_441.md)'
           '在440同一演化中，将跟随主体的真实读数接到成熟局域传播界；唯一标签位移尾界与有界读数延拓，'
           '给两个任意局部CPTP编码的完整接收者／参考迹距离上界。J＝0时主体间无信号，旧槽位却可有最大读数差。'
           '6项检查，累计2077项；634份编号科学文件、669份保护证据。'
           '初始路径和发送／接收标签位置仍为输入，后继核一般关系图；三维与GR未完成。')
for path in PATHS[:5]:
    prefix = ('research_cognition_physics/archive_231_/' if path == ROOT/'README.md'
              else ('' if path == HERE/'README.md' else 'archive_231_/'))
    old = next(line for line in texts[path].splitlines() if '**第440轮完成：**' in line)
    texts[path] = texts[path].replace(old, old+'\n\n'+('> ' if path == ROOT/'README.md' else '')+summary % prefix, 1)
    texts[path] = texts[path].replace('231—440轮', '231—441轮')
for path in (BASE/'research_direction.md', HERE/'README.md'):
    old = next(line for line in texts[path].splitlines() if line.startswith('**下一项：动态标签下的真实传播距离。**'))
    texts[path] = texts[path].replace(old,
        '**下一项：传播桥梁能否离开预置路径次序。** 441补齐初始发送与接收位置明确时的完整主体信号界，'
        '保留其他关系相干与未知参考，但没有生成初始路径。先审计一般分叉占用图及初始主体离域：'
        '已有动态交互图传播定理能否适用同一关系律，是否仍需隐藏固定槽位。'
        '不重做固定链扩散或已知LR定理；若仍必须输入固定几何，就转查440第6节允许的内部载体／虚过程，'
        '明确能削减哪项关系更新输入，不继续任意添图规则。初始分区、三维选择和新主体接入仍开放。'
        '\n\n**440后动态标签任务（已完成指定路径部门的完整接收接口）：** '+old.split('** ', 1)[1], 1)
texts[BASE/'research_direction.md'] = texts[BASE/'research_direction.md'].replace('最新科学轮次与检查数为440／2071', '最新科学轮次与检查数为441／2077')
texts[BASE/'RESEARCH_STATE.md'] = texts[BASE/'RESEARCH_STATE.md'].replace('最新440轮及累计2071项见本文开头', '最新441轮及累计2077项见本文开头')
readme = BASE/'README.md'
texts[readme] = texts[readme].replace('当前复算与冻结入口：', '当前复算与冻结入口：[441科学核验](archive_231_/verify_moving_subject_propagation_round.py)、[441整合核验](archive_231_/verify_round441_integration.py)；历史入口：', 1)
texts[readme] = texts[readme].replace('第三阶段累计2071项；旧科学证据保持原字节', '第三阶段累计2077项；旧科学证据保持原字节')
archive = HERE/'README.md'
texts[archive] = texts[archive].replace('当前完成440轮', '当前完成441轮').replace('当前440不作为预定终点', '当前441不作为预定终点')
row = next(line for line in texts[archive].splitlines() if line.startswith('| [440：'))
texts[archive] = texts[archive].replace(row, row+'\n| [441：动态标签下的主体接收与完整传播界](research_note_441.md) | 唯一标签尾界与真实主体读数；两任意CPTP编码、未知参考下的传播上界 | [代码](moving_subject_propagation_audit.py)、[结果](moving_subject_propagation_audit_results.json)、[核验](research_round_441_checks.json)；6项检查；初始路径及双方位置仍为输入 |', 1)
texts[archive] = texts[archive].replace('当前整合入口：', '当前整合入口：[441整合核验](verify_round441_integration.py)、[整合记录](round441_integration_checks.json)；历史入口：', 1)
texts[archive] = texts[archive].replace('当前2071项科学检查；666份保护证据，其中本阶段编号科学文件631份', '当前2077项科学检查；669份保护证据，其中本阶段编号科学文件634份')
audit = HERE/'spatial_premise_closure_audit.md'
assert '## 69.' not in texts[audit]
texts[audit] += '''
## 69. 第441轮：动态关系中的真实主体传播，路径部门仍是输入

2026-09-24。接440，不新增动力学项。在每槽数据C²、标签Cˢ的完整张量空间，用数据SWAP和数据＋标签复合SWAP精确补全路径部门；唯一标签、端点准备和编码有资源代价，每槽维数2s随规模增长。四、五主体全部32、192列稀疏整数交织核验，无低能或有限窗编码误差；并非免费生成固定链。

唯一码中单个标签的右移T为范数≤1的部分等距；加权生成元给与未知数据、其他标签相干和参考一致的位移概率尾界。复用Nachtergaele–Sims的局域对易子迭代，常数6(|J|＋|kappa|)不依赖网络规模和参考维数。主体读数在完整张量空间须用恰一个接收标签的有界延拓；直接求和在重复标签部门可能范数膨胀，此处已明确排除。

两个主体初始分别位于确定槽α、β，局部编码仅作用发送数据。两个任意CPTP编码后，接收者与任意保留参考的迹距离≤min(1,2S_{d−R}(6g|t|)＋q_R)，涵盖其他标签叠加与全数据纠缠。条件期望／Haar平均仅作证明，不是新增物理操作。d＝12、t＝1/100等参数给严格小于10⁻⁸的全规模有理证书；五主体64维输入及参考给非零信号见证，不误称数值计算diamond范数。

J＝0时图仍换位，所有主体数据严格不受作用，因此无主体间数据信号；四点例在t＝π/2的接收主体Z差为0、接收者原槽位Z差却为2。实际测量必须跟踪主体身份，不能用固定位置的读数变化冒充主体间通信。

下一项核一般分叉图、初始离域主体及已有动态交互图传播定理，判断是否能摆脱预置路径次序；若需重新给固定几何，则回到440的内部载体／虚过程来源。未由此选择三维、生成钟尺或证明严格相对论光锥。6项检查、13个公式；累计2077项、634份编号科学文件、669份保护证据，旧科学字节不变。
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
