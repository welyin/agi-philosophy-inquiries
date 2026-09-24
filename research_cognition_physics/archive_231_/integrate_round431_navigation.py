"""Guarded navigation update for the internal exchange reader."""
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASE, ROOT = HERE.parent, HERE.parent.parent
PATHS = [ROOT/'README.md', BASE/'README.md', BASE/'research_direction.md',
         BASE/'RESEARCH_STATE.md', HERE/'README.md', HERE/'spatial_premise_closure_audit.md']
original = {p: p.read_bytes() for p in PATHS}
texts = {p: raw.decode('utf-8').replace('\r\n', '\n') for p, raw in original.items()}
assert '**第431轮完成：**' not in texts[BASE/'RESEARCH_STATE.md']
summary = ('**第431轮完成：** [持续交换与内部读者]({prefix}research_note_431.md)'
           '在三体对象上加入两qubit读者，四条交换一直开启，关系信息进入读者联合态。'
           '整数／有理证书给时刻3／2的迹距离大于2／3，整个[17／12,19／12]窗口的等先验理论辨识成功率大于3／4。'
           '无需精确停机；初态及参与关系仍为输入，不声称永久记录或完整读出器。'
           '6项检查，累计2014项；604份编号科学文件、639份保护证据。继续核交换参与结构的内部来源，三维与GR未完成。')
for path in PATHS[:5]:
    prefix = ('research_cognition_physics/archive_231_/' if path == ROOT/'README.md'
              else ('' if path == HERE/'README.md' else 'archive_231_/'))
    old = next(line for line in texts[path].splitlines() if '**第430轮完成：**' in line)
    texts[path] = texts[path].replace(old, old+'\n\n'+('> ' if path == ROOT/'README.md' else '')+summary.format(prefix=prefix), 1)
    texts[path] = texts[path].replace('231—430轮', '231—431轮')
for path in (BASE/'research_direction.md', HERE/'README.md'):
    old = next(line for line in texts[path].splitlines() if line.startswith('**下一项：同一交换机制怎样产生内部记录。**'))
    texts[path] = texts[path].replace(old,
        '**下一项：交换参与关系和初态资源的内部来源。** 429—431已依次得到条件性交互形式、三方关系量子块、'
        '同一持续交换规则下的内部可读转移。继续此方向，先核完全均匀且无额外参与结构的多体交换所允许的关系动力学，'
        '再定位具有认知动机的最小补充；也须核singlet准备与读者复用。'
        '不继续扫描五体链长度、读取峰值或精确门编译，不把已给参与图认作已生成空间。'
        '有限可读窗口不要求永久存储；实际坐标来源仍开放。\n\n**430后内部读者任务（已完成有限窗口部分）：** '+old.split('** ', 1)[1], 1)
texts[BASE/'research_direction.md'] = texts[BASE/'research_direction.md'].replace('最新科学轮次与检查数为430／2008', '最新科学轮次与检查数为431／2014')
texts[BASE/'RESEARCH_STATE.md'] = texts[BASE/'RESEARCH_STATE.md'].replace('最新430轮及累计2008项见本文开头', '最新431轮及累计2014项见本文开头')
readme = BASE/'README.md'
texts[readme] = texts[readme].replace('当前复算与冻结入口：', '当前复算与冻结入口：[431科学核验](archive_231_/verify_exchange_reader_round.py)、[431整合核验](archive_231_/verify_round431_integration.py)；历史入口：', 1)
texts[readme] = texts[readme].replace('第三阶段累计2008项；旧科学证据保持原字节', '第三阶段累计2014项；旧科学证据保持原字节')
archive = HERE/'README.md'
texts[archive] = texts[archive].replace('当前完成430轮', '当前完成431轮').replace('当前430不作为预定终点', '当前431不作为预定终点')
row = next(line for line in texts[archive].splitlines() if line.startswith('| [430：'))
texts[archive] = texts[archive].replace(row, row+'\n| [431：持续交换中的内部读者](research_note_431.md) | 两qubit载体接收三方关系；精确有理证书与非零宽度读取窗口 | [代码](exchange_reader_audit.py)、[结果](exchange_reader_audit_results.json)、[核验](research_round_431_checks.json)；6项检查；参与关系及初态仍为输入 |', 1)
texts[archive] = texts[archive].replace('当前整合入口：', '当前整合入口：[431整合核验](verify_round431_integration.py)、[整合记录](round431_integration_checks.json)；历史入口：', 1)
texts[archive] = texts[archive].replace('当前2008项科学检查；636份保护证据，其中本阶段编号科学文件601份', '当前2014项科学检查；639份保护证据，其中本阶段编号科学文件604份')
audit = HERE/'spatial_premise_closure_audit.md'
assert '## 59.' not in texts[audit]
texts[audit] += '''
## 59. 第431轮：相同持续交换作用中的内部记录转移

2026-09-24。430已有实质关系动力学结果，本轮沿用户要求继续；没有将瞬时可测概率直接称为已生成记录器。对象为430的三个qubit，读者为内部singlet对；H=S₁₂＋S₂₃＋S₃₄＋S₄₅固定且一直开启。初态、四条参与关系及强度单位分别列为输入。

两种原三体态具有相同完整两体边缘，读者初态相同。完整持续作用使读者singlet概率不同；差的首项为−√3 t⁵／4。整数置换递推给精确嵌套对易系数，80阶有理多项式配解析余项界，在t=3／2证明读者迹距离大于2／3。独立矩阵指数给约0.6965511。由效果交换子范数√3／2，整个[17／12,19／12]窗口的迹距离大于1／2，等先验理论单次辨识成功率大于3／4，不要求精确停机或读取瞬间。

所有单qubit边缘始终为I／2，差异仅在读者两体关系中。整体未知输入与参考保持于全局酉演化，但原对象边缘会受扰；没有无扰克隆、后选择或永久存储声明。逆酉只核信息保存，不要求宇宙回滚。当前补齐内部信息载体与有限窗口，宏观放大读出、物理钟、资源自发准备和持续复用仍待实现。

6项新检查、10个编号公式。累计2014项，604份编号科学文件、639份保护证据；历史证据哈希完整，新JSON只读复现，未做图像检查。精确整数证书和独立谱指数共同核验，不声称本轮已由其他代理终审。

下一项沿同一交换分支核实际参与结构及初态来源，先检验完全均匀多体交换的关系自由度，再判断需要什么最小内部结构。不再重复指定链上的峰值优化或把给定图维数当空间推导。三维、阶段收尾和完整GR均未完成，目标保持活动。
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
