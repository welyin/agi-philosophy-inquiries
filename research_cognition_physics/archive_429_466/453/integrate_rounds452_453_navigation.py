"""Guarded integration of two independently researched rounds based on round 451."""
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASE, ROOT = HERE.parent, HERE.parent.parent
PATHS = [ROOT/'README.md', BASE/'README.md', BASE/'research_direction.md',
         BASE/'RESEARCH_STATE.md', HERE/'README.md', HERE/'spatial_premise_closure_audit.md']
original = {p: p.read_bytes() for p in PATHS}
texts = {p: raw.decode('utf-8').replace('\r\n', '\n') for p, raw in original.items()}
assert '**第452轮完成：**' not in texts[BASE/'RESEARCH_STATE.md']
summaries = [
    '**第452轮完成：** [宏观接口的关系相干与可用性边界](%sresearch_note_452.md)复用449的四主体48维有效模型，在保留现有物理读数的酉演化不变含幺星子代数合同下，证明最小接口恰在关系重组与内部处理均非零时为完整M48(C)。相同经典关系摘要的两态在原始384维系统、整个指定时间窗内产生大于1／10的读数差；同时给不同整体态的全部单次等待读数永久相同的反例，故代数完整不等于实际层析或控制。6项检查，累计2145项；667份编号科学文件、702份保护证据。没有完成宏观主体或三维。',
    '**第453轮完成：** [两个关系主体的组合与旧接口保持边界](%sresearch_note_453.md)独立从451基线研究两个三qubit关系编码，证明任意静态两体交换权重下，全程精确保留原联合编码当且仅当九条跨块权重相同，而这时没有跨块逻辑交互。给被省略内部变量改变实际读数的整窗见证，以及允许中途离开旧编码、端点精确交换两个未知逻辑态的正例。6项检查，累计2151项；670份编号科学文件、705份保护证据。结论限于固定旧编码及静态交换，未排除扩大共同接口、已有有效编码或SoCA。'
]
for path in PATHS[:5]:
    prefix = ('research_cognition_physics/archive_231_/' if path == ROOT/'README.md'
              else ('' if path == HERE/'README.md' else 'archive_231_/'))
    old = next(line for line in texts[path].splitlines() if '**第451轮完成：**' in line)
    added = '\n\n'.join(('> ' if path == ROOT/'README.md' else '')+s % prefix for s in summaries)
    texts[path] = texts[path].replace(old, old+'\n\n'+added, 1).replace('231—451轮', '231—453轮')

next_item = '**下一项：宏观主体的可用组合接口。** 452证明当前接口不能删除仍会影响未来物理读数的关系相干；453证明两个固定三qubit编码在静态交换下若全程零误差保持原编码，就没有跨块逻辑交互。下一项检验允许共同接口扩展时，哪些原有未知信息能在同一自然演化中保留并继续参与实际读出和再次组合；共同结构中的变量归属、准备来源与读取能力必须显式给出。优先复用431、434—436已有有效编码及误差界，与扩大关系接口作对照，不把旧主体永久保持独立编码新增为认知公理。若采用近似接口，须给未知参考一致的误差、时间窗和内部资源账；仅计算六自旋的矩阵分块、重做精确门序列、弱耦合调参或得到抽象复代数都不算新主体形成。SoCA的保真、按需共享和协调权限继续作功能候选，分别核充分性与必要性。仍沿429连续自然交互路线；三维及GR目标未完成。'
for path in (BASE/'research_direction.md', BASE/'RESEARCH_STATE.md', HERE/'README.md'):
    old = next(line for line in texts[path].splitlines() if line.startswith('**下一项：宏观复量子主体的组织与保持机制。**'))
    texts[path] = texts[path].replace(old, next_item+'\n\n**451提出的组织问题（452—453已给两类具体接口边界）：** '+old.split('** ', 1)[1], 1)

direction = BASE/'research_direction.md'
texts[direction] = texts[direction].replace('最新科学轮次与检查数为451／2139', '最新科学轮次与检查数为453／2151')
state = BASE/'RESEARCH_STATE.md'
texts[state] = texts[state].replace('最新451轮及累计2139项见本文开头', '最新453轮及累计2151项见本文开头')
readme = BASE/'README.md'
texts[readme] = texts[readme].replace('当前复算与冻结入口：', '当前复算与冻结入口：[452科学核验](archive_231_/verify_macro_interface_closure_round.py)、[453科学核验](archive_231_/verify_relational_subject_composition_round.py)、[453整合核验](archive_231_/verify_round453_integration.py)；历史入口：', 1)
texts[readme] = texts[readme].replace('第三阶段累计2139项；旧科学证据保持原字节', '第三阶段累计2151项；旧科学证据保持原字节')
archive = HERE/'README.md'
texts[archive] = texts[archive].replace('当前完成451轮', '当前完成453轮').replace('当前451不作为预定终点', '当前453不作为预定终点')
row = next(line for line in texts[archive].splitlines() if line.startswith('| [451：'))
texts[archive] = texts[archive].replace(row, row+
    '\n| [452：宏观接口的关系相干与可用性边界](research_note_452.md) | 最小不变星代数M48；实际关系相干整窗见证；代数与层析区别 | [代码](macro_interface_closure_audit.py)、[结果](macro_interface_closure_audit_results.json)、[核验](research_round_452_checks.json)；6项检查；科学基线451 |'+
    '\n| [453：两个关系主体的组合与旧接口保持边界](research_note_453.md) | 全权重泄漏Gram；零泄漏只许均匀跨块交换；实际读取与端点SWAP | [代码](relational_subject_composition_audit.py)、[结果](relational_subject_composition_audit_results.json)、[核验](research_round_453_checks.json)；6项检查；科学基线451，与452独立 |', 1)
texts[archive] = texts[archive].replace('当前整合入口：', '当前整合入口：[452整合核验](verify_round452_integration.py)、[453整合核验](verify_round453_integration.py)；历史入口：', 1)
texts[archive] = texts[archive].replace('当前2139项科学检查；699份保护证据，其中本阶段编号科学文件664份', '当前2151项科学检查；705份保护证据，其中本阶段编号科学文件670份')

audit = HERE/'spatial_premise_closure_audit.md'
assert '## 80.' not in texts[audit] and '## 81.' not in texts[audit]
texts[audit] += '''
## 80. 第452轮：关系相干不能删除，代数闭合不等于宏观主体可用

2026-09-24。科学基线451。接续用户由SoCA治理设计提出的组织问题，复用449物理编号下的四主体模型：三种匹配与四个内部二能级给48维有效码，现有四个接收空白效应Ei作为输出，B＝κK＋j hA＋ℓL。额外合同明确为包含这些效应并在同一演化下不变的含幺星子代数。证明它恰在κℓ非零时为M48(C)，j可为任意实数；κ零保留匹配中心，ℓ零保留激发数。证明使用单比特翻转迫使各数据块的对易子相同，再由单激发扇区三条独立Klein置换迫使匹配对易子为标量；有限星代数结构排除任何剩余不变真子代数。正时间解析延拓只用于数学证明，没有假设物理倒流。

两态(|M0,0001>±i|M1,1000>)/sqrt(2)有相同匹配概率、全部数据边缘态及匹配去相干后的完整摘要，未来物理E0读数却不同。高斯有理数20阶Taylor与余项界给有效模型间隔大于3／25；结合449已有未知态一致误差界，在ε＝1／1024及有效时间[31／250,63／500]上，原始384维系统读数差严格大于70856101／655360000，即大于1／10。不后选择，不把有效模型精确性扩大到原始模型；准备、读出及参数仍是建模输入。相同经典摘要预测器至少一态误差超过1／20。

另一方面，H2＝2B111，D＝H2²−2H2−24I构成非零、与B对易且对I和四个Ei均正交的方向。I／48±D／8064是不同合法态，却给所有等待时间的四种单次读数完全一致。故完整输出星代数不能冒充实际线性层析、全部状态准备或控制。这里没有排除所有不采用子代数的量子压缩，没有实现完整宏观主体。方法对接Grigoletto等2025精确量子模型约简及Bény–Richter有限代数结构，新增的是现有物理接口的参数分类和原始演化整窗证据，不重报一般代数定理。

6项检查、13公式，累计2145项；667份编号科学文件、702份保护证据。独立只读审核通过，历史科学文件不改。SoCA中保留原始差异与关系可作功能启发，但具体治理层级的必要性、实际宏观接口及三维仍未证明。

## 81. 第453轮：保持两个旧主体的固定编码，何时妨碍共同交互

2026-09-24。与452独立，科学基线同为451。复用430两个三qubit关系主体，每个固定码包含一个2维逻辑因子及一个2维被省略内部因子；联合原始空间64维、旧码16维。允许块内交换和任意实九条静态跨块交换权重J。证明全程精确保留旧联合码当且仅当所有J相等；此时跨块交换只作用于被省略内部因子，没有跨块逻辑交互。保持旧编码这一条件单独已经足够得出限制，并未把它当成认知公理。

对9P的九个交换对易子，整数Gram元素为2592δacδbd−720(δac＋δbd)＋192，秩8，零空间只有均匀权重。若J＝均值＋行偏差u＋列偏差v＋双零和矩阵D，则实际旧码泄漏对易子范数平方为16(||u||²＋||v||²)＋32||D||²，初始旧码均匀态的泄漏概率为t²[(||u||²＋||v||²)／2＋||D||²]＋O(t⁴)。这覆盖全部实权重，不是扫描几个例子。独立的压缩逻辑自治检查得到相同均匀条件。

实际64维演化中，单条跨块交换在相同初始逻辑态、不同被省略内部态下，使完整物理效应E的概率成为1／2±sin(2t)／12；|t−π／4|不超过1／8时差至少31／192。未使用PHP后选择，不能把内部因子任意抹去后宣称旧逻辑自治。边界正例为三条对齐交换：中途可离开旧码，π／2时恰给整体块SWAP，因而端点交换任意未知逻辑态及其参考。此正例未声称产生纠缠门。既有脉冲交换、434—436低能有效编码、扩大共同接口及特定内部态准备均未被排除。

对接Knill–Laflamme–Viola、van Meter–Knill及Fong–Wandzura，明确省略变量对编码门的影响与交换门设计已有文献。新增的是本研究候选的全权重零泄漏分类、实际效应整窗证据及适用范围。6项检查、13公式；两独立轮合计新增12项，累计2151项，670份编号科学文件、705份保护证据。复算与独立审核通过，旧科学文件保持原字节。

两轮共同推进的是可用宏观接口问题：仍影响后续认知的量子关系不能任意丢弃，而永远封闭在原有小主体的固定编码也可能阻止共同交互。下一步给允许扩大共同接口或有界误差时的正面可用构造，并核实际准备、读取、未知参考和再次组合，不将固定旧编码或全量私人信息共享新增为SoCA要求。它们尚不是从认知原则推导复结构、SoCA唯一性、物理三维或广义相对论的证明；阶段结项门槛未触发。
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
