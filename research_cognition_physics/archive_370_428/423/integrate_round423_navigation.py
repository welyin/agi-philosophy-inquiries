"""Guarded one-time navigation update; preserves each file's newline convention."""
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
ROOT = BASE.parent
PATHS = [ROOT/'README.md', BASE/'README.md', BASE/'research_direction.md',
         BASE/'RESEARCH_STATE.md', HERE/'README.md', HERE/'spatial_premise_closure_audit.md']
original = {p: p.read_bytes() for p in PATHS}
texts = {p: raw.decode('utf-8').replace('\r\n', '\n') for p, raw in original.items()}
assert '第423轮完成' not in texts[BASE/'RESEARCH_STATE.md'], 'Already integrated'

summary = ('**第423轮完成：** [有限测量来源的精确容量]({prefix}research_note_423.md)证明：'
           '在所有未知概率信息仅存于初始来源的明确合同下，N次全偏置精确独立二元记录最少需要N＋1维；'
           '足够小的内部偏置区间则恰需⌈√(N＋1)⌉维。有限有理证书、完整仪器及范围反例通过9项检查，累计1960项；'
           '580份编号科学文件、615份保护证据。此为413的一种来源实现审计，不否定其坐标算法，也未选择空间维数。')
for path in PATHS[:5]:
    prefix = 'research_cognition_physics/archive_231_/' if path == ROOT/'README.md' else ('' if path == HERE/'README.md' else 'archive_231_/')
    marker = next(line for line in texts[path].splitlines() if '**第421—422轮完成：**' in line)
    insertion = ('> ' if path == ROOT/'README.md' else '') + summary.format(prefix=prefix)
    texts[path] = texts[path].replace(marker, marker+'\n\n'+insertion, 1)
    texts[path] = texts[path].replace('231—422轮', '231—423轮')

direction = BASE/'research_direction.md'
texts[direction] = texts[direction].replace('最新科学轮次与检查数为422／1951', '最新科学轮次与检查数为423／1960')
texts[direction] = texts[direction].replace(
    '**下一项：实际关系的位移来源。** 421—422',
    '**下一项：实际关系的位移来源。** 423厘清固定初始来源的查询容量及偏置范围差异；该来源合同并非所有定位协议必须采用，不再沿N、区间宽度和噪声优化扩展支线。用户四种变化性质继续作为认知动力学候选，优先在主体之间可共同测量的关系上检验其作用及跨基点比较，不能把多方向直接解释为恰好三个方向。421—422')
state = BASE/'RESEARCH_STATE.md'
texts[state] = texts[state].replace('最新422轮及累计1951项见本文开头', '最新423轮及累计1960项见本文开头')

readme = BASE/'README.md'
texts[readme] = texts[readme].replace('当前复算与冻结入口：',
    '当前复算与冻结入口：[423科学核验](archive_231_/verify_finite_sampling_source_round.py)、[423整合核验](archive_231_/verify_round423_integration.py)；历史入口：', 1)
texts[readme] = texts[readme].replace('第三阶段累计1951项；旧科学证据保持原字节', '第三阶段累计1960项；旧科学证据保持原字节')

archive = HERE/'README.md'
texts[archive] = texts[archive].replace('当前完成422轮', '当前完成423轮').replace('当前422不作为预定终点', '当前423不作为预定终点')
row = next(line for line in texts[archive].splitlines() if line.startswith('| [422：'))
texts[archive] = texts[archive].replace(row, row+'\n| [423：有限测量来源的精确容量](research_note_423.md) | 全偏置N＋1与足够小局部区间⌈√(N＋1)⌉的锐界；有限有理证书及顺序仪器 | [代码](finite_sampling_source_audit.py)、[结果](finite_sampling_source_audit_results.json)、[核验](research_round_423_checks.json)；9项检查；[草稿](round423_drafts/)保留 |', 1)
texts[archive] = texts[archive].replace('当前整合入口：',
    '当前整合入口：[423整合核验](verify_round423_integration.py)、[整合记录](round423_integration_checks.json)；历史入口：', 1)
texts[archive] = texts[archive].replace('当前1951项科学检查；610份保护证据，其中本阶段编号科学文件577份，另含422',
    '当前1960项科学检查；615份保护证据，其中本阶段编号科学文件580份，另含423的两份草稿、422')
texts[archive] = texts[archive].replace('## 下一步\n',
    '## 下一步\n\n423已完成固定初始来源的有限查询审计；不继续扫描样本数、偏置窗口或误差。接回用户提出的认知变化性质：在同一内部交互过程中确定可共同查询的主体间关系，解释变化如何改变这些关系及跨基点比较，再使用已有坐标构造。多方向不预定方向数为3，查询来源维数不当作空间维数。\n', 1)
texts[archive] += '\n423以422为科学基线，只读复算：`python -B -X utf8 research_cognition_physics/archive_231_/verify_finite_sampling_source_round.py`。9项检查、14个编号公式；全偏置及足够小局部区间的不同下界不相互外推。\n'

audit = HERE/'spatial_premise_closure_audit.md'
texts[audit] += '''
## 50. 第423轮：实际查询的有限来源及认知变化猜想的接续

2026-09-24。先读取导航、421—422笔记与结果，再查385、412—413、343及相关早期记录。本回合有新的有限任务定理，不把一般概率查询、商空间定义或旧独立证据限制重新编号。

### 50.1 去重与实际新增

“以全部关系读数定义位置、检查操作能否下降”由221及378—379覆盖；343已有不同维关系网络的实际Born读数，392、397、411已经检查机制图及交互代数。因此取消通用响应向量候选。单一来源上适应读取仍是一个POVM、第一个结果复制多遍不等于独立样本，由92直接复用，也取消重做。

[423](research_note_423.md)新增明确有限任务的锐容量：所有未知p信息只在初始D维来源、后继过程不另带p、完整N位记录精确独立。全p∈(0,1)的最小D为N＋1；任意内部p_0周围足够小区间的最小D为⌈√(N＋1)⌉。来源可以消费，不增加催化返还或无限服务要求。固定来源合同是额外输入，不从认知原则直接推出。

全偏置下界对接Lee、Wei与de Wolf的PSD-rank定理24及推论26，给N＋1个严格内部有理参数的有限证书；不宣称下界方法原创。计数来源、覆盖全部源态的逐次仪器及Dicke等距给可达上界。局部构造有解析正性窗口，N＝2的三等分例在[1/4,3/4]用qubit实现，明确阻止全偏置结论误用。

### 50.2 来源容量不代替空间来源

423未假定413实际必须采用固定初始来源，也未使其方向交会失效。参数可由持续关系、实际耦合、新鲜来源等承载，须按各自过程列账。来源维数、独立样本数与空间维数不是同一个量。制备、空白记录及步数控制仍为明示资源，不称为已由自治认知演化生成。

用户再次提出连续变化、资源成本、幅度可调和方向可选，直接接续§41与415—416。可将四点作为认知动力学候选继续严格化；当前证据尚未将其提升为385的全部位移前提。变化的方向数量、主体间关系及跨基点比较未被四点指定。成本先按内部实现账理解，不追加每次变化必须净耗散正能量。无需引入认知之外的地点，也无需全部内部认知状态只有三个参数。

### 50.3 交付与停止扩展边界

9项检查、14个编号公式通过，第三阶段累计1960项；580份编号科学文件、615份保护证据，含423两份早期草稿。独立代理终审复算及原始文献核对完成。此前科学文件、结果、核验快照及前两阶段档案原字节保留，不做文章图像检查。

后继只在补充真实关系来源、关系改变或跨基点比较时增加轮次，不扫描N、局部窗口及噪声。实际三维及阶段结项未完成，原342—343在原形式前提内的不足性结论继续有效，强化内部实施范围未被本轮完全证明或反驳。待答的实施量词问题保留，不能当作正向空间研究的唯一阻碍。
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
