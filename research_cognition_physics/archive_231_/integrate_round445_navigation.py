"""Guarded navigation update for autonomous reciprocal-address rematching."""
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASE, ROOT = HERE.parent, HERE.parent.parent
PATHS = [ROOT/'README.md', BASE/'README.md', BASE/'research_direction.md',
         BASE/'RESEARCH_STATE.md', HERE/'README.md', HERE/'spatial_premise_closure_audit.md']
original = {p: p.read_bytes() for p in PATHS}
texts = {p: raw.decode('utf-8').replace('\r\n', '\n') for p, raw in original.items()}
assert '**第445轮完成：**' not in texts[BASE/'RESEARCH_STATE.md']
summary = ('**第445轮完成：** [互指一致性与自主伙伴重配](%sresearch_note_445.md)'
           '普通两主体地址SWAP配合明示一致性罚能，通过非一致中间态给完美匹配重组；'
           '任意偶数N的二阶项为−N²I／8−A_match／2。'
           '四主体模型对全部未知代码相干及参考给误差证书；'
           '指定初始匹配在整段窗口的实际重配概率严格大于0.78。'
           '该作用允许跨旧分量，不自动选择440已有路径条件。'
           '6项检查，累计2101项；646份编号科学文件、681份保护证据。'
           '一致性罚能、地址接触及初态仍为输入，未生成三维。')
for path in PATHS[:5]:
    prefix = ('research_cognition_physics/archive_231_/' if path == ROOT/'README.md'
              else ('' if path == HERE/'README.md' else 'archive_231_/'))
    old = next(line for line in texts[path].splitlines() if '**第444轮完成：**' in line)
    texts[path] = texts[path].replace(old, old+'\n\n'+('> ' if path == ROOT/'README.md' else '')+summary % prefix, 1)
    texts[path] = texts[path].replace('231—444轮', '231—445轮')
next_item = (
    '**下一项：完整主体交换怎样同时搬运数据与改变关系。** 445已给地址重配，'
    '不继续调整ε、窗口或精度。检验将429的交换施于主体完整地址与数据块时，'
    '同一互指模型是否直接提供实际主体数据传递；保留未知联合输入及参考，'
    '不另加按关系控制的数据SWAP。若存在精确条件置换化简，须同步变换真实读数，'
    '区分表示随地址移动和固定主体实际收到信息，不重做436通用模拟。'
    '一致性罚能、全对接触与初始代码仍为输入；空间局域选择、测距、三维及新主体接入继续开放。')
for path in (BASE/'research_direction.md', BASE/'RESEARCH_STATE.md', HERE/'README.md'):
    old = next(line for line in texts[path].splitlines() if line.startswith('**下一项：直接交换怎样维持伙伴一致性。**'))
    texts[path] = texts[path].replace(old, next_item+
        '\n\n**444后互指任务（445已补普通交换的自主重配与边界）：** '+old.split('** ', 1)[1], 1)
texts[BASE/'research_direction.md'] = texts[BASE/'research_direction.md'].replace(
    '最新科学轮次与检查数为444／2095', '最新科学轮次与检查数为445／2101')
state = BASE/'RESEARCH_STATE.md'
texts[state] = texts[state].replace('最新444轮及累计2095项见本文开头', '最新445轮及累计2101项见本文开头')
readme = BASE/'README.md'
texts[readme] = texts[readme].replace('当前复算与冻结入口：',
    '当前复算与冻结入口：[445科学核验](archive_231_/verify_reciprocal_exchange_dynamics_round.py)、[445整合核验](archive_231_/verify_round445_integration.py)；历史入口：', 1)
texts[readme] = texts[readme].replace('第三阶段累计2095项；旧科学证据保持原字节', '第三阶段累计2101项；旧科学证据保持原字节')
archive = HERE/'README.md'
texts[archive] = texts[archive].replace('当前完成444轮', '当前完成445轮').replace('当前444不作为预定终点', '当前445不作为预定终点')
row = next(line for line in texts[archive].splitlines() if line.startswith('| [444：'))
texts[archive] = texts[archive].replace(row, row+
    '\n| [445：互指一致性与自主伙伴重配](research_note_445.md) | 两主体虚过程；全N匹配flip；裸代码误差；整段窗口概率；路径选择反例 | [代码](reciprocal_exchange_dynamics_audit.py)、[结果](reciprocal_exchange_dynamics_audit_results.json)、[核验](research_round_445_checks.json)；6项检查；罚能和全对接触仍为输入 |', 1)
texts[archive] = texts[archive].replace('当前整合入口：',
    '当前整合入口：[445整合核验](verify_round445_integration.py)、[整合记录](round445_integration_checks.json)；历史入口：', 1)
texts[archive] = texts[archive].replace(
    '当前2095项科学检查；678份保护证据，其中本阶段编号科学文件643份',
    '当前2101项科学检查；681份保护证据，其中本阶段编号科学文件646份')
audit = HERE/'spatial_premise_closure_audit.md'
assert '## 73.' not in texts[audit]
texts[audit] += '''
## 73. 第445轮：连续两主体交换、互指稳定与重配范围

2026-09-24。沿444明确下一项，检验短暂非一致中间态，不再重做434—436虚过程与通用模拟存在性，也不再发现438的占用度／潜在作用区分。匹配配置图对接Cioabă、Royle、Tan（2021）的已知完全图完美匹配flip图。

偶数N、每主体一份N维地址，互指效果n_ij为明确两主体算符；H＝Δ(NI−2Σn_ij)＋gΣS_addr，初态取互指完美匹配。罚能正定，代码外能隙2Δ；全H仅含两主体作用，每地址出现次数严格守恒。零能代码在N!维置换部门，原始重复／自指态仍保留在完整张量定义，不通过后选删去。

同旧边两端交换的中间代价2Δ，跨两旧边交换的代价4Δ。一般偶数N的二阶项为(g²／Δ)(−N²I／8−A_match／2)，每种替代配对有两条虚路径。N＝2、4、6、8全部匹配逐列用实际两次交换及独立换边算法核对。初态代码的所有时间码外概率≤min(1,m²ε²)，m＝N(N−1)／2；该界只控制已准备代码，不代表冷却、修复任意不一致初态或无限规模统一稳定。

N＝4的24×3精确分数Gram给裸代码Duhamel误差B＝(2＋4|τ|)|ε|＋(13／5＋4|τ|)ε²，覆盖任意未知匹配相干／参考，无需制备修饰代码。ε＝1／128、τ∈[1.9,2.1]时B<0.082。针对任一指定初始匹配，其他两匹配的实际未后选总概率，由sin凹性和端点交错多项式给严格有理下界0.78809…>0.78。不能推广为任意未知相干态都改变；均匀匹配态为有效本征态。完整256维谱指数独立核对，未用浮点值代替证书。

C＝1的旧匹配没有任何三路径，440的flip为零；本轮重配却非零并跨原分量。因此这份一致性机制不选择旧路径局域性。每主体虽只占一边，其二阶非对角作用范数仍为g²(N−2)／(2Δ)。地址字典、共同身份比较、罚能、全对接触及初态均明确为新增输入，429没有推出H₀；数据本轮仍是旁观因子。没有光锥、三维或GR的生成结论。

下一项检验完整主体地址＋数据SWAP是否同时给伙伴变化与真实数据接口，先求精确条件置换并变换固定主体读数；不另输入受关系控制的数据项，不优化本轮小模型参数。6项检查、12公式，累计2101项、646份编号科学文件、681份保护证据。旧轮次原字节保留；阶段未结项。
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
