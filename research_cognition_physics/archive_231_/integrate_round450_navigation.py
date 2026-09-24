"""Guarded navigation update for round 450 composite agreement and fixed comparison audit."""
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASE, ROOT = HERE.parent, HERE.parent.parent
PATHS = [ROOT/'README.md', BASE/'README.md', BASE/'research_direction.md',
         BASE/'RESEARCH_STATE.md', HERE/'README.md', HERE/'spatial_premise_closure_audit.md']
original = {p: p.read_bytes() for p in PATHS}
texts = {p: raw.decode('utf-8').replace('\r\n', '\n') for p, raw in original.items()}
assert '**第450轮完成：**' not in texts[BASE/'RESEARCH_STATE.md']
summary = '**第450轮完成：** [相同状态交互规则在复合主体上的继承边界](%sresearch_note_450.md)证明弱“继续一致”等价于H与完整交换对易，并在对应子部件交换组合下闭合；独立两Bell主体经微观交换后仍一致，却各自距原态3／4，故429强不更新合同不自动继承。原互指势在数值字典下不满足弱合同，双主体改变固定识别可恢复相容；多主体识别来源仍开放。6项检查，累计2133项；661份编号科学文件、696份保护证据。未改写429或采纳新公理，三维未生成。'
for path in PATHS[:5]:
    prefix = ('research_cognition_physics/archive_231_/' if path == ROOT/'README.md'
              else ('' if path == HERE/'README.md' else 'archive_231_/'))
    old = next(line for line in texts[path].splitlines() if '**第449轮完成：**' in line)
    texts[path] = texts[path].replace(old, old+'\n\n'+('> ' if path == ROOT/'README.md' else '')+summary % prefix, 1)
    texts[path] = texts[path].replace('231—449轮', '231—450轮')
next_item = '**下一项：多主体比较识别与自然交互的整体相容性。** 450已区分强不更新与弱继续一致，给独立复合主体的继承反例及双主体扭转识别。先回查85／88／177／212和221参考、组合结果，再核当前地址载体中各对固定识别能否来自同一组内部表示，同时与完整交换和互指作用相容。不能逐对任选F后就宣称整体统一，也不重复290环传播、一般参考存在性或434—436模拟。新增轮次须有全体相容条件、反例或减少来源输入；若仍只能任选字典与势，应记录未选定之处。实际定位、三维与GR仍开放。'
for path in (BASE/'research_direction.md', BASE/'RESEARCH_STATE.md', HERE/'README.md'):
    old = next(line for line in texts[path].splitlines() if line.startswith('**下一项：组合主体的交互合同与记录作用来源。**'))
    texts[path] = texts[path].replace(old, next_item+
        '\n\n**449后组合合同任务（450已核继承范围与识别依赖）：** '+old.split('** ', 1)[1], 1)
texts[BASE/'research_direction.md'] = texts[BASE/'research_direction.md'].replace(
    '最新科学轮次与检查数为449／2127', '最新科学轮次与检查数为450／2133')
state = BASE/'RESEARCH_STATE.md'
texts[state] = texts[state].replace('最新449轮及累计2127项见本文开头', '最新450轮及累计2133项见本文开头')
readme = BASE/'README.md'
texts[readme] = texts[readme].replace('当前复算与冻结入口：',
    '当前复算与冻结入口：[450科学核验](archive_231_/verify_composite_agreement_round.py)、[450整合核验](archive_231_/verify_round450_integration.py)；历史入口：', 1)
texts[readme] = texts[readme].replace('第三阶段累计2127项；旧科学证据保持原字节', '第三阶段累计2133项；旧科学证据保持原字节')
archive = HERE/'README.md'
texts[archive] = texts[archive].replace('当前完成449轮', '当前完成450轮').replace('当前449不作为预定终点', '当前450不作为预定终点')
row = next(line for line in texts[archive].splitlines() if line.startswith('| [449：'))
texts[archive] = texts[archive].replace(row, row+
    '\n| [450：相同状态交互规则在复合主体上的继承边界](research_note_450.md) | 弱合同充要条件；独立复合关联反例；微观支撑边界；比较字典 | [代码](composite_agreement_audit.py)、[结果](composite_agreement_audit_results.json)、[核验](research_round_450_checks.json)；6项检查；未改写429或采纳弱合同 |', 1)
texts[archive] = texts[archive].replace('当前整合入口：',
    '当前整合入口：[450整合核验](verify_round450_integration.py)、[整合记录](round450_integration_checks.json)；历史入口：', 1)
texts[archive] = texts[archive].replace(
    '当前2127项科学检查；693份保护证据，其中本阶段编号科学文件658份',
    '当前2133项科学检查；696份保护证据，其中本阶段编号科学文件661份')
audit = HERE/'spatial_premise_closure_audit.md'
assert '## 78.' not in texts[audit]
texts[audit] += '\n## 78. 第450轮：完整复合主体的认同合同与比较字典\n\n2026-09-24。接续449作用来源。429强合同F要求独立同态输入的双方原态均不更新；弱合同E只要求输出边缘彼此相同。有限同维、固定识别、封闭静态H中，E恰为[H,S]=0。H的交换反对称部分A定义Φ(X)=Tr_B[A(I⊗X)]，HS反对称性与初始边缘差−2i[Φ(ρ),ρ]给Φ(I)=0、TrΦ(X)=0及[Φ(X),X]=0。逐秩一投影、任意补基的对角性迫使Φ=0，故A=0。只要求全部同态输入的初始边缘导数相等已足够。d2/3精确约束秩6/36，弱H维数10/45而非强合同2维。\n\n组合主结果：S_block=乘积S_r，Σg_rS_r与其对易，因此对任意内部关联的独立完整ρ⊗ρ保持E。两主体各含两qubit且各自为Bell，初始整体严格独立。两条微观交换产生rho_A(t)=rho_B(t)=cos²(2t)Bell＋sin²(2t)I4/4；tπ/4旧态迹距离3/4，全部单qubit边缘不变。区别于429§7的初始跨主体关联，也不恢复429取消的不含内部标记裸更新分类。微观429不自动推出复合强F。\n\n附论限同一原始空间、所有时间精确强F：H=aI+bS_block，任何少于2k微观因子的静态项之和与全无迹张量见证正交，完整SWAP却不正交，故b0。没有排除某一端点SWAP、编码辅助或虚过程；434—436保留。\n\n来源接口：N2原始h_R奖励|1,blank>_0|0,blank>_1，完整物理交换后是自指、罚能2，故原数值比较下[h_R,S]非零。两相同加态经h_R在tπ/4边缘差1/2，弱化合同本身并不生成该势。但F同时翻两寄存器0/1、保持blank时，S_F=(F⊗F)S与h_R、h_A、均匀L及G均对易；对输入ρ⊗FρF†的新比较rho_A=F†rho_BF成立。F恰为自伴 involution；81维原空间核验，不把被动重命名混为物理交换，也不声称任意识别都失败。具体F及全体主体识别的一致来源未选定。\n\n6项检查、13公式，累计2133项；661份编号科学文件、696份保护证据。独立只读审核，旧科学文件保持原字节，无图像。未改写429、未自动采纳E公理。后继核多主体固定识别能否共同来自一套内部表示并与当前交互相容，先复用既有参考及组合成果，不重做传播环路或一般模拟。实际定位、三维及GR未完成。\n'
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
