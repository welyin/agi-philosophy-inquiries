"""One-time, reversible insertion of the SR-stage navigation; no historical edits."""
from pathlib import Path
import hashlib
import json
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import relocation_original_bytes

manifest = HERE/'navigation_layer.json'
if manifest.exists():
    raise SystemExit('Already applied; use the recorded layer instead of inserting again.')

content = '''## 当前阶段：狭义相对论的实际钟尺与信号（1086起）

**三维阶段已封存为[archive_1063_1085](archive_1063_1085/README.md)，后续研究使用[archive_1086_](archive_1086_/README.md)。** 已正式采用的CO1—CO5及[三维阶段论文](认知操作五公理与三维关系空间_阶段论文.md)保持原科学结论。

当前目标是在同一认知操作模型和声明有效域中，连接实际内部钟、三维位置、真实惯性观察者、共同信号与钟尺物质记录，判断是否迫使Lorentz变换。[新目标回执](archive_1086_/_shared/goal_start.json)已保存；目标进行中。整体ROADMAP继续，当前阶段不要求同时完成引力、标准模型或无限高能模型。

四项桥为：**实际时间计量；运动状态之间原任务的双向接续；由实际信号实现的共同不变速度；钟尺、物质及记录的共同变换。** 294、575、963、1041按原范围复用，见[时钟审计](archive_1086_/_admission/clock_reuse.md)和[运动学查重](archive_1086_/_admission/kinematics_reuse.md)。平直有效域、时空齐性、方向等价、连续组合、同步与单位也明确登记，不由三维或有限传播上限直接假定Lorentz对称。

最新正式[1086](archive_1086_/research_note_1086.md)：同一有效钟尺模型的回声仍允许尺度$\\lambda$，而同一对固定单位钟的双向滴答间隔比为$\\lambda^2$。有限实例回声全过，双向滴答差仍为11/15；这给运动角色交换候选一个可验判据，尚未由CO1—CO5推出它。成熟运动学采用计0，本轮新增独立科学校准0，累计保持3860。[解析与复算](archive_1086_/1086/proof.md) · [验收与审阅](archive_1086_/1086/acceptance.json) · [下一项](archive_1086_/1086/NEXT.md)。

报告放在阶段根目录，代码、结果和证明在各轮编号子目录。数学排版采用块级 `$$`、行内 `$`；历史代码／JSON及原核验事实保留，目录与导航变化使用可逆差分。[目录迁移说明](_migration/closure_1063_1085_20261010/README.md)。

**以下是此前阶段的历史导航，旧“当前”、complete及“不自动启动”均指写入当时，不覆盖本节。**

'''
entries = []
for name in ['README.md','research_direction.md','RESEARCH_STATE.md','ROADMAP.md']:
    p = ROOT/'research_cognition_physics'/name
    before = p.read_bytes()
    enc = 'utf-8-sig' if before.startswith(b'\xef\xbb\xbf') else 'utf-8'
    text = before.decode(enc)
    pos = text.index('\n')+1
    nl = '\r\n' if '\r\n' in text else '\n'
    addition = nl+content.replace('\n',nl)
    after = (text[:pos]+addition+text[pos:]).encode(enc)
    rel = p.relative_to(ROOT).as_posix()
    e = {'original':rel,'destination':rel,
         'original_sha256':hashlib.sha256(before).hexdigest(),
         'current_sha256':hashlib.sha256(after).hexdigest(),
         'original_bytes':len(before),'current_bytes':len(after),
         'edits':[{'start':pos,'end':pos,'before':'','after':addition,'kind':'sr_current_navigation'}]}
    assert relocation_original_bytes(after,e) == before
    entries.append(e)
    p.write_bytes(after)
manifest.write_text(json.dumps({'schema':'research_sr_navigation_layer_v1','date':'2026-10-10',
    'purpose':'Current SR-stage navigation, reversible to spatial-closure bytes before math normalization.',
    'entries':entries},ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'navigation_files':len(entries),'reversible':True}))
