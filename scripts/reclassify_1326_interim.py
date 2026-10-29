"""Move the withdrawn manuscript and record reversible editorial changes once."""
from pathlib import Path
import difflib
import hashlib
import json
import os
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
RESEARCH = ROOT / 'research_cognition_physics'
STAGE = RESEARCH / 'archive_1323_'
OLD = RESEARCH / '16_未知参考交接的作用选择与动态来源_阶段论文.md'
NEW = STAGE / '未知参考交接的作用选择与动态来源_中间成果.md'
MIG = RESEARCH / '_migration/reclassify_1326_interim'
MANIFEST = MIG / 'manifest.json'
PROPOSAL = STAGE / '_shared/next_long_horizon_proposal.md'


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def rel(path):
    return path.relative_to(ROOT).as_posix()


def edits(old, new):
    a, b = old.splitlines(keepends=True), new.splitlines(keepends=True)
    offsets = [0]
    for line in a:
        offsets.append(offsets[-1] + len(line))
    return [dict(start=offsets[i], end=offsets[j], before=''.join(a[i:j]), after=''.join(b[k:l]))
            for kind, i, j, k, l in difflib.SequenceMatcher(None, a, b, autojunk=False).get_opcodes()
            if kind != 'equal']


def rebase_links(text, source, destination):
    def replace(match):
        target = match.group(1)
        if '://' in target or target.startswith('#'):
            return match.group(0)
        name, marker, anchor = target.partition('#')
        resolved = (source.parent / name).resolve()
        if resolved == OLD:
            resolved = NEW
        mapped = Path(os.path.relpath(resolved, destination.parent)).as_posix()
        return '](' + mapped + (marker + anchor if marker else '') + ')'
    return re.sub(r'\]\(([^)\s]+)\)', replace, text)


def prepare():
    assert not MANIFEST.exists() and OLD.is_file() and not NEW.exists()
    assert OLD.resolve().is_relative_to(ROOT) and NEW.resolve().is_relative_to(STAGE)
    assert PROPOSAL.is_file(), 'Write the proposed longer horizon first'
    changes = []

    def add(source, destination, text):
        raw = source.read_bytes()
        encoding = 'utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
        revised = text.encode(encoding)
        entry = dict(original=rel(source), destination=rel(destination),
                     original_sha256=sha(raw), current_sha256=sha(revised),
                     original_bytes=len(raw), current_bytes=len(revised),
                     edits=edits(raw.decode(encoding), text))
        changes.append((source, destination, raw, revised, entry))

    paper = rebase_links(OLD.read_text(encoding='utf8'), OLD, NEW)
    paper = paper.replace('# 未知参考交接的作用选择与动态来源', '# 未知参考交接的作用选择与动态来源：中间成果汇总', 1)
    first = paper.splitlines()[4]
    paper = paper.replace(first, '2026-10-19。第1323—1326轮中间成果。依用户要求撤回原“第16阶段论文”定位，取消编号，移入研究档案。已核验的具体结论保留，但不作为完整研究阶段收尾；现行状态见[研究目录](README.md)，后续范围见[长程目标建议](_shared/next_long_horizon_proposal.md)。', 1)
    paper = paper.replace('本文完成的是原阶段目标内的条件性作用选择及同对象来源连接。', '本文仅记录有限操作合同下的条件性作用选择及同对象来源连接。')
    paper = paper.replace('本阶段', '本工作包')
    paper = paper.replace('最终按[原五项联合标准]', '当时按[窄工作包五项标准]')
    paper = paper.replace('阶段完成只表示这五项在声明范围内交付，不表示完整标准模型或现实参数已经推出。', '这些检查仅验收局部工作包。用户现明确：它们不足以作为独立物理阶段的完成标准；后续须继续形成从认知机制到统一物质相互作用及可辨后果的完整链条。')
    paper = paper.replace('## 10. 下一实质问题', '## 10. 后续长程研究中的首项连接')
    paper += '\n场身份连接本身也只是下一目标的里程碑。不能在一个有限实现、一次源匹配或一项候选失败后再次宣布整个阶段结束；需按长程目标的联合成果判断。\n'
    add(OLD, NEW, paper)

    prefix = '\n## 当前定位：1323—1326为中间成果，完整研究阶段未收尾\n\n2026-10-19，依用户要求撤回原第16阶段论文，移入[中间成果汇总](archive_1323_/未知参考交接的作用选择与动态来源_中间成果.md)，取消正式论文编号。四轮报告、具体证明和复算结果保留；此前的窄工作包验收不能再作为完整物理阶段结项。\n\n[研究目录](archive_1323_/README.md) · [下一长程目标建议](archive_1323_/_shared/next_long_horizon_proposal.md)。建议围绕“从认知操作机制到统一物质相互作用的生成与选择”组织连续研究，将场身份、共同物质作用、结构／参数选择及可辨预测作为一条联合链，局部模型和技术接口仅作里程碑。\n\n目前尚未创建下一应用目标，下一可用研究编号1327。旧目标及验收文件保留为历史记录，其阶段定位由本次用户修正覆盖；没有撤销已经成立的局部数学结论。[引用迁移与历史复算](_migration/reclassify_1326_interim/README.md)。\n'
    layer = json.loads((RESEARCH / '_migration/interaction_rules_active_1323/layer.json').read_bytes())
    for entry in layer['entries']:
        path = ROOT / entry['destination']
        text = path.read_text(encoding='utf8')
        old_prefix = entry['edits'][0]['after']
        assert text.count(old_prefix) == 1
        add(path, path, text.replace(old_prefix, prefix, 1))

    readme = '# 1323起：认知操作机制与物质相互作用研究\n\n2026-10-19。研究目录保持开放。1323—1326是已核验的中间工作包，依用户要求不再独立认定为完整研究阶段；原第16论文已撤回编号并移入本目录。\n\n## 当前定位与后续范围\n\n- [中间成果汇总](未知参考交接的作用选择与动态来源_中间成果.md)\n- [下一长程目标建议](_shared/next_long_horizon_proposal.md)：尚未创建新的应用目标。\n- [共同前件与剩余自由](_shared/joint_premise_ledger_after1326.md)\n- [输入与选择账](_shared/selection_ledger.md)\n\n现有工作解决了指定有限接口中的作用选择、角色和共同来源问题，尚未形成统一物质相互作用及其物理选择的完整成果。后续应把这些结果作为材料，优先检验能够改变大方向的连接和竞争；不得把单个候选的内部技术修补当作目标。\n\n## 已保存的研究轮次\n\n|轮次|中间结果|逐轮核验|\n|---|---|---|\n|[1323](research_note_1323.md)|未知参考交付的角色、归一作用与真实组合|[验收](1323/acceptance.json)|\n|[1324](research_note_1324.md)|相干有效过程、完整后态与资源反作用|[验收](1324/acceptance.json)|\n|[1325](research_note_1325.md)|双端口角色、重数关系与跨部门预算|[验收](1325/acceptance.json)|\n|[1326](research_note_1326.md)|动态参考低带、实际来源与有限预测|[验收](1326/acceptance.json)|\n\n## 历史记录与复算\n\n[原窄工作包目标](_shared/goal.md)、[当时的联合验收](_shared/stage_acceptance.json)、[数学核查](_shared/stage_math_review.md)、[范围核查](_shared/stage_scope_review.md)、[联合结果](_shared/stage_joint_results.json)保留原证据含义。文件中原“论文16／阶段完成”属于当时记录，现不具有完整研究阶段收尾的地位。下一可用编号为1327，无新增科学轮次或公理。\n\n[迁移说明](../_migration/reclassify_1326_interim/README.md)提供当前只读复算入口；旧路径与旧字节在内存中按可逆编辑核验，不维护旧论文副本。\n'
    add(STAGE / 'README.md', STAGE / 'README.md', readme)
    names = ('joint_premise_ledger_after1326.md', 'selection_ledger.md', 'stage_math_review.md', 'stage_scope_review.md')
    for name in names:
        path = STAGE / '_shared' / name
        text = path.read_text(encoding='utf8')
        text = text.replace('../../16_未知参考交接的作用选择与动态来源_阶段论文.md', '../未知参考交接的作用选择与动态来源_中间成果.md')
        text = text.replace('[论文16', '[中间成果')
        index = text.index('\n') + 1
        note = '\n> 2026-10-19用户修正：原“论文16”已撤回，现为[中间成果](../未知参考交接的作用选择与动态来源_中间成果.md)。下文原阶段验收措辞只对应当时的窄工作包，不表示完整研究阶段收尾；已核科学结论保留。表内哈希是迁移前冻结值，通过[迁移映射](../../_migration/reclassify_1326_interim/README.md)恢复旧字节复核，不冒称当前文件哈希。\n'
        text = text[:index] + note + text[index:]
        add(path, path, text)

    scientific = json.loads((STAGE / '_shared/stage_acceptance.json').read_bytes())['scientific_sha256']
    assert len(scientific) == 40
    for name, expected in scientific.items():
        assert sha((ROOT / name).read_bytes()) == expected
    manifest = dict(schema='interim_reclassification_v1', date='2026-10-19',
                    reason='User withdrew manuscript 16 and reclassified 1323-1326 as intermediate results.',
                    stage_complete=False, new_goal_started=False, new_scientific_rounds=0,
                    next_available_round=1327, paper_number_withdrawn=16,
                    entries=[c[4] for c in changes], scientific_sha256=scientific,
                    unchanged_json_and_python=True, previous_goal_status='historically completed; no active goal now',
                    current_goal_proposal=rel(PROPOSAL), historical_layer_preserved=True)
    return changes, manifest


def main():
    changes, manifest = prepare()
    if '--apply' not in sys.argv:
        print(json.dumps(dict(status='READY', changed_markdown=len(changes), paper_destination=rel(NEW)), ensure_ascii=False))
        return
    written = []
    try:
        for source, destination, before, after, entry in changes:
            assert source.read_bytes() == before
            if source != destination:
                assert not destination.exists()
                source.rename(destination)
            written.append((source, destination, before, after))
            destination.write_bytes(after)
        MIG.mkdir(parents=True, exist_ok=True)
        with MANIFEST.open('x', encoding='utf8') as output:
            json.dump(manifest, output, ensure_ascii=False, indent=2)
            output.write('\n')
    except BaseException:
        for source, destination, before, after in reversed(written):
            if destination.is_file() and destination.read_bytes() in (before, after):
                destination.write_bytes(before)
                if source != destination and not source.exists():
                    destination.rename(source)
        raise
    print(json.dumps(dict(status='APPLIED', changed_markdown=len(changes), scientific_files_unchanged=40,
                          paper_destination=rel(NEW), stage_complete=False, new_goal_started=False), ensure_ascii=False))


if __name__ == '__main__':
    main()
