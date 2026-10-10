"""Publish 1099 through the established reversible navigation layer."""
from pathlib import Path
from difflib import SequenceMatcher
import hashlib
import json
import subprocess
import sys

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
ROOT = HERE.parents[2]
BASE = ROOT/'research_cognition_physics'
sys.path.insert(0, str(ROOT/'scripts'))
from research_layout import relocation_original_bytes


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def run(path, *args):
    p = subprocess.run([sys.executable, '-B', '-X', 'utf8', str(path), *args],
                       capture_output=True, text=True, encoding='utf8')
    assert p.returncode == 0, p.stdout+p.stderr
    return json.loads(p.stdout)


def main():
    acceptance = HERE/'acceptance.json'
    assert not acceptance.exists(), 'Already published; refusing overwrite.'
    reviews = json.loads((HERE/'review_evidence.json').read_text(encoding='utf8'))['reviews']
    assert len(reviews) == 2
    hashes = reviews[0]['source_sha256']
    assert len(hashes) == 5
    for r in reviews:
        assert r['status'] == 'PASS_CONDITIONAL_STORAGE_FRAME_BRIDGE'
        assert r['source_sha256'] == hashes and r['actual_recomputation']['groups'] == 6
    for name, value in hashes.items():
        assert sha((STAGE/name).read_bytes()) == value, name
    check = run(HERE/'check.py')
    assert check['status'] == 'PASS' and check['groups'] == 6
    lp = BASE/'_migration/active_documents/layer.json'
    original = lp.read_bytes()
    layer = json.loads(original)
    common = '''## 当前目标：验证可扩展共识与共同洛伦兹结构（1092起）

目标active，[猜想](../猜想/可扩展共识与共同洛伦兹结构猜想.md)继续检验，目标内容未改。F＋U＋C＋P、CO1—CO5、六条任务限定协议及H1—H3保留，G未采用。[六协议范围](archive_1086_/_admission/six_protocols_and_common_motion.md)与原稿物理推论箭头区分。

**最新正式1099，累计3860。** [本轮报告](archive_1086_/research_note_1099.md)得到一条绕开1093短程分解D4的条件路线：同一完整任务中的趋零时长真实保留H、同标定同背景实际参考群R、完整空间旋转与真实运动O，以及全资源非退化基线D，共同迫使正有限任务速度确界，再复用1097得到共形Lorentz形式。[证明](archive_1086_/1099/proof.md) · [6组结果](archive_1086_/1099/results.json) · [验收](archive_1086_/1099/acceptance.json)。新公理0，科学计0。

证明通过旋转的有限共轭词取消时间尺度；若参考不混合时间与位置，实际保留会造成固定距离等待下界趋零。无需中间单人恢复、D2实际输运串接或D4短程替代。H只需一列实际时长，不预设微观连续。完整群词实际资格和全域适用仍是公开的额外条件。

H/R/O/D均未由六协议证明，也未采用为新公理。随身保持是否已完成原任务必须明确；固定接收口等待界不能代替覆盖移动输出的D。[下一步](archive_1086_/1099/NEXT.md)优先认证同一个任务中的保留、空间分离及实际参考。1098的全部影响连接、实际尺度与全部物质责任仍开放，整个加强猜想及共同Lorentz物理未判定。

研究在[archive_1086_](archive_1086_/README.md)，报告在阶段根目录、资产在编号目录。[目标启动回执](archive_1086_/_shared/consensus_goal_start_20261010.json)保持；以下1086—1091是旧目标的限定反证，不覆盖后来加强的猜想。

'''
    stage = '''## 当前目标：可扩展共识与共同洛伦兹结构猜想

目标active，[启动回执](_shared/consensus_goal_start_20261010.json)保留。最新正式[1099](research_note_1099.md)、累计3860，新公理0。四原则、五公理、六协议及H1—H3保留，G未采用。

[1099证明](1099/proof.md) · [结果](1099/results.json) · [验收](1099/acceptance.json) · [下一步](1099/NEXT.md)。新增条件定理：同一任务的短时保留H、同背景同标定实际参考群R、完整旋转与非静止运动O，以及全资源基线非退化D，直接给正有限任务速度，再接1097共形Lorentz。无需1093的D2与D4。6组精确算术校准及两份只读独立审阅通过；候选条件未采用，无实际仪器数据。

同体保持与固定接收装置收件必须区分；短时完成及覆盖移动输出的D须在同一实际任务认证。矩阵可乘、内部换轴或背景整体变化不认证参考权限。H仅要求时长序列，不预设微观连续。后继优先关闭H/R/O/D的认知来源；[1098](research_note_1098.md)的全影响合同、钟尺尺度及全部物质共同运动仍开放。

'''
    addition = '''
- [1099报告](research_note_1099.md)：短时保留、非退化距离及实际参考的条件性Lorentz选择，绕开D2/D4。
- [1099证明](1099/proof.md)、[代码](1099/check.py)、[结果](1099/results.json)、[验收](1099/acceptance.json)、[独立审阅](1099/review_evidence.json)、[核验](1099/verify.py)与[下一步](1099/NEXT.md)。
'''
    roots = {f'research_cognition_physics/{n}' for n in
             ['README.md', 'research_direction.md', 'RESEARCH_STATE.md', 'ROADMAP.md']}
    changes = []
    for entry in layer['entries']:
        name = entry['destination']
        p = ROOT/name
        old = p.read_bytes()
        base = relocation_original_bytes(old, entry)
        enc = 'utf-8-sig' if old.startswith(b'\xef\xbb\xbf') else 'utf8'
        content = old.decode(enc)
        nl = '\r\n' if '\r\n' in content else '\n'
        new = content
        if name in roots:
            start = content.index('## 当前目标：')
            end = content.index('## 狭义相对论连接', start)
            new = content[:start]+common.replace('\n', nl)+content[end:]
        elif name.endswith('archive_1086_/README.md'):
            start = content.index('## 当前目标：')
            end = content.index('## 阶段目标', start)
            new = content[:start]+stage.replace('\n', nl)+content[end:]
        elif name.endswith('archive_1086_/文件索引.md'):
            pos = content.index('\n')+1
            new = content[:pos]+addition.replace('\n', nl)+content[pos:]
        if new == content:
            continue
        a = base.decode(enc).splitlines(keepends=True)
        b = new.splitlines(keepends=True)
        offsets = [0]
        for line in a:
            offsets.append(offsets[-1]+len(line))
        edits = []
        for tag, i, j, k, l in SequenceMatcher(None, a, b, autojunk=False).get_opcodes():
            if tag != 'equal':
                edits.append(dict(start=offsets[i], end=offsets[j], before=''.join(a[i:j]),
                    after=''.join(b[k:l]), kind='round1099_storage_frame_bridge'))
        raw = new.encode(enc)
        entry.update(current_sha256=sha(raw), current_bytes=len(raw), edits=edits)
        assert relocation_original_bytes(raw, entry) == base
        changes.append((p, old, raw))
    assert len(changes) == 6
    assets = dict(hashes)
    for name in ['verify.py', 'publish.py', 'review_evidence.json']:
        assets['1099/'+name] = sha((HERE/name).read_bytes())
    prior = {}
    for name in ['archive_1086_/1098/acceptance.json', 'archive_1086_/1098/proof.md',
                 'archive_1086_/1097/proof.md', 'archive_1086_/1093/proof.md',
                 'archive_370_428/research_note_402.md',
                 'archive_1086_/_admission/six_protocols_and_common_motion.md']:
        p = BASE/name
        prior[p.relative_to(ROOT).as_posix()] = sha(p.read_bytes())
    evidence = dict(round=1099, date='2026-10-10', status='PASS_CONDITIONAL_STORAGE_FRAME_BRIDGE',
        conditional_finite_speed_derived=True, conditional_conformal_lorentz=True,
        D2_required=False, D4_required=False, new_H_R_O_D_adopted=False,
        actual_task_contract_certified=False, entire_conjecture_decided=False, Lorentz_derived=False,
        all_physical_influences_certified=False, physical_clock_scale_fixed=False,
        all_matter_motion_transport_certified=False, new_adopted_axioms=0,
        scientific_count_increment=0, scientific_count_total=3860,
        assets_sha256=assets, prior_sources_sha256=prior, independent_read_only_reviews=2,
        review_delivery='Root records actual read-only agent messages; no reviewer-authored files.',
        goal_status='active')
    assert lp.read_bytes() == original
    for p, old, _ in changes:
        assert p.read_bytes() == old, str(p)
    acceptance.write_text(json.dumps(evidence, ensure_ascii=False, indent=2)+'\n', encoding='utf8')
    for p, _, raw in changes:
        p.write_bytes(raw)
    lp.write_text(json.dumps(layer, ensure_ascii=False, indent=2)+'\n', encoding='utf8')
    verification = run(HERE/'verify.py')
    layout = run(ROOT/'scripts/run_research_active.py', '--verify-layout')
    previous = run(STAGE/'1098/verify.py')
    receipt = dict(round=1099, date='2026-10-10', status='PASS', goal_status='active',
        goal_content_changed=False, entire_conjecture_decided=False,
        round_verification=verification, layout_verification=layout, previous_round_verification=previous,
        navigation_documents_updated=6, actual_new_review_groups=6, new_scientific_count=0,
        stage_paper_sha256=sha((BASE/'认知操作五公理与三维关系空间_阶段论文.md').read_bytes()))
    (HERE/'publication_checks.json').write_text(
        json.dumps(receipt, ensure_ascii=False, indent=2)+'\n', encoding='utf8')
    p = BASE/'_migration/active_documents/checks.json'
    maintenance = json.loads(p.read_text(encoding='utf8'))
    maintenance['latest_maintenance'] = dict(date='2026-10-10',
        scope='Round1099 conditional short storage and actual reference bridge',
        active_documents=7, layer_sha256=sha(lp.read_bytes()), layout_verified=True, goal_status='active',
        latest_formal_round=1099, new_formal_rounds=1,
        receipt='../../archive_1086_/1099/publication_checks.json')
    p.write_text(json.dumps(maintenance, ensure_ascii=False, indent=2)+'\n', encoding='utf8')
    print(json.dumps(dict(round=1099, status='PASS', navigation_documents_updated=6,
        independent_reviews=2, recomputation_groups=6, layout_passed=True,
        previous_round_passed=True, goal_status='active'), ensure_ascii=False))


if __name__ == '__main__':
    main()
