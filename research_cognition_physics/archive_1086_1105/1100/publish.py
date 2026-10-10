"""Publish 1100 through the established reversible navigation layer."""
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
        assert r['status'] == 'PASS_CONDITIONAL_DEADLINE_SELECTION'
        assert r['source_sha256'] == hashes and r['actual_recomputation']['groups'] == 5
    for name, value in hashes.items():
        assert sha((STAGE/name).read_bytes()) == value, name
    check = run(HERE/'check.py')
    assert check['status'] == 'PASS' and check['groups'] == 5
    lp = BASE/'_migration/active_documents/layer.json'
    original = lp.read_bytes()
    layer = json.loads(original)
    common = '''## 当前目标：验证可扩展共识与共同洛伦兹结构（1092起）

目标active，[猜想](../猜想/可扩展共识与共同洛伦兹结构猜想.md)继续检验，目标内容未改。F＋U＋C＋P、CO1—CO5、六条任务限定协议及H1—H3保留，G未采用。[六协议范围](archive_1086_/_admission/six_protocols_and_common_motion.md)与原稿物理推论箭头区分。

**最新正式1100，累计3860。** [本轮报告](archive_1086_/research_note_1100.md)严格削弱1099：在同一实际正向时间、全域仿射字典、同背景同标定参考群R及完整旋转／真实运动O下，候选FD——某有限期限内已能完成任务，但仍有原空间目标不能完成——与正有限任务速度等价，再复用1097给共形Lorentz。[证明](archive_1086_/1100/proof.md) · [5组结果](archive_1086_/1100/results.json) · [验收](archive_1086_/1100/acceptance.json)。新公理0，科学计0。

无需1099的趋零保留H、单基线D、D2或D4；允许指定任务存在固定正固有处理间隔。若参考时间完全不混合位置，1099的有限共轭词会使每个非空固定时长截面遍及全空间，违反FD。成熟轨道及锥证明直接复用。权限相容例未认证物理器件或微观最小时间。

FD仍约束全部允许有限准备的并集；固定组织分别有限和空窗口稳定都不够。CO1的静态端点分离不自动给实际期限内不可达；[下一步](archive_1086_/1100/NEXT.md)优先检验同一期限的一份成功任务与一个原关系读口的全菜单排除界。R/O、全域仿射字典、实际单位和全部物质／影响责任仍开放，未由六协议推出完整Lorentz物理。

研究在[archive_1086_](archive_1086_/README.md)，报告在阶段根目录、资产在编号目录。[目标启动回执](archive_1086_/_shared/consensus_goal_start_20261010.json)保持；以下1086—1091是旧目标的限定反证，不覆盖后来加强的猜想。

'''
    stage = '''## 当前目标：可扩展共识与共同洛伦兹结构猜想

目标active，[启动回执](_shared/consensus_goal_start_20261010.json)保留。最新正式[1100](research_note_1100.md)、累计3860，新公理0。四原则、五公理、六协议及H1—H3保留，G未采用。

[1100证明](1100/proof.md) · [结果](1100/results.json) · [验收](1100/acceptance.json) · [下一步](1100/NEXT.md)。在同一实际正向计时、全域仿射字典及实际R/O下，FD（某有限期限的完整输出位置集非空且非全空间）与正有限任务速度等价，后接1097共形Lorentz。去掉1099的趋零保留H与单基线D，也无需D2/D4；带正固有完成间隔的权限例证明弱化真实。

5组精确复算及两份只读独立审阅通过。FD尚未由原CO空间或六协议推出；固定预算有界不保证全部资源并集有界，空的短窗口也不认证FD。后继检验一个非空期限、原关系读口及同菜单空间排除机制，不继续修微观短时输出。实际参考、单位与全部物质／影响连接仍开放；没有实际仪器数据。

'''
    addition = '''
- [1100报告](research_note_1100.md)：非空有限期限的空间区分取代H＋D，容许固定交付延迟。
- [1100证明](1100/proof.md)、[代码](1100/check.py)、[结果](1100/results.json)、[验收](1100/acceptance.json)、[独立审阅](1100/review_evidence.json)、[核验](1100/verify.py)与[下一步](1100/NEXT.md)。
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
                    after=''.join(b[k:l]), kind='round1100_nonempty_deadline_selection'))
        raw = new.encode(enc)
        entry.update(current_sha256=sha(raw), current_bytes=len(raw), edits=edits)
        assert relocation_original_bytes(raw, entry) == base
        changes.append((p, old, raw))
    assert len(changes) == 6
    assets = dict(hashes)
    for name in ['verify.py', 'publish.py', 'review_evidence.json']:
        assets['1100/'+name] = sha((HERE/name).read_bytes())
    prior = {}
    for name in ['archive_1086_/1099/acceptance.json', 'archive_1086_/1099/proof.md',
                 'archive_1086_/1099/check.py', 'archive_1086_/1097/proof.md',
                 'archive_1086_/1093/proof.md', 'archive_1086_/1098/proof.md',
                 '认知操作五公理与三维关系空间_阶段论文.md',
                 'archive_1086_/_admission/six_protocols_and_common_motion.md']:
        p = BASE/name
        prior[p.relative_to(ROOT).as_posix()] = sha(p.read_bytes())
    evidence = dict(round=1100, date='2026-10-10', status='PASS_CONDITIONAL_DEADLINE_SELECTION',
        conditional_finite_speed_derived=True, conditional_conformal_lorentz=True,
        H_required=False, D_required=False, D2_required=False, D4_required=False,
        new_FD_adopted=False, fd_speed_equivalence=True,
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
    previous = run(STAGE/'1099/verify.py')
    receipt = dict(round=1100, date='2026-10-10', status='PASS', goal_status='active',
        goal_content_changed=False, entire_conjecture_decided=False,
        round_verification=verification, layout_verification=layout, previous_round_verification=previous,
        navigation_documents_updated=6, actual_new_review_groups=5, new_scientific_count=0,
        stage_paper_sha256=sha((BASE/'认知操作五公理与三维关系空间_阶段论文.md').read_bytes()))
    (HERE/'publication_checks.json').write_text(
        json.dumps(receipt, ensure_ascii=False, indent=2)+'\n', encoding='utf8')
    p = BASE/'_migration/active_documents/checks.json'
    maintenance = json.loads(p.read_text(encoding='utf8'))
    maintenance['latest_maintenance'] = dict(date='2026-10-10',
        scope='Round1100 nonempty deadline spatial distinction replaces zero-duration storage',
        active_documents=7, layer_sha256=sha(lp.read_bytes()), layout_verified=True, goal_status='active',
        latest_formal_round=1100, new_formal_rounds=1,
        receipt='../../archive_1086_/1100/publication_checks.json')
    p.write_text(json.dumps(maintenance, ensure_ascii=False, indent=2)+'\n', encoding='utf8')
    print(json.dumps(dict(round=1100, status='PASS', navigation_documents_updated=6,
        independent_reviews=2, recomputation_groups=5, layout_passed=True,
        previous_round_passed=True, goal_status='active'), ensure_ascii=False))


if __name__ == '__main__':
    main()
