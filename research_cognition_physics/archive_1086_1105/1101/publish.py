"""Publish 1101 without overwriting historical sources or navigation bases."""
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
    assert not (HERE/'acceptance.json').exists(), 'Refuse to republish.'
    reviews = json.loads((HERE/'review_evidence.json').read_text(encoding='utf8'))['reviews']
    assert len(reviews) == 2
    hashes = reviews[0]['source_sha256']
    for r in reviews:
        assert r['status'] == 'PASS_CONDITIONAL_BRIDGE'
        assert r['actual_recomputation']['groups'] == 4 and r['source_sha256'] == hashes
    for name, value in hashes.items():
        assert sha((STAGE/name).read_bytes()) == value, name
    lp = BASE/'_migration/active_documents/layer.json'
    old_layer = lp.read_bytes()
    layer = json.loads(old_layer)
    common = '''## 当前目标：验证可扩展共识与共同洛伦兹结构（1092起）

目标active，[猜想](../猜想/可扩展共识与共同洛伦兹结构猜想.md)继续检验，目标内容未改。F＋U＋C＋P、CO1—CO5、六条任务限定协议及H1—H3保留，G未采用。[六协议最新审计](archive_1086_/_admission/six_protocols_after_1100.md)保留第三方一致与实际运动互换的区别。

**最新正式1101，累计3860。** [本轮报告](archive_1086_/research_note_1101.md)补出1100条件下的剩余后果：完整实际旋转与一份非静止参考的有限词，已生成全部保取向、尺度为1的Lorentz变换；没有取闭包赠送权限。定向保持的全部参考群为D×SO⁺(1,3)，D是仍可能存在的纯伸缩子群。[证明](archive_1086_/1101/proof.md) · [4组结果](archive_1086_/1101/results.json) · [验收](archive_1086_/1101/acceptance.json)。新增采用公理0、成熟数学科学计0。

新候选S：事前定义且实际运输不变的完整任务类，若在全部允许准备下具有正端点锥间隔下界，则尺度因子为1，无需实际边界光脉冲。S尚未采用或由认知推出；固定组织分别有延迟不够，资源并集可把下界压到零。锥间隔未与任意实际钟读数识别，未引入宇宙最小尺度。

[1100](archive_1086_/research_note_1100.md)的实际正向计时、全域仿射字典、R/O与FD条件继续保留；不恢复趋零保留H或D4前置。原Q任务的完整Lorentz参考覆盖不再单列缺口，但R/O/FD来源、实际单位、全部物质和可辨影响仍有独立责任。[下一步](archive_1086_/1101/NEXT.md)先审任务资格及跨资源时标的认知来源，保持上游实际运动互换与期限区分问题，避免修补处理器替代原目标。

研究在[archive_1086_](archive_1086_/README.md)，报告在阶段根目录、资产在编号目录。4组精确复算、两份独立只读审阅和归档核验通过。以下1086—1091仅为旧九项基础下的限定反证，当前加强猜想与全部Lorentz物理仍未完成。

'''
    stage = '''## 当前目标：可扩展共识与共同洛伦兹结构猜想

目标active，[启动回执](_shared/consensus_goal_start_20261010.json)保留。最新正式[1101](research_note_1101.md)，累计3860；四原则、五公理、六协议与H1—H3保持，G未采用。

[1101证明](1101/proof.md) · [结果](1101/results.json) · [验收](1101/acceptance.json) · [下一步](1101/NEXT.md)。在1100实际R/O、FD及正向仿射字典条件下，完整保单位Lorentz子群由有限参考词生成；定向保持部分为D×SO⁺(1,3)，纯伸缩D仍可存在。新增候选S若给运输不变任务类跨全部资源的正端点锥间隔下界，就消去D；不需要实际边界光脉冲，但S尚无认知来源，也未识别实际钟读数。

4组精确复算及两份独立审阅通过。新增采用公理0，成熟数学科学计0。R/O/FD、全域仿射与实际标定来源、全部物质／影响仍开放；原Q任务在这组强条件下的完整群覆盖已补齐，不再另列任意快度权限。保留[六协议补审](_admission/six_protocols_after_1100.md)，不重复旧H/D4问题、不扩建处理器。

'''
    addition = '''
- [1101报告](research_note_1101.md)：已有R/O条件生成完整Lorentz子群；分离纯伸缩自由，正任务时标条件可消尺度。
- [证明](1101/proof.md)、[代码](1101/check.py)、[结果](1101/results.json)、[验收](1101/acceptance.json)、[独立审阅](1101/review_evidence.json)、[核验](1101/verify.py)与[下一步](1101/NEXT.md)。
'''
    roots = {f'research_cognition_physics/{n}' for n in
             ['README.md', 'research_direction.md', 'RESEARCH_STATE.md', 'ROADMAP.md']}
    changes = []
    for entry in layer['entries']:
        name = entry['destination']; p = ROOT/name; old = p.read_bytes()
        base = relocation_original_bytes(old, entry)
        enc = 'utf-8-sig' if old.startswith(b'\xef\xbb\xbf') else 'utf8'
        content = old.decode(enc); nl = '\r\n' if '\r\n' in content else '\n'
        new = content
        if name in roots:
            start = content.index('## 当前目标：')
            end = content.index('## 狭义相对论连接', start)
            new = content[:start]+common.replace('\n', nl)+content[end:]
        elif name.endswith('archive_1086_/README.md'):
            start = content.index('## 当前目标：'); end = content.index('## 阶段目标', start)
            new = content[:start]+stage.replace('\n', nl)+content[end:]
        elif name.endswith('archive_1086_/文件索引.md'):
            pos = content.index('\n')+1
            new = content[:pos]+addition.replace('\n', nl)+content[pos:]
        if new == content:
            continue
        a = base.decode(enc).splitlines(keepends=True); b = new.splitlines(keepends=True)
        offsets = [0]
        for line in a:
            offsets.append(offsets[-1]+len(line))
        edits = []
        for tag, i, j, k, l in SequenceMatcher(None, a, b, autojunk=False).get_opcodes():
            if tag != 'equal':
                edits.append(dict(start=offsets[i], end=offsets[j], before=''.join(a[i:j]),
                    after=''.join(b[k:l]), kind='round1101_full_lorentz_subgroup_and_scale'))
        raw = new.encode(enc)
        entry.update(current_sha256=sha(raw), current_bytes=len(raw), edits=edits)
        assert relocation_original_bytes(raw, entry) == base
        changes.append((p, old, raw))
    assert len(changes) == 6
    assets = dict(hashes)
    for name in ['verify.py', 'publish.py', 'review_evidence.json']:
        assets['1101/'+name] = sha((HERE/name).read_bytes())
    prior = {}
    for name in ['1100/acceptance.json', '1100/proof.md', '1099/check.py',
                 '1098/proof.md', '1097/proof.md', '1086/proof.md',
                 '_admission/six_protocols_after_1100.md']:
        p = STAGE/name
        prior[p.relative_to(ROOT).as_posix()] = sha(p.read_bytes())
    previous = json.loads((STAGE/'1100/acceptance.json').read_text(encoding='utf8'))
    for name, digest in previous['assets_sha256'].items():
        assert sha((STAGE/name).read_bytes()) == digest
    for name, digest in previous['prior_sources_sha256'].items():
        assert sha((ROOT/name).read_bytes()) == digest
    evidence = dict(round=1101, date='2026-10-10',
        status='PASS_CONDITIONAL_LORENTZ_SUBGROUP_AND_SCALE',
        conditional_full_Lorentz_subgroup=True, conditional_scale_elimination=True,
        finite_words_not_topological_closure=True, central_scale_kernel_can_remain=True,
        new_S_adopted=False, S_derived_from_cognition=False, R_O_FD_derived_from_cognition=False,
        actual_clock_equals_cone_interval=False, all_matter_Lorentz_certified=False,
        all_physical_influences_certified=False, entire_conjecture_decided=False,
        new_adopted_axioms=0, scientific_count_increment=0, scientific_count_total=3860,
        assets_sha256=assets, prior_sources_sha256=prior, independent_read_only_reviews=2,
        goal_status='active')
    assert lp.read_bytes() == old_layer
    for p, old, _ in changes:
        assert p.read_bytes() == old, str(p)
    (HERE/'acceptance.json').write_text(json.dumps(evidence, ensure_ascii=False, indent=2)+'\n', encoding='utf8')
    for p, _, raw in changes:
        p.write_bytes(raw)
    lp.write_text(json.dumps(layer, ensure_ascii=False, indent=2)+'\n', encoding='utf8')
    verification = run(HERE/'verify.py')
    layout = run(ROOT/'scripts/run_research_active.py', '--verify-layout')
    receipt = dict(round=1101, date='2026-10-10', status='PASS', goal_status='active',
        goal_content_changed=False, entire_conjecture_decided=False,
        round_verification=verification, layout_verification=layout,
        previous_round_assets_and_prior_sources_unchanged=True,
        navigation_documents_updated=6, new_scientific_count=0,
        stage_paper_sha256=sha((BASE/'认知操作五公理与三维关系空间_阶段论文.md').read_bytes()))
    (HERE/'publication_checks.json').write_text(json.dumps(receipt, ensure_ascii=False, indent=2)+'\n', encoding='utf8')
    mp = BASE/'_migration/active_documents/checks.json'
    maintenance = json.loads(mp.read_text(encoding='utf8'))
    maintenance['latest_maintenance'] = dict(date='2026-10-10', scope='1101 Lorentz finite generation and scale separation',
        active_documents=7, layer_sha256=sha(lp.read_bytes()), layout_verified=True, goal_status='active',
        latest_formal_round=1101, new_formal_rounds=1, receipt='../../archive_1086_/1101/publication_checks.json')
    mp.write_text(json.dumps(maintenance, ensure_ascii=False, indent=2)+'\n', encoding='utf8')
    print(json.dumps(dict(round=1101, status='PASS', independent_reviews=2,
        recomputation_groups=4, navigation_documents_updated=6, layout_passed=True,
        entire_conjecture_decided=False, goal_status='active'), ensure_ascii=False))


if __name__ == '__main__':
    main()
