"""Publish round 530 with retained navigation snapshots and concurrent-edit checks."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core
import verify_round530 as evidence

HERE = Path(__file__).resolve().parent
RESEARCH, ROOT = HERE.parent, HERE.parent.parent
checked = evidence.verify()
assert core.read(HERE/'research_round_530_checks.json') == checked
paths = [ROOT/'README.md', RESEARCH/'README.md', RESEARCH/'research_direction.md',
    RESEARCH/'RESEARCH_STATE.md', HERE/'README.md',
    HERE/'spatial_premise_closure_audit.md', HERE/'three_dimensional_four_conditions_review.md']
folder = HERE/'navigation_before_round530_20260930'
resuming = folder.exists()
assert not (HERE/'round530_navigation_checks.json').exists(), 'already published'
if resuming:
    saved_manifest = core.read(folder/'manifest.json')
    raws = {p:(folder/saved_manifest[str(p.relative_to(ROOT))]['snapshot']).read_bytes()
            for p in paths}
    for p, raw in raws.items():
        assert hashlib.sha256(raw).hexdigest() == saved_manifest[str(p.relative_to(ROOT))]['sha256']
else:
    raws = {p:p.read_bytes() for p in paths}
    folder.mkdir(exist_ok=False)
count = checked['cumulative_unique_protected_evidence_files']
summary = (
    '**第530轮完成：** [无永久边界支撑的有限局部参考]({p}research_note_530.md)'
    '在529同一关系门及普通闭合周期图中，以有限共同场准备和缓冲区替代永久仿射twist，'
    '证明有限时窗的均值误差、边缘能量及实际有限记录保证。'
    '所给三维例的一批81个读数用729位，全部核心坐标差同时正确的成功率至少99.1965%。'
    '几何、地址、热源及准备仍输入，未从旧度3动态图生成该度6格点，也未选择三维。'
    f'6组、12式及独立终审通过；最新530／2622，902份编号科学文件，{count}份保护证据；'
    '[核验]({p}research_round_530_checks.json)。本段为最新状态，旧接续保留为历史。')
next_steps = '''

### 第530轮后：检验实际条件记录能否形成参考来源

530在528—529同一G与普通闭合周期H中，给出有限时窗来源替代：有限支撑共同Weyl准备在核心与缓冲区仿射，外围回零；边缘梯度能全部计入。所有更大容纳盒中，核心均值偏差有整函数尾界，同一真实χ读取继续使用529非Gaussian热尾。所给27站、三份参考的81个末读有明确有限编码与至少99.1965%联合成功率。没有重算环路、同步、BL或旧IR阈值来增加轮次。

该结果取消永久twist和边界钉死，保留一次指定非均匀f的准备及其程序／地址。额外位移能量对更大盒一致，背景总热能仍随规模增长。均值界对全时窗成立，概率只保证一次预选时刻的整批读取；末读回冲后不能自动复用产品热态。三维例是给定度6格点上的条件构造，旧度3动态图到它的映射未证；有限局部任务本身也允许其他维数。

**下一项先回查既有条件化来源：** 在同一H及同一实际χ读口下，少量有限标记记录是否可以形成可用的参考剖面，从而减少逐站指定f？必须核真实条件态、来源成功概率及对称退化。χ含U与非Gaussian W，条件测量通常使两部门关联；不能直接把529产品矩生成函数移到后验，不能把随机数值不同当成正确坐标。已有Gaussian条件化、520状态学习、287身份桥和来源障碍应直接继承，未解决的接口才形成下一编号报告。

原动态图到有效几何、全局实际图册、连续方向壳、三维唯一性与量子应力—几何反馈仍开放。当前双向逻辑统一目标保持；每轮及时交付完整research_note_编号.md、代码、结果及核验，保留所有旧稿，不以未编号长报告替代正式编号。空间阶段未结项，不改应用目标、不设置定时任务、不做图像检验。
'''
planned, manifest = {}, {}
for path, raw in raws.items():
    key = ('root_' if path.parent==ROOT else 'research_' if path.parent==RESEARCH else 'archive_')+path.name
    if not resuming:
        with (folder/key).open('xb') as stream:
            stream.write(raw)
    manifest[str(path.relative_to(ROOT))] = dict(snapshot=key, sha256=hashlib.sha256(raw).hexdigest())
    encoding = 'utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline = '\r\n' if b'\r\n' in raw else '\n'
    body = raw.decode(encoding).replace('\r\n', '\n')
    assert '**第530轮完成：**' not in body
    prefix = 'research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block = summary.format(p=prefix)
    if path in paths[:5]:
        title, tail = body.split('\n\n', 1)
        body = title+'\n\n'+block+'\n\n'+tail
    elif path.name == 'spatial_premise_closure_audit.md':
        assert '\n## 175.' not in body
        body += '\n\n## 175. 有限准备替代永久边界的局部参考\n\n'+block+next_steps
    else:
        assert '\n## 80.' not in body
        body += '\n\n## 80. 第530轮：有限局部参考不选择维数\n\n'+block
    if path == HERE/'README.md':
        body += ('\n\n## 第530轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|530|[无永久边界支撑的有限局部参考](research_note_530.md)|'
            '[代码](finite_patch_reference_source.py)、[结果](finite_patch_reference_source_results.json)、'
            '[核验](research_round_530_checks.json)|\n')
    if path.name == 'RESEARCH_STATE.md':
        body += next_steps
    if path.name == 'research_direction.md':
        assert '最新科学轮次与检查数为529／2616' in body
        body = body.replace('最新科学轮次与检查数为529／2616', '最新科学轮次与检查数为530／2622')
        body = body.replace('完成231—529轮。', '完成231—530轮。')
        body += ('\n\n**530后接续（目标不变）：** 无twist的有限局部参考已有同规则、有限时窗和实际记录证书；'
            '转向实际条件记录的来源，核后验关联与身份／几何退化，不重复噪声、环路或IR实验。详见RESEARCH_STATE.md末尾。\n')
    planned[path] = body.replace('\n', newline).encode(encoding)
if not resuming:
    with (folder/'manifest.json').open('x', encoding='utf8') as stream:
        json.dump(manifest, stream, ensure_ascii=False, indent=2)
else:
    assert manifest == saved_manifest
links = 0
for path, body in planned.items():
    for link in core.link_parser()(body.decode('utf-8-sig')):
        assert (path.parent/link).resolve().exists(), (path, link)
        links += 1
for path, raw in raws.items():
    assert path.read_bytes() in (raw, planned[path]), ('concurrent navigation edit', path)
for path, body in planned.items():
    if path.read_bytes() == body:
        continue
    temporary = path.with_name(path.name+'.round530.tmp')
    if temporary.exists():
        assert temporary.read_bytes() == body
    else:
        with temporary.open('xb') as stream:
            stream.write(body)
    assert path.read_bytes() == raws[path], ('concurrent navigation edit', path)
    os.replace(temporary, path)
for group in ('new_file_hashes', 'preserved_draft_hashes'):
    for name, digest in checked[group].items():
        assert core.digest(HERE/name) == digest, name
report = dict(date='2026-09-30', latest_round=530, cumulative_tests=2622,
    numbered_scientific_files=902, unique_protected_evidence_files=count,
    navigation_files=7, navigation_links=links, broken_links=0,
    navigation_snapshots_preserved=True, recovered_partial_navigation_publish=resuming,
    complete_report_published=True, existing_results_preserved=True, next_round=531,
    active_goal_unchanged=True, stage_complete=False, all_checks_passed=True)
with (HERE/'round530_navigation_checks.json').open('x', encoding='utf8', newline='\n') as stream:
    stream.write(json.dumps(report, ensure_ascii=False, indent=2)+'\n')
print(json.dumps(report, ensure_ascii=False))
