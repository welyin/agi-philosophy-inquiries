"""Write current indexes after the evidence-preserving early-round migration."""
import json
from pathlib import Path
import organize_research_231_775 as m
from organize_research_001_230 import atomic

ROOT,BASE=m.ROOT,m.RESEARCH
FIRST,SECOND=BASE/'archive_001_222',BASE/'archive_223_230'
LAST=BASE/'archive_764_'
MAINT=LAST/'_migration/layout_001_230_20261004'
manifest=json.loads((MAINT/'manifest.json').read_text('utf8'))
entries=manifest['entries']
mapping={(ROOT/e['original']).resolve():(ROOT/e['destination']).resolve() for e in entries}
for d in (FIRST,SECOND):
    mapping[d/'README.md']=d/'README.md'
records=[]


def write(path,text):
    before=path.read_bytes() if path.exists() else None
    raw=text.encode('utf8')
    if before==raw:
        return
    if before is not None:
        backup=MAINT/'navigation_before_final_edit'/path.relative_to(ROOT)
        assert not backup.exists()
        backup.parent.mkdir(parents=True,exist_ok=True)
        backup.write_bytes(before)
    atomic(path,raw)
    records.append(dict(path=path.relative_to(ROOT).as_posix(),before_sha256=m.sha(before) if before else None,after_sha256=m.sha(raw)))


original_index=json.loads((FIRST/'222/archive_closure/ROUND_INDEX.json').read_text('utf8'))
rows=[]
for row in original_index['rounds']:
    updated=dict(row)
    updated['note']=m.relative(mapping[(FIRST/row['note']).resolve()],FIRST)
    updated['indexed_artifacts']=[m.relative(mapping[(FIRST/a).resolve()],FIRST) for a in row['indexed_artifacts']]
    updated['current_note_sha256']=m.digest(FIRST/updated['note'])
    rows.append(updated)
m.dump(MAINT/'round_index_001_222.json',dict(count=222,rounds=rows))

text='# 第1—222轮：逐轮结论与范围\n\n[阶段README](README.md) · [配套文件索引](文件索引.md)。标题与结论继承原索引，数学结论未因目录整理而改变。\n\n|轮次|报告|原索引结论|适用范围|配套材料|\n|---|---|---|---|---|\n'
for row in rows:
    clean=lambda s:s.replace('|','\\|').replace('\n',' ')
    text+=f"|{row['round']}|[{clean(row['title'])}]({row['note']})|{clean(row['indexed_conclusion'])}|{clean(row['indexed_scope'])}|[{row['round']}/]({row['round']}/)|\n"
write(FIRST/'ROUND_INDEX.md',text)

write(FIRST/'README.md', '''# 第1—222轮研究档案（已结项）

[研究总目录](../README.md) · [阶段论文](../可组合认知结构与复量子状态空间_阶段论文.md) · [逐轮结论与范围](ROUND_INDEX.md) · [配套文件索引](文件索引.md)

## 阶段范围

本阶段建立并审查从认知操作要求到复量子状态空间的条件性重建链；阶段论文保留明确的工作定义、数学前提和适用范围。不把工作定义的采用等同于证明所有认知系统或所有物理现象都满足这些前提。

正式报告为本目录的 `research_note_01.md` 至 `research_note_222.md`，保留原有编号写法。代码、结果、检查及补充材料按轮次放在 `1/` 至 `222/` 中。共用文件按首次可确认的轮次归属，后续轮次通过链接引用。

- [第219轮修订稿](research_note_219.md)及[原稿和清单](219/history/round219_v1/archive_manifest.json)均保留。
- [第222轮报告](research_note_222.md)接至[223—230轮有限维量子理论阶段](../archive_223_230/README.md)。
- [原研究导航](_history/indexes/README.md)和[原归档导航](222/archive_closure/README_before_round_layout.md)按历史语境阅读；当前方向以[研究状态](../RESEARCH_STATE.md)为准。

## 目录约定

|位置|内容|
|---|---|
|`research_note_编号.md`|正式编号研究报告|
|`1/` … `222/`|对应轮次的代码、结果和核验；联合检查放在文件名中最早一轮|
|`219/history/`|219轮修订前的冻结版本|
|`222/archive_closure/`|阶段论文维护、原索引、原清单、历史归档核验和迁移备份|
|`_shared/`|没有单一轮次归属的背景讨论及共用依赖说明|
|`_history/`|原导航和保留的运行缓存；缓存不计科研证据|

原 `research_process/` 的材料已实际迁移。独立的物理构造、信息几何路线仍在项目根目录，不恢复此前已删除的两套归档副本；旧迁移备份仅作为已有历史原件保存。

## 复算与原件保护

代码、JSON结果及历史收据保持原字节；阅读版Markdown仅调整链接目标。原件与路径变化见[本次迁移说明](../archive_764_/_migration/layout_001_230_20261004/README.md)。旧脚本使用固定目录与跨轮导入，统一通过隔离复算入口恢复原布局：

```powershell
python -B -X utf8 research_cognition_physics/archive_764_/_migration/layout_001_230_20261004/replay_early.py --stage 1
```

命令从项目根目录运行，使用既有Python/NumPy环境；可加 `--module local_calibration_contract` 只检查指定模块。[本次复算记录](../archive_764_/_migration/layout_001_230_20261004/replay_checks.json)与[原1915项回归记录](222/archive_closure/stage1_test_results.json)分别保留，不把目录整理计为新研究轮次。
''')

prior=(SECOND/'230/archive_closure/README_before_round_layout.md').read_bytes()
rewritten,_=m.rewrite_links(prior,SECOND/'README.md',SECOND/'README.md',mapping)
old=rewritten.decode('utf-8-sig')
table=old[old.index('## 逐轮证据'):old.index('## 完整性与复算')]
write(SECOND/'README.md','''# 第223—230轮：有限维量子理论（已结项）

[研究总目录](../README.md) · [整合阶段论文v1.1](../可组合认知结构与有限维量子理论_阶段论文.md) · [配套文件索引](文件索引.md)

正式 `research_note_223.md` 至 `research_note_230.md` 直接放在本目录；其它材料分别放在 `223/` 至 `230/`，阶段结项与论文维护材料在相应轮次的 `archive_closure/`。当前研究转至 [archive_764_](../archive_764_/README.md)。

结项范围为有限维量子操作规则与连续可逆动力学形式的条件性重建。229—230轮已整合进论文；保留原合同下的有限精度结论，以及增加连续种子、路径或闭群条件后的精确实现结论。具体自然Hamiltonian、物理钟尺、无限维和引力属于后续问题。

'''+table+'''## 原件、论文版本与复算

[原科学清单](228/archive_closure/STAGE2_MANIFEST.json)及[补充清单](230/archive_closure/STAGE2_CLOSURE_ADDENDUM.json)的哈希与路径含义保持原样；[旧论文v1.0](230/archive_closure/paper_versions/finite_quantum_v1.0.md)、[原补充稿](230/archive_closure/STAGE2_ADDENDUM.md)和[旧归档导航](230/archive_closure/README_before_round_layout.md)继续保留。

代码和结果保持原字节，阅读版报告只更新链接。原目录中的验证脚本依赖历史布局，现通过[隔离复算入口](../archive_764_/_migration/layout_001_230_20261004/README.md)运行八轮原单元测试，不改写结果文件：

```powershell
python -B -X utf8 research_cognition_physics/archive_764_/_migration/layout_001_230_20261004/replay_early.py --stage 2
```

从项目根目录运行。[本次测试记录](../archive_764_/_migration/layout_001_230_20261004/replay_checks.json)和[历史论文整合核验](230/archive_closure/integration_checks.json)分开保存；历史论文的整合生成与哈希合同依照冻结快照解释，不能直接套用到更新过链接的阅读版。
''')

for archive in (FIRST,SECOND):
    selected=[e for e in entries if (ROOT/e['destination']).is_relative_to(archive)]
    text='# 配套文件索引\n\n[阶段README](README.md)。正式报告在阶段根目录，其余材料按轮次存放；原路径和哈希见[迁移记录](../archive_764_/_migration/layout_001_230_20261004/manifest.json)。缓存保留但不逐项列入研究索引。\n\n|轮次|文件|原路径|\n|---|---|---|\n'
    for e in selected:
        if '/runtime_cache/' in e['destination']:
            continue
        p=ROOT/e['destination']
        text+=f"|{e['round_owner'] or '公共/历史'}|[{p.name}]({m.relative(p,archive)})|`{Path(e['original']).relative_to(archive.relative_to(ROOT)).as_posix()}`|\n"
    write(archive/'文件索引.md',text)

text=(LAST/'README.md').read_text('utf8')
text=text.replace('# 764—775轮：共同量子场与局部来源','# 第764轮起：共同量子场与局部来源（持续研究目录）',1)
text=text.replace('正式报告直接位于本目录；配套代码、结果和核验在相应轮次目录中。','目前已完成764—775轮，下一轮为776。正式报告直接位于本目录；配套代码、结果和核验在相应轮次目录中。按用户要求，本目录由 `archive_764_775` 改名为 `archive_764_`，后续研究持续写入这里。',1)
text=text.replace('恢复研究时从已保存的776入口开始；先核高次共同规范化，不重复已完成首项或旧空间证明。','从[776已有入口](776/drafts/STATUS.md)接续，先核高次共同规范化，不重复已完成首项或旧空间证明。完成后新增本目录下的 `research_note_776.md`，代码、结果与草稿放入 `776/`。本次整理未新增研究结果，也未更改应用目标状态。')
write(LAST/'README.md',text)
text=(BASE/'README.md').read_text('utf8')
text=text.replace('目标由用户暂停。231—775轮共545份编号报告，已实际分拆到以下15个并列目录。','001—775轮共775份正式编号报告，现统一采用“阶段根目录放报告、编号目录放配套材料”的格式。231轮以后的15个阶段保留既有划分，末阶段 `archive_764_` 为持续研究目录。')
text=text.replace('|15|764—775|','|15|764起（当前至775）|')
text=text.replace('恢复研究从[776已有入口](archive_764_/776/drafts/STATUS.md)开始；本次整理不会自动恢复目标。','后续在 [archive_764_](archive_764_/README.md) 从[776已有入口](archive_764_/776/drafts/STATUS.md)接续。应用目标仍保留用户暂停状态，本次整理未自动恢复目标或新增科学轮次。')
write(BASE/'README.md',text)
text=(BASE/'research_direction.md').read_text('utf8')
text=text.replace('本次只做231—775分期归档','本次完成001—230轮格式统一，并将末阶段更名为archive_764_')
text=text.replace('新增研究另建工作目录，不覆盖本次归档。','后续在 `archive_764_` 持续研究：正式报告直接放阶段根目录，配套材料放编号目录，下一编号776；已有775轮及以前的证据继续保留。')
write(BASE/'research_direction.md',text)
text=(BASE/'RESEARCH_STATE.md').read_text('utf8')
text=text.replace('## 2026-10-04：目标暂停，231—775分期归档','## 2026-10-04：001—775目录格式统一，archive_764_持续研究')
text=text.replace('- 231—775共545份正式报告，实际分拆至15个并列archive目录。','- 001—230新增完成按轮次整理；全部001—775共775份正式报告无缺号。231轮以后的15个阶段划分保持，末阶段改名为 `archive_764_`。')
text=text.replace('新增科学工作应另建目录，引用新路径，不恢复向旧archive_231_写入。','后续科学工作写入 `archive_764_`；776的草稿与材料放 `776/`，正式报告放阶段根目录。')
text+='\n## 本次整理核验\n\n[001—230迁移及复算](archive_764_/_migration/layout_001_230_20261004/README.md)；[目录与链接核验](archive_764_/_migration/layout_001_230_20261004/layout_checks.json)。旧代码、结果、历史收据和历史版本继续保留；未新增科学轮次。\n'
write(BASE/'RESEARCH_STATE.md',text)
text=(ROOT/'README.md').read_text('utf8')
text=text.replace('**231—775轮已按15个并列archive目录归档，研究目标暂停。** [阶段目录与成果](research_cognition_physics/README.md)。报告、代码和结果已实际迁移；恢复研究从776入口继续。下方保留历史进展。','**001—775轮已统一为报告在阶段根目录、材料按轮次归档。** [阶段目录与成果](research_cognition_physics/README.md)。后续在 [archive_764_](research_cognition_physics/archive_764_/README.md) 从776入口继续；应用目标仍保留暂停状态。下方保留历史进展。',1)
write(ROOT/'README.md',text)
text=(LAST/'_migration/README.md').read_text('utf8').replace('archive_764_775的_shared与_history','archive_764_的_shared与_history').replace('旧脚本直接在各code目录运行','旧脚本直接在各编号目录运行')
text+='\n## 001—230轮与持续目录\n\n001—230已使用相同轮次布局。末阶段改名为 `archive_764_`，776入口在 `776/drafts/`。[本次整理与完整早期复算](layout_001_230_20261004/README.md)提供新路径清单及原件快照；当前231—775文件清单已更新目的路径，原始大快照未改动。\n'
write(LAST/'_migration/README.md',text)
write(MAINT/'README.md','''# 001—230轮目录整理与末阶段更名

2026-10-04。本次是文件组织维护，不增加研究轮次或科学结论。

## 变更

- 001—222和223—230两个阶段保留原边界，230篇正式报告放各阶段根目录，配套材料归入对应数字目录。
- `archive_764_775` 更名为 `archive_764_`，后续仍在此研究。下一轮入口为 [776/drafts/STATUS.md](../../776/drafts/STATUS.md)。
- 共用背景与原导航单列；219修订前版本、论文旧版、原核验及迁移备份保留。

[1045份早期文件的迁移清单](manifest.json) · [迁移前快照](before_layout.zip) · [快照哈希](snapshot.json) · [全部775轮及链接检查](layout_checks.json) · [早期全量单元测试](replay_checks.json)

## 字节与路径合同

代码、结果和历史收据保持原字节。阅读版Markdown只更新链接目标；精确替换位置记入清单。原导航的完整改写另存前版本和变更收据。冻结稿、原哈希清单中的旧路径按写作时布局理解。

231—775的[当前迁移清单](../manifest.json)已同步目的路径，原[大快照](../original_workspace_231_775.zip)没有改动。本次不重新验证数学论证，也不把两条已有整链审计限制改成通过，详见[此前复算范围](../README.md)。

## 复算

从项目根目录用既有Python/NumPy运行：

```powershell
python -B -X utf8 research_cognition_physics/archive_764_/_migration/layout_001_230_20261004/replay_early.py
python -B -X utf8 research_cognition_physics/archive_764_/_migration/layout_001_230_20261004/replay_early.py --stage 1 --module local_calibration_contract
python -B -X utf8 research_cognition_physics/archive_764_/_migration/replay.py --script verify_round775.py
```

早期入口检查现存代码与冻结原件一致，在临时目录恢复原有布局，只运行原单元测试，不调用结果写入器。默认检查第一阶段原218个测试模块的1915项测试与第二阶段八轮70项测试；SciPy复用项目已有 `.research_runtime`。可用 `--stage 1` 或 `--stage 2` 限定范围，`--report` 只允许创建新的收据。

历史发布核验要求当时的完整导航和论文哈希；本次单元测试与当前路径完整性检查分别记录，不把历史发布合同默默套到新阅读版。
''')
m.dump(MAINT/'navigation_updates.json',dict(documents=records,new_scientific_rounds=0,goal_status_unchanged=True))
print('Updated navigation files:',len(records))
