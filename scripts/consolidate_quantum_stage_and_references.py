"""Consolidate reference locations and present rounds 001–230 as one major stage."""
import json
import os
from pathlib import Path
import re
import subprocess
import zipfile
import organize_research_231_775 as m
from organize_research_001_230 import compose_link_edits

ROOT,BASE=m.ROOT,m.RESEARCH
MIG=BASE/'archive_764_/_migration'
EARLY=MIG/'layout_001_230_20261004'
RECEIPT=MIG/'quantum_stage_001_230_20261004'
OLD_OVERVIEW=BASE/'archive_217_222/001—222阶段总览.md'
OVERVIEW=BASE/'archive_223_230/001—230阶段总览.md'
NAMES=['认知联合体架构-AGI工程方案.md','研究方向讨论：从认知架构到物理世界.md']
BASELINE={}


def baseline(path):
    return BASELINE[path] if path in BASELINE else path.read_bytes()


def atomic(path,raw):
    assert path.resolve().is_relative_to(ROOT.resolve())
    temporary=path.with_name(path.name+'.stage_update_tmp')
    assert not temporary.exists()
    temporary.write_bytes(raw)
    os.replace(temporary,path)


def read(path):
    return json.loads(baseline(path).decode('utf8'))


def key(path):
    return os.path.normcase(os.path.abspath(path))


def edit_targets(raw,source,mapping):
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    text=raw.decode(enc);edits=[]
    for a,b,target,local in m.links(text):
        destination=mapping.get(key(source.parent/local))
        if destination is None:
            continue
        fragment=('#'+target.split('#',1)[1]) if '#' in target else ''
        after=m.relative(destination,source.parent)+fragment
        if target!=after:
            edits.append(dict(start=a,end=b,before=target,after=after))
    for e in reversed(edits):
        text=text[:e['start']]+e['after']+text[e['end']:]
    return text.encode(enc),edits


def main():
    resuming=(RECEIPT/'before_update.zip').exists()
    assert OLD_OVERVIEW.is_file() and not OVERVIEW.exists()
    if resuming:
        assert not (RECEIPT/'updated_documents.json').exists()
        with zipfile.ZipFile(RECEIPT/'before_update.zip') as z:
            BASELINE.update({ROOT/name:z.read(name) for name in z.namelist()})
    early,late=read(EARLY/'manifest.json'),read(MIG/'manifest.json')
    phase_data=read(ROOT/'scripts/research_phases_001_222.json')
    phases=phase_data['phases']
    copies=[BASE/'archive_217_222/_shared'/name for name in NAMES]
    originals=[ROOT/'猜想'/name for name in NAMES]
    comparison=[]
    for copy,original in zip(copies,originals):
        assert copy.resolve().parent==(BASE/'archive_217_222/_shared').resolve()
        assert original.is_file() and not copy.is_symlink() and not original.is_symlink()
        raw,current=copy.read_bytes(),original.read_bytes()
        edits=compose_link_edits(raw,current)
        comparison.append(dict(deleted_copy=copy.relative_to(ROOT).as_posix(),
            canonical=original.relative_to(ROOT).as_posix(),deleted_sha256=m.sha(raw),canonical_sha256=m.sha(current),
            substantive_text_identical=True,link_only_differences=edits))
    mapping={key(a):b for a,b in zip(copies,originals)}
    mapping[key(OLD_OVERVIEW)]=OVERVIEW
    raw_frozen={key(ROOT/e['destination']) for e in early['entries'] if not e['rewrite_markdown_links']}
    raw_frozen|={key(BASE/e['destination']) for e in late['entries'] if not e['rewrite_markdown_links']}
    split=read(EARLY/'phase_split_plan.json')
    raw_frozen|={key(ROOT/p) for p in split['raw_historical_navigation']}
    files=subprocess.check_output(['rg','--files','--hidden','-g','*.md','-g','!.git/**','-g','!.codex/**','-g','!.agents/**','-g','!.research_runtime/**'],cwd=ROOT,text=True,encoding='utf8').splitlines()
    before,outputs,link_changes={},{},[]
    for relative in files:
        p=ROOT/relative
        if key(p) in raw_frozen or any(t in ('navigation_before_final_edit','navigation_before_phase_summary') for t in p.parts):
            continue
        if p in copies or p in originals:
            continue
        raw=baseline(p)
        if not any(name.encode('utf8') in raw for name in NAMES+[OLD_OVERVIEW.name]):
            continue
        changed,edits=edit_targets(raw,p,mapping)
        if raw!=changed:
            before[p]=raw;outputs[p]=changed
            link_changes.append(dict(path=p.relative_to(ROOT).as_posix(),edits=edits))

    def modify(path,fn):
        if path not in before:
            before[path]=baseline(path)
        raw=outputs.get(path,before[path])
        enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
        outputs[path]=fn(raw.decode(enc)).encode(enc)

    def prefix_only(fn):
        def change(text):
            marker='\n## 原项目根README迁入的历史记录'
            head,sep,tail=text.partition(marker)
            return fn(head)+sep+tail
        return change

    overview_link='archive_223_230/001—230阶段总览.md'
    def main_readme(text):
        text=text.replace('## 001—222阶段目录','## 第一大阶段：量子理论重建（001—230）')
        text=text.replace('|轮次|主题|','001—230视为同一个大阶段，包含状态空间研究（001—222）和操作与动力学研究（223—230），按问题转折保留以下12个子阶段目录。\n\n|轮次|主题|',1)
        anchor='|217—222|[可融入性、公理澄清与局部层析](archive_217_222/README.md)|'
        assert anchor in text
        text=text.replace(anchor,anchor+'\n|223—230|[有限维量子操作、连续动力学与结项](archive_223_230/README.md)|',1)
        text=text.replace(f'[前期阶段总览]({overview_link})。223—230仍为[有限维量子理论阶段](archive_223_230/README.md)。',f'[001—230大阶段总览]({overview_link})。两篇综合论文分别覆盖状态空间部分与操作、动力学部分，合起来构成本大阶段的成果。')
        text=text.replace('## 231—775阶段目录','## 后续物理生成与统一研究（231起）')
        text=text.replace('## 前两阶段论文','## 量子理论重建大阶段的两篇综合论文')
        text=text.replace('： [001—222分期总览]', '：[001—230大阶段总览]')
        text=text.replace('[001—222分期总览]', '[001—230大阶段总览]')
        return text
    modify(BASE/'README.md',prefix_only(main_readme))
    modify(BASE/'research_direction.md',lambda t:t.replace('本次完成001—230轮格式统一，将001—222按问题转折分成11段，并将末阶段更名为archive_764_，不改写应用中的目标。',
        '001—230统一归为“量子理论重建”大阶段，含12个子阶段；其中001—222研究状态空间，223—230补齐操作与动力学。231轮起进入后续物理生成与统一研究，当前目录为archive_764_。本次只更新组织与引用，不改写应用中的目标。'))
    modify(BASE/'RESEARCH_STATE.md',lambda t:t.replace('001—222进一步按问题转折分为11个阶段，223—230保留原阶段；',
        '001—230作为同一量子理论重建大阶段，按问题转折分为12个子阶段；').replace('001—222的[分期依据与阅读入口]', '001—230的[大阶段总览与阅读入口]'))

    for p in phases:
        d=BASE/f"archive_{p['start']:03d}_{p['end']:03d}"
        def change(text):
            text=text.replace('[001—222总览]','[001—230大阶段总览]')
            sentence='正式编号报告直接放本目录；代码、结果、检查和补充稿放相应轮次的数字目录。'
            text=text.replace(sentence,'本目录属于001—230轮“量子理论重建”大阶段，是其中一个子阶段。\n\n'+sentence,1)
            text=text.replace('## 阶段成果','## 本段成果').replace('前一阶段：','前一子阶段：').replace('[后一阶段]','[后一子阶段]')
            text=text.replace('本目录同时保存001—222阶段的公共背景', '本目录同时保存001—222轮的公共背景')
            if d.name=='archive_217_222':
                text=text.replace('223—230继续补齐有限维量子过程、连续可逆动力学与有限实验逼近。','同一大阶段的最后一段223—230继续补齐有限维量子过程、连续可逆动力学与有限实验逼近。')
                text+='\n两篇背景文档统一引用“猜想”目录：[认知联合体架构](../../猜想/认知联合体架构-AGI工程方案.md)、[从认知架构到物理世界的方向讨论](../../猜想/研究方向讨论：从认知架构到物理世界.md)。\n'
            return text
        modify(d/'README.md',prefix_only(change))
    def finish_part(text):
        text=text.replace('# 第223—230轮：有限维量子理论（已结项）','# 第223—230轮：量子理论重建大阶段的操作、动力学与收尾')
        text=text.replace('## 逐轮证据','## 在001—230大阶段中的位置\n\n001—222建立状态空间与相关操作原则，223—230是同一大阶段的最后一个子阶段，补齐量子过程、连续可逆动力学和有限实验范围；大阶段在230轮收尾。参见[001—230大阶段总览](001—230阶段总览.md)。原“第一阶段／第二阶段”是写作时的分段称呼，历史收据及版本沿用原名。\n\n## 逐轮证据',1)
        return text
    modify(BASE/'archive_223_230/README.md',prefix_only(finish_part))

    def overview(text):
        text=text.replace('# 001—222轮阶段总览','# 001—230轮：量子理论重建大阶段总览')
        text=text.replace('按问题转折划分11段；这是阅读组织，不改变原结论、轮次与历史适用范围。','001—230作为同一个大阶段，按问题转折划分12个子阶段；其中001—222研究状态空间，223—230补齐操作与动力学。此处更新组织层次，不改变数学前提、轮次与结论范围。')
        text=text.replace('|轮次|阶段|转入下一阶段的原因|','|轮次|子阶段|衔接|')
        anchor='|217—222|[可融入性、公理澄清与局部层析](../archive_217_222/README.md)|223—230继续补齐有限维量子过程、连续可逆动力学与有限实验逼近。|'
        assert anchor in text
        text=text.replace(anchor,anchor+'\n|223—230|[有限维量子操作、连续动力学与结项](../archive_223_230/README.md)|在有限维及明示正则条件下完成量子理论结构重建；231轮起进入物理生成与统一研究。|')
        text=text.replace('论文：[可组合认知结构与复量子状态空间](../可组合认知结构与复量子状态空间_阶段论文.md)。后续有限维量子过程见[223—230](../archive_223_230/README.md)，当前继续研究见[archive_764_](../archive_764_/README.md)。',
            '## 两篇论文，共同组成一大阶段\n\n- [状态空间部分：可组合认知结构与复量子状态空间](../可组合认知结构与复量子状态空间_阶段论文.md)，综合001—222。\n- [操作与动力学部分：可组合认知结构与有限维量子理论](../可组合认知结构与有限维量子理论_阶段论文.md)，综合223—230。\n\n结项范围是有限维量子理论结构的条件性重建：原合同下的仪器逼近、增加连续种子或相应正则条件后的精确操作，以及Time条件下的连续可逆动力学形式。具体自然哈密顿量、物理钟尺、无限维场论、时空与引力仍属于后续研究。\n\n后续从[231—258](../archive_231_258/README.md)接续，当前工作目录为[archive_764_](../archive_764_/README.md)。')
        return text
    modify(OLD_OVERVIEW,overview)

    paper1=BASE/'可组合认知结构与复量子状态空间_阶段论文.md'
    paper2=BASE/'可组合认知结构与有限维量子理论_阶段论文.md'
    modify(paper1,lambda t:t.replace('**阶段论文 v1.0 ·','**量子理论重建大阶段·状态空间部分论文 v1.0 ·',1).replace('## 摘要','编排说明（2026-10-04）：001—230统一归为量子理论重建大阶段；本文覆盖001—222，[操作与动力学部分](可组合认知结构与有限维量子理论_阶段论文.md)覆盖223—230。两篇论文共同构成本大阶段的成果，原证明范围不变。\n\n## 摘要',1))
    modify(paper2,lambda t:t.replace('**第二阶段论文 v1.1','**量子理论重建大阶段·操作与动力学部分论文 v1.1',1).replace('第一阶段论文','状态空间部分论文').replace('第一阶段','状态空间部分').replace('**第二阶段于第230轮完成收尾。**','**001—230轮量子理论重建大阶段于第230轮完成收尾。**').replace('## 摘要','编排说明（2026-10-04）：本文与001—222轮的状态空间研究属于同一大阶段，并构成最后一个子阶段。历史核验中的stage1／stage2仍指原测试与出版分段。\n\n## 摘要',1))
    def root_readme(text):
        text=text.replace('### 一、认知物理主线（research_cognition_physics/，230轮，第二阶段收尾完成）','### 一、认知物理主线（research_cognition_physics/）')
        text=text.replace('**阶段一成果（F+U+C+P ⟹ 复矩阵状态锥）**：','001—230轮共同构成“量子理论重建”大阶段：001—222为状态空间部分，223—230为操作与动力学部分。后续进展见[研究主目录](research_cognition_physics/README.md)。\n\n**状态空间部分成果（001—222轮，F+U+C+P ⟹ 复矩阵状态锥）**：')
        text=text.replace('**阶段二成果（223—228轮）**','**操作与动力学部分成果（223—230轮）**').replace('六轮55项检查','八轮70项检查')
        text=text.replace('# 认知物理主线（228轮，两阶段已结项）','# 认知物理主线（001—230为量子理论重建大阶段）')
        text=text.replace('# 第一阶段论文','# 状态空间部分论文').replace('# 第二阶段论文','# 操作与动力学部分论文')
        text=text.replace('│   ├── archive_001_222/                 # 222 轮完整档案\n│   │   ├── research_process/             # 每轮原始笔记与代码\n│   │   └── README.md                     # 档案入口与复算说明\n│   └── archive_223_/                    # 223—228轮、结项清单与核验入口',
            '│   ├── archive_001_020/ 等             # 001—222的11个子阶段\n│   └── archive_223_230/                 # 同一大阶段的最后子阶段及结项材料')
        return text
    modify(ROOT/'README.md',root_readme)
    def file_index(text):
        lines=[line for line in text.splitlines() if not (line.startswith('|') and any(name in line for name in NAMES))]
        return '\n'.join(lines)+'\n\n## 外置参考原件\n\n两份重复副本已按用户要求移除，统一使用[认知联合体架构](../../猜想/认知联合体架构-AGI工程方案.md)与[研究方向讨论](../../猜想/研究方向讨论：从认知架构到物理世界.md)。原版本及删除前哈希保留在迁移快照中。\n'
    modify(BASE/'archive_217_222/文件索引.md',file_index)
    modify(BASE/'archive_231_258/README.md',prefix_only(lambda t:t.replace('## 研究问题','本段接在001—230轮[量子理论重建大阶段](../archive_223_230/001—230阶段总览.md)之后，开始物理生成与统一研究。\n\n## 研究问题',1)))

    # Confirm the requested editorial update leaves displayed math untouched.
    for paper in (paper1,paper2):
        old_math=[x.group() for x in m.EXCLUDED_BLOCK.finditer(before[paper].decode('utf-8-sig'))]
        new_math=[x.group() for x in m.EXCLUDED_BLOCK.finditer(outputs[paper].decode('utf-8-sig'))]
        assert old_math==new_math,paper
    personal_head=before[ROOT/'README.md'].decode('utf8').split('## 研究进展',1)[0]
    personal_tail=before[ROOT/'README.md'].decode('utf8').split('## 给后来者的话',1)[1]
    assert outputs[ROOT/'README.md'].decode('utf8').split('## 研究进展',1)[0]==personal_head
    assert outputs[ROOT/'README.md'].decode('utf8').split('## 给后来者的话',1)[1]==personal_tail

    # Save every previous version, including the explicitly removed duplicates.
    if not resuming:
        RECEIPT.mkdir()
    snapshot={**before,**{p:baseline(p) for p in copies}}
    metadata=[EARLY/'manifest.json',MIG/'manifest.json',EARLY/'round_index_001_222.json',ROOT/'scripts/research_phases_001_222.json']
    snapshot.update({p:baseline(p) for p in metadata})
    if not resuming:
        with zipfile.ZipFile(RECEIPT/'before_update.zip','x',zipfile.ZIP_DEFLATED) as z:
            for p,raw in sorted(snapshot.items()):
                z.writestr(p.relative_to(ROOT).as_posix(),raw)
    with zipfile.ZipFile(RECEIPT/'before_update.zip') as z:
        for p,raw in snapshot.items():
            assert z.read(p.relative_to(ROOT).as_posix())==raw
    if not resuming:
        m.dump(RECEIPT/'snapshot_checks.json',dict(file='before_update.zip',sha256=m.digest(RECEIPT/'before_update.zip'),entries=[dict(path=p.relative_to(ROOT).as_posix(),sha256=m.sha(raw)) for p,raw in snapshot.items()]))
    for p,raw in outputs.items():
        actual=p.read_bytes()
        assert actual in (before[p],raw),str(p)
        if actual!=raw:
            atomic(p,raw)
    assert OLD_OVERVIEW.resolve().is_relative_to(BASE.resolve()) and OVERVIEW.resolve().parent==(BASE/'archive_223_230').resolve()
    OLD_OVERVIEW.replace(OVERVIEW)
    for p,item in zip(copies,comparison):
        assert p.resolve().parent==(BASE/'archive_217_222/_shared').resolve()
        assert m.digest(p)==item['deleted_sha256'] and m.digest(ROOT/item['canonical'])==item['canonical_sha256']
        p.unlink()

    # The two retired entries resolve to the user-owned canonical references;
    # the saved historical bytes remain in the original immutable snapshots.
    with zipfile.ZipFile(EARLY/'before_layout.zip') as z:
        for e in early['entries']:
            prior=ROOT/e['destination']
            if prior in copies:
                target=originals[copies.index(prior)]
                e['removed_duplicate_path']=e['destination']
                e['destination']=target.relative_to(ROOT).as_posix()
                e['canonical_external_reference']=True
            else:
                target=prior
            if prior not in outputs and prior not in copies:
                continue
            original=z.read(e['original']);current=target.read_bytes()
            e.update(current_sha256=m.sha(current),current_bytes=len(current),link_edits=compose_link_edits(original,current))
    with zipfile.ZipFile(MIG/'original_workspace_231_775.zip') as z:
        for e in late['entries']:
            p=BASE/e['destination']
            if p not in outputs:
                continue
            original=z.read('research_cognition_physics/archive_231_/'+e['original']);current=p.read_bytes()
            e.update(current_sha256=m.sha(current),current_bytes=len(current),link_edits=compose_link_edits(original,current))
    major=dict(title='量子理论重建',rounds=[1,230],status='closed within the stated finite-dimensional scope',
               parts=[dict(rounds=[1,222],name='状态空间'),dict(rounds=[223,230],name='操作与动力学')],
               subphases=[dict(rounds=[p['start'],p['end']],directory=f"archive_{p['start']:03d}_{p['end']:03d}",title=p['title']) for p in phases]+[dict(rounds=[223,230],directory='archive_223_230',title='有限维量子操作、连续动力学与结项')],
               overview=OVERVIEW.relative_to(BASE).as_posix(),next_major_stage_starts_at=231)
    early['version']+=1;late['version']+=1
    early['major_stage']=major
    early['canonical_reference_changes']=comparison
    late['major_stage_001_230']=major
    for p,data in [(EARLY/'manifest.json',early),(MIG/'manifest.json',late)]:
        atomic(p,(json.dumps(data,ensure_ascii=False,indent=2)+'\n').encode('utf8'))
    index=read(EARLY/'round_index_001_222.json')
    for row in index['rounds']:
        row['current_note_sha256']=m.digest(BASE/row['note'])
    atomic(EARLY/'round_index_001_222.json',(json.dumps(index,ensure_ascii=False,indent=2)+'\n').encode('utf8'))
    phase_data['major_stage']=major
    atomic(ROOT/'scripts/research_phases_001_222.json',(json.dumps(phase_data,ensure_ascii=False,indent=2)+'\n').encode('utf8'))
    m.dump(RECEIPT/'stage_structure.json',major)
    m.dump(RECEIPT/'reference_changes.json',dict(deleted_duplicates=comparison,link_changes=link_changes,overview_previous=OLD_OVERVIEW.relative_to(ROOT).as_posix(),overview_current=OVERVIEW.relative_to(ROOT).as_posix()))
    m.write(RECEIPT/'README.md','''# 参考文档去重与001—230大阶段整合

2026-10-04。两份指定背景文档以项目“猜想”目录为原件，归档内重复副本已删除，现行引用已调整。SoCA两版本的差异仅是文档索引链接；方向讨论文本完全一致。原件保持不变，旧副本原文保存在既有快照及[本次更新前快照](before_update.zip)中。

001—230统一归为量子理论重建大阶段，保留12个子阶段目录；两篇论文分别覆盖状态空间和操作、动力学。231轮起进入后续物理生成与统一研究。历史报告、出版清单及旧进展中“第一／第二阶段”的称呼按写作时语境保留，原测试入口stage1／stage2也仅代表测试套件。

[阶段结构](stage_structure.json) · [引用变更](reference_changes.json) · [快照核验](snapshot_checks.json) · [最终核验](checks.json)

本次修改只涉及组织说明与引用；论文公式、科学代码、实验结果、研究目标和当前轮次保持不变。
''')
    m.dump(RECEIPT/'updated_documents.json',dict(documents=[dict(previous=p.relative_to(ROOT).as_posix(),current=(OVERVIEW if p==OLD_OVERVIEW else p).relative_to(ROOT).as_posix(),sha256=m.digest(OVERVIEW if p==OLD_OVERVIEW else p)) for p in outputs],papers_math_unchanged=True,root_personal_sections_unchanged=True,new_scientific_rounds=0))
    print(json.dumps(dict(deleted_duplicates=2,updated_documents=len(outputs),major_stage=[1,230],subphases=12,overview=str(OVERVIEW)),ensure_ascii=False))


if __name__=='__main__':
    main()
