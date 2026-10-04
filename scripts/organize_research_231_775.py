"""Explicit, reversible phase migration. Run --plan, inspect, then --apply.

Scientific code/results retain their bytes. Markdown changes are confined to
link destinations. Historical receipts keep their original path/hash meaning;
the replay tool reconstructs that exact layout from the verified snapshot.
"""
import argparse
from collections import Counter, defaultdict
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
from urllib.parse import unquote
import zipfile

ROOT = Path(__file__).resolve().parent.parent
RESEARCH = ROOT/'research_cognition_physics'
OLD = RESEARCH/'archive_231_'
NEW = RESEARCH/'archive_231_775'
DATA = json.loads((ROOT/'scripts/research_phases_231_775.json').read_text('utf8'))
PHASES = DATA['phases']
MIG = NEW/'_migration'
LINK = re.compile(r'(?<!!)\[[^\]\n]*\]\(([^)\n]+)\)|^\[[^\]\n]+\]:\s*(\S+)\s*$', re.M)
EXCLUDED_BLOCK = re.compile(r'\\\[.*?\\\]|\$\$.*?\$\$|^```.*?^```\s*$', re.S | re.M)
# Four unmarked mathematical expressions in preserved pre-link-fix drafts.
# These are multiplication, spectral calculus and coefficient extraction, not links.
NON_LINK_MATH_TOKENS = {'7θ/3', 'H_*', 'A₀', '1−2At'}


def sha(data):
    return hashlib.sha256(data).hexdigest()


def digest(p):
    return sha(p.read_bytes())


def dump(p, obj):
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open('x', encoding='utf8', newline='\n') as f:
        json.dump(obj, f, ensure_ascii=False, indent=2)
        f.write('\n')


def write(p, text):
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open('x', encoding='utf8', newline='\n') as f:
        f.write(text)


def folder(phase):
    return f"{phase['id']}_{phase['start']}_{phase['end']}_{phase['slug']}"


def phase_of(n):
    return next((p for p in PHASES if p['start'] <= n <= p['end']), None)


def links(text):
    masked = EXCLUDED_BLOCK.sub(lambda m: ' '*len(m[0]), text)
    for match in LINK.finditer(masked):
        group = 1 if match[1] is not None else 2
        a, b = match.span(group)
        value = text[a:b]
        stripped = value.strip()
        angled = stripped.startswith('<') and '>' in stripped
        target = stripped[1:stripped.index('>')] if angled else stripped
        # The old collection predominantly has bare targets; preserve optional titles.
        if not angled:
            target = re.split(r'\s+["\']', target, 1)[0]
        local = unquote(target.split('#', 1)[0])
        if local in NON_LINK_MATH_TOKENS:
            continue
        if not local or re.match(r'^[a-zA-Z]+:', local) or local.startswith('//'):
            continue
        target_start = a+value.index(target)
        yield target_start, target_start+len(target), target, local


def relative(target, base):
    return os.path.relpath(target, base).replace('\\', '/')


def make_plan():
    assert OLD.is_dir() and not (MIG/'manifest.json').exists()
    assert [n for p in PHASES for n in range(p['start'],p['end']+1)] == list(range(231,776))
    files = sorted(p for p in OLD.rglob('*') if p.is_file())
    ownership = defaultdict(set)
    for n in range(231,776):
        f = OLD/f'research_note_{n}.md'
        assert f.is_file(), n
        for _, _, _, local in links(f.read_text('utf-8-sig')):
            linked = (f.parent/local).resolve()
            if linked.is_relative_to(OLD) and linked.exists():
                ownership[linked].add(n)
        receipt = OLD/f'research_round_{n}_checks.json'
        if receipt.exists():
            values = json.loads(receipt.read_text('utf-8-sig'))
            for name in values.get('new_file_hashes', {}):
                linked=(OLD/name).resolve()
                if linked.is_relative_to(OLD) and linked.exists():
                    ownership[linked].add(n)
    for f in files:
        if f in ownership or f.suffix!='.py':
            continue
        content=f.read_text('utf-8-sig')
        explicit=re.findall(r'research_round_(\d{3})_checks',content)
        if not explicit:
            explicit=re.findall(r'(?:^|[,{])\s*(\d{3})\s*:',content,re.M)
        candidates=[int(x) for x in explicit if 231<=int(x)<=775]
        if candidates:
            ownership[f].add(max(candidates))
    entries = []
    for f in files:
        rel = f.relative_to(OLD)
        top = rel.parts[0]
        raw = f.read_bytes()
        snapshot = top.startswith('navigation_') or 'navigation_before' in rel.as_posix()
        n = None
        reason = 'shared support; no numbered owner identified'
        candidates = re.findall(r'(?:research_note_|research_round_|round)(\d{3})(?!\d)', top)
        if not candidates:
            candidates = re.findall(r'(?:^|_)([2-7]\d{2})(?=[_.]|$)', top)
        ns = [int(x) for x in candidates if 231 <= int(x) <= 776]
        if ns:
            n, reason = min(ns), 'explicit round in original path'
        elif f in ownership:
            n, reason = min(ownership[f]), 'earliest direct numbered-report reference'
        elif len(rel.parts)>1 and (OLD/top) in ownership:
            n, reason = min(ownership[OLD/top]), 'earliest reference to containing material directory'
        if snapshot:
            dest = Path('_history/navigation_snapshots')/rel
            kind = 'historical_navigation_snapshot'
        elif top=='__pycache__':
            dest = Path('_history/runtime_cache')/rel
            kind = 'historical_runtime_cache'
        elif n==776:
            dest = Path('next_round_776')/Path(*rel.parts[1:]) if len(rel.parts)>1 else Path('next_round_776')/rel
            kind = 'unfinished_next_entry'
        elif len(rel.parts)==1 and f.name in ('README.md','spatial_premise_closure_audit.md','three_dimensional_four_conditions_review.md'):
            dest = Path('_history/indexes')/rel
            kind = 'historical_index'
        else:
            phase = phase_of(n) if n else None
            base = Path(folder(phase)) if phase else Path('_shared')
            if len(rel.parts)>1:
                kind = 'draft_or_support'
                dest = base/'drafts'/rel
            else:
                kind = 'report' if re.fullmatch(r'research_note_\d+\.md', f.name) else 'code' if f.suffix in ('.py','.mjs') else 'result_or_receipt' if f.suffix=='.json' else 'document'
                category = 'code' if kind=='code' else 'results' if kind=='result_or_receipt' else 'notes'
                dest = base/category/f.name
        assert f.resolve().is_relative_to(OLD.resolve())
        assert (NEW/dest).resolve().is_relative_to(NEW.resolve())
        entries.append(dict(original=rel.as_posix(), destination=dest.as_posix(),
                            original_sha256=sha(raw), bytes=len(raw), round_owner=n,
                            attribution=reason, kind=kind,
                            rewrite_markdown_links=f.suffix.lower()=='.md' and not snapshot))
    assert len({e['destination'] for e in entries})==len(entries)
    plan = dict(version=1, source=str(OLD), destination=str(NEW), date='2026-10-04',
                latest_completed_round=775, next_round=776, research_goal_paused=True,
                original_file_count=len(entries), original_bytes=sum(e['bytes'] for e in entries),
                phases=PHASES, entries=entries)
    dump(MIG/'plan.json', plan)
    print(json.dumps(dict(files=len(entries), reports=545,
                         owners=Counter((e['destination'].split('/')[0] for e in entries)),
                         shared=[e['original'] for e in entries if e['destination'].startswith('_shared/')]),ensure_ascii=False,indent=2))


def rewrite_links(raw, src, dst, mapping):
    enc = 'utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    text = raw.decode(enc)
    changes = []
    for a,b,target,local in links(text):
        old_target=(src.parent/local).resolve()
        if old_target in mapping:
            new_target=mapping[old_target]
        elif old_target==OLD:
            new_target=NEW
        elif old_target.is_dir() and old_target.is_relative_to(OLD):
            children=[(p,q) for p,q in mapping.items() if p.is_relative_to(old_target)]
            if not children:
                continue
            # A directory may have split across phases. Its stable index is the root map.
            parents={q.parent for p,q in children if p.parent==old_target}
            new_target=next(iter(parents)) if len(parents)==1 else NEW
        else:
            new_target=old_target
        suffix='#'+target.split('#',1)[1] if '#' in target else ''
        replacement=relative(new_target,dst.parent)+suffix
        if replacement!=target:
            changes.append(dict(start=a,end=b,before=target,after=replacement))
    for c in reversed(changes):
        assert text[c['start']:c['end']]==c['before']
        text=text[:c['start']]+c['after']+text[c['end']:]
    return text.encode(enc), changes


def context_files():
    excluded={'.git','.codex','.agents','.research_runtime','.tmp','node_modules','__pycache__'}
    selected=[]
    for base,dirs,names in os.walk(ROOT):
        base=Path(base)
        dirs[:]=[d for d in dirs if d not in excluded and (base/d)!=NEW]
        for name in names:
            p=base/name
            if not p.is_relative_to(OLD):
                selected.append(p)
    return selected


def snapshot(plan):
    archive=MIG/'original_workspace_231_775.zip'
    assert not archive.exists()
    originals=[OLD/e['original'] for e in plan['entries']]
    context=context_files()
    inventory=[]
    with zipfile.ZipFile(archive,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=3,allowZip64=True) as z:
        for p in originals+context:
            raw=p.read_bytes()
            name=p.relative_to(ROOT).as_posix()
            z.writestr(name,raw)
            inventory.append(dict(path=name,sha256=sha(raw),bytes=len(raw)))
    with zipfile.ZipFile(archive) as z:
        for e in inventory:
            assert sha(z.read(e['path']))==e['sha256'], e['path']
    proof=dict(file=archive.name,sha256=digest(archive),files=len(inventory),
               archive_files=len(originals),context_files=len(context),
               entries=inventory,all_bytes_verified=True)
    dump(MIG/'snapshot_checks.json',proof)
    print(f"Snapshot verified: {len(inventory)} files; compressed bytes {archive.stat().st_size}",flush=True)


def note_link(n, from_dir, mapping):
    return relative(mapping[OLD/f'research_note_{n}.md'],from_dir)


def document_phases(plan,mapping):
    root_rows=[]
    for phase in PHASES:
        base=NEW/folder(phase)
        artifacts=[e for e in plan['entries'] if e['destination'].startswith(folder(phase)+'/')]
        keys='、'.join(f"[{n}]({note_link(n,base,mapping)})" for n in phase['keys'])
        text=f"# 阶段{phase['id']}：{phase['slug']}\n\n轮次：{phase['start']}—{phase['end']}；{phase['end']-phase['start']+1}份正式报告。\n\n"
        text+='[返回总览](../README.md) · [逐文件索引](文件索引.md) · [迁移与复算说明](../_migration/README.md)\n\n'
        text+=f"## 研究问题\n\n{phase['question']}\n\n## 已得到的结果\n\n"+'\n'.join('- '+r for r in phase['results'])+'\n\n'
        text+=f"## 适用范围与未完成项\n\n{phase['limits']}\n\n本目录是研究过程的主题分期，不表示每一项开放问题已结项。经典受限模型、条件性定理、反例及量子构造按原报告的量词保留。\n\n## 关键报告\n\n{keys}\n\n## 与后续阶段的连接\n\n{phase['transition']}\n\n## 全部轮次\n\n|轮次|报告|\n|---|---|\n"
        for n in range(phase['start'],phase['end']+1):
            title=mapping[OLD/f'research_note_{n}.md'].read_text('utf-8-sig').splitlines()[0].lstrip('# ').replace('|','／')
            text+=f"|{n}|[{title}]({note_link(n,base,mapping)})|\n"
        write(base/'README.md',text)
        index='# 文件索引\n\n[阶段综述](README.md)。科学代码与结果保持原字节；跨阶段复算用统一入口。来源归属采用显式轮次或最早报告直接引用，不宣称最早引用就是理论首次提出。\n\n|类型|文件|原路径|\n|---|---|---|\n'
        for e in artifacts:
            index+=f"|{e['kind']}|[{Path(e['destination']).name}]({relative(NEW/e['destination'],base)})|`{e['original']}`|\n"
        write(base/'文件索引.md',index)
        root_rows.append(f"|{phase['id']}|{phase['start']}—{phase['end']}|[{phase['slug']}]({folder(phase)}/README.md)|{phase['question']}|")
    intro='''# 231—775轮：阶段研究档案

整理日期：2026-10-04。545份编号报告全部保留，按主问题转折分为15个阶段。目标由用户暂停；本次做归档，不新增研究轮次，也不宣布整个统一目标完成。

**已有阶段性成果：** 条件性三维与光滑坐标；若干明确前提不足或特定候选失败的反证；同一物质、记录、来源和几何的联合构造；最新完整线性物理态与局部一圈来源的连接。它们尚未合成从少量认知原则无条件推出全部物理的定理。

[阶段成果总览](阶段成果总览.md) · [跨阶段主题索引](跨阶段主题索引.md) · [迁移与复算说明](_migration/README.md) · [待续776入口](next_round_776/STATUS.md)

## 阶段目录

|阶段|轮次|主题|核心问题|
|---|---|---|---|
'''
    intro+='\n'.join(root_rows)+'\n\n## 目录规则\n\n每阶段的`notes`保存原始报告及专题说明，`code`保存代码，`results`保存结果与历史核验，`drafts`保留工作稿。跨阶段公共材料在`_shared`，历史索引与导航快照在`_history`。所有历史文件均有迁移映射，旧核验中的路径与哈希按原快照解释。\n\n本次没有重新证明或重跑全部545轮；迁移核验和代表性复算单独列在迁移报告中，不混入历史科学检查计数。\n'
    write(NEW/'README.md',intro)
    overview='''# 231—775轮阶段成果总览

## 一、现在是否有阶段性成果

有。最恰当的阶段定位是：**从认知操作出发，得到一组有范围的空间重构定理和失败边界，并形成一个正在接通量子态、物质、参考、记录与几何反作用的联合候选。**

前两阶段001—230的有限维量子结论继续作为已有基础。231—775增加的是时空、动力学和共同实现问题；不能因报告数量增加，就把这些不同证据等级合并成一个已经完成的“认知必然推出宇宙”的证明。

## 二、四类可交付成果

### 1. 空间的条件性重构

382—386把完整反向、可逆重定向及实际位移邻域接到三维上下界与有限尺度证书。384已经消去旧下界中的额外Lipschitz要求。425给一致半幅与成本收缩到光滑坐标的条件连接。522—523在热参考和实际方向仪器上给共同实现。

这些是可复用的定理接口。尚未证明所有认知过程都必然具有它们的全部前提；内部qubit状态的三个分量也不能直接算作物理空间维数。

### 2. 真正限制候选的反例

342—344显示已列明的操作和几何条件未唯一决定自然演化与Einstein动力学；604排除指定中心差分物种接法；649排除原独立几何动能直接下降到指定面积匹配商；699排除声明的辅助候选对所有时间参数具有物理正性；759排除对全部输入统一有效的单确定背景轨道合同。

这些反例留下可检验的替代分支。它们的量词不同，不能合成对所有认知模型或所有量子引力路线的反证。

### 3. 同一材料承担多个角色的联合候选

531以后不再只追逐单个维数，而是联立表示、反常、质量、真空、约束、参考与尺度。548—573、753给出同一物质参考与经典动态几何相容的具体构造。574、598及其后续提供有限图量子过程；625、702—704、723等连接实际记录、热态与来源。

这比“分别找几个相似数学结构”更强：同一候选中的不同角色需要共同满足方程和资源域。但四维、群/物种、经典作用及若干参数仍属于输入，不是本阶段全部独立推出的结果。

### 4. 同一背景上的量子态与首阶来源

730—741完成声明费米部门与有限阶反作用的连接。764—767在已有共同背景上建立完整线性物理代数与正Hadamard态；768—772处理同态表示和完整首阶局部来源。773—775把局部修复、共同插入与实际来源首项接起来。

775关闭的是首项相容接口。774所需的高次共同规范化N1/N2、实际相互作用准备与记录仍未全部签收。因此可以把当前成果整理成阶段综述或范围明确的专题论文，不能命名为完整量子引力或全部物理推导已经完成。

## 三、为什么按15阶段分目录

分期依据是问题和验证对象的变化，不是每隔固定轮数切一刀。阶段01—02处理事件、网络和测量几何；03—04处理引力路线及推导边界；05—07处理实际空间接口与内部组织；08—11转入联合物理候选；12专门保留已经严格关闭的辅助正性路线；13—15收束到正过程、实际关联和共同量子场。

同一旧结论可以服务于多个后续阶段，故原文件只归属一个目录，另用跨阶段主题索引连接，避免复制出相互分叉的“最新版”。

## 四、恢复研究后的真实起点

下一编号仍为776。先核原混合体系在扩大但明确的局部系数类中的高次因果规范化及量子作用原理；若接通，复用773—775，不重做其首项证明。之后仍需共同实际准备/记录、跨尺度与图到连续的连接，以及独立输入和可区分预测的进一步压缩。

空间旧接口继续复用，不重新把已消去的条件列为缺口。目标保持暂停，归档完成不会自动恢复科学研究。

## 五、证据等级和本次核验

阶段综述概括历史报告的声明，不冒充逐行独立审稿。原解析论证、有限数值校准、成熟文献输入与物理解释保持区分。本次核对全部编号与文件迁移、链接和代表性复算，不把整理当作新增科学轮次。

详细分期、关键原文和完整545轮入口见[总目录](README.md)。原始字节、迁移映射及复算记录见[迁移说明](_migration/README.md)。
'''
    write(NEW/'阶段成果总览.md',overview)
    topics=[('三维与坐标',[372,382,383,384,386,425,522,523]),
            ('引力条件与非唯一性',[301,304,324,342,344,351,358]),
            ('内生交互、SoCA与传播',[429,451,459,465,466,504]),
            ('同一物质与经典几何',[531,548,549,551,572,573,753]),
            ('量子记录、热态与来源',[574,598,625,702,704,723,730,735,741]),
            ('应保留的限定失败',[604,649,699,759]),
            ('关联、量子背景与当前入口',[756,758,762,763,767,771,773,774,775])]
    text='# 跨阶段主题索引\n\n原文件只归属一个阶段；这里按后续复用问题交叉连接。\n\n'
    for title,ns in topics:
        text+='## '+title+'\n\n'+'、'.join(f'[{n}]({note_link(n,NEW,mapping)})' for n in ns)+'\n\n'
    text+='[总目录](README.md) · [待续776](next_round_776/STATUS.md)\n'
    write(NEW/'跨阶段主题索引.md',text)
    shared='# 跨阶段公共材料\n\n[总目录](../README.md)。下列文件没有唯一编号归属，或服务于多轮审计。原字节及脚本复算方式见[迁移说明](../_migration/README.md)。\n\n|文件|原路径|\n|---|---|\n'
    for e in plan['entries']:
        if e['destination'].startswith('_shared/'):
            shared+=f"|[{Path(e['destination']).name}]({relative(NEW/e['destination'],NEW/'_shared')})|`{e['original']}`|\n"
    write(NEW/'_shared/README.md',shared)
    write(NEW/'_history/README.md','''# 原导航与历史快照

[总目录](../README.md) · [迁移与复算](../_migration/README.md)

`indexes`保存原长篇索引和空间前提审计，链接已迁移。`navigation_snapshots`保存历次导航原始快照，字节完全不改；这些快照中的相对路径是当时布局的证据，不作为当前导航使用。重建原布局请使用统一复算入口。`runtime_cache`仅保留迁移前已有缓存，不作为科学证据。

迁移前的研究README、方向和状态原文包含在原工作区快照中；当前入口使用简洁的新导航。历史安排不覆盖用户暂停和最新目标。
''')


def apply():
    plan=json.loads((MIG/'plan.json').read_text('utf8'))
    assert not (MIG/'manifest.json').exists()
    actual={p.relative_to(OLD).as_posix() for p in OLD.rglob('*') if p.is_file()}
    assert actual=={e['original'] for e in plan['entries']}
    for e in plan['entries']:
        assert digest(OLD/e['original'])==e['original_sha256'], e['original']
        assert not (NEW/e['destination']).exists(), e['destination']
    snapshot(plan)
    mapping={(OLD/e['original']).resolve():(NEW/e['destination']).resolve() for e in plan['entries']}
    records=[]
    for e in plan['entries']:
        src,dst=OLD/e['original'],NEW/e['destination']
        assert src.resolve().is_relative_to(OLD.resolve())
        assert dst.resolve().is_relative_to(NEW.resolve())
        raw=src.read_bytes()
        assert sha(raw)==e['original_sha256']
        changed,edits=rewrite_links(raw,src,dst,mapping) if e['rewrite_markdown_links'] else (raw,[])
        dst.parent.mkdir(parents=True,exist_ok=True)
        # Atomic per-file relocation inside the verified workspace, not a shell-built recursive move.
        src.replace(dst)
        if edits:
            dst.write_bytes(changed)
        records.append(dict(**e,current_sha256=sha(changed),current_bytes=len(changed),link_edits=edits))
    document_phases(plan,mapping)
    external=[]
    for src in [ROOT/'README.md',RESEARCH/'archive_223_230/README.md']:
        before=src.read_bytes()
        after,edits=rewrite_links(before,src,src,mapping)
        if src==ROOT/'README.md':
            enc='utf-8-sig' if after.startswith(b'\xef\xbb\xbf') else 'utf8'
            text=after.decode(enc).replace('\r\n','\n')
            head,rest=text.split('\n\n',1) if '\n\n' in text else (text,'')
            text=head+'\n\n**231—775轮已按15阶段归档，目标暂停：** [阶段成果与完整目录](research_cognition_physics/archive_231_775/README.md)。545份编号报告、代码、结果及历史稿件已实际迁移；恢复研究从776入口继续。下方为保留的历史进展。\n\n'+rest
            after=text.encode(enc)
        src.write_bytes(after)
        external.append(dict(path=src.relative_to(ROOT).as_posix(),before_sha256=sha(before),after_sha256=sha(after),link_edits=edits))
    nav={
        'README.md': '''# 认知物理研究：阶段论文与研究档案

## 当前状态

目标由用户暂停。231—775轮已完成阶段回顾与实体归档；545份编号报告分为15个阶段。当前有条件性定理、限定反例和联合候选的阶段性成果，完整统一目标尚未完成。

- [231—775阶段成果与目录](archive_231_775/README.md)
- [阶段成果总览](archive_231_775/阶段成果总览.md)
- [跨阶段主题索引](archive_231_775/跨阶段主题索引.md)
- [研究方向](research_direction.md)
- [当前状态与恢复入口](RESEARCH_STATE.md)

## 已有阶段论文

- [可组合认知结构与复量子状态空间](可组合认知结构与复量子状态空间_阶段论文.md)：[001—222原始研究](archive_001_222/README.md)。
- [可组合认知结构与有限维量子理论](可组合认知结构与有限维量子理论_阶段论文.md)：[223—230原始研究](archive_223_230/README.md)。

## 档案与复算

231—775的原始报告、代码、结果已实际分拆；没有删掉旧稿，也没有把迁移计作新的科学轮次。历史校验中的旧路径与哈希通过原字节快照和统一复算入口保留，详见[迁移说明](archive_231_775/_migration/README.md)。旧archive_231_仅留重定向入口。
''',
        'research_direction.md': '''# 研究方向：认知本体论与现代物理的共同模型

## 目标保持

本次只整理231—775轮；不改写应用中的研究目标。继续追求同一基础结构在明确范围内共同承载量子理论、3+1时空、广义相对论、标准模型及测量记录。允许认知到物理与物理到认知的双向构建；认知观察用于提出假说，数学检查连接，物理预测检验适用性。

## 当前阶段成果

[15阶段目录](archive_231_775/README.md)与[成果总览](archive_231_775/阶段成果总览.md)区分条件性空间定理、具体路线的失败、共同经典模型、量子过程及局部来源连接。旧全部方向记录已保存在迁移原字节快照中，不能把历史“下一步”当作当前指令。

## 恢复后的顺序

1. 从[776已有入口](archive_231_775/next_round_776/STATUS.md)继续，先核同一原混合体系的高次因果规范化N1/N2。
2. 复用773—775的局部修复与首项匹配，不再要求先建完整in-in有效作用才能证明首项接口。
3. 接实际相互作用准备与记录、图到连续和共同跨尺度映射；继续列清独立物理输入及可区分后果。

空间、连续极限和跨尺度接口须回用382—386、425及522—523；384已消去的额外Lipschitz条件不得重新列为缺口。604、649、699、759的限定反例按各自量词保留。

## 执行约束

目标当前暂停，待用户恢复。不新增应用任务或定时任务。每个实质研究单元仍保存编号报告、可复算代码和结果；目录整理不新增编号。新阶段工作另建工作目录，不向旧分期混入未完成内容。历史科学文件、文献输入、解析证明、数值校准和物理解释分别记账；不把预置Einstein作用的验证写成认知独立生成引力。

复算方式见[迁移说明](archive_231_775/_migration/README.md)。
''',
        'RESEARCH_STATE.md': '''# 研究状态

## 2026-10-04：目标暂停，完成231—775阶段归档

- 最新完成轮次：775；下一编号：776，尚未开始新的科学轮次。
- 231—775共545份正式编号报告，已实体分拆到15个阶段目录。
- 最新历史科学核验记载累计3510项检查、1637份编号科学文件、3829份受保护证据。本次迁移检查另行统计，不能冒充全部历史实验重新运行。
- [阶段目录](archive_231_775/README.md)、[成果总览](archive_231_775/阶段成果总览.md)、[跨阶段索引](archive_231_775/跨阶段主题索引.md)。

## 当前真正接通的部分

原共同背景上的完整线性物理代数与正Hadamard态已经连接；同一自由态的完整首阶来源可作局部Ward修复。775进一步直接证明该来源与局部BV反常首项匹配，因而不必先有完整闭时路径有效作用，才能使用这个首项。

## 尚未接通

高次共同规范化N1/N2；实际相互作用准备、正性与记录；严格UV/全局拼接；原有限图过程到连续模型的共同映射。维数、群/物种、经典作用和若干参数仍按输入列账，不能声称全部物理已经从认知独立推出。

## 恢复入口

[776冻结入口](archive_231_775/next_round_776/STATUS.md)保存下一项。用户恢复目标前，只完成本次档案整理；不自动启动776。新轮次应采用新的工作目录，引用迁移后原报告。

## 证据和旧状态

所有原始字节及迁移前导航均有可恢复快照，原校验可通过[统一复算入口](archive_231_775/_migration/README.md)在隔离的旧布局中重跑。迁移后阅读版只更新链接，代码与结果不改。旧状态中的历史下一步不覆盖本状态。
'''
    }
    for name,text in nav.items():
        p=RESEARCH/name
        before=p.read_bytes()
        p.write_text(text,encoding='utf8',newline='\n')
        external.append(dict(path=p.relative_to(ROOT).as_posix(),before_sha256=sha(before),after_sha256=digest(p),navigation_compacted=True))
    # Only empty directories are removed; the original bytes are in the checked snapshot.
    for d in sorted((p for p in OLD.rglob('*') if p.is_dir()),key=lambda p:len(p.parts),reverse=True):
        assert d.resolve().is_relative_to(OLD.resolve())
        d.rmdir()
    write(OLD/'README.md','''# 研究档案已迁移

231—775轮已实际分拆到[15阶段档案](../archive_231_775/README.md)。此目录仅为旧入口提供定位，不再写入研究文件。

[阶段成果](../archive_231_775/阶段成果总览.md) · [复算与原始快照](../archive_231_775/_migration/README.md) · [待续776](../archive_231_775/next_round_776/STATUS.md)
''')
    write(MIG/'README.md','''# 迁移、冻结证据与复算

## 这次怎样迁移

原archive_231_中的全部文件已实际移动至15阶段、公共材料、历史快照及待续入口。原文件名保留。545份正式报告无缺号且各有唯一归属；草稿和导航历史保留。

科学代码、JSON结果和历史核验文件不改字节。Markdown阅读版只修改本地链接目标；数学与文字内容不改。每个变化的精确链接替换及迁移前后哈希记录在[manifest.json](manifest.json)。历史核验中的旧路径和哈希仍指迁移前原件，不能用当前阅读版哈希替代。

## 原始快照

[原工作区快照](original_workspace_231_775.zip)包含旧档案全部原始字节，以及重建其相对路径所需的工作区文本和前期档案。[快照核验](snapshot_checks.json)逐文件记录哈希。快照是恢复证据，不是日常重复维护的第二份研究。

## 统一复算入口

在此文件所在目录，用原Python运行时执行：

```powershell
& 'C:\\Users\\admin\\.cache\\codex-runtimes\\codex-primary-runtime\\dependencies\\python\\python.exe' -B -X utf8 .\\replay.py --script verify_round775.py
```

`--script`接受旧archive_231_中的脚本相对路径，也接受manifest中的新路径。额外参数在`--`之后提供。旧脚本默认路径相互依赖，请通过此入口执行，不在分拆后的code目录直接运行历史脚本。运行器在隔离临时目录恢复原布局，验证快照和所选脚本与当前代码一致，使用同一个Python与NumPy；结束清理自己的临时目录，不改正式档案。

历史报告中的旧运行命令作为历史记录保留。没有向1858份Python文件逐个注入兼容代码，也没有改科学算法以适配目录。以后新增代码应使用新路径或明确调用复算入口。

## 核验范围

[迁移核验](migration_checks.json)检查全部文件、允许的链接变化及545份报告；[代表性复算](replay_checks.json)单列实际重跑的脚本。未重跑全部545轮，未重审全部数学证明，未做图像检查。原本存在的历史断链另列，不伪称由迁移修好了旧事实。
''')
    manifest=dict(version=1,date='2026-10-04',rounds=[231,775],phases=PHASES,
                  original_file_count=len(records),entries=records,external_navigation_changes=external,
                  scientific_code_and_json_unchanged=True,mathematical_prose_unchanged=True,
                  goal_unchanged_and_paused=True,new_scientific_rounds=0)
    dump(MIG/'manifest.json',manifest)
    shutil.copyfile(ROOT/'scripts/research_phases_231_775.json',MIG/'phase_plan.json')
    print(json.dumps(dict(moved=len(records),phases=len(PHASES),reports=545,goal='paused'),ensure_ascii=False))


if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--plan',action='store_true')
    ap.add_argument('--apply',action='store_true')
    a=ap.parse_args()
    assert a.plan != a.apply
    make_plan() if a.plan else apply()
