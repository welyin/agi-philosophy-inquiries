"""Merge the two completed SR stages and number five selected stage papers.

Only Markdown paths/navigation change. Frozen scientific files keep their bytes.
Every move is checked against the workspace and recorded with reversible edits.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess

from organize_research_231_775 import links
from organize_research_1044_1062 import text_edits
from research_layout import relocation_original_bytes

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'research_cognition_physics'
NEW=BASE/'archive_1106_1119'
MIG=BASE/'_migration/merge_1106_1119'
OLD={'archive_1106_':'1106_1115','archive_1116_':'1116_1119'}
PAPERS=[
    '可组合认知结构与复量子状态空间_阶段论文.md',
    '可组合认知结构与有限维量子理论_阶段论文.md',
    '认知操作五公理与三维关系空间_阶段论文.md',
    '认知时空操作公理与洛伦兹运动学_阶段论文.md',
    '认知动力学公理与狭义相对论基础动力学_阶段论文.md',
]
NUMBERED={name:f'{n:02}_{name}' for n,name in enumerate(PAPERS,1)}

def sha(raw):return hashlib.sha256(raw).hexdigest()
def rel(p):return p.relative_to(ROOT).as_posix()
def safe(p):
    p=p.resolve()
    assert p.is_relative_to(ROOT) and p!=ROOT,('Outside workspace',str(p))
    assert not set(p.relative_to(ROOT).parts)&{'.git','.codex','.agents','.research_runtime'}
    return p

def mapped(p):
    p=Path(os.path.abspath(p))
    if p.parent==BASE and p.name in NUMBERED:return BASE/NUMBERED[p.name]
    for old,segment in OLD.items():
        old_root=BASE/old
        if p==old_root:return NEW
        if p.is_relative_to(old_root):
            sub=p.relative_to(old_root)
            if sub.parts[0]=='_shared':return NEW/'_shared'/segment/Path(*sub.parts[1:])
            if sub.as_posix() in ('README.md','文件索引.md'):
                return NEW/'_shared'/segment/sub
            return NEW/sub
    return p

def plain_paths(text):
    for old,segment in OLD.items():
        for sep in ('/', '\\'):
            for name in ('README.md','文件索引.md'):
                text=text.replace(old+sep+name,'archive_1106_1119'+sep+'_shared'+sep+segment+sep+name)
            text=text.replace(old+sep+'_shared','archive_1106_1119'+sep+'_shared'+sep+segment)
        text=re.sub(re.escape(old)+r'(?!\d)', 'archive_1106_1119', text)
    for old,new in NUMBERED.items():
        text=re.sub(r'(?<!'+re.escape(new[:3])+r')'+re.escape(old),lambda _:new,text)
    return text

def paper_table(prefix=''):
    labels=['复量子状态空间（001—222）','有限维量子理论（223—230）',
            '三维关系空间（1063—1085）','Lorentz运动学（1086—1105）',
            '狭义相对论基础动力学（1116—1119）']
    out=['|序号|阶段论文|','|---|---|']
    for n,(old,new) in enumerate(NUMBERED.items(),1):
        out.append(f'|{n:02}|[{labels[n-1]}]({prefix}{new})|')
    return '\n'.join(out)

def notice():
    return ('## 阅读入口：阶段论文顺序与1106—1119合并归档\n\n'
      '[1106—1119统一档案](archive_1106_1119/README.md)已合并实际相对论连接与基础动力学两段研究。'
      '14份编号报告中，13轮正式验收、1113为保留的候选；本次整理不新增研究结论。\n\n'
      +paper_table()+'\n\n编号按理论推导与阅读层次排列；原研究轮次及各文采用的独立前提保持不变。'
      '另保留[条件生成链与物理选择的剩余自由](条件生成链与物理选择的剩余自由_阶段论文.md)，按用户要求暂不参与阶段论文编号，文件位置和内容保留。历史工作稿和冻结版本保留原名。'
      '旧路径及哈希收据保留原研究时点，通过[只读复算说明](_migration/merge_1106_1119/README.md)映射到唯一现存文件。\n\n---\n\n')

def rewrite(p,raw):
    if p.suffix.lower()!='.md':return raw,[]
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf-8'
    before=raw.decode(enc); after=before; dst=mapped(p)
    changes=[]
    for start,end,target,local in links(before):
        if re.match(r'^[A-Za-z]:',local):continue
        old_target=Path(os.path.abspath(p.parent/local))
        target_new=mapped(old_target)
        if target_new!=old_target or dst.parent!=p.parent:
            updated=os.path.relpath(target_new,dst.parent).replace('\\','/')
            if '#' in target:updated+='#'+target.split('#',1)[1]
            if target.startswith('<'):updated='<'+updated+'>'
            changes.append((start,end,updated))
    for start,end,updated in reversed(changes):after=after[:start]+updated+after[end:]
    after=plain_paths(after)
    if p.parent==BASE and p.name in {'README.md','RESEARCH_STATE.md','research_direction.md','ROADMAP.md'}:
        nl='\r\n' if '\r\n' in after else '\n'
        at=after.index('\n')+1
        after=after[:at]+nl+notice().replace('\n',nl)+after[at:]
    return after.encode(enc),text_edits(before,after)

def inventory():
    paths=set()
    for old in OLD:
        stage=BASE/old;assert stage.is_dir(),str(stage)
        paths.update(p for p in stage.rglob('*') if p.is_file())
    for name in PAPERS:
        p=BASE/name;assert p.is_file(),str(p);paths.add(p)
    files=subprocess.run(['rg','--files','-g','*.md','-g','!research_cognition_physics/_migration/**',
        '-g','!**/.git/**','-g','!**/.research_runtime/**'],cwd=ROOT,check=True,capture_output=True,encoding='utf-8').stdout.splitlines()
    for name in files:
        p=ROOT/name
        raw=p.read_bytes()
        if any(token.encode() in raw for token in [*OLD,*PAPERS]):paths.add(p)
    return sorted(paths)

def stage_readme():
    return '''# 1106—1119：实际狭义相对论连接与基础动力学

2026-10-11。按用户要求合并1106—1115与1116—1119两个已完成目标的档案。正式报告在本目录，证明、代码、结果与审阅在对应编号目录；跨轮材料按原段落各保留一份。

## 研究范围与阅读顺序

1. **1106—1115：同一有限有效域的实际连接。** 实际钟、原CO位置、未知Q及参考、完整物质反冲和已纳入影响使用共同Lorentz结构。[该段结项](_shared/1106_1115/stage_completion.md)。其中1113是保留的可选加强候选，不计正式验收轮次。
2. **1116—1119：认知动力学、惯性质量与质能关系。** 正式采用D-P，连接完整平移、静止内能、真实速度—动量、低速惯性与完整有效边界反作用。[正式采用](_shared/1116_1119/adoption_and_dependencies.md) · [该段结项](_shared/1116_1119/stage_completion.md)。

合计14份编号报告、13轮正式验收、81组构造复算（59＋22）；1113候选地位不变。新增科学校准0、累计3860。两个原目标均已完成，本次目录整理不新增轮次、公理或经验验证。

结论仍限各自公开的有限有效域和材料部门。绝对质量、材料谱、具体相互作用与更广物质覆盖仍有自由；不是全部物理理论完成。

## 入口

- [逐轮索引](文件索引.md)
- [05 狭义相对论基础动力学阶段论文](../05_认知动力学公理与狭义相对论基础动力学_阶段论文.md)
- [04 Lorentz运动学阶段论文](../04_认知时空操作公理与洛伦兹运动学_阶段论文.md)
- [1106—1115原段说明](_shared/1106_1115/README.md)及[1116—1119原段说明](_shared/1116_1119/README.md)
- [迁移与历史复算说明](../_migration/merge_1106_1119/README.md)

原结项、目标、依赖账分别保存在`_shared/1106_1115/`与`_shared/1116_1119/`，不覆盖同名文件。历史“当前／待审／未采用”保留原时点，当前状态以已完成结项及正式采用为准。

## 历史复算

从项目根目录使用`python -B -X utf8 scripts/run_research_merged.py --script research_cognition_physics/archive_1106_1119/1119/check.py`。统一入口把旧逻辑路径及原哈希与当前文件对应；不依赖旧目录、ZIP或文件副本。
'''

def stage_index():
    out=['# 1106—1119轮索引','','[阶段说明](README.md) · [研究总目录](../README.md)','',
         '|轮次|报告|状态|材料|','|---|---|---|---|']
    for n in range(1106,1120):
        old='archive_1106_' if n<=1115 else 'archive_1116_'
        lines=(BASE/old/f'research_note_{n}.md').read_text(encoding='utf-8-sig').splitlines()
        title=next(x[2:] for x in lines if x.startswith('# ')).replace('|','／')
        status='候选，未计正式验收' if n==1113 else '已验收'
        out.append(f'|{n}|[{title}](research_note_{n}.md)|{status}|[{n}/]({n}/)|')
    out+=['','## 分段材料','',
      '- [1106—1115原索引](_shared/1106_1115/文件索引.md)与[依赖账](_shared/1106_1115/dependency_ledger.json)',
      '- [1116—1119原索引](_shared/1116_1119/文件索引.md)与[依赖账](_shared/1116_1119/dependency_ledger.json)',
      '- [完整迁移清单](../_migration/merge_1106_1119/manifest.json)与[迁移说明](../_migration/merge_1106_1119/README.md)',
      '', '历史代码、结果和签审收据保持原字节；Markdown只修引用及维护入口。复算入口见阶段说明。','']
    return '\n'.join(out)

def build():
    assert not NEW.exists(),'Merged directory already exists'
    assert not (MIG/'manifest.json').exists(),'Migration already recorded'
    changes={};entries=[];destinations=set()
    for p in inventory():
        safe(p);dst=safe(mapped(p));raw=p.read_bytes();new,edits=rewrite(p,raw)
        if p==dst and raw==new:continue
        assert dst not in destinations,('Destination collision',str(dst));destinations.add(dst)
        assert p==dst or not dst.exists(),('Would overwrite',str(dst))
        entry=dict(original=rel(p),destination=rel(dst),original_sha256=sha(raw),current_sha256=sha(new),
            original_bytes=len(raw),current_bytes=len(new),edits=edits)
        assert relocation_original_bytes(new,entry)==raw
        if p.suffix.lower() in {'.py','.json'}:assert raw==new
        entries.append(entry);changes[p]=(dst,raw,new)
    created={NEW/'README.md':stage_readme().encode('utf-8'),NEW/'文件索引.md':stage_index().encode('utf-8')}
    support={rel(p):sha(p.read_bytes()) for p in (ROOT/'scripts').glob('*.py')
        if p.name not in {Path(__file__).name,'run_research_merged.py'}}
    support['research_cognition_physics/_migration/motion_active_1106/layer.json']=sha((BASE/'_migration/motion_active_1106/layer.json').read_bytes())
    manifest=dict(schema='research_merge_1106_1119_v1',date='2026-10-11',entries=entries,
        created_files=[dict(path=rel(p),sha256=sha(b),bytes=len(b)) for p,b in created.items()],unchanged_support_sha256=support,
        old_stages=[rel(BASE/x) for x in OLD],merged_stage=rel(NEW),
        paper_order=[dict(number=n,original=name,current=NUMBERED[name]) for n,name in enumerate(PAPERS,1)],
        scientific_files_byte_identical=True,zip_used=False,directory_copy_used=False,
        formal_rounds=[n for n in range(1106,1120) if n!=1113],candidate_rounds=[1113])
    return manifest,changes,created

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--apply',action='store_true');args=ap.parse_args()
    manifest,changes,created=build()
    # Compare each pre-existing usable local link by destination, not just text.
    checked_links=0
    for p,(dst,old,new) in changes.items():
        if p.suffix.lower()!='.md':continue
        old_links=list(links(old.decode('utf-8-sig')))
        new_links=list(links(new.decode('utf-8-sig')))
        # Root navigation prepends a new, separately checked reading notice.
        if p.parent==BASE and p.name in {'README.md','RESEARCH_STATE.md','research_direction.md','ROADMAP.md'}:
            new_links=new_links[len(list(links(notice()))):]
        assert len(old_links)==len(new_links),(rel(p),'Link inventory changed')
        for before,after in zip(old_links,new_links):
            original=Path(os.path.abspath(p.parent/before[3]))
            if not original.exists():continue
            expected=mapped(original)
            actual=Path(os.path.abspath(dst.parent/after[3]))
            assert actual==expected,(rel(p),before[2],after[2],str(expected),str(actual))
            checked_links+=1
    summary=dict(entries=len(changes),moved_files=sum(p!=d for p,(d,_,_) in changes.items()),
        markdown_updated=sum(old!=new for _,old,new in changes.values()),
        frozen_python_json=sum(p.suffix.lower() in {'.py','.json'} for p in changes),
        created_navigation=len(created),preserved_local_link_targets=checked_links,paper_order=manifest['paper_order'])
    print(json.dumps(summary,ensure_ascii=False))
    if not args.apply:return
    # Final concurrent-edit / absolute-path checks before any mutation.
    for p,(dst,raw,_) in changes.items():
        assert p.read_bytes()==raw,('Concurrent edit',rel(p));safe(p);safe(dst)
        if p!=dst:assert not dst.exists()
    for p in created:safe(p);assert not p.exists()
    MIG.mkdir(parents=True,exist_ok=True)
    with (MIG/'manifest.json').open('x',encoding='utf-8',newline='\n') as f:
        json.dump(manifest,f,ensure_ascii=False,indent=2);f.write('\n')
    for p,(dst,raw,new) in changes.items():
        dst.parent.mkdir(parents=True,exist_ok=True)
        if p!=dst:p.rename(dst)
        if raw!=new:dst.write_bytes(new)
    for p,data in created.items():p.write_bytes(data)
    for old in OLD:
        directory=safe(BASE/old)
        for sub in sorted(directory.rglob('*'),key=lambda p:len(p.parts),reverse=True):
            assert sub.is_dir(),('Unmoved file',str(sub));safe(sub).rmdir()
        directory.rmdir()
    print('APPLIED: moved original files; old stages absent; no copies or ZIPs')

if __name__=='__main__':main()
