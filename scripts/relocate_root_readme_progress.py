"""Keep the user's exact README passage and relocate the surrounding logs."""
from collections import Counter,defaultdict
import argparse
import json
from pathlib import Path
import re
import zipfile
import organize_research_231_775 as m

ROOT,BASE=m.ROOT,m.RESEARCH
MAINT=BASE/'archive_764_/_migration/root_readme_simplification_20261004'
EARLY=json.loads((ROOT/'scripts/research_phases_001_222.json').read_text('utf8'))['phases']
LATE=json.loads((ROOT/'scripts/research_phases_231_775.json').read_text('utf8'))['phases']


def archive(n):
    if 1<=n<=222:
        p=next(p for p in EARLY if p['start']<=n<=p['end'])
        return BASE/f"archive_{p['start']:03d}_{p['end']:03d}"
    if 223<=n<=230:
        return BASE/'archive_223_230'
    if n>=764:
        return BASE/'archive_764_'
    if 231<=n<=763:
        p=next(p for p in LATE if p['start']<=n<=p['end'])
        return BASE/f"archive_{p['start']}_{p['end']}"
    return None


def make_plan():
    with zipfile.ZipFile(MAINT/'before_simplification.zip') as z:
        raw=z.read('README.md')
    source=raw.decode('utf-8-sig')
    start=source.index('## 写在前面')
    signature='—— 一个曾经在这里想过这些问题的人'
    end=source.index(signature,start)+len(signature)
    retained=source[start:end]+'\n'
    entries=[]
    for side,chunk,offset in [('before',source[:start],0),('after',source[end:],end)]:
        previous=None
        cursor=0
        for part in re.split(r'(\r?\n[ \t]*\r?\n)',chunk):
            part_start=offset+cursor
            cursor+=len(part)
            if not part.strip():
                continue
            # Leading/trailing newlines are formatting; retain every non-space
            # character and record its exact range in the original README.
            left=len(part)-len(part.lstrip())
            text=part.strip()
            a,b=part_start+left,part_start+left+len(text)
            assert source[a:b]==text
            header=re.sub(r'^[>\s*#]+','',text).split('：',1)[0].split(':',1)[0]
            match=re.match(r'第\s*(\d{1,3})\s*轮',header)
            if not match:
                match=re.match(r'(\d{3})(?=\D|$)',header)
            if not match:
                match=re.search(r'当前执行顺序[（(](\d{3})后',header)
            if re.match(r'\d{1,3}[—–-]\d{1,3}',header):
                match=None
            n=int(match[1]) if match else None
            target=archive(n) if n else None
            reason='numbered paragraph heading' if target else None
            if header.startswith(('001—775轮已统一','运行方式','当前阶段','最新范围')):
                target,reason=BASE,'project-wide historical direction or navigation'
            if header.startswith('核心定理'):
                target,reason=archive(222),'early-stage theorem summary'
            if target is None:
                destinations=[]
                for _,_,_,local in m.links(text):
                    resolved=(ROOT/local).resolve()
                    if resolved.is_relative_to(BASE):
                        rel=resolved.relative_to(BASE)
                        if len(rel.parts)>1 and rel.parts[0].startswith('archive_'):
                            destinations.append(BASE/rel.parts[0])
                    elif resolved.is_relative_to(ROOT/'research_physics_construction'):
                        destinations.append(ROOT/'research_physics_construction')
                    elif resolved.is_relative_to(ROOT/'research_information_geometry'):
                        destinations.append(ROOT/'research_information_geometry')
                if destinations:
                    target=destinations[0]
                    reason='first research artifact linked by the original paragraph'
                elif previous and not text.startswith('#'):
                    target=previous
                    reason='continuation of the preceding historical entry'
                else:
                    target=BASE
                    reason='project-wide historical navigation'
            previous=target
            path=target/'README.md'
            assert path.is_file(),path
            entries.append(dict(side=side,start=a,end=b,original_sha256=m.sha(text.encode('utf8')),
                                destination=path.relative_to(ROOT).as_posix(),attribution=reason,text=text))
    # Nothing outside the selected passage is dropped apart from whitespace.
    for side,chunk in [('before',source[:start]),('after',source[end:])]:
        original=re.sub(r'\s','',chunk)
        moved=''.join(re.sub(r'\s','',e['text']) for e in entries if e['side']==side)
        assert original==moved
    return raw,source,retained,entries


def apply():
    original,source,retained,entries=make_plan()
    groups=defaultdict(list)
    for entry in entries:
        groups[ROOT/entry['destination']].append(entry)
    existing={p:p.read_bytes() for p in groups}
    existing[ROOT/'README.md']=(ROOT/'README.md').read_bytes()
    existing[BASE/'README.md']=(BASE/'README.md').read_bytes()
    existing[MAINT/'README.md']=(MAINT/'README.md').read_bytes()
    backup=MAINT/'before_boundary_correction.zip'
    with zipfile.ZipFile(backup,'x',zipfile.ZIP_DEFLATED) as z:
        for p,raw in sorted(existing.items()):
            z.writestr(p.relative_to(ROOT).as_posix(),raw)
    with zipfile.ZipFile(backup) as z:
        for p,raw in existing.items():
            assert z.read(p.relative_to(ROOT).as_posix())==raw
    outputs={}
    for p,items in groups.items():
        before=existing[p].decode('utf-8-sig')
        if p==BASE/'README.md':
            before=before.replace('[项目根README](../README.md)：项目简介与一级目录导航。',
                                  '[项目根README](../README.md)：用户指定保留的项目说明与作者寄语。')
            before=before.replace('原项目根README中的逐轮进展和作者文字已保存在[历史备份](archive_764_/_migration/root_readme_simplification_20261004/README.md)。详细研究介绍由所属目录维护，后续不再向项目根README追加轮次日志。',
                                  '项目根README保留从“写在前面”到作者署名的完整原文。此前及此后的逐轮进展已移入各所属阶段README的历史记录区；[迁移与原文备份](archive_764_/_migration/root_readme_simplification_20261004/README.md)保留逐段来源。后续研究进展继续写入所属目录。')
        assert '## 原项目根README迁入的历史记录' not in before
        addition='\n\n## 原项目根README迁入的历史记录\n\n以下保留当时的进展、核验与后续安排，按原文出现顺序整理；它们是历史记录，当前状态以本页前文及[研究状态]('+m.relative(BASE/'RESEARCH_STATE.md',p.parent)+')为准。\n\n<details>\n<summary>展开历史记录</summary>\n\n'
        for entry in items:
            raw=entry['text'].encode('utf8')
            rewritten,edits=m.rewrite_links(raw,ROOT/'README.md',p,{})
            content=rewritten.decode('utf8')
            addition+='<!-- 原项目README字符区间 '+str(entry['start'])+':'+str(entry['end'])+' -->\n\n'
            offset=len(before+addition)
            addition+=content+'\n\n'
            entry.update(destination_start=offset,destination_end=offset+len(content),
                         current_sha256=m.sha(rewritten),link_edits=edits)
        addition+='</details>\n'
        outputs[p]=before.rstrip()+addition
        # Offsets above include the original trailing whitespace; normalize it
        # only after extracting exact markers below, so checks remain explicit.
        for entry in items:
            marker='<!-- 原项目README字符区间 '+str(entry['start'])+':'+str(entry['end'])+' -->\n\n'
            start=outputs[p].index(marker)+len(marker)
            changed,_=m.rewrite_links(entry['text'].encode('utf8'),ROOT/'README.md',p,{})
            entry['destination_start']=start
            entry['destination_end']=start+len(changed.decode('utf8'))
    outputs[ROOT/'README.md']=retained
    outputs[MAINT/'README.md']='''# 项目根README整理记录

2026-10-04。按用户最新指定，根README完整保留从“写在前面”至“一个曾经在这里想过这些问题的人”的原文。边界之前及之后的内容已逐段移到所属研究目录README的历史记录区。

[最初原文快照](before_simplification.zip)保存整理前完整README；[边界修正前快照](before_boundary_correction.zip)保留中间导航版本。历史记录只调整相对链接，不改研究正文；当时的“下一步”不覆盖当前研究状态。

[逐段迁移与核验](boundary_correction_checks.json)记录原文区间、归属、目标区间和链接变化。先前29行精简版本的核验仍作为中间记录保留，不代表最终根README格式。科学代码、结果及正式编号报告未改写。

[研究主目录](../../../README.md) · [项目根README](../../../../README.md)
'''
    for p,text in outputs.items():
        assert p.read_bytes()==existing[p]
        p.write_text(text,encoding='utf8',newline='\n')
    verify_relocation()


def verify_relocation():
    original,source,retained,entries=make_plan()
    assert (ROOT/'README.md').read_bytes().decode('utf8')==retained
    paths={ROOT/e['destination'] for e in entries}|{ROOT/'README.md',BASE/'README.md',MAINT/'README.md'}
    outputs={p:p.read_bytes().decode('utf8') for p in paths}
    for e in entries:
        p=ROOT/e['destination']
        marker='<!-- 原项目README字符区间 '+str(e['start'])+':'+str(e['end'])+' -->\n\n'
        assert outputs[p].count(marker)==1
        a=outputs[p].index(marker)+len(marker)
        expected,edits=m.rewrite_links(e['text'].encode('utf8'),ROOT/'README.md',p,{})
        e.update(destination_start=a,destination_end=a+len(expected.decode('utf8')),
                 current_sha256=m.sha(expected),link_edits=edits)
        current=outputs[p][e['destination_start']:e['destination_end']]
        # Explicitly reproduce each migrated paragraph from its original text.
        origin=source[e['start']:e['end']]
        for edit in reversed(e['link_edits']):
            assert origin[edit['start']:edit['end']]==edit['before']
            origin=origin[:edit['start']]+edit['after']+origin[edit['end']:]
        assert origin==current and m.sha(current.encode('utf8'))==e['current_sha256']
        del e['text']
    broken=[];count=0
    check_path=MAINT/'boundary_correction_checks.json'
    for p,text in outputs.items():
        for _,_,target,local in m.links(text):
            count+=1
            destination=(p.parent/local).resolve()
            if destination==check_path:
                continue
            if not destination.exists():
                broken.append(dict(file=p.relative_to(ROOT).as_posix(),target=target))
    result=dict(date='2026-10-04',user_selected_passage_preserved_verbatim=True,
                retained_start='## 写在前面',retained_end='—— 一个曾经在这里想过这些问题的人',
                retained_lines=len(retained.splitlines()),retained_sha256=m.sha(retained.encode('utf8')),
                moved_paragraphs=len(entries),all_outside_nonwhitespace_content_preserved=True,
                destinations=dict(Counter(e['destination'] for e in entries)),entries=entries,
                modified_readme_sha256={p.relative_to(ROOT).as_posix():m.digest(p) for p in outputs},
                local_links_checked=count,broken_links=broken,scientific_code_results_reports_changed=False,
                new_scientific_rounds=0,all_checks_passed=not broken)
    m.dump(check_path,result)
    print(json.dumps({k:result[k] for k in ('retained_lines','moved_paragraphs','destinations','local_links_checked','broken_links','all_checks_passed')},ensure_ascii=False,indent=2))
    assert not broken


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--apply',action='store_true')
    parser.add_argument('--verify',action='store_true')
    args=parser.parse_args()
    if args.verify:
        verify_relocation()
    elif args.apply:
        apply()
    else:
        _,_,retained,entries=make_plan()
        print(json.dumps(dict(retained_lines=len(retained.splitlines()),moved_paragraphs=len(entries),destinations=dict(Counter(e['destination'] for e in entries)),global_entries=[e['text'][:160] for e in entries if e['destination']=='research_cognition_physics/README.md']),ensure_ascii=False,indent=2))
