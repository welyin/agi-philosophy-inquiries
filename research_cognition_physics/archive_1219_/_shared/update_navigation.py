"""Publish one reversible 1219 navigation notice above the frozen closure."""
from pathlib import Path
import argparse
import hashlib
import json
import sys

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / 'research_cognition_physics'
HERE = Path(__file__).resolve().parent
NOTICE = HERE / 'navigation_notice.md'
LAYER = BASE / '_migration/thermal_geometry_active_1219/layer.json'
FROZEN = BASE / '_migration/closure_1161_1218/manifest.json'
NAMES = ('README.md','research_direction.md','RESEARCH_STATE.md','ROADMAP.md')
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import relocation_original_bytes

def sha(raw): return hashlib.sha256(raw).hexdigest()

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group()
    group.add_argument('--apply',action='store_true')
    group.add_argument('--check',action='store_true')
    args=parser.parse_args()
    frozen_raw=FROZEN.read_bytes()
    frozen=json.loads(frozen_raw)
    assert frozen['schema']=='research_closure_1161_1218_v1'
    boundary={e['destination']:e for e in frozen['entries']}
    notice_raw=NOTICE.read_bytes()
    notice=notice_raw.decode('utf-8-sig').replace('\r\n','\n').strip()
    assert notice.startswith('## ') and notice.count('$$')%2==0
    old_raw=LAYER.read_bytes() if LAYER.exists() else None
    old=json.loads(old_raw) if old_raw else None
    if old: assert old['schema']=='thermal_geometry_active_navigation_v1'
    existing={e['original']:e for e in old['entries']} if old else {}
    keys={(BASE/n).relative_to(ROOT).as_posix() for n in NAMES}
    if old: assert set(existing)==keys and len(old['entries'])==4
    entries,changes=[],[]
    for name in NAMES:
        p=BASE/name; key=p.relative_to(ROOT).as_posix(); current=p.read_bytes()
        original=relocation_original_bytes(current,existing[key]) if old else current
        assert (sha(original),len(original))==(boundary[key]['current_sha256'],boundary[key]['current_bytes']),name
        enc='utf-8-sig' if original.startswith(b'\xef\xbb\xbf') else 'utf-8'
        text=original.decode(enc); assert text.startswith('# ') and '\n' in text
        newline='\r\n' if '\r\n' in text else '\n'; split=text.index('\n')+1
        insert=newline+notice.replace('\n',newline)+newline+newline
        after=(text[:split]+insert+text[split:]).encode(enc)
        entry=dict(original=key,destination=key,original_sha256=sha(original),current_sha256=sha(after),
                   original_bytes=len(original),current_bytes=len(after),
                   edits=[dict(start=split,end=split,before='',after=insert,kind='active_stage_navigation_notice')])
        assert relocation_original_bytes(after,entry)==original
        entries.append(entry);changes.append((p,current,after))
    supports=['scripts/run_research_thermodynamics_closed.py',FROZEN.relative_to(ROOT).as_posix(),
              'research_cognition_physics/_migration/gr_active_1161/layer.json']
    support={n:sha((ROOT/n).read_bytes()) for n in supports}
    if old: assert old['unchanged_support_sha256']==support,'Frozen support changed'
    layer=dict(schema='thermal_geometry_active_navigation_v1',date='2026-10-15',
               notice_path=NOTICE.relative_to(ROOT).as_posix(),notice_sha256=sha(notice_raw),
               entries=entries,unchanged_support_sha256=support)
    if args.apply:
        assert FROZEN.read_bytes()==frozen_raw and NOTICE.read_bytes()==notice_raw
        assert (LAYER.read_bytes() if LAYER.exists() else None)==old_raw
        assert all(p.read_bytes()==current for p,current,after in changes)
        LAYER.parent.mkdir(parents=True,exist_ok=True)
        written=[]
        try:
            for p,current,after in changes:
                if current!=after:
                    p.write_bytes(after);written.append((p,current,after))
            LAYER.write_text(json.dumps(layer,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
        except Exception:
            for p,current,after in reversed(written):
                if p.read_bytes()==after:p.write_bytes(current)
            if old_raw is not None:LAYER.write_bytes(old_raw)
            raise
    print(json.dumps(dict(status='APPLIED' if args.apply else 'CHECKED',files=4,
        changed_files=sum(a!=b for p,a,b in changes),closed_boundary_unchanged=True,
        original_navigation_recoverable=True,scientific_files_modified=False),ensure_ascii=False))

if __name__=='__main__':main()
