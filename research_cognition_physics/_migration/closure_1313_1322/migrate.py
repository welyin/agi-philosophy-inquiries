"""Reversible closure of 1313--1322; only PowerShell performs directory rename.
Scientific Python and JSON bytes are never edited. New navigation stays outside.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import re
import sys
import uuid

ROOT=Path(__file__).resolve().parents[3]
BASE=ROOT/"research_cognition_physics"
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/"scripts"))
from research_layout import relocation_original_bytes

SOURCE="research_cognition_physics/archive_1313_"
DEST="research_cognition_physics/archive_1313_1322"
STAGES=[dict(source=SOURCE,destination=DEST,first_round=1313,last_round=1322)]
NAV=[BASE/name for name in ("README.md","research_direction.md","RESEARCH_STATE.md","ROADMAP.md")]
PAPER=BASE/"15_认知交接机制的物质角色选择与参数剩余自由_阶段论文.md"
PATTERN=re.compile(r"archive_1313_(?![0-9])")


def sha(raw):return hashlib.sha256(raw).hexdigest()
def dump(path,data):
    with path.open("x",encoding="utf-8",newline="\n") as stream:
        json.dump(data,stream,ensure_ascii=False,indent=2)
        stream.write("\n")
def rel(path):return path.relative_to(ROOT).as_posix()
def edited(raw,edits):
    if not edits:return raw
    enc="utf-8-sig" if raw.startswith(b"\xef\xbb\xbf") else "utf-8"
    text=raw.decode(enc)
    for item in reversed(edits):
        assert text[item["start"]:item["end"]]==item["before"]
        text=text[:item["start"]]+item["after"]+text[item["end"]:]
    return text.encode(enc)


def inventory():
    source=ROOT/SOURCE
    assert source.resolve().parent==BASE.resolve() and source.is_dir()
    assert not (ROOT/DEST).exists()
    files=[]
    for path in sorted(source.rglob("*")):
        assert not path.is_symlink(),path
        if path.is_file():
            raw=path.read_bytes()
            files.append(dict(path=rel(path),bytes=len(raw),sha256=sha(raw)))
    binding_count=0
    for number in range(1313,1323):
        acc=json.loads((source/str(number)/"acceptance.json").read_text(encoding="utf-8"))
        assert acc["status"].startswith("ACCEPTED"),number
        for name,expected in acc["sha256"].items():
            assert sha((ROOT/name).read_bytes())==expected,(number,name)
            binding_count+=1
    stage_path=source/"_shared/stage_acceptance_after1322.json"
    stage=json.loads(stage_path.read_text(encoding="utf-8"))
    assert stage["status"]=="ACCEPTED_DECLARED_STAGE_SCOPE"
    assert stage["scientific_rounds"]==list(range(1313,1323))
    for name,expected in stage["sha256"].items():
        assert sha((ROOT/name).read_bytes())==expected,("stage",name)
    before=dict(schema="closure_1313_1322_inventory_v1",date="2026-10-19",
                source=SOURCE,destination=DEST,files=files,
                scientific_python_json=sum(Path(x["path"]).suffix.lower() in (".py",".json") for x in files),
                round_bindings=binding_count,stage_bindings=len(stage["sha256"]),
                stage_acceptance_sha256=sha(stage_path.read_bytes()),
                directory_copy_used=False,zip_used=False)
    dump(HERE/"before.json",before)
    print(json.dumps({key:before[key] for key in ("scientific_python_json","round_bindings","stage_bindings")}|{"files":len(files),"status":"INVENTORY_PASS"}))


def prepare():
    before=json.loads((HERE/"before.json").read_text(encoding="utf-8"))
    files={ROOT/x["path"] for x in before["files"]}
    assert (ROOT/SOURCE).is_dir() and not (ROOT/DEST).exists()
    actual={p for p in (ROOT/SOURCE).rglob("*") if p.is_file()}
    assert actual==files,"Inventory changed"
    for item in before["files"]:
        raw=(ROOT/item["path"]).read_bytes()
        assert sha(raw)==item["sha256"] and len(raw)==item["bytes"],item["path"]
    external=set(NAV)|{PAPER}
    entries=[]
    for path in sorted(files|external):
        raw=path.read_bytes();name=rel(path)
        destination=DEST+name[len(SOURCE):] if name.startswith(SOURCE+"/") else name
        edits=[]
        if path.suffix==".md":
            enc="utf-8-sig" if raw.startswith(b"\xef\xbb\xbf") else "utf-8"
            text=raw.decode(enc)
            for match in PATTERN.finditer(text):
                edits.append(dict(start=match.start(),end=match.end(),before=match.group(),after="archive_1313_1322"))
            if name==SOURCE+"/README.md":
                at=text.index("\n")+1
                title=text[:at].rstrip("\r\n")
                if not title.startswith("# 1313—1322"):
                    edits.append(dict(start=0,end=len(title),before=title,after="# 1313—1322：已结项研究档案"))
                addition=("\n## 正式归档与后继\n\n"
                    "2026-10-19。十轮已在既定机制范围内通过逐轮和阶段验收，目录收尾为 archive_1313_1322。"
                    "原报告、代码、结果、审阅及冻结材料全部保留；科学 Python 与 JSON 逐字节不变。"
                    "Markdown路径与本段导航修改可逆。"
                    "[迁移清单及历史复算](../_migration/closure_1313_1322/README.md)。\n\n"
                    "后继从1323起按新阶段独立目标继续；本次归档不新建科学轮次、公理或扩大已验收结论。\n")
                edits.append(dict(start=at,end=at,before="",after=addition))
        edits.sort(key=lambda x:(x["start"],x["end"]))
        current=edited(raw,edits)
        entry=dict(original=name,destination=destination,original_sha256=sha(raw),
                   current_sha256=sha(current),original_bytes=len(raw),current_bytes=len(current),edits=edits)
        assert relocation_original_bytes(current,entry)==raw
        if path.suffix.lower() in (".py",".json"):
            assert current==raw
        entries.append(entry)
    layer=BASE/"_migration/interaction_generation_active_1313/layer.json"
    support=dict(json.loads(layer.read_text(encoding="utf-8"))["unchanged_support_sha256"])
    for p in (layer,ROOT/"scripts/run_research_interaction_generation_active.py"):
        support[rel(p)]=sha(p.read_bytes())
    for name,expected in support.items():assert sha((ROOT/name).read_bytes())==expected,name
    plan=dict(schema="research_closure_1313_1322_v1",date="2026-10-19",stages=STAGES,
              entries=entries,unchanged_support_sha256=support,moved_files=len(files),
              external_markdown_files=len(external),directory_copy_used=False,zip_used=False,
              scientific_results_rewritten=False)
    dump(HERE/"plan.json",plan)
    print(json.dumps(dict(status="PLAN_READY",moved_files=len(files),external_files=len(external),
                         edited_markdown=sum(bool(e["edits"]) for e in entries))))


def apply():
    plan=json.loads((HERE/"plan.json").read_text(encoding="utf-8"))
    assert not (ROOT/SOURCE).exists() and (ROOT/DEST).is_dir()
    for entry in plan["entries"]:
        raw=(ROOT/entry["destination"]).read_bytes()
        assert sha(raw) in (entry["original_sha256"],entry["current_sha256"]),entry["destination"]
    for entry in plan["entries"]:
        path=ROOT/entry["destination"];raw=path.read_bytes()
        if sha(raw)!=entry["current_sha256"]:
            current=edited(raw,entry["edits"])
            temporary=path.with_name(".closure-write-"+uuid.uuid4().hex)
            with temporary.open("xb") as stream:stream.write(current)
            os.replace(temporary,path)
        relocation_original_bytes(path.read_bytes(),entry)
    dump(HERE/"manifest.json",plan)
    assert (HERE/"plan.json").read_bytes()==(HERE/"manifest.json").read_bytes()
    (HERE/"plan.json").unlink()
    print(json.dumps(dict(status="MIGRATED",files=len(plan["entries"]),
                         manifest_sha256=sha((HERE/"manifest.json").read_bytes()))))


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("mode",choices=["inventory","prepare","apply"])
    {"inventory":inventory,"prepare":prepare,"apply":apply}[parser.parse_args().mode]()

