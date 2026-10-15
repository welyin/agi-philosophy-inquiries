"""Publish an explicitly accepted round; startup receipts remain historical.

Preview: --acceptance ../1219/acceptance.json --plan publication_1219.json
Publish: add --apply. Verify current publication: --verify progress_checks_1219.json
Relative CLI paths resolve beside this script. No scientific check is promoted
or counted without the supplied acceptance and explicit publication plan.
"""
from pathlib import Path
import argparse
import hashlib
import importlib
import json
import re
import subprocess
import sys

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
ROOT = HERE.parents[2]
BASE = ROOT / 'research_cognition_physics'
STARTUP = HERE / 'startup_checks.json'
LEDGER = HERE / 'dependency_ledger.json'
NOTICE = HERE / 'navigation_notice.md'
NAV_NAMES = ('README.md','research_direction.md','RESEARCH_STATE.md','ROADMAP.md')
LAYER = BASE / '_migration/thermal_geometry_active_1219/layer.json'
RUNNER = ROOT / 'scripts/run_research_thermal_geometry_active.py'
EXTERNAL_PAPER = 'research_cognition_physics/08_认知热学角色与热力学四律_阶段论文.md'
MUTABLE_STARTUP = {p.relative_to(ROOT).as_posix() for p in [LEDGER,NOTICE,
    HERE/'usage_reset_policy.json',STAGE/'README.md',STAGE/'文件索引.md']}

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p): return json.loads(p.read_text(encoding='utf-8-sig'))
def rel(p): return p.relative_to(ROOT).as_posix()
def local(value): return (HERE / value).resolve()
def encode(value): return (json.dumps(value,ensure_ascii=False,indent=2)+'\n').encode('utf-8')
def validate_map(bindings):
    for name, expected in bindings.items():
        assert sha(ROOT/name)==expected, 'Changed: '+name

def run(*args):
    result=subprocess.run([sys.executable,'-B','-X','utf8',*map(str,args)],
                          capture_output=True,text=True,encoding='utf-8')
    assert result.returncode==0,result.stderr+'\n'+result.stdout
    return json.loads(result.stdout)

def frozen_startup():
    data=read(STARTUP)
    assert data['accepted_scientific_rounds']==[] and data['new_scientific_rounds']==0
    frozen={name:h for name,h in data['source_files_sha256'].items() if name not in MUTABLE_STARTUP}
    validate_map(frozen)
    return data,frozen

def acceptance(path):
    data=read(path); number=data['round']
    assert type(number) is int and number>=1219
    assert path==(STAGE/str(number)/'acceptance.json').resolve()
    assert data['status'].startswith('ACCEPTED'), 'Acceptance is not signed'
    assert type(data['check_groups']) is int and data['check_groups']>0
    assert isinstance(data['files_sha256'],dict) and data['files_sha256']
    files={}
    for name,h in data['files_sha256'].items():
        p=(path.parent/name).resolve()
        assert p.is_relative_to(STAGE) and sha(p)==h,(name,'acceptance binding')
        files[rel(p)]=h
    note=STAGE/f'research_note_{number}.md'
    assert rel(note) in files,'Numbered note must be bound by acceptance'
    return data,files

def audit_bindings(plan):
    bindings={}
    for item in plan.get('management_audit_documents',[]):
        p=(ROOT/item['path']).resolve()
        assert p.is_relative_to(STAGE/'_admission') and sha(p)==item['sha256']
        bindings[rel(p)]=item['sha256']
    return bindings

def closure_check(plan,scientific_bindings):
    """Keep whole-closure failure explicit; never substitute old paper bytes.

    The only allowed deviation is the specifically registered current root
    paper. Every archive entry, other external entry, inventory and frozen
    support hash is still mandatory. The frozen runner remains unchanged.
    """
    allowed=plan.get('external_current_document_edits',[])
    if not allowed:
        result=run(RUNNER,'--verify-closure')
        assert result['status']=='PASS'
        return result
    assert len(allowed)==1 and allowed[0]['path']==EXTERNAL_PAPER
    permission=allowed[0]
    assert permission['not_a_current_executable_or_byte_bound_source'] is True
    assert EXTERNAL_PAPER not in scientific_bindings
    audit=(ROOT/permission['audit_path']).resolve()
    assert audit.is_relative_to(STAGE/'_admission') and sha(audit)==permission['audit_sha256']
    assert permission['audit_path'] in audit_bindings(plan)
    registered=read(audit)
    assert registered['paper']==EXTERNAL_PAPER
    assert registered['current_sha256']==permission['observed_sha256']
    assert registered['current_bytes']==permission['observed_bytes']
    assert registered['expected_closure_current_sha256']==permission['expected_closure_current_sha256']
    validate_map(registered['dependency_bindings_sha256'])
    # Results contain the direct scientific byte-bound source menu. Do not
    # allow a later round that names this paper to inherit the exception.
    source_bindings={}
    for name in scientific_bindings:
        if name.endswith('.json'):
            for source,h in read(ROOT/name).get('source_sha256',{}).items():
                full=source if source.startswith('research_cognition_physics/') else 'research_cognition_physics/'+source
                assert full!=EXTERNAL_PAPER,'External edit is now a scientific source'
                assert full not in source_bindings or source_bindings[full]==h
                source_bindings[full]=h
    validate_map(source_bindings)
    sys.path.insert(0,str(ROOT/'scripts'))
    try: module=importlib.import_module('run_research_thermal_geometry_active')
    finally: sys.path.pop(0)
    store=module.ThermalGeometryLayout()
    store.verify_active_navigation()
    closed=module.previous
    assert sha(closed.MANIFEST)==store.thermal_sha256
    expected=set();verified=scientific=markdown=0;external=[]
    for name,entry in store.thermal_entries.items():
        if closed.moved(name):expected.add(entry['destination'])
        if name==EXTERNAL_PAPER:
            assert not closed.moved(name) and entry['destination']==name and name.endswith('.md')
            assert entry['current_sha256']==permission['expected_closure_current_sha256']
            current=(ROOT/name).read_bytes()
            assert hashlib.sha256(current).hexdigest()==permission['observed_sha256']
            assert len(current)==permission['observed_bytes']
            assert permission['observed_sha256']!=entry['current_sha256']
            external.append(dict(path=name,status='EXTERNAL_EDIT',
                expected_sha256=entry['current_sha256'],expected_bytes=entry['current_bytes'],
                observed_sha256=permission['observed_sha256'],observed_bytes=len(current),
                original_bytes_substituted=False,audit_path=permission['audit_path'],audit_sha256=permission['audit_sha256']))
            continue
        prior=store.before_thermal(ROOT/name)
        verified+=1;markdown+=bool(entry['edits'])
        if Path(name).suffix.lower() in ('.py','.json'):
            assert hashlib.sha256(prior).hexdigest()==entry['current_sha256']
            scientific+=1
    assert len(external)==1
    inventory=set()
    for source,(destination,_,_) in closed.STAGES.items():
        assert not (ROOT/source).exists() and (ROOT/destination).is_dir()
        inventory.update(rel(p) for p in (ROOT/destination).rglob('*') if p.is_file())
    assert inventory==expected
    support=store.thermal_manifest['unchanged_support_sha256']
    assert isinstance(support,dict) and support
    validate_map(support)
    validate_map(store.thermal_manifest.get('new_support_sha256',{}))
    return dict(status='EXTERNAL_EDIT',whole_closure_passed=False,
        scope='Unchanged archived evidence verified; the current root paper is an explicitly unresolved external edit',
        manifest_sha256=store.thermal_sha256,total_manifest_entries=len(store.thermal_entries),
        unchanged_entries_verified=verified,moved_archive_files_verified=len(expected),
        scientific_python_json_byte_identical=scientific,reversible_markdown_files_verified=markdown,
        frozen_support_files_verified=len(support),scientific_source_bindings_verified=len(source_bindings),
        external_edits=external,old_directories_absent=True,archive_inventory_matches=True,
        historical_complete_replay_status='NOT_PASSED_EXTERNAL_DOCUMENT_EDIT',
        current_round_publication_allowed_by_explicit_plan=True,old_manifest_modified=False,
        old_scientific_files_modified=False,current_paper_modified=False,overlay_or_copy_created=False)

def replace_section(text,heading,body):
    pattern=re.escape(heading)+r'\n.*?(?=\n## |\Z)'
    assert re.search(pattern,text,re.S),heading
    return re.sub(pattern,lambda m:heading+'\n\n'+body.rstrip()+'\n',text,count=1,flags=re.S)

def prepare(acceptance_path,plan_path):
    plan=read(plan_path); ac,science=acceptance(acceptance_path); number=ac['round']
    out=HERE/f'progress_checks_{number}.json'
    assert not out.exists(),'Publication already exists; use --verify'
    initial,frozen=frozen_startup(); ledger=read(LEDGER)
    assert ledger['status']=='active' and not ledger['stage_requirement_status']['stage_deliverables_all_completed']
    assert sha(LEDGER)==plan['expected_ledger_sha256']
    assert sha(NOTICE)==plan['expected_notice_sha256']
    assert sha(LAYER)==plan['expected_navigation_layer_sha256']
    assert sha(acceptance_path)==plan['acceptance_sha256']
    previous=HERE/ledger['latest_management_receipt']
    assert previous.is_file() and sha(previous)==plan['previous_receipt_sha256']
    old=read(previous)
    if previous!=STARTUP:
        validate_map(old['management_sha256'])
        validate_map(old['science_files_sha256'])
        for previous_number,h in old['acceptance_sha256'].items():
            p=STAGE/previous_number/'acceptance.json'
            assert sha(p)==h
            acceptance(p)
    rounds=ledger['accepted_stage_rounds']
    assert number==(rounds[-1]+1 if rounds else ledger['first_round'])
    assert plan['accepted_stage_rounds']==rounds+[number]
    assert plan['stage_check_groups']==ledger['stage_check_groups']+ac['check_groups']
    assert ledger['goal_objective_sha256']==initial['goal_objective_sha256']
    requirements=dict(ledger['stage_requirement_status'])
    patch=plan.get('stage_requirement_status',{})
    assert set(patch)<=set(requirements) and patch.get('stage_deliverables_all_completed',False) is False
    assert not any(v is True or str(v).lower() in {'complete','completed','pass','passed'} for v in patch.values()),'Stage closure needs a separate goal audit'
    requirements.update(patch)
    ledger.update(latest_accepted_round=number,accepted_stage_rounds=plan['accepted_stage_rounds'],
        stage_check_groups=plan['stage_check_groups'],current_work=plan['current_work'],
        next_direction=plan['next_direction'],stage_requirement_status=requirements,
        current_direction_document=f'../{number}/NEXT.md',latest_management_receipt=out.name)
    status=(f"{plan['date']}。应用目标active；本阶段已验收{len(ledger['accepted_stage_rounds'])}轮、"
            f"{ledger['stage_check_groups']}组，最新为{number}。阶段尚未结项。")
    readme=(STAGE/'README.md').read_text(encoding='utf-8')
    first=readme.split('\n\n',2);assert len(first)==3
    readme=first[0]+'\n\n'+status+'\n\n'+first[2]
    readme=replace_section(readme,'## 当前任务',plan['current_work']+'\n\n下一步：'+plan['next_direction'])
    note_link=f'[第{number}轮报告](research_note_{number}.md)'
    current=f'{note_link} · [本轮验收]({number}/acceptance.json) · [当前进展核验](_shared/{out.name})'
    history=('启动检查与原导航验收是启动时点的历史证据，原字节保持；其中零轮状态及管理文件哈希不代表后续当前状态。'
             '后续管理状态由对应progress_checks核验，不覆盖startup_checks或反复复制启动文件。')
    if '## 当前验收与历史启动' in readme:
        readme=replace_section(readme,'## 当前验收与历史启动',current+'\n\n'+history)
    else:readme+='\n## 当前验收与历史启动\n\n'+current+'\n\n'+history+'\n'
    index=(STAGE/'文件索引.md').read_text(encoding='utf-8');first=index.split('\n\n',2)
    index=first[0]+'\n\n'+status+'\n\n'+first[2]
    index=re.sub(r'尚无已验收编号报告。[^\n]*','',index)
    index+='\n- '+current+'\n'
    notice=(f'## 当前阶段（1219起）：认知操作中的热力学—几何对应\n\n{status}\n\n'
        f'[阶段目录](archive_1219_/README.md) · [正式目标](archive_1219_/_shared/goal.md) · '
        f'[第{number}轮](archive_1219_/research_note_{number}.md) · '
        f'[进展核验](archive_1219_/_shared/{out.name}) · [索引](archive_1219_/文件索引.md)。\n\n'
        +plan['current_work']+'\n\n下一步：'+plan['next_direction']+'\n\n'
        '仍须在同一实际过程中完成对象/区域、时间—温度、独立边界—熵和完整来源/反作用的连接、有限预测、阶段论文与逐项审计。'
        '具体候选的条件结论不等于阶段目标完成；复用382—386、425及522—523，不重加已消去假设。\n\n'
        '[1161—1200](archive_1161_1200/README.md)的GR目标保持未完成；'
        '[1201—1218](archive_1201_1218/README.md)四律阶段保持声明范围内完成。'
        '历史复算入口为`scripts/run_research_thermal_geometry_active.py`，旧档案和迁移证据不回写。\n\n'
        '以下为前阶段归档及历史记录，其中“当前／下一步”属于原时点；现行工作以本栏和archive_1219_为准。\n')
    updates={LEDGER:encode(ledger),STAGE/'README.md':readme.encode('utf-8'),
             STAGE/'文件索引.md':index.encode('utf-8'),NOTICE:notice.encode('utf-8')}
    bindings=dict(old.get('science_files_sha256',{}));bindings.update(science)
    acceptances=dict(old.get('acceptance_sha256',{}));acceptances[str(number)]=sha(acceptance_path)
    audits=audit_bindings(plan)
    if audits:
        index+='\n## 当前管理审计\n\n'
        for item in plan['management_audit_documents']:
            target=(ROOT/item['path']).relative_to(STAGE).as_posix()
            index+=f"- [{item['label']}]({target})\n"
        updates[STAGE/'文件索引.md']=index.encode('utf-8')
    closure=closure_check(plan,bindings)
    data=dict(schema='thermal_geometry_round_progress_v1',date=plan['date'],status='PASS',
        status_scope='current accepted round and management publication, not unrestricted historical replay',
        application_goal_status='active',stage_complete=False,latest_round=number,
        accepted_stage_rounds=ledger['accepted_stage_rounds'],stage_check_groups=ledger['stage_check_groups'],
        new_scientific_rounds=1,goal_objective_sha256=initial['goal_objective_sha256'],
        startup_receipt_path=rel(STARTUP),startup_receipt_sha256=sha(STARTUP),
        startup_receipt_is_historical_not_current=True,
        startup_mutable_bindings=list(sorted(MUTABLE_STARTUP)),startup_immutable_bindings_sha256=frozen,
        previous_receipt_path=rel(previous),previous_receipt_sha256=sha(previous),
        publication_plan_sha256=sha(plan_path),publication_plan_path=rel(plan_path),
        publisher_sha256=sha(Path(__file__)),science_files_sha256=bindings,acceptance_sha256=acceptances,
        audit_files_sha256=audits,external_current_document_edits=plan.get('external_current_document_edits',[]),
        frozen_closure_verify=closure,historical_closure_status=closure['status'],
        stage_requirement_status=requirements,old_GR_goal_completed=False,
        four_laws_previous_goal_status='complete_in_declared_scope',new_goal_created=False)
    return out,updates,data

def verify(path):
    data=read(path)
    assert sha(Path(__file__))==data['publisher_sha256']
    assert sha(STARTUP)==data['startup_receipt_sha256']
    validate_map(data['startup_immutable_bindings_sha256'])
    validate_map(data['science_files_sha256']);validate_map(data['management_sha256'])
    validate_map(data.get('audit_files_sha256',{}))
    assert sha(ROOT/data['previous_receipt_path'])==data['previous_receipt_sha256']
    assert sha(ROOT/data['publication_plan_path'])==data['publication_plan_sha256']
    for number,h in data['acceptance_sha256'].items():
        p=STAGE/number/'acceptance.json';assert sha(p)==h;acceptance(p)
    ledger=read(LEDGER)
    assert ledger['accepted_stage_rounds']==data['accepted_stage_rounds']
    assert ledger['stage_check_groups']==data['stage_check_groups'] and ledger['latest_management_receipt']==path.name
    nav=run(HERE/'update_navigation.py','--check');assert nav['changed_files']==0
    assert run(RUNNER,'--verify-navigation')['status']=='PASS'
    closure=closure_check(read(ROOT/data['publication_plan_path']),data['science_files_sha256'])
    assert closure==data['frozen_closure_verify']
    return dict(status='PASS',latest_round=data['latest_round'],stage_rounds=len(data['accepted_stage_rounds']),
                stage_check_groups=data['stage_check_groups'],stage_complete=False,startup_preserved_as_historical=True,
                historical_closure_status=closure['status'])

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--acceptance');parser.add_argument('--plan')
    parser.add_argument('--apply',action='store_true');parser.add_argument('--verify')
    args=parser.parse_args()
    if args.verify:
        assert not args.apply and not args.acceptance and not args.plan
        print(json.dumps(verify(local(args.verify)),ensure_ascii=False));return
    assert args.acceptance and args.plan
    out,updates,data=prepare(local(args.acceptance),local(args.plan))
    if not args.apply:
        print(json.dumps(dict(status='READY',round=data['latest_round'],stage_check_groups=data['stage_check_groups'],
                              files_to_update=[rel(p) for p in updates],receipt=rel(out),
                              historical_closure_status=data['historical_closure_status']),ensure_ascii=False));return
    snapshot={p:p.read_bytes() for p in [*updates,*[BASE/n for n in NAV_NAMES],LAYER]}
    known=dict(snapshot)
    try:
        for p,value in updates.items():
            assert p.read_bytes()==snapshot[p],str(p)
            p.write_bytes(value);known[p]=value
        applied=run(HERE/'update_navigation.py','--apply')
        for p in [*[BASE/n for n in NAV_NAMES],LAYER]:known[p]=p.read_bytes()
        checked=run(HERE/'update_navigation.py','--check');assert checked['changed_files']==0
        nav=run(RUNNER,'--verify-navigation');assert nav['status']=='PASS'
        closure=closure_check(read(local(args.plan)),data['science_files_sha256'])
        data['management_sha256']={rel(p):sha(p) for p in snapshot}
        data['navigation_apply']=applied;data['navigation_check']=checked
        data['navigation_verify']=nav;data['frozen_closure_verify']=closure
        with out.open('x',encoding='utf-8',newline='\n') as stream:
            json.dump(data,stream,ensure_ascii=False,indent=2);stream.write('\n')
    except Exception:
        if out.exists():raise RuntimeError('Receipt exists: inspect it before any rollback')
        for p,value in reversed(list(snapshot.items())):
            assert p.read_bytes()==known[p],'Concurrent management change blocks rollback: '+str(p)
            if known[p]!=value:p.write_bytes(value)
        raise
    print(json.dumps(dict(status='PUBLISHED',round=data['latest_round'],receipt=rel(out),stage_complete=False,
                         historical_closure_status=closure['status']),ensure_ascii=False))

if __name__=='__main__':main()
