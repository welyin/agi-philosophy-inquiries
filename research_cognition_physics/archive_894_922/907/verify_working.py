"""907 working checkpoint; no formal round or scientific group is added."""
from pathlib import Path
import argparse,ast,hashlib,json,re
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
TARGET=HERE/'drafts/working_checks.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(writing=False):
    receipt=json.loads((STAGE/'906/research_round_906_checks.json').read_text('utf-8'))
    prior={}
    for key in ('frozen_inputs','new_scientific_and_entry_files'):
        for name,d in receipt[key].items():assert sha(ROOT/name)==d,name;prior[name]=d
    raw=json.loads((HERE/'material_reference_jets_results.json').read_text('utf-8'))
    cross=json.loads((HERE/'reference_chart_crosscheck_results.json').read_text('utf-8'))
    assert raw['formal_previous']==cross['formal_previous']==906
    assert raw['rows'][0]['data']['new']['determinant']>0
    assert cross['new_determinant_sign_consistent_across_refinements'] and cross['old_determinant_sign_consistent_across_refinements']
    assert [(r['N'],r['steps']) for r in cross['rows']]==[(16,4),(24,2),(24,4),(32,4)]
    assert not (STAGE/'research_note_907.md').exists()
    paths=[HERE/n for n in ('material_reference_jets.py','material_reference_jets_results.json','reference_chart_crosscheck.py','reference_chart_crosscheck_results.json','drafts/STATUS.md','drafts/material_source_implementation_audit.md')]+[Path(__file__),STAGE/'906/research_round_906_checks.json',STAGE/'773/joint_material_local_brst.py',STAGE/'863/smeared_loop_current_probe.py',STAGE/'863/physical_relational_loop_bridge.py',STAGE/'904/coupled_boson_tangent.py']
    for p in paths:
        if p.suffix=='.py':ast.parse(p.read_text('utf-8'),filename=str(p))
    docs=[HERE/'drafts/material_source_implementation_audit.md',HERE/'drafts/STATUS.md',STAGE/'_shared/notes/unified_physics_condition_ledger_current.md']
    docs += [STAGE/p for p in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md')]
    docs += [STAGE.parent/p for p in ('README.md','research_direction.md','RESEARCH_STATE.md')]
    count=0
    for p in docs:
        text=re.sub(r'\$\$.*?\$\$','',p.read_text('utf-8-sig'),flags=re.S)
        for target in re.findall(r'\]\(([^)]+)\)',text):
            if re.match(r'^[a-zA-Z]+://',target) or target.startswith('#'):continue
            q=(p.parent/target.split('#')[0].strip('<>')).resolve()
            assert q.exists() or (writing and q==TARGET.resolve()),(p,target)
            count+=1
    return dict(date='2026-10-06',status='working',all_checkpoint_checks_passed=True,formal_rounds=906,cumulative_numbered_groups=3691,new_scientific_groups=0,
        prior906_frozen_files_verified=len(prior),local_links_checked=count,
        source_time_derivative_checks_saved=True,actual_diffeomorphism_checks_saved=True,independent_grid_time_crosschecks_saved=True,
        continuous_chart_failure_proved=False,whole_common_model_failed=False,physical_source_implementation_complete=False,
        files={str(p.relative_to(ROOT)):sha(p) for p in paths},full_goal_completed=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args();out=run(a.write)
    if a.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:
        old=json.loads(TARGET.read_text('utf-8'));assert old['files']==out['files']
    print(json.dumps({k:v for k,v in out.items() if k!='files'},ensure_ascii=False,indent=2))
