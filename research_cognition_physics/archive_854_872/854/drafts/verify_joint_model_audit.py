"""Check the non-numbered 853 stage audit and its exact source files."""
from pathlib import Path
import argparse,hashlib,json,re
HERE=Path(__file__).resolve().parent;STAGE=HERE.parents[1];ROOT=HERE.parents[3]
REPORT=HERE/'共同模型阶段报告_截至853.md';RECEIPT=HERE/'joint_model_audit_853_receipt.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(writing=False):
    body=REPORT.read_text('utf-8');assert body.count('$$')==2 and '\t' not in body
    assert '不算完成854' in body and '完整持续目标不结项' in body
    refs={};links=0
    for target in re.findall(r'\]\(([^)]+)\)',body):
        if re.match(r'^[a-zA-Z]+://',target) or target.startswith('#'):continue
        path=(REPORT.parent/target.split('#')[0].strip('<>')).resolve()
        assert path.exists() or (writing and path==RECEIPT.resolve()),target
        links+=1
        if path!=RECEIPT.resolve() and path.is_file() and path.name!='unified_physics_condition_ledger_current.md':
            refs[str(path.relative_to(ROOT))]=sha(path)
    return dict(kind='non_numbered_joint_model_stage_audit',date='2026-10-05',formal_reports=853,
        cumulative_numbered_test_groups=3638,new_scientific_test_groups=0,
        report_path=str(REPORT.relative_to(ROOT)),report_sha256=sha(REPORT),local_links_checked=links,
        reference_files=refs,review_method='Current ledger plus directly read key reports and preserved verification evidence, not a new exhaustive proof review of all historical notes.',
        scoped_local_formal_joint_model_supported=True,
        full_graph_continuum_mapping_or_physical_parameter_fit_completed=False,
        full_goal_completed=False,app_goal_changed=False,visual_checks_performed=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    if args.write:assert not RECEIPT.exists()
    result=run(args.write)
    if args.write:RECEIPT.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert result==json.loads(RECEIPT.read_text('utf-8'))
    print(json.dumps({k:v for k,v in result.items() if k!='reference_files'},ensure_ascii=False,indent=2))
