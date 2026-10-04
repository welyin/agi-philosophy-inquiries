"""Read-only checks for rounds 333-334 and the unnumbered GR dependency review."""
import argparse
import json
import verify_effective_gravity_rounds as previous

CONFIG={333:('scalar_spatial_response_audit',12),
        334:('disformal_probe_geometry_audit',10)}
core=previous.core
core.CONFIG.update(CONFIG)
REVIEW='gr_assumption_dependency_review.md'


def verify(number,pending=False):
    result=previous.verify(number,pending)
    path=core.HERE/REVIEW
    checked=core.text_checks(path)
    assert checked['display_formulas']==1
    for link in core.link_parser()(path.read_text(encoding='utf-8')):
        assert (path.parent/link).resolve().exists(),link
    if number==333:
        result['new_file_hashes'][REVIEW]=core.digest(path)
    result['gr_dependency_review_hash_verified']=core.digest(path)
    result['dependency_review_is_not_a_scientific_test']=True
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('round',type=int,choices=CONFIG)
    parser.add_argument('--write-checks',action='store_true')
    args=parser.parse_args()
    result=verify(args.round,args.write_checks)
    if args.write_checks:
        path=core.HERE/f'research_round_{args.round}_checks.json'
        if path.exists():
            raise RuntimeError('Report already exists; use read-only mode.')
        path.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
