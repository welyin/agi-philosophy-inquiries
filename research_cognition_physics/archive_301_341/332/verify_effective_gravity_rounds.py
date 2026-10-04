"""Scientific and text audit for rounds 331-332, including additive round-330 errata."""
import argparse
import json
import re
import verify_heavy_scalar_rounds as previous

CONFIG={331:('heavy_scalar_geometry_audit',10),
        332:('occupied_scalar_gravity_audit',12)}
core=previous.core
core.CONFIG.update(CONFIG)


def maintenance():
    manifest=core.HERE/'round330_text_correction.json'
    record=core.read(manifest)
    original=core.HERE/record['original']; corrected=core.HERE/record['corrected']
    assert core.digest(original)==record['original_sha256']
    assert core.digest(corrected)==record['corrected_sha256']
    raw=original.read_bytes(); bad=bytes([13])+b'm '
    assert raw.count(bad)==record['recovered_backslashes']==4
    assert corrected.read_bytes()==record['banner'].encode('utf-8')+raw.replace(bad,b'\\rm ')
    assert not re.search(rb'\r(?!\n)',corrected.read_bytes())
    assert core.text_checks(corrected)['display_formulas']==6
    return {p.name:core.digest(p) for p in (manifest,corrected)}


def verify(number,pending=False):
    hashes=maintenance()
    result=previous.verify(number,pending)
    assert not re.search(rb'\r(?!\n)',(core.HERE/f'research_note_{number}.md').read_bytes())
    if number==331:
        result['new_file_hashes'].update(hashes)
    result['additive_text_maintenance_hashes_verified']=hashes
    result['text_maintenance_is_not_a_scientific_test']=True
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
