"""Recompute round232 and check source, formula and local-link evidence."""
import argparse
import ast
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import unittest

HERE=Path(__file__).resolve().parent
BASE=HERE.parent
REPORT=HERE/'research_round_232_checks.json'
def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p): return json.loads(p.read_text(encoding='utf-8'))


def verify(pending_report=False):
    import record_composition_audit as science
    primary=['record_composition_audit.py','record_composition_audit_results.json',
             'research_note_232.md']
    if REPORT.exists():
        previous=read(REPORT)
        for name,sha in previous['new_file_hashes'].items():
            assert digest(HERE/name)==sha,name
    for p in HERE.glob('*.py'): ast.parse(p.read_text(encoding='utf-8'))
    saved=read(HERE/'record_composition_audit_results.json')
    computed=json.loads(json.dumps(science.report()))
    # Runtime changes do not change the scientific certificate.
    assert {k:v for k,v in saved.items() if k not in ('checks','runtime')} == {
        k:v for k,v in computed.items() if k!='runtime'}
    output=io.StringIO()
    result=unittest.TextTestRunner(stream=output).run(
        unittest.defaultTestLoader.loadTestsFromModule(science))
    assert result.wasSuccessful() and result.testsRun==8,output.getvalue()
    assert saved['checks']=={'run':8,'failures':0,'errors':0}
    math=read(HERE/'round232_math_checks.json')
    assert math['all_passed']
    for item in math['documents']: assert digest(HERE/item['path'])==item['sha256']
    spec=importlib.util.spec_from_file_location('old_link_parser',BASE/'archive_001_222/verify_stage1.py')
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    links=0
    for p in [HERE/'research_note_232.md', HERE/'README.md']:
        source=p.read_text(encoding='utf-8')
        assert not any(ord(c)<32 and c not in '\r\n\t' for c in source)
        for target in mod.local_links(source):
            dest=(p.parent/target).resolve()
            if pending_report and dest==REPORT and not dest.exists(): continue
            assert dest.exists(),(p.name,target)
            links+=1
    return {'date':'2026-09-22','round':232,'new_tests':saved['checks'],
            'fresh_tests_run':result.testsRun,'saved_scientific_results_reproduced':True,
            'scientific_results_rewritten':False,'math_documents_checked':len(math['documents']),
            'local_links_checked':links,'broken_links':0,
            'new_file_hashes':{name:digest(HERE/name) for name in primary},
            'explicit_extra_input':'Known single-use Boolean environment plus explicit independent-copy and port-regrouping interface.',
            'not_claimed':'Spacetime, physical clocks, a universal speed, mass, gauge group or gravity.',
            'all_reported_checks_passed':True}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-checks',action='store_true')
    args=parser.parse_args()
    report=verify(pending_report=args.write_checks)
    if args.write_checks:
        REPORT.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,ensure_ascii=False,indent=2))
