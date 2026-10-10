"""Read-only verification of round 1044; --record writes its first acceptance.

Numerical and rational checks verify the declared calculations, not the
semantics or all cognitive axioms. The independent prose review remains needed.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import runpy
import subprocess
import sys

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
ROOT = STAGE.parent
ASSETS = [
    STAGE / 'research_note_1044.md',
    HERE / 'scattering_instrument_proof.md',
    HERE / 'scattering_family.py',
    HERE / 'scattering_family_results.json',
    HERE / 'finite_window_proof.md',
    HERE / 'finite_window_bounds.py',
    HERE / 'finite_window_bounds_results.json',
    HERE / 'independent_review.md',
    HERE / 'dependency_delta_1044.json',
    Path(__file__).resolve(),
]
HISTORY = [
    ROOT / 'archive_1009_1043/research_note_1043.md',
    ROOT / 'archive_956_989/research_note_970.md',
    ROOT / 'archive_956_989/research_note_971.md',
    ROOT / 'archive_956_989/research_note_956.md',
    ROOT / 'archive_956_989/research_note_959.md',
    ROOT / 'archive_956_989/research_note_967.md',
    ROOT / 'archive_956_989/research_note_985.md',
    ROOT / 'archive_923_934/research_note_929.md',
    ROOT / 'archive_990_1008/research_note_1002.md',
]
RECEIPT = HERE / 'research_round_1044_checks.json'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def manifest(paths):
    return {p.relative_to(ROOT).as_posix():digest(p) for p in paths}


def compare(a,b,path='root'):
    if isinstance(a,dict):
        assert isinstance(b,dict) and a.keys()==b.keys(),path
        for k in a: compare(a[k],b[k],path+'.'+k)
    elif isinstance(a,list):
        assert isinstance(b,list) and len(a)==len(b),path
        for i,(x,y) in enumerate(zip(a,b)): compare(x,y,path+f'[{i}]')
    elif isinstance(a,float):
        assert isinstance(b,(int,float)) and math.isclose(a,b,rel_tol=1e-10,abs_tol=1e-12),(path,a,b)
    else:
        assert a==b,(path,a,b)


def run_check(script, result):
    env=os.environ.copy()
    env['OPENBLAS_NUM_THREADS']='1'
    cp=subprocess.run([sys.executable,'-B','-X','utf8',str(HERE/script)],
                      check=True,capture_output=True,text=True,encoding='utf-8',env=env)
    actual=json.loads(cp.stdout)
    if script=='scattering_family.py':
        assert actual['passed']
        os.environ['OPENBLAS_NUM_THREADS']='1'
        namespace=runpy.run_path(str(HERE/script))
        actual=namespace['run']()
        refinement=namespace['run'](384)
        actual['quadrature_refinement_probability_error']=max(
            abs(a['reflection_probability_data1']-b['reflection_probability_data1'])
            for a,b in zip(actual['rows'],refinement['rows']))
    compare(actual,json.loads((HERE/result).read_text(encoding='utf-8')))
    return actual


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--record',action='store_true')
    args=parser.parse_args()
    before=manifest(ASSETS+HISTORY)
    scattering=run_check('scattering_family.py','scattering_family_results.json')
    finite=run_check('finite_window_bounds.py','finite_window_bounds_results.json')
    assert finite['all_passed']
    assert finite['finite_probability_gap_lower']['decimal'] > .317
    assert finite['common_record_probability_lower']['decimal'] > .55
    assert scattering['g4_minus_g2_probability_gap'] > .32
    assert scattering['max_boundary_match_error'] < 1e-12
    links=0
    for doc in [p for p in ASSETS if p.suffix=='.md']:
        body=doc.read_text(encoding='utf-8')
        body=re.sub(chr(96)*3+r'.*?'+chr(96)*3,'',body,flags=re.S)
        body=re.sub(r'\$\$.*?\$\$|\\\[.*?\\\]','',body,flags=re.S)
        body=re.sub(r'(?<!\\)\$(?:\\.|[^$\n])*(?<!\\)\$','',body)
        for target in re.findall(r'(?<![_\w])\[[^\]\n]+\]\(([^)]+)\)',body):
            target=target.strip('<>').split('#')[0]
            if not target or re.match(r'^\w+://',target): continue
            resolved=(doc.parent/target).resolve()
            if not (args.record and resolved==RECEIPT):
                assert resolved.exists(),(doc,target)
            links+=1
    after=manifest(ASSETS+HISTORY)
    assert before==after,'read-only verification changed assets'
    report={
        'schema':'round1044_joint_scattering_acceptance_v1',
        'round':1044,'new_scientific_groups':1,
        'cumulative_before':3819,'cumulative_after':3820,
        'new_cognitive_axioms':0,
        'scientific_assets':manifest(ASSETS),
        'history_inputs':manifest(HISTORY),
        'read_only_recomputations_passed':True,
        'rational_checks_passed':len(finite['checks']),
        'local_links_checked':links,
        'finite_probability_gap_lower':finite['finite_probability_gap_lower'],
        'same_preparation_time_and_readout':True,
        'all_unknown_data_and_passive_reference':True,
        'physical_poststates_and_recoil_retained':True,
        'terminal_instrument_energy_closed':False,
        'all_A1_A7_or_GR_SM_realized':False,
        'scope':'Adopted bounded two-body effective interaction and fixed finite task; not full cognitive non-implication.',
        'independent_review_file':'independent_review.md',
        'semantic_completion_not_proved_by_this_script':True
    }
    if args.record:
        assert not RECEIPT.exists(),'Do not overwrite a frozen receipt'
        RECEIPT.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:
        assert RECEIPT.exists(),'First acceptance has not been recorded'
        compare(report,json.loads(RECEIPT.read_text(encoding='utf-8')))
    print(json.dumps({'all_passed':True,'round':1044,'rational_checks':len(finite['checks']),
                      'local_links':links,'mode':'record' if args.record else 'read_only',
                      'finite_gap_lower':finite['finite_probability_gap_lower']['decimal']},
                     ensure_ascii=False))


if __name__=='__main__':
    main()
