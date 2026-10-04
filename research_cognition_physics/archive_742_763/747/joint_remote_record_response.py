"""747: certified sterile-hop contribution, not the total remote jet."""
import argparse,hashlib,json,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'round747_drafts'))
import remote_population_entry as entry
import remote_hop_sign_certificate as sign
TARGET=HERE/'joint_remote_record_response_results.json'
def run():
    a=entry.run();assert a==json.loads(entry.TARGET.read_text('utf8'))
    b=sign.run();assert b==json.loads(sign.TARGET.read_text('utf8'))
    deps=('research_note_577.md','research_note_598.md','research_note_746.md',
        'round745_drafts/exact_history_density.py',
        'round746_drafts/higgs_readout_certificate.py',
        'round747_drafts/remote_population_entry.py',
        'round747_drafts/remote_hop_coefficient.py',
        'round747_drafts/remote_hop_sign_certificate.py')
    return dict(round=747,tests_run=2,failures=0,errors=0,
        checks=['original_full_CAR_remote_occupation_entry',
                'exact_remote_hop_squared_jet_and_certified_native_readout_integral'],
        entry=a,certificate=b,
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope='The sixth time-derivative coefficient with two original sterile hopping, two mass and two kinetic factors is exactly extracted and its Gaussian native T-readout integral is strictly negative. The original scalar-edge contribution remains unevaluated. No nonzero total remote readout, usable time window, ideal occupation measurement, autonomous detector or gravity closure is claimed. The hopping marker is bookkeeping, not a changed physical coupling.')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true')
    args=p.parse_args();r=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8') as f:json.dump(r,f,ensure_ascii=False,indent=2)
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(round=747,tests_run=2,total_remote_signal_still_open=True)))
