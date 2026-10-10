"""Read-only replay: formula inverse -> SR navigation inverse -> frozen spatial/older layers.
Existing runtime source, archived Python/JSON and their hash receipts are untouched.
"""
from __future__ import annotations
import argparse, hashlib, json, os, sys
from pathlib import Path
import run_research_spatial_archive as spatial
import normalize_markdown_math as formatting
from research_layout import ROOT, relocation_original_bytes, no_archive_or_research_writes, Layout
THIS=Path(__file__).resolve()
NAV=ROOT/'research_cognition_physics/archive_1086_/_shared/navigation_layer.json'

class FormattedLayout(spatial.SpatialLayout):
    def __init__(self):
        self.format_layer=json.loads(formatting.MANIFEST.read_text('utf8'))
        if self.format_layer['schema']!='markdown_math_format_v1': raise ValueError('Unknown formula layer')
        self.formats={}
        for e in self.format_layer['entries']:
            path=ROOT/e['path']
            if not path.resolve().is_relative_to(ROOT) or path.suffix.lower()!='.md':
                raise ValueError('Formula layer must affect only project Markdown')
            self.formats[os.path.normcase(os.path.abspath(path))]=e
        self.sr_navigation=json.loads(NAV.read_text('utf8'))
        if self.sr_navigation['schema']!='research_sr_navigation_layer_v1': raise ValueError('Unknown navigation layer')
        self.sr={}
        for e in self.sr_navigation['entries']:
            if e['original']!=e['destination'] or e['original'] not in spatial.live_runtime.NAVIGATION:
                raise ValueError('Unexpected SR navigation file')
            self.sr[os.path.normcase(os.path.abspath(ROOT/e['destination']))]=e
        if len(self.sr)!=4: raise ValueError('Expected four SR navigation files')
        super().__init__()

    def physical_bytes(self,path):
        actual=Path(str(path)); raw=actual.read_bytes()
        key=os.path.normcase(os.path.abspath(actual))
        if key in self.formats: raw=formatting.restore_bytes(raw,self.formats[key])
        if key in self.sr: raw=relocation_original_bytes(raw,self.sr[key])
        return raw

    def verify_layout(self):
        formatted=formatting.verify(self.format_layer)
        prior=super().verify_layout()
        # Every older source-entry inverse is checked through the inherited read_pre_split chain.
        old=Layout.verify(self)
        return {'format_layer':formatted,'sr_navigation':{'files':len(self.sr),'inverse_hashes_verified':True},
                'spatial_and_prior_layout':prior,'older_source_layout':old,'all_checks_passed':True}

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    g=ap.add_mutually_exclusive_group(required=True)
    g.add_argument('--script'); g.add_argument('--verify-layout',action='store_true')
    ap.add_argument('arguments',nargs=argparse.REMAINDER)
    args=ap.parse_args()
    os.environ['PYTHONDONTWRITEBYTECODE']='1'; sys.dont_write_bytecode=True
    store=FormattedLayout(); runtime=spatial.SpatialRuntime(store)
    sys.path.append(str(ROOT/'.research_runtime'))
    # The inherited subprocess router must start this outer reader in child checks as well.
    spatial.THIS=THIS
    sys.addaudithook(no_archive_or_research_writes)
    with runtime.installed(),spatial.execution_routes(store,runtime):
        if args.verify_layout: print(json.dumps(store.verify_layout(),ensure_ascii=False))
        else:
            arguments=args.arguments[1:] if args.arguments[:1]==['--'] else args.arguments
            runtime.run_script(args.script,arguments)

if __name__=='__main__': main()
