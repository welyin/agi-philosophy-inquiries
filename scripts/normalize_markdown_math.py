"""Reversible delimiter-only Markdown maintenance. No source backups are created.
Plan is read-only; apply creates one manifest then changes only planned Markdown.
Verify checks exact inverses, protected regions, formula bodies and idempotence.
"""
from __future__ import annotations
import argparse, hashlib, json, os, re, subprocess, sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'research_cognition_physics/_migration/markdown_math_20261010'
MANIFEST=BASE/'manifest.json'
ACTIVE='research_cognition_physics/archive_1086_/'
EXCLUDE_COMPONENTS={'.git','.research_runtime','node_modules','site-packages','vendor','third_party'}
EXCLUDE_PREFIXES=('research_cognition_physics/_history/navigation_snapshots/',
                  'research_cognition_physics/_migration/markdown_math_20261010/')
TOKEN=re.compile(r'(?<!\\)(?:\\\\)*(\\[\[\]\(\)])')

def sha(b): return hashlib.sha256(b).hexdigest()
def decode(b):
    enc='utf-8-sig' if b.startswith(b'\xef\xbb\xbf') else 'utf-8'
    return b.decode(enc),enc

def file_list():
    p=subprocess.run(['rg','--files','--hidden','--no-ignore','-g','!**/.git/**','-g','!**/.research_runtime/**'],cwd=ROOT,capture_output=True,text=True,encoding='utf-8',check=True)
    return sorted(set(x.replace('\\','/') for x in p.stdout.splitlines()),key=str.casefold)

def excluded(rel):
    return bool(set(rel.split('/')) & EXCLUDE_COMPONENTS) or rel.startswith(EXCLUDE_PREFIXES)

def protected(text):
    """Conservative Markdown code/comment mask; preserves source offsets/newlines."""
    spans=[]; lines=text.splitlines(keepends=True); offset=0; fence=None; start=None
    indented=None; previous_blank=True
    for line in lines:
        clean=re.sub(r'^(?:[ \t]{0,3}> ?)+','',line)
        if fence:
            if re.match(r'^[ \t]*'+re.escape(fence[0])+'{'+str(fence[1])+r',}[ \t]*(?:\r?\n)?$',clean):
                spans.append((start,offset+len(line),'fenced_code')); fence=None
            offset+=len(line); previous_blank=not line.strip(); continue
        fm=re.match(r'^[ \t]*(?:(?:[-+*]|\d+[.)])[ \t]+)?(`{3,}|~{3,})([^\r\n]*)',clean)
        if fm and not(fm[1][0]=='`' and '`' in fm[2]):
            if indented is not None:
                spans.append((indented,offset,'indented_code')); indented=None
            fence=(fm[1][0],len(fm[1])); start=offset
        else:
            is_indent=bool(re.match(r'^(?: {4}|\t)',line))
            if indented is not None and line.strip() and not is_indent:
                spans.append((indented,offset,'indented_code')); indented=None
            if indented is None and previous_blank and is_indent and line.strip():
                indented=offset
        offset+=len(line); previous_blank=not line.strip()
    if fence: raise ValueError('Unclosed fenced code block')
    if indented is not None: spans.append((indented,len(text),'indented_code'))
    for m in re.finditer(r'<!--[\s\S]*?(?:-->|\Z)',text):
        if not m[0].endswith('-->'): raise ValueError('Unclosed HTML comment')
        spans.append((m.start(),m.end(),'html_comment'))
    mask=bytearray(len(text))
    for a,b,_ in spans: mask[a:b]=b'1'*(b-a)
    ticks=list(re.finditer(r'(?<!\\)(?:\\\\)*(`+)',text)); i=0
    while i<len(ticks):
        first=ticks[i]; a=first.start(1)
        if mask[a]: i+=1; continue
        n=len(first[1]); j=i+1
        while j<len(ticks):
            last=ticks[j]; z=last.start(1)
            if not mask[z] and len(last[1])==n:
                b=last.end(1); spans.append((a,b,'inline_or_multiline_code'))
                mask[a:b]=b'1'*(b-a); i=j+1; break
            j+=1
        else: i+=1
    return mask,sorted(spans)

def active_pairs(text,mask):
    stack=[]; pairs=[]
    for m in TOKEN.finditer(text):
        a=m.start(1); token=m[1]
        if mask[a]: continue
        if token in ('\\[','\\('):
            if stack: raise ValueError('Nested old math delimiters at '+str(a))
            stack.append((a,token))
        else:
            if not stack: raise ValueError('Unmatched closing '+token+' at '+str(a))
            a0,opening=stack.pop()
            if (opening,token) not in (('\\[','\\]'),('\\(','\\)')):
                raise ValueError('Crossed math delimiters')
            pairs.append((a0,a,'display' if opening=='\\[' else 'inline'))
    if stack: raise ValueError('Unmatched math opening '+str(stack))
    return pairs

def dollar_pairs(text,mask):
    tokens=[]; mentions=[]
    for m in re.finditer(r'(?<!\\)(?:\\\\)*(\$\$)',text):
        a=m.start(1)
        if mask[a]: continue
        # This exact plain-language notation mention is not a formula delimiter.
        if text[max(0,a-3):a]=='独占行':
            mentions.append({'line':text.count('\n',0,a)+1,'kind':'literal_notation_mention'}); continue
        tokens.append(a)
    if len(tokens)%2: raise ValueError('Unpaired existing dollar display delimiter')
    return [(a,b,'existing_display') for a,b in zip(tokens[::2],tokens[1::2])],mentions

def plan_text(text):
    mask,regions=protected(text); pairs=active_pairs(text,mask)
    dollars,_=dollar_pairs(text,mask)
    pairs=sorted(pairs+dollars)

    edits=[]; formulas=[]; nl='\r\n' if '\r\n' in text else '\n'
    for a,b,kind in pairs:
        body=text[a+2:b]
        if kind=='inline' and re.search(r'(?<!\\)(?:\\\\)*\$',body):
            raise ValueError('Unescaped dollar within inline formula at '+str(a))
        formula={'kind':kind,'body_sha256':sha(body.encode('utf8')),'body_characters':len(body)}
        if kind=='inline':
            formulas.append(formula)
            edits.extend([{'start':a,'end':a+2,'before':'\\(','after':'$','kind':'inline_open'},
                          {'start':b,'end':b+2,'before':'\\)','after':'$','kind':'inline_close'}])
            continue
        # Only delimiter replacements and boundary newline insertions are recorded.
        left=text[:a]; right=text[b+2:]
        prefix=''
        if left:
            if not left.endswith(('\n','\r')):
                prefix=nl+nl
            elif not re.search(r'(?:\r?\n)[ \t]*(?:\r?\n)$',left): prefix=nl
        suffix=''
        if right:
            if not right.startswith(('\n','\r')): suffix=nl+nl
            elif not re.match(r'^\r?\n[ \t]*\r?\n',right): suffix=nl
        inside_left='' if body.startswith(('\n','\r')) else nl
        inside_right='' if body.endswith(('\n','\r')) else nl
        proposed=[{'start':a,'end':a+2,'before':text[a:a+2],'after':prefix+'$$'+inside_left,'kind':kind+'_open'},
                  {'start':b,'end':b+2,'before':text[b:b+2],'after':inside_right+'$$'+suffix,'kind':kind+'_close'}]
        changed_edits=[e for e in proposed if e['before']!=e['after']]
        if changed_edits:
            formulas.append(formula)
            edits.extend(changed_edits)
    edits.sort(key=lambda e:e['start'])
    current=apply_edits(text,edits)
    # Body/nonmath bytes are copied directly; disjoint edits only replace two-character delimiters.
    for e in edits:
        assert e['end']-e['start']==2 and text[e['start']:e['end']]==e['before']
        assert e['after'].strip() in ('$','$$')
        assert not any(mask[e['start']:e['end']])
    assert inverse_text(current,edits)==text
    if edits:
        again_mask,_=protected(current)
        assert not active_pairs(current,again_mask), 'Conversion is not idempotent'
    return current,edits,formulas,regions

def apply_edits(text,edits):
    pieces=[]; prior=0
    for e in edits:
        assert prior<=e['start'] and text[e['start']:e['end']]==e['before']
        pieces.extend((text[prior:e['start']],e['after'])); prior=e['end']
    pieces.append(text[prior:]); return ''.join(pieces)

def inverse_text(text,edits):
    pieces=[]; offset=prior=0
    for e in edits:
        a=e['start']+offset; b=a+len(e['after'])
        if text[a:b]!=e['after']: raise ValueError('Formula reverse edit mismatch')
        pieces.extend((text[prior:a],e['before'])); prior=b
        offset+=len(e['after'])-len(e['before'])
    pieces.append(text[prior:]); return ''.join(pieces)

def restore_bytes(raw,entry):
    if sha(raw)!=entry['current_sha256'] or len(raw)!=entry['current_bytes']:
        raise ValueError('Formatted file changed: '+entry['path'])
    text,enc=decode(raw); prior=inverse_text(text,entry['edits']).encode(enc)
    if sha(prior)!=entry['original_sha256'] or len(prior)!=entry['original_bytes']:
        raise ValueError('Pre-format bytes not recovered: '+entry['path'])
    return prior

def scan():
    names=file_list(); entries=[]; skipped=[]; errors=[]; inventory={}; assets={}; notation_mentions=[]
    for rel in names:
        if excluded(rel):
            if rel.lower().endswith('.md'): skipped.append(rel)
            continue
        path=ROOT/rel
        if path.suffix.lower() in ('.py','.json') and not rel.startswith(ACTIVE):
            raw=path.read_bytes(); assets[rel]={'sha256':sha(raw),'bytes':len(raw)}
        if path.suffix.lower()!='.md': continue
        raw=path.read_bytes(); inventory[rel]={'sha256':sha(raw),'bytes':len(raw)}
        try:
            text,enc=decode(raw); changed,edits,formulas,regions=plan_text(text)
            pmask,_=protected(text)
            _,mentions=dollar_pairs(text,pmask)
            notation_mentions.extend({'path':rel,**item} for item in mentions)
        except Exception as ex:
            errors.append({'path':rel,'error':str(ex)}); continue
        if not edits: continue
        after=changed.encode(enc)
        entries.append({'path':rel,'original_sha256':sha(raw),'current_sha256':sha(after),
                        'original_bytes':len(raw),'current_bytes':len(after),'encoding':enc,
                        'edits':edits,'formulas':formulas,
                        'protected_regions':len(regions)})
    return {'schema':'markdown_math_format_v1','scope':'project Markdown except explicit snapshots/runtime/external trees',
            'root':str(ROOT),'entries':entries,'markdown_inventory':inventory,'old_python_json_inventory':assets,
            'excluded_markdown':skipped,'errors':errors,'literal_notation_mentions':notation_mentions,
            'summary':{'markdown_scanned':len(inventory),'files_changed':len(entries),
                       'display_formulas':sum(f['kind']=='display' for e in entries for f in e['formulas']),
                       'inline_formulas':sum(f['kind']=='inline' for e in entries for f in e['formulas']),
                       'existing_dollar_blocks_boundary_adjusted':sum(f['kind']=='existing_display' for e in entries for f in e['formulas']),
                       'old_python_json_files':len(assets),'errors':len(errors)},
            'protections':['fenced and indented code','inline and multiline code','HTML comments','even backslash runs including TeX line-break spacing','bare brackets unchanged'],
            'limitations':['Conservative scanner, not full CommonMark AST; protected ambiguities are skipped rather than rewritten.',
                           'Existing dollar display boundaries are also normalized; their mathematical bodies remain byte-identical.',
                           'Exact reversal and body hashes certify format-only edits, not visual rendering in every Markdown engine.'],
            'tool_sha256':sha(Path(__file__).read_bytes())}

def verify(data):
    bypath={e['path']:e for e in data['entries']}; total_bodies=0
    for rel,old in data['markdown_inventory'].items():
        raw=(ROOT/rel).read_bytes()
        if rel in bypath:
            e=bypath[rel]; prior=restore_bytes(raw,e)
            text,enc=decode(prior); after,edits,formulas,_=plan_text(text)
            assert after.encode(enc)==raw and edits==e['edits'] and formulas==e['formulas'],rel
            current,_=decode(raw); again,redo,_,_=plan_text(current)
            assert again==current and not redo,rel
            total_bodies+=len(formulas)
        elif sha(raw)!=old['sha256'] or len(raw)!=old['bytes']:
            raise ValueError('Noncandidate Markdown changed: '+rel)
    for rel,old in data['old_python_json_inventory'].items():
        raw=(ROOT/rel).read_bytes()
        if sha(raw)!=old['sha256'] or len(raw)!=old['bytes']:
            raise ValueError('Existing Python/JSON changed: '+rel)
    return {'status':'pass','files_changed':len(bypath),'formula_bodies_verified':total_bodies,
            'all_old_python_json_byte_identical':len(data['old_python_json_inventory']),
            'all_scanned_markdown_accounted':len(data['markdown_inventory']),
            'inverse_hashes':True,'idempotent':True,'nonmath_and_code_preserved':True}

def selftest():
    cases=[r'alpha \(x+1\) omega',r'\[x+y\]',
           'before\n\\[\nx\\\\[6pt]y\n\\]\nafter',
           '```python\n\\[notmath\\]\n```\n\\(z\\)',
           '`one \\[ two \\]`\n\\(z\\)',
           '`multiline\n\\[literal\\]\ncode`\n\\(z\\)',
           '<!-- \\[comment\\] -->\n\\(z\\)',
           '    \\[indented\\]\n\n\\(z\\)',
           '> ~~~\n> \\[code\\]\n> ~~~\n\\(z\\)',
           '$$\nx\\\\[6pt]y\n$$\n\\(z\\)', '$$x+y$$', 'before\n$$\nx\n$$\nafter']
    for text in cases:
        out,edits,_,_=plan_text(text)
        assert inverse_text(out,edits)==text
        assert not plan_text(out)[1]
    assert len(plan_text(cases[5])[2])==1
    assert len(plan_text(cases[2])[2])==1
    return {'cases':len(cases),'pass':True}

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('mode',choices=['plan','apply','verify','selftest'])
    args=ap.parse_args()
    if args.mode=='selftest': print(json.dumps(selftest())); return
    if args.mode=='verify': print(json.dumps(verify(json.loads(MANIFEST.read_text('utf8'))),ensure_ascii=False)); return
    data=scan()
    if args.mode=='plan':
        print(json.dumps({'summary':data['summary'],'errors':data['errors'],'excluded_markdown':data['excluded_markdown']},ensure_ascii=False,indent=2)); return
    if data['errors']: raise ValueError(json.dumps(data['errors'],ensure_ascii=False))
    if MANIFEST.exists(): raise FileExistsError('Format layer exists; use verify, never overwrite')
    selftest(); BASE.mkdir(parents=True,exist_ok=True)
    # Verify all old inputs again before committing the single reversible manifest.
    for e in data['entries']:
        if sha((ROOT/e['path']).read_bytes())!=e['original_sha256']: raise ValueError('Concurrent Markdown change: '+e['path'])
    MANIFEST.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    for e in data['entries']:
        path=ROOT/e['path']; raw=path.read_bytes()
        if sha(raw)!=e['original_sha256']: raise ValueError('Concurrent edit stopped apply: '+e['path'])
        text,enc=decode(raw); after=apply_edits(text,e['edits']).encode(enc)
        assert sha(after)==e['current_sha256']
        path.write_bytes(after)
    receipt=verify(data)
    (BASE/'format_verification.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    print(json.dumps({'summary':data['summary'],'verification':receipt,'manifest_sha256':sha(MANIFEST.read_bytes())},ensure_ascii=False))

if __name__=='__main__': main()
