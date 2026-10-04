"""Preserve precheck drafts, then fix prose accidentally parsed as an MD link."""
from pathlib import Path
import hashlib
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
note=ROOT/'research_note_637.md';review=HERE/'final_review.txt'
for src,dest in ((note,HERE/'pre_text_fix_note.md'),(review,HERE/'pre_text_fix_review.txt')):
    with dest.open('xb') as f:f.write(src.read_bytes())
old=note.read_text('utf8')
new=old.replace('P_N=1_[1,N](A₀)，p_N=TrρβP_N',
                'P_N为A₀在区间[1,N]上的谱投影，p_N=TrρβP_N')
assert new!=old
oldhash=hashlib.sha256(note.read_bytes()).hexdigest()
for p in (note,HERE/'research_note_637_draft.md'):
    assert p.read_text('utf8')==old
    p.write_text(new,encoding='utf8',newline='\n')
r=review.read_text('utf8').replace(oldhash,hashlib.sha256(note.read_bytes()).hexdigest())
r+='Prepublication text check parsed the prose interval as a Markdown link.\n'
r+='Reworded only that interval; preserved precheck note and review as two extra evidence files.\n'
review.write_text(r,encoding='utf8',newline='\n')
v=ROOT/'verify_round637.py'
s=v.read_text('utf8').replace(
    "'round637_drafts/final_review.txt','round638_drafts/STATUS.md')",
    "'round637_drafts/final_review.txt','round638_drafts/STATUS.md',\n"
    "           'round637_drafts/pre_text_fix_note.md','round637_drafts/pre_text_fix_review.txt')")
s=s.replace('2057','2059')
v.write_text(s,encoding='utf8',newline='\n')
p=ROOT/'publish_round637.py'
p.write_text(p.read_text('utf8').replace('2057','2059'),encoding='utf8',newline='\n')
print('Preserved two precheck files; fixed interval text; final evidence count 2059')
