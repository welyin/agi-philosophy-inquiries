"""Preserve and correct one escaped-backspace typo before743 publication."""
from pathlib import Path
import hashlib
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
assert not (ROOT/'research_round_743_checks.json').exists()
p=ROOT/'research_note_743.md';raw=p.read_bytes()
with (HERE/'research_note_743_before_text_repair.md').open('xb') as f:f.write(raw)
text=raw.decode('utf8');assert text.count(chr(8))==1
text=text.replace(chr(8)+'ar',chr(92)+'bar')
p.write_text(text,encoding='utf8',newline='\n')
(HERE/'research_note_743_draft.md').write_text(text,encoding='utf8',newline='\n')
review=HERE/'final_review.txt'
with (HERE/'final_review_before_text_repair.txt').open('xb') as f:f.write(review.read_bytes())
names=('research_note_743.md','joint_native_quantum_record.py','joint_native_quantum_record_results.json','unified_physics_condition_ledger_743.md')
head=review.read_text('utf8').split('\n',1)[0]
review.write_text(head+'\n'+'\n'.join(n+' '+hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names)+'\n',encoding='utf8',newline='\n')
print('Preserved initial743 draft; repaired one formula text escape.')
