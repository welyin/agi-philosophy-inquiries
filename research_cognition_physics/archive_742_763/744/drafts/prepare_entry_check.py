"""Prepare744 entry publication, preserving completed743."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent
s=(HERE.parent/'round743_drafts/check_and_publish_entry.py').read_text('utf8')
mapping={'743':'744','742':'743','3422':'3424','3387':'3401',
         'quantum_scalar_record_entry':'native_instrument_selection_entry'}
s=re.sub('|'.join(re.escape(k) for k in sorted(mapping,key=len,reverse=True)),lambda m:mapping[m.group()],s)
s=s.replace("checks['display_formulas']==10","checks['display_formulas']==9")
start=s.index('summary=heading+');end=s.index('\nsummary=summary.replace',start)
s=s[:start]+"summary=heading+' [原未知态读口与精确选择规则]( {p}round744_drafts/research_note_744_working.md)原全H的联合反射／费米相位使对称准备下的单次sin s效果精确只读编码相干，任意等待均不分空态与配对态。实际CP输出已列出；下一项核多次读取的联合历史与来源，不扩大单次限制。'"+s[end:]
s=s.replace('native-scalar quantum-jet evidence','native-instrument selection evidence')
with (HERE/'check_and_publish_entry.py').open('x',encoding='utf8',newline='\n') as f:f.write(s)
