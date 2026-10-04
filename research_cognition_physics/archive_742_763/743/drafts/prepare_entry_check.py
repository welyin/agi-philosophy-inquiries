"""Prepare743 entry publication, preserving completed742."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent
s=(HERE.parent/'round742_drafts/check_and_publish_entry.py').read_text('utf8')
mapping={'742':'743','741':'742','3421':'3422','3373':'3387',
         'gaussian_record_entry':'quantum_scalar_record_entry'}
s=re.sub('|'.join(re.escape(k) for k in sorted(mapping,key=len,reverse=True)),lambda m:mapping[m.group()],s)
s=s.replace("checks['display_formulas']==12","checks['display_formulas']==10")
start=s.index('summary=heading+');end=s.index('\nsummary=summary.replace',start)
s=s[:start]+"summary=heading+' [原量子标量与实际读口]( {p}round743_drafts/research_note_743_working.md)原Yukawa／Majorana的量子系数可在短时产生非高斯关联，解除742限制的资源已在原作用内；同一涨落也扰动占据数，不能直接认成理想记录。原128维对象核验通过，下一项接实际仪器、Gauss准备与共同来源。'"+s[end:]
s=s.replace('Gaussian-implementation evidence','native-scalar quantum-jet evidence')
with (HERE/'check_and_publish_entry.py').open('x',encoding='utf8',newline='\n') as f:f.write(s)
