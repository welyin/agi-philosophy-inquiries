"""Prepare747 entry publication using the existing conflict-checked workflow."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
text=(ROOT/'round744_drafts/check_and_publish_entry.py').read_text('utf8')
mapping={'744':'747','743':'746','3424':'3430','3401':'3449',
         'native_instrument_selection_entry':'remote_population_entry'}
text=re.sub('|'.join(re.escape(k) for k in sorted(mapping,key=len,reverse=True)),lambda m:mapping[m[0]],text)
text=text.replace("checks['display_formulas']==9","checks['display_formulas']==4")
summary="summary=heading+' [原占据转移与实际异地读口]({p}round747_drafts/research_note_747_working.md)原完整H首步使邻端sterile占据在二阶出现严格输入差，原质量背景同时保留；这尚非原T读口的实际记录。下一项核真实异地读口及两条传播通道的共同时间次序，不插入理想中间占据测量。'"
text=re.sub(r'^summary=heading\+.*$',lambda m:summary,text,flags=re.M)
text=text.replace('native-instrument selection evidence','native remote-population entry evidence')
with (HERE/'check_and_publish_entry.py').open('x',encoding='utf8') as f:f.write(text)
print('Prepared747 evidence and navigation check.')
