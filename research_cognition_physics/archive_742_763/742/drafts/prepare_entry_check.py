"""Prepare742 entry publication; no prior evidence edits."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent
s=(HERE.parent/'round741_drafts/check_and_publish_entry.py').read_text('utf8')
mapping={'741':'742','740':'741','3419':'3421','3359':'3373',
         'absolute_source_order_entry':'gaussian_record_entry'}
s=re.sub('|'.join(re.escape(k) for k in sorted(mapping,key=len,reverse=True)),lambda m:mapping[m.group()],s)
s=s.replace("checks['display_formulas']==7","checks['display_formulas']==12")
start=s.index('summary=heading+');end=s.index('\nsummary=summary.replace',start)
s=s[:start]+"summary=heading+' [原记录与二次实现限制]({p}round742_drafts/research_note_742_working.md)原实际纯参考的非选择记录有非零四点Wick缺陷，确定二次演化加独立高斯探针无法任意逼近，给出正迹距离下界。仅限该实现类，741双线性来源保持；下一项核原量子玻色—费米相互作用。'"+s[end:]
s=s.replace('source-order evidence','Gaussian-implementation evidence')
with (HERE/'check_and_publish_entry.py').open('x',encoding='utf8',newline='\n') as f:f.write(s)
