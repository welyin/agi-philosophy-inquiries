"""Register703 bounded-source entry, not a completed science round."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent
s=(HERE.parent/'round702_drafts/check_and_publish_entry.py').read_text('utf8')
s=re.sub(r'(?<!\d)(701|702)(?!\d)',lambda m:str(int(m.group())+1),s)
s=s.replace('source_scale_entry','bounded_source_entry')
s=s.replace('normalized source scale corollary of699 negative certificate',
            'bounded original-record source derivatives and finite time-cut identity')
s=s.replace("['display_formulas']==2","['display_formulas']==6")
a=s.index('    entry=');b=s.index('    head,rest=',a)
s=s[:a]+'''    entry=(marker+' [原记录来源的热响应与时间切口]('+prefix+'round703_drafts/bounded_source_entry.md)'
        '702共同热极限可扩展到原有界读口的来源导数；'
        '有限步配分函数得分与原测量切口期望不能直接混同，二者差随同一极限消失。'
        '解析控制及原耦合诊断已核；一般几何来源仍待接。正式科学轮次保持702，目标不变。')
'''+s[b:]
with (HERE/'check_and_publish_entry.py').open('x',encoding='utf8',newline='\n') as f:f.write(s)
print('703 entry publication prepared')
