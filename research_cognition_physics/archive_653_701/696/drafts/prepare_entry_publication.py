"""Prepare696 entry check; retain695 formal count and all frozen evidence."""
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent
src=(HERE.parent/'round695_drafts/check_and_publish_entry.py').read_text('utf8')
mapping={'695':'696','694':'695','self_complement_entry':'exterior_average_entry'}
pattern='|'.join(r'(?<!\d)'+k+r'(?!\d)' if k.isdecimal() else re.escape(k) for k in mapping)
src=re.sub(pattern,lambda m:mapping[m.group()],src)
src=src.replace('critical-source diagnostic','local exterior-average screen')
a=src.index('    entry=');b=src.index('    head,rest=',a)
src=src[:a]+'''    entry=(marker+' [原局部球面张量的外代数表示]('+prefix+'round696_drafts/exterior_average_entry.md)'
        '精确120维块具有正负2/5特征值，不能把这一直接重排当作正算符或投影。'
        '这只筛查表示，未计算原完整物理平均；正式科学轮次仍为695，'
        '旧空间接口逐项复用，目标不变。')
'''+src[b:]
src=src.replace(",'primary_source_audit.json'",'')
with (HERE/'check_and_publish_entry.py').open('x',encoding='utf8',newline='\n') as f:f.write(src)
print('696 entry publication check prepared')
