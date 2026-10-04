"""Prepare697 entry check, preserving696 formal count and frozen evidence."""
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent
src=(HERE.parent/'round696_drafts/check_and_publish_entry.py').read_text('utf8')
mapping={'696':'697','695':'696','exterior_average_entry':'center_holonomy_entry'}
pattern='|'.join(r'(?<!\d)'+k+r'(?!\d)' if k.isdecimal() else re.escape(k) for k in mapping)
src=re.sub(pattern,lambda m:mapping[m.group()],src)
src=src.replace('local exterior-average screen','original centre-holonomy sewing')
a=src.index('    entry=');b=src.index('    head,rest=',a)
src=src[:a]+'''    entry=(marker+' [中心半区的全群闭合与四角积分]('+prefix+'round697_drafts/center_holonomy_entry.md)'
        '原256维来源的同步输送已复算，完整辅助平均后的最终holonomy可用全原群Weyl积分处理。'
        '剩余积分尚未求值，不能由固定接缝负核判定物理RP；正式科学轮次仍为696，目标不变。')
'''+src[b:]
src=src.replace("'check_and_publish_entry.py')","'check_and_publish_entry.py','primary_source_audit.json')")
with (HERE/'check_and_publish_entry.py').open('x',encoding='utf8',newline='\n') as f:f.write(src)
print('697 entry publication check prepared')
