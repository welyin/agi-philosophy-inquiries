"""Prepare695 entry checks using the existing non-round publication protocol."""
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent
src=(HERE.parent/'round694_drafts/check_and_publish_entry.py').read_text('utf8')
mapping={'694':'695','693':'694','critical_source_entry':'self_complement_entry'}
pattern='|'.join(r'(?<!\d)'+k+r'(?!\d)' if k.isdecimal() else re.escape(k) for k in mapping)
src=re.sub(pattern,lambda m:mapping[m.group()],src)
a=src.index('    entry=');b=src.index('    head,rest=',a)
src=src[:a]+'''    entry=(marker+' [原自旋互补与辅助相位]('+prefix+'round695_drafts/self_complement_entry.md)'
        '复用674互补Pfaffian，同一原两传播方向模型的辅助权重只剩实符号抵消问题。'
        '原矩阵相位检验已复算，全辅助非负与完整积分仍待证；'
        '正式科学轮次仍为694，旧空间接口及目标不变。')
'''+src[b:]
src=src.replace("'check_and_publish_entry.py')","'check_and_publish_entry.py','primary_source_audit.json')")
with (HERE/'check_and_publish_entry.py').open('x',encoding='utf8',newline='\n') as f:f.write(src)
print('695 entry publication check prepared')
