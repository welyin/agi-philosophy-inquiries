"""Prepare694 entry checks from the existing publication protocol."""
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent
src=(HERE.parent/'round693_drafts/check_and_publish_entry.py').read_text('utf8')
mapping={'693':'694','692':'693','finite_character_entry':'critical_source_entry'}
pattern='|'.join(r'(?<!\d)'+k+r'(?!\d)' if k.isdecimal() else re.escape(k) for k in mapping)
src=re.sub(pattern,lambda m:mapping[m.group()],src)
src=src.replace('route screen','critical-source diagnostic')
a=src.index('    entry=');b=src.index('    head,rest=',a)
src=src[:a]+'''    entry=(marker+' [原通量路径的临界来源诊断]('+prefix+'round694_drafts/critical_source_entry.md)'
        '保留全部原通道和512维来源，观察到首次谱穿越附近的单侧非零数值与零分区；'
        '严格根包围及非零极限尚待认证，不宣称跳变定理或RP反例。'
        '正式科学轮次仍为693，旧空间接口直接复用，目标不变。')
'''+src[b:]
with (HERE/'check_and_publish_entry.py').open('x',encoding='utf8',newline='\n') as f:f.write(src)
print('694 entry publication check prepared')
