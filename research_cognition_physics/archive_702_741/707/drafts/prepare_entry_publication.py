"""Register707 spatial map audit without a new formal scientific round."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent
s=(HERE.parent/'round706_drafts/check_and_publish_entry.py').read_text('utf8')
s=re.sub(r'(?<!\d)(705|706)(?!\d)',lambda m:str(int(m.group())+1),s)
s=s.replace('local_domain_entry','spatial_map_entry')
s=s.replace('graded local Stinespring and two-sided comparison-energy domain',
            'spatial cut/refinement contract and original character flux allocation')
s=s.replace("['display_formulas']==6","['display_formulas']==4")
s=s.replace(",'parity_completion_706.md'","")
a=s.index('    entry=');b=s.index('    head,rest=',a)
s=s[:a]+'''    entry=(marker+' [区域切口与真实空间细化的系数合同]('+prefix+'round707_drafts/spatial_map_entry.md)'
        '复用617／637，核原各向同性串联细化的电能因子4；'
        '原字符四路径通量比较另有整数分配条件。'
        '仅为指定电部门的入口审计，不当完整物理极限或新整轮；'
        '正式科学轮次保持706，旧空间接口及目标不变。')
'''+s[b:]
with (HERE/'check_and_publish_entry.py').open('x',encoding='utf8',newline='\n') as f:f.write(s)
print('707 entry publication prepared')
