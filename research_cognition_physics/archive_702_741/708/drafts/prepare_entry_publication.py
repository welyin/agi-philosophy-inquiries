"""Register708 nested-amplitude entry without a new formal round."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent
s=(HERE.parent/'round707_drafts/check_and_publish_entry.py').read_text('utf8')
s=re.sub(r'(?<!\d)(706|707)(?!\d)',lambda m:str(int(m.group())+1),s)
s=s.replace('spatial_map_entry','nested_mass_entry')
s=s.replace('spatial cut/refinement contract and original character flux allocation',
            'nested hyperboloid amplitude and original collective CAR masses')
s=s.replace("['display_formulas']==4","['display_formulas']==6")
a=s.index('    entry=');b=s.index('    head,rest=',a)
s=s[:a]+'''    entry=(marker+' [嵌套中点与共同质量幅度]('+prefix+'round708_drafts/nested_mass_entry.md)'
        '原双曲向量和的同一幅度同时支撑嵌套组合与全部Dirac／Majorana均匀质量块。'
        '707径向偏迹不能保更强的幅度合同，保幅度的较大块压缩下一步检验；'
        '此为旧式推论及原128模式校准，不计整轮，正式科学轮次保持707，目标不变。')
'''+s[b:]
with (HERE/'check_and_publish_entry.py').open('x',encoding='utf8',newline='\n') as f:f.write(s)
print('708 entry publication prepared')
