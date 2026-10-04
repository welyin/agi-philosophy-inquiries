"""Register704 own-Gibbs compression entry without a completed new round."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent
s=(HERE.parent/'round703_drafts/check_and_publish_entry.py').read_text('utf8')
s=re.sub(r'(?<!\d)(702|703)(?!\d)',lambda m:str(int(m.group())+1),s)
s=s.replace('bounded_source_entry','own_gibbs_compression_entry')
s=s.replace('bounded original-record source derivatives and finite time-cut identity',
            'own finite Gibbs initial-state corollary for original625 real histories')
s=s.replace("['display_formulas']==6","['display_formulas']==4")
a=s.index('    entry=');b=s.index('    head,rest=',a)
s=s[:a]+'''    entry=(marker+' [自身热初态与原真实记录响应]('+prefix+'round704_drafts/own_gibbs_compression_entry.md)'
        '复用625完整CP及来源定理，以A²热尾界证明有限过程自身基点Gibbs态也保真实记录二阶响应。'
        '有限初态不相同，替代族不冒充703对数转移；准备几何与动态来源混合项下一步检验。'
        '精确初态身份及实际历史已核，正式科学轮次保持703，目标不变。')
'''+s[b:]
with (HERE/'check_and_publish_entry.py').open('x',encoding='utf8',newline='\n') as f:f.write(s)
print('704 entry publication prepared')
