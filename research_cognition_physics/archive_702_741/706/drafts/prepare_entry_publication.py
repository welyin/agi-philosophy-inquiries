"""Register706 domain entry and705 graded completion without a new formal round."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent
s=(HERE.parent/'round705_drafts/check_and_publish_entry.py').read_text('utf8')
s=re.sub(r'(?<!\d)(704|705)(?!\d)',lambda m:str(int(m.group())+1),s)
s=s.replace('regional_compression_entry','local_domain_entry')
s=s.replace('regional compression algebra and inherited thermal tail corollary',
            'graded local Stinespring and two-sided comparison-energy domain')
s=s.replace("saved['dependency_hashes']","saved['dependencies']")
s=s.replace("'check_and_publish_entry.py','prepare_entry_publication.py')",
            "'check_and_publish_entry.py','prepare_entry_publication.py','parity_completion_705.md')")
a=s.index('    entry=');b=s.index('    head,rest=',a)
s=s[:a]+'''    entry=(marker+' [保边界通道的双向能量域]('+prefix+'round706_drafts/local_domain_entry.md)'
        '以固定辅助能量标签接通705局部通道的正向及伴随图范数，623原共同域直接复用。'
        '[705的CAR宇称分块补充]('+prefix+'round706_drafts/parity_completion_705.md)'
        '已核验，旧冻结文件保留；完整真实来源历史下一步检验。'
        '正式科学轮次保持705，目标不变。')
'''+s[b:]
with (HERE/'check_and_publish_entry.py').open('x',encoding='utf8',newline='\n') as f:f.write(s)
print('706 entry publication prepared')

