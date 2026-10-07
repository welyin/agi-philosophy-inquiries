"""Prepare publication of the scope-only 775 working report."""
from pathlib import Path
import re

HERE = Path(__file__).resolve().parent
mapping = {'773': '775', '772': '774', '3502': '3508', '3789': '3815',
           'local_insertion_entry': 'causal_source_entry'}
rx = re.compile('|'.join(map(re.escape, sorted(mapping, key=len, reverse=True))))
convert = lambda text: rx.sub(lambda m: mapping[m[0]], text)
source = convert((HERE/'publish_round773_entry.py').read_text('utf8'))
summary = '**775共同来源入口已推进，正式仍774／3508：** [工作报告]({p}round775_drafts/research_note_775_working.md)保772原态，以代数期望检查闭时路径入口；不另要求全局迹类密度矩阵。初态最低边界次数与同一反常/来源映射尚待核。[条件增量]({p}round775_drafts/entry_condition_ledger.md)、[入口核验]({p}round775_drafts/causal_source_entry_checks.json)。本项不增加科学轮次，目标保持。'
start, stop = source.index('summary = '), source.index('planned = ')
source = source[:start]+'summary = '+repr(summary)+'\n'+source[stop:]
source = source.replace('**775共同插入入口已核', '**775共同来源入口已推进')
post = convert((HERE/'postcheck_round773_entry.py').read_text('utf8'))
post = post.replace('**775共同插入入口已核', '**775共同来源入口已推进')
for name, value in [('publish_round775_entry.py', source), ('postcheck_round775_entry.py', post)]:
    with (HERE/name).open('x', encoding='utf8', newline='\n') as stream:
        stream.write(value)
print('Prepared 775 working-entry publication; scientific counters unchanged.')
