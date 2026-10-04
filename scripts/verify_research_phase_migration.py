"""Legacy command name for the current ZIP-free layout verification."""
from datetime import datetime
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
arguments = sys.argv[1:]
if not arguments:
    folder = ROOT / '.research_runtime'
    folder.mkdir(exist_ok=True)
    arguments = ['--report', str(folder / ('migration_check_' + datetime.now().strftime('%Y%m%d_%H%M%S') + '.json'))]
raise SystemExit(subprocess.call([sys.executable, '-B', '-X', 'utf8',
    str(ROOT / 'scripts/verify_research_final_layout.py'), *arguments], cwd=ROOT))
