from pathlib import Path

root = Path(__file__).parent.parent

html = str(root / 'templates' / 'html.tpl')
csv = str(root / 'templates' / 'csv.tpl')
tex = str(root / 'templates' / 'tex.tpl')
