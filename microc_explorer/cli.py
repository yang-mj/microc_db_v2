"""Optional one-region export, suitable for shell batch loops."""
import argparse
from pathlib import Path
from .core import load_config
from .service import run

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config',default='config.json')
    p.add_argument('--source',default='merged_1k')
    group=p.add_mutually_exclusive_group(required=True)
    group.add_argument('--gene'); group.add_argument('--region')
    p.add_argument('--contact',help='Contact id from config')
    p.add_argument('--output',required=True,help='Output prefix, without extension')
    a=p.parse_args()
    data=run(load_config(a.config),a.source,a.gene or a.region,'gene' if a.gene else 'region',contact_id=a.contact)
    prefix=Path(a.output); prefix.parent.mkdir(parents=True,exist_ok=True)
    Path(str(prefix)+'_loops.csv').write_bytes(data['result'].csv())
    for fmt,content in data['exports'].items(): Path(str(prefix)+'.'+fmt).write_bytes(content)
    for warning in data['warnings']: print('NOTICE:',warning)
    if not data['exports']: raise SystemExit(2)
if __name__=='__main__': main()
