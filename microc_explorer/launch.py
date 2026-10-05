"""Launch the installed Streamlit app without needing a source checkout."""
import argparse
from importlib.resources import files
import os
from pathlib import Path
import subprocess
import sys

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',default=os.environ.get('MICROC_CONFIG','config.json'))
    parser.add_argument('--host',default='127.0.0.1')
    parser.add_argument('--port',type=int,default=8501)
    parser.add_argument('--init-config',metavar='PATH',help='Write the supplied configuration template and exit; never overwrite.')
    args=parser.parse_args()
    if args.init_config:
        target=Path(args.init_config)
        target.parent.mkdir(parents=True,exist_ok=True)
        with target.open('x') as handle:
            handle.write(files('microc_explorer').joinpath('assets/config.example.json').read_text())
        print('Created',target);return
    if not 1<=args.port<=65535:parser.error('Port must be 1–65535.')
    config=Path(args.config).resolve()
    if not config.is_file():parser.error('Configuration missing. Use --init-config config.json, then edit data paths.')
    env=dict(os.environ,MICROC_CONFIG=str(config))
    app=str(files('microc_explorer').joinpath('webapp.py'))
    command=[sys.executable,'-m','streamlit','run',app,
             '--server.address',args.host,'--server.port',str(args.port),
             '--server.headless=true','--server.enableCORS=true','--server.enableXsrfProtection=true',
             '--server.enableStaticServing=false','--server.disconnectedSessionTTL=60',
             '--server.fileWatcherType=none','--server.maxUploadSize=1','--server.maxMessageSize=64',
             '--client.showErrorDetails=none','--client.toolbarMode=viewer',
             '--runner.fastReruns=false','--browser.gatherUsageStats=false']
    try:raise SystemExit(subprocess.call(command,env=env))
    except KeyboardInterrupt:raise SystemExit(130)

if __name__=='__main__':main()
