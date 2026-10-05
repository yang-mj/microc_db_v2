"""SciLifeLab Serve entry point: streamlit run app.py."""
import os
from pathlib import Path
import uuid
import time
import pandas as pd
import streamlit as st
from microc_explorer.core import load_config,available_sources,metadata,fingerprint
from microc_explorer.readers import read_bigwig,read_contact
from microc_explorer.service import run,confirm
from microc_explorer.safety import RequestError, report_failure
from microc_explorer.store import STORE, start_janitor
from microc_explorer.manual import render_manual
from microc_explorer import __version__

st.set_page_config(page_title='Micro-C Explorer',page_icon='🧬',layout='wide')

# Cache numeric data, not file handles, figures, mutable connections, or user results.
@st.cache_data(max_entries=128,ttl=900,show_spinner=False)
def cached_signal(track,region,bins,statistic,exact,file_version):
    return read_bigwig(track,region,bins,statistic,exact)

@st.cache_data(max_entries=8,ttl=300,show_spinner=False)
def cached_contact(track,region,resolution,file_version):
    return read_contact(track,region,resolution)

def signal_reader(track,region,bins,statistic,exact):
    return cached_signal(track,region,bins,statistic,exact,fingerprint(track['path']))

def contact_reader(track,region,resolution):
    return cached_contact(track,region,resolution,fingerprint(track['path']))

@st.cache_data(max_entries=8,ttl=300,show_spinner=False)
def cached_catalog(database,file_version):
    return metadata(database),available_sources(database)

start_janitor()  # Module-level guard starts exactly one thread per process.
if 'export_owner' not in st.session_state:
    st.session_state.export_owner=uuid.uuid4().hex
session_id=st.session_state.export_owner

@st.fragment(run_every=30)
def keep_exports_alive():
    current=st.session_state.get('plot_result')
    if current and current.get('has_exports') and not current.get('export_expired'):
        if STORE.get(st.session_state.export_owner) is None:
            current['export_expired']=True
            st.rerun()
keep_exports_alive()

st.title('Micro-C Explorer')
st.caption('Pediatric BCP-ALL · promoter-linked chromatin interactions · hg38')
page=st.sidebar.radio('Page',['Explore','User manual','About / dataset'])
if page=='User manual':
    render_manual()
    st.stop()

config_path=os.environ.get('MICROC_CONFIG',str(Path.cwd()/'config.json'))
try:
    cfg=load_config(config_path)
    meta,sources=cached_catalog(cfg['database'],fingerprint(cfg['database']))
except Exception as e:
    st.error(report_failure('Cannot load the configured dataset',e))
    st.info('Open User manual in the sidebar for installation guidance. Prepare the database and set MICROC_CONFIG as described there.')
    st.stop()
if not sources:
    st.error('No supported loop sources exist in this database. Check loopSource values in the input TSV.')
    st.stop()
if meta.get('sample_dataset'):
    st.warning('Example dataset: this package contains only the uploaded loop excerpts. Missing subtypes/results do not indicate biological absence. Prepare your full concat_loops_v2.tab before publication.')

if page=='About / dataset':
    st.header('About / dataset')
    st.write(f'Micro-C Explorer {__version__}')
    st.markdown('Curated promoter-linked CRE and noCRE chromatin loops in hg38. Open **User manual** for usage and interpretation.')
    st.markdown('''This database visualizes loop-calling results from **32 BCP-ALL samples**, integrated with multi-omics data as described in the paper.

The app supports gene-centric and region-centric exploration of regulatory landscapes, with general and subtype-specific queries. See the **User manual** for usage and content.

**Reference:** [bioRxiv preprint, version 1](https://www.biorxiv.org/content/10.1101/2025.07.22.666073v1).

**Contact:**  
Aneuploidy in Cancer Lab, Division of Clinical Genetics, Lund University Faculty of Medicine  
BMC, D14  
221 84 Lund, Sweden
''')
    st.stop()

labels={t['id']:t['label'] for t in cfg.get('bigwigs',[])}
contacts={t['id']:t['label'] for t in cfg.get('contacts',[])}
# Streamlit combines the Enter hint and character counter in this element.
# Scope the override to this form; retain the native 128-character limit and
# describe it in the accessible help tooltip instead.
st.markdown("""
<style>
[data-testid="stForm"] [data-testid="InputInstructions"] {
    display: none;
}
</style>
""", unsafe_allow_html=True)
with st.sidebar.form('query_form',enter_to_submit=False):
    source=st.selectbox('Dataset',list(sources),format_func=lambda k:sources[k])
    mode=st.radio('Search by',['Gene','Coordinates'],horizontal=True)
    value=st.text_input('Gene symbol or genomic interval',value='KRAS',max_chars=128,help='Enter a gene symbol (for example, KRAS) or a GENCODE ENSG gene ID (with or without a version suffix) or coordinates such as chr1:25,500,000-25,850,000. Maximum input length: 128 characters. All input coordinates are 0-based, end-exclusive.')
    submitted=st.form_submit_button('Plot region',type='primary')
    confirmed=st.form_submit_button('Confirm gene / interval')
    track_ids=st.multiselect('Signal tracks',list(labels),default=list(labels),format_func=labels.get)
    contact_id=st.selectbox('Contact map',['None']+list(contacts),format_func=lambda k:contacts.get(k,k))
    with st.expander('Plot settings'):
        bins=st.select_slider('Display bins',options=[600,1200,1800,2400,3600],value=1800)
        statistic=st.selectbox('Signal per bin',['max','mean'])
        exact=st.checkbox('Exact BigWig summaries',value=True,help='Turn off for faster precomputed zoom summaries; values may differ.')
        shared=st.checkbox('Match y-axis scales within each assay',value=True)
        highlight=st.checkbox('Highlight searched gene or interval',value=True)
        resolution=st.number_input('Hi-C resolution (bp; 0 = automatic)',min_value=0,max_value=10_000_000,value=0,step=1000)

version=(fingerprint(cfg['database']),fingerprint(cfg['_config_path']),
         fingerprint(cfg['cytoband']) if cfg.get('cytoband') else 'bundled-hg38',
         tuple(fingerprint(t['path']) for t in cfg.get('bigwigs',[])+cfg.get('contacts',[])))
if confirmed or submitted:
    now=time.monotonic()
    if now-st.session_state.get('last_submit',-float('inf'))<2:
        st.warning('Please wait two seconds between requests.')
        st.stop()
    st.session_state.last_submit=now
    st.session_state.pop('confirm_msg',None)

if confirmed:
    st.session_state.confirm_key=(version,source,mode,value)
    # Retain the checked records for CSV even when the expanded view exceeds 2 Mb.
    st.session_state.pop('plot_result',None)
    STORE.drop(session_id)
    try:
        r=confirm(cfg,source,value,'gene' if mode=='Gene' else 'region')
        what=f'{r.gene} at {r.seed.label}' if r.gene else r.seed.label
        text=(f'{what} · {len(r.loops):,} gene-associated loop records in {sources[source]} · '
              f'expanded view {r.region.size/1e6:.3f} Mb'+('' if r.can_plot else ' (over the 2 Mb plot limit; CSV results remain available)'))
        st.session_state.confirm_msg=('success' if r.can_plot else 'warning',text)
        st.session_state.plot_result=dict(result=r,genes=[],warnings=list(r.notices),
            has_exports=False,query_label=f'{sources[source]} · {value}',version=version)
    except RequestError as e: st.session_state.confirm_msg=('error',str(e))
    except Exception as e: st.session_state.confirm_msg=('error',report_failure('The input could not be checked',e))
if st.session_state.get('confirm_key')!=(version,source,mode,value):
    st.session_state.pop('confirm_msg',None)
if 'confirm_msg' in st.session_state:
    kind,text=st.session_state.confirm_msg
    getattr(st.sidebar,kind)(text)

if submitted:
    # Clear old results before a new request, including failures.
    st.session_state.pop('plot_result',None); STORE.drop(session_id); st.session_state.pop('confirm_msg',None)
    try:
        with st.spinner('Finding loops and drawing aligned tracks…'):
            payload=run(cfg,source,value,'gene' if mode=='Gene' else 'region',track_ids,
                        None if contact_id=='None' else contact_id,bins=bins,statistic=statistic,
                        exact=exact,shared_scale=shared,resolution=int(resolution) or None,highlight_query=highlight,
                        signal_reader=signal_reader,contact_reader=contact_reader)
        # Numeric arrays are cached; session state retains only results/metadata.
        payload.pop('signals',None)
        # Exports (MBs) live in the capped, janitor-cleaned store, not in session state.
        exports=payload.pop('exports')
        payload['has_exports']=bool(exports)
        if exports and not STORE.put(session_id,exports):
            payload['has_exports']=False
            payload['warnings'].append('The plot is larger than the server export budget and cannot be shown. Use a smaller region.')
        del exports
        payload['query_label']=f'{sources[source]} · {value}'
        payload['version']=version
        st.session_state.plot_result=payload
    except RequestError as e: st.error(str(e))
    except Exception as e: st.error(report_failure('The request could not be completed',e))

payload=st.session_state.get('plot_result')
if payload and payload['version']!=version:
    st.info('The dataset, track files or configuration changed. Submit your query again.')
    st.session_state.pop('plot_result',None)
    st.session_state.pop('confirm_msg',None)
    STORE.drop(session_id)
    payload=None
if payload and st.button('Clear current result'):
    STORE.drop(session_id)
    st.session_state.pop('plot_result',None)
    st.session_state.pop('confirm_msg',None)
    payload=None


if not payload:
    st.info('Choose a dataset and enter a gene symbol or coordinates, then select Plot region.')
else:
    result=payload['result']
    st.subheader(payload['query_label'])
    st.caption(f'Expanded view: {result.region.label} · {result.region.size/1e6:.3f} Mb · {len(result.loops):,} matching records')
    for warning in payload['warnings']: st.warning(warning)
    filename=f'chr{result.region.chrom}_{result.region.start}_{result.region.end}'
    # Always offer results, including blocked (>2 Mb) queries.
    st.download_button('Download loop results (CSV)',result.csv(),filename+'_loops.csv','text/csv')
    if payload['has_exports']:
        exports=STORE.get(session_id)
        if exports is None:
            st.info('This plot was released from server memory (idle connection or memory limit). Select Plot region to draw it again.')
        else:
            # SVG is vector: sharp at any width. Fall back to PNG if it is unusually large.
            if len(exports['svg'])<=8_000_000: st.image(exports['svg'].decode('utf-8'),width='stretch')
            else:
                st.image(exports['png'],width='stretch')
                st.caption('Large plot: showing a 300 dpi PNG preview. Download SVG for vector detail.')
            for col,fmt in zip(st.columns(3),('svg','pdf','png')):
                mime={'png':'image/png','pdf':'application/pdf','svg':'image/svg+xml'}[fmt]
                col.download_button(f'Download {fmt.upper()}',exports[fmt],filename+'.'+fmt,mime)
    with st.expander('Loop table'):
        st.dataframe(result.loops.drop(columns='id',errors='ignore'),hide_index=True,width='stretch')
    if payload['genes']:
        with st.expander('Genes and representative transcripts'):
            genes=pd.DataFrame(payload['genes'])
            st.dataframe(genes.drop(columns=['exons','cds']),hide_index=True,width='stretch')
            st.download_button('Download genes (CSV)',genes.to_csv(index=False),filename+'_genes.csv','text/csv')

# Streamlit retains script globals for fragment callbacks; do not let a preview
# dictionary keep a second application reference after store eviction/expiry.
exports = None
