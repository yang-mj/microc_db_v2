"""SciLifeLab Serve entry point: streamlit run app.py."""
import os
from pathlib import Path
import re
import time
import pandas as pd
import streamlit as st
from microc_explorer.core import load_config,available_sources,metadata,fingerprint
from microc_explorer.readers import read_bigwig,read_contact
from microc_explorer.service import confirm,get_plot,warmup
from microc_explorer.safety import RequestError, report_failure
from microc_explorer.store import PLOT_CACHE, start_janitor
from microc_explorer.manual import render_manual
from microc_explorer import __version__

st.set_page_config(page_title='Micro-C Explorer',page_icon='🧬',layout='wide',initial_sidebar_state='expanded')

# Static stylesheet only: no user-supplied text is ever placed in HTML.
st.markdown('''<style>
.block-container{padding-top:2.2rem;padding-bottom:3rem;max-width:1500px}
[data-testid="stSidebar"]{border-right:1px solid #e3e8ef}
[data-testid="stSidebar"] [data-testid="stForm"]{padding:0}
[data-testid="stMetric"]{background:#f6f8fb;border:1px solid #e3e8ef;border-radius:12px;padding:.7rem 1rem}
[data-testid="stMetricLabel"]{color:#667085}
h1{letter-spacing:-.02em;margin-bottom:0}
.mc-sub{color:#667085;margin:.15rem 0 1.1rem 0;font-size:.98rem}
.mc-step{color:#4775a5;font-weight:700;font-size:.8rem;letter-spacing:.06em;text-transform:uppercase}
[data-testid="stImage"] img{border:1px solid #e3e8ef;border-radius:10px;background:#fff}
</style>''',unsafe_allow_html=True)

FORMAT_CHOICES={'SVG · vector, sharp at any zoom':'svg','PNG · high resolution (300 dpi)':'png'}
FORMAT_NAMES={'svg':'SVG (vector)','png':'PNG (300 dpi)','pdf':'PDF (vector)'}
MIME={'png':'image/png','pdf':'application/pdf','svg':'image/svg+xml'}
SVG_PREVIEW_LIMIT=8_000_000

def md_escape(text):
    """Visitor text is shown as plain text: neutralise Markdown/LaTeX/emoji-shortcode syntax."""
    text=''.join(ch for ch in str(text) if ch.isprintable())
    return re.sub(r'([!-/:-@\[-`{-~])',r'\\\1',text)

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

@st.cache_data(max_entries=4,ttl=10,show_spinner=False)
def dataset_version(paths):
    """File stats are re-checked at most every 10 s instead of on every widget rerun."""
    return tuple(fingerprint(p) if p!='bundled-hg38' else p for p in paths)

def fetch_format(cfg,recipe,version,fmt):
    """Runs when a download is clicked (worker thread, no Streamlit calls). Uses the shared cache."""
    return get_plot(cfg,recipe,(fmt,),version=version)['exports'][fmt]

start_janitor()  # Module-level guard starts exactly one thread per process.

st.title('Micro-C Explorer')
st.markdown('<div class="mc-sub">Pediatric BCP-ALL · promoter-linked chromatin interactions · hg38</div>',unsafe_allow_html=True)
page=st.sidebar.radio('Page',['Explore','User manual','About / dataset'],label_visibility='collapsed')
if page=='User manual':
    render_manual()
    st.stop()

config_path=os.environ.get('MICROC_CONFIG',str(Path.cwd()/'config.json'))
try:
    cfg=load_config(config_path)
    file_paths=(cfg['database'],cfg['_config_path'],cfg['cytoband'] if cfg.get('cytoband') else 'bundled-hg38',
                *[t['path'] for t in cfg.get('bigwigs',[])+cfg.get('contacts',[])])
    version=dataset_version(file_paths)
    meta,sources=cached_catalog(cfg['database'],version[0])
except Exception as e:
    st.error(report_failure('Cannot load the configured dataset',e))
    st.info('Open User manual in the sidebar for installation guidance. Prepare the database and set MICROC_CONFIG as described there.')
    st.stop()
if not sources:
    st.error('No supported loop sources exist in this database. Check loopSource values in the input TSV.')
    st.stop()
warmup(cfg)  # once per process, in the background
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

# ---------------------------------------------------------------- query form
# Pressing Enter in the text box submits the form's FIRST submit button, i.e. "Plot region".
with st.sidebar.form('query_form',enter_to_submit=True,border=False):
    st.markdown('<span class="mc-step">Query</span>',unsafe_allow_html=True)
    source=st.selectbox('Dataset',list(sources),format_func=lambda k:sources[k])
    mode=st.radio('Search by',['Gene','Coordinates'],horizontal=True)
    value=st.text_input('Gene symbol or genomic interval',value='KRAS',max_chars=128,
        help='Enter a gene symbol (for example, KRAS) or a GENCODE ENSG gene ID (with or without a version suffix) or coordinates such as chr1:25,500,000-25,850,000. Maximum input length: 128 characters. All input coordinates are 0-based, end-exclusive.')
    c_plot,c_check=st.columns([3,2])
    submitted=c_plot.form_submit_button('Plot region',type='primary',width='stretch',help='Press Enter in the text box as a shortcut.')
    confirmed=c_check.form_submit_button('Confirm',width='stretch',help='Confirm gene / interval: shows the resolved location and loop count without reading any tracks.')
    st.markdown('<span class="mc-step">Tracks</span>',unsafe_allow_html=True)
    track_ids=st.multiselect('Signal tracks',list(labels),default=list(labels),format_func=labels.get)
    contact_id='None'
    if contacts: contact_id=st.selectbox('Contact map',['None']+list(contacts),format_func=lambda k:contacts.get(k,k))
    with st.expander('Plot settings'):
        fmt_label=st.radio('Figure format',list(FORMAT_CHOICES),help='The format shown on the page and offered first for download. The other formats are produced on demand when you download them.')
        bins=st.select_slider('Display bins',options=[600,1200,1800,2400,3600],value=1800)
        statistic=st.selectbox('Signal per bin',['max','mean'])
        exact=st.checkbox('Exact BigWig summaries',value=True,help='Base-level maximum or mean per bin. Turn off to use the precomputed zoom summaries; values may differ slightly.')
        shared=st.checkbox('Match y-axis scales within each assay',value=True)
        highlight=st.checkbox('Highlight searched gene or interval',value=True)
        resolution=st.number_input('Hi-C resolution (bp; 0 = automatic)',min_value=0,max_value=10_000_000,value=0,step=1000) if contacts else 0
fmt=FORMAT_CHOICES[fmt_label]
st.sidebar.caption(f'Micro-C Explorer v{__version__}')

if confirmed or submitted:
    now=time.monotonic()
    if now-st.session_state.get('last_submit',-float('inf'))<2:
        st.warning('Please wait two seconds between requests.')
        confirmed=submitted=False   # keep showing the current result
    else:
        st.session_state.last_submit=now
        st.session_state.pop('confirm_msg',None)

if confirmed:
    st.session_state.confirm_key=(version,source,mode,value)
    # Retain the checked records for CSV even when the expanded view exceeds 2 Mb.
    st.session_state.pop('plot_result',None)
    try:
        r=confirm(cfg,source,value,'gene' if mode=='Gene' else 'region')
        what=f'{r.gene} at {r.seed.label}' if r.gene else r.seed.label
        text=(f'{what} · {len(r.loops):,} gene-associated loop records in {sources[source]} · '
              f'expanded view {r.region.size/1e6:.3f} Mb'+('' if r.can_plot else ' (over the 2 Mb plot limit; CSV results remain available)'))
        st.session_state.confirm_msg=('success' if r.can_plot else 'warning',text)
        st.session_state.plot_result=dict(result=r,genes=[],warnings=list(r.notices),has_exports=False,key=None,
            query_label=f'{sources[source]} · {value}',version=version)
    except RequestError as e: st.session_state.confirm_msg=('error',str(e))
    except Exception as e: st.session_state.confirm_msg=('error',report_failure('The input could not be checked',e))
if st.session_state.get('confirm_key')!=(version,source,mode,value):
    st.session_state.pop('confirm_msg',None)
if 'confirm_msg' in st.session_state:
    kind,text=st.session_state.confirm_msg
    getattr(st.sidebar,kind)(text)

if submitted:
    # Clear old results before a new request, including failures.
    st.session_state.pop('plot_result',None); st.session_state.pop('confirm_msg',None)
    recipe=dict(source=source,value=value,mode='gene' if mode=='Gene' else 'region',tracks=list(track_ids),
                contact_id=None if contact_id=='None' else contact_id,bins=bins,statistic=statistic,exact=exact,
                shared_scale=shared,resolution=int(resolution) or None,highlight_query=highlight)
    try:
        t0=time.perf_counter()
        with st.spinner('Finding loops and drawing aligned tracks…'):
            out=get_plot(cfg,recipe,(fmt,),version=version,signal_reader=signal_reader,contact_reader=contact_reader)
            if fmt=='svg' and out['stored'] and len(out['exports']['svg'])>SVG_PREVIEW_LIMIT:
                # Very large vector drawing: also prepare a PNG for the on-page preview.
                out=get_plot(cfg,recipe,('svg','png'),version=version,signal_reader=signal_reader,contact_reader=contact_reader)
        payload=out['payload']
        payload.update(has_exports=bool(out['exports']),key=out['key'],fmt=fmt,recipe=recipe,version=version,
                       query_label=f'{sources[source]} · {value}',elapsed=time.perf_counter()-t0,from_cache=out['from_cache'])
        st.session_state.plot_result=payload
    except RequestError as e: st.error(str(e))
    except Exception as e: st.error(report_failure('The request could not be completed',e))

payload=st.session_state.get('plot_result')
if payload and payload['version']!=version:
    st.info('The dataset, track files or configuration changed. Submit your query again.')
    st.session_state.pop('plot_result',None)
    st.session_state.pop('confirm_msg',None)
    payload=None

# ---------------------------------------------------------------- results
@st.fragment
def loop_table(df):
    if st.toggle('Show loop table',key='show_loops',help='Loaded only when switched on, so the page stays light.'):
        st.dataframe(df.drop(columns='id',errors='ignore'),hide_index=True,width='stretch')

@st.fragment
def gene_table(genes,filename):
    if st.toggle('Show genes and representative transcripts',key='show_genes'):
        table=pd.DataFrame(genes)
        st.dataframe(table.drop(columns=['exons','cds']),hide_index=True,width='stretch')
        st.download_button('Download genes (CSV)',lambda:table.to_csv(index=False),filename+'_genes.csv','text/csv',on_click='ignore')

if not payload:
    st.info('Choose a dataset, enter a gene symbol or coordinates in the sidebar, then press **Enter** or click **Plot region**.')
    cols=st.columns(3)
    for col,(title,text) in zip(cols,[
        ('1 · Search','Gene symbol (e.g. KRAS), ENSG ID, or hg38 coordinates (0-based, end-exclusive).'),
        ('2 · Choose tracks','Pick ChIP signal tracks and plot settings, including SVG or high-resolution PNG output.'),
        ('3 · Explore & export','Inspect aligned genes, loops and signal; download the figure and loop table.')]):
        with col.container(border=True):
            st.markdown(f'**{title}**'); st.caption(text)
else:
    result=payload['result']
    st.subheader(md_escape(payload['query_label']))
    cre=int((result.loops.canon_annot=='CRE').sum()); no=int((result.loops.canon_annot=='noCRE').sum())
    m1,m2,m3=st.columns(3)
    m1.metric('Loop records',f'{len(result.loops):,}')
    m2.metric('CRE · noCRE',f'{cre:,} · {no:,}')
    m3.metric('Expanded view',f'{result.region.size/1e6:.3f} Mb')
    timing=''
    if payload.get('elapsed') is not None and payload.get('has_exports'):
        timing=f' · ready in {payload["elapsed"]:.1f} s'+(' (shared cache)' if payload.get('from_cache') else '')
    st.caption(f'{result.region.label}{timing}')
    for warning in payload['warnings']: st.warning(warning)
    filename=f'chr{result.region.chrom}_{result.region.start}_{result.region.end}'
    entry=PLOT_CACHE.get(payload['key']) if payload.get('has_exports') else None
    if payload.get('has_exports') and entry is None:
        st.info('This plot was released from server memory (idle time or memory limit). Select Plot region to draw it again; it takes about a second.')
    shown=entry['exports'] if entry else {}
    chosen=payload.get('fmt','svg')
    if shown:
        if chosen=='svg' and len(shown['svg'])<=SVG_PREVIEW_LIMIT: st.image(shown['svg'].decode('utf-8'),width='stretch')
        elif chosen=='png' and 'png' in shown: st.image(shown['png'],width='stretch')
        elif 'png' in shown:
            st.image(shown['png'],width='stretch')
            st.caption('Large plot: showing a 300 dpi PNG preview. Download SVG for vector detail.')
    # Always offer results, including blocked (>2 Mb) queries. Files are produced when clicked.
    d1,d2,d3,_=st.columns([2,2,2,3])
    if shown and chosen in shown:
        d1.download_button(f'Download {FORMAT_NAMES[chosen]}',shown[chosen],filename+'.'+chosen,MIME[chosen],type='primary',on_click='ignore',width='stretch')
    d2.download_button('Loop results (CSV)',result.csv,filename+'_loops.csv','text/csv',on_click='ignore',width='stretch')
    if shown:
        with d3.popover('Other formats',width='stretch'):
            for other in ('svg','png','pdf'):
                if other==chosen: continue
                st.download_button(FORMAT_NAMES[other],(shown[other] if other in shown else
                    (lambda f=other:fetch_format(cfg,payload['recipe'],payload['version'],f))),
                    filename+'.'+other,MIME[other],on_click='ignore',key=f'dl_{other}',width='stretch')
    loop_table(result.loops)
    if payload['genes']: gene_table(payload['genes'],filename)
