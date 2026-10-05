"""Shared orchestration with process-local admission and per-request state."""
from .core import query,region_genes
from .readers import read_bigwig,read_contact
from .plotting import render_exports
from .ideogram import load_cytobands
from .safety import REQUEST_GATE,MAX_PLOT_LOOPS,validate_request,report_failure

def run(cfg,source,value,mode,tracks=None,contact_id=None,*,bins=1800,statistic='max',
        exact=True,shared_scale=True,resolution=None,highlight_query=True,
        signal_reader=read_bigwig,contact_reader=read_contact):
    validate_request(source,value,mode,bins,statistic,resolution,cfg,tracks,contact_id)
    with REQUEST_GATE.acquire():
        return _run(cfg,source,value,mode,tracks,contact_id,bins,statistic,exact,shared_scale,
                    resolution,highlight_query,signal_reader,contact_reader)

def confirm(cfg,source,value,mode):
    """Same admission gate as full plots; no track I/O or rendering."""
    with REQUEST_GATE.acquire():
        return query(cfg['database'],source,value,mode)

def _run(cfg,source,value,mode,tracks,contact_id,bins,statistic,exact,shared_scale,
         resolution,highlight_query,signal_reader,contact_reader):
    result=query(cfg['database'],source,value,mode)
    payload={'result':result,'genes':[],'exports':{},'warnings':list(result.notices),'signals':[]}
    if not result.can_plot:return payload
    if len(result.loops)>MAX_PLOT_LOOPS:
        payload['warnings'].append(f'This result has more than {MAX_PLOT_LOOPS:,} loops. Plotting was stopped to protect server resources; all matched records remain downloadable.')
        return payload
    genes=region_genes(cfg['database'],result.region)
    payload['genes']=genes
    chosen=set(tracks) if tracks is not None else {t['id'] for t in cfg.get('bigwigs',[])}
    for track in cfg.get('bigwigs',[]):
        if track['id'] not in chosen:continue
        try:payload['signals'].append(signal_reader(track,result.region,bins,statistic,exact))
        except Exception as e:payload['warnings'].append(report_failure(f'{track["label"]} could not be plotted',e))
    contact=None
    if contact_id:
        track=next(t for t in cfg.get('contacts',[]) if t['id']==contact_id)
        try:contact=contact_reader(track,result.region,resolution)
        except Exception as e:payload['warnings'].append(report_failure(f'{track["label"]} heatmap could not be plotted',e))
    bands=None
    try: bands=load_cytobands(cfg.get('cytoband')).get(result.region.chrom)
    except Exception as e: payload['warnings'].append(report_failure('The chromosome ideogram could not be loaded',e))
    try:
        payload['exports']=render_exports(result,genes,payload['signals'],contact,
            title=cfg.get('title','Micro-C regulatory landscape'),shared_scale=shared_scale,
            signal_summary=f'{statistic} ({"exact" if exact else "zoom summaries"})',
            headroom=cfg.get('headroom',1.5),highlight_query=highlight_query,
            cytobands=bands,track_height=cfg.get('track_height',.42),
            ideogram_height=cfg.get('ideogram_height',.55),ideogram_width=cfg.get('ideogram_width',1.0),
            loops_height=cfg.get('loops_height',.8),loops_width=cfg.get('loops_width',1.0))
    except Exception as e:
        # Keep valid query results even when rendering fails.
        payload['warnings'].append(report_failure('The plot could not be rendered; query results are still available',e))
    return payload
