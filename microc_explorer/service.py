"""Shared orchestration with process-local admission, a shared plot cache and per-request state."""
from __future__ import annotations
import hashlib
import logging
import threading
from io import BytesIO
from .core import query,region_genes
from .readers import read_bigwig,read_contact
from .plotting import render_exports,FORMATS
from .ideogram import load_cytobands
from .safety import REQUEST_GATE,GLOBAL_LIMITER,MAX_PLOT_LOOPS,validate_request,report_failure
from .store import PLOT_CACHE

def run(cfg,source,value,mode,tracks=None,contact_id=None,*,bins=1800,statistic='max',
        exact=True,shared_scale=True,resolution=None,highlight_query=True,formats=FORMATS,
        signal_reader=read_bigwig,contact_reader=read_contact):
    validate_request(source,value,mode,bins,statistic,resolution,cfg,tracks,contact_id,formats)
    with REQUEST_GATE.acquire():
        return _run(cfg,source,value,mode,tracks,contact_id,bins,statistic,exact,shared_scale,
                    resolution,highlight_query,tuple(formats),signal_reader,contact_reader)

def confirm(cfg,source,value,mode):
    """Same admission gate as full plots; no track I/O or rendering."""
    GLOBAL_LIMITER.check()
    with REQUEST_GATE.acquire():
        return query(cfg['database'],source,value,mode)

def _run(cfg,source,value,mode,tracks,contact_id,bins,statistic,exact,shared_scale,
         resolution,highlight_query,formats,signal_reader,contact_reader):
    result=query(cfg['database'],source,value,mode)
    payload={'result':result,'genes':[],'exports':{},'warnings':list(result.notices),'signals':[],'render':None}
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
    kwargs=dict(title=cfg.get('title','Micro-C regulatory landscape'),shared_scale=shared_scale,
        signal_summary=f'{statistic} ({"exact" if exact else "zoom summaries"})',
        headroom=cfg.get('headroom',1.5),highlight_query=highlight_query,
        cytobands=bands,track_height=cfg.get('track_height',.42),
        ideogram_height=cfg.get('ideogram_height',.55),ideogram_width=cfg.get('ideogram_width',1.0),
        loops_height=cfg.get('loops_height',.8),loops_width=cfg.get('loops_width',1.0))
    # Kept (small numeric arrays) so other formats can be produced later without re-reading tracks.
    payload['render']=dict(genes=genes,signals=payload['signals'],contact=contact,kwargs=kwargs)
    try:
        payload['exports']=render_exports(result,genes,payload['signals'],contact,formats=formats,**kwargs)
    except Exception as e:
        # Keep valid query results even when rendering fails.
        payload['render']=None
        payload['warnings'].append(report_failure('The plot could not be rendered; query results are still available',e))
    return payload

def render_more(result,render,formats):
    """Produce additional formats from the inputs of an earlier plot (no query, no track I/O)."""
    with REQUEST_GATE.acquire():
        return render_exports(result,render['genes'],render['signals'],render['contact'],
                              formats=tuple(formats),**render['kwargs'])

def cache_key(version,source,value,mode,tracks,contact_id,bins,statistic,exact,shared_scale,resolution,highlight_query):
    """Identity of a plot: dataset files + validated parameters. Output format is NOT part of it."""
    parts=(version,source,mode,value.strip().upper(),tuple(sorted(tracks)) if tracks is not None else None,
           contact_id,bins,statistic,bool(exact),bool(shared_scale),resolution,bool(highlight_query))
    return hashlib.sha256(repr(parts).encode()).hexdigest()

def get_plot(cfg,req,formats,*,version,signal_reader=read_bigwig,contact_reader=read_contact,cache=PLOT_CACHE):
    """Return {'key','payload','exports','from_cache','stored'} for the request.

    `req` holds source, value, mode, tracks, contact_id, bins, statistic, exact, shared_scale,
    resolution and highlight_query. A repeated identical request is served from the shared cache;
    a request that differs only in output format re-renders from cached numeric inputs.
    """
    formats=tuple(formats)
    validate_request(req['source'],req['value'],req['mode'],req['bins'],req['statistic'],req['resolution'],
                     cfg,req['tracks'],req['contact_id'],formats)
    key=cache_key(version,**req)
    entry=cache.get(key)
    if entry is not None and entry['render'] is not None:
        missing=[f for f in formats if f not in entry['exports']]
        if missing:
            new=render_more(entry['payload']['result'],entry['render'],missing)
            cache.add_exports(key,new); entry['exports'].update(new)
        payload=dict(entry['payload']); payload['warnings']=list(payload['warnings'])
        return dict(key=key,payload=payload,exports=entry['exports'],from_cache=True,stored=True)
    GLOBAL_LIMITER.check()
    payload=run(cfg,req['source'],req['value'],req['mode'],req['tracks'],req['contact_id'],bins=req['bins'],
                statistic=req['statistic'],exact=req['exact'],shared_scale=req['shared_scale'],
                resolution=req['resolution'],highlight_query=req['highlight_query'],formats=formats,
                signal_reader=signal_reader,contact_reader=contact_reader)
    payload.pop('signals',None)
    render=payload.pop('render',None); exports=payload.pop('exports')
    stored=False
    if exports and render is not None:
        stored=cache.put(key,payload,render,exports)
        if not stored: payload['warnings'].append('The plot is larger than the server export budget and cannot be shown. Use a smaller region.')
    return dict(key=key,payload=payload,exports=exports if stored else {},from_cache=False,stored=stored)

_warm_lock=threading.Lock(); _warmed=False

def warmup(cfg):
    """Once per process, off the request path: build Matplotlib's font cache and load cytobands so the
    first visitor does not pay for them (a cold container can otherwise add seconds)."""
    global _warmed
    with _warm_lock:
        if _warmed: return
        _warmed=True
    def work():
        try:
            from matplotlib.figure import Figure
            from matplotlib.backends.backend_agg import FigureCanvasAgg
            fig=Figure(figsize=(2,1)); FigureCanvasAgg(fig); fig.text(.1,.5,'warm-up 0123 ABC kb Mb',fontsize=9,weight='bold')
            fig.savefig(BytesIO(),format='svg'); fig.savefig(BytesIO(),format='png',dpi=50)
            load_cytobands(cfg.get('cytoband'))
        except Exception as e: logging.getLogger('microc_explorer').warning('warm-up skipped: %s',e)
    threading.Thread(target=work,name='microc-warmup',daemon=True).start()
