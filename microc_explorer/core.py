from __future__ import annotations
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
import json
import re
import sqlite3
import time
import math
from urllib.parse import urlsplit
from .safety import RequestError, MAX_ROWS, QUERY_SECONDS
import pandas as pd

# GRCh38 primary chromosome lengths. Internal coordinates: zero-based, half-open.
_LENGTHS = [248956422,242193529,198295559,190214555,181538259,170805979,
            159345973,145138636,138394717,133797422,135086622,133275309,
            114364328,107043718,101991189,90338345,83257441,80373285,
            58617616,64444167,46709983,50818468,156040895,57227415,16569]
CHROMS = dict(zip([str(i) for i in range(1,23)]+['X','Y','MT'], _LENGTHS))
SOURCES = {'merged_1k':'All cases · merged 1 kb', 'HeH_10k':'High hyperdiploidy',
           'ER_10k':'ETV6::RUNX1', 'BA_10k':'BCR::ABL1', 'DUX4r_10k':'DUX4-r',
           'TP_10k':'TCF3::PBX1', 'KMT2Ar_10k':'KMT2A-r',
           'iAMP_10k':'iAMP21', 'nearHaploid_10k':'Near haploid'}
LOOP_COLUMNS = ['chr','start','end','targetChr','targetStart','targetEnd',
                'targetGeneID_canon','targetGene_canon','canon_annot','loopSource']

def normalize_gene_id(value):
    value = str(value).strip().upper()
    return re.sub(r'\.\d+$', '', value)


def chromosome(value):
    value = str(value).strip()
    if value.lower().startswith('chr'): value = value[3:]
    return 'MT' if value.upper() in ('M','MT') else value.upper()

@dataclass(frozen=True)
class Region:
    chrom: str
    start: int
    end: int
    def __post_init__(self):
        object.__setattr__(self, 'chrom', chromosome(self.chrom))
        if self.chrom not in CHROMS:
            raise RequestError('Use an hg38 primary chromosome: 1–22, X, Y or MT.')
        if not 0 <= self.start < self.end <= CHROMS[self.chrom]:
            raise RequestError(f'Invalid hg38 interval on chr{self.chrom}: {self.start:,}–{self.end:,}.')
    @property
    def size(self): return self.end - self.start
    @property
    def label(self): return f'chr{self.chrom}:{self.start:,}-{self.end:,}'

def parse_region(value):
    m = re.fullmatch(r'(?:chr)?([\w]+):([\d,]+)-([\d,]+)', value.strip(), re.I)
    if not m: raise RequestError('Enter a region such as chr13:28,000,000-28,400,000 (0-based, end-exclusive).')
    return Region(m[1], int(m[2].replace(',','')), int(m[3].replace(',','')))

@contextmanager
def connect(path, query_seconds=QUERY_SECONDS):
    uri = Path(path).resolve().as_uri() + '?mode=ro'
    con = sqlite3.connect(uri, uri=True, timeout=2)
    con.row_factory = sqlite3.Row
    con.execute('PRAGMA query_only=ON')
    con.execute('PRAGMA trusted_schema=OFF')
    con.execute('PRAGMA cache_size=-8192')  # Target page-cache size; not a hard process-memory cap.
    schema = con.execute("SELECT value FROM metadata WHERE key='schema_version'").fetchone()
    if not schema or json.loads(schema[0]) != 3:
        con.close()
        raise RequestError('This version requires a rebuilt database. Run microc-prepare with the new gene-ID loop table and GENCODE annotation.')
    con.execute('BEGIN')  # One consistent read snapshot for metadata and rows.
    deadline=time.monotonic()+query_seconds
    con.set_progress_handler(lambda: int(time.monotonic()>deadline),1000)
    try:
        yield con
    except sqlite3.OperationalError as exc:
        if 'interrupted' in str(exc): raise RequestError('The query exceeded the server time limit. Please use a smaller region.') from exc
        raise
    finally: con.close()

def metadata(path):
    with connect(path) as con:
        return {r['key']:json.loads(r['value']) for r in con.execute('SELECT * FROM metadata')}

def available_sources(path):
    with connect(path) as con:
        present = {r[0] for r in con.execute('SELECT DISTINCT loopSource FROM loops')}
    return {k:v for k,v in SOURCES.items() if k in present}

@dataclass
class QueryResult:
    seed: Region
    region: Region
    loops: pd.DataFrame
    gene: str | None
    notices: list[str]
    @property
    def can_plot(self): return self.region.size <= 2_000_000
    def csv(self): return self.loops.drop(columns=['id'],errors='ignore').to_csv(index=False).encode()


def _read_limited(con,sql,params,max_rows):
    cursor=con.execute(sql+' LIMIT ?',(*params,max_rows+1))
    rows=cursor.fetchall()
    if len(rows)>max_rows:
        raise RequestError(f'This query exceeds {max_rows:,} result records. Narrow the input; no partial results were returned.')
    return pd.DataFrame.from_records([tuple(r) for r in rows],columns=[c[0] for c in cursor.description])

def region_query_sql(con,source,seed,use_bounds=True):
    clauses=[];params=[]
    bounded=use_bounds and con.execute("SELECT 1 FROM sqlite_master WHERE name='anchor_bounds' AND type='table'").fetchone()
    for anchor,chrom,start,end in [(1,'chr','start','end'),(2,'targetChr','targetStart','targetEnd')]:
        bounds=con.execute('SELECT max_length FROM anchor_bounds WHERE loopSource=? AND chrom=? AND anchor=?',
            (source,seed.chrom,anchor)).fetchone() if bounded else None
        clause=f"SELECT id FROM loops WHERE loopSource=? AND canon_annot IN ('CRE','noCRE') AND {chrom}=? AND {start}<? AND {end}>?"
        params.extend((source,seed.chrom,seed.end,seed.start))
        if bounds:
            # end = start + length, length <= L. For end > query.start,
            # start must be >= query.start-L+1. The exact end predicate stays.
            clause+=f' AND {start}>=?'
            params.append(max(0,seed.start-int(bounds[0])+1))
        clauses.append(clause)
    return 'SELECT * FROM loops WHERE id IN ('+' UNION '.join(clauses)+') ORDER BY id',tuple(params)

def query(path, source, value, mode='gene', max_rows=MAX_ROWS):
    """Gene-associated CRE and noCRE loops. One expansion from the original query, never recursive."""
    if not isinstance(value,str) or not 1<=len(value.strip())<=128:
        raise RequestError('Enter a query of at most 128 characters.')
    if not isinstance(max_rows,int) or not 1<=max_rows<=MAX_ROWS: raise RequestError('Invalid result limit.')
    if source not in SOURCES: raise RequestError(f'Unsupported loopSource: {source}')
    notices = []
    gene = value.strip().upper() if mode == 'gene' else None
    with connect(path) as con:
        if mode == 'gene':
            if not re.fullmatch(r'[A-Z0-9_.-]+',gene): raise RequestError('Enter a gene symbol using letters, numbers, hyphens, periods or underscores.')
            if gene.startswith('ENSG'):
                gid = normalize_gene_id(gene)
                matches = con.execute('SELECT DISTINCT gene_id,symbol FROM genes WHERE gene_id=?', (gid,)).fetchall()
            else:
                matches = con.execute('SELECT DISTINCT gene_id,symbol FROM genes WHERE symbol=?', (gene,)).fetchall()
            if not matches: raise RequestError(f'{gene} is not present in the prepared protein-coding annotation.')
            if len(matches)!=1: raise RequestError(f'{gene} resolves to multiple gene IDs. Search using an ENSG ID.')
            gid,gene = matches[0]
            loci = con.execute('SELECT chrom,MIN(start),MAX(end) FROM genes WHERE gene_id=? GROUP BY chrom', (gid,)).fetchall()
            if not loci: raise RequestError(f'{gene} is not present in the prepared gene annotation.')
            if len(loci) != 1:
                raise RequestError(f'{gene} has multiple chromosome locations. Please search by coordinates.')
            seed = Region(*loci[0])
            loops = _read_limited(con,'''SELECT l.* FROM loops l JOIN loop_genes g ON l.id=g.loop_id
                WHERE g.gene_id=? AND l.loopSource=? AND l.canon_annot IN ('CRE','noCRE') ORDER BY l.id''', (gid,source),max_rows)
        elif mode == 'region':
            seed = parse_region(value)
            sql,params=region_query_sql(con,source,seed)
            loops = _read_limited(con,sql,params,max_rows)
        else: raise RequestError('mode must be gene or region')
    a,b = seed.start,seed.end
    if len(loops):
        cis = (loops.chr == seed.chrom) & (loops.targetChr == seed.chrom)
        if (~cis).any():
            notices.append(f'{int((~cis).sum())} cross-chromosome/off-locus records are included in the table, but cannot be drawn on this single-chromosome axis.')
        visible = loops[cis]
        if len(visible):
            a = min(a,int(visible[['start','targetStart']].min().min()))
            b = max(b,int(visible[['end','targetEnd']].max().max()))
    else: notices.append('No matching gene-associated CRE/noCRE loops; genes and available signal tracks can still be plotted.')
    region = Region(seed.chrom,max(0,a-200_000),min(CHROMS[seed.chrom],b+200_000))
    result = QueryResult(seed,region,loops,gene,notices)
    if not result.can_plot:
        notices.append(f'Expanded region is {region.size/1e6:.3f} Mb, exceeding the 2 Mb limit. Plotting is stopped; all query results remain downloadable.')
    return result

def region_genes(path, region):
    with connect(path) as con:
        return [dict(r) for r in con.execute('SELECT * FROM genes WHERE chrom=? AND start<? AND end>? ORDER BY start,symbol',
                                           (region.chrom,region.end,region.start))]

def load_config(path):
    path = Path(path).resolve()
    cfg = json.loads(path.read_text())
    if not isinstance(cfg,dict) or not isinstance(cfg.get('database'),str) or not cfg['database'].strip(): raise ValueError('Invalid database configuration.')
    if not 1<=float(cfg.get('headroom',1.5))<=10: raise ValueError('headroom must be 1–10.')
    if not .2<=float(cfg.get('track_height',.42))<=3: raise ValueError('track_height must be 0.2–3 inches.')
    if not .2<=float(cfg.get('ideogram_height',.55))<=1.5: raise ValueError('ideogram_height must be 0.2–1.5 inches.')
    if not .2<=float(cfg.get('ideogram_width',1.0))<=1: raise ValueError('ideogram_width must be 0.2–1 (fraction of plot width).')
    if not .2<=float(cfg.get('loops_height',.8))<=4: raise ValueError('loops_height must be 0.2–4 inches.')
    if float(cfg.get('loops_width',1.0))!=1: raise ValueError('loops_width must be 1 to keep loops aligned with genomic tracks; use loops_height to make loops shorter.')
    for key,default in [('headroom',1.5),('track_height',.42),('ideogram_height',.55),('ideogram_width',1.),('loops_height',.8),('loops_width',1.)]:
        cfg[key]=float(cfg.get(key,default))
    if cfg.get('cytoband') is not None and (not isinstance(cfg['cytoband'],str) or '://' in cfg['cytoband']):
        raise ValueError('cytoband must be a local file path.')
    for key in ('bigwigs','contacts'):
        entries=cfg.get(key,[])
        if not isinstance(entries,list) or len(entries)>16: raise ValueError('At most 16 curated entries per track type.')
        ids=[]
        for entry in entries:
            if not isinstance(entry,dict): raise ValueError('Invalid track entry.')
            if not re.fullmatch(r'[A-Za-z0-9_-]{1,64}',str(entry.get('id',''))): raise ValueError('Invalid track id.')
            if not isinstance(entry.get('label'),str) or not 1<=len(entry['label'])<=100: raise ValueError('Invalid track label.')
            if not isinstance(entry.get('path'),str) or not entry['path']: raise ValueError('Invalid track path.')
            url=urlsplit(entry['path'].split('::',1)[0] if key=='contacts' else entry['path'])
            if url.scheme:
                if url.scheme!='https' or not url.hostname or url.username or url.password:
                    raise ValueError('Remote curated tracks must use HTTPS without embedded credentials.')
                if key=='contacts': raise ValueError('Contact matrices must be local files.')
            lo,hi=entry.get('min_value'),entry.get('max_value')
            if (lo is None)!=(hi is None): raise ValueError('Set both min_value and max_value, or neither.')
            if lo is not None and not (math.isfinite(float(lo)) and math.isfinite(float(hi)) and float(lo)<float(hi)):
                raise ValueError('Fixed signal limits must be finite and increasing.')
            ids.append(entry['id'])
        if len(ids)!=len(set(ids)): raise ValueError('Track ids must be unique within each track type.')
    def resolve(p):
        if not p or str(p).startswith(('https://','http://')): return p
        base,separator,group=str(p).partition('::')
        return str((path.parent / base).resolve())+separator+group
    cfg['database'] = resolve(cfg['database'])
    for t in cfg.get('bigwigs',[]): t['path'] = resolve(t['path'])
    for t in cfg.get('contacts',[]): t['path'] = resolve(t['path'])
    # Explicit local override; otherwise use bundled hg38 bands (no runtime network).
    cfg['cytoband'] = resolve(cfg['cytoband']) if cfg.get('cytoband') else None
    cfg['_config_path'] = str(path)
    return cfg

def fingerprint(path):
    """Cache identity; remote URL must be versioned when remote content changes."""
    if path.startswith(('https://','http://')): return path
    p = Path(path.split('::',1)[0])
    if not p.exists(): return (str(p),None)
    stat = p.stat()
    return (path,stat.st_size,stat.st_mtime_ns,stat.st_ino)
