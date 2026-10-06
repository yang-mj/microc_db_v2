"""Regional I/O. Handles are opened per call and never shared across sessions."""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import math
import threading

_HIC_LOCK = threading.RLock()
import numpy as np
from .core import Region, chromosome
from .safety import MAX_SIGNAL_BINS, RequestError

@dataclass
class Signal:
    edges: np.ndarray
    values: np.ndarray
    label: str
    color: str
    group: str
    min_value: float | None = None
    max_value: float | None = None

@dataclass
class Contact:
    edges: np.ndarray
    matrix: np.ndarray
    label: str
    normalization: str
    resolution: int

def matching_chrom(names, chrom):
    matches=[n for n in names if chromosome(n)==chromosome(chrom)]
    if len(matches)!=1: raise ValueError(f'Expected one match for chromosome {chrom}; found {matches}.')
    return matches[0]

def _binned_exact(bw, chrom, start, end, n, statistic):
    """Exact per-bin max/mean from one base-level read, binned with NumPy.

    Equivalent to ``bw.stats(..., exact=True)`` (verified value-for-value on dense and gapped
    tracks) but one read of the region instead of one lookup per display bin: ~10-30x faster.
    Uncovered bases are missing (NaN), never zero; a bin without coverage stays NaN.
    """
    v = np.asarray(bw.values(chrom, start, end, numpy=True), dtype=np.float64)
    v[~np.isfinite(v)] = np.nan
    first = (np.arange(n, dtype=np.int64) * (end - start)) // n     # first base of each bin
    if statistic == 'max':
        out = np.fmax.reduceat(v, first)                            # fmax ignores NaN
    else:
        ok = np.isfinite(v)
        total = np.add.reduceat(np.where(ok, v, 0.0), first)
        count = np.add.reduceat(ok.astype(np.int64), first)
        out = np.full(n, np.nan)
        np.divide(total, count, out=out, where=count > 0)
    out[~np.isfinite(out)] = np.nan
    return out

def read_bigwig(track, region, bins=1800, statistic='max', exact=True):
    if not isinstance(bins,int) or not 1<=bins<=MAX_SIGNAL_BINS: raise RequestError('Invalid signal bin count.')
    if region.size>2_000_000: raise RequestError('Signal region exceeds 2 Mb.')
    import pyBigWig
    if statistic not in ('mean','max'): raise ValueError('BigWig summary must be mean or max.')
    path=track['path']
    if path.startswith(('https://','http://')) and not pyBigWig.remote:
        raise RuntimeError('This pyBigWig build lacks remote support; use local files or a libcurl-enabled build.')
    bw=pyBigWig.open(path)
    if bw is None: raise ValueError(f'Cannot open {track["label"]}.')
    try:
        chrom=matching_chrom(bw.chroms(),region.chrom)
        if region.end>bw.chroms(chrom): raise ValueError('BigWig chromosome is shorter than this hg38 region; check assembly.')
        n=min(int(bins),region.size)
        if exact and getattr(pyBigWig,'numpy',0):
            values=_binned_exact(bw,chrom,region.start,region.end,n,statistic)
        else:
            values=bw.stats(chrom,region.start,region.end,nBins=n,type=statistic,exact=exact)
            # Missing values remain NaN; absence of coverage is not a measured zero.
            values=np.array([np.nan if v is None else v for v in values],dtype=float)
            values[~np.isfinite(values)]=np.nan
    finally: bw.close()
    edges=np.linspace(region.start,region.end,n+1)
    return Signal(edges,values,track['label'],track.get('color','#3878a8'),track.get('group',track['label']),track.get('min_value'),track.get('max_value'))

def _choose_resolution(resolutions, region, requested, max_bins):
    resolutions=sorted(int(x) for x in resolutions)
    if requested is not None:
        if int(requested) not in resolutions: raise ValueError(f'Resolution {requested} unavailable; choices: {resolutions}')
        chosen=int(requested)
    else:
        options=[r for r in resolutions if math.ceil(region.end/r)-region.start//r<=max_bins]
        if not options: raise ValueError('No resolution satisfies the contact-matrix memory limit.')
        chosen=options[0]
    if math.ceil(region.end/chosen)-region.start//chosen>max_bins:
        raise ValueError(f'Contact matrix exceeds {max_bins} bins. Choose a coarser resolution.')
    return chosen

def read_contact(track,region,requested=None,max_bins=600):
    """No fallback to raw counts if requested normalization is absent."""
    if region.size>2_000_000: raise RequestError('Contact region exceeds 2 Mb.')
    if not isinstance(max_bins,int) or not 1<=max_bins<=600: raise RequestError('Invalid contact bin limit.')
    path=track['path']
    kind=track.get('format',Path(path.split('::',1)[0]).suffix.lstrip('.')).lower()
    if kind in ('mcool','cool'):
        import cooler
        if '::' in path:
            c=cooler.Cooler(path)
            resolutions=[c.binsize]
            uri=lambda r:path
        elif kind=='mcool':
            groups=cooler.fileops.list_coolers(path)
            resolutions=[int(g.rsplit('/',1)[-1]) for g in groups if g.startswith('/resolutions/')]
            uri=lambda r:f'{path}::/resolutions/{r}'
        else:
            c=cooler.Cooler(path)
            resolutions=[c.binsize]
            uri=lambda r:path
        if None in resolutions: raise ValueError('Variable-width Cooler bins are not supported.')
        resolution=_choose_resolution(resolutions,region,requested,max_bins)
        c=cooler.Cooler(uri(resolution))
        chrom=matching_chrom(c.chromnames,region.chrom)
        if region.end>c.chromsizes[chrom]: raise ValueError('Contact file chromosome is shorter than the query.')
        weight=track.get('weight_column','weight')
        if weight not in c.bins().columns: raise ValueError(f'ICE weight column {weight!r} is missing. Set weight_column to the existing ICE weights.')
        bins=c.bins().fetch((chrom,region.start,region.end))
        weights=bins[weight].to_numpy(dtype=float)
        if not np.isfinite(weights).any(): raise ValueError('No finite ICE weights in the selected region.')
        matrix=c.matrix(balance=weight,divisive_weights=bool(track.get('divisive_weights',False))).fetch((chrom,region.start,region.end))
        edges=np.r_[bins.start.to_numpy(),int(bins.end.iloc[-1])]
        normalization=f'ICE ({weight})'
    elif kind=='hic':
        # Serialize native .hic calls conservatively; do not share native handles.
        with _HIC_LOCK:
            return _read_hic(track,region,requested,max_bins)
    else: raise ValueError('Contact format must be cool, mcool or hic.')
    return _validated_contact(edges,matrix,track['label'],normalization,resolution)

def _read_hic(track,region,requested,max_bins):
    path=track['path']
    import hicstraw
    hic=hicstraw.HiCFile(path)
    chroms={c.name:c.length for c in hic.getChromosomes() if c.name.lower()!='all'}
    chrom=matching_chrom(chroms,region.chrom)
    if region.end>chroms[chrom]: raise ValueError('Contact file chromosome is shorter than the query.')
    resolution=_choose_resolution(hic.getResolutions(),region,requested,max_bins)
    normalization=track.get('normalization','SCALE')
    mzd=hic.getMatrixZoomData(chrom,chrom,'observed',normalization,'BP',resolution)
    lo=region.start//resolution*resolution
    hi=min(math.ceil(region.end/resolution)*resolution,chroms[chrom])
    n=math.ceil((hi-lo)/resolution)
    matrix=np.zeros((n,n),dtype=float)
    # Sparse records avoid ambiguous inclusive end-bin behavior of getRecordsAsMatrix.
    records=mzd.getRecords(lo,hi-1,lo,hi-1)
    if not records: raise ValueError(f'No {normalization} contact records returned; check normalization and region.')
    for rec in records:
        i,j=(rec.binX-lo)//resolution,(rec.binY-lo)//resolution
        if 0<=i<n and 0<=j<n:
            matrix[i,j]=rec.counts; matrix[j,i]=rec.counts
    edges=np.minimum(lo+np.arange(n+1)*resolution,chroms[chrom])
    return _validated_contact(edges,matrix,track['label'],normalization,resolution)

def _validated_contact(edges,matrix,label,normalization,resolution):
    matrix=np.asarray(matrix,dtype=float)
    matrix[~np.isfinite(matrix)]=np.nan
    if not np.isfinite(matrix).any(): raise ValueError('Contact matrix has no finite normalized values.')
    if matrix.shape!=(len(edges)-1,len(edges)-1): raise ValueError('Contact matrix and genomic bins do not align.')
    return Contact(edges,matrix,label,normalization,resolution)
