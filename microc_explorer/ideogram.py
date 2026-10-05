"""UCSC-style chromosome ideogram drawn from a UCSC cytoBand table (hg38).

Bands are flat rectangles with the band name inside (when it fits), the
centromere is a maroon bow-tie, and the plotted region is a red marker.
"""
from __future__ import annotations
import gzip
from functools import lru_cache
from pathlib import Path
from importlib.resources import files
from matplotlib.patches import Rectangle, Polygon
from .core import chromosome, CHROMS

# Approximate UCSC Giemsa colours.
STAIN_COLORS = {'gneg': '#ffffff', 'gpos25': '#c8c8c8', 'gpos50': '#969696',
                'gpos75': '#646464', 'gpos100': '#000000', 'gvar': '#dcdcdc',
                'stalk': '#647fa4', 'acen': '#8c2d2d'}
DARK_STAINS = {'gpos75', 'gpos100', 'stalk', 'acen'}
MARKER_COLOR = '#e01b24'

@lru_cache(maxsize=4)
def _load(path, mtime_ns, size, inode):
    opener = gzip.open if path.endswith('.gz') else open
    bands = {}
    with opener(path, 'rt') as f:
        for line in f:
            c = line.rstrip('\n').split('\t')
            if len(c) < 5 or line.startswith('#'): continue
            chrom = chromosome(c[0])
            if chrom not in CHROMS: continue  # skips alt/unplaced contigs
            start,end=int(c[1]),int(c[2])
            if not 0<=start<end<=CHROMS[chrom]: raise ValueError('Cytobands must use hg38 half-open coordinates.')
            if c[4].strip() not in STAIN_COLORS: raise ValueError('Unknown cytoband stain.')
            bands.setdefault(chrom, []).append((start,end,c[3],c[4].strip()))
    for v in bands.values():
        v.sort()
        if any(a[1]>b[0] for a,b in zip(v,v[1:])): raise ValueError('Overlapping cytobands.')
    return {key:tuple(value) for key,value in bands.items()}

def load_cytobands(path=None):
    """{chrom: [(start, end, name, stain), ...]}; cached until the file changes."""
    p = Path(path) if path else Path(str(files('microc_explorer').joinpath('assets/hg38_cytoBandIdeo.tsv')))
    stat=p.stat()
    return dict(_load(str(p),stat.st_mtime_ns,stat.st_size,stat.st_ino))

def region_band_label(bands, region):
    """'q12.2' or 'q12.2–q13.1' for the bands overlapped by the region."""
    hit = [b[2] for b in bands if b[0] < region.end and b[1] > region.start]
    if not hit: return ''
    return hit[0] if hit[0] == hit[-1] else f'{hit[0]}–{hit[-1]}'

def draw_ideogram(ax, bands, region, width_in, fontsize=7):
    """Draw the whole chromosome on `ax` (x in bp) with the region marked."""
    length = CHROMS[region.chrom]
    ax.set_xlim(0, length); ax.set_ylim(-.2, 1.2)
    ax.axis('off')
    for s, e, name, stain in bands:
        color = STAIN_COLORS.get(stain, '#ffffff')
        if stain == 'acen':
            # p-side triangle points right, q-side points left: a bow-tie.
            pts = [(s, 1), (e, .5), (s, 0)] if name.startswith('p') else [(e, 1), (s, .5), (e, 0)]
            ax.add_patch(Polygon(pts, closed=True, facecolor=color, edgecolor='black', lw=.6))
            continue
        ax.add_patch(Rectangle((s, 0), e - s, 1, facecolor=color, edgecolor='black', lw=.6))
        # Label only when the text fits inside the band (as UCSC does).
        if (e - s) / length * width_in > len(name) * fontsize * .58 / 72 + .05:
            ax.text((s + e) / 2, .5, name, ha='center', va='center', fontsize=fontsize,
                    color='white' if stain in DARK_STAINS else 'black', clip_on=True)
    # Exact region marker plus a visible pointer for very narrow windows.
    w = region.size
    x0 = region.start
    # Preserve true interval width; a pointer makes subpixel regions visible.
    ax.plot((region.start+region.end)/2,-.18,marker='^',ms=3,color=MARKER_COLOR,clip_on=False)
    ax.add_patch(Rectangle((x0, -.14), w, 1.28, facecolor=MARKER_COLOR, alpha=.25, edgecolor='none', zorder=5))
    ax.add_patch(Rectangle((x0, -.14), w, 1.28, fill=False, edgecolor=MARKER_COLOR, lw=1.3, zorder=6))
