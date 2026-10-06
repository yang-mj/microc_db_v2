"""Static renderer; plain data objects also serve a future interactive renderer."""
from __future__ import annotations
from io import BytesIO
import json
import threading
import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle, PathPatch
from matplotlib.lines import Line2D
from matplotlib.path import Path
from matplotlib.ticker import FuncFormatter, MaxNLocator
from .ideogram import draw_ideogram, region_band_label

_RENDER_LOCK=threading.RLock()

def _pack_genes(genes,region,max_rows=12):
    ends=[]; packed=[]; omitted=0
    for g in sorted(genes,key=lambda g:(g['tx_start'],g['tx_end'],g['symbol'])):
        start=max(region.start,g['tx_start']); end=min(region.end,g['tx_end'])
        if end<=start: continue
        # Reserve label width in genomic units; conservative for a 12-inch panel.
        reserved=max(end,start+len(g['symbol'])*region.size*0.009)
        row=next((i for i,e in enumerate(ends) if e+region.size*.012<start),len(ends))
        if row>=max_rows: omitted+=1; continue
        if row==len(ends): ends.append(reserved)
        else: ends[row]=reserved
        packed.append((g,row,start,end))
    return packed,max(1,len(ends)),omitted

def _label(ax,text,color='#344054',size=10):
    ax.text(-.025,.5,text,transform=ax.transAxes,ha='right',va='center',fontsize=size,color=color,linespacing=1.15)

FIG_WIDTH=11.0           # inches; ~85-100 px/inch in a wide browser pane, so fonts below are sized for that
LEFT,RIGHT=.18,.90       # horizontal axes extent as figure fractions

def build_figure(result,genes,signals,contact=None,*,title='',shared_scale=True,headroom=1.5,
                 signal_summary='max',max_gene_rows=12,highlight_query=True,
                 cytobands=None,track_height=.42,ideogram_height=.55,ideogram_width=1.0,loops_height=.8,loops_width=1.0):
    """Caller must serialize rendering if used concurrently. See render_exports().

    Panel heights are real inches (no gridspec shrinkage), so track_height means what it says.
    `cytobands` is the list of (start,end,name,stain) for this chromosome, or None for no ideogram.
    """
    region=result.region
    if not result.can_plot: raise ValueError('Expanded region exceeds the 2 Mb limit.')
    packed,rows,omitted=_pack_genes(genes,region,max_gene_rows)
    gap=.12
    panels=[max(.9,.3*rows+.4)]+[float(track_height)]*len(signals)+[float(loops_height)]
    if contact is not None: panels.append(2.2)
    ideo_h=float(ideogram_height) if cytobands else 0.
    # Vertical layout in inches from the top: title block, ideogram, room for scale bar, panels.
    cursor=.95; ideo_top=cursor
    if ideo_h: cursor+=ideo_h+.1
    cursor+=.32
    tops=[]
    for h in panels: tops.append(cursor); cursor+=h+gap
    figure_height=cursor-gap+.95
    fig=Figure(figsize=(FIG_WIDTH,figure_height),facecolor='white')
    FigureCanvasAgg(fig)
    width=RIGHT-LEFT
    def add(top,h): return fig.add_axes([LEFT,1-(top+h)/figure_height,width,h/figure_height])
    axes=[add(t,h) for t,h in zip(tops,panels)]
    for ax in axes:
        ax.set_xlim(region.start,region.end)
        ax.spines[['left','right','top','bottom']].set_visible(False)
        ax.tick_params(axis='x',bottom=False,labelbottom=False)
        ax.set_yticks([])
        if highlight_query:
            ax.axvspan(result.seed.start,result.seed.end,color='#aab8cc',alpha=.14,zorder=-1)
    fig.text(LEFT,1-.22/figure_height,title or 'Micro-C regulatory landscape',fontsize=15,weight='bold',color='#182230',va='top')
    fig.text(LEFT,1-.58/figure_height,f'hg38  |  {region.label}  |  {region.size/1e3:,.0f} kb  |  canonical promoters',fontsize=9.5,color='#667085',va='top')
    if ideo_h:
        iw=width*min(max(float(ideogram_width),.2),1.)   # fraction of the plot width, left-aligned
        iax=fig.add_axes([LEFT,1-(ideo_top+ideo_h)/figure_height,iw,ideo_h/figure_height])
        draw_ideogram(iax,cytobands,region,width_in=FIG_WIDTH*iw)
        _label(iax,f'chr{region.chrom}\n{region_band_label(cytobands,region)}',size=10)
    # A scale no larger than 20% of the displayed genomic span, entirely inside the panel.
    target=region.size*.2
    magnitude=10**int(np.floor(np.log10(max(1,target))))
    scale=max(v*magnitude for v in (.1,.2,.5,1,2,5,10) if v*magnitude<=target)
    bar_start=region.end-region.size*.025-scale
    gene_h=panels[0]
    y_bar=1+.12/gene_h; y_txt=1+.17/gene_h
    axes[0].plot([bar_start,bar_start+scale],[y_bar,y_bar],transform=axes[0].get_xaxis_transform(),color='#344054',lw=1.5,clip_on=False)
    label=f'{scale/1e6:g} Mb' if scale>=1e6 else f'{scale/1e3:g} kb'
    axes[0].text(bar_start+scale/2,y_txt,label,transform=axes[0].get_xaxis_transform(),ha='center',fontsize=8.5)
    ax=axes[0]
    _label(ax,'Genes')
    for g,row,start,end in packed:
        y=rows-row
        color='#b42318' if g['symbol']==result.gene else '#475b8c'
        ax.plot([start,end],[y,y],color=color,lw=1)
        xs=np.arange(start+region.size*.009,end,region.size*.023)
        if len(xs): ax.plot(xs,np.full(len(xs),y),ls='none',marker='>' if g['strand']=='+' else '<',ms=3.5,color=color)
        blocks=[]
        for field,height in [('exons',.18),('cds',.32)]:
            for a,b in json.loads(g[field]):
                a,b=max(a,region.start),min(b,region.end)
                if b>a: blocks.append(Rectangle((a,y-height/2),b-a,height))
        if blocks: ax.add_collection(PatchCollection(blocks,facecolor=color,edgecolor='none'))
        # BED6 is a body line with arrows, never invented exons.
        ax.text(start,y+.22,g['symbol'],fontsize=8.5,color=color,clip_on=True)
    ax.set_ylim(.45,rows+.7)
    if omitted:
        ax.text(0,1+.08/gene_h,f'{omitted} crowded transcript labels/models omitted; gene table contains all genes',
                transform=ax.transAxes,ha='left',fontsize=7.5,color='#667085')
    ranges={}
    for index,s in enumerate(signals):
        vals=s.values[np.isfinite(s.values)]
        low=min(0,float(vals.min())) if len(vals) else 0
        high=max(0,float(vals.max())) if len(vals) else 0
        key=s.group if shared_scale else index
        a,b=ranges.get(key,(0,0)); ranges[key]=(min(a,low),max(b,high))
    for i,s in enumerate(signals,1):
        ax=axes[i]
        low,high=ranges[s.group if shared_scale else i-1]
        low=low*headroom; high=high*headroom
        if low==high: high=1
        if s.min_value is not None and s.max_value is not None:
            low,high=float(s.min_value),float(s.max_value)
            clipped=np.any(np.isfinite(s.values)&((s.values<low)|(s.values>high)))
            ax.text(1,.94,'Fixed scale'+(' · signal clipped' if clipped else ''),transform=ax.transAxes,
                    ha='right',va='top',fontsize=7.5,color='#667085')
        # Step traces preserve a common genomic x-axis and bin boundaries.
        xx=np.repeat(s.edges,2)[1:-1]; yy=np.repeat(s.values,2)
        ax.fill_between(xx,0,yy,color=s.color,alpha=.9,lw=0)
        ax.axhline(0,color='#d0d5dd',lw=.7)
        ax.set_ylim(low,high)
        # IGV-style range tag inside the track: tick labels would overlap the neighbouring track.
        ax.text(.004,.94,f'[{low:.3g}–{high:.3g}]',transform=ax.transAxes,ha='left',va='top',fontsize=7.5,
                color='#667085',bbox=dict(facecolor='white',alpha=.75,edgecolor='none',pad=1))
        _label(ax,s.label,s.color,size=9.5)
        if not np.isfinite(s.values).any():
            ax.text(.5,.5,'No signal coverage',transform=ax.transAxes,ha='center',va='center',fontsize=9,color='#98a2b3')
    ax=axes[1+len(signals)]
    if float(loops_width)!=1: raise ValueError('loops_width must be 1 to preserve genomic alignment.')
    _label(ax,'Promoter-linked\nloops')
    loops=result.loops
    cis=loops[(loops.chr==region.chrom)&(loops.targetChr==region.chrom)]
    seen=set()
    patches=[]
    arc_colors=[]
    for r in cis.itertuples():
        anchors=tuple(sorted(((r.start,r.end),(r.targetStart,r.targetEnd))))
        key=(anchors,r.canon_annot)
        if key in seen: continue
        seen.add(key)
        a,b=[(x+y)/2 for x,y in anchors]
        height=max(.06,((b-a)/region.size)**.65)
        path=Path([(a,0),(a,-height),(b,-height),(b,0)],
                  [Path.MOVETO,Path.CURVE4,Path.CURVE4,Path.CURVE4])
        patches.append(PathPatch(path))
        arc_colors.append('#4775a5' if r.canon_annot=='CRE' else '#d36b35')
    if patches: ax.add_collection(PatchCollection(patches,facecolor='none',edgecolor=arc_colors,linewidth=1,alpha=.75))
    ax.legend(handles=[Line2D([0],[0],color='#4775a5',label='CRE'),Line2D([0],[0],color='#d36b35',label='noCRE')],loc='lower right',fontsize=7,frameon=False,ncol=2)
    ax.axhline(0,color='#eaecf0',lw=.8)
    ax.set_ylim(-.82,.04)
    if not patches: ax.text(.5,.5,'No matching loops',transform=ax.transAxes,ha='center',fontsize=9.5,color='#98a2b3')
    if contact is not None:
        ax=axes[-1]
        edges=contact.edges.astype(float)
        row,col=np.meshgrid(edges,edges,indexing='ij')
        # Rotation maps bin pairs to their genomic midpoint and separation / 2.
        x=(row+col)/2; y=(col-row)/2
        z=contact.matrix.copy()
        z[np.tril_indices(len(z),k=-1)]=np.nan
        if np.nanmin(z)<0: raise ValueError('Negative contact counts cannot be log1p-transformed.')
        z=np.log1p(z)
        off=z[np.triu_indices(len(z),k=1)]
        finite=off[np.isfinite(off)&(off>0)]
        vmax=float(np.percentile(finite,98)) if len(finite) else 1
        mesh=ax.pcolormesh(x,y,np.ma.masked_invalid(z),cmap='Reds',vmin=0,vmax=vmax,shading='flat',rasterized=True)
        ax.set_ylim(0,region.size/2)
        _label(ax,contact.label+'\n'+contact.normalization+f'\n{contact.resolution/1000:g} kb')
        cbax=ax.inset_axes([1.025,.08,.018,.7])
        cb=fig.colorbar(mesh,cax=cbax)
        cb.ax.tick_params(labelsize=8,length=2)
        cb.set_label('log1p(contact)',fontsize=8)
    axes[-1].spines['bottom'].set_visible(True)
    axes[-1].spines['bottom'].set_color('#d0d5dd')
    axes[-1].tick_params(axis='x',bottom=True,labelbottom=True,labelsize=9,color='#d0d5dd')
    axes[-1].xaxis.set_major_locator(MaxNLocator(6))
    axes[-1].xaxis.set_major_formatter(FuncFormatter(lambda x,p:f'{x/1e6:.3f}'))
    axes[-1].set_xlabel(f'Chromosome {region.chrom} position (Mb)',fontsize=10,color='#475467')
    return fig

PNG_DPI=300
FORMATS=('svg','png','pdf')

def render_exports(result,genes,signals,contact=None,*,formats=FORMATS,**kwargs):
    """Render only the requested formats (the figure is built once).

    Matplotlib font/rendering state is not thread-safe, so rendering is serialized; query and
    file reads stay outside this lock. SVG is the sharp on-screen format (text becomes paths, so
    it looks identical in every browser). PNG is 300 dpi; its zlib level is lowered because the
    default level costs ~40% more time for ~20% smaller files, which does not pay off for web use.
    """
    bad=[f for f in formats if f not in FORMATS]
    if bad: raise ValueError(f'Unsupported export format: {bad}')
    with _RENDER_LOCK:
        fig=None
        try:
            fig=build_figure(result,genes,signals,contact,**kwargs)
            outputs={}
            for fmt in formats:
                buf=BytesIO()
                opts={'dpi':PNG_DPI,'pil_kwargs':{'compress_level':3}} if fmt=='png' else {}
                fig.savefig(buf,format=fmt,facecolor='white',**opts)
                outputs[fmt]=buf.getvalue()
            return outputs
        finally:
            if fig is not None: fig.clear()

def interactive_payload(result,genes,signals):
    """Extension point: JSON-ready tracks for a future browser renderer; no UI dependency."""
    return dict(region=dict(chrom=result.region.chrom,start=result.region.start,end=result.region.end),
                genes=genes,loops=result.loops.drop(columns=['id'],errors='ignore').to_dict('records'),
                signals=[dict(label=s.label,color=s.color,edges=s.edges.tolist(),
                              values=[float(v) if np.isfinite(v) else None for v in s.values]) for s in signals])
