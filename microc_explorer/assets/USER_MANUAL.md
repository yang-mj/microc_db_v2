# Micro-C Explorer user manual

This manual explains how to explore the administrator's curated hg38 datasets.
For installation, system requirements and server configuration, choose
**Installation guide (README)** in the Documentation selector.

## Quick start

1. Select **Explore** in the sidebar.
2. Select a **Dataset** and choose **Gene** or **Coordinates**.
3. Enter a symbol such as `KRAS`, a GENCODE ENSG gene ID, or an interval such as
   `chr1:25,500,000-25,850,000`.
4. Optionally select **Confirm gene / interval** to check the input before loading tracks.
5. Select the signal tracks and optional contact map, then click **Plot region**.
6. Inspect the notices, figure and result tables. Download the files you need.

The example query only returns gene results if that symbol is in your prepared
annotation. You can deselect all signal tracks to examine genes and loops alone.
At least two seconds must pass between confirmation/plot submissions.

## Choose a dataset

**All cases · merged 1 kb** uses `merged_1k`. Subtype datasets use their configured
10 kb source identifiers. Only sources actually present in the prepared database
appear in the menu. All queries include gene-associated `CRE` and `noCRE` records; changing the
source can therefore change the loop count and expanded plotting interval.

A source missing from the menu does not demonstrate that a biological interaction
is absent. Ask the administrator whether the complete table has been deployed.
An example-data banner indicates an incomplete preparation marked with `--sample`.

## Search by gene

Select **Gene** and enter an exact gene symbol, for example `KRAS`. Matching is
case-insensitive. ENSG IDs are also accepted, with or without a numeric version suffix (for example, `ENSG00000000003.15` is searched as `ENSG00000000003`). Gene symbols and IDs are resolved against the prepared protein-coding annotation; aliases,
free-text descriptions and fuzzy matches are not supported.

The input region is the annotated gene span. Loops are selected from canonical
gene associations identified by normalized GENCODE gene IDs in the loop table, which need not coincide exactly with the
representative transcript drawn in the gene track. If the symbol maps to multiple
chromosomes, search by coordinates. A gene can still be plotted when it has no
matching gene-associated loops in the chosen dataset.

## Search by coordinates

Select **Coordinates** and enter `chr1:25,500,000-25,850,000` or
`1:25500000-25850000`. Commas are optional. Supported primary chromosomes are
1–22, X, Y and mitochondrial M/MT.

Coordinates must use **GRCh38/hg38, 0-based, end-exclusive** intervals. The start
is included and the end is excluded. For one base, use `chr1:100000-100001`.
Coordinates copied from another browser may use a different convention; convert
1-based closed starts by subtracting one while keeping the end unchanged.

A loop matches if either anchor overlaps the original query interval. Merely
touching an interval endpoint is not overlap. Coordinates must lie within the
specified chromosome and the start must be smaller than the end.

## Confirm and plot

**Confirm gene / interval** reports the resolved region, gene-associated loop count and
expanded view. It does not read BigWigs/contact matrices or render a figure. It
retains matching loop records, so **Download loop results (CSV)** is available
after successful confirmation, including views exceeding 2 Mb.

**Plot region** performs the query and reads the selected tracks. Confirmation
is optional. Editing sidebar settings takes effect when you submit the form;
the displayed result remains labeled with its original submitted query.
A new confirmation or plotting request replaces the previous result.

## Automatic region expansion and the 2 Mb rule

The app expands the original gene/coordinate interval once to include directly
selected same-chromosome loop anchors, then adds **200 kb on each side**. It clips
the interval at chromosome boundaries. It does not recursively follow new links.

If the final interval exceeds **2,000,000 bp**, images are not generated. Accepted
matching loop records remain available as CSV; distant links are not silently
removed to fit the limit. Searching a smaller input does not necessarily produce
a smaller expanded interval if the same distant loop is selected.

Cross-chromosome/off-locus records remain in the loop table but cannot be drawn or
used to expand a single-chromosome figure. A notice reports their presence.

## Signal tracks and plot settings

| Control | Interpretation |
| --- | --- |
| Signal tracks | Select curated BigWigs to display. Deselect all for genes/loops only. |
| Display bins | 600–3,600 bins summarize the displayed interval; default 1,800. More bins preserve finer displayed detail but may cost more I/O/rendering. |
| Signal per bin: max | Maximum covered signal in each bin. Useful for displaying narrow peaks. |
| Signal per bin: mean | Mean over covered bases; missing bases are not filled with zero. |
| Exact BigWig summaries | Enabled by default. Disabling it uses available zoom summaries, which can be faster but may yield different values. |
| Match y-axis scales within each assay | Tracks sharing a configured assay group use a common automatic scale. Disable for independent track scales. |
| Highlight searched gene or interval | Shades the original query interval. The full plotted interval usually extends beyond this highlight. |
| Hi-C resolution | `0` selects an automatic available resolution; otherwise enter an available positive bin size in bp. |

Track heights, colors, order and optional fixed limits are controlled by the
administrator. Fixed limits override automatic scaling for that track and are
identified in the figure. Signal outside fixed limits is labeled as clipped.
Changing display settings does not renormalize the underlying BigWig values.

## Read the figure

**Chromosome ideogram:** the top strip displays the whole chromosome. The red
marker and pointer locate the expanded plotted interval. Cytoband labels are
shown where they fit; this overview has a different scale from the regional tracks.

**Genes:** lines show transcript bodies, arrows show strand direction, and exon/CDS
blocks appear when the database was prepared from GTF. BED-only annotations have
body lines/arrows. The searched gene is highlighted in red. Representative GTF
transcripts prioritize MANE Select, Ensembl canonical, coding transcript type and
then genomic span. Crowded models can be omitted after 12 rows; the gene table
retains all overlapping genes returned by the annotation query.

**Signals:** each colored panel represents a selected BigWig. Its range is shown
inside the panel. Gaps are missing coverage, not measured zero. Comparing heights
requires attention to the scale tags and the source data's normalization.

**Loops:** blue arcs represent `CRE` loops and orange arcs represent `noCRE` loops.
Both are promoter-associated chromatin loops. `noCRE` means a loop anchor contains
a gene promoter but the available epigenetic data do not support its CRE annotation.
This may reflect lower loop confidence or the need for additional epigenetic data;
it does not demonstrate absence of regulatory function. Both types are included
in gene searches, coordinate searches and automatic region expansion.

Loops without an assigned gene are not displayed. Support for non-gene-associated
chromatin loops is planned for the next version of the web server. Rows whose gene
IDs are absent from the prepared protein-coding GENCODE annotation are excluded
when the database is built.

Arcs join the centers of linked anchors. Arc depth is a drawing choice
based on genomic distance, not a measurement of interaction strength. Duplicate
geometric arcs within the same annotation category may be drawn once even when distinct annotations remain in CSV.

**Contact map:** the optional triangle shares the regional genomic axis. Color
represents `log1p` of normalized contact values, using an off-diagonal percentile
for the display color limit. It is not a direct linear scale of raw read counts.

All regional gene/signal/loop/contact panels align horizontally to one genomic
coordinate axis. These curated associations are not new enhancer predictions or
proof of a causal regulatory effect.

## Optional Hi-C contact maps

Contact maps appear only when the administrator configures files and installs the
optional dependencies. Select a map or leave **None** selected.

`.cool`/`.mcool` files use stored ICE weights, usually a column named `weight`.
`.hic` files use the configured stored normalization, normally `SCALE`.
Missing normalization does not trigger a silent fallback to raw counts.
The automatic resolution chooses the finest available size that stays within
600 bins per side; too-fine manual resolutions are rejected. Ask the administrator
if the requested normalization/resolution is unavailable.

## Downloads and tables

| Download | Contents and use |
| --- | --- |
| Loop CSV | All accepted matching records, including records that cannot be drawn. Available after successful confirmation or plotting. |
| SVG | Vector genes, labels, signals and arcs; suitable for enlargement and figure editing. Contact heatmaps are rasterized. |
| PDF | Vector figure suitable for publication workflows; contact heatmaps are rasterized. |
| PNG | Raster figure at 300 dpi, convenient for slides and image software. |
| Gene CSV | Overlapping genes and representative transcript information loaded for a plotted view. |

The web preview normally uses SVG for sharp text and curves. Very large SVGs
use a labeled PNG preview; the SVG file is still downloadable. Open **Loop table**
and **Genes and representative transcripts** to inspect tabular results.
The gene table/download is not populated by confirmation alone or a blocked plot.

Use the image fullscreen control or download SVG/PDF when examining fine details.
Keep the dataset/source, query, genome build, annotation and plot settings with
figures used in publications. Pixel width and display-bin count limit the detail
of a signal summary even when its drawing is vector.

## Sessions, storage and cleanup

The web app does not create a growing folder of plot files. It keeps the current
exports in memory, with at most one set per session and a **100 MB total export
budget per server process**. Least recently accessed sets can be evicted to admit
a new set. The limit concerns retained export bytes, not total server memory.

Connected sessions refresh their current exports every 30 seconds without
replotting. After 120 seconds without access, exports expire; a background worker
checks every 15 seconds. Closing the browser stops the refresh. Suspended tabs
can also expire. **Clear current result** releases the current application-held
exports immediately. A released plot can be regenerated; retained query CSV
results stay available. Files already downloaded to your laptop are unaffected.

Opening the manual keeps the current plot's refresh active while you stay connected.
Results are not a permanent analysis history. Download files you want to keep.

## Common messages

| Message or symptom | Meaning and next step |
| --- | --- |
| Gene not present | Verify the exact prepared symbol or use coordinates. |
| No matching loops | No gene-associated match in this source/query; genes and signal can still be shown. |
| Expanded region exceeds 2 Mb | Download loop CSV; use another query if a drawable view is needed. |
| More than 3,000 loops | Images are skipped to limit rendering cost; accepted loop CSV remains available. |
| More than 50,000 result records | The whole query is rejected; narrow it. No partial CSV is returned. |
| Server processing other plots | Retry when a slot becomes available; the default is two compute requests per process. |
| Plot released from memory | Select Plot region again. |
| No signal coverage | The file has no measured coverage in those bins; it is not a measured zero. |
| Track or matrix could not be plotted | Other valid tracks can remain; give the error reference to the administrator. |
| Dataset/configuration changed | Resubmit the query so results match the new data. |
| Cannot load configured dataset | An administrator must check paths/permissions; the manual remains available. |

## Reproducibility and support

Use **About / dataset** for the study preprint and laboratory contact details.
Administrators can inspect preparation metadata, counts and input checksums in
the preparation summary and SQLite metadata table. Record the software version, dataset source,
query interval, selected tracks, bin count/statistic, exact-summary setting,
signal scales and contact normalization/resolution with your analysis.

Send the deployment administrator the query, selected source, error reference and
software version when reporting a problem. The browser does not expose detailed
server paths or credentials. See the installation guide for server-side diagnosis.


## Database version and search help

Version 1.4 requires a database rebuilt from the new gene-ID loop table.
Old databases cannot be upgraded by adding indexes alone. Ask the administrator
to prepare the complete updated table with GENCODE v50 protein-coding annotation.

The search field accepts up to 128 characters, described in its **?** tooltip.
Click **Plot region** immediately below the search field to submit; Enter does
not submit the form. Plot settings and track selections apply to that submission.
