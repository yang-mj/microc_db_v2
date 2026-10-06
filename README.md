# Micro-C Explorer

**Version 1.5.0** | Human GRCh38/hg38 | Linux | Streamlit

Micro-C Explorer is a web application for querying and visualizing curated
chromatin interactions in pediatric B-cell precursor acute lymphoblastic leukemia
(BCP-ALL). It combines gene-associated CRE and noCRE chromatin loops, gene models,
BigWig signals and optional normalized Hi-C contact maps in an aligned genomic
figure. Queries accept a gene symbol or genomic coordinates.

The application supports Linux pip installations, Docker, Singularity/Apptainer
and container-based Streamlit hosting. A command-line plotting interface is also
available. It does not require pyGenomeTracks or a separate database server.

## Contents

- [Features](#features)
- [System requirements](#system-requirements)
- [Install with pip](#install-with-pip)
- [Server and laptop access](#server-and-laptop-access)
- [Docker](#docker)
- [Singularity / Apptainer](#singularity--apptainer)
- [SciLifeLab Serve / other Streamlit hosting](#scilifelab-serve--other-streamlit-hosting)
- [Database preparation and optimization](#database-preparation-and-optimization)
- [Configuration reference](#configuration-reference)
- [Using the web application](#using-the-web-application)
- [Plot and query behavior](#plot-and-query-behavior)
- [Generated outputs, concurrency and security](#generated-outputs-concurrency-and-security)
- [Troubleshooting](#troubleshooting)
- [Upgrading](#upgrading)
- [Validation and limitations](#validation-and-limitations)
- [Attribution and support](#attribution-and-support)

## Features

- Indexed SQLite lookup of gene-associated loops using normalized GENCODE gene IDs.
- One-step expansion to directly linked anchors, with 200 kb flanks and a 2 Mb plot limit.
- Chromosome ideograms with a marker for the displayed region.
- Gene bodies, strand arrows and exon/CDS models when GTF annotation is available.
- Compact BigWig tracks with shared or independent signal scales.
- Optional triangular `.cool`/`.mcool` ICE and `.hic` SCALE contact maps.
- SVG, PDF, 300 dpi PNG, loop CSV and gene CSV downloads.
- A searchable **User manual** page, including this installation guide.
- Session-specific results, bounded in-memory exports and background cleanup.

The distribution includes runtime code, hg38 cytobands, configuration templates,
dependency constraints and deployment files. Scientific datasets are supplied
separately: no SQLite database, BigWig, contact matrix or GTF is included. The
supplied configuration contains six ER/HeH tracks for CTCF, H3K27ac and H3K4me3.
Additional curated tracks can be configured by the administrator.

## System requirements

### Software

| Component | Requirement or support scope |
| --- | --- |
| Operating system | Linux server/workstation; native Windows/macOS execution is outside the validated scope. Clients can use a web browser on any OS. |
| Python | Python 3.11 or later; Python 3.12 was used for validation and is the recommended deployment version. Other interpreter versions require dependency compatibility checks. |
| Browser | A current browser with JavaScript, WebSocket and SVG support. No Python installation is needed on client laptops. |
| Database | Readable prepared SQLite file; Python's built-in `sqlite3` module is used. No SQLite service or `pip install sqlite3` is needed. |
| Core Python libraries | Streamlit, NumPy, pandas, Matplotlib and pyBigWig; installed automatically by pip. Version ranges are in `pyproject.toml`; tested primary pins are in `constraints.txt`. |
| Optional contact libraries | Cooler for `.cool`/`.mcool`; hic-straw for `.hic`. Install the `contact` extra. |
| Container installation | Docker or Apptainer/Singularity on the host. Container users do not need a separate host Python environment. |
| GPU | Not required or used. |
| Network | Needed to install dependencies and retrieve external data. Runtime can be offline with local datasets; remote BigWigs require working HTTPS byte-range access. |

### Hardware and storage planning

These are **starting allocations, not measured minimum requirements or capacity
guarantees**. Measure memory and response times using representative datasets.

| Workload | Suggested starting allocation |
| --- | --- |
| Small deployment, mostly genes/loops/BigWigs | 2 CPU cores/vCPUs and 4 GB RAM |
| Multiple active users or optional contact maps | 4 CPU cores/vCPUs and 8 GB RAM; increase if measured usage requires it |
| Full-data database preparation or container builds | Prefer 8 GB RAM and additional temporary disk space; profile the actual input/build |
| Disk | Space for the Python environment/container image plus all scientific data and temporary build files; optimization additionally requires space for another database copy |

Contact matrices may dominate storage. Use actual dataset sizes rather than
estimating disk demand from the small application ZIP. A local SSD and read-only
local track mounts usually provide more predictable access than remote files.
The 100 MB plot-export budget is not a limit on total application RAM.

The default compute limit is two concurrent requests **per process**, with
serialized Matplotlib rendering. Increasing the limit does not imply a matching
increase in rendering throughput. Client hardware only needs to display the
resulting web page; computation runs on the server.

## Install with pip

Extract this ZIP into a new directory and enter it. Keep your production database
and tracks outside the extracted source tree.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -c constraints.txt .
```

Copy/adapt your existing `config.json`: paths are relative to that configuration
file. Set `database` to your full `explorer.sqlite` or optimized copy, and set
BigWig paths to local files or curated direct HTTPS URLs supporting byte ranges.
`cytoband: null` uses bundled UCSC hg38 bands; a local .txt/.tsv or .gz path overrides
those bands. The public UI does not accept arbitrary paths or remote URLs.

```bash
microc-web --config /absolute/path/config.json
```

Open http://localhost:8501 on that computer. `microc-web` also works from an
installed wheel without keeping the source folder. For a wheel-only installation:

```bash
python -m pip install /path/to/microc_explorer-1.5.0-py3-none-any.whl
microc-web --init-config config.json
# Edit config.json with the production database and track paths.
microc-web --config config.json
```

The template writer refuses to overwrite an existing file. If console commands
are not on PATH, use `python -m microc_explorer.launch` instead of `microc-web`.
Python's sqlite3 module is sufficient; no separate SQLite server is needed.

Optional contact readers:

```bash
python -m pip install -c constraints.txt '.[contact]'
```

For a wheel-only installation, install its optional extra using:

```bash
python -m pip install '/path/to/microc_explorer-1.5.0-py3-none-any.whl[contact]'
```

On Rocky Linux 9, building hic-straw may require `gcc-c++`, `libcurl-devel`,
`zlib-devel` and the development headers matching your Python. On Debian/Ubuntu,
use `build-essential libcurl4-openssl-dev zlib1g-dev` and Python headers if needed.
These system packages may need administrator installation. With conda, create a
Python 3.12 environment and install this package through pip inside it:

```bash
conda create -n microc-explorer python=3.12 pip
conda activate microc-explorer
python -m pip install -c constraints.txt .
```

This is a pip installation in a conda environment; no published conda package is
assumed. Install the contact extra only if contact maps are required.

## Server and laptop access

For private use, start on the server bound to localhost:

```bash
microc-web --config /datasets/config.json --host 127.0.0.1 --port 8501
```

From your laptop, open an SSH tunnel:

```bash
ssh -N -L 8501:127.0.0.1:8501 USER@SERVER
```

Then browse http://localhost:8501 on the laptop. No remote Firefox or X11 is needed.
For a managed public service, bind to `0.0.0.0` behind the platform's HTTPS reverse
proxy with WebSocket support. Keep CORS/XSRF enabled. Run under your platform's
service manager or container restart policy so it survives logout/reboots.

A direct Streamlit launch from this directory is also supported:

```bash
MICROC_CONFIG=/datasets/config.json python -m streamlit run app.py \
  --server.address=0.0.0.0 --server.port=8501
```

Keep `.streamlit/config.toml` alongside this checkout. The `microc-web` launcher
passes the essential security/session settings explicitly, including when launched
outside the source directory.

## Docker

From the clean extracted directory:

```bash
docker build --platform linux/amd64 -t microc-explorer:1.5.0 .
docker run --rm -p 127.0.0.1:8501:8501 \
  -e MICROC_CONFIG=/datasets/config.json \
  -v /absolute/path/to/datasets:/datasets:ro \
  microc-explorer:1.5.0
```

The `datasets` directory should contain your config, database and tracks, or mount
additional directories and reference their container paths. They must be readable
by UID 1000. The image runs as a non-root user. Runtime data are not baked into it.
For the SSH-tunnel workflow, leave the host port bound to 127.0.0.1 as above.

Build optional contact readers into the image:

```bash
docker build --platform linux/amd64 --build-arg WITH_CONTACT=1 \
  -t microc-explorer:1.5.0-contact .
```

For local persistent operation add an appropriate restart policy and omit `--rm`.
Public access should go through the hosting service/reverse proxy, not an unprotected
internet-facing HTTP port.

## Singularity / Apptainer

Build **from the clean extracted package**, before copying large datasets into it:

```bash
apptainer build --fakeroot microc-explorer.sif Singularity.def
apptainer run --bind /absolute/path/to/datasets:/datasets:ro \
  --env MICROC_CONFIG=/datasets/config.json microc-explorer.sif
```

The definition installs optional contact readers. Use `singularity` instead of
`apptainer` where appropriate. `--fakeroot` depends on local configuration; use your
administrator's build workflow if unavailable. No port mapping is needed for
Apptainer; use the SSH tunnel above and follow your HPC site's web-server rules.
An existing Docker image can alternatively be converted into a SIF.

## SciLifeLab Serve / other Streamlit hosting

Build and push the Docker image to your chosen registry, then configure the host
to run it on port 8501. The image contains `/app/app.py` for platforms that require
that entry point. Mount your persistent data directory and set
`MICROC_CONFIG=/datasets/config.json` (adjust the mount path to match your host).
A local Docker bind mount does not transfer your files to the hosting platform.
For a pip-based host, install the wheel or `pip install .` and start `microc-web`
with `--host 0.0.0.0 --port 8501`, or use the direct Streamlit command above.

Reference: https://serve.scilifelab.se/docs/application-hosting/streamlit/

## Database preparation and optimization

Use your existing full prepared database. If you need to build one from scratch:

```bash
microc-prepare --loops /datasets/concat_loops_v2.tab \
  --gtf /datasets/gencode.v50.primary_assembly.annotation.gtf.gz \
  --output /datasets/explorer.sqlite
```

GENCODE GTF is recommended and gives exon/CDS models. The optional `--genes-bed`
input now requires a curated protein-coding BED7: chromosome, start, end, symbol,
score, strand, ENSG gene ID. A previous BED6 without gene IDs is insufficient. The GENCODE URL is:
https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.primary_assembly.annotation.gtf.gz

Input loop columns (tab-separated header):
`chr start end targetChr targetStart targetEnd targetGeneID_canon targetGene_canon canon_annot loopSource`.
Loop/BED coordinates must be hg38 0-based, end-exclusive. GTF 1-based closed
coordinates are converted internally. Retained sources are `merged_1k` and
`HeH_10k`, `ER_10k`, `BA_10k`, `DUX4r_10k`, `TP_10k`, `KMT2Ar_10k`, `iAMP_10k`,
`nearHaploid_10k`; plots include gene-associated `CRE` and `noCRE` records. Use `--sample` when
preparing incomplete example input. Preparation refuses existing output unless
explicitly passed `--force`; never replace a live database while serving requests.

Inspect/optimize an existing database offline:

```bash
microc-index --database /datasets/explorer.sqlite
microc-index --database /datasets/explorer.sqlite \
  --output /datasets/explorer.optimized.sqlite
```

Point the config to the optimized copy and restart. Optimization does not alter
the source and refuses an existing output. Existing version 1 databases work
without optimization; version 1.1+ databases use the bounded coordinate indexes.

## Configuration reference

The administrator edits `config.json`; browser users select only configured datasets.
Relative paths are resolved from the directory containing the configuration file.
Use absolute paths inside containers when data are mounted externally.

| Top-level key | Meaning |
| --- | --- |
| `title` | Title displayed in exported figures. |
| `database` | Path to the prepared SQLite database. Required. |
| `cytoband` | `null` for bundled hg38 bands, or a local tab-separated plain/gzipped cytoband file. |
| `headroom` | Automatic signal-scale multiplier, allowed range 1–10; supplied value 1.5. |
| `track_height` | Signal-panel height in inches, allowed range 0.2–3; supplied value 0.42. |
| `loops_height` | Loop-panel height in inches, allowed range 0.2–4; supplied value 0.5. |
| `loops_width` | Must be 1 to preserve horizontal genomic alignment. |
| `ideogram_height` | Ideogram height in inches, allowed range 0.2–1.5; supplied value 0.35. |
| `ideogram_width` | Fraction of plot width, allowed range 0.2–1; supplied value 1. |
| `bigwigs` | Array of up to 16 curated signal entries. An empty array permits genes/loops-only plots. |
| `contacts` | Array of configured local contact maps; an empty array disables this data option. |

A signal entry has this form (within the `bigwigs` array):

```json
{
  "id": "CTCF_ER",
  "label": "CTCF · ER",
  "group": "CTCF",
  "color": "#4775a5",
  "path": "/datasets/tracks/CTCF_ER.perCopy.bw"
}
```

IDs must be unique within each track type and contain only letters, numbers,
underscores or hyphens. `group` determines shared assay scales. Optional
`min_value` and `max_value` must both be present, finite and increasing; they
replace automatic limits for that track. Clipped signal is indicated on the plot.
Remote BigWigs must be trusted direct HTTPS URLs without embedded user/password
credentials. Contact matrices and cytoband overrides must be local files.

| Environment variable | Default | Meaning |
| --- | --- | --- |
| `MICROC_CONFIG` | `config.json` in the launch directory | Configuration path. Explicit launcher `--config` takes precedence. |
| `MICROC_MAX_CONCURRENT` | `2` | Active compute requests per process; integer 1–8. |
| `MICROC_EXPORT_CACHE_MB` | `100` | Retained export budget per process in decimal MB; integer 1–100. |
| `MPLCONFIGDIR` | Matplotlib default | Writable Matplotlib cache directory; containers set a directory under `/tmp`. |

## Using the web application

1. Open the app URL and select **Explore** in the sidebar.
2. Choose a dataset. Only sources present in the prepared database are offered.
3. Select **Gene** or **Coordinates**. Enter a gene symbol such as `KRAS`, or
   a region such as `chr1:25,500,000-25,850,000`.
4. Optionally click **Confirm** (confirm gene / interval) to inspect the resolved location,
   loop count and expanded span. This reads no tracks and already allows loop CSV download.
5. Choose signal tracks and, if configured, a contact map. Adjust **Plot settings**.
6. Press **Enter** in the search box or click **Plot region**. Allow at least two seconds between submitted requests.
   Choose SVG or high-resolution PNG under **Plot settings → Figure format**.
7. Read notices, inspect the figure and tables, then download the desired files.

Input intervals use **hg38, 0-based, end-exclusive** coordinates. A single-base
interval is written as `chr1:100000-100001`, not `chr1:100000-100000`.
Gene search uses exact prepared symbols, case-insensitively; aliases and fuzzy
matches are not resolved. Available results depend on the curated source data.

Open **User manual** for a searchable explanation of each control, scientific
interpretation, exports and common messages. The manual and README can be downloaded
from that page. It is also accessible before a database has been configured.
**About / dataset** shows the loaded dataset's preparation metadata.

## Plot and query behavior

- The confirmation button validates the gene/interval and retains matching CSV
  results without reading tracks. The plot button adds tracks and images. Both
  share the concurrency limit and per-session two-second click cooldown.
- Expand once to directly linked same-chromosome anchors, then add 200 kb on each
  side. A final span above 2 Mb stops plotting while preserving accepted CSV results.
- `track_height`, `loops_height`, and `ideogram_height` specify physical panel
  heights in inches. `ideogram_width` controls the chromosome overview width.
  `loops_width` must remain 1 so loop anchors stay aligned with signal/gene positions.
- Gene/loop/signal panels share a genomic axis. The ideogram has its own full
  chromosome scale and marks the exact expanded interval, with a pointer for tiny
  windows. Labels appear only where band names fit.
- SVG is the sharp browser preview. Very large SVGs (>8 MB) use a labeled PNG
  preview instead; SVG remains downloadable. PNG exports are 300 dpi. PDF/SVG
  contain vector labels/signals/loops; contact matrices are rasterized.
- Matching assay groups share signal scales by default. Fixed `min_value` and
  `max_value` together override a track's scale. Missing signal stays missing.
- Limits: 50,000 result records (larger queries explicitly rejected); 3,000 loops
  for rendering (larger accepted results keep CSV); 3,600 signal bins; 600 contact
  bins per side; 10-second SQLite execution budget.

Optional local contact config example:

```json
"contacts": [
  {"id":"ER_ice","label":"ER Micro-C","path":"/datasets/ER.mcool","format":"mcool","weight_column":"weight","divisive_weights":false},
  {"id":"HeH_scale","label":"HeH Hi-C","path":"/datasets/HeH.hic","format":"hic","normalization":"SCALE"}
]
```

ICE weights and SCALE normalization must already exist; no silent fallback to raw
counts. Explicit Cooler URIs such as `file.mcool::/resolutions/10000` are supported.
These optional libraries/data are unnecessary for BigWig-and-loop-only hosting.

## Generated outputs, concurrency and security

No web plot files are written to disk. Finished plots are kept in a **shared in-memory
cache**: because the datasets are curated and public, one visitor's plot also serves the next
visitor who submits the identical query, tracks and settings (seconds saved on popular genes).
The cache key is a hash of validated request parameters and the dataset file fingerprints; no
session data are stored. The cache is capped at 100 MB (decimal), evicting the least recently used
plot first; `MICROC_EXPORT_CACHE_MB` can lower the cap to 1–100 MB. A single oversized plot is
rejected; CSV remains. Entries expire after 15 minutes without access (swept every 30 seconds),
so memory is released when the server is idle. The cap excludes Streamlit media/download copies,
temporary rendering buffers, tables, numeric caches and Python overhead; use hosting memory
limits for total RSS. Files downloaded to a laptop are never deleted.

Only the format chosen under **Plot settings → Figure format** (SVG or 300 dpi PNG) is rendered
when you plot. The other formats (SVG, PNG, PDF) are produced when you click their download
button, from the cached numeric inputs and without re-reading tracks or the database.

Uncached plots and input checks are limited to `MICROC_RATE_PER_MINUTE` per minute for the whole
process (default 240; `0` disables). This protects the server from floods that open many
browser sessions, which the per-session two-second cooldown cannot do. **Per-visitor** limits,
TLS, security headers (CSP, X-Frame-Options, HSTS) and authentication for private data belong
at the hosting layer/reverse proxy: Streamlit does not let an application set response headers.

SQLite connections are read-only and per request; query values are parameterized.
Native file handles are per call. Matplotlib exports and native hic reads are
serialized. The export store uses a lock and random server-generated owner IDs.
Two active compute requests per process are allowed by default; set
`MICROC_MAX_CONCURRENT=1` through `8` after measuring capacity. Extra requests get
a retry message. Rate/connection limits, HTTPS and private-data authentication
belong at the hosting layer. Global caches assume curated **public** data.

CORS/XSRF protections stay enabled. Public error messages omit diagnostic paths;
admin logs may include them and should be restricted. Curated HTTPS URLs are
administrator-controlled; this is not a public URL-fetch service. Remote native
I/O can stall independently of the SQLite timeout. Use process supervision,
resource limits and trusted local read-only data for predictable service.

Optional CLI: `microc-plot --config config.json --gene KRAS --output results/KRAS`.
CLI files are deliberately written to that prefix and are not auto-deleted; use
different prefixes for simultaneous jobs. Exit status 2 means no plot was produced.

## Troubleshooting

| Symptom | Action |
| --- | --- |
| Cannot load the configured dataset | Check `MICROC_CONFIG`, the database path, mount location and read permissions. The manual is still available. |
| Gene is not present | Check the prepared symbol/annotation and assembly. Use coordinates if the symbol is an alias or has multiple chromosome loci. |
| No supported loop sources | Check `loopSource` values against the supported identifiers listed above. |
| Expanded region exceeds 2 Mb | Download accepted loop CSV results. Select another query if a plot is needed; no links are dropped to force a smaller view. |
| Plot released from server memory | Regenerate it. Export eviction/expiry does not invalidate retained query CSV results. |
| Server busy | Retry after an active request finishes; administrators should measure capacity before increasing concurrency. |
| Missing signal/contact track | Check the admin log reference, track path, chromosome names, assembly, and stored normalization. Other valid tracks can still render. |
| Remote BigWig does not open | Confirm the endpoint serves the file with byte-range support and `python -c "import pyBigWig; print(pyBigWig.remote)"` prints `1`. Prefer local mounting when remote access is unreliable. |
| Hi-C dependency cannot import | Install the contact extra in the same Python environment; confirm compiler/libcurl/zlib dependencies for hic-straw builds. |
| Console command not found | Activate the environment or use the corresponding `python -m microc_explorer...` command. |
| Web page unavailable from laptop | Check server state and the SSH tunnel or reverse-proxy configuration; use the laptop's local browser. |
| Documentation seems out of date | Upgrade the installed package and restart Streamlit. The manual is bundled with that package version. |

## Upgrading

Stop the app, extract a release into a new directory and retain the production
config/database/tracks separately. Install the new source package or wheel into
the intended environment, then restart using the existing config. For example:

```bash
python -m pip install --upgrade -c constraints.txt .
microc-web --config /datasets/config.json
```

Versions 1.4.0 to 1.5.0 do not change the database schema (schema version 3) or
require rebuilding an existing optimized database. Rebuild the container image
when upgrading a container deployment. Keep the previous release for rollback.

## Validation and limitations

The 1.3 release passed 26 functional tests against the installed wheel, including
SQLite queries, concurrent requests, export cleanup, synthetic BigWig/Cooler reads,
and Streamlit confirmation/plot/download behavior. The installed launcher returned
HTTP 200 / `ok`, and a synthetic six-track layout was inspected visually.
The 1.4.0 update added packaged documentation and in-app manual navigation. On
2026-10-02, all 28 functional tests passed against its installed wheel, including
manual access without a configured database, keyword/topic filtering, and returning
to the current plot after visiting the manual.

Actual production datasets, real SCALE-normalized `.hic` files, container builds,
public hosting, and performance/security testing under production traffic remain
deployment validation tasks. Container runtimes were not available during package
validation. Hardware allocations above are planning suggestions, not benchmarked
minimum specifications. Pinned dependencies are not a vulnerability certification.

## Attribution and support

Gene annotations can be prepared from [GENCODE human release 50](https://www.gencodegenes.org/human/release_50.html).
Bundled cytobands come from [UCSC hg38 cytoBandIdeo](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBandIdeo.txt.gz);
the source checksum/date are recorded in `microc_explorer/assets/SOURCE.json`.
The supplied original application code credited Efe Aydın. No new software license,
publication DOI, institutional endorsement or support address is asserted here.
Confirm the original code/data redistribution terms before public publication.

When reporting a problem to the deployment administrator, include the application
version, selected dataset, query, error reference and relevant configuration details
with credentials or signed URLs removed. Cite the application version, genome build,
annotation source, scientific datasets and plot settings when reporting results.

Additional documentation:

- [Streamlit documentation](https://docs.streamlit.io/)
- [SciLifeLab Serve Streamlit hosting](https://serve.scilifelab.se/docs/application-hosting/streamlit/)
- [pyBigWig](https://github.com/deeptools/pyBigWig)
- [Cooler](https://cooler.readthedocs.io/en/latest/)
- [hic-straw](https://github.com/aidenlab/straw/tree/master/pybind11_python)


## Version 1.5: speed, security and interface

- BigWig signals are read once per region and binned with NumPy: values are identical to the
  previous exact summaries (verified value-for-value, including gapped tracks) at roughly 10–30× the speed.
- Only the selected figure format is rendered up front; others are made on demand.
- Shared plot cache, startup warm-up of Matplotlib's font cache, lazy tables, deferred downloads
  that no longer rerun the page.
- Enter in the search box submits **Plot region**; **Figure format** (SVG or 300 dpi PNG) is under Plot settings.
- Process-wide rate limiter, Markdown-escaped visitor text, `.streamlit/config.toml` with
  hardened server settings, and a cleaner layout. No database rebuild is needed.

## Version 1.4: gene-ID preparation and loop annotations

Rebuild `explorer.sqlite` from the full new-format loop table and GENCODE v50
GTF before using this version. Existing schema-2 databases are rejected with a
rebuild message; indexing alone cannot migrate their missing gene associations.
Prepare a new output filename offline, update `config.json`, and restart the app.

Each loop row must contain exactly one `targetGeneID_canon` ENSG identifier.
Numeric version suffixes are removed in both GTF and loop inputs. Only IDs
present in the prepared protein-coding annotation are retained; missing IDs,
non-protein-coding IDs and invalid/empty IDs are excluded without symbol fallback.
The displayed/exported symbol is resolved from GENCODE by ID, rather than trusted
from the input symbol column. Coordinates stay 0-based, end-exclusive.

Both `CRE` and `noCRE` are imported, indexed and queried. Preparation prints
`annotation_counts`, `excluded_gene_rows`, `unmatched_gene_ids` and duplicate
counts in its JSON summary. Input checksums remain stored in SQLite metadata
for administrator reproducibility, but are not displayed on the About page.

Search accepts exact gene symbols, unversioned ENSG IDs and versioned ENSG IDs.
CRE loops are blue; noCRE loops are orange. noCRE denotes a promoter-associated
loop lacking supporting epigenetic evidence in the available data, potentially
reflecting lower loop confidence or incomplete epigenetic coverage. It does not
establish that the loop has no regulatory function. Non-gene-associated loops
are not displayed; support is planned for the next version.

The About page links to the study preprint:
https://www.biorxiv.org/content/10.1101/2025.07.22.666073v1
