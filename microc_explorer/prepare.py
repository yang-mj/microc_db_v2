"""One-time conversion of loop TSV + BED6/GTF into a compact, indexed SQLite database."""
from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import sqlite3
import tempfile
from decimal import Decimal, InvalidOperation
from .core import chromosome, CHROMS, LOOP_COLUMNS, SOURCES, normalize_gene_id
from .index_db import install_indexes

GTF_URL = 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.primary_assembly.annotation.gtf.gz'

def integer(value):
    try: n = Decimal(str(value))
    except InvalidOperation as e: raise ValueError(f'Invalid coordinate: {value}') from e
    if not n.is_finite() or n != n.to_integral_value(): raise ValueError(f'Noninteger coordinate: {value}')
    return int(n)

def validate(chrom,start,end):
    if chrom not in CHROMS or not 0 <= start < end <= CHROMS[chrom]:
        raise ValueError(f'Invalid hg38 coordinates: {chrom}:{start}-{end}')

def checksum(path):
    h = hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(1024*1024), b''): h.update(block)
    return h.hexdigest()

def gtf_records(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path,'rt') as f:
        for line in f:
            if line.startswith('#'): continue
            c = line.rstrip('\n').split('\t')
            if len(c) != 9: raise ValueError('GTF requires nine tab-separated columns.')
            if c[2] not in ('gene','transcript','exon','CDS'): continue
            chrom = chromosome(c[0])
            if chrom not in CHROMS: continue
            attrs = {}
            for key,value in re.findall(r'(\S+)\s+"([^"]*)"',c[8]):
                attrs.setdefault(key,[]).append(value)
            yield chrom,c,attrs

def prepare_genes(con, bed=None, gtf=None):
    if gtf:
        # First pass: rank transcripts. GENCODE GTF is 1-based closed: subtract one from start only.
        genes, chosen = {}, {}
        for chrom,c,a in gtf_records(gtf):
            gid = normalize_gene_id(a.get('gene_id',[''])[0])
            if a.get('gene_type',[''])[0] != 'protein_coding': continue
            start,end = int(c[3])-1,int(c[4])
            validate(chrom,start,end)
            key = (chrom,gid)
            if c[2] == 'gene':
                genes[key] = dict(chrom=chrom,start=start,end=end,symbol=a.get('gene_name',[gid])[0].upper(),strand=c[6])
            elif c[2] == 'transcript':
                tid = a['transcript_id'][0]
                tags = a.get('tag',[])
                score = (int('MANE_Select' in tags), int('Ensembl_canonical' in tags),
                         int(a.get('transcript_type',[''])[0]=='protein_coding'), end-start,tid)
                if key not in chosen or score > chosen[key]['score']:
                    chosen[key] = dict(score=score,tid=tid,start=start,end=end,exons=[],cds=[])
        selected = {(chrom,v['tid']):v for (chrom,_),v in chosen.items()}
        for chrom,c,a in gtf_records(gtf):
            if c[2] not in ('exon','CDS'): continue
            t = selected.get((chrom,a.get('transcript_id',[''])[0]))
            if t is not None: t['exons' if c[2]=='exon' else 'cds'].append([int(c[3])-1,int(c[4])])
        for key,g in genes.items():
            t = chosen.get(key)
            con.execute('INSERT INTO genes VALUES (?,?,?,?,?,?,?,?,?,?,?)',
                (*g.values(),key[1],t['tid'] if t else '', t['start'] if t else g['start'],
                 t['end'] if t else g['end'],json.dumps(sorted(t['exons'])) if t else '[]',
                 json.dumps(sorted(t['cds'])) if t else '[]'))
    elif bed:
        with open(bed) as f:
            for line in f:
                if not line.strip() or line.startswith(('#','track','browser')): continue
                c = line.rstrip().split('\t')
                if len(c)<7: raise ValueError('Gene BED requires seven columns: chrom,start,end,symbol,score,strand,ENSG_gene_id. Use GENCODE GTF for standard annotations.')
                gid=normalize_gene_id(c[6])
                if not re.fullmatch(r'ENSG\d{11}',gid): raise ValueError('Gene BED column 7 must contain an ENSG gene ID.')
                chrom,start,end = chromosome(c[0]),integer(c[1]),integer(c[2])
                validate(chrom,start,end)
                if c[5] not in ('+','-'): raise ValueError('Gene strand must be + or -.')
                con.execute('INSERT INTO genes VALUES (?,?,?,?,?,?,?,?,?,?,?)',
                            (chrom,start,end,c[3].upper(),c[5],gid,'',start,end,'[]','[]'))
    else: raise ValueError('Provide --gtf or --genes-bed.')

def build_database(loops, output, bed=None, gtf=None, sample=False, force=False):
    output = Path(output).resolve()
    if any(output==Path(p).resolve() for p in (loops,bed,gtf) if p): raise ValueError('Output must not replace an input file.')
    if output.exists() and not force: raise FileExistsError('Output exists. Choose a new filename or explicitly use --force offline.')
    output.parent.mkdir(parents=True,exist_ok=True)
    fd,tmp = tempfile.mkstemp(prefix='.prepare-',suffix='.sqlite',dir=output.parent)
    os.close(fd)
    con = sqlite3.connect(tmp)
    try:
        con.executescript('''
          CREATE TABLE loops(id INTEGER PRIMARY KEY, chr TEXT,start INTEGER,end INTEGER,
            targetChr TEXT,targetStart INTEGER,targetEnd INTEGER,targetGeneID_canon TEXT,
            targetGene_canon TEXT,canon_annot TEXT,loopSource TEXT,
            UNIQUE(chr,start,end,targetChr,targetStart,targetEnd,targetGeneID_canon,targetGene_canon,canon_annot,loopSource));
          CREATE TABLE loop_genes(gene_id TEXT,loop_id INTEGER,PRIMARY KEY(gene_id,loop_id)) WITHOUT ROWID;
          CREATE TABLE genes(chrom TEXT,start INTEGER,end INTEGER,symbol TEXT,strand TEXT,gene_id TEXT,
            transcript_id TEXT,tx_start INTEGER,tx_end INTEGER,exons TEXT,cds TEXT);
          CREATE TABLE metadata(key TEXT PRIMARY KEY,value TEXT);
        ''')
        prepare_genes(con,bed,gtf)
        annotation = {}
        for gid,symbol in con.execute('SELECT gene_id,symbol FROM genes'):
            if gid in annotation and annotation[gid] != symbol:
                raise ValueError(f'Conflicting symbols for gene ID {gid} in annotation.')
            annotation[gid] = symbol
        n=kept=excluded=duplicates=0
        unmatched=set()
        counts={'CRE':0,'noCRE':0}
        with open(loops,newline='') as f:
            reader = csv.DictReader(f,delimiter='\t')
            if not set(LOOP_COLUMNS) <= set(reader.fieldnames or []):
                raise ValueError(f'Loop TSV must contain columns: {LOOP_COLUMNS}')
            for row in reader:
                n+=1
                if row['loopSource'].strip() not in SOURCES: continue
                values = [row[k] for k in LOOP_COLUMNS]
                values[0],values[3] = chromosome(values[0]),chromosome(values[3])
                for i in (1,2,4,5): values[i] = integer(values[i])
                validate(*values[:3]); validate(*values[3:6])
                for i in (6,7,8,9): values[i] = values[i].strip()
                gid = normalize_gene_id(values[6])
                if not re.fullmatch(r'ENSG\d{11}',gid) or gid not in annotation:
                    excluded+=1
                    unmatched.add(gid or '(empty)')
                    continue
                if values[8] not in ('CRE','noCRE'):
                    raise ValueError(f"Row {n}: canon_annot must be CRE or noCRE.")
                values[6],values[7] = gid,annotation[gid]
                cur = con.execute('INSERT OR IGNORE INTO loops('+','.join(LOOP_COLUMNS)+') VALUES ('+','.join('?'*10)+')',values)
                if not cur.rowcount:
                    duplicates+=1
                    continue
                kept+=1
                counts[values[8]]+=1
                con.execute('INSERT INTO loop_genes VALUES (?,?)',(gid,cur.lastrowid))
        con.executescript('''
          CREATE INDEX anchor1 ON loops(loopSource,canon_annot,chr,start,end);
          CREATE INDEX anchor2 ON loops(loopSource,canon_annot,targetChr,targetStart,targetEnd);
          CREATE INDEX gene_symbol ON genes(symbol);
          CREATE INDEX gene_region ON genes(chrom,start,end);
        ''')
        meta = dict(schema_version=3,input_rows=n,stored_rows=kept,sample_dataset=sample,
                    excluded_gene_rows=excluded,unmatched_gene_ids=sorted(unmatched),duplicate_rows=duplicates,annotation_counts=counts,
                    annotation='GENCODE GTF (protein-coding genes)' if gtf else 'Curated protein-coding BED7 gene bodies (no exons)',
                    loops_sha256=checksum(loops),annotation_sha256=checksum(gtf or bed),
                    available_sources=[r[0] for r in con.execute('SELECT DISTINCT loopSource FROM loops')])
        con.executemany('INSERT INTO metadata VALUES (?,?)',[(k,json.dumps(v)) for k,v in meta.items()])
        install_indexes(con)
        con.commit(); con.close()
        os.chmod(tmp,0o644)  # Curated public dataset; readable by the non-root container user.
        if force: os.replace(tmp,output)
        else: os.link(tmp,output); os.unlink(tmp)
        return meta
    except Exception:
        con.close()
        Path(tmp).unlink(missing_ok=True)
        raise

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--loops',required=True)
    p.add_argument('--output',required=True)
    p.add_argument('--genes-bed')
    p.add_argument('--gtf',help='Plain or gzipped GENCODE GTF; takes precedence over BED6.')
    p.add_argument('--sample',action='store_true',help='Mark incomplete/example data visibly in the app.')
    p.add_argument('--force',action='store_true',help='Explicitly replace an existing database offline.')
    a=p.parse_args()
    print(json.dumps(build_database(a.loops,a.output,a.genes_bed,a.gtf,a.sample,a.force),indent=2))
if __name__=='__main__': main()
