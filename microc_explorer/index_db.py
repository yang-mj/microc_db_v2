"""Inspect indexes or create an optimized COPY of an existing v1 database offline."""
import argparse
import json
import os
from pathlib import Path
import sqlite3
import tempfile

INDEX_SQL = (
    'CREATE INDEX IF NOT EXISTS anchor1 ON loops(loopSource,canon_annot,chr,start,end)',
    'CREATE INDEX IF NOT EXISTS anchor2 ON loops(loopSource,canon_annot,targetChr,targetStart,targetEnd)',
    'CREATE INDEX IF NOT EXISTS gene_id ON genes(gene_id)',
    'CREATE INDEX IF NOT EXISTS gene_symbol ON genes(symbol)',
    'CREATE INDEX IF NOT EXISTS gene_region ON genes(chrom,start,end)',
)

def install_indexes(con):
    """Run during offline preparation; callers commit. Input loop data is not altered."""
    for sql in INDEX_SQL: con.execute(sql)
    con.execute('''CREATE TABLE IF NOT EXISTS anchor_bounds(
        loopSource TEXT,chrom TEXT,anchor INTEGER,max_length INTEGER NOT NULL CHECK(max_length>0),
        PRIMARY KEY(loopSource,chrom,anchor)) WITHOUT ROWID''')
    con.execute('DELETE FROM anchor_bounds')
    for anchor,chrom,start,end in [(1,'chr','start','end'),(2,'targetChr','targetStart','targetEnd')]:
        # Identifiers are constants from this tuple, never user input.
        con.execute(f'''INSERT INTO anchor_bounds SELECT loopSource,{chrom},?,MAX({end}-{start})
            FROM loops WHERE canon_annot IN ('CRE','noCRE') GROUP BY loopSource,{chrom}''',(anchor,))
    # If an administrator later edits loops, invalidate the optimization rather
    # than use stale length bounds and accidentally omit matches. Queries fall back.
    for action in ('INSERT','UPDATE','DELETE'):
        con.execute(f'''CREATE TRIGGER IF NOT EXISTS invalidate_bounds_{action.lower()}
            AFTER {action} ON loops BEGIN DELETE FROM anchor_bounds; END''')
    # Preserve the database schema version; indexing must not migrate columns.
    con.execute("INSERT OR REPLACE INTO metadata VALUES ('index_strategy','\"bounded_start_range_v1\"')")
    con.execute('ANALYZE')
    con.execute('PRAGMA optimize')

def inspect_database(path):
    p=Path(path).resolve()
    con=sqlite3.connect(p.as_uri()+'?mode=ro',uri=True)
    try:
        indexes=[dict(name=n,sql=s) for n,s in con.execute("SELECT name,sql FROM sqlite_master WHERE type='index'")]
        bounded=bool(con.execute("SELECT 1 FROM sqlite_master WHERE name='anchor_bounds' AND type='table'").fetchone())
        groups=con.execute('SELECT COUNT(*) FROM anchor_bounds').fetchone()[0] if bounded else 0
        return {'database':str(p),'indexes':indexes,'bounded_anchor_groups':groups,
                'note':'loop_genes also has PRIMARY KEY(gene_id,loop_id) in its WITHOUT ROWID table.'}
    finally:con.close()

def optimize_copy(database,output):
    source=Path(database).resolve();output=Path(output).resolve()
    if source==output: raise ValueError('Use a different output filename; never update the live database in place.')
    if output.exists(): raise FileExistsError('Output already exists. Choose a new filename.')
    output.parent.mkdir(parents=True,exist_ok=True)
    fd,tmp=tempfile.mkstemp(prefix='.index-',suffix='.sqlite',dir=output.parent);os.close(fd)
    src=dst=None
    try:
        src=sqlite3.connect(source.as_uri()+'?mode=ro',uri=True)
        dst=sqlite3.connect(tmp)
        src.backup(dst)  # Consistent SQLite backup, including any committed WAL contents.
        install_indexes(dst);dst.commit()
        if dst.execute('PRAGMA quick_check').fetchone()[0]!='ok': raise ValueError('Database integrity check failed.')
        dst.close();src.close()
        os.chmod(tmp,0o644)
        # Exclusive creation: do not overwrite a concurrently created output.
        os.link(tmp,output);os.unlink(tmp)
    except BaseException:
        if dst is not None: dst.close()
        if src is not None: src.close()
        Path(tmp).unlink(missing_ok=True)
        raise
    return inspect_database(output)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--database',required=True)
    p.add_argument('--output',help='Write optimized copy here; omit for read-only inspection.')
    a=p.parse_args()
    print(json.dumps(optimize_copy(a.database,a.output) if a.output else inspect_database(a.database),indent=2))
if __name__=='__main__':main()
