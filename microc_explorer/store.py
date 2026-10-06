"""Bounded, thread-safe, in-memory cache of finished plots shared by all sessions.

The data are curated and public, so one visitor's plot can safely serve the next visitor who asks
the identical question (same dataset files, query, tracks and settings). The cache key is a hash of
validated request parameters only; nothing session-specific is ever stored. No plot files touch disk.

Memory is bounded by a hard byte budget (least recently used entries are evicted first) and by an
idle time-to-live; a background thread sweeps expired entries even when the app is idle.
"""
from collections import OrderedDict
import os
import threading
import time

def _nbytes(render):
    """Approximate memory held by the inputs kept for on-demand re-rendering."""
    if not render: return 0
    total=0
    for s in render.get('signals',()): total+=s.edges.nbytes+s.values.nbytes
    c=render.get('contact')
    if c is not None: total+=c.edges.nbytes+c.matrix.nbytes
    return total

class PlotCache:
    def __init__(self,max_bytes=100_000_000,ttl=900,clock=time.monotonic):
        self.max_bytes=max_bytes; self.ttl=ttl; self.clock=clock
        self._entries=OrderedDict()   # key -> dict(payload, render, exports, size, seen)
        self._bytes=0
        self._lock=threading.RLock()

    def _drop(self,key):
        entry=self._entries.pop(key,None)
        if entry is not None: self._bytes-=entry['size']

    def sweep(self):
        with self._lock:
            now=self.clock()
            for key,entry in list(self._entries.items()):
                if now-entry['seen']>=self.ttl: self._drop(key)

    def _make_room(self,size):
        while self._entries and self._bytes+size>self.max_bytes:
            self._drop(next(iter(self._entries)))

    def put(self,key,payload,render,exports):
        """Store a finished plot. Returns False when it alone exceeds the budget."""
        size=sum(len(v) for v in exports.values())+_nbytes(render)+int(payload['result'].loops.memory_usage(deep=True).sum())
        with self._lock:
            self.sweep(); self._drop(key)
            if size>self.max_bytes: return False
            self._make_room(size)
            self._entries[key]=dict(payload=payload,render=render,exports=dict(exports),size=size,seen=self.clock())
            self._bytes+=size
            return True

    def get(self,key):
        """Snapshot of an entry (refreshes its recency) or None."""
        with self._lock:
            self.sweep()
            entry=self._entries.get(key)
            if entry is None: return None
            entry['seen']=self.clock(); self._entries.move_to_end(key)
            return dict(payload=entry['payload'],render=entry['render'],exports=dict(entry['exports']))

    def add_exports(self,key,exports):
        """Attach extra formats to an existing entry. Returns False if it is gone or would not fit."""
        extra=sum(len(v) for v in exports.values())
        with self._lock:
            entry=self._entries.get(key)
            if entry is None: return False
            if entry['size']+extra>self.max_bytes: return False
            self._bytes-=entry['size']; entry['size']+=extra
            self._bytes+=entry['size']
            entry['exports'].update(exports); entry['seen']=self.clock()
            self._entries.move_to_end(key)
            # Evict others (never this entry) until within budget.
            while self._bytes>self.max_bytes and len(self._entries)>1:
                oldest=next(iter(self._entries))
                if oldest==key: self._entries.move_to_end(key); continue
                self._drop(oldest)
            return True

    def drop(self,key):
        with self._lock: self._drop(key)

    @property
    def size_bytes(self):
        with self._lock:
            self.sweep(); return self._bytes

def _budget():
    value=int(os.environ.get('MICROC_EXPORT_CACHE_MB','100'))
    if not 1<=value<=100: raise ValueError('MICROC_EXPORT_CACHE_MB must be 1–100 (decimal MB).')
    return value*1_000_000

PLOT_CACHE=PlotCache(max_bytes=_budget())
_started=False
_start_lock=threading.Lock()

def start_janitor():
    """One daemon per process; expired plots are released even without new requests."""
    global _started
    with _start_lock:
        if _started: return
        def worker():
            while True:
                time.sleep(30)
                PLOT_CACHE.sweep()
        threading.Thread(target=worker,name='microc-cache-cleanup',daemon=True).start()
        _started=True
