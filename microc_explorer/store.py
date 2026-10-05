"""Bounded, thread-safe in-memory web exports; no plot files are written to disk."""
from collections import OrderedDict
import threading
import time
import os

class ExportStore:
    def __init__(self,max_bytes=100_000_000,ttl=120,clock=time.monotonic):
        self.max_bytes=max_bytes
        self.ttl=ttl
        self.clock=clock
        self._entries=OrderedDict()
        self._bytes=0
        self._lock=threading.RLock()

    def _drop(self,owner):
        entry=self._entries.pop(owner,None)
        if entry is not None:self._bytes-=entry[2]

    def sweep(self):
        with self._lock:
            now=self.clock()
            for owner,(_,seen,_) in list(self._entries.items()):
                if now-seen>=self.ttl:self._drop(owner)

    def put(self,owner,exports):
        """Replace this session's exports. Evict least recently accessed if needed."""
        size=sum(len(value) for value in exports.values())
        with self._lock:
            self.sweep()
            self._drop(owner)
            if size>self.max_bytes:return False
            if not exports:return True
            while self._bytes+size>self.max_bytes:
                self._drop(next(iter(self._entries)))
            self._entries[owner]=(dict(exports),self.clock(),size)
            self._bytes+=size
            return True

    def get(self,owner):
        with self._lock:
            self.sweep()
            entry=self._entries.get(owner)
            if entry is None:return None
            data,_,size=entry
            self._entries[owner]=(data,self.clock(),size)
            self._entries.move_to_end(owner)
            return dict(data)

    def drop(self,owner):
        with self._lock:self._drop(owner)

    @property
    def size_bytes(self):
        with self._lock:
            self.sweep()
            return self._bytes

def _budget():
    value=int(os.environ.get('MICROC_EXPORT_CACHE_MB','100'))
    if not 1<=value<=100: raise ValueError('MICROC_EXPORT_CACHE_MB must be 1–100 (decimal MB).')
    return value*1_000_000

STORE=ExportStore(max_bytes=_budget())
_started=False
_start_lock=threading.Lock()

def start_janitor():
    """One daemon per process; bounded background expiry even without new requests."""
    global _started
    with _start_lock:
        if _started:return
        def worker():
            while True:
                time.sleep(15)
                STORE.sweep()
        threading.Thread(target=worker,name='microc-export-cleanup',daemon=True).start()
        _started=True
