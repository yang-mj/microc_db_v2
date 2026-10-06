"""Application limits and process-local concurrency admission; no shared file handles."""
from contextlib import contextmanager
import logging
import os
import threading
import time
import uuid
from collections import deque

class RequestError(ValueError):
    """A controlled message safe to present to public users."""

class BusyError(RequestError): pass

MAX_ROWS = 50_000
MAX_PLOT_LOOPS = 3_000
MAX_SIGNAL_BINS = 3_600
MAX_TRACKS = 16
QUERY_SECONDS = 10

class AdmissionGate:
    def __init__(self, capacity):
        if not 1 <= capacity <= 8: raise ValueError('Concurrency capacity must be 1–8.')
        self._slots = threading.BoundedSemaphore(capacity)
    @contextmanager
    def acquire(self):
        if not self._slots.acquire(blocking=False):
            raise BusyError('The server is processing other plots. Please try again shortly.')
        try: yield
        finally: self._slots.release()

class RateLimiter:
    """Process-wide sliding-window limit on expensive requests (0 disables it).

    Protects CPU and memory from floods that open many sessions, which a per-session cooldown cannot
    do. Per-visitor limits need the visitor's real address, which only the hosting proxy knows
    reliably; configure those at the ingress.
    """
    def __init__(self, max_events, window=60.0, clock=time.monotonic):
        self.max_events=max_events; self.window=window; self.clock=clock
        self._events=deque(); self._lock=threading.Lock()
    def allow(self):
        if self.max_events<=0: return True
        now=self.clock()
        with self._lock:
            while self._events and now-self._events[0]>=self.window: self._events.popleft()
            if len(self._events)>=self.max_events: return False
            self._events.append(now); return True
    def check(self):
        if not self.allow(): raise BusyError('The server is receiving many requests right now. Please try again in a minute.')

def _env_int(name, default, low, high):
    value=int(os.environ.get(name, default))
    if not low<=value<=high: raise ValueError(f'{name} must be {low}–{high}.')
    return value

# Initialized once when this module is imported, not on each Streamlit rerun.
REQUEST_GATE = AdmissionGate(int(os.environ.get('MICROC_MAX_CONCURRENT', '2')))
# Uncached plots and input checks per minute, all visitors combined (0 = unlimited).
GLOBAL_LIMITER = RateLimiter(_env_int('MICROC_RATE_PER_MINUTE', 240, 0, 100_000))

def report_failure(context, exc):
    """Opaque reference for users; full diagnostics stay in restricted server logs."""
    reference = uuid.uuid4().hex[:12]
    logging.getLogger('microc_explorer').error('%s [%s]', context, reference,
        exc_info=(type(exc), exc, exc.__traceback__))
    return f'{context}. Please contact the app administrator with reference {reference}.'

def validate_request(source, value, mode, bins, statistic, resolution, cfg, tracks, contact_id, formats=('svg','png','pdf')):
    if not formats or not set(formats)<={'svg','png','pdf'}: raise RequestError('Choose SVG, PNG or PDF.')
    if not isinstance(value,str) or not 1 <= len(value.strip()) <= 128:
        raise RequestError('Enter a gene symbol or genomic interval of at most 128 characters.')
    if mode not in ('gene','region'): raise RequestError('Choose gene or coordinate search.')
    if isinstance(bins,bool) or not isinstance(bins,int) or not 1<=bins<=MAX_SIGNAL_BINS:
        raise RequestError(f'Display bins must be 1–{MAX_SIGNAL_BINS}.')
    if statistic not in ('max','mean'): raise RequestError('Select maximum or mean signal.')
    if resolution is not None and (isinstance(resolution,bool) or not isinstance(resolution,int) or not 0<resolution<=10_000_000):
        raise RequestError('Choose a positive contact-map resolution, or automatic.')
    configured = {t['id'] for t in cfg.get('bigwigs',[])}
    chosen = list(configured) if tracks is None else list(tracks)
    if len(chosen)>MAX_TRACKS or len(chosen)!=len(set(chosen)) or not set(chosen)<=configured:
        raise RequestError('Select only the configured signal tracks, without duplicates.')
    if contact_id is not None and contact_id not in {t['id'] for t in cfg.get('contacts',[])}:
        raise RequestError('Select a configured contact map.')
