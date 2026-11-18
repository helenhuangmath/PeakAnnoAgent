"""Small stdlib HTTP helpers with retries, so the core package has no dependencies."""

from __future__ import annotations

import json
import os
import time
import urllib.error
import urllib.parse
import urllib.request

from . import __version__
from .engine.common import _ssl_contexts

USER_AGENT = f"PeakAnnoAgent/{__version__} (+https://github.com/helenhuangmath/PeakAnnoAgent)"


class HttpError(RuntimeError):
    pass


_context = [None]      # the certificate store that worked, reused for later requests


def _open(req, timeout: int) -> bytes:
    """Open a request, trying the next certificate store only when verification fails because this
    Python has no root certificates (python.org builds on macOS). Verification is never turned off."""
    import ssl
    contexts = [_context[0]] if _context[0] is not None else _ssl_contexts()
    error = None
    for context in contexts:
        try:
            with urllib.request.urlopen(req, timeout=timeout, context=context) as resp:
                data = resp.read()
            _context[0] = context
            return data
        except urllib.error.HTTPError:
            raise
        except urllib.error.URLError as err:
            error = err
            if not isinstance(err.reason, ssl.SSLCertVerificationError):
                raise
    raise error


def _request(url: str, data: bytes = None, headers: dict = None, timeout: int = 60, retries: int = 2) -> bytes:
    hdrs = {"User-Agent": USER_AGENT}
    hdrs.update(headers or {})
    last = None
    for attempt in range(retries + 1):
        try:
            req = urllib.request.Request(url, data=data, headers=hdrs)
            return _open(req, timeout)
        except urllib.error.HTTPError as err:
            last = err
            if err.code not in (429, 500, 502, 503, 504):
                break
        except (urllib.error.URLError, TimeoutError, ConnectionError) as err:
            last = err
        time.sleep(2 * (attempt + 1))
    raise HttpError(f"Request failed for {url.split('?')[0]}: {last}")


def get_text(url: str, params: dict = None, timeout: int = 60) -> str:
    if params:
        url = f"{url}?{urllib.parse.urlencode(params)}"
    return _request(url, timeout=timeout).decode("utf-8", errors="replace")


def get_json(url: str, params: dict = None, timeout: int = 60):
    return json.loads(get_text(url, params, timeout))


def post_form(url: str, fields: dict, timeout: int = 120) -> str:
    body = urllib.parse.urlencode(fields).encode()
    return _request(url, data=body, timeout=timeout, retries=1,
                    headers={"Content-Type": "application/x-www-form-urlencoded"}).decode("utf-8", "replace")


_last_ncbi = [0.0]


def ncbi_get_json(endpoint: str, params: dict):
    """NCBI E-utilities with the documented rate limit (3/s, or 10/s with NCBI_API_KEY)."""
    params = dict(params, retmode="json", tool="PeakAnnoAgent")
    key = os.environ.get("NCBI_API_KEY")
    if key:
        params["api_key"] = key
    gap = 0.11 if key else 0.35
    wait = _last_ncbi[0] + gap - time.time()
    if wait > 0:
        time.sleep(wait)
    _last_ncbi[0] = time.time()
    return get_json(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/{endpoint}.fcgi", params)
