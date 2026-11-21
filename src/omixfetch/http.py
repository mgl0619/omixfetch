"""A shared HTTP session with retries and polite pacing.

NCBI throttles aggressively: a burst of requests returns an HTML error page
with a 200 status rather than a 429, so retry logic has to inspect the body,
not just the status code.
"""

from __future__ import annotations

import time
from typing import Any

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from .errors import ResolveError

USER_AGENT = "omixfetch/0.2.0 (+https://github.com/mgl0619/omixfetch)"

DEFAULT_TIMEOUT = 120


def make_session(total_retries: int = 5, backoff: float = 1.0) -> requests.Session:
    """Build a session that retries idempotent requests with exponential backoff."""
    session = requests.Session()
    session.headers["User-Agent"] = USER_AGENT
    retry = Retry(
        total=total_retries,
        backoff_factor=backoff,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=frozenset(["GET", "HEAD"]),
        raise_on_status=False,
    )
    adapter = HTTPAdapter(max_retries=retry, pool_maxsize=32)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    return session


def get_text(
    session: requests.Session,
    url: str,
    *,
    expect_text: bool = True,
    tries: int = 4,
    pause: float = 2.0,
    **kwargs: Any,
) -> str:
    """GET *url* and return the body, retrying when an HTML error page comes back.

    Args:
        expect_text: when True, a body that begins with ``<`` is treated as an
            error page (NCBI returns these with a 200 status under load) and
            the request is retried.

    Raises:
        ResolveError: if every attempt fails.
    """
    kwargs.setdefault("timeout", DEFAULT_TIMEOUT)
    last = ""
    for attempt in range(tries):
        try:
            response = session.get(url, **kwargs)
        except requests.RequestException as exc:
            last = str(exc)
        else:
            # NCBI serves UTF-8 text/plain without a charset parameter, and
            # requests then falls back to ISO-8859-1, which mangles quotation
            # marks and accented author names in SOFT records.
            if "charset" not in response.headers.get("Content-Type", "").lower():
                response.encoding = "utf-8"
            body = response.text
            if response.ok and not (expect_text and body.lstrip().startswith("<")):
                return body
            last = f"HTTP {response.status_code}"
        if attempt < tries - 1:
            time.sleep(pause * (attempt + 1))
    raise ResolveError(f"Could not retrieve {url} after {tries} attempts ({last}).")


def get_json(session: requests.Session, url: str, **kwargs: Any) -> Any:
    """GET *url* and parse the body as JSON."""
    kwargs.setdefault("timeout", DEFAULT_TIMEOUT)
    try:
        response = session.get(url, **kwargs)
        response.raise_for_status()
        return response.json()
    except (requests.RequestException, ValueError) as exc:
        raise ResolveError(f"Could not retrieve JSON from {url}: {exc}") from exc


def head_size(session: requests.Session, url: str) -> int | None:
    """Return the ``Content-Length`` of *url*, or ``None`` if unavailable."""
    try:
        response = session.head(url, timeout=DEFAULT_TIMEOUT, allow_redirects=True)
        if response.ok:
            length = response.headers.get("Content-Length")
            return int(length) if length is not None else None
    except (requests.RequestException, ValueError):
        pass
    return None
