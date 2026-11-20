"""Dispatch accessions to the right archive resolver."""

from __future__ import annotations

import requests

from .accession import Archive, parse_accession
from .http import make_session
from .models import Dataset
from .sources.geo import resolve_geo

_RESOLVERS = {
    Archive.GEO: resolve_geo,
}


def resolve(
    accession: str, session: requests.Session | None = None, **kwargs
) -> Dataset:
    """Resolve a supported accession into a :class:`Dataset`.

    Args:
        accession: a GEO accession (GSE, GSM or GPL).
        session: optional shared HTTP session; one is created if omitted.
        **kwargs: forwarded to the archive-specific resolver.

    Raises:
        AccessionError: if the accession belongs to an unsupported archive.
        ResolveError: if the archive cannot be queried.
    """
    archive, normalised = parse_accession(accession)
    return _RESOLVERS[archive](normalised, session or make_session(), **kwargs)
