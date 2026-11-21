"""Dispatch accessions to the right archive resolver."""

from __future__ import annotations

import requests

from .accession import Archive, parse_accession
from .errors import AccessionError
from .http import make_session
from .models import Dataset
from .sources.geo import resolve_geo

_RESOLVERS = {
    Archive.GEO: resolve_geo,
}


def resolve(
    accession: str,
    session: requests.Session | None = None,
    resolver: str = "native",
    **kwargs,
) -> Dataset:
    """Resolve a supported accession into a :class:`Dataset`.

    Args:
        accession: a GEO accession (GSE, GSM or GPL).
        session: optional shared HTTP session; one is created if omitted.
        resolver: ``"native"`` for the built-in SOFT parser, or ``"geofetch"``
            to delegate metadata resolution to geofetch (optional extra).
        **kwargs: forwarded to the archive-specific resolver.

    Raises:
        AccessionError: if the accession belongs to an unsupported archive.
        ResolveError: if the archive cannot be queried.
    """
    archive, normalised = parse_accession(accession)
    http = session or make_session()

    if resolver == "geofetch":
        from .sources.geofetch_resolver import resolve_geo_via_geofetch

        kwargs.pop("fetch_sizes", None)
        kwargs.pop("include_archives", None)
        kwargs.pop("include_matrix", None)
        kwargs.pop("size_workers", None)
        return resolve_geo_via_geofetch(normalised, http, **kwargs)
    if resolver != "native":
        raise AccessionError(
            f"Unknown resolver {resolver!r}. Use 'native' or 'geofetch'."
        )

    return _RESOLVERS[archive](normalised, http, **kwargs)
