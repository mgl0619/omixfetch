"""Dispatch accessions to the right archive resolver."""

from __future__ import annotations

import requests

from .accession import GEO_RE, Archive, parse_accession
from .errors import AccessionError
from .http import make_session
from .models import Dataset
from .sources.geo import resolve_geo
from .sources.project import resolve_project

_RESOLVERS = {
    Archive.GEO: resolve_geo,
    Archive.SRA: resolve_project,
}

#: Arguments the geofetch backend does not accept.
_NATIVE_ONLY_KWARGS = ("fetch_sizes", "include_archives", "include_matrix", "size_workers")

#: Arguments only the geofetch backend accepts.
_GEOFETCH_ONLY_KWARGS = ("verbose",)

VALID_RESOLVERS = ("auto", "native", "geofetch")


def choose_resolver(
    accession: str, requested: str = "auto", include_reads: bool = False
) -> str:
    """Decide which resolver to use for *accession*.

    ``auto`` prefers geofetch, the more mature and widely exercised GEO
    parser, but only where it is actually the better source:

    * GSM and GPL accessions always use the built-in resolver. geofetch's
      project API returns an empty result for the former and raises on the
      latter.
    * Raw-read workflows always use the built-in resolver. geofetch's sample
      table carries no ``relation`` columns, so the BioSample identifiers that
      link an ENA run back to its GEO sample are unavailable -- runs would be
      left ungrouped. It also returns no sample table at all for series that
      have no processed files, which is common for reads-only submissions.
    """
    if not GEO_RE.match(accession):
        # Only GEO has a choice of parser; read-archive projects are always
        # resolved from ENA.
        return "native"

    if requested != "auto":
        return requested

    if include_reads:
        return "native"

    match = GEO_RE.match(accession)
    is_series = bool(match) and match.group(1).upper() == "GSE"
    if not is_series:
        return "native"

    from .sources.geofetch_resolver import geofetch_available

    return "geofetch" if geofetch_available() else "native"


def resolve(
    accession: str,
    session: requests.Session | None = None,
    resolver: str = "auto",
    **kwargs,
) -> Dataset:
    """Resolve a supported accession into a :class:`Dataset`.

    Args:
        accession: a GEO accession (GSE, GSM or GPL), or a read-archive
            accession (PRJNA..., SRP..., SRR..., SAMN... and their ENA/DDBJ
            equivalents).
        session: optional shared HTTP session; one is created if omitted.
        resolver: ``"auto"`` (geofetch for series when installed, else
            built-in), ``"native"`` or ``"geofetch"``.
        **kwargs: forwarded to the archive-specific resolver.

    Raises:
        AccessionError: for an unsupported archive or unknown resolver.
        ResolveError: if the archive cannot be queried.
    """
    if resolver not in VALID_RESOLVERS:
        raise AccessionError(
            f"Unknown resolver {resolver!r}. Use one of: {', '.join(VALID_RESOLVERS)}."
        )

    archive, normalised = parse_accession(accession)
    http = session or make_session()
    chosen = choose_resolver(
        normalised, resolver, include_reads=bool(kwargs.get("include_reads"))
    )

    if chosen == "geofetch":
        from .sources.geofetch_resolver import resolve_geo_via_geofetch

        geofetch_kwargs = {
            k: v for k, v in kwargs.items() if k not in _NATIVE_ONLY_KWARGS
        }
        dataset = resolve_geo_via_geofetch(normalised, http, **geofetch_kwargs)

        # geofetch returns only a series table for submissions with no
        # processed files; fall back rather than report a dataset with no
        # samples.
        if kwargs.get("include_samples", True) and not dataset.samples:
            native_kwargs = {
                k: v for k, v in kwargs.items() if k not in _GEOFETCH_ONLY_KWARGS
            }
            dataset = _RESOLVERS[archive](normalised, http, **native_kwargs)
            chosen = "native"
    else:
        native_kwargs = {
            k: v for k, v in kwargs.items() if k not in _GEOFETCH_ONLY_KWARGS
        }
        dataset = _RESOLVERS[archive](normalised, http, **native_kwargs)

    dataset.metadata.resolved_by = chosen
    return dataset
