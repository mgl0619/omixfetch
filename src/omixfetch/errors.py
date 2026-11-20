"""Exception hierarchy for omixfetch."""


class OmixfetchError(Exception):
    """Base class for every error raised by omixfetch."""


class AccessionError(OmixfetchError):
    """The accession is malformed or belongs to an unsupported archive."""


class ResolveError(OmixfetchError):
    """The archive could not be queried, or returned an unusable record."""


class BackendError(OmixfetchError):
    """The download backend failed or is unavailable."""


class ChecksumError(OmixfetchError):
    """A downloaded file did not match its expected size or checksum."""
