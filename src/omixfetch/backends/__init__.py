"""Transfer backends."""

from .aria2 import Aria2Backend, aria2_available, aria2_version
from .aspera import AsperaBackend, ascp_available
from .fallback import RequestsBackend

__all__ = [
    "Aria2Backend",
    "AsperaBackend",
    "RequestsBackend",
    "aria2_available",
    "aria2_version",
    "ascp_available",
]
