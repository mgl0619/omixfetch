"""Transfer backends."""

from .aria2 import Aria2Backend, aria2_available, aria2_version
from .fallback import RequestsBackend

__all__ = ["Aria2Backend", "RequestsBackend", "aria2_available", "aria2_version"]
