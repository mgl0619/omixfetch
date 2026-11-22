"""Transfer backends."""

from .aria2 import Aria2Backend, aria2_available, aria2_version
from .aspera import AsperaBackend, ascp_available
from .fallback import RequestsBackend
from .s5cmd import S5cmdBackend, s3_uri, s5cmd_available, s5cmd_version

__all__ = [
    "Aria2Backend",
    "AsperaBackend",
    "RequestsBackend",
    "S5cmdBackend",
    "aria2_available",
    "aria2_version",
    "ascp_available",
    "s3_uri",
    "s5cmd_available",
    "s5cmd_version",
]
