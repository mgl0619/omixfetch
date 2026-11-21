"""Backend command construction and the aria2 input-file format."""

from pathlib import Path

import pytest

from omixfetch.backends.aria2 import (
    Aria2Backend,
    aria2_available,
    ca_bundle,
    render_input_file,
    subprocess_env,
)
from omixfetch.backends.fallback import RequestsBackend, file_md5
from omixfetch.errors import ChecksumError
from omixfetch.models import RemoteFile


def test_render_input_file_sets_dir_and_out(tmp_path):
    files = [
        RemoteFile(url="https://example.org/a.tar", filename="a.tar"),
        RemoteFile(url="https://example.org/b.gz", filename="b.gz", md5="d41d8cd9"),
    ]
    text = render_input_file(files, tmp_path)
    lines = text.splitlines()

    assert lines[0] == "https://example.org/a.tar"
    assert lines[1] == f"  dir={tmp_path}"
    assert lines[2] == "  out=a.tar"
    # Only the file with a published checksum gets a checksum line.
    assert "  checksum=md5=d41d8cd9" in lines
    assert sum(1 for line in lines if line.startswith("  checksum=")) == 1


def test_render_input_file_empty():
    assert render_input_file([], Path("/tmp")) == ""


@pytest.mark.skipif(not aria2_available(), reason="aria2c not installed")
def test_build_command_carries_parallelism(tmp_path):
    backend = Aria2Backend(connections=16, concurrent=8)
    command = backend.build_command(tmp_path / "in.txt")
    assert "--max-connection-per-server=16" in command
    assert "--split=16" in command
    assert "--max-concurrent-downloads=8" in command
    assert "--continue=true" in command


def test_subprocess_env_drops_socks_all_proxy(monkeypatch):
    """aria2c aborts on a SOCKS all_proxy, so it must not be inherited."""
    monkeypatch.setenv("all_proxy", "socks5h://localhost:49954")
    monkeypatch.setenv("http_proxy", "http://localhost:49953")
    env = subprocess_env()
    assert "all_proxy" not in env
    # HTTP proxies are still honoured.
    assert env["http_proxy"] == "http://localhost:49953"


def test_subprocess_env_keeps_http_all_proxy(monkeypatch):
    monkeypatch.setenv("ALL_PROXY", "http://proxy.example.org:3128")
    assert subprocess_env()["ALL_PROXY"] == "http://proxy.example.org:3128"


def test_ca_bundle_prefers_environment(monkeypatch, tmp_path):
    bundle = tmp_path / "ca.pem"
    bundle.write_text("")
    monkeypatch.setenv("SSL_CERT_FILE", str(bundle))
    assert ca_bundle() == str(bundle)


def test_ca_bundle_ignores_missing_path(monkeypatch, tmp_path):
    monkeypatch.setenv("SSL_CERT_FILE", str(tmp_path / "absent.pem"))
    monkeypatch.delenv("REQUESTS_CA_BUNDLE", raising=False)
    monkeypatch.delenv("CURL_CA_BUNDLE", raising=False)
    # Falls through to certifi, which is installed alongside requests.
    result = ca_bundle()
    assert result is None or result.endswith(".pem")


def test_file_md5(tmp_path):
    p = tmp_path / "x.bin"
    p.write_bytes(b"hello")
    assert file_md5(p) == "5d41402abc4b2a76b9719d911017c592"


class _FakeResponse:
    status_code = 200

    def __init__(self, payload: bytes):
        self._payload = payload

    def raise_for_status(self):
        return None

    def iter_content(self, _size):
        yield self._payload


def _backend_with(payload: bytes):
    backend = RequestsBackend(concurrent=1, max_tries=2, quiet=True)
    backend.session = type(
        "S", (), {"get": staticmethod(lambda *a, **k: _FakeResponse(payload))}
    )()
    return backend


def test_fallback_accepts_matching_md5(tmp_path):
    backend = _backend_with(b"hello")
    f = RemoteFile(
        url="https://example.org/x.bin",
        filename="x.bin",
        size=5,
        md5="5d41402abc4b2a76b9719d911017c592",
    )
    backend.download([f], tmp_path)
    assert (tmp_path / "x.bin").read_bytes() == b"hello"


def test_fallback_rejects_bad_md5_and_removes_partial(tmp_path):
    backend = _backend_with(b"hello")
    f = RemoteFile(
        url="https://example.org/x.bin", filename="x.bin", size=5, md5="0" * 32
    )
    with pytest.raises(Exception) as excinfo:
        backend.download([f], tmp_path)
    assert "md5 mismatch" in str(excinfo.value)
    # Neither the corrupt partial nor a misleading final file is left behind.
    assert not (tmp_path / "x.bin").exists()
    assert not (tmp_path / "x.bin.part").exists()
