"""Finding and building the conversion toolkit without user configuration."""

import pytest

from omixfetch import convert
from omixfetch.convert import (
    DEFAULT_CONVERT_ENV,
    ENV_NAME_VAR,
    PATH_ENV_VAR,
    build_command,
    create_convert_env,
    discover,
    env_prefix,
)
from omixfetch.errors import BackendError


@pytest.fixture
def clean(monkeypatch, tmp_path):
    """No toolkit anywhere, no configuration -- including no bundled copy."""
    monkeypatch.setenv("OMIXFETCH_HOME", str(tmp_path / "empty-home"))
    monkeypatch.delenv(PATH_ENV_VAR, raising=False)
    monkeypatch.delenv(ENV_NAME_VAR, raising=False)
    monkeypatch.setattr(convert.shutil, "which", lambda name: None)
    monkeypatch.setattr(convert, "env_exists", lambda name: False)
    return monkeypatch


def _binary(tmp_path):
    path = tmp_path / "fasterq-dump"
    path.write_text("#!/bin/sh\n")
    path.chmod(0o755)
    return path


def test_discover_finds_nothing_when_nothing_exists(clean):
    assert discover() == (None, None)


def test_discover_uses_the_setup_environment_automatically(clean):
    """After `omixfetch setup`, no flags or exports are needed."""
    clean.setattr(convert, "env_exists", lambda name: name == DEFAULT_CONVERT_ENV)
    assert discover() == (None, DEFAULT_CONVERT_ENV)


def test_discover_prefers_path_over_the_setup_environment(clean, tmp_path):
    """On Linux sra-tools sits beside aria2; no second environment needed."""
    binary = _binary(tmp_path)
    clean.setattr(convert.shutil, "which",
                  lambda name: str(binary) if name == "fasterq-dump" else None)
    clean.setattr(convert, "env_exists", lambda name: True)
    assert discover() == (str(binary), None)


def test_discover_honours_an_explicit_binary_first(clean, tmp_path):
    binary = _binary(tmp_path)
    clean.setattr(convert, "env_exists", lambda name: True)
    assert discover(str(binary), "other-env") == (str(binary), None)


def test_discover_honours_a_named_environment(clean):
    assert discover(None, "my-sra") == (None, "my-sra")
    clean.setenv(ENV_NAME_VAR, "from-env-var")
    assert discover() == (None, "from-env-var")


def test_create_convert_env_command(clean):
    clean.setattr(convert.shutil, "which",
                  lambda name: "/usr/bin/mamba" if name == "mamba" else None)
    cmd = create_convert_env()
    assert cmd[:4] == ["/usr/bin/mamba", "create", "-y", "-n"]
    assert DEFAULT_CONVERT_ENV in cmd
    # conda-forge must precede bioconda or the solver may pick sra-tools 3.2.
    assert cmd.index("conda-forge") < cmd.index("bioconda")
    assert "sra-tools>=3.4" in cmd


def test_create_convert_env_needs_a_launcher(clean):
    with pytest.raises(BackendError, match="micromamba, mamba or conda"):
        create_convert_env()


def test_env_prefix_runs_in_another_environment(clean):
    clean.setattr(convert.shutil, "which",
                  lambda name: "/opt/micromamba" if name == "micromamba" else None)
    assert env_prefix("omixfetch-convert") == [
        "/opt/micromamba", "run", "-n", "omixfetch-convert"
    ]
    assert env_prefix(None) == []


def test_build_command_wraps_fasterq_dump_in_the_environment(clean, tmp_path):
    clean.setattr(convert.shutil, "which",
                  lambda name: "/opt/mamba" if name == "mamba" else None)
    cmd = build_command(tmp_path / "a.sra", tmp_path, env_name="omixfetch-convert")
    assert cmd[:4] == ["/opt/mamba", "run", "-n", "omixfetch-convert"]
    assert cmd[4] == "fasterq-dump"
