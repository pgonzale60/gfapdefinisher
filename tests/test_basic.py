"""Basic smoke tests for gfapdefinisher."""

import pytest

from gfapdefinisher import hello
from gfapdefinisher.cli import main


def test_hello() -> None:
    assert hello() == "Hello from gfapdefinisher!"


def test_cli(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as exc_info:
        main([])
    assert exc_info.value.code == 0
    captured = capsys.readouterr()
    assert "gfapdefinisher" in captured.out
