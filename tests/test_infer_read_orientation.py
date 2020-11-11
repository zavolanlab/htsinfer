"""Unit tests for read orientation inference module."""

import os

import pytest

from htsinfer import infer_read_orientation

test_files_dir = os.path.dirname(__file__)
test_files_dir = os.path.join(test_files_dir, "test_files")
file_1 = os.path.join(test_files_dir, "first_mate.fastq")
file_2 = os.path.join(test_files_dir, "second_mate.fastq")
file_na = os.path.join(test_files_dir, "not_available")


def _raise(exception) -> None:
    """General purpose exception raiser."""
    raise exception


class TestInfer:
    """Tests for the main function `infer()`."""

    def test_single_file(self):
        assert infer_read_orientation.infer(file_1=file_1) == "U"

    def test_cannot_create_tmp_dir(self, monkeypatch):
        monkeypatch.setattr(
            'tempfile.mkdtemp',
            lambda *args, **kwargs: _raise(OSError),
        )
        with pytest.raises(OSError):
            infer_read_orientation.infer(file_1=file_1)

    def test_cannot_delete_tmp_dir(self, monkeypatch):
        monkeypatch.setattr(
            'shutil.rmtree',
            lambda *args, **kwargs: _raise(OSError),
        )
        with pytest.raises(OSError):
            infer_read_orientation.infer(file_1=file_1)
