"""Unit tests for read orientation inference module."""

import os
from pathlib import Path

import pytest

from htsinfer import infer_read_orientation

dir_path = os.path.abspath(os.path.join(__file__, "../.."))
transcript_fasta_path = os.path.join(
    dir_path, "htsinfer", "transcripts.fasta.zip"
    )
test_files_dir = Path(__file__).parent.absolute() / "test_files"
file_1 = str(test_files_dir / "first_mate.fastq")
file_2 = str(test_files_dir / "second_mate.fastq")
fasta_human = str(test_files_dir / "transcripts_human.fa")
fasta_mixed = str(test_files_dir / "transcripts_mixed.fa")
fasta_no_orgs = str(test_files_dir / "transcripts_no_orgs.fa")


def _raise(exception) -> None:
    """General purpose exception raiser."""
    raise exception


def file_len(fname):
    """Count lines in file."""
    return sum(1 for line in open(fname))


class TestInfer:
    """Tests for the main function `infer()`."""

    def test_single_file(self):
        """Function returns without errors."""
        assert infer_read_orientation.infer(
            fasta=transcript_fasta_path,
            file_1=file_1,
        ) == "U"

    def test_cannot_create_tmp_dir(self, monkeypatch):
        """Fails to create temporary directory."""
        monkeypatch.setattr(
            'tempfile.mkdtemp',
            lambda *args, **kwargs: _raise(OSError),
        )
        with pytest.raises(OSError):
            infer_read_orientation.infer(
                fasta=fasta_human,
                file_1=file_1,
            )

    def test_cannot_delete_tmp_dir(self, monkeypatch, tmp_path):
        """Fails to deleted temporary directory."""
        monkeypatch.setattr(
            'tempfile.mkdtemp',
            lambda *args, **kwargs: str(tmp_path),
        )
        monkeypatch.setattr(
            'shutil.rmtree',
            lambda *args, **kwargs: _raise(OSError),
        )
        with pytest.raises(OSError):
            infer_read_orientation.infer(
                fasta=transcript_fasta_path,
                file_1=file_1,
            )


class TestSubsetFastaByOrgansim:
    """Test for function `subset_fasta_by_organism()`."""

    def test_fasta_subsetting(self, tmp_path):
        """Writes FASTA records of specified organism only."""
        fasta_out = str(tmp_path / "out.fa")
        infer_read_orientation.subset_fasta_by_organism(
            fasta_in=fasta_mixed,
            fasta_out=fasta_out,
            organism="hsapiens",
        )
        assert file_len(fasta_out) == file_len(fasta_human)

    def test_no_orgs(self, tmp_path):
        """All FASTA records are skipped because organism information is
        absent.
        """
        fasta_out = str(tmp_path / "out.fa")
        infer_read_orientation.subset_fasta_by_organism(
            fasta_in=fasta_no_orgs,
            fasta_out=str(tmp_path / "out.fa"),
            organism="hsapiens",
        )
        assert file_len(fasta_out) == 0

    def test_invalid_fasta(self, tmp_path):
        """No FASTA records are written because input file is not of FASTA
        format.
        """
        fasta_out = str(tmp_path / "out.fa")
        infer_read_orientation.subset_fasta_by_organism(
            fasta_in=file_1,
            fasta_out=str(tmp_path / "out.fa"),
            organism="hsapiens",
        )
        assert file_len(fasta_out) == 0

    def test_fasta_na(self, monkeypatch, tmp_path):
        """Input FASTA file cannot be opened."""
        fasta_out = str(tmp_path / "out.fa")
        monkeypatch.setattr(
            'Bio.SeqIO.parse',
            lambda *args, **kwargs: _raise(OSError),
        )
        with pytest.raises(OSError):
            infer_read_orientation.subset_fasta_by_organism(
                fasta_in=fasta_human,
                fasta_out=fasta_out,
                organism="hsapiens",
            )

    def test_out_path_not_writable(self, monkeypatch, tmp_path):
        """Output FASTA file cannot be written."""
        fasta_out = str(tmp_path / "out.fa")
        monkeypatch.setattr(
            'Bio.SeqIO.write',
            lambda *args, **kwargs: _raise(OSError),
        )
        with pytest.raises(OSError):
            infer_read_orientation.subset_fasta_by_organism(
                fasta_in=fasta_human,
                fasta_out=fasta_out,
                organism="hsapiens",
            )
