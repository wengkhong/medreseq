"""Tests for ReferenceGenome sequence retrieval (requires pysam)."""

from __future__ import annotations

from pathlib import Path

import pytest

from medreseq.sequences import ReferenceGenome
from medreseq.variants import Variant

pytestmark = pytest.mark.integration


def test_fetch_around_snv_correct_length(synthetic_fasta: tuple[Path, str]) -> None:
    fasta_path, _ = synthetic_fasta
    v = Variant('chr1', 5000, 'A', 'T')
    with ReferenceGenome(fasta_path) as ref:
        seq, offset = ref.fetch_around(v, flank=100)
    # 100 left + 1 bp ref + 100 right = 201
    assert len(seq) == 201
    assert offset == 100


def test_fetch_around_snv_matches_raw_seq(synthetic_fasta: tuple[Path, str]) -> None:
    fasta_path, raw = synthetic_fasta
    v = Variant('chr1', 5000, 'A', 'T')
    with ReferenceGenome(fasta_path) as ref:
        seq, _ = ref.fetch_around(v, flank=100)
    # raw is 0-based; variant pos 5000 → 0-based index 4999
    expected = raw[4900:5101].upper()
    assert seq == expected


def test_fetch_around_deletion_wider_target(synthetic_fasta: tuple[Path, str]) -> None:
    fasta_path, raw = synthetic_fasta
    v = Variant('chr1', 5000, 'AC', 'A')   # 2-bp deletion
    with ReferenceGenome(fasta_path) as ref:
        seq, offset = ref.fetch_around(v, flank=100)
    # 100 + len('AC') + 100 = 202
    assert len(seq) == 202
    assert offset == 100


def test_fetch_clamps_at_chromosome_start(synthetic_fasta: tuple[Path, str]) -> None:
    fasta_path, _ = synthetic_fasta
    v = Variant('chr1', 10, 'A', 'T')   # near start; flank > pos
    with ReferenceGenome(fasta_path) as ref:
        seq, offset = ref.fetch_around(v, flank=200)
    # Cannot go before position 0
    assert offset < 200
    assert len(seq) < 401   # clamped


def test_fetch_is_uppercase(synthetic_fasta: tuple[Path, str]) -> None:
    fasta_path, _ = synthetic_fasta
    v = Variant('chr1', 5000, 'A', 'T')
    with ReferenceGenome(fasta_path) as ref:
        seq, _ = ref.fetch_around(v, flank=50)
    assert seq == seq.upper()


def test_missing_fai_raises(tmp_path: Path) -> None:
    fasta = tmp_path / 'ref.fa'
    fasta.write_text('>chr1\nACGT\n')
    # No .fai created
    with pytest.raises(FileNotFoundError, match="fai"):
        ReferenceGenome(fasta)


def test_missing_fasta_raises(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError, match="FASTA not found"):
        ReferenceGenome(tmp_path / 'nonexistent.fa')


def test_context_manager(synthetic_fasta: tuple[Path, str]) -> None:
    fasta_path, _ = synthetic_fasta
    with ReferenceGenome(fasta_path) as ref:
        assert 'chr1' in ref.references
