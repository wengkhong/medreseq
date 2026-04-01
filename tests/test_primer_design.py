"""Tests for primer design via primer3-py (requires primer3-py)."""

from __future__ import annotations

from pathlib import Path

import pytest

from medreseq.primer_design import PrimerPair, design_primers
from medreseq.sequences import ReferenceGenome
from medreseq.variants import Variant

pytestmark = pytest.mark.integration


# A ~600 bp sequence with ~50 % GC and no long repeats, known to yield primers.
# Generated from synthetic_fasta seed-42 bases 4700..5300.
_GOOD_SEQ_OFFSET = 200   # target centred at 200 in a 401 bp window


def test_design_returns_primer_pairs(synthetic_fasta: tuple[Path, str]) -> None:
    _, raw = synthetic_fasta
    window = raw[4800:5201].upper()     # 401 bp; target at offset 200, len 1
    pairs = design_primers(window, target_offset=200, target_len=1, num_return=3)
    assert isinstance(pairs, list)
    assert len(pairs) >= 1


def test_primer_sequences_are_nonempty(synthetic_fasta: tuple[Path, str]) -> None:
    _, raw = synthetic_fasta
    window = raw[4800:5201].upper()
    pairs = design_primers(window, target_offset=200, target_len=1)
    for p in pairs:
        assert len(p.left_seq) >= 18
        assert len(p.right_seq) >= 18


def test_tm_within_requested_range(synthetic_fasta: tuple[Path, str]) -> None:
    _, raw = synthetic_fasta
    window = raw[4800:5201].upper()
    pairs = design_primers(window, target_offset=200, target_len=1,
                           min_tm=57.0, opt_tm=60.0, max_tm=63.0)
    for p in pairs:
        assert 57.0 <= p.left_tm <= 63.0
        assert 57.0 <= p.right_tm <= 63.0


def test_product_size_within_range(synthetic_fasta: tuple[Path, str]) -> None:
    _, raw = synthetic_fasta
    window = raw[4800:5201].upper()
    pairs = design_primers(window, target_offset=200, target_len=1,
                           product_size_range=(200, 500))
    for p in pairs:
        assert 200 <= p.product_size <= 500


def test_ranks_are_sequential(synthetic_fasta: tuple[Path, str]) -> None:
    _, raw = synthetic_fasta
    window = raw[4800:5201].upper()
    pairs = design_primers(window, target_offset=200, target_len=1, num_return=3)
    for i, p in enumerate(pairs):
        assert p.rank == i


def test_homopolymer_sequence_yields_no_primers() -> None:
    """primer3 cannot design primers in a degenerate sequence."""
    bad = 'A' * 500
    pairs = design_primers(bad, target_offset=200, target_len=1)
    assert pairs == []


def test_design_primers_for_deletion(synthetic_fasta: tuple[Path, str]) -> None:
    """A 2-bp deletion target (len=2) should still yield primers."""
    _, raw = synthetic_fasta
    window = raw[4800:5202].upper()   # 402 bp; ref=2bp at offset 200
    pairs = design_primers(window, target_offset=200, target_len=2)
    assert len(pairs) >= 1


def test_integration_with_reference_genome(synthetic_fasta: tuple[Path, str]) -> None:
    """End-to-end: fetch sequence from FASTA, then design primers."""
    fasta_path, _ = synthetic_fasta
    v = Variant('chr1', 5000, 'A', 'T')
    with ReferenceGenome(fasta_path) as ref:
        seq, offset = ref.fetch_around(v, flank=200)
    pairs = design_primers(seq, offset, len(v.ref), num_return=3)
    assert len(pairs) >= 1
