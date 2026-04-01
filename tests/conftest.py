"""Shared pytest fixtures."""

from __future__ import annotations

import random
from pathlib import Path

import pytest

from medreseq.variants import Variant


@pytest.fixture(scope='session')
def synthetic_fasta(tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, str]:
    """Write a 10 kb synthetic chr1 FASTA and index it with pysam.

    Uses a fixed random seed so the sequence is deterministic across runs.
    The sequence has ~50 % GC content and no long homopolymer runs, which
    gives primer3 a reasonable chance of finding primers anywhere in it.
    """
    import pysam

    tmp_dir = tmp_path_factory.mktemp('fixtures')
    rng = random.Random(42)
    seq = ''.join(rng.choices('ACGT', k=10_000))

    fasta_path = tmp_dir / 'synthetic.fa'
    with open(fasta_path, 'w') as fh:
        fh.write('>chr1\n')
        for i in range(0, len(seq), 60):
            fh.write(seq[i:i + 60] + '\n')

    pysam.faidx(str(fasta_path))
    return fasta_path, seq


@pytest.fixture(scope='session')
def snv_variant() -> Variant:
    """SNV at position 5000 on chr1 of the synthetic reference."""
    return Variant(chrom='chr1', pos=5000, ref='A', alt='T')


@pytest.fixture(scope='session')
def del_variant() -> Variant:
    """2-bp deletion at position 5000 on chr1 of the synthetic reference."""
    return Variant(chrom='chr1', pos=5000, ref='AC', alt='A')
