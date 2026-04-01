"""Tests for VCF and TSV variant parsing."""

from __future__ import annotations

import textwrap
from pathlib import Path

import pytest

from medreseq.variants import Variant, parse_variants


# ── VCF parsing ───────────────────────────────────────────────────────────────

def test_parse_vcf_basic(tmp_path: Path) -> None:
    vcf = tmp_path / 'test.vcf'
    vcf.write_text(textwrap.dedent("""\
        ##fileformat=VCFv4.2
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
        chr1\t1000\t.\tA\tT\t100\tPASS\t.
        chr2\t2000\t.\tATCG\tA\t100\tPASS\t.
    """))
    variants = list(parse_variants(vcf))
    assert variants == [
        Variant('chr1', 1000, 'A', 'T'),
        Variant('chr2', 2000, 'ATCG', 'A'),
    ]


def test_parse_vcf_multiallelic_takes_first(tmp_path: Path) -> None:
    vcf = tmp_path / 'test.vcf'
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t1000\t.\tC\tG,T\t.\t.\t.\n"
    )
    variants = list(parse_variants(vcf))
    assert len(variants) == 1
    assert variants[0].alt == 'G'


def test_parse_vcf_skips_missing_alt(tmp_path: Path) -> None:
    vcf = tmp_path / 'test.vcf'
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t100\t.\tA\t.\t.\t.\t.\n"
        "chr1\t200\t.\tA\t*\t.\t.\t.\n"
        "chr1\t300\t.\tA\tT\t.\t.\t.\n"
    )
    variants = list(parse_variants(vcf))
    assert len(variants) == 1
    assert variants[0].pos == 300


def test_parse_vcf_ref_alt_uppercased(tmp_path: Path) -> None:
    vcf = tmp_path / 'test.vcf'
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t500\t.\ta\tt\t.\t.\t.\n"
    )
    v = list(parse_variants(vcf))[0]
    assert v.ref == 'A'
    assert v.alt == 'T'


# ── TSV parsing ───────────────────────────────────────────────────────────────

def test_parse_tsv_basic(tmp_path: Path) -> None:
    tsv = tmp_path / 'test.tsv'
    tsv.write_text("chrom\tpos\tref\talt\nchr1\t1000\tA\tT\nchr2\t2000\tATCG\tA\n")
    variants = list(parse_variants(tsv))
    assert variants == [
        Variant('chr1', 1000, 'A', 'T'),
        Variant('chr2', 2000, 'ATCG', 'A'),
    ]


def test_parse_tsv_flexible_column_names(tmp_path: Path) -> None:
    tsv = tmp_path / 'test.tsv'
    tsv.write_text("chr\tposition\treference\talternate\nchr1\t1000\tA\tT\n")
    variants = list(parse_variants(tsv))
    assert len(variants) == 1
    assert variants[0] == Variant('chr1', 1000, 'A', 'T')


def test_parse_tsv_missing_columns_raises(tmp_path: Path) -> None:
    tsv = tmp_path / 'test.tsv'
    tsv.write_text("chromosome\tposition\tsequence\nchr1\t1000\tA\n")
    with pytest.raises(ValueError, match="missing required columns"):
        list(parse_variants(tsv))


# ── Auto-detection ────────────────────────────────────────────────────────────

def test_autodetect_vcf_by_extension(tmp_path: Path) -> None:
    vcf = tmp_path / 'vars.vcf'
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t100\t.\tA\tT\t.\t.\t.\n"
    )
    assert list(parse_variants(vcf))[0].chrom == 'chr1'


def test_fixture_vcf(tmp_path: Path) -> None:
    """Smoke-test the bundled sample VCF fixture."""
    fixture = Path(__file__).parent / 'fixtures' / 'sample.vcf'
    variants = list(parse_variants(fixture))
    assert len(variants) == 3


def test_fixture_tsv(tmp_path: Path) -> None:
    """Smoke-test the bundled sample TSV fixture."""
    fixture = Path(__file__).parent / 'fixtures' / 'sample.tsv'
    variants = list(parse_variants(fixture))
    assert len(variants) == 3
