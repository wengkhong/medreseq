"""Command-line interface for MedReseq."""

from __future__ import annotations

import sys
from pathlib import Path

import click

from .variants import parse_variants
from .sequences import ReferenceGenome
from .primer_design import design_primers
from .insilico_pcr import run_insilico_pcr
from .output import write_results


@click.command(context_settings={'help_option_names': ['-h', '--help']})
@click.option(
    '--fasta', '-f', required=True,
    type=click.Path(exists=True, path_type=Path),
    help='Reference FASTA file (requires a .fai index alongside it).',
)
@click.option(
    '--input', '-i', 'input_file', required=True,
    type=click.Path(exists=True, path_type=Path),
    help='Variants file: VCF (.vcf) or TSV with chrom/pos/ref/alt columns.',
)
@click.option(
    '--output', '-o',
    type=click.Path(path_type=Path), default=Path('output.tsv'), show_default=True,
    help='Output TSV file.',
)
@click.option(
    '--flank', default=200, show_default=True,
    help='Flanking bases around each variant for primer design.',
)
@click.option(
    '--pcr-window', default=5000, show_default=True,
    help='Search window (±bp) around each variant for in silico PCR.',
)
@click.option(
    '--num-primers', default=3, show_default=True,
    help='Number of primer pairs to design per variant.',
)
@click.option(
    '--min-product', default=200, show_default=True,
    help='Minimum PCR product size (bp).',
)
@click.option(
    '--max-product', default=800, show_default=True,
    help='Maximum PCR product size (bp).',
)
@click.option('--verbose', '-v', is_flag=True, help='Print per-variant progress.')
def main(
    fasta: Path,
    input_file: Path,
    output: Path,
    flank: int,
    pcr_window: int,
    num_primers: int,
    min_product: int,
    max_product: int,
    verbose: bool,
) -> None:
    """Design Sanger sequencing validation primers for genomic variants.

    Accepts VCF files or TSV files with columns: chrom, pos, ref, alt.
    Each primer pair is validated by in silico PCR against the reference.

    \b
    Example:
        medreseq -f hg38.fa -i variants.vcf -o primers.tsv
        medreseq -f hg38.fa -i variants.tsv --num-primers 5 --min-product 150
    """
    try:
        ref = ReferenceGenome(fasta)
    except FileNotFoundError as exc:
        click.echo(f"ERROR: {exc}", err=True)
        sys.exit(1)

    all_results: list[tuple] = []

    with ref:
        for variant in parse_variants(input_file):
            if verbose:
                click.echo(f"Processing {variant}")

            # Primer design: fetch a window around the variant
            seq, offset = ref.fetch_around(variant, flank)
            pairs = design_primers(
                seq, offset, len(variant.ref),
                num_return=num_primers,
                product_size_range=(min_product, max_product),
            )

            if not pairs:
                click.echo(f"  WARNING: no primers designed for {variant}", err=True)
                all_results.append((variant, []))
                continue

            # In silico PCR: fetch a larger window for primer binding search
            pcr_seq, _ = ref.fetch_around(variant, pcr_window)
            pair_results = []
            for pair in pairs:
                pcr = run_insilico_pcr(pcr_seq, pair, min_product, max_product)
                if verbose:
                    if pcr.found:
                        status = (
                            f"PCR OK  product={pcr.product_size}bp"
                            + (f"  off-targets={pcr.off_target_count}" if pcr.off_target_count else "")
                        )
                    else:
                        status = "PCR FAIL  " + "; ".join(pcr.notes)
                    click.echo(
                        f"  rank{pair.rank + 1}  "
                        f"{pair.left_seq} / {pair.right_seq}  {status}"
                    )
                pair_results.append((pair, pcr))

            all_results.append((variant, pair_results))

    write_results(output, all_results)
    click.echo(f"Wrote {len(all_results)} variant(s) → {output}")
