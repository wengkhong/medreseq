"""Write primer design results to a TSV file."""

from __future__ import annotations

import csv
from pathlib import Path

from .variants import Variant
from .primer_design import PrimerPair
from .insilico_pcr import PCRResult

HEADER = [
    'chrom', 'pos', 'ref', 'alt',
    'primer_rank',
    'left_primer', 'right_primer',
    'left_tm', 'right_tm',
    'left_gc_pct', 'right_gc_pct',
    'product_size',
    'primer3_penalty',
    'pcr_validated', 'pcr_product_size', 'pcr_off_targets', 'pcr_notes',
]


def write_results(
    output_path: Path,
    results: list[tuple[Variant, list[tuple[PrimerPair, PCRResult]]]],
) -> None:
    """Write all results to a tab-separated file.

    Each primer pair for each variant is one row. Variants with no primers get
    a single row with empty primer columns and a note.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w', newline='') as fh:
        writer = csv.writer(fh, delimiter='\t')
        writer.writerow(HEADER)
        for variant, pairs in results:
            if not pairs:
                writer.writerow([
                    variant.chrom, variant.pos, variant.ref, variant.alt,
                    '', '', '', '', '', '', '', '', '',
                    False, '', 0, 'no primers designed',
                ])
                continue
            for pp, pcr in pairs:
                writer.writerow([
                    variant.chrom, variant.pos, variant.ref, variant.alt,
                    pp.rank + 1,
                    pp.left_seq, pp.right_seq,
                    f'{pp.left_tm:.1f}', f'{pp.right_tm:.1f}',
                    f'{pp.left_gc:.1f}', f'{pp.right_gc:.1f}',
                    pp.product_size,
                    f'{pp.penalty:.3f}',
                    pcr.found,
                    pcr.product_size if pcr.product_size is not None else '',
                    pcr.off_target_count,
                    '; '.join(pcr.notes),
                ])
