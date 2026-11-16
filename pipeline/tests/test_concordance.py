"""concordance.py: agreeing studies give ok, a flipped sign or swapped alleles are caught, a swapped REF/ALT with the sign
changed to match is still ok."""
import gzip
import tempfile
import unittest
from pathlib import Path

import numpy as np

import testlib  # noqa: F401  (puts pipeline/bin on the path)
import concordance


def write_vcf(path: Path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as handle:
        handle.write("##fileformat=VCFv4.2\n##contig=<ID=1>\n" + "".join(f'##FORMAT=<ID={k},Number=A,Type=Float,Description="{k}">\n' for k in ("ES", "SE", "LP", "AF", "AFS")) + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tx\n")
        for chrom, pos, ref, alt, es, se, lp, af in rows:
            handle.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tES:SE:LP:AF:AFS\t{es}:{se}:{lp}:{af}:1\n")


def variants(n=200):
    rng = np.random.default_rng(1)
    out = []
    for i in range(n):
        es = rng.normal(0, 0.3) if i > 5 else 0.8 * (1 if i % 2 else -1)
        se = 0.05
        lp = min(50.0, abs(es / se) ** 2 / 4.6)
        out.append((1, 1000 + i, "A", "C" if i % 3 else "G", es, se, lp, round(rng.uniform(0.1, 0.9), 3)))
    return out


class ConcordanceTest(unittest.TestCase):
    def test_concordance(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            base = variants()
            write_vcf(tmp / "a/standardise/P1_1-1/P1_1-1.vcf.gz", base)
            write_vcf(tmp / "same/standardise/P1_1-1/P1_1-1.vcf.gz", [(c, p, r, a, e * 1.1, s, l, f) for c, p, r, a, e, s, l, f in base])
            write_vcf(tmp / "flipped/standardise/P1_1-1/P1_1-1.vcf.gz", [(c, p, r, a, -e, s, l, f) for c, p, r, a, e, s, l, f in base])
            write_vcf(tmp / "swapped/standardise/P1_1-1/P1_1-1.vcf.gz", [(c, p, a, r, -e, s, l, 1 - f) for c, p, r, a, e, s, l, f in base])
            studies = {name: tmp / name for name in ("a", "same", "flipped", "swapped")}
            found = concordance.find_files(studies)
            self.assertEqual(list(found), ["P1_1-1"])
            verdict = {(r["study_a"], r["study_b"]): r["verdict"] for r in concordance.protein_pairs(("P1_1-1", found["P1_1-1"]))}
            self.assertEqual(verdict[("a", "same")], "ok")
            self.assertEqual(verdict[("a", "swapped")], "ok")  # REF/ALT swapped and the sign changed: the same effect
            self.assertEqual(verdict[("a", "flipped")], "SIGN OR ALLELE FLIPPED")
            self.assertEqual(verdict[("flipped", "same")], "SIGN OR ALLELE FLIPPED")


if __name__ == "__main__":
    unittest.main()
