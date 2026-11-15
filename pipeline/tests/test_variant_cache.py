"""The variant cache (gwas_cache.py): a GWAS standardised with a cold, a warm or no cache gives the same records, the
second GWAS of a variant set reads its alignment and indel results from the cache, and a damaged cache is only a miss."""

import os
import tempfile
import unittest
from pathlib import Path

from testlib import example_reference, generic_grch38_table, make_reference, run_standardise, write_table
from test_step09_alignment import HEADER, sequence

import gwas_cache


def records(result) -> list[str]:
    return [line for line in result.text.splitlines() if not line.startswith("##")]


class VariantCacheTest(unittest.TestCase):
    def test_variant_cache(self):
        tmp = Path(tempfile.mkdtemp())
        genome = sequence(3000, {100: "G", 200: "G", 300: "G", 400: "G", 500: "G", 1001: "GCACACATT", 2001: "GCACACATT",
                                 2501: "TAG", 2700: "CCC", 600: "N", 700: "GNN"})
        reference = make_reference(tmp / "reference", lookup_rows=[], fasta={"1": genome})
        rows = [[1, 100, "keep", "A", "G", 0.3, 0.1, 0.02, 1e-6, 5000], [1, 200, "swap", "G", "A", 0.3, 0.2, 0.02, 1e-6, 5000],
                [1, 300, "comp", "T", "C", 0.3, 0.3, 0.02, 1e-6, 5000], [1, 400, "compswap", "C", "T", 0.3, 0.4, 0.02, 1e-6, 5000],
                [1, 500, "neither", "A", "T", 0.3, 0.5, 0.02, 1e-6, 5000], [1, 1005, "del", "A", "ACA", 0.3, 0.6, 0.02, 1e-6, 5000],
                [1, 2003, "ins", "ACAC", "AC", 0.3, 0.7, 0.02, 1e-6, 5000], [1, 2502, "indelswap", "A", "AT", 0.3, 0.8, 0.02, 1e-6, 5000],
                [1, 2700, "indelbad", "A", "GTT", 0.3, 0.9, 0.02, 1e-6, 5000], [1, 5000, "beyond", "A", "G", 0.3, 1.0, 0.02, 1e-6, 5000],
                [2, 100, "nocontig", "A", "G", 0.3, 1.1, 0.02, 1e-6, 5000]]
        raw = write_table(tmp / "gwas.tsv", HEADER, rows)
        other = write_table(tmp / "other.tsv", HEADER, [[*r[:6], r[6] * 2, *r[7:]] for r in rows])  # same variants, other statistics
        cache_dir = tmp / "cache"

        def run(name, source, cache, chunk="100000"):
            env = {**os.environ, "GWAS_VARIANT_CACHE": str(cache), "GWAS_TABLE_CHUNK": chunk}
            result = run_standardise(tmp, source, reference, name=name, env=env, force=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            return result

        off = run("off", raw, "off")
        cold = run("cold", raw, cache_dir)
        warm = run("warm", raw, cache_dir)
        pieces = run("pieces", raw, tmp / "cache_pieces", chunk="3")
        pieces_warm = run("pieces_warm", raw, tmp / "cache_pieces", chunk="4")

        with self.subTest("cold, warm and no cache write the same records"):
            self.assertEqual(records(cold), records(off))
            self.assertEqual(records(warm), records(off))
            self.assertEqual(records(pieces), records(off))
            self.assertEqual(records(pieces_warm), records(off))
            self.assertEqual(cold.qc["reference_alignment"], off.qc["reference_alignment"])
            self.assertEqual(cold.qc["indels_left_aligned"], off.qc["indels_left_aligned"])
            self.assertEqual(warm.qc["indels_left_aligned"], off.qc["indels_left_aligned"])

        with self.subTest("the first GWAS fills the cache and the second reads it"):
            self.assertNotIn("variant_cache", off.qc)
            self.assertEqual(cold.qc["variant_cache"]["alignment"], {"hits": 0, "misses": len(rows)})
            self.assertEqual(warm.qc["variant_cache"]["alignment"], {"hits": len(rows), "misses": 0})
            self.assertEqual(warm.qc["variant_cache"]["indels"]["misses"], 0)
            self.assertGreater(warm.qc["variant_cache"]["indels"]["hits"], 0)

        with self.subTest("another GWAS of the same variants takes the cached results but its own statistics"):
            fresh = run("other_off", other, "off")
            reused = run("other_cached", other, cache_dir)
            self.assertEqual(reused.qc["variant_cache"]["alignment"]["misses"], 0)
            self.assertEqual(records(reused), records(fresh))
            self.assertNotEqual(records(reused), records(warm))

        with self.subTest("the 1000 Genomes annotation is cached too (frequencies fill EAF, concordance is unchanged)"):
            example = example_reference()
            table = generic_grch38_table(tmp / "generic.tsv")
            runs = {}
            for label, cache in (("a_off", "off"), ("a_cold", tmp / "cache_a"), ("a_warm", tmp / "cache_a")):
                env = {**os.environ, "GWAS_VARIANT_CACHE": str(cache)}
                runs[label] = run_standardise(tmp, table, example, name=label, env=env, force=True, population="ALL")
                self.assertEqual(runs[label].returncode, 0, runs[label].stderr)
            self.assertEqual(records(runs["a_cold"]), records(runs["a_off"]))
            self.assertEqual(records(runs["a_warm"]), records(runs["a_off"]))
            self.assertEqual(runs["a_warm"].qc["variant_cache"]["annotation"], {"hits": 2, "misses": 0})
            self.assertEqual(runs["a_cold"].qc["af_concordance"], runs["a_off"].qc["af_concordance"])
            self.assertEqual(runs["a_cold"].qc["reference_fills"], runs["a_off"].qc["reference_fills"])

        with self.subTest("a damaged or vanished segment is a miss, not an error"):
            segments = sorted(cache_dir.rglob("seg-*.arrow"))
            self.assertTrue(segments)
            segments[0].write_bytes(b"not an arrow file")
            segments[-1].unlink()
            damaged = run("damaged", raw, cache_dir)
            self.assertEqual(records(damaged), records(off))
            self.assertGreater(damaged.qc["variant_cache"]["alignment"]["misses"], 0)

        with self.subTest("compaction merges a bucket's segments and keeps every key"):
            caches = sorted(p for p in cache_dir.iterdir() if p.is_dir())
            before = sum(len(list(p.rglob("seg-*.arrow"))) for p in caches)
            self.assertGreater(before, 0)
            for directory in caches:
                for bucket in directory.iterdir():
                    if bucket.is_dir():  # copy a segment to push the bucket above the threshold
                        for segment in list(bucket.glob("seg-*.arrow")):
                            for k in range(10):
                                (bucket / f"seg-{k}-0-0.arrow").write_bytes(segment.read_bytes())
                            break
            compacted = 0
            for directory in caches:
                cache = gwas_cache.VariantCache(directory, {})
                compacted += cache.compact()
            self.assertGreater(compacted, 0)
            again = run("after_compaction", raw, cache_dir)
            self.assertEqual(records(again), records(off))
            self.assertEqual(again.qc["variant_cache"]["alignment"]["misses"], 0)


if __name__ == "__main__":
    unittest.main()
