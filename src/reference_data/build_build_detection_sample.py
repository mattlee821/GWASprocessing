#!/usr/bin/env python3
"""Build the paired GRCh37/GRCh38 SNP sample used to detect an input's genome build.

Usage:
    python3 src/reference_data/build_build_detection_sample.py --reference-root <root> [--size 2000000]

Takes the common biallelic SNPs of the GRCh37 1000 Genomes EUR lookup (<root>/1000genomes/phase3/lookup/EUR.tsv.gz,
with rsID, GRCh37 position and alleles), places each in GRCh38 through the dbSNP index (<root>/dbsnp/GRCh38), and
keeps those found in both builds: the same SNPs in each build, so neither build is favoured. A fixed-seed (821)
sample of --size SNPs is written to <root>/build_detection/ as grch37.{key.i64,ref.u1,alt.u1} and
grch38.{key.i64,ref.u1,alt.u1}, each sorted by key (chromosome code * 1e9 + position).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "pipeline" / "bin"))
import gwas_alleles as alleles  # noqa: E402
import gwas_reference as ref  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--reference-root", required=True)
    parser.add_argument("--size", type=int, default=2_000_000)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    root = Path(args.reference_root)
    out = root / "build_detection"
    if (out / "meta.json").exists() and not args.force:
        print(f"Skip: {out} already built")
        return
    dbsnp = ref.Dbsnp.find(root)
    if dbsnp is None:
        raise SystemExit(f"no dbSNP index under {root}/dbsnp/GRCh38; run build_dbsnp_index.py first")
    lookup, used = alleles.read_reference_lookup(root, "EUR")
    if not used:
        raise SystemExit(f"no 1000 Genomes GRCh37 lookup under {root}/1000genomes/phase3/lookup")
    lookup = lookup[(lookup["RS"].to_numpy() > 0) & np.isfinite(lookup["POS19"].to_numpy())]
    code38, pos38 = dbsnp.positions(lookup["RS"].to_numpy(dtype="int64"))
    code37 = ref.chrom_codes(lookup["CHR"].astype(str).to_numpy())
    both = (code38 > 0) & (code38 == code37)
    sample = np.where(both)[0]
    if len(sample) > args.size:
        sample = np.sort(np.random.default_rng(821).choice(sample, size=args.size, replace=False))
    alt = ref.base_codes(lookup["A1"].astype(str).to_numpy()[sample])  # A1 is ALT in the lookup, A2 is REF
    refb = ref.base_codes(lookup["A2"].astype(str).to_numpy()[sample])
    out.mkdir(parents=True, exist_ok=True)
    for name, keys in (("grch37", code37[sample] * ref.KEY_SCALE + lookup["POS19"].to_numpy(dtype="int64")[sample]),
                       ("grch38", code38[sample] * ref.KEY_SCALE + pos38[sample])):
        order = np.argsort(keys, kind="stable")
        keys[order].astype("<i8").tofile(out / f"{name}.key.i64")
        refb[order].astype("u1").tofile(out / f"{name}.ref.u1")
        alt[order].astype("u1").tofile(out / f"{name}.alt.u1")
    meta = {"snps": int(len(sample)), "source": f"1000 Genomes phase 3 {used} lookup (GRCh37) placed in GRCh38 by {dbsnp.label}",
            "seed": 821, "version": 1}
    (out / "meta.json").write_text(json.dumps(meta, indent=1), encoding="utf-8")
    print(json.dumps(meta))


if __name__ == "__main__":
    main()
