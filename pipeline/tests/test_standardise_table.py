"""End-to-end tests for the table path of standardise_gwas.py (summary-statistics tables -> GWAS-VCF).

Each test runs the real script on a small input. Expected values are hand-derived; derivation
literals were computed in R (log(1.2); (log(1.3091)-log(1.1))/(2*qnorm(0.975)); 0.1/qnorm(1-1e-6/2)).
"""

import os
import tempfile
import unittest
from pathlib import Path

import numpy as np

from testlib import EXAMPLES, ROOT, chain_block, make_reference, run_standardise, write_chain, write_table

HEADER = ["CHR", "POS", "SNP", "EA", "OA", "EAF", "BETA", "SE", "P", "N"]


class AdaptedBehaviourTests(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def test_positional_match_supplies_rsid_and_effect_allele_frequency(self):
        result = run_standardise(self.tmp, EXAMPLES / "tables/snpid_missing_eaf.tsv", EXAMPLES / "reference",
                                 ancestry="all", input_build="hg19")
        self.assertEqual(result.returncode, 0, result.stderr)
        rows = result.rows
        self.assertEqual([r["ID"] for r in rows], ["rs3", "."])
        self.assertEqual((rows[0]["FORMAT"]["AF"], rows[0]["FORMAT"]["AFS"][:6]), ("0.33", "1000G_"))
        self.assertEqual(rows[1]["FORMAT"]["AF"], ".")
        self.assertEqual((result.qc["reference_fills"]["rsid"], result.qc["reference_fills"]["EAF"]), (1, 1))

    def test_rsid_only_table_gets_coordinates_from_the_reference(self):
        result = run_standardise(self.tmp, EXAMPLES / "tables/locke_direct.tsv", EXAMPLES / "reference", sex="female")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.qc["source_columns"]["EA"], "A1")
        self.assertEqual(result.qc["source_columns"]["OA"], "A2")
        self.assertEqual(result.qc["source_columns"]["EAF"], "Freq1.Hapmap")
        self.assertEqual((result.qc["reference_fills"]["CHR"], result.qc["reference_fills"]["POS19"]), (2, 2))
        self.assertEqual(result.qc["coordinate_source"], "population_lookup")
        self.assertEqual((result.metadata["sex"], result.metadata["coordinate_build"]), ("female", "GRCh38"))
        self.assertEqual(result.sample["Sex"], "female")
        self.assertIn("##reference=GRCh38", result.header)
        rs1 = result.by_id()["rs1"]
        self.assertEqual((rs1["CHROM"], rs1["POS"], rs1["REF"], rs1["ALT"]), ("1", 1100, "G", "A"))

    def test_plink_table_derives_other_allele_from_ref_and_alt(self):
        result = run_standardise(self.tmp, EXAMPLES / "tables/zhang_plink.tsv", EXAMPLES / "reference",
                                 config_root=ROOT / "manifests", phenotype="protein", ancestry="EUR",
                                 author="Zhang", year="2022", PMID="35501419", study_id="zhang_2022_35501419")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(result.qc["derived_oa_from_ref_alt"])
        rows = result.by_id()
        self.assertEqual((rows["rs5"]["POS"], rows["rs5"]["REF"], rows["rs5"]["ALT"]), (5100, "G", "A"))
        self.assertEqual((rows["rs6"]["REF"], rows["rs6"]["ALT"]), ("C", "T"))
        self.assertEqual(rows["rs5"]["INFO"]["ERRCODE"], ".")
        self.assertEqual(result.qc["build_detection"]["build"], "hg38")

    def test_missing_ancestry_fails_loudly(self):
        result = run_standardise(self.tmp, EXAMPLES / "tables/generic.tsv", EXAMPLES / "reference", ancestry="")
        self.assertEqual(result.returncode, 78)
        self.assertIn("ancestry", result.stderr.lower())

    def test_standard_steps_keep_rare_variants_and_indels_and_remove_only_duplicates(self):
        raw = write_table(self.tmp / "qc.tsv", ["CHR", "POS38", "SNP", "EA", "OA", "EAF", "BETA", "SE", "P", "N"], [
            [1, 1000, "rs1", "A", "G", 0.30, 0.1, 0.01, 0.2, 1000],
            [1, 2000, "rs2", "AT", "G", 0.30, 0.1, 0.01, 0.2, 1000],
            [1, 3000, "rs3", "A", "T", 0.30, 0.1, 0.01, 0.2, 1000],
            [1, 4000, "rs4", "A", "G", 0.005, 0.1, 0.01, 0.2, 1000],
            [1, 1000, "rs5", "A", "G", 0.30, 0.1, 0.01, 0.2, 1000],
            [1, 5000, "rs6", "D", "I", 0.30, 0.1, 0.01, 0.2, 1000],
        ])
        result = run_standardise(self.tmp, raw, EXAMPLES / "reference")
        self.assertEqual(sorted(r["ID"] for r in result.rows), ["rs1", "rs2", "rs3", "rs4"])
        self.assertEqual(result.qc["duplicates_removed"], 1)
        self.assertEqual(result.qc["incomplete_by_cause"]["invalid_allele"]["rows"], 1)
        self.assertTrue(result.by_id()["rs3"]["INFO"]["HM"].startswith("palindrome_"))


class ReviewRegressionTests(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def reference(self, rows, **kwargs):
        return make_reference(self.tmp / "ref", lookup_rows=rows, **kwargs)

    def test_empty_tab_cell_keeps_later_columns_in_place(self):
        raw = write_table(self.tmp / "gwas.tsv", HEADER, [[1, 3000, "rs30", "A", "G", None, 0.10, 0.02, 1e-6, 5000],
                                                          [1, 3100, "rs31", "C", "T", 0.4, 0.05, 0.01, 0.2, 5000]])
        result = run_standardise(self.tmp, raw, self.reference([]), input_build="hg19")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.by_id()["rs30"]["FORMAT"], {"ES": "0.1", "SE": "0.02", "LP": "6", "AF": ".", "AFS": ".", "SS": "5000",
                                                          "ID": "rs30"})

    def test_filled_frequency_is_that_of_the_effect_allele(self):
        ref = self.reference([("1", "rs10", 1000, "A", "G", 0.8), ("1", "rs11", 1100, "A", "G", 0.2)])
        raw = write_table(self.tmp / "gwas.tsv", ["CHR", "POS", "SNP", "EA", "OA", "BETA", "SE", "P", "N"],
                          [[1, 1000, "rs10", "A", "G", 0.1, 0.02, 1e-6, 5000], [1, 1100, "rs11", "G", "A", 0.1, 0.02, 1e-6, 5000]])
        rows = run_standardise(self.tmp, raw, ref, input_build="hg19").by_id()
        self.assertEqual((rows["rs10"]["FORMAT"]["AF"], rows["rs10"]["FORMAT"]["AFS"][:6]), ("0.8", "1000G_"))
        self.assertEqual(rows["rs11"]["FORMAT"]["AF"], "0.8")

    def test_palindrome_with_reference_filled_frequency_is_kept_ambiguous_not_flipped(self):
        ref = self.reference([("1", "rs20", 2000, "A", "T", 0.8), ("1", "rs21", 2100, "A", "G", 0.3)])
        raw = write_table(self.tmp / "gwas.tsv", ["CHR", "POS", "SNP", "EA", "OA", "BETA", "SE", "P", "N"],
                          [[1, 2000, "rs20", "A", "T", 0.1, 0.02, 1e-6, 5000], [1, 2100, "rs21", "A", "G", 0.1, 0.02, 1e-6, 5000]])
        result = run_standardise(self.tmp, raw, ref, input_build="hg19")
        self.assertEqual(list(result.by_id()), ["rs20", "rs21"])
        self.assertTrue(result.by_id()["rs20"]["INFO"]["HM"].startswith("palindrome_ambiguous"))
        self.assertEqual(result.by_id()["rs20"]["FORMAT"]["AFS"][:6], "1000G_")
        detail = result.qc["palindromic_resolution"]
        self.assertEqual((detail["ambiguous_no_study_eaf"], detail["flipped"]), (1, 0))

    def detection_reference(self):
        return self.reference([("1", f"rs6{i:02d}", 6000 + 1000 * i, "A", "G", 0.3) for i in range(30)])

    def detection_table(self, positions, ids=None):
        ids = ids or [f"1:{p}:A:G" for p in positions]
        return write_table(self.tmp / "gwas.tsv", ["CHR", "POS", "MarkerName", "EA", "OA", "EAF", "BETA", "SE", "P", "N"],
                           [[1, p, i, "A", "G", 0.3, 0.1, 0.02, 1e-6, 5000] for p, i in zip(positions, ids)])

    def test_undeclared_build_is_detected_from_positions(self):
        result = run_standardise(self.tmp, self.detection_table([6000 + 1000 * i for i in range(30)]), self.detection_reference())
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual((result.qc["build_detection"]["build"], result.qc["build_detection"]["method"]), ("hg19", "detected"))
        first = result.rows[0]
        self.assertEqual((first["POS"], first["ID"]), (6100, "rs600"))

    def test_build_that_cannot_be_determined_fails_loudly(self):
        result = run_standardise(self.tmp, self.detection_table([123_457 + 7 * i for i in range(30)]), self.detection_reference())
        self.assertEqual(result.returncode, 78)
        self.assertIn("input_build", result.stderr)

    def test_build_detected_from_the_data_wins_over_a_contradicting_declaration(self):
        result = run_standardise(self.tmp, self.detection_table([6000 + 1000 * i for i in range(30)]), self.detection_reference(),
                                 input_build="hg38")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual((result.qc["build_detection"]["build"], result.qc["build_detection"]["method"]), ("hg19", "detected"))
        self.assertTrue(any("disagrees" in w for w in result.qc["warnings"]))

    def test_grch38_chr_pos_ids_never_borrow_grch37_annotations(self):
        # 30 genuine GRCh38 variants plus 3 whose GRCh38 position coincides with an unrelated GRCh37 variant.
        true_rows = [("1", f"rs{9000 + i}", 9000 + 1000 * i, "A", "G", 0.3) for i in range(30)]
        decoys = [("1", f"rs{7000 + j}", 200000 + 1000 * j, "A", "G", 0.3) for j in range(3)]
        ref = self.reference(true_rows + decoys)
        positions = [9100 + 1000 * i for i in range(30)] + [200000 + 1000 * j for j in range(3)]
        result = run_standardise(self.tmp, self.detection_table(positions, ids=[f"1:{p}" for p in positions]), ref, input_build="hg38")
        self.assertEqual(result.returncode, 0, result.stderr)
        by_pos = {r["POS"]: r for r in result.rows}
        self.assertEqual(by_pos[200000]["ID"], ".")
        self.assertEqual(by_pos[9100]["ID"], "rs9000")

    def test_p_value_of_zero_is_recovered_from_the_z_score(self):
        raw = write_table(self.tmp / "gwas.tsv", HEADER, [[1, 7000, "rs70", "A", "G", 0.3, 0.5, 0.01, 0, 5000]])
        rows = run_standardise(self.tmp, raw, self.reference([]), input_build="hg19").rows
        self.assertEqual(rows[0]["FORMAT"]["LP"], "544.66531")

    def test_log10p_column_keeps_values_beyond_double_precision(self):
        raw = write_table(self.tmp / "gwas.tsv", ["CHR", "POS", "SNP", "EA", "OA", "EAF", "BETA", "SE", "LOG10P", "N"],
                          [[1, 7000, "rs70", "A", "G", 0.3, 0.5, 0.01, 450, 5000], [1, 7100, "rs71", "A", "G", 0.3, 0.1, 0.02, 12.5, 5000]])
        rows = run_standardise(self.tmp, raw, self.reference([]), input_build="hg19").by_id()
        self.assertEqual((rows["rs70"]["FORMAT"]["LP"], rows["rs71"]["FORMAT"]["LP"]), ("450", "12.5"))

    def test_minus_strand_liftover_complements_and_other_chromosome_is_dropped(self):
        root = self.tmp / "ref"
        write_chain(root / "genome" / "hg19ToHg38.over.chain.gz", [
            chain_block("chr1", "chr2", 249250621, 242193529, 5000, 1000, 90000, 1),
            chain_block("chr1", "chr1", 249250621, 248956422, 7000, 1000, 10000, 2, dst_strand="-"),
        ])
        raw = write_table(self.tmp / "gwas.tsv", HEADER, [[1, 5500, "rs50", "A", "G", 0.3, 0.1, 0.02, 1e-6, 5000],
                                                          [1, 7500, "rs51", "A", "C", 0.3, 0.1, 0.02, 1e-6, 5000]])
        result = run_standardise(self.tmp, raw, root, input_build="hg19")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(list(result.by_id()), ["rs51"])
        rs51 = result.by_id()["rs51"]
        self.assertEqual((rs51["POS"], rs51["REF"], rs51["ALT"]), (248945923, "G", "T"))
        self.assertEqual(result.qc["liftover"]["hg19_to_hg38"]["other_chrom"], 1)

    def test_exclusions_without_joins_keep_rows_aligned(self):
        ref = self.reference([("1", "rs40", 4000, "A", "G", 0.3), ("1", "rs41", 4100, "C", "T", 0.4), ("1", "rs42", 4200, "A", "C", 0.2)])
        manifests = self.tmp / "manifests"
        manifests.mkdir()
        (manifests / "exclude.tsv").write_text("SNP\nrs40\n", encoding="utf-8")
        (manifests / "test_2026_1.yaml").write_text("exclude:\n  - file: exclude.tsv\n    key: SNP\n", encoding="utf-8")
        raw = write_table(self.tmp / "gwas.tsv", HEADER, [[1, 4000, "rs40", "A", "G", 0.3, 0.1, 0.02, 1e-6, 5000],
                                                          [1, 4100, "1:4100", "C", "T", 0.4, 0.05, 0.01, 0.2, 5000],
                                                          [1, 4200, "rs42", "A", "C", 0.2, 0.01, 0.01, 0.5, 5000]])
        result = run_standardise(self.tmp, raw, ref, config_root=manifests, input_build="hg19")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(sorted(result.by_id()), ["rs41", "rs42"])
        self.assertEqual(result.by_id()["rs41"]["FORMAT"]["ES"], "0.05")
        self.assertEqual(result.qc["exclusions"][0]["removed"], 1)

    def test_no_usable_rows_is_a_data_error(self):
        # "Estimate" is not a recognised effect-size column, so no row has BETA.
        raw = write_table(self.tmp / "gwas.tsv", ["CHR", "POS", "SNP", "EA", "OA", "EAF", "Estimate", "SE", "P", "N"],
                          [[1, 9000, "rs90", "A", "G", 0.3, 0.1, 0.02, 1e-6, 5000]])
        result = run_standardise(self.tmp, raw, self.reference([]), input_build="hg19")
        self.assertEqual(result.returncode, 65)
        self.assertIn("BETA", result.stderr)

    def test_low_retention_is_reported_as_a_warning(self):
        raw = write_table(self.tmp / "gwas.tsv", HEADER, [[1, 9000 + i, f"rs9{i}", "A", "G", 0.3, 0.1 if i == 0 else None, 0.02, 1e-6, 5000]
                                                          for i in range(4)])
        result = run_standardise(self.tmp, raw, self.reference([]), input_build="hg19")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(any("retained 1 of 4" in w for w in result.qc["warnings"]), result.qc["warnings"])
        self.assertEqual(result.qc["missing_required_by_field"]["BETA"], 3)

    def test_reruns_follow_manifest_and_setting_changes(self):
        ref = self.reference([("1", "rs80", 8000, "A", "G", 0.3), ("1", "rs81", 8100, "A", "T", 0.3)])
        raw = write_table(self.tmp / "gwas.tsv", HEADER, [[1, 8000, "rs80", "A", "G", 0.3, 0.1, 0.02, 1e-6, 5000],
                                                          [1, 8100, "rs81", "A", "T", 0.45, 0.1, 0.02, 1e-6, 5000]])
        first = run_standardise(self.tmp, raw, ref, input_build="hg19", ancestry="EUR", name="row")
        self.assertEqual(first.sample["Ancestry"], "EUR")
        second = run_standardise(self.tmp, raw, ref, input_build="hg19", ancestry="EAS", phenotype="Body mass index", name="row")
        self.assertEqual((second.sample["Ancestry"], second.sample["Phenotype"]), ("EAS", "Body mass index"))
        self.assertFalse(second.record.get("reused_existing_outputs", False))
        third = run_standardise(self.tmp, raw, ref, input_build="hg19", ancestry="EAS", phenotype="Body mass index", name="row")
        self.assertTrue(third.record["reused_existing_outputs"])
        self.assertEqual(third.qc["completed_at"], second.qc["completed_at"])

    def test_unusable_raw_gwas_catalog_file_falls_back_to_the_harmonised_file(self):
        from testlib import LocalHttpServer
        served = self.tmp / "served"
        served.mkdir()
        rows = [[1, 1100 + i, f"rs{i}", "A", "G", 0.3, 0.1, 0.02, 1e-6] for i in range(3)]
        write_table(served / "x.h.tsv", ["chromosome", "base_pair_location", "rsid", "effect_allele", "other_allele",
                                         "effect_allele_frequency", "beta", "standard_error", "p_value"], rows)
        raw = write_table(self.tmp / "raw.tsv", ["SNP", "A1", "A2", "P"], [["rs1", "A", "G", 0.1]])
        with LocalHttpServer(served) as server:
            metadata = {"source_type": "gwas_catalog", "file_kind": "raw", "harmonised_url": f"{server.url}/x.h.tsv",
                        "genome_assembly": "GRCh38"}
            result = run_standardise(self.tmp, raw, EXAMPLES / "reference", metadata=metadata)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.qc["gwas_catalog_fallback"]["harmonised_file"], "x.h.tsv")
        self.assertEqual(len(result.rows), 3)
        self.assertIn("harmonised file x.h.tsv was used", result.stderr)

    def test_id_coded_indels_with_grch37_positions_count_as_invalid_alleles_not_liftover_failures(self):
        raw = write_table(self.tmp / "id.tsv", HEADER, [[1, 1000 + i, f"rs{i}", "A", "G", 0.3, 0.1, 0.02, 1e-6, 5000] for i in range(3)]
                          + [[1, 2000, "chr1:2000:D", "D", "I", 0.3, 0.1, 0.02, 1e-6, 5000]])
        result = run_standardise(self.tmp, raw, self.reference([]), input_build="hg19")
        self.assertEqual(result.returncode, 0, result.stderr)
        causes = result.qc["incomplete_by_cause"]
        self.assertEqual(causes["invalid_allele"]["rows"], 1)
        self.assertNotIn("liftover_to_grch38_failed", causes)

    def test_yengo_style_frequency_column_is_read_as_the_study_eaf(self):
        raw = write_table(self.tmp / "yengo.tsv", ["CHR", "POS", "SNP", "Tested_Allele", "Other_Allele", "Freq_Tested_Allele_in_HRS",
                                                   "BETA", "SE", "P", "N"], [[1, 1000, "rs1", "A", "G", 0.31, 0.1, 0.02, 1e-6, 5000]])
        result = run_standardise(self.tmp, raw, self.reference([]), input_build="hg19")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.qc["source_columns"]["EAF"], "Freq_Tested_Allele_in_HRS")
        self.assertEqual(result.rows[0]["FORMAT"]["AFS"], "study")

    def test_odds_ratios_and_their_confidence_interval_are_converted(self):
        raw = write_table(self.tmp / "gwas.tsv", ["CHR", "POS", "SNP", "EA", "OA", "EAF", "OR", "L95", "U95", "P", "N"],
                          [[1, 7000, "rs70", "A", "G", 0.3, 1.2, 1.1, 1.3091, 1e-4, 5000]])
        result = run_standardise(self.tmp, raw, self.reference([]), input_build="hg19")
        self.assertEqual(result.returncode, 0, result.stderr)
        fmt = result.rows[0]["FORMAT"]
        self.assertAlmostEqual(float(fmt["ES"]), 0.182321556794, places=7)
        self.assertAlmostEqual(float(fmt["SE"]), 0.044396147014, places=7)
        self.assertEqual((result.qc["derivations"]["beta_from_or"], result.qc["derivations"]["se_from_ci"]), (1, 1))

    def test_missing_standard_error_is_derived_from_p(self):
        raw = write_table(self.tmp / "gwas.tsv", ["CHR", "POS", "SNP", "EA", "OA", "EAF", "BETA", "P", "N"],
                          [[1, 7000, "rs70", "A", "G", 0.3, 0.1, 1e-6, 5000]])
        result = run_standardise(self.tmp, raw, self.reference([]), input_build="hg19")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertAlmostEqual(float(result.rows[0]["FORMAT"]["SE"]), 0.0204430479678, places=8)
        self.assertEqual(result.qc["derivations"]["se_from_p"], 1)

    def test_harmonised_and_original_allele_columns_are_never_mixed(self):
        # The original beta is missing, so the harmonised family (alleles + beta) must be used together.
        raw = write_table(self.tmp / "gwas.tsv",
                          ["chromosome", "base_pair_location", "variant_id", "effect_allele", "other_allele", "beta",
                           "hm_effect_allele", "hm_other_allele", "hm_beta", "standard_error", "p_value"],
                          [[1, 1000, "rs1", "A", "G", None, "G", "A", -0.1, 0.02, 1e-6]])
        result = run_standardise(self.tmp, raw, self.reference([]), input_build="hg19")
        self.assertEqual(result.returncode, 0, result.stderr)
        row = result.rows[0]
        self.assertEqual((row["REF"], row["ALT"], row["FORMAT"]["ES"]), ("A", "G", "-0.1"))
        self.assertEqual(result.qc["source_columns"]["BETA"], "hm_beta")

    def test_missing_bgzip_is_a_configuration_error(self):
        raw = write_table(self.tmp / "gwas.tsv", HEADER, [[1, 1100, "rs11", "A", "G", 0.3, 0.1, 0.02, 1e-6, 5000]])
        env = {**os.environ, "PATH": os.pathsep.join(p for p in os.environ["PATH"].split(os.pathsep) if not Path(p, "bgzip").exists())}
        result = run_standardise(self.tmp, raw, self.reference([]), input_build="hg19", env=env)
        self.assertEqual(result.returncode, 78)
        self.assertIn("bgzip", result.stderr)

    def test_effect_allele_frequencies_anticorrelated_with_reference_fail(self):
        rng = np.random.default_rng(821)
        freqs = rng.uniform(0.05, 0.95, 1200)
        ref = self.reference([("1", f"rs{i}", 10000 + 10 * i, "A", "G", round(f, 4)) for i, f in enumerate(freqs)])
        # Effect allele labelled G but the frequency column describes A: the classic swapped-columns error.
        raw = write_table(self.tmp / "gwas.tsv", HEADER, [[1, 10000 + 10 * i, f"rs{i}", "G", "A", round(f, 4), 0.1, 0.02, 0.5, 5000]
                                                          for i, f in enumerate(freqs)])
        result = run_standardise(self.tmp, raw, ref, input_build="hg19")
        self.assertEqual(result.returncode, 65)
        self.assertIn("anti-correlate", result.stderr)

    def test_alleles_are_oriented_to_the_grch38_reference_when_a_fasta_is_present(self):
        sequence = ["N"] * 2500
        sequence[1100 - 1], sequence[1200 - 1] = "G", "A"
        ref = self.reference([], fasta={"1": "".join(sequence)})
        raw = write_table(self.tmp / "gwas.tsv", HEADER, [[1, 1000, "rs1", "A", "G", 0.3, 0.1, 0.02, 1e-6, 5000],
                                                          [1, 1100, "rs2", "A", "G", 0.3, 0.1, 0.02, 1e-6, 5000]])
        result = run_standardise(self.tmp, raw, ref, input_build="hg19")
        self.assertEqual(result.returncode, 0, result.stderr)
        rows = result.by_id()
        self.assertEqual((rows["rs1"]["REF"], rows["rs1"]["ALT"], rows["rs1"]["FORMAT"]["ES"]), ("G", "A", "0.1"))
        self.assertEqual((rows["rs2"]["REF"], rows["rs2"]["ALT"], rows["rs2"]["FORMAT"]["ES"], rows["rs2"]["FORMAT"]["AF"]),
                         ("A", "G", "-0.1", "0.7"))
        self.assertEqual((result.qc["reference_alignment"]["performed"], result.qc["reference_alignment"]["swapped"]), (True, 1))

    def test_reference_alignment_is_reported_as_skipped_without_a_fasta(self):
        raw = write_table(self.tmp / "gwas.tsv", HEADER, [[1, 1000, "rs1", "A", "G", 0.3, 0.1, 0.02, 1e-6, 5000]])
        result = run_standardise(self.tmp, raw, self.reference([]), input_build="hg19")
        self.assertFalse(result.qc["reference_alignment"]["performed"])

    def test_non_primary_contigs_are_dropped(self):
        raw = write_table(self.tmp / "gwas.tsv", HEADER, [[1, 1000, "rs1", "A", "G", 0.3, 0.1, 0.02, 1e-6, 5000],
                                                          ["6_ssto_hap7", 1000, "rs2", "A", "G", 0.3, 0.1, 0.02, 1e-6, 5000]])
        result = run_standardise(self.tmp, raw, self.reference([]), input_build="hg19")
        self.assertEqual(list(result.by_id()), ["rs1"])
        self.assertEqual(result.qc["missing_required_by_field"]["CHR"], 1)

    def test_no_intermediate_files_are_left_in_the_output_or_temporary_folders(self):
        scratch = self.tmp / "scratch"
        scratch.mkdir()
        raw = write_table(self.tmp / "gwas.tsv", HEADER, [[1, 1000, "rs1", "A", "G", 0.3, 0.1, 0.02, 1e-6, 5000]])
        result = run_standardise(self.tmp, raw, self.reference([]), input_build="hg19", env={**os.environ, "GWAS_TMPDIR": str(scratch)})
        self.assertEqual(result.returncode, 0, result.stderr)
        outputs = sorted(p.name for p in Path(result.record["standard_vcf"]).parent.iterdir())
        self.assertEqual(outputs, ["standardised.gwas.vcf.gz", "standardised.gwas.vcf.gz.tbi", "standardised.metadata.json",
                                   "standardised.qc.json"])
        self.assertEqual(list(scratch.iterdir()), [])

    def test_genomic_inflation_and_provenance_are_recorded(self):
        raw = write_table(self.tmp / "gwas.tsv", HEADER, [[1, 1000 + i, f"rs{i}", "A", "G", 0.3, 0.1, 0.02, p, 5000]
                                                          for i, p in enumerate([0.01, 0.05, 0.2])])
        result = run_standardise(self.tmp, raw, self.reference([]), input_build="hg19")
        self.assertAlmostEqual(result.qc["lambda_gc"], 8.44394650653, places=5)
        self.assertEqual(result.qc["provenance"]["schema_version"], result.metadata["pipeline_version"])
        self.assertEqual(len(result.qc["settings_hash"]), 16)


if __name__ == "__main__":
    unittest.main()
