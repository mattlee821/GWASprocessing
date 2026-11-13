"""Running the pipeline: src/standardise.sh and the Nextflow workflow it launches.

1 test; each README rule is a subTest. The wrapper's options and settings are checked inside a temporary copy of the
repository with a stub nextflow that records its arguments; rows running, failing, skipping and rerunning are checked
with the real Nextflow workflow on a 2-row synthetic manifest (skipped when nextflow, java, bgzip or tabix is missing).
No network: downloaders and installers are replaced by stubs that fail.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from testlib import ROOT, SRC, example_reference, generic_grch38_table, read_vcf_text, write_stub, write_table

import gwas_pipeline_common as common

BASH = shutil.which("bash")
BLOCKED_TOOLS = ("wget", "curl", "brew", "apt-get", "yum", "sudo", "conda", "mamba", "aws", "gsutil")
STUB_NEXTFLOW = r'''
printf '%s\n' "$@" > "${STUB_ARGS_FILE}"
log_root=""
while [[ $# -gt 0 ]]; do
  if [[ "$1" == "--log_root" ]]; then log_root="$2"; fi
  shift
done
if [[ -n "${STUB_SUMMARY_DIR:-}" ]]; then cp "${STUB_SUMMARY_DIR}"/* "${log_root}/"; fi
exit "${STUB_EXIT:-0}"
'''


def sandbox_env(tmp: Path, *bin_dirs: Path, **extra: str) -> dict:
    """PATH with the given stub folders first, then stubs that block network and installers, bash, and python3 only
    from the test's interpreter; GWAS_* and related settings from the caller's environment are removed."""
    guard = tmp / "guard_bin"
    for tool in BLOCKED_TOOLS:
        if not (guard / tool).exists():
            write_stub(guard, tool, f'echo "BLOCKED {tool} $*" >&2\nexit 97\n')
    bash_dir = tmp / "bash_bin"
    bash_dir.mkdir(exist_ok=True)
    if not (bash_dir / "bash").exists():
        (bash_dir / "bash").symlink_to(BASH)
    python_dir = tmp / "python_bin"
    python_dir.mkdir(exist_ok=True)
    if not (python_dir / "python3").exists():
        write_stub(python_dir, "python3", f'exec "{sys.executable}" "$@"\n')
    path = [*map(str, bin_dirs), str(guard), str(bash_dir), str(python_dir), "/usr/bin", "/bin", "/usr/sbin", "/sbin"]
    env = {k: v for k, v in os.environ.items()
           if not k.startswith(("GWAS_", "SLURM_", "NEXTFLOW", "GWASPROCESSING_", "PLINK2", "BCFTOOLS", "MAMBA"))}
    env.update({"PATH": os.pathsep.join(path)})
    env.update(extra)
    return env


def fake_repo(tmp: Path, dotenv: str = "") -> Path:
    """A copy of the parts of the repository the wrapper uses, with a placeholder main.nf."""
    repo = tmp / "repo"
    (repo / "src").mkdir(parents=True)
    (repo / "pipeline").mkdir()
    (repo / "pipeline" / "main.nf").write_text("// placeholder\n", encoding="utf-8")
    shutil.copytree(ROOT / "pipeline" / "bin", repo / "pipeline" / "bin", ignore=shutil.ignore_patterns("__pycache__"))
    shutil.copy(SRC / "standardise.sh", repo / "src" / "standardise.sh")
    shutil.copytree(SRC / "lib", repo / "src" / "lib")
    shutil.copy(ROOT / ".env.example", repo / ".env.example")
    if dotenv:
        (repo / ".env").write_text(dotenv, encoding="utf-8")
    return repo


class Wrapper:
    """src/standardise.sh in a fresh fake repository with a stub nextflow on PATH."""

    def __init__(self, tmp: Path, dotenv: str = ""):
        self.tmp = tmp
        tmp.mkdir(parents=True)
        self.repo = fake_repo(tmp, dotenv)
        self.bin = tmp / "bin"
        write_stub(self.bin, "nextflow", STUB_NEXTFLOW)
        self.args_file = tmp / "nextflow_args.txt"

    def env(self, **extra: str) -> dict:
        return sandbox_env(self.tmp, self.bin, STUB_ARGS_FILE=str(self.args_file), **extra)

    def run(self, *args: str, cwd: Path | None = None, env: dict | None = None, script: Path | None = None):
        script = script or self.repo / "src" / "standardise.sh"
        return subprocess.run([BASH, str(script), *args], cwd=cwd or self.tmp, env=env or self.env(), capture_output=True,
                              text=True, timeout=120)

    def value(self, flag: str) -> str:
        lines = self.args_file.read_text(encoding="utf-8").splitlines()
        return lines[lines.index(flag) + 1]

    def has(self, flag: str) -> bool:
        return flag in self.args_file.read_text(encoding="utf-8").splitlines()


class RunTest(unittest.TestCase):
    def test_run(self):
        base = Path(tempfile.mkdtemp(prefix="gwas_test_run_"))
        counter = iter(range(1000))

        def wrapper(dotenv: str = "") -> Wrapper:
            return Wrapper(base / f"wrapper{next(counter)}", dotenv)

        # ---- The wrapper and its settings --------------------------------------------------------------------------
        with self.subTest("every shell script parses"):
            scripts = sorted(SRC.rglob("*.sh"))
            self.assertGreater(len(scripts), 3)
            for script in scripts:
                proc = subprocess.run([BASH, "-n", str(script)], capture_output=True, text=True)
                self.assertEqual(proc.returncode, 0, f"{script.name}: {proc.stderr}")

        with self.subTest("options reach Nextflow with lowercase booleans and the local defaults"):
            w = wrapper()
            proc = w.run("--report", "False")
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertEqual((w.value("--report"), w.value("--force"), w.value("-profile"), w.value("--row_limit")),
                             ("false", "false", "local", "1"))
            self.assertEqual(w.value("run"), str(w.repo / "pipeline" / "main.nf"))
            self.assertEqual(w.value("--manifest"), str(w.repo / "manifests" / "test.tsv"))
            self.assertEqual(w.value("--work_root"), str(w.repo / "GWAS"))
            self.assertFalse(w.has("--only_study") or w.has("--study_yaml"))
            proc = w.run("--manifest", "manifests/other.tsv", "--row-limit", "0", "--only-study", "locke_2015_25673413",
                         "--study-yaml", "manifests/x.yaml", "--force", "--report", "TRUE")
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertEqual((w.value("--manifest"), w.value("--row_limit"), w.value("--only_study"), w.value("--study_yaml"),
                              w.value("--force"), w.value("--report")),
                             (str(w.repo / "manifests" / "other.tsv"), "0", "locke_2015_25673413",
                              str(w.repo / "manifests" / "x.yaml"), "true", "true"))

        with self.subTest("an invalid boolean or row limit is rejected"):
            w = wrapper()
            proc = w.run("--report", "maybe")
            self.assertEqual(proc.returncode, 2)
            self.assertIn("--report", proc.stderr)
            proc = w.run("--row-limit", "all")
            self.assertEqual(proc.returncode, 2)
            self.assertIn("--row-limit", proc.stderr)
            self.assertFalse(w.args_file.exists(), "Nextflow must not start")

        with self.subTest("an unknown option is rejected"):
            w = wrapper()
            proc = w.run("--bogus")
            self.assertEqual(proc.returncode, 2)
            self.assertIn("Unknown option: --bogus", proc.stderr)

        with self.subTest("values in the environment override the defaults in .env"):
            w = wrapper(dotenv="GWAS_WORK_ROOT=GWAS\nGWAS_REPORT=false\n# comment\nGWAS_ROW_LIMIT=3\n")
            work = base / "elsewhere" / "GWAS"
            proc = w.run(env=w.env(GWAS_WORK_ROOT=str(work)))
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertEqual((w.value("--work_root"), w.value("--report"), w.value("--row_limit")), (str(work), "false", "3"))
            proc = w.run("--row-limit", "5", env=w.env(GWAS_ROW_LIMIT="4"))
            self.assertEqual(w.value("--row_limit"), "5", "a command-line option wins over the environment")

        with self.subTest("the values in .env.example are never read"):
            w = wrapper()  # .env.example only: its /WHERE/... placeholders must not reach the run
            proc = w.run()
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertEqual(w.value("--work_root"), str(w.repo / "GWAS"))
            self.assertNotIn("/WHERE/", w.args_file.read_text(encoding="utf-8"))
            self.assertIn("src/source.sh --env-only", proc.stderr)

        with self.subTest("a path still holding a /WHERE/ placeholder in .env stops the run"):
            w = wrapper(dotenv="GWAS_WORK_ROOT=/WHERE/THE/GWAS/FILES/GO\n")
            proc = w.run()
            self.assertEqual(proc.returncode, 2)
            self.assertIn("placeholder", proc.stderr)
            self.assertFalse(w.args_file.exists())

        with self.subTest("sbatch's spooled copy finds the repository and uses the slurm profile for every row"):
            w = wrapper(dotenv="GWAS_PROFILE=local\n")
            spool = w.tmp / "spool" / "job1"
            spool.mkdir(parents=True)
            shutil.copy(w.repo / "src" / "standardise.sh", spool / "slurm_script")
            proc = w.run(script=spool / "slurm_script", cwd=w.repo, env=w.env(SLURM_JOB_ID="1", SLURM_SUBMIT_DIR=str(w.repo)))
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertEqual((w.value("run"), w.value("-profile"), w.value("--row_limit")),
                             (str(w.repo / "pipeline" / "main.nf"), "slurm", "0"))

        with self.subTest("tools come from the environment named in .env"):
            w = wrapper(dotenv="GWAS_WORK_ROOT=GWAS\nGWAS_ENV=.tools/env\n")
            write_stub(w.repo / ".tools" / "env" / "bin", "nextflow", STUB_NEXTFLOW)
            (w.bin / "nextflow").unlink()
            proc = w.run()
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertEqual(w.args_file.read_text(encoding="utf-8").splitlines()[0], "-log")

        with self.subTest("nextflow is found in the tools folder when it is not on PATH"):
            w = wrapper()
            write_stub(w.repo / ".tools" / "bin", "nextflow", STUB_NEXTFLOW)
            (w.bin / "nextflow").unlink()
            proc = w.run()
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertTrue(w.args_file.exists())

        with self.subTest("without nextflow the run stops and says how to install it"):
            w = wrapper()
            (w.bin / "nextflow").unlink()
            proc = w.run()
            self.assertEqual(proc.returncode, 2)
            self.assertIn("src/source.sh --tools-only", proc.stderr)

        with self.subTest("a container adds its profile and keeps the local row limit"):
            for runtime, profile in (("apptainer", "local,apptainer"), ("docker", "local,docker")):
                w = wrapper(dotenv=f"GWAS_WORK_ROOT=GWAS\nGWAS_CONTAINER={runtime}\n")
                proc = w.run()
                self.assertEqual(proc.returncode, 0, proc.stderr)
                self.assertEqual((w.value("-profile"), w.value("--row_limit")), (profile, "1"))
            w = wrapper(dotenv="GWAS_CONTAINER=podman\n")
            self.assertEqual(w.run().returncode, 2)

        with self.subTest("data folders in a cloud-synced folder trigger a warning"):
            w = wrapper()
            synced = base / "Library" / "CloudStorage" / "OneDrive-Org" / "GWAS"
            proc = w.run(env=w.env(GWAS_WORK_ROOT=str(synced)))
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertIn("cloud-synced", proc.stderr)
            w = wrapper()
            proc = w.run()
            self.assertNotIn("cloud-synced", proc.stderr)

        with self.subTest("failed rows are listed at the end and make the exit status non-zero"):
            w = wrapper()
            summaries = w.tmp / "summaries"
            summaries.mkdir()
            (summaries / "run_manifest_summary.json").write_text(json.dumps({"rows_to_run": 2, "rows_eligible": 2}), encoding="utf-8")
            (summaries / "run_manifest.tsv").write_text("study_id\toutput_id\trow_hash\ns1\tok\th1\ns2\tbroken\th2\n", encoding="utf-8")
            (summaries / "state_update_summary.json").write_text(json.dumps({"updated": 1, "updated_row_hashes": ["h1"]}),
                                                                  encoding="utf-8")
            proc = w.run(env=w.env(STUB_SUMMARY_DIR=str(summaries)))
            self.assertEqual(proc.returncode, 1)
            self.assertIn("1 row(s) failed", proc.stdout)
            self.assertIn("s2/broken", proc.stdout)
            self.assertNotIn("s1/ok", proc.stdout)
            self.assertTrue((w.repo / "GWAS" / "manifest.html").exists(), "manifest.html is updated at the end of every run")

        with self.subTest("summaries from an earlier run are not reported"):
            w = wrapper()
            logs = w.repo / "temp" / "logs"
            logs.mkdir(parents=True)
            (logs / "run_manifest_summary.json").write_text(json.dumps({"rows_to_run": 5}), encoding="utf-8")
            proc = w.run(env=w.env(STUB_EXIT="1"))
            self.assertEqual(proc.returncode, 1)
            self.assertNotIn("Rows to run: 5", proc.stdout)
            self.assertIn("stopped before the manifest was prepared", proc.stdout)

        # ---- The real workflow -------------------------------------------------------------------------------------
        tools = ("nextflow", "java", "bgzip", "tabix")
        missing = [t for t in tools if not shutil.which(t)]
        tmp = base / "nextflow"
        manifests = tmp / "manifests"
        manifests.mkdir(parents=True)
        good = generic_grch38_table(tmp / "good.tsv")  # GRCh38 positions: no liftover needed
        bad = write_table(tmp / "bad.tsv", ["CHR", "POS", "SNP", "EA", "OA", "Estimate", "SE", "P"],
                          [[1, 1000, "rs1", "A", "G", 0.1, 0.02, 1e-6]])  # no effect column: this row fails
        manifest = manifests / "manifest.tsv"
        manifest.write_text("GWAS_location\tphenotype\tancestry\tauthor\tyear\tPMID\tinput_build\n"
                            f'{good}\tBMI "adjusted" $HOME\tEUR\tLocke\t2015\t1\thg38\n'
                            f"{bad}\tWHR\tEUR\tBad\t2020\t2\thg38\n", encoding="utf-8")
        env = {**os.environ, "PATH": f"{Path(sys.executable).parent}{os.pathsep}{os.environ['PATH']}",
               "GWAS_WORK_ROOT": str(tmp / "GWAS"), "GWAS_REFERENCE_ROOT": str(example_reference()),
               "GWAS_CONFIG_ROOT": str(manifests), "GWAS_TEMP_ROOT": str(tmp / "temp"),
               "GWAS_LOG_ROOT": str(tmp / "temp" / "logs"), "GWAS_NEXTFLOW_WORK": str(tmp / "work"),
               "GWAS_NEXTFLOW_LOG": str(tmp / "temp" / ".nextflow.log"), "GWAS_REPORT": "true", "NXF_ANSI_LOG": "false"}

        def pipeline(*extra: str):
            cmd = [BASH, str(SRC / "standardise.sh"), "--manifest", str(manifest), "--profile", "local", "--row-limit", "0", *extra]
            return subprocess.run(cmd, env=env, capture_output=True, text=True, timeout=900)

        def state() -> dict:
            return common.read_state(tmp / "GWAS" / ".standardise_state.tsv")

        def summary() -> dict:
            return json.loads((tmp / "temp" / "logs" / "run_manifest_summary.json").read_text(encoding="utf-8"))

        first_hash = ""
        with self.subTest("rows run and fail independently; failed rows are listed and the exit status is non-zero"):
            if missing:
                self.skipTest(f"{', '.join(missing)} not found")
            first = pipeline()
            self.assertEqual(first.returncode, 1, first.stdout + first.stderr)
            self.assertIn("1 row(s) failed", first.stdout)
            self.assertIn("bad_2020_2/", first.stdout)
            (row,) = state().values()
            self.assertEqual(row["study_id"], "locke_2015_1")
            for key in ("standard_vcf", "standard_index", "metadata_json", "qc_json", "report_html"):
                self.assertGreater(Path(row[key]).stat().st_size, 0, key)
            self.assertIn('Phenotype="BMI \\"adjusted\\" $HOME"', read_vcf_text(Path(row["standard_vcf"])),
                          "manifest text reaches the VCF unchanged, never through a shell")
            page = (tmp / "GWAS" / "manifest.html").read_text(encoding="utf-8")
            self.assertIn("1 complete, 1 failed", page)
            first_hash = row["settings_hash"]

        with self.subTest("completed rows are skipped on the next run"):
            if missing:
                self.skipTest(f"{', '.join(missing)} not found")
            second = pipeline()
            self.assertEqual(second.returncode, 1)
            self.assertEqual((summary()["rows_skipped_complete"], summary()["rows_to_run"]), (1, 1))
            only = pipeline("--only-study", "locke_2015_1")
            self.assertEqual(only.returncode, 0, only.stdout + only.stderr)
            self.assertIn("Nothing launched: all selected rows already have complete outputs", only.stdout)

        with self.subTest("a change of settings reruns a completed row"):
            if missing:
                self.skipTest(f"{', '.join(missing)} not found")
            study_yaml = tmp / "tolerance.yaml"
            study_yaml.write_text("palindromic_tolerance: 0.1\n", encoding="utf-8")
            third = pipeline("--study-yaml", str(study_yaml), "--only-study", "locke_2015_1")
            self.assertEqual(third.returncode, 0, third.stdout + third.stderr)
            self.assertEqual(summary()["rows_to_run"], 1)
            (row,) = state().values()
            self.assertNotEqual(row["settings_hash"], first_hash)


if __name__ == "__main__":
    unittest.main()
