"""Tests for the bash entry points: src/standardise.sh, src/source.sh and src/reference_data/*.sh.

Scripts run for real inside a temporary copy of the repository layout, with stub executables
(nextflow, wget, plink2, ...) recording how they were called. The real repository is never modified.
"""

import gzip
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from testlib import ROOT, SRC, write_stub

BASH = shutil.which("bash")
# Commands that could reach the network or change the machine. Every shell test runs with these replaced
# by stubs that fail immediately, so no test can download data or install software whatever the script does.
BLOCKED_TOOLS = ("wget", "curl", "brew", "apt-get", "yum", "sudo", "conda", "mamba", "aws", "gsutil")


def sandbox_env(tmp: Path, *bin_dirs: Path, **extra: str) -> dict:
    guard = tmp / "guard_bin"
    for tool in BLOCKED_TOOLS:
        if not (guard / tool).exists():
            write_stub(guard, tool, f'echo "BLOCKED {tool} $*" >&2\nexit 97\n')
    bash_dir = tmp / "bash_bin"
    bash_dir.mkdir(exist_ok=True)
    if not (bash_dir / "bash").exists():
        (bash_dir / "bash").symlink_to(BASH)
    path = [*map(str, bin_dirs), str(guard), str(bash_dir), str(Path(sys.executable).parent), "/usr/bin", "/bin", "/usr/sbin", "/sbin"]
    env = {k: v for k, v in os.environ.items() if not k.startswith(("GWAS_", "SLURM_", "NEXTFLOW", "GWASPROCESSING_", "PLINK2"))}
    env.update({"PATH": os.pathsep.join(path), "PIP_NO_INDEX": "1", "PIP_DISABLE_PIP_VERSION_CHECK": "1"})
    env.update(extra)
    return env

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


def fake_repo(tmp: Path, dotenv: str = "") -> Path:
    repo = tmp / "repo"
    (repo / "src").mkdir(parents=True)
    (repo / "pipeline").mkdir()
    (repo / "pipeline" / "main.nf").write_text("// placeholder\n", encoding="utf-8")
    shutil.copytree(ROOT / "pipeline" / "bin", repo / "pipeline" / "bin", ignore=shutil.ignore_patterns("__pycache__"))
    shutil.copy(SRC / "standardise.sh", repo / "src" / "standardise.sh")
    shutil.copy(SRC / "source.sh", repo / "src" / "source.sh")
    shutil.copytree(SRC / "reference_data", repo / "src" / "reference_data")
    shutil.copytree(SRC / "lib", repo / "src" / "lib")
    shutil.copy(ROOT / "requirements.txt", repo / "requirements.txt")
    if dotenv:
        (repo / ".env").write_text(dotenv, encoding="utf-8")
    return repo


class ScriptSyntaxTests(unittest.TestCase):
    def test_every_shell_script_parses(self):
        scripts = sorted(SRC.rglob("*.sh"))
        self.assertGreater(len(scripts), 3)
        for script in scripts:
            with self.subTest(script=script.name):
                proc = subprocess.run([BASH, "-n", str(script)], capture_output=True, text=True)
                self.assertEqual(proc.returncode, 0, proc.stderr)


class WrapperTests(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.bin = self.tmp / "bin"
        write_stub(self.bin, "nextflow", STUB_NEXTFLOW)
        self.args_file = self.tmp / "nextflow_args.txt"

    def env(self, **extra):
        return sandbox_env(self.tmp, self.bin, STUB_ARGS_FILE=str(self.args_file), **extra)

    def run_wrapper(self, repo, *args, cwd=None, env=None, script=None):
        script = script or repo / "src" / "standardise.sh"
        return subprocess.run([BASH, str(script), *args], cwd=cwd or self.tmp, env=env or self.env(), capture_output=True, text=True,
                              timeout=120)

    def value(self, flag):
        lines = self.args_file.read_text(encoding="utf-8").splitlines()
        return lines[lines.index(flag) + 1]

    def test_invalid_boolean_option_is_rejected(self):
        proc = self.run_wrapper(fake_repo(self.tmp), "--qcplot", "maybe")
        self.assertEqual(proc.returncode, 2)
        self.assertIn("--qcplot", proc.stderr)

    def test_options_reach_nextflow_with_lowercase_booleans_and_local_defaults(self):
        repo = fake_repo(self.tmp)
        proc = self.run_wrapper(repo, "--qcplot", "False")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual((self.value("--qcplot"), self.value("-profile"), self.value("--row_limit")), ("false", "local", "1"))
        self.assertEqual(self.value("run"), str(repo / "pipeline" / "main.nf"))

    def test_environment_values_override_dotenv_defaults(self):
        repo = fake_repo(self.tmp, dotenv="GWAS_WORK_ROOT=GWAS\nGWAS_QCPLOT=true\n# comment\n")
        work = self.tmp / "elsewhere" / "GWAS"
        proc = self.run_wrapper(repo, env=self.env(GWAS_WORK_ROOT=str(work)))
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual((self.value("--work_root"), self.value("--qcplot")), (str(work), "true"))

    def test_sbatch_spooled_copy_finds_the_repo_and_uses_the_slurm_profile(self):
        repo = fake_repo(self.tmp, dotenv="GWAS_PROFILE=local\n")
        spool = self.tmp / "spool" / "job1"
        spool.mkdir(parents=True)
        shutil.copy(repo / "src" / "standardise.sh", spool / "slurm_script")
        proc = self.run_wrapper(repo, script=spool / "slurm_script", cwd=repo,
                                env=self.env(SLURM_JOB_ID="1", SLURM_SUBMIT_DIR=str(repo)))
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual((self.value("run"), self.value("-profile"), self.value("--row_limit")),
                         (str(repo / "pipeline" / "main.nf"), "slurm", "0"))

    def test_failed_rows_are_listed_and_make_the_exit_status_non_zero(self):
        repo = fake_repo(self.tmp)
        summaries = self.tmp / "summaries"
        summaries.mkdir()
        (summaries / "run_manifest_summary.json").write_text(json.dumps({"rows_to_run": 2, "rows_eligible": 2}), encoding="utf-8")
        (summaries / "run_manifest.tsv").write_text("study_id\toutput_id\trow_hash\ns1\tok\th1\ns2\tbroken\th2\n", encoding="utf-8")
        (summaries / "state_update_summary.json").write_text(json.dumps({"updated": 1, "updated_row_hashes": ["h1"]}), encoding="utf-8")
        proc = self.run_wrapper(repo, env=self.env(STUB_SUMMARY_DIR=str(summaries)))
        self.assertEqual(proc.returncode, 1)
        self.assertIn("s2/broken", proc.stdout)
        self.assertIn("1 row(s) failed", proc.stdout)

    def test_stale_summaries_from_an_earlier_run_are_not_reported(self):
        repo = fake_repo(self.tmp)
        logs = repo / "temp" / "logs"
        logs.mkdir(parents=True)
        (logs / "run_manifest_summary.json").write_text(json.dumps({"rows_to_run": 5}), encoding="utf-8")
        proc = self.run_wrapper(repo, env=self.env(STUB_EXIT="1"))
        self.assertEqual(proc.returncode, 1)
        self.assertNotIn("Rows to run: 5", proc.stdout)

    def test_cloud_synced_data_folders_trigger_a_warning(self):
        repo = fake_repo(self.tmp)
        synced = self.tmp / "Library" / "CloudStorage" / "OneDrive-Org" / "GWAS"
        proc = self.run_wrapper(repo, env=self.env(GWAS_WORK_ROOT=str(synced)))
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("cloud-synced", proc.stderr)

    def test_nextflow_is_found_in_the_tools_folder_when_not_on_path(self):
        repo = fake_repo(self.tmp)
        write_stub(repo / ".tools" / "bin", "nextflow", STUB_NEXTFLOW)
        (self.bin / "nextflow").unlink()
        proc = self.run_wrapper(repo)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertTrue(self.args_file.exists())


class SetupScriptTests(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.repo = fake_repo(self.tmp, dotenv="EXISTING=1\nGWAS_WORK_ROOT=GWAS\n")

    def test_set_env_var_adds_and_replaces_keys_without_touching_other_lines(self):
        (self.repo / "server.env").write_text("EXISTING=1\nGWAS_WORK_ROOT=GWAS\n", encoding="utf-8")
        dotenv = (self.repo / ".env").read_text(encoding="utf-8")
        script = (f'source "{self.repo / "src" / "source.sh"}"; set_env_var NEW_KEY first; '
                  'set_env_var NEW_KEY "second=value"; set_env_var GWAS_WORK_ROOT /data/GWAS')
        proc = subprocess.run([BASH, "-c", script], capture_output=True, text=True, cwd=self.tmp, env=sandbox_env(self.tmp), timeout=60)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual((self.repo / "server.env").read_text(encoding="utf-8").splitlines(),
                         ["EXISTING=1", "GWAS_WORK_ROOT=/data/GWAS", "NEW_KEY=second=value"])
        self.assertEqual((self.repo / ".env").read_text(encoding="utf-8"), dotenv)

    def test_unknown_option_is_rejected(self):
        proc = subprocess.run([BASH, str(self.repo / "src" / "source.sh"), "--bogus"], capture_output=True, text=True,
                              env=sandbox_env(self.tmp), timeout=60)
        self.assertEqual(proc.returncode, 2)

    def test_tools_only_run_installs_nothing_when_everything_is_present(self):
        bin_dir = self.tmp / "bin"
        log = self.tmp / "calls.log"
        for tool in ("plink2", "nextflow", "bgzip", "tabix"):
            write_stub(bin_dir, tool, f'echo "{tool} $*" >> "{log}"\n')
        write_stub(bin_dir, "bcftools", 'if [[ "$1 $2" == "plugin -l" ]]; then echo liftover; fi\n')
        proc = subprocess.run([BASH, str(self.repo / "src" / "source.sh"), "--tools-only"], capture_output=True, text=True,
                              env=sandbox_env(self.tmp, bin_dir), timeout=120)
        self.assertEqual(proc.returncode, 0, proc.stdout + proc.stderr)
        self.assertNotIn("BLOCKED", proc.stderr)
        self.assertIn("python packages: already available", proc.stdout)


class ReferenceDownloadScriptTests(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.repo = fake_repo(self.tmp, dotenv="GWAS_REFERENCE_ROOT=references\n")
        self.bin = self.tmp / "bin"
        self.log = self.tmp / "wget.log"
        fasta = gzip.compress(b">1\nACGTACGTAC\n>MT\nGATTACA\n")
        (self.tmp / "fake.fa.gz").write_bytes(fasta)
        write_stub(self.bin, "wget", f'''echo "$*" >> "{self.log}"
out=""; prev=""
for arg in "$@"; do if [[ "$prev" == "-O" ]]; then out="$arg"; fi; prev="$arg"; done
case "$out" in *.fa.gz*) cp "{self.tmp / 'fake.fa.gz'}" "$out" ;; *) echo payload > "$out" ;; esac
''')
        self.reference = self.tmp / "reference"

    def run_script(self, name):
        env = sandbox_env(self.tmp, self.bin, GWAS_REFERENCE_ROOT=str(self.reference))
        return subprocess.run([BASH, str(self.repo / "src" / "reference_data" / name)], capture_output=True, text=True, env=env, timeout=120)

    def test_existing_files_are_not_downloaded_again(self):
        raw = self.reference / "1000genomes" / "phase3" / "raw"
        raw.mkdir(parents=True)
        (raw / "all_phase3.psam").write_text("present\n", encoding="utf-8")
        proc = self.run_script("download_1000genomes_phase3.sh")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        calls = self.log.read_text(encoding="utf-8")
        self.assertEqual(calls.count("\n"), 5)
        self.assertNotIn("all_phase3.psam", calls)
        self.assertEqual((raw / "all_phase3.psam").read_text(encoding="utf-8"), "present\n")

    def test_chains_and_grch38_fasta_are_fetched_and_the_fasta_is_indexed(self):
        proc = self.run_script("download_fasta_and_chains.sh")
        self.assertEqual(proc.returncode, 0, proc.stdout + proc.stderr)
        genome = self.reference / "genome"
        self.assertTrue((genome / "hg19ToHg38.over.chain.gz").exists())
        self.assertTrue((genome / "hg38ToHg19.over.chain.gz").exists())
        self.assertEqual((genome / "GRCh38.fa").read_text(encoding="utf-8"), ">1\nACGTACGTAC\n>MT\nGATTACA\n")
        self.assertEqual((genome / "GRCh38.fa.fai").read_text(encoding="utf-8").splitlines()[0].split("\t")[:2], ["1", "10"])
        before = self.log.read_text(encoding="utf-8")
        self.assertEqual(self.run_script("download_fasta_and_chains.sh").returncode, 0)
        self.assertEqual(self.log.read_text(encoding="utf-8"), before, "a second run must not download anything")


class BuildScriptTests(unittest.TestCase):
    STUB_PLINK2 = r'''
echo "$*" >> "${STUB_PLINK_LOG}"
out=""; prev=""; make_pgen=0; make_bed=0; freq=0
for arg in "$@"; do
  [[ "$prev" == "--out" ]] && out="$arg"
  [[ "$arg" == "--make-pgen" ]] && make_pgen=1
  [[ "$arg" == "--make-bed" ]] && make_bed=1
  [[ "$arg" == "--freq" ]] && freq=1
  prev="$arg"
done
if [[ "$1" == "--zst-decompress" ]]; then echo "decompressed $2"; exit 0; fi
if [[ "$make_pgen" == 1 ]]; then echo pgen > "$out.pgen"; printf "#IID\tPAT\tMAT\tSEX\tSuperPop\tPopulation\nS1\t0\t0\t1\tEUR\tGBR\nS2\t0\t0\t2\tAFR\tYRI\n" > "$out.psam"; fi
if [[ "$make_bed" == 1 ]]; then echo bim > "$out.bim"; fi
if [[ "$freq" == 1 ]]; then echo afreq > "$out.afreq"; fi
'''

    def test_each_plink2_step_runs_once_and_population_builds_use_the_documented_filters(self):
        tmp = Path(tempfile.mkdtemp())
        repo = fake_repo(tmp)
        reference = tmp / "reference"
        raw = reference / "1000genomes" / "phase3" / "raw"
        raw.mkdir(parents=True)
        for name in ("all_phase3.pgen.zst", "all_phase3.pvar.zst", "all_phase3.psam", "deg2_phase3.king.cutoff.out.id"):
            (raw / name).write_text("x\n", encoding="utf-8")
        log = tmp / "plink.log"
        bin_dir = tmp / "bin"
        write_stub(bin_dir, "plink2", self.STUB_PLINK2)
        env = sandbox_env(tmp, bin_dir, GWAS_REFERENCE_ROOT=str(reference), STUB_PLINK_LOG=str(log))
        script = repo / "src" / "reference_data" / "build_1000genomes_phase3.sh"
        first = subprocess.run([BASH, str(script)], capture_output=True, text=True, env=env, timeout=120)
        self.assertEqual(first.returncode, 0, first.stdout + first.stderr)
        calls = log.read_text(encoding="utf-8").splitlines()
        self.assertEqual(len(calls), 15)
        bed_calls = [c for c in calls if "--make-bed" in c]
        self.assertEqual(len(bed_calls), 6)
        self.assertTrue(all("--chr 1-22,X" in c and "--maf 0.01" in c for c in bed_calls))
        self.assertEqual([("--keep" in c) for c in bed_calls], [True] * 5 + [False])
        second = subprocess.run([BASH, str(script)], capture_output=True, text=True, env=env, timeout=120)
        self.assertEqual(second.returncode, 0, second.stderr)
        self.assertEqual(len(log.read_text(encoding="utf-8").splitlines()), 15, "a second run must not call plink2 again")


if __name__ == "__main__":
    unittest.main()
