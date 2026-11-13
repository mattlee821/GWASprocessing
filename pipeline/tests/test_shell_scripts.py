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
    # Only python3 from the test's interpreter folder: in a tool environment that folder also holds nextflow, bcftools
    # and the rest, which would hide the stubs.
    python_dir = tmp / "python_bin"
    python_dir.mkdir(exist_ok=True)
    if not (python_dir / "python3").exists():  # a wrapper, not a link: a link would lose a virtualenv's packages
        write_stub(python_dir, "python3", f'exec "{sys.executable}" "$@"\n')
    path = [*map(str, bin_dirs), str(guard), str(bash_dir), str(python_dir), "/usr/bin", "/bin", "/usr/sbin", "/sbin"]
    env = {k: v for k, v in os.environ.items()
           if not k.startswith(("GWAS_", "SLURM_", "NEXTFLOW", "GWASPROCESSING_", "PLINK2", "BCFTOOLS", "MAMBA"))}
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
    shutil.copy(ROOT / "environment.yml", repo / "environment.yml")
    shutil.copy(ROOT / ".env.example", repo / ".env.example")
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
        proc = self.run_wrapper(fake_repo(self.tmp), "--report", "maybe")
        self.assertEqual(proc.returncode, 2)
        self.assertIn("--report", proc.stderr)

    def test_options_reach_nextflow_with_lowercase_booleans_and_local_defaults(self):
        repo = fake_repo(self.tmp)
        proc = self.run_wrapper(repo, "--report", "False")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual((self.value("--report"), self.value("-profile"), self.value("--row_limit")), ("false", "local", "1"))
        self.assertEqual(self.value("run"), str(repo / "pipeline" / "main.nf"))

    def test_environment_values_override_dotenv_defaults(self):
        repo = fake_repo(self.tmp, dotenv="GWAS_WORK_ROOT=GWAS\nGWAS_REPORT=true\n# comment\n")
        work = self.tmp / "elsewhere" / "GWAS"
        proc = self.run_wrapper(repo, env=self.env(GWAS_WORK_ROOT=str(work)))
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual((self.value("--work_root"), self.value("--report")), (str(work), "true"))

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

    def test_template_values_are_never_read(self):
        repo = fake_repo(self.tmp)  # .env.example only: its /WHERE/... placeholders must not reach the run
        proc = self.run_wrapper(repo)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(self.value("--work_root"), str(repo / "GWAS"))
        self.assertIn("src/source.sh --env-only", proc.stderr)

    def test_old_server_env_is_read_with_a_warning_when_there_is_no_dotenv(self):
        repo = fake_repo(self.tmp)
        work = self.tmp / "old" / "GWAS"
        (repo / "server.env").write_text(f"GWAS_WORK_ROOT={work}\n", encoding="utf-8")
        proc = self.run_wrapper(repo)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(self.value("--work_root"), str(work))
        self.assertIn("server.env, the old name", proc.stderr)

    def test_dotenv_wins_over_an_old_server_env(self):
        repo = fake_repo(self.tmp, dotenv=f"GWAS_WORK_ROOT={self.tmp / 'new'}\n")
        (repo / "server.env").write_text(f"GWAS_WORK_ROOT={self.tmp / 'old'}\n", encoding="utf-8")
        proc = self.run_wrapper(repo)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(self.value("--work_root"), str(self.tmp / "new"))
        self.assertNotIn("server.env", proc.stderr)

    def test_nextflow_is_found_in_the_tools_folder_when_not_on_path(self):
        repo = fake_repo(self.tmp)
        write_stub(repo / ".tools" / "bin", "nextflow", STUB_NEXTFLOW)
        (self.bin / "nextflow").unlink()
        proc = self.run_wrapper(repo)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertTrue(self.args_file.exists())

    def test_tools_come_from_the_environment_named_in_dotenv(self):
        repo = fake_repo(self.tmp, dotenv="GWAS_WORK_ROOT=GWAS\nGWAS_ENV=.tools/env\n")
        write_stub(repo / ".tools" / "env" / "bin", "nextflow", STUB_NEXTFLOW)
        (self.bin / "nextflow").unlink()
        proc = self.run_wrapper(repo)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertTrue(self.args_file.exists())

    def test_a_container_adds_its_profile_and_keeps_the_local_row_limit(self):
        repo = fake_repo(self.tmp, dotenv="GWAS_WORK_ROOT=GWAS\nGWAS_CONTAINER=apptainer\n")
        proc = self.run_wrapper(repo)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual((self.value("-profile"), self.value("--row_limit")), ("local,apptainer", "1"))


class SetupScriptTests(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.repo = fake_repo(self.tmp, dotenv="EXISTING=1\nGWAS_WORK_ROOT=GWAS\n")

    def test_set_env_var_adds_and_replaces_keys_without_touching_other_lines(self):
        template = (self.repo / ".env.example").read_text(encoding="utf-8")
        script = (f'source "{self.repo / "src" / "source.sh"}"; set_env_var NEW_KEY first; '
                  'set_env_var NEW_KEY "second=value"; set_env_var GWAS_WORK_ROOT /data/GWAS')
        proc = subprocess.run([BASH, "-c", script], capture_output=True, text=True, cwd=self.tmp, env=sandbox_env(self.tmp), timeout=60)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual((self.repo / ".env").read_text(encoding="utf-8").splitlines(),
                         ["EXISTING=1", "GWAS_WORK_ROOT=/data/GWAS", "NEW_KEY=second=value"])
        self.assertEqual((self.repo / ".env.example").read_text(encoding="utf-8"), template)

    def env_only(self, **extra):
        return subprocess.run([BASH, str(self.repo / "src" / "source.sh"), "--env-only"], capture_output=True, text=True,
                              env=sandbox_env(self.tmp, **extra), timeout=60)

    def settings(self):
        lines = (self.repo / ".env").read_text(encoding="utf-8").splitlines()
        return dict(line.split("=", 1) for line in lines if "=" in line and not line.lstrip().startswith("#"))

    def test_setup_creates_dotenv_from_the_template_and_fills_the_paths(self):
        (self.repo / ".env").unlink()
        home = self.tmp / "data"
        proc = self.env_only(GWAS_DATA_HOME=str(home))
        self.assertEqual(proc.returncode, 0, proc.stderr)
        settings = self.settings()
        self.assertEqual((settings["GWAS_WORK_ROOT"], settings["GWAS_REFERENCE_ROOT"], settings["GWAS_NEXTFLOW_LOG"]),
                         (str(home / "GWAS"), str(home / "references"), str(home / "temp" / ".nextflow.log")))
        self.assertFalse(any(v.startswith("/WHERE/") for v in settings.values()))
        self.assertEqual(settings["OPENGWAS_JWT"], "")
        self.assertIn("OPENGWAS_JWT", proc.stdout)
        self.assertTrue((self.repo / ".env").read_text(encoding="utf-8").startswith("# GWASprocessing settings. .env.example is the tracked template"))

    def test_setup_keeps_values_already_in_dotenv(self):
        (self.repo / ".env").write_text("GWAS_WORK_ROOT=/mine/GWAS\nOPENGWAS_JWT=abc\n", encoding="utf-8")
        proc = self.env_only(GWAS_DATA_HOME=str(self.tmp / "data"))
        self.assertEqual(proc.returncode, 0, proc.stderr)
        settings = self.settings()
        self.assertEqual((settings["GWAS_WORK_ROOT"], settings["OPENGWAS_JWT"]), ("/mine/GWAS", "abc"))
        self.assertEqual(settings["GWAS_TEMP_ROOT"], str(self.tmp / "data" / "temp"))  # missing paths are added
        self.assertNotIn("OPENGWAS_JWT", proc.stdout)
        again = self.env_only(GWAS_DATA_HOME=str(self.tmp / "elsewhere"))
        self.assertEqual(again.returncode, 0, again.stderr)
        self.assertEqual(self.settings(), settings)

    def test_setup_moves_an_old_server_env_into_dotenv(self):
        (self.repo / ".env").unlink()
        (self.repo / "server.env").write_text("GWAS_WORK_ROOT=/old/GWAS\nGWAS_BCFTOOLS=/old/bcftools\nOPENGWAS_JWT=tok\n", encoding="utf-8")
        proc = self.env_only(GWAS_DATA_HOME=str(self.tmp / "data"))
        self.assertEqual(proc.returncode, 0, proc.stderr)
        settings = self.settings()
        self.assertEqual((settings["GWAS_WORK_ROOT"], settings["GWAS_BCFTOOLS"], settings["OPENGWAS_JWT"]),
                         ("/old/GWAS", "/old/bcftools", "tok"))
        self.assertEqual(settings["GWAS_TEMP_ROOT"], str(self.tmp / "data" / "temp"))
        self.assertIn("Moved 3 settings from server.env", proc.stdout)
        self.assertNotIn("tok", proc.stdout + proc.stderr)

    def test_unknown_option_is_rejected(self):
        proc = subprocess.run([BASH, str(self.repo / "src" / "source.sh"), "--bogus"], capture_output=True, text=True,
                              env=sandbox_env(self.tmp), timeout=60)
        self.assertEqual(proc.returncode, 2)

    def install_stubs(self):
        """A stub micromamba that records its calls and creates an environment with stub tools (bcftools lists the
        liftover plugin, so nothing is built)."""
        log = self.tmp / "mamba.log"
        write_stub(self.repo / ".tools" / "bin", "micromamba", f'''echo "micromamba $*" >> "{log}"
prefix=""; while [[ $# -gt 0 ]]; do [[ "$1" == "-p" ]] && prefix="$2"; shift; done
mkdir -p "$prefix/conda-meta" "$prefix/bin" "$prefix/libexec/bcftools"
for tool in nextflow plink2 python3 bgzip tabix; do printf '#!/usr/bin/env bash\\nexit 0\\n' > "$prefix/bin/$tool"; chmod +x "$prefix/bin/$tool"; done
printf '#!/usr/bin/env bash\\n[[ "$1 $2" == "plugin -l" ]] && echo liftover\\nexit 0\\n' > "$prefix/bin/bcftools"; chmod +x "$prefix/bin/bcftools"
''')
        return log

    def tools_only(self, *extra_args, bin_dirs=()):
        return subprocess.run([BASH, str(self.repo / "src" / "source.sh"), "--tools-only", *extra_args], capture_output=True,
                              text=True, env=sandbox_env(self.tmp, *bin_dirs, GWAS_DATA_HOME=str(self.tmp / "data")), timeout=120)

    def test_tools_are_installed_into_1_environment_from_environment_yml(self):
        log = self.install_stubs()
        proc = self.tools_only()
        self.assertEqual(proc.returncode, 0, proc.stdout + proc.stderr)
        self.assertNotIn("BLOCKED", proc.stderr)
        env = self.repo / ".tools" / "env"
        calls = log.read_text(encoding="utf-8")
        self.assertIn(f"create -y -q -p {env} -f {self.repo / 'environment.yml'}", calls)
        settings = self.settings()
        self.assertEqual((settings["GWAS_ENV"], settings["NEXTFLOW"], settings["PLINK2"], settings["GWAS_BCFTOOLS"],
                          settings["BCFTOOLS_PLUGINS"]),
                         (str(env), str(env / "bin" / "nextflow"), str(env / "bin" / "plink2"), str(env / "bin" / "bcftools"),
                          str(env / "libexec" / "bcftools")))
        self.assertIn("bcftools +liftover: found", proc.stdout)

    def test_an_unchanged_environment_is_not_reinstalled_and_a_changed_one_is_updated(self):
        log = self.install_stubs()
        self.assertEqual(self.tools_only().returncode, 0)
        log.write_text("", encoding="utf-8")
        again = self.tools_only()
        self.assertEqual(again.returncode, 0, again.stdout + again.stderr)
        self.assertIn("tool environment: up to date", again.stdout)
        self.assertEqual(log.read_text(encoding="utf-8"), "")
        with open(self.repo / "environment.yml", "a", encoding="utf-8") as handle:
            handle.write("  - tabulate\n")
        changed = self.tools_only()
        self.assertEqual(changed.returncode, 0, changed.stdout + changed.stderr)
        self.assertIn("micromamba install -y -q -p", log.read_text(encoding="utf-8"))

    def test_a_tool_the_user_set_in_dotenv_is_kept(self):
        self.install_stubs()
        mine = write_stub(self.tmp / "mine", "nextflow", "exit 0\n")
        with open(self.repo / ".env", "a", encoding="utf-8") as handle:
            handle.write(f"NEXTFLOW={mine}\n")
        proc = self.tools_only()
        self.assertEqual(proc.returncode, 0, proc.stdout + proc.stderr)
        self.assertEqual(self.settings()["NEXTFLOW"], str(mine))
        self.assertIn(f"NEXTFLOW: keeping {mine}", proc.stdout)

    def test_the_container_route_pulls_the_image_and_builds_nothing_on_the_host(self):
        self.install_stubs()
        bin_dir = self.tmp / "bin"
        pulls = self.tmp / "apptainer.log"
        write_stub(bin_dir, "apptainer", f'echo "apptainer $*" >> "{pulls}"; [[ "$1" == pull ]] && touch "$3"; exit 0\n')
        proc = self.tools_only("--container", "apptainer", bin_dirs=(bin_dir,))
        self.assertEqual(proc.returncode, 0, proc.stdout + proc.stderr)
        sif = self.repo / ".tools" / "gwasprocessing.sif"
        self.assertIn(f"apptainer pull --force {sif} docker://ghcr.io/mattlee821/gwasprocessing:latest", pulls.read_text(encoding="utf-8"))
        settings = self.settings()
        self.assertEqual((settings["GWAS_CONTAINER"], settings["GWAS_CONTAINER_IMAGE"]), ("apptainer", str(sif)))
        self.assertNotIn("GWAS_BCFTOOLS", settings)
        self.assertNotIn("bcftools", proc.stdout)

    def test_an_unknown_container_runtime_is_rejected(self):
        proc = self.tools_only("--container", "podman")
        self.assertEqual(proc.returncode, 2)


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


if __name__ == "__main__":
    unittest.main()
