"""Per-study logs (copy_logs.py) and the scratch-folder free-space check (make_work_dir)."""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import testlib  # noqa: F401  (puts pipeline/bin on the path)
from testlib import BIN
import gwas_pipeline_common as common


class CopyLogsTest(unittest.TestCase):
    def run_script(self, cwd, *args):
        return subprocess.run([sys.executable, str(BIN / "copy_logs.py"), *args], cwd=cwd, check=True)

    def test_row_and_run_logs_go_to_the_study_logs_folder(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            task = tmp / "task"
            task.mkdir()
            (task / ".command.log").write_text("out\n")
            (task / ".command.err").write_text("err\n")
            (task / "row.json").write_text(json.dumps({"study_id": "s_1", "output_id": "P1_1-2"}))
            work = tmp / "GWAS"
            self.run_script(task, "--work-root", str(work), "--row-json", "row.json")
            logs = work / "s_1" / "logs"
            self.assertEqual(len(list(logs.glob("P1_1-2.*.log"))), 1)
            self.assertEqual(next(logs.glob("P1_1-2.*.err")).read_text(), "err\n")
            (tmp / "run_manifest.tsv").write_text("study_id\toutput_id\ns_1\tP1_1-2\n")
            (tmp / "nf.log").write_text("nextflow\n")
            (tmp / "slurm.out").write_text("slurm\n")
            self.run_script(tmp, "--work-root", str(work), "--run-manifest", "run_manifest.tsv",
                            "--nextflow-log", "nf.log", "--slurm-out", "slurm.out")
            self.assertEqual(len(list(logs.glob("run.*.nextflow.log"))), 1)
            self.assertEqual(len(list(logs.glob("run.*.slurm.out"))), 1)

    def test_a_missing_task_log_is_not_an_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            (tmp / "row.json").write_text(json.dumps({"study_id": "s_1", "output_id": "o"}))
            self.run_script(tmp, "--work-root", str(tmp / "GWAS"), "--row-json", "row.json")


class ScratchTest(unittest.TestCase):
    def test_free_space_check(self):
        with tempfile.TemporaryDirectory() as tmp:
            scratch = Path(tmp) / "new" / "scratch"
            with mock.patch.dict(os.environ, {"GWAS_TMPDIR": str(scratch), "GWAS_MIN_FREE_GB": "1e9"}):
                with self.assertRaises(SystemExit) as raised:
                    common.make_work_dir()
                self.assertEqual(raised.exception.code, 1)
            with mock.patch.dict(os.environ, {"GWAS_TMPDIR": str(scratch), "GWAS_MIN_FREE_GB": "0"}):
                made = common.make_work_dir()
                self.assertTrue(made.is_dir() and made.parent == scratch)


if __name__ == "__main__":
    unittest.main()
