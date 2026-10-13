"""Exercise cleanup boundaries in disposable repositories, never real run data."""

import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest


SCRIPT = Path(__file__).resolve().parents[1] / "bin/cleanup_run.sh"


@unittest.skipUnless(sys.platform.startswith("linux"), "Cleanup uses Linux /proc")
class CleanupRunTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="silac-cleanup-test-")
        self.addCleanup(self.temp.cleanup)
        self.base = Path(self.temp.name)
        self.repo = self.base / "repo with spaces"
        (self.repo / "bin").mkdir(parents=True)
        self.script = self.repo / "bin/cleanup_run.sh"
        shutil.copyfile(SCRIPT, self.script)
        for name in ("main.nf", "nextflow.config", "README.md"):
            (self.repo / name).write_text("preserve\n")
        self.disposable = (
            "work/task/intermediate.mzML", ".nextflow/cache/state",
            ".nextflow.log", ".nextflow.log.1", "bin/__pycache__/helper.pyc",
            "tests/__pycache__/test.pyc",
        )
        self.preserved = (
            "results/run/experiment.mzTab", "results/run/nextflow.log",
            "assets/metadata/pilot.sdrf.tsv", "data/raw/input.raw",
            "misc/PROJECT_STATUS_AND_CONTEXT.md", "bin/helper.py",
        )
        for name in self.disposable + self.preserved:
            path = self.repo / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("fixture\n")

    def run_script(self, *args, env=None):
        return subprocess.run(["bash", str(self.script), *args], cwd=self.base,
                              text=True, capture_output=True, timeout=15, env=env)

    def assert_preserved(self):
        for name in self.preserved:
            self.assertEqual((self.repo / name).read_text(), "fixture\n")

    def test_default_is_preview(self):
        result = self.run_script()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("Preview only", result.stdout)
        self.assertTrue(all((self.repo / name).exists() for name in self.disposable))
        self.assert_preserved()

    def test_delete_from_another_directory_preserves_scientific_files(self):
        result = self.run_script("--delete")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(all(not (self.repo / name).exists() for name in self.disposable))
        self.assert_preserved()
        repeated = self.run_script("--delete")
        self.assertEqual(repeated.returncode, 0, repeated.stderr)
        self.assertIn("Nothing to clean", repeated.stdout)

    def test_symlink_cleanup_does_not_follow_external_target(self):
        shutil.rmtree(self.repo / "work")
        external = self.base / "external data"
        external.mkdir()
        (external / "keep.raw").write_text("original\n")
        (self.repo / "work").symlink_to(external, target_is_directory=True)
        result = self.run_script("--delete")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertFalse((self.repo / "work").is_symlink())
        self.assertEqual((external / "keep.raw").read_text(), "original\n")
        self.assert_preserved()

    def test_live_work_process_blocks_all_deletion(self):
        child = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(60)"],
                                 cwd=self.repo / "work")
        try:
            result = self.run_script("--delete")
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("still using work/", result.stderr)
            self.assertTrue(all((self.repo / name).exists() for name in self.disposable))
        finally:
            child.terminate()
            child.wait(timeout=5)

    def test_nextflow_launcher_blocks_deletion(self):
        child = subprocess.Popen(
            [sys.executable, "-c", "import time; time.sleep(60)", "nextflow.cli.Launcher"],
            cwd=self.repo)
        try:
            result = self.run_script("--delete")
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Nextflow PID", result.stderr)
            self.assertTrue((self.repo / ".nextflow.log").exists())
        finally:
            child.terminate()
            child.wait(timeout=5)

    def test_missing_repository_marker_blocks_deletion(self):
        (self.repo / "main.nf").unlink()
        result = self.run_script("--delete")
        self.assertNotEqual(result.returncode, 0)
        self.assertTrue((self.repo / ".nextflow.log").exists())

    def test_bad_flags_and_sudo_without_delete_are_rejected(self):
        for args in (("--all",), ("--sudo",)):
            with self.subTest(args=args):
                result = self.run_script(*args)
                self.assertNotEqual(result.returncode, 0)
                self.assertTrue((self.repo / ".nextflow.log").exists())

    @unittest.skipIf(os.geteuid() == 0, "Root bypasses directory permissions")
    def test_unwritable_container_directory_blocks_partial_cleanup(self):
        blocked = self.repo / "work/task"
        blocked.chmod(0o555)
        try:
            result = self.run_script("--delete")
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("--delete --sudo", result.stderr)
            self.assertTrue((self.repo / ".nextflow.log").exists())
            self.assertTrue((self.repo / "work/task/intermediate.mzML").exists())
        finally:
            blocked.chmod(0o755)

    def test_failed_sudo_authorization_does_not_delete_anything(self):
        mock_bin = self.base / "mock bin"
        mock_bin.mkdir()
        sudo = mock_bin / "sudo"
        sudo.write_text("#!/bin/sh\nexit 1\n")
        sudo.chmod(0o755)
        env = dict(os.environ, PATH=str(mock_bin) + os.pathsep + os.environ["PATH"])
        result = self.run_script("--delete", "--sudo", env=env)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("authorization failed", result.stderr)
        self.assertTrue(all((self.repo / name).exists() for name in self.disposable))


if __name__ == "__main__":
    unittest.main()
