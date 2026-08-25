"""Persistent run reporting for pipelign."""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
import tempfile
import time
import zipfile
from pathlib import Path


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def tool_metadata(tool_groups) -> list[dict[str, str | None]]:
    tools = []
    for names in tool_groups:
        resolved = next((shutil.which(name) for name in names if shutil.which(name)), None)
        version = None
        if resolved:
            for flag in ("--version", "-version", "-h"):
                try:
                    result = subprocess.run(
                        [resolved, flag],
                        capture_output=True,
                        text=True,
                        timeout=5,
                        check=False,
                    )
                except (OSError, subprocess.TimeoutExpired):
                    continue
                output = (result.stdout or result.stderr).strip().splitlines()
                if output:
                    version = output[0][:500]
                    break
        tools.append(
            {
                "requested": " or ".join(names),
                "path": resolved,
                "version": version,
            }
        )
    return tools


class RunReport:
    """Atomically persist run metadata and completed checkpoints."""

    def __init__(
        self,
        output_dir: Path,
        input_path: Path,
        output_path: Path,
        arguments: dict,
        tools,
        *,
        archive: bool = False,
    ) -> None:
        output_dir = Path(output_dir)
        input_path = Path(input_path)
        output_path = Path(output_path)
        self.output_dir = output_dir.resolve()
        self.path = self.output_dir / "pipelign-run.json"
        self.archive = archive
        self.data = {
            "schema_version": 1,
            "status": "running",
            "started_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "finished_at": None,
            "duration_seconds": None,
            "input": {
                "path": str(input_path.resolve()),
                "sha256": file_sha256(input_path),
            },
            "output": str(output_path.resolve()),
            "output_directory": str(self.output_dir),
            "arguments": arguments,
            "tools": tool_metadata(tools),
            "checkpoints": {},
            "summary": {},
            "error": None,
            "archive": None,
        }
        self._started = time.monotonic()
        self._write()

    @classmethod
    def resume(
        cls,
        output_dir: Path,
        input_path: Path,
        output_path: Path,
        arguments: dict,
        *,
        archive: bool = False,
    ) -> "RunReport":
        output_dir = Path(output_dir)
        input_path = Path(input_path)
        output_path = Path(output_path)
        path = output_dir.resolve() / "pipelign-run.json"
        if not path.is_file():
            raise ValueError("cannot resume: missing <%s>" % path)
        with path.open(encoding="utf-8") as handle:
            data = json.load(handle)
        if data.get("schema_version") != 1:
            raise ValueError("cannot resume an unsupported run-report schema")
        if data.get("input", {}).get("sha256") != file_sha256(input_path):
            raise ValueError("cannot resume: input file content has changed")
        previous = data.get("arguments", {})
        for key in (
            "alphabet",
            "lenThr",
            "simPer",
            "ambigPer",
            "merge",
            "mIterateLong",
            "mIterateMerge",
            "stage",
            "keepOrphans",
            "keepBadSeqs",
            "excludeClusters",
        ):
            if previous.get(key) != arguments.get(key):
                raise ValueError("cannot resume: argument %s has changed" % key)
        if data.get("output") != str(output_path.resolve()):
            raise ValueError("cannot resume: output path has changed")

        report = cls.__new__(cls)
        report.output_dir = output_dir.resolve()
        report.path = path
        report.archive = archive
        report.data = data
        report.data["status"] = "running"
        report.data["error"] = None
        report.data["finished_at"] = None
        report.data["duration_seconds"] = None
        report.data["output"] = str(output_path.resolve())
        report.data["arguments"] = arguments
        report.data["archive"] = None
        report._started = time.monotonic()
        report._write()
        return report

    def has_checkpoint(self, name: str) -> bool:
        return name in self.data["checkpoints"]

    def checkpoint_details(self, name: str) -> dict:
        if not self.has_checkpoint(name):
            raise KeyError(name)
        return dict(self.data["checkpoints"][name])

    def checkpoint(self, name: str, **details) -> None:
        artifacts = details.pop("artifacts", ())
        artifact_hashes = {}
        for artifact in artifacts:
            artifact_path = (self.output_dir / artifact).resolve()
            if artifact_path.is_file():
                artifact_hashes[str(artifact_path.relative_to(self.output_dir))] = (
                    file_sha256(artifact_path)
                )
        self.data["checkpoints"][name] = {
            "completed_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "artifacts": artifact_hashes,
            **details,
        }
        self._write()

    def update_summary(self, **details) -> None:
        self.data["summary"].update(details)
        self._write()

    def artifact_checksums(self) -> dict[str, str]:
        checksums = {}
        for checkpoint in self.data["checkpoints"].values():
            checksums.update(checkpoint.get("artifacts", {}))
        return checksums

    def finish(self, status: str, error: str | None = None) -> None:
        if self.data["status"] != "running":
            return
        self.data["status"] = status
        self.data["error"] = error
        self.data["finished_at"] = time.strftime(
            "%Y-%m-%dT%H:%M:%SZ", time.gmtime()
        )
        self.data["duration_seconds"] = round(time.monotonic() - self._started, 3)
        if self.archive and status == "success":
            self.data["archive"] = str(self.output_dir) + ".zip"
            self.data["status"] = "finalizing"
            self._write()
            try:
                self._publish_archive()
            except Exception as exc:
                self.data["status"] = "failed"
                self.data["error"] = "archive creation failed: %s" % exc
                self.data["archive"] = None
                self._write()
                raise
            self.data["status"] = "success"
            self._write()
        else:
            self.data["archive"] = None
            self._write()

    def _publish_archive(self) -> None:
        destination = Path(str(self.output_dir) + ".zip")
        descriptor, temporary_name = tempfile.mkstemp(
            prefix=".%s-" % self.output_dir.name,
            suffix=".zip",
            dir=destination.parent,
        )
        os.close(descriptor)
        temporary = Path(temporary_name)
        archived_data = dict(self.data)
        archived_data["status"] = "success"
        try:
            with zipfile.ZipFile(
                temporary, "w", compression=zipfile.ZIP_DEFLATED
            ) as archive:
                for path in self.output_dir.rglob("*"):
                    if not path.is_file():
                        continue
                    relative = path.relative_to(self.output_dir)
                    if path == self.path:
                        archive.writestr(
                            str(relative),
                            json.dumps(archived_data, indent=2, sort_keys=True) + "\n",
                        )
                    else:
                        archive.write(path, str(relative))
            with zipfile.ZipFile(temporary) as archive:
                corrupt = archive.testzip()
                if corrupt is not None:
                    raise RuntimeError("archive contains a corrupt entry: %s" % corrupt)
            with temporary.open("rb") as handle:
                os.fsync(handle.fileno())
            os.replace(temporary, destination)
            try:
                directory_fd = os.open(destination.parent, os.O_RDONLY)
            except OSError:
                pass
            else:
                try:
                    os.fsync(directory_fd)
                finally:
                    os.close(directory_fd)
        finally:
            if temporary.exists():
                temporary.unlink()

    def _write(self) -> None:
        temporary = self.path.with_name(self.path.name + ".tmp")
        with temporary.open("w", encoding="utf-8") as handle:
            json.dump(self.data, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, self.path)
