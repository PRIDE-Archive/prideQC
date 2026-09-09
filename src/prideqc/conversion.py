"""Optional external vendor conversion, isolated from Python dependencies."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import Protocol


class RawConverter(Protocol):
    def convert(self, source: Path, output_directory: Path) -> Path: ...


class ExternalConverter:
    """Invoke an explicitly chosen installed converter with an argument list."""

    def __init__(self, tool: str, executable: str | None = None, timeout: float = 3600) -> None:
        if tool not in {"thermorawfileparser", "msconvert"}:
            raise ValueError(f"Unknown converter: {tool}")
        if timeout <= 0:
            raise ValueError("Conversion timeout must be positive.")
        self.tool = tool
        self.executable = executable or ("ThermoRawFileParser" if tool == "thermorawfileparser" else "msconvert")
        self.timeout = timeout

    def command(self, source: Path, output_directory: Path) -> list[str]:
        if self.tool == "thermorawfileparser":
            if source.suffix.lower() != ".raw" or not source.is_file():
                raise ValueError("ThermoRawFileParser requires a Thermo .raw file.")
            return [self.executable, f"-i={source}", f"-o={output_directory}", "-f=2"]
        return [self.executable, str(source), "--mzML", "--outdir", str(output_directory)]

    def convert(self, source: Path, output_directory: Path) -> Path:
        if shutil.which(self.executable) is None:
            raise RuntimeError(f"Converter executable not found: {self.executable}")
        # A fresh per-source directory prevents stale mzML being accepted as output.
        output_directory.mkdir(parents=True, exist_ok=False)
        command = self.command(source, output_directory)
        with (output_directory / "conversion.log").open("w", encoding="utf-8") as log:
            subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT,
                           timeout=self.timeout, shell=False)
        matches = [p for p in output_directory.iterdir() if p.name.lower().endswith(".mzml")]
        if len(matches) != 1 or matches[0].stat().st_size == 0:
            raise RuntimeError(
                "Converter must produce exactly one nonempty mzML; inspect conversion.log.",
            )
        return matches[0]
