#!/usr/bin/env python3
"""Archive lane: QC every run member of one PRIDE ZIP, one member on scratch at a time.

The task downloads the ZIP once (pridepy, checksum checked) to node-local scratch, lists the
run members with the same rules as inspect_archives.py, then for each run member:
extract it (SHA-256 while writing), analyze it, keep the receipt and mzQC, delete the member.
Scratch use is the archive plus the largest member. Each member gets its own receipt
(result.json) with the archive name, the member path, the member size and SHA-256. A ZIP
without a run member gives one fail-closed receipt, never a silent empty result.
"""
from __future__ import annotations

import hashlib
import json
import shutil
from collections.abc import Callable
from pathlib import Path

import archive_members as AM
from analyze_unit import analysis_timeout_s, analyze
from prepare_archive import _download_detail, _sha256, _transfer_records, acquire


def _write(path: Path, data: dict) -> None:
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def unit_id_for(task_id: int, index: int, member: str) -> str:
    fingerprint = hashlib.sha256(member.encode()).hexdigest()[:12]
    return f"archive-{task_id:09d}-member-{index:05d}-{fingerprint}"


def run_members(task: dict, protocol: str, scratch: Path, ready_units: Path, *, source: Path | None = None,
                checksum: bool = True, analyzer: Callable = analyze, plan_file: Path | None = None) -> dict:
    """Return the prepare-status summary; writes one result folder per member under ready_units.

    plan_file gets the member plan before the first member starts, so the supervisor can
    write fail-closed receipts for the members that a worker kill left without one.
    """
    pxd, task_id, archive_name = task["accession"], int(task["task_id"]), task["archive_file"]
    base = {"schema_version": 2, "accession": pxd, "task_id": task_id, "source_archive": archive_name,
            "expected_file_tasks": task["expected_runs"], "source_kind": "archive_member"}
    evidence: dict = {"checksum_check_requested": bool(checksum),
                      "transfer": "local_source" if source is not None else protocol}
    summary = dict(base, status="success", error_type=None, emitted_units=0, members_success=0, members_killed=0,
                   source_manifest=task.get("source_manifest"),
                   source_manifest_task_count=task["expected_runs"], source_evidence=evidence)
    transfer = scratch / "download"
    try:
        archive = acquire(pxd, archive_name, protocol, transfer, source=source, checksum=checksum)
        evidence["size_bytes"] = archive.stat().st_size
        evidence["sha256"] = _sha256(archive)
        plan = AM.list_zip(archive)
        if not plan["runs"]:
            raise AM.ArchiveMemberError("Archive contains no run member")
    except Exception as exc:  # noqa: BLE001 - every failure becomes one fail-closed receipt
        detail = _download_detail(_transfer_records(transfer))
        records = _transfer_records(transfer)
        if records is not None:
            evidence["transfer_records"] = records
        shutil.rmtree(scratch, ignore_errors=True)
        stage = "archive_preparation" if not evidence.get("sha256") else "archive_listing"
        unit_id = f"archive-{task_id:09d}-failure"
        unit = scratch.parent / f"{scratch.name}-failure-unit"
        unit.mkdir(parents=True, exist_ok=True)
        _write(unit / "unit.json", dict(base, unit_id=unit_id, discovery_status="failed", failure_stage=stage,
                                        error_type=type(exc).__name__, **detail))
        analyzer(unit, ready_units / unit_id)
        shutil.rmtree(unit, ignore_errors=True)
        summary.update(status="failed", error_type=type(exc).__name__, emitted_units=1, **detail)
        return summary
    records = _transfer_records(transfer)
    if records is not None:
        evidence["transfer_records"] = records
    archive_evidence = {"source_archive_sha256": evidence["sha256"], "source_archive_size_bytes": evidence["size_bytes"]}
    summary.update(member_plan={"runs": len(plan["runs"]), "unsupported": plan["unsupported"],
                                "ignored": plan["ignored"], "mzml_skipped_for_raw": plan["skipped_mzml_with_raw"]})
    members = [{"unit_id": unit_id_for(task_id, i, m["member"]), "member": m["member"], "member_index": i}
               for i, m in enumerate(plan["runs"])]
    if plan_file is not None:
        _write(plan_file, {"base": base, "members": members})
    for index, member in enumerate(plan["runs"]):
        unit_id = members[index]["unit_id"]
        unit = scratch / "member"
        shutil.rmtree(unit, ignore_errors=True)
        info = dict(base, **archive_evidence, unit_id=unit_id, source_member=member["member"],
                    input_basename=member["basename"], format=member["format"], member_index=index,
                    member_count=len(plan["runs"]), discovery_status="success")
        try:
            out = AM.extract_member(archive, member["member"], unit)
            info.update(source_member_sha256=out["sha256"], source_member_size_bytes=out["size"])
        except Exception as exc:  # noqa: BLE001
            shutil.rmtree(unit, ignore_errors=True)
            unit.mkdir(parents=True)
            info.update(discovery_status="failed", failure_stage="archive_member_extraction",
                        error_type=type(exc).__name__)
        _write(unit / "unit.json", info)
        timeout = analysis_timeout_s(int(info.get("source_member_size_bytes") or member["size"]))
        result = analyzer(unit, ready_units / unit_id, timeout_s=timeout)
        # analyze() copies the unit fields it knows; the member provenance is added here.
        for key in ("source_kind", "source_member_sha256", "source_member_size_bytes", "member_index",
                    "member_count", "source_archive_sha256", "source_archive_size_bytes"):
            if key in info:
                result[key] = info[key]
        _write(ready_units / unit_id / "result.json", result)
        summary["emitted_units"] += 1
        summary["members_success"] += result.get("status") == "success"
        summary["members_killed"] += result.get("error_type") == "AnalysisKilledResourceLimit"
        shutil.rmtree(unit, ignore_errors=True)
    shutil.rmtree(scratch, ignore_errors=True)
    return summary
