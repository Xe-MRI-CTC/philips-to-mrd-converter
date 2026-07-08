#!/usr/bin/env python3
"""
mrd_compare.py

Summarize one ISMRMRD/MRD file, or compare two files to produce a clean
debugging/troubleshooting report.

Features
--------
1) Single-file summary
   - file metadata
   - HDF5 layout summary
   - XML header presence / hash / flattened highlights
   - acquisition count
   - aggregate acquisition stats
   - sampled acquisition header summaries
   - image groups / waveforms overview

2) Two-file comparison
   - file metadata differences
   - HDF5 layout differences
   - XML header differences
   - aggregate acquisition differences
   - sampled per-acquisition differences
   - image and waveform structure differences

Dependencies
------------
    pip install ismrmrd h5py numpy

Usage
-----
Single-file summary:
    python mrd_compare.py file1.mrd

Two-file compare:
    python mrd_compare.py file1.mrd file2.mrd

Human readable output:
    python mrd_compare.py file1.mrd file2.mrd --nojson

Increase acquisition sampling depth:
    python mrd_compare.py file1.mrd file2.mrd --sample-acqs 25

Include sampled acquisition data hashes/statistics:
    python mrd_compare.py file1.mrd file2.mrd --compare-data
"""

from __future__ import annotations

import argparse
import dataclasses
import hashlib
import json
import math
import pathlib
import sys
import traceback
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from typing import Any, Dict, List, Optional, Sequence


# =============================================================================
# Helpers
# =============================================================================

def _lazy_import_numpy():
    import numpy as np
    return np


def _lazy_import_h5py():
    import h5py
    return h5py


def _lazy_import_ismrmrd():
    import ismrmrd
    return ismrmrd


def _is_number(x: Any) -> bool:
    return isinstance(x, (int, float)) and not isinstance(x, bool)


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _sha256_file(path: str, chunk_size: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            chunk = f.read(chunk_size)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def _format_bytes(n: int) -> str:
    if n < 1024:
        return f"{n} B"
    units = ["KB", "MB", "GB", "TB", "PB"]
    v = float(n)
    for u in units:
        v /= 1024.0
        if abs(v) < 1024.0:
            return f"{v:.2f} {u}"
    return f"{v:.2f} EB"


def _jsonable(value: Any) -> Any:
    """
    Convert objects to JSON-serializable structures.
    """
    np = None
    try:
        np = _lazy_import_numpy()
    except Exception:
        pass

    if value is None:
        return None
    if isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, bytes):
        try:
            return value.decode("utf-8", errors="replace")
        except Exception:
            return repr(value)
    if isinstance(value, pathlib.Path):
        return str(value)
    if dataclasses.is_dataclass(value):
        return {k: _jsonable(v) for k, v in dataclasses.asdict(value).items()}
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_jsonable(v) for v in value]
    if np is not None:
        if isinstance(value, np.ndarray):
            return value.tolist()
        if isinstance(value, np.generic):
            return value.item()
    try:
        return repr(value)
    except Exception:
        return "<unserializable>"


def _truncate(value: Any, max_len: int = 160) -> str:
    s = str(value)
    if len(s) <= max_len:
        return s
    return s[: max_len - 3] + "..."


def _normalize_scalar(x: Any) -> Any:
    """
    Normalize common scalar-ish values for cleaner comparison/reporting.
    """
    np = None
    try:
        np = _lazy_import_numpy()
    except Exception:
        pass

    if np is not None and isinstance(x, np.generic):
        x = x.item()

    if isinstance(x, bytes):
        try:
            x = x.decode("utf-8", errors="replace")
        except Exception:
            x = repr(x)

    if isinstance(x, float):
        if math.isfinite(x):
            return round(x, 6)
        return x
    return x


def _sample_indices(n: int, k: int) -> List[int]:
    if n <= 0 or k <= 0:
        return []
    if n <= k:
        return list(range(n))
    if k == 1:
        return [0]
    indices = sorted(set(round(i * (n - 1) / (k - 1)) for i in range(k)))
    return [int(i) for i in indices]


def _default_output_path(file_a: str, file_b: Optional[str] = None) -> pathlib.Path:
    """
    Build the default JSON output filename.

    Single-file mode:
        fileA_ismrmrd_summary.json

    Compare mode:
        fileA_vs_fileB_ismrmrd_report.json
    """
    a = pathlib.Path(file_a)
    if file_b:
        b = pathlib.Path(file_b)
        return a.with_name(f"{a.stem}_vs_{b.stem}_ismrmrd_report.json")
    return a.with_name(f"{a.stem}_ismrmrd_summary.json")


def _write_json_report(report_obj: Any, output_path: str, indent: int = 2) -> pathlib.Path:
    """
    Write a JSON report to disk and return the resolved output path.
    """
    output = pathlib.Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w", encoding="utf-8") as f:
        json.dump(_jsonable(report_obj), f, indent=indent, sort_keys=False)
    return output.resolve()


def _maybe_ctypes_array_to_list(obj: Any) -> Optional[Any]:
    """
    Convert ctypes/array-like fixed-size objects from ismrmrd headers into
    plain Python lists when possible.
    """
    try:
        return [_normalize_scalar(x) for x in obj]
    except Exception:
        return None

# =============================================================================
# XML flattening / comparison
# =============================================================================


def _strip_ns(tag: str) -> str:
    if "}" in tag:
        return tag.split("}", 1)[1]
    return tag


def _flatten_xml_element(elem: ET.Element, base_path: str = "") -> Dict[str, Any]:
    """
    Flatten an XML tree to path->value entries.

    Repeated sibling tags are indexed:
        root/encoding[0]/matrixSize/x
        root/encoding[1]/matrixSize/x

    Attributes become:
        root/node/@attr
    Text becomes:
        root/node = "value"
    """
    out: Dict[str, Any] = {}
    tag = _strip_ns(elem.tag)
    path = f"{base_path}/{tag}" if base_path else tag

    text = (elem.text or "").strip()
    if text:
        out[path] = text

    for attr_name, attr_value in elem.attrib.items():
        out[f"{path}/@{attr_name}"] = attr_value

    children = list(elem)
    if not children:
        return out

    counts = Counter(_strip_ns(c.tag) for c in children)
    seen = defaultdict(int)

    for child in children:
        child_tag = _strip_ns(child.tag)
        if counts[child_tag] > 1:
            child_path = f"{path}/{child_tag}[{seen[child_tag]}]"
            seen[child_tag] += 1
            # child_path currently includes child tag; recurse using parent path
            # and then patch the first segment so the child gets the indexed path.
            child_flat = _flatten_xml_element(child, base_path=path)
            for k, v in child_flat.items():
                prefix = f"{path}/{child_tag}"
                if k == prefix:
                    out[child_path] = v
                elif k.startswith(prefix + "/"):
                    out[child_path + k[len(prefix):]] = v
                else:
                    out[k] = v
        else:
            out.update(_flatten_xml_element(child, base_path=path))
    return out


def flatten_xml_string(xml_text: str) -> Dict[str, Any]:
    xml_text = xml_text or ""
    if not xml_text.strip():
        return {}
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError:
        return {"__parse_error__": "Unable to parse XML header"}
    return _flatten_xml_element(root)


def diff_flat_maps(a: Dict[str, Any], b: Dict[str, Any]) -> Dict[str, Any]:
    keys = sorted(set(a) | set(b))
    only_a = {}
    only_b = {}
    changed = {}
    same_count = 0
    for k in keys:
        in_a = k in a
        in_b = k in b
        if in_a and not in_b:
            only_a[k] = a[k]
        elif in_b and not in_a:
            only_b[k] = b[k]
        else:
            av = _normalize_scalar(a[k])
            bv = _normalize_scalar(b[k])
            if av == bv:
                same_count += 1
            else:
                changed[k] = {"a": av, "b": bv}
    return {
        "same_count": same_count,
        "only_a_count": len(only_a),
        "only_b_count": len(only_b),
        "changed_count": len(changed),
        "only_a": only_a,
        "only_b": only_b,
        "changed": changed,
    }


# =============================================================================
# HDF5 structure inspection
# =============================================================================

def inspect_hdf5_structure(path: str) -> Dict[str, Any]:
    h5py = _lazy_import_h5py()
    result: Dict[str, Any] = {
        "top_level_keys": [],
        "paths": {},
        "image_groups": [],
        "waveforms_present": False,
        "dataset_group_present": False,
    }

    with h5py.File(path, "r") as f:
        result["top_level_keys"] = sorted(list(f.keys()))

        def visit(name: str, obj: Any):
            entry: Dict[str, Any] = {
                "kind": "group" if hasattr(obj, "keys") else "dataset"
            }
            if hasattr(obj, "shape"):
                entry["shape"] = tuple(int(x) for x in obj.shape)
            if hasattr(obj, "dtype"):
                try:
                    entry["dtype"] = str(obj.dtype)
                except Exception:
                    entry["dtype"] = "<unknown>"
            result["paths"]["/" + name] = entry

        f.visititems(visit)

        if "dataset" in f:
            result["dataset_group_present"] = True
            dsgrp = f["dataset"]
            result["waveforms_present"] = "waveforms" in dsgrp
            for k in dsgrp.keys():
                if k.startswith("image_"):
                    result["image_groups"].append(k)
            result["image_groups"] = sorted(result["image_groups"])

    return result


# =============================================================================
# ISMRMRD acquisition helpers
# =============================================================================

def _struct_like_to_dict(obj: Any, depth: int = 3) -> Any:
    """
    Best-effort conversion of ISMRMRD objects/headers into plain Python objects.
    """
    np = None
    try:
        np = _lazy_import_numpy()
    except Exception:
        pass

    if depth < 0:
        return repr(obj)

    if obj is None or isinstance(obj, (str, int, float, bool)):
        return _normalize_scalar(obj)

    if np is not None:
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, np.generic):
            return _normalize_scalar(obj.item())

    if isinstance(obj, bytes):
        try:
            return obj.decode("utf-8", errors="replace")
        except Exception:
            return repr(obj)

    if isinstance(obj, (list, tuple)):
        return [_struct_like_to_dict(x, depth - 1) for x in obj]

    if isinstance(obj, dict):
        return {str(k): _struct_like_to_dict(v, depth - 1) for k, v in obj.items()}

    # Try attrs/properties
    names = []
    for name in dir(obj):
        if name.startswith("_"):
            continue
        if name in ("FLAGS",):
            # huge static field collections can be noisy
            continue
        try:
            value = getattr(obj, name)
        except Exception:
            continue
        if callable(value):
            continue
        names.append(name)

    if names:
        out = {}
        for name in names:
            try:
                out[name] = _struct_like_to_dict(getattr(obj, name), depth - 1)
            except Exception:
                out[name] = "<error>"
        return out

    # Try converting ctypes / fixed-size array-like objects to lists
    as_list = _maybe_ctypes_array_to_list(obj)
    if as_list is not None:
        return as_list

    return repr(obj)


def _extract_acquisition_head(acq: Any) -> Dict[str, Any]:
    """
    Convert acquisition header to a plain dict.
    """
    head = None
    if hasattr(acq, "getHead"):
        try:
            head = acq.getHead()
        except Exception:
            head = None
    if head is None and hasattr(acq, "head"):
        head = acq.head
    if head is None:
        return {}
    d = _struct_like_to_dict(head, depth=3)
    if not isinstance(d, dict):
        return {"repr": repr(d)}
    return d


def _extract_acquisition_data_stats(acq: Any) -> Dict[str, Any]:
    np = _lazy_import_numpy()
    stats: Dict[str, Any] = {}

    data = getattr(acq, "data", None)
    if data is not None:
        arr = np.asarray(data)
        stats["data_shape"] = tuple(int(x) for x in arr.shape)
        stats["data_dtype"] = str(arr.dtype)
        if arr.size > 0:
            mag = np.abs(arr)
            stats["data_abs_min"] = float(np.min(mag))
            stats["data_abs_max"] = float(np.max(mag))
            stats["data_abs_mean"] = float(np.mean(mag))
            stats["data_hash_sha256"] = _sha256_bytes(np.ascontiguousarray(arr).tobytes())
        else:
            stats["data_abs_min"] = None
            stats["data_abs_max"] = None
            stats["data_abs_mean"] = None
            stats["data_hash_sha256"] = _sha256_bytes(b"")

    traj = getattr(acq, "traj", None)
    if traj is not None:
        try:
            traj_arr = np.asarray(traj)
            stats["traj_shape"] = tuple(int(x) for x in traj_arr.shape)
            stats["traj_dtype"] = str(traj_arr.dtype)
            if traj_arr.size > 0:
                stats["traj_hash_sha256"] = _sha256_bytes(np.ascontiguousarray(traj_arr).tobytes())
        except Exception:
            stats["traj_shape"] = "<unavailable>"

    return stats


def _read_xml_header_via_h5(path: str) -> str:
    """
    Use h5py to read /dataset/xml if present.
    """
    h5py = _lazy_import_h5py()
    with h5py.File(path, "r") as f:
        if "dataset" not in f:
            return ""
        dsgrp = f["dataset"]
        if "xml" not in dsgrp:
            return ""
        node = dsgrp["xml"]
        try:
            value = node[()]
        except Exception:
            return ""

    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    if isinstance(value, str):
        return value
    try:
        return bytes(value).decode("utf-8", errors="replace")
    except Exception:
        return str(value)


def _scan_acquisitions(
    path: str,
    sample_acqs: int = 10,
    max_acqs_to_scan: Optional[int] = None,
    compare_data: bool = False,
) -> Dict[str, Any]:
    """
    Read acquisition count and summarize a scan of acquisitions.

    If max_acqs_to_scan is None, scan all acquisitions.
    Otherwise, scan up to that many acquisitions (evenly sampled).
    """
    ismrmrd = _lazy_import_ismrmrd()
    dset = ismrmrd.Dataset(path)

    result: Dict[str, Any] = {
        "number_of_acquisitions": None,
        "scan_mode": "all" if max_acqs_to_scan is None else f"sampled<= {max_acqs_to_scan}",
        "scanned_count": 0,
        "sampled_acquisitions": [],
        "aggregates": {},
        "errors": [],
    }

    try:
        nacq = int(dset.number_of_acquisitions())
        result["number_of_acquisitions"] = nacq

        if max_acqs_to_scan is None or nacq <= max_acqs_to_scan:
            scan_indices = list(range(nacq))
        else:
            scan_indices = _sample_indices(nacq, max_acqs_to_scan)

        sample_indices = set(_sample_indices(len(scan_indices), sample_acqs))
        actual_sampled_positions = set(scan_indices[i] for i in sample_indices) if scan_indices else set()

        result["scanned_count"] = len(scan_indices)

        # Aggregate counters
        number_of_samples = []
        active_channels = []
        trajectory_dimensions = []
        encoding_space_refs = []
        flags_counter = Counter()
        data_shapes = Counter()
        traj_shapes = Counter()
        idx_ranges: Dict[str, set] = defaultdict(set)

        for i in scan_indices:
            try:
                acq = dset.read_acquisition(i)
                head = _extract_acquisition_head(acq)

                # Best effort field extraction
                ns = head.get("number_of_samples")
                ac = head.get("active_channels")
                td = head.get("trajectory_dimensions")
                esr = head.get("encoding_space_ref")
                flg = head.get("flags")
                idx = head.get("idx", {})

                if _is_number(ns):
                    number_of_samples.append(int(ns))
                if _is_number(ac):
                    active_channels.append(int(ac))
                if _is_number(td):
                    trajectory_dimensions.append(int(td))
                if _is_number(esr):
                    encoding_space_refs.append(int(esr))
                if _is_number(flg):
                    flags_counter[int(flg)] += 1

                if isinstance(idx, dict):
                    for k, v in idx.items():
                        if isinstance(v, list):
                            # user array inside idx; skip
                            continue
                        if _is_number(v):
                            idx_ranges[k].add(int(v))

                # Shape summaries
                try:
                    data = getattr(acq, "data", None)
                    if data is not None:
                        data_shapes[tuple(int(x) for x in data.shape)] += 1
                except Exception:
                    pass

                try:
                    traj = getattr(acq, "traj", None)
                    if traj is not None:
                        traj_shapes[tuple(int(x) for x in traj.shape)] += 1
                except Exception:
                    pass

                # Sampled acquisition details
                if i in actual_sampled_positions:
                    sample_entry = {
                        "index": i,
                        "header": head,
                    }
                    if compare_data:
                        try:
                            sample_entry["data_stats"] = _extract_acquisition_data_stats(acq)
                        except Exception as ex:
                            sample_entry["data_stats_error"] = repr(ex)
                    else:
                        try:
                            data = getattr(acq, "data", None)
                            if data is not None:
                                sample_entry["data_shape"] = tuple(int(x) for x in data.shape)
                                sample_entry["data_dtype"] = str(data.dtype)
                        except Exception:
                            pass
                        try:
                            traj = getattr(acq, "traj", None)
                            if traj is not None:
                                sample_entry["traj_shape"] = tuple(int(x) for x in traj.shape)
                        except Exception:
                            pass

                    result["sampled_acquisitions"].append(sample_entry)

            except Exception as ex:
                result["errors"].append({
                    "acquisition_index": i,
                    "error": repr(ex),
                })

        def summarize_numeric(values: List[int]) -> Dict[str, Any]:
            if not values:
                return {}
            return {
                "count": len(values),
                "min": min(values),
                "max": max(values),
                "mean": round(sum(values) / len(values), 6),
                "unique": sorted(set(values)),
            }

        result["aggregates"] = {
            "number_of_samples": summarize_numeric(number_of_samples),
            "active_channels": summarize_numeric(active_channels),
            "trajectory_dimensions": summarize_numeric(trajectory_dimensions),
            "encoding_space_ref": summarize_numeric(encoding_space_refs),
            "flags_histogram": dict(sorted(flags_counter.items())),
            "data_shapes": {str(k): v for k, v in sorted(data_shapes.items(), key=lambda kv: str(kv[0]))},
            "traj_shapes": {str(k): v for k, v in sorted(traj_shapes.items(), key=lambda kv: str(kv[0]))},
            "idx_ranges": {
                k: {
                    "min": min(v) if v else None,
                    "max": max(v) if v else None,
                    "unique_count": len(v),
                    # "sample_values": sorted(v)[:20],
                }
                for k, v in sorted(idx_ranges.items())
                if v
            },
        }
        return result
    finally:
        try:
            dset.close()
        except Exception:
            pass


# =============================================================================
# Image and waveform summaries via HDF5
# =============================================================================

def summarize_images_and_waveforms(path: str) -> Dict[str, Any]:
    h5py = _lazy_import_h5py()
    result: Dict[str, Any] = {
        "image_groups": {},
        "waveforms": {},
    }

    with h5py.File(path, "r") as f:
        dsgrp = f.get("dataset")
        if dsgrp is None:
            return result

        # Image groups
        for key in sorted(dsgrp.keys()):
            if not key.startswith("image_"):
                continue
            grp = dsgrp[key]
            entry = {
                "keys": sorted(list(grp.keys())),
                "data_shape": None,
                "data_dtype": None,
                "header_shape": None,
                "header_dtype": None,
                "attributes_shape": None,
                "attributes_dtype": None,
            }
            if "data" in grp:
                entry["data_shape"] = tuple(int(x) for x in grp["data"].shape)
                entry["data_dtype"] = str(grp["data"].dtype)
            if "header" in grp:
                entry["header_shape"] = tuple(int(x) for x in grp["header"].shape)
                entry["header_dtype"] = str(grp["header"].dtype)
            if "attributes" in grp:
                entry["attributes_shape"] = tuple(int(x) for x in grp["attributes"].shape)
                entry["attributes_dtype"] = str(grp["attributes"].dtype)
            result["image_groups"][key] = entry

        # Waveforms
        if "waveforms" in dsgrp:
            wf = dsgrp["waveforms"]
            result["waveforms"] = {
                "shape": tuple(int(x) for x in wf.shape),
                "dtype": str(wf.dtype),
            }

    return result


# =============================================================================
# Summary / compare API
# =============================================================================

def summarize_ismrmrd_file(
    path: str,
    sample_acqs: int = 2,
    max_acqs_to_scan: Optional[int] = None,
    compare_data: bool = False,
) -> Dict[str, Any]:
    """
    Build a clean summary for one ISMRMRD/MRD file.
    """
    p = pathlib.Path(path)
    if not p.exists():
        raise FileNotFoundError(path)

    summary: Dict[str, Any] = {
        "path": str(p.resolve()),
        "file_name": p.name,
        "suffix": p.suffix,
        "file_size_bytes": p.stat().st_size,
        "file_size_human": _format_bytes(p.stat().st_size),
        "file_sha256": _sha256_file(str(p)),
    }
    # Read in ISMRMRD File
    ismrmrd = _lazy_import_ismrmrd()
    dset = ismrmrd.Dataset(path)

    # HDF5 structure
    try:
        summary["hdf5"] = inspect_hdf5_structure(str(p))
    except Exception as ex:
        summary["hdf5_error"] = "".join(traceback.format_exception_only(type(ex), ex)).strip()

    # XML header
    try:
        header = ismrmrd.xsd.CreateFromDocument(dset.read_xml_header())
        header_xml = ismrmrd.xsd.ToXML(header)

        summary["xml_header"] = {
            "present": bool(header_xml.strip()),
            "length": len(header_xml),
            "sha256": _sha256_bytes(header_xml.encode("utf-8")) if header_xml else None,
            "flattened": flatten_xml_string(header_xml),
        }
    except Exception as ex:
        summary["xml_header_error"] = "".join(traceback.format_exception_only(type(ex), ex)).strip()

    # Acquisition-level summary
    try:
        summary["acquisitions"] = _scan_acquisitions(
            str(p),
            sample_acqs=sample_acqs,
            max_acqs_to_scan=max_acqs_to_scan,
            compare_data=compare_data,
        )
    except Exception as ex:
        summary["acquisitions_error"] = "".join(traceback.format_exception_only(type(ex), ex)).strip()

    # Images / waveforms
    try:
        summary["images_and_waveforms"] = summarize_images_and_waveforms(str(p))
    except Exception as ex:
        summary["images_and_waveforms_error"] = "".join(traceback.format_exception_only(type(ex), ex)).strip()

    return summary


def _deep_diff(a: Any, b: Any, path: str = "") -> List[Dict[str, Any]]:
    """
    Generic recursive diff for JSON-like structures.
    """
    diffs: List[Dict[str, Any]] = []

    if type(a) is not type(b):
        diffs.append({
            "path": path or "/",
            "type": "type_mismatch",
            "a_type": type(a).__name__,
            "b_type": type(b).__name__,
            "a": _jsonable(a),
            "b": _jsonable(b),
        })
        return diffs

    if isinstance(a, dict):
        keys = sorted(set(a) | set(b))
        for k in keys:
            kp = f"{path}/{k}" if path else f"/{k}"
            if k not in a:
                diffs.append({"path": kp, "type": "only_in_b", "b": _jsonable(b[k])})
            elif k not in b:
                diffs.append({"path": kp, "type": "only_in_a", "a": _jsonable(a[k])})
            else:
                diffs.extend(_deep_diff(a[k], b[k], kp))
        return diffs

    if isinstance(a, list):
        if len(a) != len(b):
            diffs.append({
                "path": path or "/",
                "type": "list_length_mismatch",
                "a_len": len(a),
                "b_len": len(b),
            })
        for i, (av, bv) in enumerate(zip(a, b)):
            diffs.extend(_deep_diff(av, bv, f"{path}[{i}]"))
        return diffs

    av = _normalize_scalar(a)
    bv = _normalize_scalar(b)
    if av != bv:
        diffs.append({
            "path": path or "/",
            "type": "value_mismatch",
            "a": _jsonable(av),
            "b": _jsonable(bv),
        })
    return diffs


def compare_ismrmrd_files(
    path_a: str,
    path_b: str,
    sample_acqs: int = 2,
    max_acqs_to_scan: Optional[int] = None,
    compare_data: bool = False,
) -> Dict[str, Any]:
    """
    Compare two ISMRMRD/MRD files and return a structured diff report.
    """
    a = summarize_ismrmrd_file(
        path_a,
        sample_acqs=sample_acqs,
        max_acqs_to_scan=max_acqs_to_scan,
        compare_data=compare_data,
    )
    b = summarize_ismrmrd_file(
        path_b,
        sample_acqs=sample_acqs,
        max_acqs_to_scan=max_acqs_to_scan,
        compare_data=compare_data,
    )

    report: Dict[str, Any] = {
        "file_a": a["path"],
        "file_b": b["path"],
        "same_file_sha256": a.get("file_sha256") == b.get("file_sha256"),
        "file_level": {},
        "xml_header_diff": {},
        "hdf5_diff": {},
        "acquisition_diff": {},
        "images_and_waveforms_diff": {},
        "all_diffs": {},
    }

    # File-level compare (exclude heavy nested details)
    file_level_keys = [
        "file_name",
        "suffix",
        "file_size_bytes",
        "file_sha256",
    ]
    file_level = {}
    for k in file_level_keys:
        av = a.get(k)
        bv = b.get(k)
        if av != bv:
            file_level[k] = {"a": av, "b": bv}
    report["file_level"] = file_level

    # XML header compare
    ax = a.get("xml_header", {}).get("flattened", {})
    bx = b.get("xml_header", {}).get("flattened", {})
    report["xml_header_diff"] = diff_flat_maps(ax, bx)

    # HDF5 compare
    ah = a.get("hdf5", {})
    bh = b.get("hdf5", {})
    report["hdf5_diff"] = _deep_diff(ah, bh)

    # Acquisition compare
    aa = a.get("acquisitions", {})
    ba = b.get("acquisitions", {})

    acq_focus = {
        "number_of_acquisitions": {
            "a": aa.get("number_of_acquisitions"),
            "b": ba.get("number_of_acquisitions"),
        },
        "aggregates_diff": _deep_diff(aa.get("aggregates", {}), ba.get("aggregates", {})),
        "sampled_acquisitions_diff": _deep_diff(
            aa.get("sampled_acquisitions", []),
            ba.get("sampled_acquisitions", []),
        ),
        "scan_errors_diff": _deep_diff(
            aa.get("errors", []),
            ba.get("errors", []),
        ),
    }
    report["acquisition_diff"] = acq_focus

    # Images/waveforms compare
    aiw = a.get("images_and_waveforms", {})
    biw = b.get("images_and_waveforms", {})
    report["images_and_waveforms_diff"] = _deep_diff(aiw, biw)

    # Broad recursive diff (excluding some large raw maps to avoid duplication)
    a_small = {
        k: v for k, v in a.items()
        if k not in {"xml_header"}
    }
    b_small = {
        k: v for k, v in b.items()
        if k not in {"xml_header"}
    }
    report["all_diffs"] = {
        "summary_diff": _deep_diff(a_small, b_small)
    }

    return report


# =============================================================================
# Human-readable formatting
# =============================================================================

def _render_kv_lines(d: Dict[str, Any], indent: int = 0) -> List[str]:
    lines = []
    prefix = " " * indent
    for k in sorted(d.keys()):
        v = d[k]
        if isinstance(v, dict):
            lines.append(f"{prefix}{k}:")
            lines.extend(_render_kv_lines(v, indent + 2))
        elif isinstance(v, list):
            lines.append(f"{prefix}{k}:")
            if not v:
                lines.append(f"{prefix}  []")
            else:
                for item in v:
                    if isinstance(item, (dict, list)):
                        lines.append(f"{prefix}  -")
                        if isinstance(item, dict):
                            lines.extend(_render_kv_lines(item, indent + 4))
                        else:
                            for sub in item:
                                lines.append(f"{prefix}    - {_truncate(sub)}")
                    else:
                        lines.append(f"{prefix}  - {_truncate(item)}")
        else:
            lines.append(f"{prefix}{k}: {_truncate(v)}")
    return lines


def format_summary_text(summary: Dict[str, Any], show_xml_keys: int = 50) -> str:
    lines: List[str] = []
    lines.append("ISMRMRD FILE SUMMARY")
    lines.append("=" * 80)
    lines.append(f"Path:             {summary.get('path')}")
    lines.append(f"File name:        {summary.get('file_name')}")
    lines.append(f"Suffix:           {summary.get('suffix')}")
    lines.append(f"File size:        {summary.get('file_size_bytes')} bytes ({summary.get('file_size_human')})")
    lines.append(f"File SHA-256:     {summary.get('file_sha256')}")
    lines.append("")

    # HDF5
    if "hdf5" in summary:
        h = summary["hdf5"]
        lines.append("HDF5 STRUCTURE")
        lines.append("-" * 80)
        lines.append(f"Top-level keys:   {h.get('top_level_keys', [])}")
        lines.append(f"Has /dataset:     {h.get('dataset_group_present')}")
        lines.append(f"Image groups:     {h.get('image_groups', [])}")
        lines.append(f"Waveforms present:{h.get('waveforms_present')}")
        path_items = sorted(h.get("paths", {}).items())
        lines.append(f"Total paths:      {len(path_items)}")
        for path, entry in path_items[:40]:
            kind = entry.get("kind")
            shape = entry.get("shape")
            dtype = entry.get("dtype")
            lines.append(f"  {path} [{kind}] shape={shape} dtype={dtype}")
        if len(path_items) > 40:
            lines.append(f"  ... ({len(path_items) - 40} more paths)")
        lines.append("")
    elif "hdf5_error" in summary:
        lines.append("HDF5 STRUCTURE")
        lines.append("-" * 80)
        lines.append(f"Error: {summary['hdf5_error']}")
        lines.append("")

    # XML
    if "xml_header" in summary:
        x = summary["xml_header"]
        lines.append("XML HEADER")
        lines.append("-" * 80)
        lines.append(f"Present:          {x.get('present')}")
        lines.append(f"Length:           {x.get('length')}")
        lines.append(f"SHA-256:          {x.get('sha256')}")
        flat = x.get("flattened", {})
        lines.append(f"Flattened keys:   {len(flat)}")
        for k in sorted(flat.keys())[:show_xml_keys]:
            lines.append(f"  {k}: {_truncate(flat[k])}")
        if len(flat) > show_xml_keys:
            lines.append(f"  ... ({len(flat) - show_xml_keys} more XML keys)")
        lines.append("")
    elif "xml_header_error" in summary:
        lines.append("XML HEADER")
        lines.append("-" * 80)
        lines.append(f"Error: {summary['xml_header_error']}")
        lines.append("")

    # Acquisitions
    if "acquisitions" in summary:
        a = summary["acquisitions"]
        lines.append("ACQUISITIONS")
        lines.append("-" * 80)
        lines.append(f"Count:            {a.get('number_of_acquisitions')}")
        lines.append(f"Scan mode:        {a.get('scan_mode')}")
        lines.append(f"Scanned count:    {a.get('scanned_count', 0)}")
        agg = a.get("aggregates", {})
        lines.append("Aggregates:")
        lines.extend(_render_kv_lines(agg, indent=2))
        sampled = a.get("sampled_acquisitions", [])
        lines.append(f"Sampled acquisitions: {len(sampled)}")
        for item in sampled:
            lines.append(f"  Acquisition index {item.get('index')}:")
            header = item.get("header", {})
            for hk in [
                "flags",
                "measurement_uid",
                "scan_counter",
                "number_of_samples",
                "available_channels",
                "active_channels",
                "discard_pre",
                "discard_post",
                "center_sample",
                "encoding_space_ref",
                "trajectory_dimensions",
                "sample_time_us",
            ]:
                if hk in header:
                    lines.append(f"    {hk}: {header[hk]}")
            idx = header.get("idx")
            if isinstance(idx, dict):
                lines.append(f"    idx: {idx}")
            if "data_shape" in item:
                lines.append(f"    data_shape: {item['data_shape']}")
            if "data_dtype" in item:
                lines.append(f"    data_dtype: {item['data_dtype']}")
            if "traj_shape" in item:
                lines.append(f"    traj_shape: {item['traj_shape']}")
            if "data_stats" in item:
                lines.append(f"    data_stats: {item['data_stats']}")
        errs = a.get("errors", [])
        if errs:
            lines.append("Acquisition scan errors:")
            for err in errs[:20]:
                lines.append(f"  idx={err.get('acquisition_index')} error={err.get('error')}")
            if len(errs) > 20:
                lines.append(f"  ... ({len(errs) - 20} more)")
        lines.append("")
    elif "acquisitions_error" in summary:
        lines.append("ACQUISITIONS")
        lines.append("-" * 80)
        lines.append(f"Error: {summary['acquisitions_error']}")
        lines.append("")

    # Images / waveforms
    if "images_and_waveforms" in summary:
        iw = summary["images_and_waveforms"]
        lines.append("IMAGES AND WAVEFORMS")
        lines.append("-" * 80)
        lines.append("Image groups:")
        if iw.get("image_groups"):
            for k, v in sorted(iw["image_groups"].items()):
                lines.append(f"  {k}:")
                lines.extend(_render_kv_lines(v, indent=4))
        else:
            lines.append("  None")
        lines.append("Waveforms:")
        if iw.get("waveforms"):
            lines.extend(_render_kv_lines(iw["waveforms"], indent=2))
        else:
            lines.append("  None")
        lines.append("")
    elif "images_and_waveforms_error" in summary:
        lines.append("IMAGES AND WAVEFORMS")
        lines.append("-" * 80)
        lines.append(f"Error: {summary['images_and_waveforms_error']}")
        lines.append("")

    return "\n".join(lines)


def format_compare_text(report: Dict[str, Any], max_items: int = 100) -> str:
    lines: List[str] = []
    lines.append("ISMRMRD FILE COMPARISON")
    lines.append("=" * 80)
    lines.append(f"File A:           {report.get('file_a')}")
    lines.append(f"File B:           {report.get('file_b')}")
    lines.append(f"Same file SHA-256:{report.get('same_file_sha256')}")
    lines.append("")

    # File-level
    lines.append("FILE-LEVEL DIFFERENCES")
    lines.append("-" * 80)
    file_level = report.get("file_level", {})
    if file_level:
        lines.extend(_render_kv_lines(file_level, indent=2))
    else:
        lines.append("  None")
    lines.append("")

    # XML diff
    x = report.get("xml_header_diff", {})
    lines.append("XML HEADER DIFFERENCES")
    lines.append("-" * 80)
    lines.append(f"same_count:       {x.get('same_count')}")
    lines.append(f"only_a_count:     {x.get('only_a_count')}")
    lines.append(f"only_b_count:     {x.get('only_b_count')}")
    lines.append(f"changed_count:    {x.get('changed_count')}")
    if x.get("only_a"):
        lines.append("Only in A:")
        for i, (k, v) in enumerate(sorted(x["only_a"].items())):
            if i >= max_items:
                lines.append(f"  ... ({len(x['only_a']) - max_items} more)")
                break
            lines.append(f"  {k}: {_truncate(v)}")
    if x.get("only_b"):
        lines.append("Only in B:")
        for i, (k, v) in enumerate(sorted(x["only_b"].items())):
            if i >= max_items:
                lines.append(f"  ... ({len(x['only_b']) - max_items} more)")
                break
            lines.append(f"  {k}: {_truncate(v)}")
    if x.get("changed"):
        lines.append("Changed:")
        for i, (k, v) in enumerate(sorted(x["changed"].items())):
            if i >= max_items:
                lines.append(f"  ... ({len(x['changed']) - max_items} more)")
                break
            lines.append(f"  {k}:")
            lines.append(f"    A: {_truncate(v.get('a'))}")
            lines.append(f"    B: {_truncate(v.get('b'))}")
    if not x:
        lines.append("  None")
    lines.append("")

    # HDF5 diff
    lines.append("HDF5 STRUCTURE DIFFERENCES")
    lines.append("-" * 80)
    hdf5_diffs = report.get("hdf5_diff", [])
    if hdf5_diffs:
        for item in hdf5_diffs[:max_items]:
            lines.append(f"  {item.get('type')} at {item.get('path')}: "
                         f"A={_truncate(item.get('a'))} B={_truncate(item.get('b'))}")
        if len(hdf5_diffs) > max_items:
            lines.append(f"  ... ({len(hdf5_diffs) - max_items} more)")
    else:
        lines.append("  None")
    lines.append("")

    # Acquisition diff
    lines.append("ACQUISITION DIFFERENCES")
    lines.append("-" * 80)
    acq = report.get("acquisition_diff", {})
    if acq:
        lines.append(f"number_of_acquisitions: {acq.get('number_of_acquisitions')}")
        agd = acq.get("aggregates_diff", [])
        lines.append(f"aggregate_diffs:       {len(agd)}")
        for item in agd[:max_items]:
            lines.append(f"  {item.get('type')} at {item.get('path')}: "
                         f"A={_truncate(item.get('a'))} B={_truncate(item.get('b'))}")
        if len(agd) > max_items:
            lines.append(f"  ... ({len(agd) - max_items} more)")

        sad = acq.get("sampled_acquisitions_diff", [])
        lines.append(f"sampled_acq_diffs:     {len(sad)}")
        for item in sad[:max_items]:
            lines.append(f"  {item.get('type')} at {item.get('path')}: "
                         f"A={_truncate(item.get('a'))} B={_truncate(item.get('b'))}")
        if len(sad) > max_items:
            lines.append(f"  ... ({len(sad) - max_items} more)")

        sed = acq.get("scan_errors_diff", [])
        lines.append(f"scan_error_diffs:      {len(sed)}")
        for item in sed[:max_items]:
            lines.append(f"  {item.get('type')} at {item.get('path')}: "
                         f"A={_truncate(item.get('a'))} B={_truncate(item.get('b'))}")
        if len(sed) > max_items:
            lines.append(f"  ... ({len(sed) - max_items} more)")
    else:
        lines.append("  None")
    lines.append("")

    # Images/waveforms diff
    lines.append("IMAGES / WAVEFORMS DIFFERENCES")
    lines.append("-" * 80)
    iwd = report.get("images_and_waveforms_diff", [])
    if iwd:
        for item in iwd[:max_items]:
            lines.append(f"  {item.get('type')} at {item.get('path')}: "
                         f"A={_truncate(item.get('a'))} B={_truncate(item.get('b'))}")
        if len(iwd) > max_items:
            lines.append(f"  ... ({len(iwd) - max_items} more)")
    else:
        lines.append("  None")
    lines.append("")

    # Broad diff count
    summary_diff = report.get("all_diffs", {}).get("summary_diff", [])
    lines.append("OVERALL DIFF SUMMARY")
    lines.append("-" * 80)
    lines.append(f"total_recursive_diffs: {len(summary_diff)}")
    for item in summary_diff[:max_items]:
        lines.append(f"  {item.get('type')} at {item.get('path')}: "
                     f"A={_truncate(item.get('a'))} B={_truncate(item.get('b'))}")
    if len(summary_diff) > max_items:
        lines.append(f"  ... ({len(summary_diff) - max_items} more)")
    lines.append("")

    return "\n".join(lines)


# =============================================================================
# CLI
# =============================================================================

def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Summarize one ISMRMRD/MRD file or compare two files."
    )
    p.add_argument(
        "file_a",
        help="Path to the first .mrd / ISMRMRD file",
    )
    p.add_argument(
        "file_b",
        nargs="?",
        help="Optional second .mrd / ISMRMRD file to compare against file_a",
    )
    p.add_argument(
        "--nojson",
        action="store_false",
        help="Emit human-readable text instead of JSON",
    )
    p.add_argument(
        "--sample-acqs",
        type=int,
        default=1,
        help="Number of acquisitions to include in sampled details (default: 1)",
    )
    p.add_argument(
        "--max-acqs-to-scan",
        type=int,
        default=None,
        help=(
            "Maximum number of acquisitions to scan when building aggregates. "
            "Default: scan all acquisitions."
        ),
    )
    p.add_argument(
        "--compare-data",
        action="store_true",
        help=(
            "For sampled acquisitions, include data/traj hashes and simple magnitude stats. "
            "Useful for deeper debugging, but can be slower."
        ),
    )
    p.add_argument(
        "--indent",
        type=int,
        default=2,
        help="JSON indentation (default: 2)",
    )
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_arg_parser().parse_args(argv)

    try:
        if args.file_b:
            report = compare_ismrmrd_files(
                args.file_a,
                args.file_b,
                sample_acqs=args.sample_acqs,
                max_acqs_to_scan=args.max_acqs_to_scan,
                compare_data=args.compare_data,
            )
            if args.nojson:
                print(format_compare_text(report))
            else:
                print(json.dumps(_jsonable(report), indent=args.indent, sort_keys=False))

            written = _write_json_report(
                report,
                str(_default_output_path(args.file_a, args.file_b)),
                indent=args.indent,
            )
            print(f"JSON report written to: {written}", file=sys.stderr)

        else:
            summary = summarize_ismrmrd_file(
                args.file_a,
                sample_acqs=args.sample_acqs,
                max_acqs_to_scan=args.max_acqs_to_scan,
                compare_data=args.compare_data,
            )
            if args.nojson:
                print(format_summary_text(summary))
            else:
                print(json.dumps(_jsonable(summary), indent=args.indent, sort_keys=False))

            written = _write_json_report(
                summary,
                str(_default_output_path(args.file_a)),
                indent=args.indent,
            )
            print(f"JSON report written to: {written}", file=sys.stderr)

        return 0

    except KeyboardInterrupt:
        print("Interrupted.", file=sys.stderr)
        return 130
    except Exception as ex:
        print(f"ERROR: {ex}", file=sys.stderr)
        tb = traceback.format_exc()
        print(tb, file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
