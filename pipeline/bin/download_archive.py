#!/usr/bin/env python3
"""Stage a single member file from a multi-GWAS remote archive URL (<archive_url>#<member>).

Many-trait archives (for example Zhang et al. 2022 ARIC EA.zip and AA.zip, 4,657 proteins each) are downloaded once per
study into a shared cache folder (<study>/raw/_archives/) under an exclusive file lock, then each row extracts only its
own member in parallel with CRC/integrity verification.

CLI helper for manifest builders:
    download_archive.py --list-members https://jh-pwas.s3.amazonaws.com/results/EA.zip --out members.txt
"""

from __future__ import annotations

import argparse
import shutil
import struct
import tarfile
import tempfile
import urllib.request
import zipfile
from pathlib import Path
from typing import Any

from gwas_pipeline_common import (
    atomic_write_text,
    data_error,
    download_url,
    ensure_dir,
    file_lock,
    guess_filename_from_url,
    sanitize_filename,
    short_hash,
    split_archive_url,
    strip_url_secrets,
    urlopen,
    write_json,
)


def _extract_zip_member(archive_path: Path, member: str, dest: Path, archive_label: str) -> dict[str, Any]:
    try:
        zf_ctx = zipfile.ZipFile(archive_path, "r")
    except zipfile.BadZipFile as exc:
        raise IOError(f"Corrupt ZIP archive {archive_label}: {exc}") from exc
    with zf_ctx as zf:
        info = None
        try:
            info = zf.getinfo(member)
        except KeyError:
            target = member.lstrip("/")
            for candidate in zf.infolist():
                if candidate.filename.lstrip("/") == target:
                    info = candidate
                    break
        if info is None or info.is_dir():
            raise data_error(f"Archive member '{member}' not found in {archive_label}")
        part = dest.with_name(dest.name + ".part")
        part.unlink(missing_ok=True)
        try:
            with zf.open(info, "r") as src, open(part, "wb") as dst:
                shutil.copyfileobj(src, dst)
            if part.stat().st_size == 0:
                raise data_error(f"Archive member '{member}' in {archive_label} is empty")
            part.replace(dest)
        except zipfile.BadZipFile as exc:
            part.unlink(missing_ok=True)
            raise IOError(f"CRC/decompression check failed for '{member}' in {archive_label}: {exc}") from exc
        except Exception:
            part.unlink(missing_ok=True)
            raise
        return {"crc32": f"{info.CRC:08x}", "uncompressed_bytes": info.file_size}


def _extract_tar_member(archive_path: Path, member: str, dest: Path, archive_label: str) -> dict[str, Any]:
    try:
        tf_ctx = tarfile.open(archive_path, "r:*")
    except tarfile.TarError as exc:
        raise IOError(f"Corrupt tar archive {archive_label}: {exc}") from exc
    with tf_ctx as tf:
        info = None
        target = member.lstrip("./")
        for candidate in tf.getmembers():
            if candidate.isfile() and candidate.name.lstrip("./") == target:
                info = candidate
                break
        if info is None:
            raise data_error(f"Archive member '{member}' not found in {archive_label}")
        part = dest.with_name(dest.name + ".part")
        part.unlink(missing_ok=True)
        try:
            src = tf.extractfile(info)
            if src is None:
                raise data_error(f"Archive member '{member}' in {archive_label} is not a regular file")
            with src, open(part, "wb") as dst:
                shutil.copyfileobj(src, dst)
            if part.stat().st_size == 0:
                raise data_error(f"Archive member '{member}' in {archive_label} is empty")
            part.replace(dest)
        except Exception:
            part.unlink(missing_ok=True)
            raise
        return {"uncompressed_bytes": info.size}


def download_archive_member(location: str, out_dir: str | Path, cache_dir: str | Path | None = None) -> dict[str, Any]:
    """Ensure the parent archive is cached once under a file lock, then extract `member` into `out_dir`."""
    archive_url, member = split_archive_url(location)
    if not member:
        raise ValueError(f"Location is not an '<archive_url>#<member>' reference: {strip_url_secrets(location)}")

    out_dir = ensure_dir(out_dir)
    clean_archive_url = strip_url_secrets(archive_url)
    archive_base = sanitize_filename(guess_filename_from_url(archive_url, "archive.zip"))
    cached_name = f"{short_hash(clean_archive_url, 8)}_{archive_base}"
    archives_dir = ensure_dir(cache_dir if cache_dir is not None else out_dir.parent / "_archives")
    archive_path = archives_dir / cached_name
    lock_path = archives_dir / f"{cached_name}.lock"

    lower_name = archive_base.lower()
    is_zip = lower_name.endswith(".zip")

    with file_lock(lock_path):
        if not archive_path.exists() or archive_path.stat().st_size == 0:
            download_url(archive_url, archives_dir, filename=cached_name)
        valid = zipfile.is_zipfile(archive_path) if is_zip else tarfile.is_tarfile(archive_path)
        if not valid:
            archive_path.unlink(missing_ok=True)
            raise IOError(f"Downloaded file from {clean_archive_url} is not a valid archive")

    dest = out_dir / sanitize_filename(Path(member).name)
    if not dest.exists() or dest.stat().st_size == 0:
        if is_zip:
            check_info = _extract_zip_member(archive_path, member, dest, clean_archive_url)
        else:
            check_info = _extract_tar_member(archive_path, member, dest, clean_archive_url)
    else:
        check_info = {"uncompressed_bytes": dest.stat().st_size}

    return {
        "raw_path": str(dest),
        "archive_url": clean_archive_url,
        "archive_path": str(archive_path),
        "archive_member": member,
        **check_info,
    }


def _list_remote_zip_via_range(url: str) -> list[str] | None:
    """Read a remote .zip Central Directory using HTTP Range requests (~few hundred KB instead of the whole archive).
    Returns None if the server does not support Range or uses ZIP64."""
    head_req = urllib.request.Request(url, method="HEAD")
    with urlopen(head_req, timeout=30) as resp:
        cl = resp.headers.get("Content-Length")
        if not cl:
            return None
        size = int(cl)
    if size < 22:
        return None
    tail_len = min(size, 2 * 1024 * 1024)
    tail_req = urllib.request.Request(url, headers={"Range": f"bytes={size - tail_len}-{size - 1}"})
    with urlopen(tail_req, timeout=60) as resp:
        if getattr(resp, "status", 200) != 206 and tail_len < size:
            return None
        tail = resp.read()
    eocd_pos = tail.rfind(b"PK\x05\x06")
    if eocd_pos < 0 or len(tail) - eocd_pos < 22:
        return None
    _, _, _, _, total_entries, cd_size, cd_offset, _ = struct.unpack_from("<4sHHHHIIH", tail, eocd_pos)
    if total_entries == 0xFFFF or cd_size == 0xFFFFFFFF or cd_offset == 0xFFFFFFFF:
        return None
    if cd_offset >= size - tail_len:
        rel = cd_offset - (size - tail_len)
        cd = tail[rel : rel + cd_size]
    else:
        cd_req = urllib.request.Request(url, headers={"Range": f"bytes={cd_offset}-{cd_offset + cd_size - 1}"})
        with urlopen(cd_req, timeout=60) as resp:
            if getattr(resp, "status", 200) != 206:
                return None
            cd = resp.read()
    if len(cd) < cd_size:
        return None
    members: list[str] = []
    pos = 0
    while pos + 46 <= len(cd):
        if cd[pos : pos + 4] != b"PK\x01\x02":
            return None
        fname_len, extra_len, comment_len = struct.unpack_from("<HHH", cd, pos + 28)
        name = cd[pos + 46 : pos + 46 + fname_len].decode("utf-8", errors="replace")
        if not name.endswith("/"):
            members.append(name)
        pos += 46 + fname_len + extra_len + comment_len
    return members


def list_zip_members(url_or_path: str) -> list[str]:
    """Return all non-directory member paths inside a local or HTTP(S) .zip archive."""
    text = str(url_or_path).strip()
    local_candidate = Path(text.replace("file://", ""))
    if not text.lower().startswith(("http://", "https://", "ftp://")) and local_candidate.exists():
        with zipfile.ZipFile(local_candidate, "r") as zf:
            return [info.filename for info in zf.infolist() if not info.is_dir()]

    try:
        ranged = _list_remote_zip_via_range(text)
        if ranged is not None:
            return ranged
    except Exception:
        pass

    with tempfile.TemporaryDirectory(prefix="gwas_zip_list_") as tmp:
        downloaded = download_url(text, tmp, filename="listing.zip")
        with zipfile.ZipFile(downloaded, "r") as zf:
            return [info.filename for info in zf.infolist() if not info.is_dir()]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--list-members", default="", help="List member paths inside a local or remote .zip archive")
    parser.add_argument("--gwas-location", default="", help="Stage a single '<archive_url>#<member>' location")
    parser.add_argument("--out-dir", default="", help="Target raw directory when staging a member")
    parser.add_argument("--metadata", default="", help="JSON metadata output path when staging a member")
    parser.add_argument("--out", default="", help="Output newline-separated file list for --list-members")
    args = parser.parse_args()

    if args.list_members:
        if not args.out:
            parser.error("--out is required with --list-members")
        members = list_zip_members(args.list_members)
        atomic_write_text(args.out, "".join(f"{m}\n" for m in members))
        return

    if not args.gwas_location or not args.out_dir or not args.metadata:
        parser.error("--gwas-location, --out-dir and --metadata are required unless --list-members is used")
    write_json(args.metadata, download_archive_member(args.gwas_location, args.out_dir))


if __name__ == "__main__":
    main()
