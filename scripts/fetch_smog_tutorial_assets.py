#!/usr/bin/env python3
"""Fetch public SMOG tutorial pages and linked model-generation assets.

The crawler starts at https://smog-server.org/tutorials/, follows tutorial-local
HTML/directory pages, and downloads small linked input/template assets.  It is
deliberately conservative: generated manifests are committed, but downloaded
public data lives under validation/tutorials/assets/ and is ignored by git.
"""

from __future__ import annotations

import argparse
import json
import mimetypes
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import deque
from dataclasses import asdict, dataclass
from html.parser import HTMLParser
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_START = "https://smog-server.org/tutorials/"
DEFAULT_ASSET_ROOT = ROOT / "validation" / "tutorials" / "assets"
USER_AGENT = "SMOG3 tutorial asset fetcher/0.1 (+https://smog-server.org/tutorials/)"

DOWNLOAD_EXTENSIONS = {
    ".pdb",
    ".xml",
    ".bif",
    ".sif",
    ".b",
    ".nb",
    ".map",
    ".contacts",
    ".top",
    ".gro",
    ".g96",
    ".ndx",
    ".tar.gz",
    ".tgz",
    ".zip",
    ".dat",
    ".txt",
    ".py",
    ".xvg",
    ".gz",
}
SKIP_QUERY_PREFIXES = ("C=",)


@dataclass
class FetchRecord:
    """One crawled URL and the local file, status, and source page for it."""

    url: str
    source_page: str | None
    local_path: str | None
    http_status: int | None
    content_type: str | None
    file_type: str
    size: int | None
    status: str
    note: str = ""


class LinkParser(HTMLParser):
    """Small HTML link extractor for href/src attributes."""

    def __init__(self) -> None:
        super().__init__()
        self.links: list[str] = []

    def handle_starttag(self, _tag: str, attrs: list[tuple[str, str | None]]) -> None:
        for name, value in attrs:
            if name.lower() in {"href", "src"} and value:
                self.links.append(value)


def _request(url: str) -> urllib.request.Request:
    return urllib.request.Request(url, headers={"User-Agent": USER_AGENT})


def _clean_url(url: str) -> str:
    parsed = urllib.parse.urlsplit(url)
    return urllib.parse.urlunsplit((parsed.scheme, parsed.netloc, parsed.path, "", ""))


def _is_sort_query(url: str) -> bool:
    query = urllib.parse.urlsplit(url).query
    return any(query.startswith(prefix) for prefix in SKIP_QUERY_PREFIXES)


def _under_tutorials(url: str, start_host: str) -> bool:
    parsed = urllib.parse.urlsplit(url)
    return parsed.netloc == start_host and parsed.path.startswith("/tutorials")


def _html_base_url(url: str) -> str:
    parsed = urllib.parse.urlsplit(url)
    if url.endswith("/") or _compound_suffix(parsed.path) in DOWNLOAD_EXTENSIONS:
        return url
    return url + "/"


def _compound_suffix(path: str) -> str:
    lowered = path.lower()
    for suffix in (".tar.gz",):
        if lowered.endswith(suffix):
            return suffix
    return Path(urllib.parse.urlsplit(path).path).suffix.lower()


def _safe_local_path(root: Path, url: str, directory: str, default_name: str) -> Path:
    parsed = urllib.parse.urlsplit(url)
    raw_parts = [part for part in parsed.path.split("/") if part]
    parts = [urllib.parse.unquote(part) for part in raw_parts]
    if not parts or url.endswith("/"):
        parts.append(default_name)
    is_page = directory == "pages"
    if is_page and _compound_suffix(parts[-1]) not in DOWNLOAD_EXTENSIONS:
        parts[-1] = f"{parts[-1]}.html"
    if is_page and url.endswith("/"):
        parts[-1] = "index.html"
    safe_parts = [re.sub(r"[^A-Za-z0-9._+@=-]+", "_", part) for part in parts]
    return root / directory / parsed.netloc / Path(*safe_parts)


def _looks_like_html(content_type: str | None, data: bytes) -> bool:
    ctype = (content_type or "").lower()
    if "text/html" in ctype:
        return True
    return data[:200].lstrip().lower().startswith((b"<!doctype html", b"<html"))


def _file_type(url: str, content_type: str | None, is_html: bool) -> str:
    if is_html:
        return "html"
    suffix = _compound_suffix(url)
    if suffix:
        return suffix.lstrip(".")
    guessed = mimetypes.guess_extension(content_type or "")
    return (guessed or "unknown").lstrip(".")


def _fetch_url(url: str, source_page: str | None, asset_root: Path, max_bytes: int) -> tuple[FetchRecord, str | None, list[str]]:
    try:
        with urllib.request.urlopen(_request(url), timeout=30) as response:
            status = getattr(response, "status", None)
            content_type = response.headers.get("content-type")
            content_length = response.headers.get("content-length")
            expected_size = int(content_length) if content_length and content_length.isdigit() else None
            if expected_size is not None and expected_size > max_bytes:
                return (
                    FetchRecord(url, source_page, None, status, content_type, "unknown", expected_size, "LARGE_FILE_SKIPPED"),
                    None,
                    [],
                )
            data = response.read(max_bytes + 1)
    except urllib.error.HTTPError as exc:
        return FetchRecord(url, source_page, None, exc.code, exc.headers.get("content-type"), "unknown", None, "FAILED", str(exc)), None, []
    except Exception as exc:  # pragma: no cover - exercised by real network runs
        return FetchRecord(url, source_page, None, None, None, "unknown", None, "FAILED", str(exc)), None, []

    if len(data) > max_bytes:
        return FetchRecord(url, source_page, None, status, content_type, "unknown", len(data), "LARGE_FILE_SKIPPED"), None, []

    is_html = _looks_like_html(content_type, data)
    kind = "page" if is_html else "download"
    local_path = _safe_local_path(asset_root, url, "pages" if is_html else "downloads", "index.html")
    local_path.parent.mkdir(parents=True, exist_ok=True)
    local_path.write_bytes(data)

    links: list[str] = []
    if is_html:
        parser = LinkParser()
        parser.feed(data.decode("utf-8", errors="replace"))
        links = parser.links

    record = FetchRecord(
        url=url,
        source_page=source_page,
        local_path=str(local_path.relative_to(ROOT)),
        http_status=status,
        content_type=content_type,
        file_type=_file_type(url, content_type, is_html),
        size=len(data),
        status="DOWNLOADED" if not is_html else "PAGE_SAVED",
    )
    return record, kind, links


def crawl(start_url: str, asset_root: Path, max_mb: int, max_pages: int) -> dict[str, object]:
    """Crawl tutorial pages and download small linked files."""

    start_url = _clean_url(start_url)
    start_host = urllib.parse.urlsplit(start_url).netloc
    max_bytes = max_mb * 1024 * 1024
    queue: deque[tuple[str, str | None]] = deque([(start_url, None)])
    seen: set[str] = set()
    records: list[FetchRecord] = []
    page_count = 0

    while queue and page_count < max_pages:
        url, source_page = queue.popleft()
        url = _clean_url(url)
        if url in seen or not _under_tutorials(url, start_host):
            continue
        seen.add(url)

        record, kind, links = _fetch_url(url, source_page, asset_root, max_bytes)
        records.append(record)
        if kind != "page":
            continue
        page_count += 1

        base_url = _html_base_url(url)
        for raw in links:
            resolved = urllib.parse.urljoin(base_url, raw)
            if _is_sort_query(resolved):
                continue
            clean = _clean_url(resolved)
            if not _under_tutorials(clean, start_host):
                continue
            queue.append((clean, url))

    status_counts: dict[str, int] = {}
    for record in records:
        status_counts[record.status] = status_counts.get(record.status, 0) + 1
    return {
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "start_url": start_url,
        "asset_root": str(asset_root.relative_to(ROOT)),
        "max_mb": max_mb,
        "max_pages": max_pages,
        "status_counts": status_counts,
        "page_count": sum(1 for record in records if record.status == "PAGE_SAVED"),
        "download_count": sum(1 for record in records if record.status == "DOWNLOADED"),
        "records": [asdict(record) for record in records],
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Download public SMOG tutorial pages and linked input assets.")
    parser.add_argument("--start-url", default=DEFAULT_START)
    parser.add_argument("--asset-root", default=str(DEFAULT_ASSET_ROOT))
    parser.add_argument("--max-mb", type=int, default=100)
    parser.add_argument("--max-pages", type=int, default=250)
    ns = parser.parse_args(argv)

    asset_root = Path(ns.asset_root)
    asset_root.mkdir(parents=True, exist_ok=True)
    manifest = crawl(ns.start_url, asset_root, ns.max_mb, ns.max_pages)
    manifest_path = asset_root / "manifest_raw.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"Indexed {manifest['page_count']} tutorial pages")
    print(f"Downloaded {manifest['download_count']} assets")
    print(f"Wrote {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
