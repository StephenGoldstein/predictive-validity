"""download(): segmented range requests, per-segment resume, fallback without ranges, legacy partial files."""
import os
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

from pipeline.gwas_catalog import download

DATA = bytes(range(256)) * 4096   # 1 MiB


def serve(ranges=True, honour_get=True):
    """A range-capable file server on a free port; (server, url, [(method, Range header)...])."""
    seen = []

    class H(BaseHTTPRequestHandler):
        def log_message(self, *a):
            pass

        def _bounds(self):
            h = self.headers.get("Range")
            seen.append((self.command, h))
            if h and ranges and (honour_get or self.command == "HEAD"):
                a, b = h.split("=")[1].split("-")
                return int(a), (int(b) + 1 if b else len(DATA)), True
            return 0, len(DATA), False

        def _head(self, start, end, partial):
            self.send_response(206 if partial else 200)
            self.send_header("Content-Length", str(end - start))
            if ranges:
                self.send_header("Accept-Ranges", "bytes")
            if partial:
                self.send_header("Content-Range", f"bytes {start}-{end - 1}/{len(DATA)}")
            self.end_headers()

        def do_HEAD(self):
            self._head(*self._bounds())

        def do_GET(self):
            start, end, partial = self._bounds()
            self._head(start, end, partial)
            self.wfile.write(DATA[start:end])

    srv = ThreadingHTTPServer(("127.0.0.1", 0), H)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    return srv, f"http://127.0.0.1:{srv.server_address[1]}/f.gz", seen


def test_segmented_download_assembles_the_file(tmp_path):
    srv, url, seen = serve()
    dest = str(tmp_path / "f.gz")
    assert download(url, dest, streams=4) == dest
    assert open(dest, "rb").read() == DATA
    gets = [r for m, r in seen if m == "GET"]
    assert len(gets) == 4 and all(gets)
    assert not any(p.startswith("f.gz.part") for p in os.listdir(tmp_path))
    srv.shutdown()


def test_segmented_download_resumes_a_segment(tmp_path):
    srv, url, seen = serve()
    dest = str(tmp_path / "f.gz")
    n = len(DATA) // 4
    with open(dest + ".part1", "wb") as f:
        f.write(DATA[n:n + 1000])
    download(url, dest, streams=4)
    assert open(dest, "rb").read() == DATA
    assert ("GET", f"bytes={n + 1000}-{2 * n - 1}") in seen
    srv.shutdown()


def test_download_falls_back_to_one_stream_without_ranges(tmp_path):
    srv, url, seen = serve(ranges=False)
    dest = str(tmp_path / "f.gz")
    download(url, dest, streams=4)
    assert open(dest, "rb").read() == DATA
    assert [m for m, _ in seen].count("GET") == 1
    srv.shutdown()


def test_legacy_partial_file_finishes_on_one_stream(tmp_path):
    srv, url, seen = serve()
    dest = str(tmp_path / "f.gz")
    with open(dest + ".part", "wb") as f:
        f.write(DATA[:5000])
    download(url, dest, streams=4)
    assert open(dest, "rb").read() == DATA
    assert [r for m, r in seen if m == "GET"] == ["bytes=5000-"]
    srv.shutdown()


def test_segment_whose_range_is_ignored_raises(tmp_path):
    import pytest
    import requests
    srv, url, seen = serve(honour_get=False)
    dest = str(tmp_path / "f.gz")
    with pytest.raises(requests.HTTPError):
        download(url, dest, attempts=1, streams=4)
    assert not os.path.exists(dest)
    srv.shutdown()
