"""Get the demo data and view the tables of an SQLite database in a browser."""

import configparser
import html
import os
import sqlite3
import webbrowser
from http.server import BaseHTTPRequestHandler, HTTPServer
from pathlib import Path
from urllib.parse import parse_qs

from ..tools import get_default_cache_dir

REPO_DIR = Path(__file__).parent.parent.parent.parent
DEMO_DATA_ENV_VAR = "IXDAT_DEMO_DATA_DIR"
DEMO_DATA_CONFIG = Path(__file__).parent / "demo_data.ini"


def read_archive_config():
    """Return the [archive] section of demo_data.ini as a dict."""
    config = configparser.ConfigParser()
    config.read(DEMO_DATA_CONFIG)
    return dict(config["archive"])


def get_downloaded_data_dir(cache_dir=None):
    """Return where download_demo_data() puts the demo data, as a Path.

    Args:
        cache_dir (str or Path): The folder to download into. Defaults to the
            ``demo_data`` folder in ixdat's cache folder.
    """
    archive = read_archive_config()
    cache_dir = Path(cache_dir or get_default_cache_dir("ixdat") / "demo_data")
    return cache_dir / f"ixdat-demo-data-{archive['version']}" / archive["folder"]


def download_demo_data(cache_dir=None):
    """Download and unzip the demo data archive given in demo_data.ini.

    This requires the `pooch` package. Pooch checks the archive against its sha256
    checksum and skips the download if the archive is already there.

    Args:
        cache_dir (str or Path): The folder to download into. Defaults to the
            ``demo_data`` folder in ixdat's cache folder.

    Returns:
        Path: the folder with the demo data
    """
    import pooch  # not a requirement of ixdat, so we import it here.

    archive = read_archive_config()
    cache_dir = Path(cache_dir or get_default_cache_dir("ixdat") / "demo_data")
    name = f"ixdat-demo-data-{archive['version']}"
    pooch.retrieve(
        url=archive["url"],
        known_hash="sha256:" + archive["sha256"],
        fname=name + ".zip",
        path=cache_dir,
        processor=pooch.Unzip(extract_dir=name),
        progressbar=False,
    )
    return get_downloaded_data_dir(cache_dir)


def get_demo_data_dir(data_dir=None):
    """Return the folder with the demo data used by the demos, as a Path.

    Args:
        data_dir (str or Path): The folder to use. By default, the first of these
            which is set or exists: the environment variable IXDAT_DEMO_DATA_DIR,
            ``demo_data/`` in the root of the ixdat repository, and the folder made
            by :func:`download_demo_data`.

    Raises:
        FileNotFoundError: if the folder does not exist.
    """
    if not data_dir:
        data_dir = os.environ.get(DEMO_DATA_ENV_VAR)
    if not data_dir and (REPO_DIR / "demo_data").is_dir():
        data_dir = REPO_DIR / "demo_data"
    data_dir = Path(data_dir or get_downloaded_data_dir())
    if not data_dir.is_dir():
        raise FileNotFoundError(
            f"No demo data found at {data_dir}. Download it with "
            "`python -m ixdat.demos.download` (requires `pip install pooch`), or set "
            f"the environment variable {DEMO_DATA_ENV_VAR} to a demo data folder."
        )
    return data_dir


def get_test_data_dir(data_dir=None):
    """Return `data_dir`, by default the ``test_data/`` folder of the repo, as a Path."""
    return Path(data_dir or REPO_DIR / "test_data")


def make_page(selected, sqlite_file):
    db = sqlite3.connect(Path(sqlite_file).resolve().as_uri() + "?mode=ro", uri=True)
    tables = [
        row[0]
        for row in db.execute(
            "SELECT name FROM sqlite_master "
            "WHERE type='table' AND name NOT LIKE 'sqlite_%' ORDER BY name"
        )
    ]
    table = selected if selected in tables else tables[0]
    cursor = db.execute('SELECT * FROM "{}"'.format(table))
    columns = [column[0] for column in cursor.description]
    rows = cursor.fetchall()
    db.close()

    links = " ".join(
        '<a href="/?table={0}">{0}</a>'.format(html.escape(name)) for name in tables
    )
    headings = "".join("<th>{}</th>".format(html.escape(name)) for name in columns)
    body = "".join(
        "<tr>{}</tr>".format("".join("<td>{}</td>".format(show(value)) for value in row))
        for row in rows
    )
    return (
        "<!doctype html><title>ixdat SQLite tables</title>"
        "<style>body{{font:14px sans-serif}}a{{margin-right:1em}}"
        "table{{border-collapse:collapse}}th,td{{border:1px solid;padding:4px}}"
        "th{{background:#eee}}</style><h1>ixdat SQLite tables</h1>"
        "<nav>{}</nav><h2>{}</h2><table><tr>{}</tr>{}</table>"
    ).format(links, html.escape(table), headings, body)


def show(value):
    return (
        "{} bytes".format(len(value))
        if isinstance(value, bytes)
        else html.escape(str(value))
    )


class TableHandler(BaseHTTPRequestHandler):
    """Serve one page per table of the SQLite file given by ``sqlite_file``."""

    sqlite_file = None

    def do_GET(self):
        query = parse_qs(self.path.partition("?")[2])
        page = make_page(
            query.get("table", [None])[0], sqlite_file=self.sqlite_file
        ).encode()
        self.send_response(200)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.end_headers()
        self.wfile.write(page)


def view_tables(sqlite_file):
    """Open the tables of `sqlite_file` in a browser. Ctrl+C stops the server."""
    # The server makes a new handler for each request, so the file goes on a
    # handler class made for this call.
    handler = type("TableHandler", (TableHandler,), {"sqlite_file": sqlite_file})
    server = HTTPServer(("127.0.0.1", 0), handler)
    url = "http://127.0.0.1:{}".format(server.server_port)
    print("Open {} (Ctrl+C stops the server)".format(url))
    webbrowser.open(url)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        server.server_close()
