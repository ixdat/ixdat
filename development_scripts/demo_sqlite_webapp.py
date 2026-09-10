"""Save ixdat data and view the SQLite tables in a browser."""

import html
import sqlite3
import tempfile
import webbrowser
from http.server import BaseHTTPRequestHandler, HTTPServer
from pathlib import Path
from urllib.parse import parse_qs

from ixdat import Measurement, Spectrum
from ixdat.db import DB, change_database


DATA = Path(__file__).resolve().parents[1] / "test_data"
DATABASE = Path(tempfile.gettempdir()) / "ixdat_web_demo.sqlite"


def save_data():
    if DATABASE.exists():
        DATABASE.unlink()
    change_database("sqlite", db_path=DATABASE)
    Measurement.read(DATA / "biologic/Pt_poly_cv_CUT.mpt", reader="biologic").save()
    Measurement.read(DATA / "biologic/Pt_poly_cv.mpt", reader="biologic").save()
    Spectrum.read(DATA / "xrd/twotheta_2th.xy", reader="xrdxy").save()
    DB.backend.close()


def make_page(selected):
    db = sqlite3.connect(DATABASE.resolve().as_uri() + "?mode=ro", uri=True)
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


class Handler(BaseHTTPRequestHandler):
    def do_GET(self):
        query = parse_qs(self.path.partition("?")[2])
        page = make_page(query.get("table", [None])[0]).encode()
        self.send_response(200)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.end_headers()
        self.wfile.write(page)


def main():
    save_data()
    server = HTTPServer(("127.0.0.1", 0), Handler)
    url = "http://127.0.0.1:{}".format(server.server_port)
    print("Open {} (Ctrl+C stops the server)".format(url))
    webbrowser.open(url)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        server.server_close()


if __name__ == "__main__":
    main()
