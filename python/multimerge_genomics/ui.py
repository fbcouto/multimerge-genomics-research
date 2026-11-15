"""MMGEN's web interface: the marimo app that comes with the package.

    pip install "multimerge-genomics[ui]"
    mmgen-ui                         # opens it in the browser of this computer
    mmgen-ui --folder D:/my/data     # the screens list (and write) the files of that folder
    python -m multimerge_genomics.ui --headless      # on a server or on Colab

The screens list the files of the folder it is started in (or of --folder)
and write their results there. By default it answers only on this computer
(127.0.0.1) and needs no password; with --host 0.0.0.0 other machines can
reach it, and then a password is required (one is made up and printed when
none is given).
"""
from __future__ import annotations

import argparse
import importlib.util
import os
import secrets
import subprocess
import sys

APP = os.path.join(os.path.dirname(os.path.abspath(__file__)), "ui_main.py")
_LOCAL = ("127.0.0.1", "localhost", "::1")


def app_path() -> str:
    """The path of the interface's marimo file."""
    return APP


def command(port: int = 2718, host: str = "127.0.0.1", headless: bool = False, password: str | None = None) -> list:
    """The `marimo run` command line that serves the interface."""
    cmd = [sys.executable, "-m", "marimo", "run", APP, "--port", str(port), "--host", host]
    if headless:
        cmd.append("--headless")
    cmd += ["--token", "--token-password", password] if password else ["--no-token"]
    return cmd


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="mmgen-ui", description="MMGEN's web interface (29 screens, with ready-made examples).")
    ap.add_argument("--folder", default=os.getcwd(),
                    help="the folder whose files the screens list and where results are written (default: the current one)")
    ap.add_argument("--port", type=int, default=2718)
    ap.add_argument("--host", default="127.0.0.1",
                    help="127.0.0.1: this computer only (the default); 0.0.0.0: other machines too (a password is then required)")
    ap.add_argument("--headless", action="store_true", help="do not open a browser (on a server, or on Colab)")
    ap.add_argument("--password", help="require this password (made up and printed when --host is not this computer)")
    args = ap.parse_args(argv)

    if importlib.util.find_spec("marimo") is None:
        sys.exit('The web interface needs marimo and a few other packages: pip install "multimerge-genomics[ui]"')
    folder = os.path.abspath(args.folder)
    if not os.path.isdir(folder):
        sys.exit(f"--folder {args.folder!r} is not a folder")
    password = args.password
    if args.host not in _LOCAL and not password:
        password = secrets.token_urlsafe(16)
    shown = "127.0.0.1" if args.host in ("0.0.0.0", "::") else args.host
    url = f"http://{shown}:{args.port}/" + (f"?access_token={password}" if password else "")
    print(f"MMGEN web interface: {url}")
    print(f"Files: {folder}")
    if password and not args.password:
        print(f"Password (made up for this run): {password}")
    print("Ctrl+C stops it.", flush=True)
    try:
        return subprocess.call(command(args.port, args.host, args.headless, password), cwd=folder)
    except KeyboardInterrupt:
        return 0


if __name__ == "__main__":
    sys.exit(main())
