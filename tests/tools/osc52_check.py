#!/usr/bin/env python3
"""End-to-end clipboard check in a real pseudo-terminal (development tool).

Runs Neovim with the plugin in a pty, yanks a region with each strand
orientation, and reads what Neovim writes to the terminal. With no clipboard
provider installed, each plain yank must arrive as an OSC 52 sequence
(ESC ] 52 ; c ; <base64> ESC \\) whose payload equals the expected strand,
computed here independently with Biopython from the same GenBank file.

    python tests/tools/osc52_check.py   # needs biopython; exits 1 on mismatch (GenBank yanks are upper case, as displayed)
"""
import base64
import os
import pathlib
import pty
import re
import select
import shutil
import sys
import time

from Bio import SeqIO

ROOT = pathlib.Path(__file__).resolve().parents[2]
FIXTURE = ROOT / "tests/fixtures/ncbi/J01749.gb"
OSC52 = re.compile(rb"\x1b\]52;c;([A-Za-z0-9+/=]*)(?:\x1b\\|\x07)")


def main():
    if any(shutil.which(t) for t in ("xclip", "xsel", "wl-copy", "pbcopy", "win32yank.exe")):
        print("a clipboard provider is installed; this check covers the OSC 52 fallback only")
        return 0
    seq = str(SeqIO.read(FIXTURE, "genbank").seq).upper()
    plus = seq[9:40]  # 10..40
    comp = plus.translate(str.maketrans("ACGT", "TGCA"))
    expected = [
        ("top", plus),
        ("top35", plus[::-1]),
        ("bottom", comp[::-1]),
        ("bottom35", comp),
    ]
    pid, fd = pty.fork()
    if pid == 0:
        os.environ.update(TERM="xterm-256color", LANG="C.UTF-8")
        os.execvp("nvim", ["nvim", "-u", "NONE", "-i", "NONE", "-n",
                           "--cmd", f"set rtp^={ROOT}", "--cmd", "runtime plugin/plasmid.lua", str(FIXTURE)])
    os.set_blocking(fd, False)
    import fcntl, struct, termios
    fcntl.ioctl(fd, termios.TIOCSWINSZ, struct.pack("HHHH", 40, 130, 0, 0))
    raw = bytearray()

    def pump(seconds):
        end = time.time() + seconds
        while time.time() < end:
            r, _, _ = select.select([fd], [], [], 0.05)
            if r:
                try:
                    data = os.read(fd, 65536)
                except OSError:
                    return
                raw.extend(data)
                if b"\x1b[c" in data or b"\x1b[0c" in data:
                    os.write(fd, b"\x1b[?62;22c")

    pump(2.0)
    os.write(fd, b":DnaSelect 10..40\r")
    pump(0.5)
    got = []
    for key, _ in expected:
        start = len(raw)
        os.write(fd, f":DnaYank {key}\r".encode())
        pump(1.0)
        m = OSC52.findall(bytes(raw[start:]))
        got.append(base64.b64decode(m[-1]).decode() if m else None)
    os.write(fd, b":qa!\r")
    pump(0.5)

    ok = True
    for (key, want), have in zip(expected, got):
        status = "ok" if have == want else "MISMATCH"
        ok = ok and status == "ok"
        print(f"{key:9} {status:8} {have!r}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
