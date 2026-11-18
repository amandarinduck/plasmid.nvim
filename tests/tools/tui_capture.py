"""Drive real Neovim in a pseudo-terminal and capture rendered screens.

End-to-end check of what unit tests cannot see: the decoration provider's
highlights, the winbar, floats, and key handling, as a terminal displays them.

    python tests/tools/tui_capture.py FILE OUTPREFIX [--cols 130 --rows 46]

Steps are (name, keys) pairs; after each step the screen is saved as
OUTPREFIX_<name>.txt and OUTPREFIX_<name>.png.
"""
import argparse
import os
import pathlib
import pty
import select
import subprocess
import time

import pyte
from PIL import Image, ImageDraw, ImageFont

ROOT = pathlib.Path(__file__).resolve().parents[2]
FONT = "/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf"
FONT_BOLD = "/usr/share/fonts/truetype/dejavu/DejaVuSansMono-Bold.ttf"

NAMED = {
    "black": "000000", "red": "cd3131", "green": "0dbc79", "brown": "e5e510", "blue": "2472c8",
    "magenta": "bc3fbc", "cyan": "11a8cd", "white": "e5e5e5", "default": None,
}


def color(c, fallback):
    if c in (None, "default"):
        return fallback
    if c in NAMED:
        return "#" + NAMED[c]
    if len(c) == 6:
        return "#" + c
    return fallback


def render_png(screen, path, cell_w=9, cell_h=18, size=15):
    font = ImageFont.truetype(FONT, size)
    bold = ImageFont.truetype(FONT_BOLD, size) if os.path.exists(FONT_BOLD) else font
    img = Image.new("RGB", (screen.columns * cell_w, screen.lines * cell_h), "#14161b")
    d = ImageDraw.Draw(img)
    for y in range(screen.lines):
        row = screen.buffer[y]
        for x in range(screen.columns):
            ch = row[x]
            fg = color(ch.fg, "#e0e2ea")
            bg = color(ch.bg, "#14161b")
            if ch.reverse:
                fg, bg = bg, fg
            x0, y0 = x * cell_w, y * cell_h
            if bg != "#14161b":
                d.rectangle([x0, y0, x0 + cell_w - 1, y0 + cell_h - 1], fill=bg)
            cp = ord(ch.data[0]) if ch.data else 32
            if 0x2800 < cp <= 0x28FF:
                # most monospace fonts lack braille; draw the dots from the bit pattern
                bits = cp - 0x2800
                for bit, (dx, dy) in enumerate([(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (0, 3), (1, 3)]):
                    if bits & (1 << bit):
                        px = x0 + 2 + dx * (cell_w - 4) // 1.6
                        py = y0 + 2 + dy * (cell_h - 4) / 3.4
                        d.ellipse([px, py, px + 2.4, py + 2.4], fill=fg)
            elif ch.data.strip():
                d.text((x0, y0 + 1), ch.data, font=bold if ch.bold else font, fill=fg)
    # cursor
    cx, cy = screen.cursor.x, screen.cursor.y
    d.rectangle([cx * cell_w, cy * cell_h, cx * cell_w + cell_w - 1, cy * cell_h + cell_h - 1],
                outline="#ffcc00")
    img.save(path)


def run(file, outprefix, steps, cols, rows, pre_cmds=()):
    screen = pyte.Screen(cols, rows)
    stream = pyte.ByteStream(screen)
    pid, fd = pty.fork()
    if pid == 0:
        os.environ.update(TERM="xterm-256color", COLORTERM="truecolor", LANG="C.UTF-8")
        os.execvp("nvim", [
            "nvim", "-u", "NONE", "-i", "NONE", "-n",
            "--cmd", f"set rtp^={ROOT}", "--cmd", "set termguicolors noshowmode laststatus=2",
            *[x for c in pre_cmds for x in ("--cmd", c)],
            "--cmd", "runtime plugin/plasmid.lua",
            file,
        ])
    os.set_blocking(fd, False)
    import fcntl, struct, termios
    fcntl.ioctl(fd, termios.TIOCSWINSZ, struct.pack("HHHH", rows, cols, 0, 0))

    def pump(seconds):
        end = time.time() + seconds
        while time.time() < end:
            r, _, _ = select.select([fd], [], [], 0.05)
            if r:
                try:
                    data = os.read(fd, 65536)
                except OSError:
                    return
                stream.feed(data)
                # answer terminal queries nvim sends at startup
                if b"\x1b[c" in data or b"\x1b[0c" in data:
                    os.write(fd, b"\x1b[?62;22c")

    pump(2.5)
    saved = []
    for name, keys in steps:
        for chunk in keys:
            if chunk == "<wait>":
                pump(3.0)  # let asynchronous work (network requests) finish
                continue
            os.write(fd, chunk.encode())
            pump(0.4)
        pump(1.0)
        if name.startswith("_"):
            continue  # setup step, no screenshot
        txt = "\n".join(screen.display)
        pathlib.Path(f"{outprefix}_{name}.txt").write_text(txt)
        render_png(screen, f"{outprefix}_{name}.png")
        saved.append(name)
    os.write(fd, b"\x1b:qa!\r")
    pump(0.5)
    return saved


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("file")
    ap.add_argument("outprefix")
    ap.add_argument("--cols", type=int, default=130)
    ap.add_argument("--rows", type=int, default=46)
    ap.add_argument("--preset", default="tour", choices=["tour", "map", "linear"])
    ap.add_argument("--script", help="file of 'name|keys' lines; \\r = Enter, \\e = Esc, \\t = Tab, <wait> = pause 3 s")
    ap.add_argument("--pre", action="append", default=[], help="extra --cmd before the plugin loads")
    a = ap.parse_args()
    steps = {
        "tour": [
            ("linear", [":DnaGoto 100\r", "j", "j"]),
            ("info", ["K"]),
            ("circular", ["\x1b", "\t"]),
        ],
        "map": [("circular", [":DnaView circular\r"]), ("linear", ["\t", "]f", "]f"])],
        "linear": [("linear", [":set cmdheight=1\r", "]f"])],
    }[a.preset]
    if a.script:
        steps = []
        for line in pathlib.Path(a.script).read_text().splitlines():
            if line.strip() and not line.startswith("#"):
                name, keys = line.split("|", 1)
                keys = keys.replace("\\r", "\r").replace("\\e", "\x1b").replace("\\t", "\t").replace("\\u", "\x15")
                parts = []
                for k in keys.split("\x00"):
                    for i, piece in enumerate(k.split("<wait>")):
                        if i > 0:
                            parts.append("<wait>")
                        if piece:
                            parts.append(piece)
                steps.append((name, parts))
    print(run(a.file, a.outprefix, steps, a.cols, a.rows, a.pre))
