"""Generate the SpatialScribe favicon: a dark-field 'cell neighbourhood' mark in brand colors.

A dark rounded tile with a tight cluster of fluorescence-colored cells (one larger cyan cell
with a nucleus, plus magenta / green / amber / rose satellites) - reads as spatial cell types on
a dark-field view even at 16px. Regenerate with:  python scripts/make_favicon.py
"""
from __future__ import annotations

import pathlib

from PIL import Image, ImageDraw

SIZE = 256
BG = (10, 14, 22, 255)  # #0a0e16
OUT = pathlib.Path(__file__).resolve().parents[1] / "src" / "spatialscribe" / "app" / "favicon.png"

# (cx, cy, r, rgb) - a compact neighbourhood, cyan cell dominant.
CELLS = [
    (128, 122, 46, (34, 211, 238)),   # cyan (large, with nucleus)
    (80, 182, 26, (232, 121, 249)),   # magenta
    (182, 172, 24, (52, 211, 153)),   # green
    (178, 86, 20, (251, 191, 36)),    # amber
    (86, 78, 18, (244, 113, 133)),    # rose
]


def main() -> None:
    img = Image.new("RGBA", (SIZE, SIZE), (0, 0, 0, 0))
    d = ImageDraw.Draw(img)
    d.rounded_rectangle([0, 0, SIZE - 1, SIZE - 1], radius=58, fill=BG)
    d.rounded_rectangle([0, 0, SIZE - 1, SIZE - 1], radius=58, outline=(29, 37, 49, 255), width=3)
    for cx, cy, r, rgb in CELLS:
        d.ellipse([cx - r, cy - r, cx + r, cy + r], fill=rgb + (255,))
    # nucleus of the dominant cell
    d.ellipse([128 - 17, 122 - 17, 128 + 17, 122 + 17], fill=(9, 26, 33, 255))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    img.save(OUT)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
