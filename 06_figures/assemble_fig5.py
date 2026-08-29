#!/usr/bin/env python3
"""Place the four Figure 5 panels on one page.

Panels are written at 252 x 226 pt by their own scripts, so they drop in without
rescaling and text keeps its intended size:

    a  make_fig5a_track.R      Fig5a_track.pdf
    b  make_fig5b_recomb.R     Fig5b_recomb.pdf
    c  assemble_fig4_fig5.R    Fig4_panel_gynae.pdf
    d  make_fig5d_maletract.R  Fig5d_maletract.pdf

Usage: python assemble_fig5.py <panel_dir> <output.pdf>
"""
import sys, os
import fitz  # PyMuPDF

PW, PH = 252.0, 226.0          # panel size, points
W = 517.06                     # double-column width
GX = W - 2 * PW                # horizontal gap
GY = 16.0                      # vertical gap
PANELS = [("a", "Fig5a_track.pdf"), ("b", "Fig5b_recomb.pdf"),
          ("c", "Fig4_panel_gynae.pdf"), ("d", "Fig5d_maletract.pdf")]

def main(panel_dir, out_path):
    out = fitz.open()
    page = out.new_page(width=W, height=2 * PH + GY)
    for i, (label, fname) in enumerate(PANELS):
        src = fitz.open(os.path.join(panel_dir, fname))
        col, row = i % 2, i // 2
        x0, y0 = col * (PW + GX), row * (PH + GY)
        page.show_pdf_page(fitz.Rect(x0, y0, x0 + PW, y0 + PH), src, 0)
        src.close()
        page.insert_text(fitz.Point(x0 + 2, y0 + 10), label,
                         fontname="hebo", fontsize=11)
    out.save(out_path)
    out.close()
    print("wrote", out_path)

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "figures_final",
         sys.argv[2] if len(sys.argv) > 2 else "Fig5_4panel_final.pdf")
