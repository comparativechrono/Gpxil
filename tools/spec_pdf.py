"""Render a Markdown specification as a PDF: pandoc makes the HTML, Chromium prints it.

    python tools/spec_pdf.py SPEC.md gpxil-spec.pdf --title "gpxil 0.1: Genome Pixel Image format"

Needs pandoc and Playwright with Chromium
(pip install playwright; python -m playwright install chromium).
"""

from __future__ import annotations

import argparse
import html
import pathlib
import subprocess
import tempfile

from playwright.sync_api import sync_playwright

CSS = """
@page { size: A4; margin: 22mm 20mm 24mm 20mm; }
html { font: 10.5pt/1.5 "Inter", "DejaVu Sans", sans-serif; color: #15151f; }
body { max-width: none; margin: 0; }
h1 { font-size: 21pt; line-height: 1.2; margin: 0 0 6pt; color: #000; }
h2 { font-size: 14pt; margin: 20pt 0 6pt; border-bottom: 1.5pt solid #e0306f; padding-bottom: 2pt; break-after: avoid; }
h3 { font-size: 11.5pt; margin: 14pt 0 4pt; break-after: avoid; }
p, li { orphans: 3; widows: 3; }
a { color: #1a55c4; text-decoration: none; }
code { font: 9pt "DejaVu Sans Mono", monospace; background: #f1f1f6; padding: 0 2pt; border-radius: 2pt; }
pre { background: #f1f1f6; padding: 7pt 9pt; border-left: 3pt solid #20b3d6; overflow-wrap: anywhere; white-space: pre-wrap; break-inside: avoid; }
pre code { background: none; padding: 0; }
table { border-collapse: collapse; margin: 6pt 0 10pt; width: 100%; font-size: 9.5pt; break-inside: auto; }
th, td { border: 0.75pt solid #c9c9d6; padding: 3pt 6pt; text-align: left; vertical-align: top; }
th { background: #ececf4; }
tr { break-inside: avoid; }
blockquote { margin: 0; padding-left: 10pt; border-left: 3pt solid #ccc; color: #444; }
"""


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("markdown")
    ap.add_argument("pdf")
    ap.add_argument("--title", required=True, help="shown in the page footers")
    a = ap.parse_args(argv)
    with tempfile.TemporaryDirectory() as tmp:
        head = pathlib.Path(tmp, "head.html")
        head.write_text(f"<style>{CSS}</style>")
        page_html = pathlib.Path(tmp, "spec.html")
        subprocess.run(
            ["pandoc", a.markdown, "--from", "gfm", "--to", "html5", "--standalone", "--embed-resources", "--resource-path", str(pathlib.Path(a.markdown).resolve().parent),
             "--metadata", f"pagetitle={a.title}", "--include-in-header", str(head), "--output", str(page_html)],
            check=True,
        )
        footer = (
            '<div style="width:100%;font:7.5pt Inter,DejaVu Sans,sans-serif;color:#666;padding:0 20mm;display:flex;justify-content:space-between">'
            f"<span>{html.escape(a.title)}</span><span><span class=pageNumber></span> / <span class=totalPages></span></span></div>"
        )
        with sync_playwright() as p:
            browser = p.chromium.launch()
            page = browser.new_page()
            page.goto(page_html.resolve().as_uri())
            page.pdf(path=a.pdf, format="A4", print_background=True, display_header_footer=True, header_template="<span></span>", footer_template=footer, margin={"top": "20mm", "bottom": "22mm", "left": "20mm", "right": "20mm"})
            browser.close()
    print(f"wrote {a.pdf}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
