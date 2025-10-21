#!/usr/bin/env python3
"""ARDaP HTML report parser

Simple script to parse per-sample ARDaP HTML reports and extract a few
fields into a TSV. Intended to live in `bin/` so it can be used by other
scripts or run manually.

Usage:
  python bin/parse_ardap.py <input_folder> <output_tsv>

Dependencies: beautifulsoup4
"""

import os
import sys
import csv
import re
from bs4 import BeautifulSoup


def parse_report(html_path):
    """
    Parse a single NGS HTML report and return a dict with:
      - SampleID
      - SummaryLine1
      - SummaryLine2
      - ResistancePredict (the selected option)
      - AntimicrobialDeterminantDetails
    """
    with open(html_path, "r", encoding="utf-8") as f:
        soup = BeautifulSoup(f, "html.parser")

    # Field placeholders
    parsed = {
        "SampleID": "",
        "SummaryLine1": "",
        "SummaryLine2": "",
        "ResistancePredict": "",
        "AntimicrobialDeterminantDetails": ""
    }

    # 1) SAMPLE ID: prefer an explicit "Sample ID" label in the report
    #    (robust, case-insensitive search in <td> or <th>); fall back to
    #    the HTML <title>, then to the filename.
    sample_id = None

    sample_id_label = soup.find(lambda tag: tag.name in ["td", "th"] and 'sample id' in tag.get_text(strip=True).lower())
    if sample_id_label:
        # Try to find the value in the next sibling <td>
        value_cell = sample_id_label.find_next_sibling(lambda t: t.name == 'td')
        if not value_cell:
            parent_tr = sample_id_label.find_parent('tr')
            if parent_tr:
                cells = parent_tr.find_all(['td', 'th'])
                for i,cell in enumerate(cells):
                    if cell is sample_id_label and i+1 < len(cells):
                        value_cell = cells[i+1]
                        break
        if value_cell:
            txt = value_cell.get_text(strip=True)
            if txt:
                sample_id = txt.rstrip(':').strip()

    # fallback: use <title> if present
    if not sample_id and soup.title and soup.title.get_text(strip=True):
        sample_id = soup.title.get_text(strip=True)

    # final fallback: basename of the HTML file
    if not sample_id:
        sample_id = os.path.splitext(os.path.basename(html_path))[0]

    parsed["SampleID"] = sample_id

    # If the detected sample_id looks like a site title (common fixed banner),
    # try to find a sample-like token elsewhere in the page and prefer that.
    text = soup.get_text(" ", strip=True)
    # patterns: SAMPLE_1_S2_L001 style, SRRxxxx, ERRxxxx
    patterns = [r"\b[A-Za-z0-9\-]+_S\d+_L\d+\b", r"\bSRR\d+\b", r"\bERR\d+\b"]
    for pat in patterns:
        m = re.search(pat, text)
        if m:
            candidate = m.group(0)
            if candidate and candidate != parsed["SampleID"]:
                parsed["SampleID"] = candidate
                break

    # 2) SUMMARY LINES: locate the table whose <thead> or <th> says "Summary"
    summary_table = None
    all_detail_tables = soup.find_all("table", class_="detail_table")
    for tbl in all_detail_tables:
        thead = tbl.find("thead")
        if thead and "Summary" in thead.get_text(strip=True):
            summary_table = tbl
            break

    if summary_table:
        tbody = summary_table.find("tbody")
        if tbody:
            rows = tbody.find_all("tr", recursive=False)
            if len(rows) >= 1:
                parsed["SummaryLine1"] = rows[0].get_text(strip=True)
            if len(rows) >= 2:
                parsed["SummaryLine2"] = rows[1].get_text(strip=True)

    # 3) RESISTANCE PREDICTION: prefer explicit checked <input> elements (type=checkbox)
    #    and their associated <label>. If not present, look for known option
    #    strings (e.g. "No drug resistance predicted") that are adjacent to a
    #    checkbox glyph. Avoid returning the entire page text by constraining
    #    matches to short labelled options.
    def find_checked_option(soup, page_text):
        # 1) Look for an actual <input type="checkbox" checked> element
        inp = soup.find(lambda t: t.name == 'input' and t.get('type') == 'checkbox' and (t.has_attr('checked') or t.get('aria-checked') == 'true'))
        if inp:
            # Try label association via for/id
            if inp.has_attr('id'):
                lbl = soup.find('label', attrs={'for': inp['id']})
                if lbl:
                    return lbl.get_text(strip=True)
            # sibling or parent label
            sib = inp.find_next_sibling('label')
            if sib:
                return sib.get_text(strip=True)
            parent_lbl = inp.find_parent('label')
            if parent_lbl:
                return parent_lbl.get_text(strip=True)
            # as a last resort return nearby short text (limit to 200 chars)
            container = inp.find_parent()
            if container:
                txt = container.get_text(' ', strip=True)
                return txt[:200].strip()

        # 2) Known option strings and their checked glyphs
        box_chars = r'[\u2610\u2611\u2713\u2714☑☒☐✔✖]'
        options = [
            "No drug resistance predicted",
            "Mono-resistance predicted",
            "Multi-drug resistance predicted",
            "Extensive drug resistance predicted",
        ]
        # Search for option preceded or followed by a checkbox glyph within a small window
        for opt in options:
            m = re.search(rf'({box_chars}).{{0,50}}{re.escape(opt)}', page_text, flags=re.I)
            if m:
                return opt
            m2 = re.search(rf'{re.escape(opt)}.{{0,50}}({box_chars})', page_text, flags=re.I)
            if m2:
                return opt

        # 3) If none explicitly marked, return the first option that appears on the page
        for opt in options:
            if re.search(re.escape(opt), page_text, flags=re.I):
                return opt

        return ""

    parsed['ResistancePredict'] = find_checked_option(soup, text)

    # 4) ANTIMICROBIAL DETERMINANT DETAILS
    antidet_table = None
    for tbl in all_detail_tables:
        thead = tbl.find("thead")
        if thead and "Antimicrobial determinant details" in thead.get_text(strip=True):
            antidet_table = tbl
            break

    if antidet_table:
        tbody = antidet_table.find("tbody")
        if tbody:
            raw = tbody.get_text(" | ", strip=True)
            txt = re.sub(r"\s+", " ", raw)
            txt = re.sub(r"\s*\|\s*", " | ", txt)
            txt = txt.strip()
            parsed["AntimicrobialDeterminantDetails"] = txt

    return parsed


def parse_folder_to_tsv(input_folder, output_tsv):
    """
    Loop over all .html files in `input_folder`, parse them, and write TSV output.
    """
    rows = []
    for fname in sorted(os.listdir(input_folder)):
        if fname.lower().endswith(".html"):
            path = os.path.join(input_folder, fname)
            row_data = parse_report(path)
            rows.append(row_data)

    fieldnames = [
        "SampleID",
        "SummaryLine1",
        "SummaryLine2",
        "ResistancePredict",
        "AntimicrobialDeterminantDetails",
    ]

    with open(output_tsv, "w", newline="", encoding="utf-8") as out_f:
        writer = csv.DictWriter(out_f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main():
    if len(sys.argv) < 3:
        print("Usage: python bin/parse_ardap.py <input_folder> <output_tsv>")
        sys.exit(1)

    input_folder = sys.argv[1]
    output_tsv = sys.argv[2]

    parse_folder_to_tsv(input_folder, output_tsv)


if __name__ == "__main__":
    main()
