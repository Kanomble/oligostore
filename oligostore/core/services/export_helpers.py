from openpyxl.styles import Font

PRIMER_HEADERS = [
    "Name",
    "Sequence",
    "5' Overhang",
    "Restriction Sites",
    "Length",
    "GC Content",
    "Temperature",
    "Hairpin",
    "Self Dimer",
    "Creator",
    "Created",
]


def style_header_row(sheet):
    header_font = Font(bold=True)
    for cell in sheet[1]:
        cell.font = header_font


def apply_column_widths(sheet, min_width=12, max_width=50, padding=2):
    for column_cells in sheet.columns:
        max_length = 0
        for cell in column_cells:
            cell_value = cell.value
            if cell_value is None:
                continue
            max_length = max(max_length, len(str(cell_value)))
        adjusted_width = min(max_length + padding, max_width)
        column_letter = column_cells[0].column_letter
        sheet.column_dimensions[column_letter].width = max(adjusted_width, min_width)


def build_primer_worksheet(workbook, primers, title="Primers"):
    sheet = workbook.active
    sheet.title = title

    sheet.append(PRIMER_HEADERS)
    style_header_row(sheet)

    for primer in primers:
        sheet.append(
            [
                primer.primer_name,
                primer.sequence,
                primer.overhang_sequence,
                primer.restriction_site_summary,
                primer.length,
                primer.gc_content,
                primer.tm,
                primer.hairpin_dg,
                primer.self_dimer_dg,
                str(primer.creator),
                primer.created_at.strftime("%Y-%m-%d %H:%M"),
            ]
        )

    apply_column_widths(sheet)
    return sheet


def build_primerpair_worksheet(workbook, primer_pairs, title="Primer Pairs"):
    sheet = workbook.active
    sheet.title = title

    headers = ["Pair Name"]
    headers.extend([f"Forward {header}" for header in PRIMER_HEADERS])
    headers.extend([f"Reverse {header}" for header in PRIMER_HEADERS])
    sheet.append(headers)
    style_header_row(sheet)

    for pair in primer_pairs:
        forward = pair.forward_primer
        reverse = pair.reverse_primer
        sheet.append(
            [
                pair.name,
                forward.primer_name,
                forward.sequence,
                forward.overhang_sequence,
                forward.restriction_site_summary,
                forward.length,
                forward.gc_content,
                forward.tm,
                forward.hairpin_dg,
                forward.self_dimer_dg,
                str(forward.creator),
                forward.created_at.strftime("%Y-%m-%d %H:%M"),
                reverse.primer_name,
                reverse.sequence,
                reverse.overhang_sequence,
                reverse.restriction_site_summary,
                reverse.length,
                reverse.gc_content,
                reverse.tm,
                reverse.hairpin_dg,
                reverse.self_dimer_dg,
                str(reverse.creator),
                reverse.created_at.strftime("%Y-%m-%d %H:%M"),
            ]
        )

    apply_column_widths(sheet)
    return sheet