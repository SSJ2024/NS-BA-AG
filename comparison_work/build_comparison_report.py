from pathlib import Path

import pandas as pd
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.platypus import (
    PageBreak,
    Paragraph,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)


BASE = Path(r"C:\Users\a3188798\OneDrive - Adelaide University\Desktop\NS BA AG")
WORK = BASE / "comparison_work"
OUT = BASE / "Bile_Acid_Report_Comparison_Audit.pdf"
PAGE_SIZE = landscape(A4)


def fmt_p(x):
    if pd.isna(x):
        return "NA"
    x = float(x)
    if x < 0.001:
        return "<0.001"
    return f"{x:.3f}"


def fmt_num(x, n=3):
    if pd.isna(x):
        return "NA"
    return f"{float(x):.{n}f}"


def style_table(table, header_bg="#1F4E79", font_size=7, align="CENTER"):
    table.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor(header_bg)),
                ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
                ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                ("FONTSIZE", (0, 0), (-1, -1), font_size),
                ("LEADING", (0, 0), (-1, -1), font_size + 2),
                ("ALIGN", (0, 0), (-1, -1), align),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("GRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#D9E2EC")),
                ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#F7FAFC")]),
                ("LEFTPADDING", (0, 0), (-1, -1), 4),
                ("RIGHTPADDING", (0, 0), (-1, -1), 4),
                ("TOPPADDING", (0, 0), (-1, -1), 4),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
            ]
        )
    )
    return table


def p(text, style):
    return Paragraph(text, style)


def footer(canvas, doc):
    canvas.saveState()
    canvas.setFont("Helvetica", 8)
    canvas.setFillColor(colors.HexColor("#66788A"))
    canvas.drawString(18 * mm, 10 * mm, "Bile acid analysis comparison audit")
    canvas.drawRightString(doc.pagesize[0] - 18 * mm, 10 * mm, f"Page {doc.page}")
    canvas.restoreState()


styles = getSampleStyleSheet()
styles.add(
    ParagraphStyle(
        name="TitleCenter",
        parent=styles["Title"],
        fontSize=20,
        leading=24,
        alignment=TA_CENTER,
        textColor=colors.HexColor("#17324D"),
        spaceAfter=6,
    )
)
styles.add(
    ParagraphStyle(
        name="Subtle",
        parent=styles["BodyText"],
        fontSize=8,
        leading=11,
        textColor=colors.HexColor("#546A7B"),
        alignment=TA_CENTER,
        spaceAfter=12,
    )
)
styles.add(
    ParagraphStyle(
        name="H1Blue",
        parent=styles["Heading1"],
        fontSize=14,
        leading=18,
        textColor=colors.HexColor("#17324D"),
        spaceBefore=10,
        spaceAfter=6,
    )
)
styles.add(
    ParagraphStyle(
        name="H2Blue",
        parent=styles["Heading2"],
        fontSize=11,
        leading=14,
        textColor=colors.HexColor("#1F4E79"),
        spaceBefore=8,
        spaceAfter=4,
    )
)
styles.add(
    ParagraphStyle(
        name="Small",
        parent=styles["BodyText"],
        fontSize=8,
        leading=10,
    )
)
styles.add(
    ParagraphStyle(
        name="Body",
        parent=styles["BodyText"],
        fontSize=9,
        leading=12,
        spaceAfter=5,
    )
)
styles.add(
    ParagraphStyle(
        name="Callout",
        parent=styles["BodyText"],
        fontSize=9,
        leading=12,
        leftIndent=8,
        rightIndent=8,
        spaceBefore=4,
        spaceAfter=8,
        textColor=colors.HexColor("#17324D"),
    )
)


lmm = pd.read_csv(WORK / "lmm_composite_check.csv")
pair_report = pd.read_csv(WORK / "r_pairwise_report_style_check.csv")
pair_corr = pd.read_csv(WORK / "r_pairwise_corrected_missing_check.csv")

doc = SimpleDocTemplate(
    str(OUT),
    pagesize=PAGE_SIZE,
    rightMargin=15 * mm,
    leftMargin=15 * mm,
    topMargin=15 * mm,
    bottomMargin=16 * mm,
)

story = []
story.append(p("Bile Acid Report Comparison and Verification Audit", styles["TitleCenter"]))
story.append(
    p(
        "Inputs reviewed: Bile_Acid_LMM_Analysis_Report.pdf, Bile_Acid_Analysis_Report.pdf, "
        "Bile_Acid_LMM_Results.xlsx, IRMA23930_MS-BA_Results.csv, and the local R scripts. "
        "Generated May 20, 2026.",
        styles["Subtle"],
    )
)

story.append(p("Executive Verdict", styles["H1Blue"]))
story.append(
    p(
        "<b>Most reliable answer for the main post-meal question:</b> the audited finite-sample "
        "LMM results in <i>Bile_Acid_LMM_Analysis_Report.pdf</i> are the stronger result set. "
        "The original XLSX p-values are reproducible, but they are asymptotic z-test p-values; "
        "the finite-sample Satterthwaite/Kenward-Roger correction is more appropriate for 32 subjects.",
        styles["Body"],
    )
)
story.append(
    p(
        "<b>Use the CSV/R PDF with caution.</b> Its individual bile-acid pairwise model estimates are reproducible, "
        "but the PDF is not a clean final statistical report: it reports raw p-value significance rather than FDR q-values, "
        "answers a different contrast structure, contains aggregate-class missing-data problems, and includes narrative "
        "claims that contradict its own tables.",
        styles["Body"],
    )
)
story.append(
    p(
        "<b>Not an apples-to-apples disagreement:</b> the LMM report tests timepoint 23 against the average of the fasting "
        "states (2, 18, 26), while the CSV/R report tests three separate pairwise contrasts: 23 vs 2, 23 vs 18, and 23 vs 26.",
        styles["Body"],
    )
)

story.append(p("Data Lineage Checks", styles["H1Blue"]))
data_rows = [
    ["Check", "Verification result", "Implication"],
    [
        "Raw CSV vs XLSX raw sheet",
        "Same 203 x 19 assay table; maximum numeric difference approximately 1.4e-14.",
        "Both reports start from the same raw data.",
    ],
    [
        "XLSX batch-adjusted sheet",
        "Exactly matches LTR mean-based plate factors (max difference approximately 9.1e-13).",
        "The XLSX method is mean-based, not median-based.",
    ],
    [
        "CSV/R script batch adjustment",
        "Uses LTR median-based plate factors.",
        "This is defensible but different from the XLSX pipeline; results are not expected to be identical.",
    ],
    [
        "Timepoints analyzed",
        "Both exclude timepoint 10 and use 2, 18, 23, and 26.",
        "The core longitudinal subset is aligned.",
    ],
]
table = Table(data_rows, colWidths=[38 * mm, 78 * mm, 62 * mm], repeatRows=1)
story.append(style_table(table, font_size=7))

story.append(p("Major Method Differences", styles["H1Blue"]))
method_rows = [
    ["Area", "LMM/XLSX audit report", "CSV/R analysis report"],
    ["Batch correction", "LTR mean-based adjustment in the workbook.", "LTR median-based adjustment in run_analysis.R."],
    ["Transform", "log10(concentration).", "Natural log of concentration; adds 1 if any zero is present."],
    ["Primary endpoint", "23 vs average fasting state (2, 18, 26).", "Separate 23 vs 2, 23 vs 18, and 23 vs 26 contrasts."],
    ["Multiple testing", "BH FDR over 15 bile acids; PDF explicitly audits asymptotic vs finite-sample p-values.", "Existing PDF shows raw p-value stars. Current R script adds FDR, but the PDF in the folder is stale."],
    ["Degrees of freedom", "Audit identifies finite-sample correction as best practice.", "Text says Satterthwaite; current emmeans default in this environment is Kenward-Roger unless lmer.df is set explicitly."],
    ["Aggregate classes", "Not the focus of the workbook LMM results.", "Aggregate classes have a missing-data bug: all-missing component rows become zero via rowSums(..., na.rm=TRUE)."],
]
table = Table(method_rows, colWidths=[33 * mm, 72 * mm, 73 * mm], repeatRows=1)
story.append(style_table(table, font_size=6.7))

story.append(PageBreak())
story.append(p("Verified Key Results", styles["H1Blue"]))
story.append(
    p(
        "For the LMM/XLSX endpoint (23 vs average fasting), the workbook p-values match R only when "
        "using asymptotic degrees of freedom. Recomputing the same model with finite-sample correction changes the borderline calls.",
        styles["Body"],
    )
)

sig = lmm[lmm["R_Satterthwaite_FDR_q"] < 0.05].copy()
sig_rows = [["Compound", "Fold change 23 / fasting avg", "Satterthwaite p", "Satterthwaite FDR q", "Sig"]]
for _, r in sig.iterrows():
    sig_rows.append(
        [
            r["Compound"],
            fmt_num(r["R_Satterthwaite_Fold"], 2),
            fmt_p(r["R_Satterthwaite_Post_p"]),
            fmt_p(r["R_Satterthwaite_FDR_q"]),
            r["R_Satterthwaite_Sig"],
        ]
    )
table = Table(sig_rows, colWidths=[30 * mm, 46 * mm, 34 * mm, 38 * mm, 20 * mm], repeatRows=1)
story.append(style_table(table, font_size=7.2))

story.append(p("Borderline And Changed Calls", styles["H2Blue"]))
border = lmm[(lmm["Compound"].isin(["GCDCA", "LCA", "TLCA", "TCA"]))].copy()
border_rows = [["Compound", "XLSX/asymptotic FDR q", "XLSX sig", "Finite-sample FDR q", "Finite-sample sig", "Interpretation"]]
notes = {
    "GCDCA": "Borderline, remains non-significant after FDR.",
    "LCA": "False positive in asymptotic XLSX call; finite-sample q is just above 0.05.",
    "TLCA": "Still significant, but one star weaker.",
    "TCA": "Still significant, but one star weaker.",
}
for _, r in border.iterrows():
    border_rows.append(
        [
            r["Compound"],
            fmt_p(r["Xlsx_Post_FDR_q"]),
            r["Xlsx_Sig"],
            fmt_p(r["R_Satterthwaite_FDR_q"]),
            r["R_Satterthwaite_Sig"],
            notes[r["Compound"]],
        ]
    )
table = Table(border_rows, colWidths=[22 * mm, 34 * mm, 22 * mm, 34 * mm, 26 * mm, 49 * mm], repeatRows=1)
story.append(style_table(table, font_size=6.8))

story.append(p("CSV/R Pairwise Results After FDR", styles["H2Blue"]))
ba_cols = ["CDCA", "GCDCA", "TCDCA", "UDCA", "GUDCA", "TUDCA", "DCA", "GDCA", "LCA", "GLCA", "TLCA", "CA", "GCA", "TCA", "TDCA"]
pair_rows = [["Contrast", "FDR-significant individual bile acids", "Important correction to PDF interpretation"]]
for label, qcol, note in [
    ("23 vs baseline 2", "Q_23_vs_2", "Consistent strong post-meal signal in conjugated acids."),
    ("23 vs NS4 fasting 18", "Q_23_vs_18", "GCA is raw-p significant in the PDF but not FDR significant."),
    ("23 vs RTDS fasting 26", "Q_23_vs_26", "Only LCA and GLCA survive FDR."),
]:
    comps = pair_report[pair_report["Compound"].isin(ba_cols) & (pair_report[qcol] < 0.05)]["Compound"].tolist()
    pair_rows.append([label, ", ".join(comps), note])
table = Table(pair_rows, colWidths=[38 * mm, 88 * mm, 54 * mm], repeatRows=1)
story.append(style_table(table, font_size=7))

story.append(PageBreak())
story.append(p("Aggregate-Class Issue In The CSV/R Report", styles["H1Blue"]))
story.append(
    p(
        "The aggregate-class section is the least reliable part of the CSV/R PDF. The script calculates aggregate pools with "
        "rowSums(..., na.rm = TRUE). When every component in an aggregate is missing, R returns 0, not NA. The script then includes "
        "that row and models log(value + 1). This created one extra aggregate observation (N=127 instead of N=126) and materially changed p-values.",
        styles["Body"],
    )
)
merged = pair_report.merge(pair_corr, on="Compound", suffixes=("_pdf", "_corrected"))
agg = merged[merged["Compound"].isin(["Total_BA", "Glycine_Conjugated", "Taurine_Conjugated", "Unconjugated"])]
agg_rows = [
    [
        "Aggregate",
        "PDF-style N",
        "PDF-style q 23 vs 2",
        "PDF-style sig",
        "Corrected N",
        "Corrected q 23 vs 2",
        "Corrected sig",
        "Corrected q 23 vs 18",
        "Corrected sig",
    ]
]
for _, r in agg.iterrows():
    agg_rows.append(
        [
            r["Compound"],
            str(int(r["N_Obs_pdf"])),
            fmt_p(r["Q_23_vs_2_pdf"]),
            r["Sig_23_vs_2_pdf"],
            str(int(r["N_Obs_corrected"])),
            fmt_p(r["Q_23_vs_2_corrected"]),
            r["Sig_23_vs_2_corrected"],
            fmt_p(r["Q_23_vs_18_corrected"]),
            r["Sig_23_vs_18_corrected"],
        ]
    )
table = Table(agg_rows, colWidths=[35 * mm, 18 * mm, 27 * mm, 18 * mm, 20 * mm, 30 * mm, 20 * mm, 30 * mm, 20 * mm], repeatRows=1)
story.append(style_table(table, font_size=6.2))

story.append(p("Internal Narrative Inconsistencies In The CSV/R PDF", styles["H2Blue"]))
issue_rows = [
    ["Location", "What the table says", "What the narrative says", "Why it matters"],
    [
        "Glycine-conjugated aggregate",
        "23 vs baseline p=0.0556, non-significant in the PDF table.",
        "Claims a highly significant spike compared to all fasting timepoints.",
        "Conclusion overstates the reported table result.",
    ],
    [
        "Taurine-conjugated aggregate",
        "23 vs RTDS fasting p=0.0985, non-significant.",
        "Claims significance relative to all fasting points.",
        "The RTDS comparison does not support that statement.",
    ],
    [
        "Unconjugated aggregate",
        "All three contrasts are non-significant in the table.",
        "Claims p<0.001 vs baseline and fasting 18.",
        "This is a direct contradiction and should not be used.",
    ],
]
table = Table(issue_rows, colWidths=[36 * mm, 50 * mm, 50 * mm, 42 * mm], repeatRows=1)
story.append(style_table(table, font_size=6.4, header_bg="#7A2E2E"))

story.append(p("Recommended Final Position", styles["H1Blue"]))
rec_rows = [
    ["Question", "Recommended answer"],
    [
        "Which report is more likely correct?",
        "Use the LMM audit report's finite-sample composite results as the primary answer for post-meal divergence.",
    ],
    [
        "Can the CSV/R report be used?",
        "Use only as a pairwise exploratory analysis after applying FDR q-values and regenerating the PDF from the current script.",
    ],
    [
        "What should be fixed before publication?",
        "Choose one LTR normalization rule, set lmer.df explicitly, report FDR q-values, fix aggregate missingness, and rewrite narrative text from the tables.",
    ],
    [
        "Biological conclusion that survives both approaches",
        "The strongest and most consistent postprandial signals are concentrated in conjugated bile acids, especially TCDCA, GLCA, TLCA, GCA, TCA, and TDCA.",
    ],
]
table = Table(rec_rows, colWidths=[52 * mm, 126 * mm], repeatRows=1)
story.append(style_table(table, font_size=7.2, header_bg="#2F5D50"))

story.append(Spacer(1, 6))
story.append(
    p(
        "Audit artifacts generated in comparison_work: lmm_composite_check.csv, "
        "r_pairwise_report_style_check.csv, and r_pairwise_corrected_missing_check.csv.",
        styles["Small"],
    )
)

doc.build(story, onFirstPage=footer, onLaterPages=footer)
print(OUT)
