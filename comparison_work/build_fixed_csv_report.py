from pathlib import Path

import pandas as pd
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.platypus import (
    Image,
    PageBreak,
    Paragraph,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)


BASE = Path(r"C:\Users\a3188798\OneDrive - Adelaide University\Desktop\NS BA AG")
OUT_DIR = BASE / "comparison_work" / "fixed_csv_report_outputs"
PLOT_DIR = OUT_DIR / "plots"
PDF_PATH = BASE / "Bile_Acid_Analysis_Report_FIXED_Verified.pdf"
PAGE_SIZE = landscape(A4)


def fmt_p(x):
    if pd.isna(x):
        return "NA"
    x = float(x)
    if x < 0.001:
        return "<0.001"
    return f"{x:.3f}"


def fmt_num(x, digits=3):
    if pd.isna(x):
        return "NA"
    return f"{float(x):.{digits}f}"


def paragraph(text, style):
    return Paragraph(str(text), style)


def make_table(rows, widths, font_size=6.5, header="#17324D", align="CENTER"):
    tbl = Table(rows, colWidths=widths, repeatRows=1)
    tbl.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor(header)),
                ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
                ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                ("FONTSIZE", (0, 0), (-1, -1), font_size),
                ("LEADING", (0, 0), (-1, -1), font_size + 2),
                ("ALIGN", (0, 0), (-1, -1), align),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("GRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#D9E2EC")),
                ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#F7FAFC")]),
                ("LEFTPADDING", (0, 0), (-1, -1), 3),
                ("RIGHTPADDING", (0, 0), (-1, -1), 3),
                ("TOPPADDING", (0, 0), (-1, -1), 3),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
            ]
        )
    )
    return tbl


def footer(canvas, doc):
    canvas.saveState()
    canvas.setFont("Helvetica", 8)
    canvas.setFillColor(colors.HexColor("#66788A"))
    canvas.drawString(15 * mm, 9 * mm, "Corrected CSV/R bile acid LMM report")
    canvas.drawRightString(doc.pagesize[0] - 15 * mm, 9 * mm, f"Page {doc.page}")
    canvas.restoreState()


def result_sentence(row):
    sig_bits = []
    checks = [
        ("baseline fasting (2)", row["LFC_23_vs_2"], row["Q_23_vs_2"], row["Sig_23_vs_2"]),
        ("NS4 fasting (18)", row["LFC_23_vs_18"], row["Q_23_vs_18"], row["Sig_23_vs_18"]),
        ("RTDS fasting (26)", row["LFC_23_vs_26"], row["Q_23_vs_26"], row["Sig_23_vs_26"]),
    ]
    for label, lfc, q, sig in checks:
        if float(q) < 0.05:
            direction = "higher" if float(lfc) > 0 else "lower"
            sig_bits.append(f"{direction} than {label} (FDR q={fmt_p(q)}, {sig})")
    if not sig_bits:
        return (
            "After BH FDR correction, timepoint 23 was not significantly different from "
            "baseline fasting, NS4 fasting, or RTDS fasting for this outcome."
        )
    return "After BH FDR correction, timepoint 23 was " + "; ".join(sig_bits) + "."


styles = getSampleStyleSheet()
styles.add(
    ParagraphStyle(
        "TitleCenter",
        parent=styles["Title"],
        alignment=TA_CENTER,
        fontSize=20,
        leading=24,
        textColor=colors.HexColor("#17324D"),
        spaceAfter=8,
    )
)
styles.add(
    ParagraphStyle(
        "Subtitle",
        parent=styles["BodyText"],
        alignment=TA_CENTER,
        fontSize=9,
        leading=12,
        textColor=colors.HexColor("#546A7B"),
        spaceAfter=12,
    )
)
styles.add(
    ParagraphStyle(
        "H1",
        parent=styles["Heading1"],
        fontSize=14,
        leading=17,
        textColor=colors.HexColor("#17324D"),
        spaceBefore=8,
        spaceAfter=6,
    )
)
styles.add(
    ParagraphStyle(
        "H2",
        parent=styles["Heading2"],
        fontSize=11,
        leading=14,
        textColor=colors.HexColor("#1F4E79"),
        spaceBefore=7,
        spaceAfter=4,
    )
)
styles.add(
    ParagraphStyle(
        "Body",
        parent=styles["BodyText"],
        fontSize=9,
        leading=12,
        spaceAfter=5,
    )
)
styles.add(
    ParagraphStyle(
        "Small",
        parent=styles["BodyText"],
        fontSize=7.2,
        leading=9,
    )
)
styles.add(
    ParagraphStyle(
        "TableText",
        parent=styles["BodyText"],
        fontSize=6.5,
        leading=8.2,
        alignment=TA_LEFT,
    )
)

results = pd.read_csv(OUT_DIR / "corrected_lmm_results.csv")
verification = pd.read_csv(OUT_DIR / "verification_checks.csv")

doc = SimpleDocTemplate(
    str(PDF_PATH),
    pagesize=PAGE_SIZE,
    rightMargin=14 * mm,
    leftMargin=14 * mm,
    topMargin=14 * mm,
    bottomMargin=14 * mm,
)

story = []
story.append(paragraph("Corrected Bile Acid Analysis Report", styles["TitleCenter"]))
story.append(
    paragraph(
        "Fixed and verified CSV/R-derived longitudinal mixed-model analysis of fasting and postprandial bile-acid dynamics. "
        "Source data: IRMA23930_MS-BA_Results.csv. Generated from corrected scripts in comparison_work.",
        styles["Subtitle"],
    )
)

story.append(paragraph("Statistical Methods Guide", styles["H1"]))
story.append(
    paragraph(
        "This report uses a longitudinal linear mixed-effects model because each subject contributes repeated bile-acid "
        "measurements across multiple timepoints. A standard linear model or independent t-test would ignore within-subject "
        "correlation. The model used for each bile acid or aggregate class was <b>log(concentration) ~ timepoint + (1 | subject)</b>, "
        "where timepoint is the fixed effect of interest and subject is a random intercept that accounts for each participant's "
        "baseline level.",
        styles["Body"],
    )
)

method_guide_rows = [
    ["Term", "Meaning in this report"],
    [
        "Log transform",
        paragraph(
            "Concentrations were analyzed on the natural-log scale. This stabilizes variance and makes skewed concentration "
            "data more suitable for mixed-model inference.",
            styles["TableText"],
        ),
    ],
    [
        "LFC",
        paragraph(
            "Log fold change. Here it is the model-estimated difference on the natural-log scale, such as log(23) - log(2). "
            "Positive LFC means the post-meal timepoint is higher; negative LFC means it is lower.",
            styles["TableText"],
        ),
    ],
    [
        "Fold",
        paragraph(
            "Back-transformed LFC: Fold = exp(LFC). A fold of 1.50 means 50% higher; 0.75 means 25% lower.",
            styles["TableText"],
        ),
    ],
    [
        "EMM",
        paragraph(
            "Estimated marginal mean from the LMM. These are model-adjusted means for each timepoint, back-transformed to "
            "the concentration scale for plots.",
            styles["TableText"],
        ),
    ],
    [
        "Satterthwaite",
        paragraph(
            "A finite-sample degrees-of-freedom approximation for mixed models. It gives more appropriate p-values for "
            "moderate repeated-measures datasets than treating test statistics as if sample size were infinite.",
            styles["TableText"],
        ),
    ],
    [
        "BH FDR",
        paragraph(
            "Benjamini-Hochberg false discovery rate correction. Because many contrasts are tested, raw p-values are adjusted "
            "to q-values to control the expected proportion of false discoveries among significant results.",
            styles["TableText"],
        ),
    ],
    [
        "q-value",
        paragraph(
            "The FDR-adjusted p-value. This report calls results significant using q < 0.05, not raw p < 0.05.",
            styles["TableText"],
        ),
    ],
]
story.append(make_table(method_guide_rows, [42 * mm, 196 * mm], font_size=7.0, header="#1F4E79", align="LEFT"))

story.append(paragraph("Contrasts And Significance Rules", styles["H2"]))
contrast_guide_rows = [
    ["Analysis", "How it was done"],
    [
        "Main contrasts",
        paragraph(
            "Timepoint 23 was compared separately against baseline fasting (2), NS4 fasting (18), and RTDS fasting (26). "
            "These are shown as 23 vs 2, 23 vs 18, and 23 vs 26.",
            styles["TableText"],
        ),
    ],
    [
        "Multiple testing families",
        paragraph(
            "FDR correction was applied separately to the 45 individual bile-acid contrasts (15 bile acids x 3 contrasts) "
            "and the 12 aggregate-class contrasts (4 aggregate classes x 3 contrasts).",
            styles["TableText"],
        ),
    ],
    [
        "Significance stars",
        paragraph(
            "Stars are based on FDR q-values: *** q<0.001, ** q<0.01, * q<0.05, ns = not significant.",
            styles["TableText"],
        ),
    ],
    [
        "Plot layout",
        paragraph(
            "Plots use the requested divergence layout: blue branch 2->18->26 and orange branch 18->23->26. "
            "Only the timepoint 23 marker is orange because it is the meal branch point.",
            styles["TableText"],
        ),
    ],
    [
        "Missing values",
        paragraph(
            "Missing values were not imputed. The all-missing sample 1623-26 was kept missing for aggregate classes rather "
            "than converted to zero.",
            styles["TableText"],
        ),
    ],
]
story.append(make_table(contrast_guide_rows, [50 * mm, 188 * mm], font_size=7.0, header="#2F5D50", align="LEFT"))

story.append(PageBreak())

story.append(paragraph("What Was Fixed", styles["H1"]))
fix_rows = [
    ["Problem in prior PDF", "Correction in this report"],
    ["Raw p-value stars were used for interpretation.", "Tables and text use BH FDR q-values for significance."],
    ["Aggregate rowSums(..., na.rm=TRUE) converted an all-missing BA row to zero.", "All-missing aggregate rows are preserved as NA and excluded compound-by-compound."],
    ["Model df method was not explicit.", "emmeans is called with lmer.df='satterthwaite' for finite-sample inference."],
    ["Narrative overclaimed significance for several aggregate classes.", "Interpretation text is generated from the corrected q-values."],
]
story.append(make_table(fix_rows, [68 * mm, 170 * mm], font_size=7.2, header="#2F5D50", align="LEFT"))

story.append(paragraph("Corrected Statistical Method", styles["H1"]))
story.append(
    paragraph(
        "For each individual bile acid and each aggregate class, adjusted positive concentrations were modeled as "
        "<b>log(concentration) ~ timepoint + (1 | subject)</b>. Timepoint 23 was compared with baseline fasting (2), "
        "NS4 fasting (18), and RTDS fasting (26). P-values use Satterthwaite finite-sample degrees of freedom. "
        "BH FDR correction was applied separately to the 45 individual-bile-acid contrasts and the 12 aggregate-class contrasts.",
        styles["Body"],
    )
)

story.append(paragraph("Verification Checks", styles["H1"]))
ver_rows = [["Check", "Result"]]
for _, r in verification.iterrows():
    ver_rows.append([r["Check"], r["Result"]])
story.append(make_table(ver_rows, [86 * mm, 152 * mm], font_size=7.2, header="#1F4E79", align="LEFT"))

story.append(paragraph("High-Level Corrected Result", styles["H1"]))
story.append(
    paragraph(
        "The strongest corrected postprandial signals remain concentrated in conjugated bile acids. "
        "Against baseline fasting (2), GCDCA, TCDCA, GLCA, TLCA, GCA, TCA, and TDCA are FDR-significant. "
        "Against NS4 fasting (18), TCDCA, LCA, GLCA, TLCA, TCA, and TDCA are FDR-significant. "
        "Against RTDS fasting (26), only LCA and GLCA survive FDR among individual bile acids.",
        styles["Body"],
    )
)

story.append(PageBreak())

contrast_specs = [
    ("23 vs Baseline Fasting (2)", "LFC_23_vs_2", "Fold_23_vs_2", "P_23_vs_2", "Q_23_vs_2", "Sig_23_vs_2"),
    ("23 vs NS4 Fasting (18)", "LFC_23_vs_18", "Fold_23_vs_18", "P_23_vs_18", "Q_23_vs_18", "Sig_23_vs_18"),
    ("23 vs RTDS Fasting (26)", "LFC_23_vs_26", "Fold_23_vs_26", "P_23_vs_26", "Q_23_vs_26", "Sig_23_vs_26"),
]

for title, lfc_col, fold_col, p_col, q_col, sig_col in contrast_specs:
    story.append(paragraph(f"Master Results Table: {title}", styles["H1"]))
    rows = [["Compound", "Class", "Family", "N", "LFC", "Fold", "raw p", "FDR q", "Sig"]]
    for _, r in results.iterrows():
        rows.append(
            [
                r["Compound"],
                paragraph(r["Class"], styles["TableText"]),
                r["Result_Family"],
                str(int(r["N_Obs"])),
                fmt_num(r[lfc_col]),
                fmt_num(r[fold_col], 2),
                fmt_p(r[p_col]),
                fmt_p(r[q_col]),
                r[sig_col],
            ]
        )
    story.append(
        make_table(
            rows,
            [22 * mm, 46 * mm, 32 * mm, 14 * mm, 22 * mm, 20 * mm, 22 * mm, 22 * mm, 14 * mm],
            font_size=6.6,
            header="#17324D",
        )
    )
    story.append(Spacer(1, 4))
    if title != contrast_specs[-1][0]:
        story.append(PageBreak())

story.append(PageBreak())
story.append(paragraph("Trajectory Overview", styles["H1"]))
story.append(
    paragraph(
        "Plots show LMM estimated marginal means back-transformed to the concentration scale. Error bars are 95% confidence intervals. "
        "All trajectories use a divergence layout: the fasting branch runs from baseline (2) through NS4 fasting (18) to "
        "RTDS fasting (26), while the meal branch starts at NS4 fasting (18), passes through NS4 after-meal (23), and "
        "merges again at RTDS fasting (26). "
        "The y-axis is log10-scaled for readability; significance labels in the detailed plots use FDR q-values.",
        styles["Body"],
    )
)
grid = Image(str(PLOT_DIR / "all_15_corrected_trajectories_grid.png"), width=126 * mm, height=159 * mm)
agg = Image(str(PLOT_DIR / "aggregate_corrected_trajectories_grid.png"), width=110 * mm, height=80 * mm)
story.append(Table([[grid, agg]], colWidths=[132 * mm, 116 * mm]))

for _, r in results.iterrows():
    story.append(PageBreak())
    story.append(paragraph(f"{r['Compound']}: {r['Compound_Name']}", styles["H1"]))
    story.append(paragraph(f"{r['Class']} | {r['Result_Family']} | N={int(r['N_Obs'])}, subjects={int(r['N_Subjects'])}", styles["Body"]))
    img = Image(str(PLOT_DIR / f"{r['Compound']}_corrected_trajectory.png"), width=142 * mm, height=96 * mm)
    detail_rows = [
        ["Contrast", "LFC", "Fold", "raw p", "FDR q", "Sig"],
        ["23 vs 2", fmt_num(r["LFC_23_vs_2"]), fmt_num(r["Fold_23_vs_2"], 2), fmt_p(r["P_23_vs_2"]), fmt_p(r["Q_23_vs_2"]), r["Sig_23_vs_2"]],
        ["23 vs 18", fmt_num(r["LFC_23_vs_18"]), fmt_num(r["Fold_23_vs_18"], 2), fmt_p(r["P_23_vs_18"]), fmt_p(r["Q_23_vs_18"]), r["Sig_23_vs_18"]],
        ["23 vs 26", fmt_num(r["LFC_23_vs_26"]), fmt_num(r["Fold_23_vs_26"], 2), fmt_p(r["P_23_vs_26"]), fmt_p(r["Q_23_vs_26"]), r["Sig_23_vs_26"]],
    ]
    detail = make_table(
        detail_rows,
        [25 * mm, 19 * mm, 18 * mm, 20 * mm, 20 * mm, 14 * mm],
        font_size=7.0,
        header="#1F4E79",
    )
    side = Table([[img, detail]], colWidths=[148 * mm, 118 * mm])
    side.setStyle(TableStyle([("VALIGN", (0, 0), (-1, -1), "TOP")]))
    story.append(side)
    story.append(Spacer(1, 6))
    story.append(paragraph(f"<b>FDR-based interpretation:</b> {result_sentence(r)}", styles["Body"]))
    if r["Result_Family"] == "Aggregate class":
        story.append(
            paragraph(
                "Aggregate-class inference is derived from summed component bile acids after preserving all-missing rows as missing. "
                "These aggregate tests should be interpreted as class-level summaries, not replacements for the individual bile-acid results.",
                styles["Small"],
            )
        )

story.append(PageBreak())
story.append(paragraph("Reproducibility Files", styles["H1"]))
files = [
    ["File", "Purpose"],
    ["comparison_work/fixed_csv_lmm_analysis.R", "Recomputes corrected LMM results, q-values, verification checks, and plots."],
    ["comparison_work/fixed_csv_report_outputs/corrected_lmm_results.csv", "Primary corrected result table used for this PDF."],
    ["comparison_work/fixed_csv_report_outputs/corrected_emm_results.csv", "Back-transformed estimated marginal means and confidence intervals."],
    ["comparison_work/fixed_csv_report_outputs/verification_checks.csv", "Audit checks printed in this report."],
    ["comparison_work/fixed_csv_report_outputs/session_info.txt", "R version and package session details."],
]
story.append(make_table(files, [92 * mm, 150 * mm], font_size=7.2, header="#2F5D50", align="LEFT"))

doc.build(story, onFirstPage=footer, onLaterPages=footer)
print(PDF_PATH)
