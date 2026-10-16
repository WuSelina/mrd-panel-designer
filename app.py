"""
Tumor-Informed MRD Panel Designer — TCGA KIRC
Teaching demo: designs a ctDNA panel from somatic mutations in clear cell RCC.
Run: streamlit run app.py --server.port 8081
"""

from typing import Optional

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st
from scipy.stats import poisson

# ── Page config — must be the very first Streamlit call ──────────────────────
st.set_page_config(
    page_title="MRD Panel Designer — KIRC",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ── Constants ─────────────────────────────────────────────────────────────────
CHIP_GENES = {"DNMT3A", "TET2", "ASXL1", "PPM1D", "JAK2"}
DATA_PATH = "data/data_mutations.txt"


# ── Data loader ───────────────────────────────────────────────────────────────
@st.cache_data(show_spinner="Loading mutations…")
def load_mutations(path: str = DATA_PATH) -> pd.DataFrame:
    """Load the TCGA KIRC MAF, keep only the columns we need, compute VAF."""
    cols = [
        "Tumor_Sample_Barcode",
        "Hugo_Symbol",
        "Variant_Classification",
        "t_alt_count",
        "t_depth",
    ]
    df = pd.read_csv(path, sep="\t", low_memory=False, usecols=cols)
    df["t_alt_count"] = pd.to_numeric(df["t_alt_count"], errors="coerce")
    df["t_depth"] = pd.to_numeric(df["t_depth"], errors="coerce")
    # Drop rows with missing or zero depth (can't compute VAF)
    df = df[df["t_depth"].notna() & (df["t_depth"] > 0)].copy()
    df["VAF"] = df["t_alt_count"] / df["t_depth"]
    return df


def lod_tf(ccf_sum: float, D: int, k: int, target_p: float) -> Optional[float]:
    """
    Binary-search for the plasma tumor fraction (TF) where P(detect) = target_p.
    Returns None if the target is not reached within the search range.
    """
    lo, hi = 1e-7, 1.0
    for _ in range(60):
        mid = (lo + hi) / 2
        lam = D * mid * ccf_sum / 2.0
        p = 1.0 - poisson.cdf(k - 1, lam)
        if p < target_p:
            lo = mid
        else:
            hi = mid
    tf = (lo + hi) / 2
    # Return None if P never reaches target (e.g. panel has zero CCF sum)
    lam_max = D * 1.0 * ccf_sum / 2.0
    if 1.0 - poisson.cdf(k - 1, lam_max) < target_p:
        return None
    return tf


# ── Load data ─────────────────────────────────────────────────────────────────
try:
    df = load_mutations()
except FileNotFoundError:
    st.error(
        f"**Data file not found:** `{DATA_PATH}`\n\n"
        "Please follow the README instructions to download the TCGA KIRC MAF "
        "into `data/data_mutations.txt`."
    )
    st.stop()

all_samples = sorted(df["Tumor_Sample_Barcode"].unique())

# ── Sidebar ───────────────────────────────────────────────────────────────────
st.sidebar.title("⚙️ Panel Controls")

# ── Sample selector ───────────────────────────────────────────────────────────
sample = st.sidebar.selectbox("Tumor sample", all_samples)

# ── Filter for selected sample ────────────────────────────────────────────────
sample_df = df[df["Tumor_Sample_Barcode"] == sample].copy()

n_chip = int(sample_df["Hugo_Symbol"].isin(CHIP_GENES).sum())
n_low_depth = int((sample_df["t_depth"] < 30).sum())

filtered = sample_df[
    (~sample_df["Hugo_Symbol"].isin(CHIP_GENES)) & (sample_df["t_depth"] >= 30)
].copy()

st.sidebar.metric("Mutations (after filter)", len(filtered))

# CHIP exclusion caption with ⓘ popover
chip_col, chip_info_col = st.sidebar.columns([5, 1])
with chip_col:
    st.caption(f"Excluded: {n_chip} CHIP, {n_low_depth} low-depth (<30×)")
with chip_info_col:
    with st.popover("ⓘ"):
        st.markdown(
            "**CHIP — Clonal Haematopoiesis of Indeterminate Potential**\n\n"
            "CHIP mutations arise in blood stem cells with age. Most cfDNA in "
            "plasma comes from blood cells, so a CHIP mutation shows up in plasma "
            "even when no tumor is present, which looks like a false-positive MRD "
            "call. We exclude ***DNMT3A*, *TET2*, *ASXL1*, *PPM1D*** and ***JAK2*** "
            "as a simple gene-level filter. Real assays instead sequence the "
            "patient's own white blood cells and remove any variant found there.\n\n"
            "MRD = Minimal Residual Disease; the detection of tumor-derived DNA "
            "in blood after treatment."
        )

# ── Purity inference (must happen before the slider is rendered) ──────────────
vhl_muts = filtered[filtered["Hugo_Symbol"] == "VHL"]
if len(vhl_muts) > 0:
    vhl_vaf = float(vhl_muts["VAF"].median())
    purity_inferred = 2.0 * vhl_vaf / (1.0 + vhl_vaf)
    purity_method = "vhl"
else:
    vhl_vaf = None
    purity_inferred = 0.6
    purity_method = "slider"

# ── Purity slider (sidebar) with ⓘ popover ───────────────────────────────────
st.sidebar.markdown("---")
pur_label_col, pur_info_col = st.sidebar.columns([5, 1])
with pur_label_col:
    st.markdown("**Purity override**")
with pur_info_col:
    with st.popover("ⓘ"):
        st.markdown(
            "**Tumor purity** is the fraction of cells in the biopsy that are tumor "
            "cells.\n\n"
            "In ccRCC (clear cell renal cell carcinoma), *VHL* is almost always "
            "truncally mutated on a background of **3p LOH** (loss of heterozygosity "
            "— one copy of chromosome 3p is deleted), making tumor cells hemizygous "
            "at that locus. For a hemizygous clonal mutation:\n\n"
            "$$\\text{VAF} = \\frac{p}{2 - p}$$\n\n"
            "Rearranging: $p = \\dfrac{2 \\cdot \\text{VAF}}{1 + \\text{VAF}}$\n\n"
            "This is more accurate than the naive diploid estimate "
            "$p = 2 \\cdot \\text{VAF}$, which would overestimate purity when one "
            "allele is lost."
        )

purity = st.sidebar.slider(
    "Purity",
    min_value=0.05,
    max_value=1.0,
    value=float(round(purity_inferred, 2)),
    step=0.01,
    label_visibility="collapsed",
)

# ── Panel size ────────────────────────────────────────────────────────────────
st.sidebar.markdown("---")
max_n = max(1, len(filtered))
n_panel = int(
    st.sidebar.number_input(
        "Panel size (N)",
        min_value=1,
        max_value=min(200, max_n),
        value=min(16, max_n),
        step=1,
    )
)

# ── Detection model sliders ───────────────────────────────────────────────────
st.sidebar.markdown("---")
st.sidebar.markdown("**Detection model**")

d_label_col, d_info_col = st.sidebar.columns([5, 1])
with d_label_col:
    st.markdown("**Unique molecules per locus (genome equivalents, GE)**")
with d_info_col:
    with st.popover("ⓘ"):
        st.markdown(
            "**D = genome equivalents (GE) per locus — not sequencing read depth.**\n\n"
            "cfDNA (cell-free DNA) input is limited by plasma volume: a standard "
            "10 mL blood draw yields roughly 4–5 mL of plasma, which contains on "
            "the order of 3,000–10,000 haploid GE per locus after cfDNA extraction "
            "and capture losses. Sequencing 100,000 reads per locus does not help "
            "if there are only 5,000 unique DNA molecules — extra reads just "
            "re-sequence the same molecules.\n\n"
            "After UMI (Unique Molecular Identifier) deduplication, **D is set by "
            "how much cfDNA you extracted**, not how deeply you sequenced. "
            "More plasma volume → higher D → better sensitivity."
        )
D = st.sidebar.slider(
    "Unique molecules per locus (GE)",
    min_value=500,
    max_value=50_000,
    value=5_000,
    step=500,
    label_visibility="collapsed",
)
k = st.sidebar.slider("Min mutant molecules to call (k)", min_value=1, max_value=10, value=2)
tf_ref_pct = st.sidebar.slider(
    "Reference TF line (%)",
    min_value=0.001,
    max_value=1.0,
    value=0.01,          # default 0.01%
    step=0.001,
    format="%.3f",
)
tf_ref = tf_ref_pct / 100.0

# ══════════════════════════════════════════════════════════════════════════════
# Main body
# ══════════════════════════════════════════════════════════════════════════════
st.title("🧬 Tumor-Informed MRD Panel Designer")
st.caption(
    "TCGA KIRC · PanCancer Atlas 2018 · Teaching Demo  "
    f"· {len(all_samples)} of 512 cohort samples have WXS (whole-exome sequencing) "
    "mutation calls in this MAF"
)

# ── Purity banner ─────────────────────────────────────────────────────────────
banner_col, banner_info_col = st.columns([11, 1])
with banner_col:
    if purity_method == "vhl":
        st.info(
            f"🔬 **Purity estimated from *VHL* VAF** ({vhl_vaf:.3f}) using the "
            f"hemizygous 3p-loss formula  "
            f"**p = 2·VAF / (1 + VAF) = {purity_inferred:.2f}**. "
            f"The slider on the left overrides this value (currently **{purity:.2f}**)."
        )
    else:
        st.warning(
            "⚠️ No *VHL* mutation found in this sample after filtering — "
            f"purity defaults to slider value (**{purity:.2f}**)."
        )
with banner_info_col:
    with st.popover("ⓘ"):
        st.markdown(
            "**VAF — Variant Allele Fraction**\n\n"
            "VAF = `t_alt_count ÷ t_depth`: the fraction of sequencing reads "
            "that carry the mutation at this locus.\n\n"
            "In a *pure* tumor, a clonal heterozygous (one mutant copy out of "
            "two) mutation on a diploid background has VAF ≈ 0.50. In a mixed "
            "biopsy (tumor + stroma + immune cells), VAF is lower because normal "
            "cells dilute the signal. With 3p LOH (one allele lost), a clonal "
            "mutation has VAF = p/(2−p), so we use the corrected purity formula."
        )

# ── CCF computation ───────────────────────────────────────────────────────────
filtered["CCF"] = (2.0 * filtered["VAF"] / purity).clip(upper=1.0)

# ── Variant classification filter ────────────────────────────────────────────
all_classes = sorted(filtered["Variant_Classification"].dropna().unique().tolist())

vc_col, vc_info_col = st.columns([11, 1])
with vc_col:
    selected_classes = st.multiselect(
        "Filter by Variant Classification (applied before top-N selection)",
        options=all_classes,
        default=all_classes,
    )
with vc_info_col:
    with st.popover("ⓘ"):
        st.markdown(
            "**Variant Classification filter**\n\n"
            "This filter controls which mutation types are eligible for the panel "
            "before the top-N selection by CCF (Cancer Cell Fraction).\n\n"
            "Silent (synonymous) mutations are **valid MRD markers** — we are "
            "tracking a tumor-specific DNA sequence, not a functional protein "
            "change. However, SNVs (single nucleotide variants) are generally more "
            "reliable than indels (insertions/deletions) in cfDNA assays because "
            "indels are harder to distinguish from PCR (polymerase chain reaction) "
            "and sequencing artefacts at low tumor fractions."
        )

if not selected_classes:
    st.warning("No variant classes selected — showing all mutations.")
    vc_filtered = filtered
else:
    vc_filtered = filtered[filtered["Variant_Classification"].isin(selected_classes)]

# ── Panel selection (top N by CCF, ties broken by t_depth descending) ─────────
panel_cols = [
    "Hugo_Symbol",
    "Variant_Classification",
    "t_alt_count",
    "t_depth",
    "VAF",
    "CCF",
]
panel = (
    vc_filtered
    .sort_values(["CCF", "t_depth"], ascending=[False, False])
    .head(n_panel)[panel_cols]
    .reset_index(drop=True)
)
panel.index += 1  # 1-based row numbers

# ── Panel table section ───────────────────────────────────────────────────────
tbl_col, tbl_info_col = st.columns([11, 1])
with tbl_col:
    st.subheader(f"Top {n_panel} of {len(vc_filtered)} Panel Mutations")
    st.caption(
        f"Showing the top {n_panel} mutations by CCF out of "
        f"{len(vc_filtered)} that passed all filters. "
        f"Adjust **Panel size (N)** in the sidebar to see more or fewer."
    )
with tbl_info_col:
    with st.popover("ⓘ"):
        st.markdown(
            "**CCF — Cancer Cell Fraction**\n\n"
            "CCF = `2 × VAF ÷ purity` (diploid assumption, capped at 1.0). "
            "It estimates the fraction of tumor cells that carry a given "
            "mutation.\n\n"
            "A CCF near **1.0** means the mutation is *clonal* — present in "
            "virtually all tumor cells — making it a more reliable MRD target "
            "than a subclonal mutation (CCF < 1.0), which may be absent in some "
            "tumor regions or metastases.\n\n"
            "Clonal mutations cluster around CCF = 1 at any purity, and sampling "
            "noise pushes about half of them above 1, so many hit the cap and tie. "
            "Ties are broken by sequencing depth (t_depth, descending) to prefer "
            "better-evidenced mutations."
        )

display = panel.copy()
display["VAF"] = display["VAF"].map("{:.3f}".format)
display["CCF"] = display["CCF"].map("{:.3f}".format)
st.dataframe(display, use_container_width=True)

st.download_button(
    label="⬇️ Download panel TSV",
    data=panel.to_csv(sep="\t", index_label="rank"),
    file_name=f"{sample}_mrd_panel.tsv",
    mime="text/tab-separated-values",
)

st.markdown("---")

# ── Detection model ───────────────────────────────────────────────────────────
det_col, det_info_col = st.columns([11, 1])
with det_col:
    st.subheader("Detection Probability Model")
with det_info_col:
    with st.popover("ⓘ"):
        st.markdown(
            "**CCF-weighted Poisson detection model**\n\n"
            "Each panel mutation contributes independently. The expected total "
            "mutant molecules across all N loci:\n\n"
            "$$\\lambda = \\sum_{i=1}^{N} D \\cdot TF \\cdot \\frac{CCF_i}{2}$$\n\n"
            "where **D** is the number of unique cfDNA molecules (GE — genome "
            "equivalents) per locus — set by plasma input, not sequencing depth — "
            "**TF** is the plasma tumor fraction, and **CCF_i / 2** is the "
            "expected mutant allele fraction for mutation *i* in plasma.\n\n"
            "**P(detect)** = P(Poisson(λ) ≥ k) = probability of observing at "
            "least k mutant molecules in total across the panel.\n\n"
            "**LoD** (Limit of Detection) is the TF at which P(detect) reaches "
            "50% (LoD50) or 95% (LoD95).\n\n"
            "*Simplifications:* this model ignores sequencing error rates, "
            "inter-mutation correlation, and panel capture efficiency."
        )

# Compute detection curve
tf_values = np.logspace(-5, -2, 300)          # 0.001% → 1%
ccf_sum = float(panel["CCF"].values.sum())
lam_vec = D * tf_values * ccf_sum / 2.0       # CCF-weighted lambda
p_detect = 1.0 - poisson.cdf(k - 1, lam_vec)

# P at the user-chosen reference TF
p_at_ref = float(1.0 - poisson.cdf(k - 1, D * tf_ref * ccf_sum / 2.0))

# LoD50 and LoD95
lod50 = lod_tf(ccf_sum, D, k, 0.50)
lod95 = lod_tf(ccf_sum, D, k, 0.95)

def fmt_lod(tf: Optional[float]) -> str:
    if tf is None:
        return "not reached"
    pct = tf * 100
    if pct < 0.01:
        return f"{pct:.4f}%"
    return f"{pct:.3f}%"

# Plotly figure
fig = go.Figure()
fig.add_trace(
    go.Scatter(
        x=tf_values * 100,
        y=p_detect,
        mode="lines",
        line=dict(color="royalblue", width=2.5),
        name="P(detect)",
        hovertemplate="TF = %{x:.4f}%<br>P(detect) = %{y:.3f}<extra></extra>",
    )
)

# User-adjustable reference TF line
fig.add_vline(
    x=tf_ref * 100,
    line_dash="dash",
    line_color="tomato",
    line_width=1.5,
    annotation_text=f"Ref TF = {tf_ref_pct:.3f}%",
    annotation_position="top right",
    annotation_font_color="tomato",
)

# LoD50 marker
if lod50 is not None:
    fig.add_scatter(
        x=[lod50 * 100], y=[0.50],
        mode="markers", marker=dict(color="green", size=10, symbol="diamond"),
        name=f"LoD50 = {fmt_lod(lod50)}",
        hovertemplate=f"LoD50: TF = {fmt_lod(lod50)}<extra></extra>",
    )

# LoD95 marker
if lod95 is not None:
    fig.add_scatter(
        x=[lod95 * 100], y=[0.95],
        mode="markers", marker=dict(color="darkorange", size=10, symbol="diamond"),
        name=f"LoD95 = {fmt_lod(lod95)}",
        hovertemplate=f"LoD95: TF = {fmt_lod(lod95)}<extra></extra>",
    )

# Explicit log-scale x-axis ticks at 0.001%, 0.01%, 0.1%, 1%
tick_vals = [0.001, 0.01, 0.1, 1.0]
tick_text = ["0.001%", "0.01%", "0.1%", "1%"]

fig.update_xaxes(
    type="log",
    title="Plasma Tumor Fraction (TF)",
    tickmode="array",
    tickvals=tick_vals,
    ticktext=tick_text,
)
fig.update_yaxes(title="P(detect \u2265 k mutant molecules)", range=[0, 1.05])
fig.update_layout(
    title=f"Detection Probability | N = {n_panel} mutations, D = {D:,} GE/locus, k \u2265 {k}",
    height=440,
    margin=dict(l=60, r=30, t=60, b=60),
    legend=dict(x=0.01, y=0.99),
)
st.plotly_chart(fig, use_container_width=True)

# Metrics row
# Metrics row: P(detect) | LoD50 | LoD95 | ⓘ
m1, m2, m3, m4 = st.columns([3, 3, 3, 1])
m1.metric(
    label=f"P(detect) at ref TF ({tf_ref_pct:.3f}%)",
    value=f"{p_at_ref:.1%}",
)
m2.metric(label="LoD50 (P = 50%)", value=fmt_lod(lod50))
m3.metric(label="LoD95 (P = 95%)", value=fmt_lod(lod95))
with m4:
    with st.popover("ⓘ"):
        st.markdown(
            "**LoD — Limit of Detection**\n\n"
            "LoD is the plasma tumor fraction (TF) at which the panel achieves a "
            "given detection probability under the Poisson model.\n\n"
            "- **LoD50**: the TF at which half of samples with that true TF would "
            "be called positive (50% sensitivity). Useful as a central estimate.\n"
            "- **LoD95**: the TF at which 95% of samples would be detected. "
            "This is the more clinically relevant threshold — below it, most "
            "positive samples would be missed.\n\n"
            "Both assume perfect specificity (no false positives), which requires "
            "error suppression in practice."
        )

st.markdown("---")

# ── Summary ───────────────────────────────────────────────────────────────────
st.subheader("Summary")

mean_ccf = float(panel["CCF"].mean())
n_at_cap = int((panel["CCF"] >= 1.0).sum())

purity_desc = (
    f"purity {purity:.2f} (*VHL* VAF {vhl_vaf:.3f} → {purity_inferred:.2f}, overridden via slider)"
    if purity_method == "vhl" and abs(purity - purity_inferred) > 0.01
    else (
        f"purity {purity:.2f} (from *VHL* VAF {vhl_vaf:.3f})"
        if purity_method == "vhl"
        else f"purity {purity:.2f} (slider default; no *VHL* mutation)"
    )
)

lod95_str = fmt_lod(lod95)
lod50_str = fmt_lod(lod50)

cap_warning = (
    " ⚠️ All mutations hit the CCF cap — purity may be underestimated; "
    "see caveats."
    if n_at_cap == n_panel else ""
)

st.markdown(
    f"**Technical summary — {sample}:** "
    f"{purity_desc[0].upper() + purity_desc[1:]}; {n_panel}-mutation panel; "
    f"mean CCF {mean_ccf:.2f} ({n_at_cap}/{n_panel} at cap){cap_warning}. "
    f"LoD50 = **{lod50_str}**, LoD95 = **{lod95_str}** "
    f"(D = {D:,} GE/locus, k ≥ {k})."
)

# Plain-language summary
if lod95 is not None:
    lod95_pct = lod95 * 100
    if lod95_pct <= 0.01:
        sensitivity_plain = (
            f"At its modeled sensitivity limit of {lod95_pct:.4f}% tumor DNA "
            f"in plasma, the panel would detect a positive sample in 95 out of "
            f"100 tests — a level that may be achievable with deep cfDNA sequencing."
        )
    elif lod95_pct <= 0.1:
        sensitivity_plain = (
            f"At its modeled sensitivity limit of {lod95_pct:.3f}% tumor DNA "
            f"in plasma, the panel would detect a positive sample in 95 out of "
            f"100 tests — within range of a well-optimized cfDNA assay."
        )
    else:
        sensitivity_plain = (
            f"At its modeled sensitivity limit of {lod95_pct:.2f}% tumor DNA "
            f"in plasma, the panel would detect a positive sample in 95 out of "
            f"100 tests — a relatively high threshold; the panel or input material "
            f"may need improvement."
        )
else:
    sensitivity_plain = (
        "The panel does not reach 95% detection probability within the modeled "
        "range (0.001%–1% tumor fraction). Consider increasing panel size or "
        "genome equivalents."
    )

clonal_desc = (
    f"All {n_panel} panel mutations are at the CCF = 1.0 cap — they are either "
    "truly clonal or the purity estimate is too low to distinguish clonal from "
    "subclonal (see caveats)"
    if n_at_cap == n_panel
    else (
        "Most panel mutations appear clonal (CCF near 1.0), making them stable "
        "targets expected to be present in any tumor cell"
        if mean_ccf >= 0.8
        else (
            "Panel mutations are a mix of clonal and subclonal variants — subclonal "
            "targets may be absent from metastatic or recurrent disease"
            if mean_ccf >= 0.5
            else "Many panel mutations appear subclonal, which may limit their "
            "reliability as MRD markers across tumor sites"
        )
    )
)

st.markdown(
    f"**Plain-language summary — {sample}:** "
    f"This panel tracks {n_panel} tumor-specific DNA changes in a blood test. "
    f"{clonal_desc}. "
    f"{sensitivity_plain}"
)

st.markdown("---")

# ── Caveats ───────────────────────────────────────────────────────────────────
with st.expander("⚠️ Caveats & Limitations", expanded=True):
    st.markdown(
        """
- **Diploid CCF assumption.** The formula CCF = 2 × VAF / purity assumes two
  chromosomal copies at every locus. In ccRCC, 3p loss makes tumor cells
  *hemizygous* at the *VHL* locus — this is why we use the corrected purity formula.
  Other loci may be amplified or carry additional LOH (loss of heterozygosity), so
  CCF values are approximations. Values are capped at 1.0 but the underlying
  biology can produce apparent CCF > 1 at amplified sites.

- **3p genes are also hemizygous.** *PBRM1*, *SETD2* and *BAP1* sit on 3p, like
  *VHL*. After 3p loss, their mutations behave like *VHL* (VAF = p/(2−p)), so the
  diploid CCF formula overestimates their CCF.

- **CCF cap ties.** Clonal mutations cluster around CCF = 1 at any purity, and
  sampling noise pushes about half of them above 1, so many tie at the cap and
  ranking among them is arbitrary. Two things make this worse: low purity (small
  VAFs have larger relative noise) and an underestimated purity (for example from
  a subclonal *VHL* mutation), which inflates every CCF. Check the *VHL* VAF and
  consider overriding purity with the slider.

- **Single-biopsy clonal inference.**
  [TRACERx Renal](https://doi.org/10.1016/j.cell.2018.03.043) demonstrated
  substantial inter- and intra-tumor heterogeneity in ccRCC. A single biopsy may
  miss subclonal mutations or mis-estimate clonal architecture; the strategy of
  selecting top-CCF mutations may not generalise across all tumor regions or
  metastatic sites.

- **Low ctDNA shedding in ccRCC.** Renal cell carcinoma is among the
  lowest ctDNA (circulating tumor DNA)-shedding tumor types. Achieving the
  plasma tumor fractions shown on the left of the plot may require very deep
  sequencing or may be practically undetectable even for clonal mutations,
  particularly in localized disease.

- **No sequencing-error model.** The Poisson model treats any mutant molecule as
  signal, which assumes perfect specificity. In practice, raw sequencing and PCR
  errors occur at rates comparable to or above the TFs on this plot, so
  error-suppression strategies such as UMI (Unique Molecular Identifier) or duplex
  sequencing are needed for low-TF detection, and are not modeled here.
        """
    )
