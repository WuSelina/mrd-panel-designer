# DEVLOG — Tumor-Informed MRD Panel Designer

Design decisions, prompt history, and corrections made during development.

---

## Session overview

The app was built interactively in a single session via an AI-assisted coding
workflow. The prompts below are paraphrased; the key design decisions and corrections
that followed each review round are documented in detail.

---

## Prompts and iterations

### Prompt 1 — Data download

> "Please download the cBioPortal TCGA KIRC PanCancer Atlas study from the cBioPortal
> datahub. Extract only data_mutations.txt and data_clinical_sample.txt into a data/
> folder. Tell me the columns available for read counts and purity."

**Outcome:** The cBioPortal S3 bucket returned 403; the GitHub raw URL serves git-LFS
pointer files (131 bytes), not real data. Resolved by installing `git-lfs`, cloning
with `--filter=blob:none --no-checkout`, configuring sparse checkout, then calling
the git-LFS batch API directly to get signed S3 download URLs for the two files.

**Key finding:** `data_clinical_sample.txt` contains no purity column. Tumor purity
must be inferred from the MAF itself (via *VHL* VAF) or set manually.

---

### Prompt 2 — Initial app plan

> "Build a Streamlit app (port 8081) that designs a tumor-informed MRD panel..."

Specified: VAF computation, VHL purity inference, CCF, CHIP exclusion, top-N panel,
Poisson detection model with sliders for D/k/TF, panel table with TSV download, and
a caveats box.

**Planning decisions:**
- Single-file app (`app.py`), `@st.cache_data` for the MAF load.
- Sidebar for all controls; main area for output.
- Use `plotly` for the detection curve (log x-axis, interactive).

---

### Prompt 3 — Revised design before implementation

> "The plan looks good with these changes: VHL purity using p = 2·VAF/(1+VAF),
> user-adjustable reference TF line, CCF-weighted λ, and ⓘ popovers."

**Key corrections at this stage:**

1. **VHL purity formula.** The original spec used "purity ≈ VAF", which is only
   exact at 100% purity and underestimates purity otherwise (VAF 0.43 actually
   means purity 0.60). Corrected to the hemizygous formula: since tumor cells at the
   *VHL* locus are hemizygous (3p LOH), VAF = p/(2−p), so p = 2·VAF/(1+VAF). The
   naive diploid formula p = 2·VAF would instead overestimate purity.

2. **CCF-weighted detection model.** Original spec used λ = N × D × TF / 2 (uniform
   CCF = 1). Corrected to λ = Σᵢ (D × TF × CCFᵢ / 2), which weights each mutation
   by its clonal fraction. This is a meaningful difference for panels with subclonal
   mutations.

3. **Reference TF line.** Original spec placed a fixed dashed line at TF = 0.01%
   labelled as a "ccRCC floor." Removed the biological label (no published source was
   cited) and made it a user-adjustable slider instead.

4. **ⓘ popovers.** Added for VAF, purity, CHIP, CCF, and Poisson model using
   `st.popover()` (available in Streamlit ≥ 1.31).

---

### Prompt 4 — Implementation (Tasks 1–6)

All six tasks implemented in one pass. The app loaded 29,473 mutations from 356
samples. Syntax and logic verified with `python3 -c "import ast; ast.parse(...)"` and
manual numerical checks.

**Post-implementation corrections:**

- **jinja2 error:** `panel.style.format(...)` requires `jinja2`, which was not
  installed. Replaced with direct column formatting:
  `display["VAF"] = display["VAF"].map("{:.3f}".format)`.

---

### Prompt 5 — Depth parameter correction

> "69% detection at 0.001% TF is far too optimistic. D isn't the number of reads you
> sequence, it's the number of unique DNA molecules in the blood sample."

**Key correction — genome equivalents vs. sequencing depth:**

The original default D = 30,000 was labelled "depth per locus" and interpreted as
sequencing reads. This is wrong: after UMI deduplication, additional reads from the
same molecule add no information. D is set by plasma input volume (cfDNA molecules
available), not sequencing depth. A standard 10 mL blood draw yields ~4–5 mL plasma,
providing ~3,000–10,000 GE/locus after extraction and capture losses.

**Changes:** D relabelled "Unique molecules per locus (genome equivalents, GE)";
default changed from 30,000 to 5,000; range changed from 1,000–200,000 to 500–50,000.
The GE popover also corrected a factual error: the draft said "10 mL blood yields
~30 mL plasma" — corrected to "~4–5 mL plasma."

With D = 5,000 and N = 16 fully clonal mutations at TF = 0.001%, λ = 0.4 → P(≥2) ≈
6% (realistic). The previous D = 30,000 gave λ = 2.4 → P(≥2) ≈ 69% (inflated).

**Additional correction — CCF cap and low purity:**

A new caveat was added: at low purity (< ~0.3), the CCF formula 2×VAF/purity
amplifies noise and pushes many mutations to the cap of 1.0 simultaneously. This
makes CCF ranking uninformative. The caveat notes that this may reflect either a
genuinely low-purity biopsy or a subclonal *VHL* mutation that underestimated purity.

---

### Prompt 6 — Feature additions

> "Add variant classification multiselect, CCF tie-breaking by t_depth, LoD50/LoD95
> metrics with ⓘ, explicit x-axis tick labels, summary section, default ref TF 0.01%,
> expanded D label, WXS not WGS."

**Changes:**
- `st.multiselect` for `Variant_Classification` above the panel table, with popover
  noting that silent mutations are valid MRD markers (DNA tracking, not function) but
  SNVs are more reliable than indels.
- Tie-breaking: replaced `nlargest()` with `sort_values(["CCF","t_depth"], ascending=[False,False]).head(n)`.
- LoD50/LoD95 computed by 60-iteration binary search; displayed as diamond markers on
  the detection curve and as `st.metric` values.
- X-axis ticks explicitly set to `["0.001%", "0.01%", "0.1%", "1%"]`.
- Reference TF default changed from 0.1% to 0.01%.
- Caption updated: "356 of 512 cohort samples have WXS mutation calls in this MAF"
  (156 samples are in the clinical file but not present in the MAF; reason not
  verified).

---

### Prompt 7 — Python 3.9 compatibility

> "TypeError: unsupported operand type(s) for |: 'type' and 'NoneType'"

**Fix:** Python 3.9 does not support `float | None` union syntax in type hints.
Replaced with `Optional[float]` from `typing`.

---

### Prompt 8 — Spelling, gene nomenclature, summary style

> "I use 'tumor' not 'tumour'. Gene names should follow proper nomenclature.
> Technical summary too long; both summaries should show sample ID up front."

**Changes:**
- All British spellings replaced globally (`tumour→tumor`, `modelled→modeled`,
  `optimised→optimized`, `localised→localized`).
- Human gene symbols italicised in all prose markdown: *VHL*, *DNMT3A*, *TET2*,
  *ASXL1*, *PPM1D*, *JAK2* (HGNC convention).
- Technical summary rewritten from three verbose sentences to one compact line:
  `purity; panel size; mean CCF (N at cap). LoD50, LoD95.`
- Both summary headers now bold with colon and sample ID: **Technical summary —
  SAMPLE:** and **Plain-language summary — SAMPLE:**
- `.capitalize()` bug fixed: it was lowercasing `*VHL* VAF` to `*vhl* vaf`, breaking
  both the gene case and the markdown italics. Replaced with `s[0].upper() + s[1:]`.
- CCF-cap warning added inline to technical summary when all N mutations hit the cap.

---

### Prompt 9 — Panel size context

> "For this sample, it's not clear why there are only 16 mutations reported."

**Fix:** Subheader changed from "Top 16 Panel Mutations" to "Top 16 of 87 Panel
Mutations" with a caption: "Showing the top N mutations by CCF... Adjust Panel size
(N) in the sidebar."

---

### Prompt 10 — Plain-language summary ambiguity

> "'This level' in 'if a patient has this level of tumor DNA' is unclear."

**Fix:** Restructured `sensitivity_plain` to be a self-contained sentence leading
with the explicit LoD95 value: "At its modeled sensitivity limit of 0.012% tumor DNA
in plasma, the panel would detect a positive sample in 95 out of 100 tests."

---

### Prompt 11 — CCF cap description in plain-language summary

> "For this sample aren't all 16/16 genes CCF = 1? The text should say 'all' not 'most'."

**Fix:** Added a fourth case to `clonal_desc`: when `n_at_cap == n_panel`, the text
reads "All N panel mutations are at the CCF = 1.0 cap — they are either truly clonal
or the purity estimate is too low to distinguish clonal from subclonal (see caveats)."
This is more honest than "most appear clonal" and explicitly flags the artefact.

---

## Summary of key design decisions

| Decision | Rationale |
|---|---|
| VHL purity via p = 2·VAF/(1+VAF) | 3p LOH makes tumor cells hemizygous; naive p = 2·VAF overestimates |
| D = genome equivalents, not reads | Sequencing depth beyond molecule count adds no information post-UMI dedup |
| D default = 5,000 GE | Realistic for 10 mL blood draw; 30,000 reads was misleadingly optimistic |
| CCF-weighted λ | Mutations with CCF < 1.0 contribute proportionally less signal |
| Tie-break by t_depth | Prefer better-evidenced mutations when CCF is tied at the cap |
| Reference TF line is user-adjustable | No defensible published ccRCC-specific floor; don't assert one |
| LoD via binary search | Direct, exact, no closed-form needed |
| CHIP exclusion | DNMT3A/TET2/ASXL1/PPM1D/JAK2 mutations originate from blood, not tumor |
| Silent mutations included | MRD tracks DNA sequence, not protein function |
| 356 of 512 samples have WXS | 156 absent from MAF; reason not verified |

---

### Review round — post-workshop code and documentation review

- **get_data.sh bug:** `DATADIR` was a relative path, but the script `cd`s into the
  temp clone before copying, so running `bash get_data.sh` wrote the data into the
  temp folder and then deleted it. Fixed with an absolute path, a cleanup `trap`, a
  skip-if-present check, and a check that the downloaded file is not an LFS pointer.
- **CCF cap wording:** the earlier caveat blamed cap ties on low purity alone. Clonal
  mutations cluster at CCF ≈ 1 at any purity and noise pushes about half above the
  cap; low or underestimated purity only makes it worse. Text and popover corrected
  (the popover had said the opposite, "high purity").
- **3p genes:** added a caveat that *PBRM1*, *SETD2* and *BAP1* are also on 3p, so
  the diploid CCF formula overestimates their CCF.
- **CHIP rationale:** clarified that CHIP is a plasma false-positive problem (most
  cfDNA comes from blood cells) and that production assays remove CHIP by
  sequencing matched white blood cells, not by gene lists.
- **Unsourced numbers:** removed specific per-base error-rate and "below ~0.1% TF"
  thresholds; reworded VHL frequency claims (mutated in ~half of TCGA ccRCC; 3p loss
  in >90%) and the MC3 exclusion guess.
- **Minor:** P(detect) at the reference TF is now computed exactly instead of from
  the nearest plot grid point; TSV export labels the index column `rank`.
