# Tumor-Informed MRD Panel Designer — TCGA KIRC

## Purpose

This is a Streamlit teaching app that designs a tumor-informed Minimal Residual
Disease (MRD) ctDNA panel from somatic mutations in clear cell renal cell carcinoma
(ccRCC). Given a TCGA KIRC tumor sample, the app infers tumor purity from the *VHL*
variant allele fraction (exploiting 3p loss, which occurs in >90% of ccRCC), computes
Cancer Cell Fractions (CCF) for each somatic mutation, selects the top-N clonal
mutations as panel targets, and models ctDNA detection probability as a function of
plasma tumor fraction (TF) using a CCF-weighted Poisson model. The goal is to make
the key quantitative steps in tumor-informed MRD panel design — purity inference,
clonality estimation, and detection modeling — transparent and interactive for
educational purposes.

---

## Data

**Source:** TCGA KIRC PanCancer Atlas 2018 (`kirc_tcga_pan_can_atlas_2018`),
downloaded from the [cBioPortal Datahub](https://github.com/cBioPortal/datahub).

The cohort contains 512 tumor samples; 356 have whole-exome sequencing (WXS) mutation
calls in this MAF. The remaining 156 samples are in the clinical file but not present
in this MAF (reason not verified).

The `data/` directory is **gitignored**. Use `get_data.sh` to download the two
required files:

```bash
bash get_data.sh
```

This requires `git-lfs`:

```bash
# macOS
brew install git-lfs

# Ubuntu/Debian
sudo apt install git-lfs

# RHEL/Amazon Linux
sudo yum install git-lfs
```

The script downloads `data_mutations.txt` (~31 MB) and `data_clinical_sample.txt`
(~100 KB) from the cBioPortal datahub git-LFS store into `data/`.

---

## Install

```bash
pip install -r requirements.txt
```

Or manually:

```bash
pip install streamlit==1.50.0 pandas==2.3.3 numpy==2.0.2 scipy==1.13.1 plotly==7.1.0
```

---

## Run

```bash
streamlit run app.py --server.port 8081
```

Then open http://localhost:8081 in your browser.

---

## Methods

### VAF

Variant Allele Fraction is computed per mutation as:

```
VAF = t_alt_count / t_depth
```

where `t_alt_count` and `t_depth` are columns from the MAF file. Sites with
`t_depth < 30` are excluded as low-confidence.

### Purity estimation

In ccRCC, *VHL* is inactivated in most tumors: it is mutated in roughly half of TCGA
cases, and 3p loss of heterozygosity (LOH) occurs in >90%. Both events are typically
truncal. When a *VHL* mutation sits on a 3p-loss background, tumor cells are
hemizygous at the locus: one copy is lost, and the remaining copy carries the
mutation.
For a hemizygous clonal mutation, the expected VAF in a mixed biopsy is:

```
VAF = p / (2 - p)
```

Rearranging to solve for purity *p*:

```
p = 2 * VAF / (1 + VAF)
```

If the sample carries a *VHL* mutation after filtering, purity is estimated from the
median *VHL* VAF using this formula. Otherwise the app defaults to a slider value
(default 0.6). The user can override purity at any time via the sidebar slider.

Note: the naive diploid formula `p = 2 * VAF` overestimates purity when one allele
is lost, because it ignores the denominator change from LOH.

### CCF

Cancer Cell Fraction is estimated using the standard diploid correction:

```
CCF = min(2 * VAF / purity, 1.0)
```

A CCF near 1.0 indicates a clonal mutation (present in nearly all tumor cells).
Values are capped at 1.0; apparent CCF > 1 before capping reflects copy-number
variation or purity underestimation. When all panel mutations hit the CCF = 1.0
cap simultaneously, ranking is uninformative and the panel selection is essentially
arbitrary. Clonal mutations cluster at CCF ≈ 1 regardless of purity, and sampling noise
pushes about half of them above the cap; low or underestimated purity makes this worse.

CCF ties are broken by sequencing depth (`t_depth`, descending) to prefer
better-evidenced mutations.

### CHIP exclusion

Mutations in common CHIP (Clonal Hematopoiesis of Indeterminate Potential) genes
(*DNMT3A*, *TET2*, *ASXL1*, *PPM1D*, *JAK2*) are excluded before panel selection.
Most plasma cfDNA comes from blood cells, so CHIP mutations appear in plasma without
any tumor present and would produce false-positive MRD calls. This gene-level filter
is a simplification: production assays sequence the patient's own white blood cells
and remove any variant found there.

### Poisson detection model

The expected number of mutant cfDNA molecules observed across the panel is:

```
λ = Σᵢ (D × TF × CCFᵢ / 2)
  = D × TF × Σ CCFᵢ / 2
```

where:
- **D** = unique cfDNA molecules per locus (genome equivalents, GE) — set by plasma
  input volume, not sequencing read depth. A standard 10 mL blood draw yields ~4–5 mL
  plasma, providing roughly 3,000–10,000 GE/locus after extraction and capture losses.
- **TF** = plasma tumor fraction (proportion of cfDNA derived from tumor cells)
- **CCFᵢ / 2** = expected mutant allele fraction for mutation *i* in plasma
  (diploid assumption: one mutant copy out of two)

Detection probability:

```
P(detect) = P(Poisson(λ) ≥ k)  =  1 - CDF(k-1, λ)
```

where **k** is the minimum number of mutant molecules required to call a sample
positive (default k = 2).

### Limit of Detection (LoD)

LoD50 and LoD95 are the plasma tumor fractions at which P(detect) = 50% and 95%,
respectively, computed by binary search over TF. LoD95 is the clinically relevant
threshold: below it, 5% or more of truly positive samples would be missed.

---

## Assumptions

1. **Diploid genome at all loci.** The CCF formula assumes two chromosomal copies
   everywhere. In reality, 3p loss makes tumor cells hemizygous at *VHL* (handled by
   the purity formula) but other loci may have amplifications or additional LOH.

2. **3p LOH is universal and clonal.** The *VHL* purity formula assumes the sample
   has undergone 3p LOH. If *VHL* is mutated on a diploid background (rare in ccRCC),
   purity will be underestimated.

3. **Single-biopsy clonal inference.** CCF is estimated from a single tumor biopsy.
   TRACERx Renal demonstrated substantial inter- and intra-tumor heterogeneity in
   ccRCC; a single biopsy may miss subclonal mutations, or call a mutation clonal
   that is subclonal tumor-wide.

4. **Other 3p genes are hemizygous too.** *PBRM1*, *SETD2* and *BAP1* are on 3p;
   after 3p loss, the diploid CCF formula overestimates their CCF.

5. **D is limited by cfDNA input, not sequencing depth.** The model parameter D
   represents unique DNA molecules (GE), not sequencing reads. After UMI
   deduplication, additional reads from the same molecule add no information.

6. **No sequencing-error model (perfect specificity assumed).** The Poisson model
   treats all mutant molecules as true signal. In practice, raw sequencing and PCR
   error rates are comparable to or higher than the TFs modeled, so UMI or duplex
   error suppression is needed for low-TF detection.

---

## Limitations

- Teaching demo only; not validated for clinical use.
- Uses TCGA WXS data, which has lower sensitivity than targeted deep sequencing.
- No copy-number correction for CCF.
- No multi-region sampling; clonal architecture is inferred from a single biopsy.
- ccRCC is among the lowest ctDNA-shedding tumor types; the TF values on the left
  of the detection plot may be practically unachievable without large plasma volumes.

---

## Screenshots

Manual screenshots can be saved to `docs/screenshots/`. No automated screenshot
tooling is included.

---

## Project structure

```
mrd_panel/
├── app.py                  # Streamlit app (single file, ~600 lines)
├── get_data.sh             # Download script for MAF and clinical data
├── requirements.txt        # Pinned Python dependencies
├── README.md               # This file
├── DEVLOG.md               # Design decisions and review history
├── .gitignore
└── data/                   # Gitignored — populated by get_data.sh
    ├── data_mutations.txt
    └── data_clinical_sample.txt
```
