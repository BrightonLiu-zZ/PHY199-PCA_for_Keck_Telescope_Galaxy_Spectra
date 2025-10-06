# PHY199 — DEIMOS Spectral PCA & Line-Ratio Analysis

This repo contains the full, reproducible pipeline I used in my PHY199 research project to go from **raw DEIMOS spectra** → **PCA analysis & filtering** → **paper-ready figures** (ggplot in R).
The final figure bundle is **`plots/plot_collection.pdf`**.

---

## Repo layout (current)

```
.
├─ .gitattributes
├─ environment.yml                  # Python env (conda/miniforge)
├─ data/
│  ├─ list_file_9.txt              # list of usable spectra with λ_min/λ_max
│  └─ pc_scores_9.csv              # PCA scores + (line ratios if available)
├─ notebooks/
│  ├─ main_jupyter_notebook.ipynb  # PCA EDA, filtering, elbow, 3D proj, cornerplot
│  └─ plot_collection.Rmd          # optional R Markdown gallery
├─ plots/
│  ├─ plot_collection.pdf          # final bundle of 5 plots (ggplot)
│  └─ *.png                        # example anomaly spectra
├─ scripts/
│  ├─ get_flux_matrix.py           # build/align flux matrix for PCA
│  ├─ get_list_file.py             # make list_file_9.txt from FITS headers
│  ├─ histogram_lambda_min.py      # explore λ_min distribution
│  ├─ histogram_lambda_max.py      # explore λ_max distribution
│  ├─ pc_reconstruction.py         # reconstruct spectra with first n PCs
│  └─ plotting.R                   # generates plot_collection.pdf (ggplot)
└─ src/
   └─ pyzeutil.py                  # pyze helper (incl. unsharp-mask continuum)
```

> **Data note:** The large **flux matrix** is hosted on Google Drive (not in the repo).
> You can (a) **rebuild** it locally from FITS using `scripts/get_flux_matrix.py`, or (b) **download** the precomputed matrix from Drive and place it under `data/` (see Step 3 below).

---

## What the pipeline produces

**`plots/plot_collection.pdf`** contains:

1. **Stacked PC1, PC2, PC3** (eigenspectra)
2. **PC3 vs log10([O II]/[O III])**
3. **Linear combo of standardized PC1–3 vs F(Hβ)**
4. **Linear combo of standardized PC1–2 vs F(Hβ)**
5. **Example spectra (stacked, Gaussian-smoothed)** near representative Hβ fluxes (~25, 400, 850)

It also saves example anomaly spectra as individual PNGs in `plots/`.

---

## Quickstart (one-time setup)

### 0) Python environment (conda / miniforge)

```bat
# From repo root:
conda env create -f environment.yml
conda activate phy199
# (optional) make Jupyter kernel visible
python -m ipykernel install --user --name phy199 --display-name "Python (phy199)"
```

### 1) R packages for plotting

Open R (or Rterm) in the repo root and run:

```r
install.packages(c("ggplot2", "readr", "dplyr", "scales", "patchwork", "rmarkdown"))
```

(If you prefer pinned R deps, initialize **renv** later with `renv::init()` and commit `renv.lock`.)

---

## End-to-end pipeline (detailed)

> You can run only the steps you need. If you already have `data/pc_scores_9.csv`, jump to **Step 5**.

### Step 1 — Build a list of usable spectra

Create `data/list_file_9.txt` by scanning your FITS (skipping ones without `PYZELMIN/PYZELMAX`):

```bat
conda activate phy199
python scripts/get_list_file.py
```

* **Input:** a folder or zip of your DEIMOS **spec1d** FITS files (configure the path at the top of `get_list_file.py` if needed).
* **Output:** `data/list_file_9.txt` with lines like:

  ```
  Spec1D File: <name>, Redshift: <z>, Best Matching Template: <tmpl>, lambda_min: <…>, lambda_max: <…>
  ```
* Files with missing/NaN `lambda_min` or `lambda_max` are **skipped**.

### Step 2 — Choose a wavelength window (sanity check)

Inspect the wavelength coverage to pick a safe common range for alignment (e.g., 4900–7200 Å if that matches your dataset; adjust to your data):

```bat
python scripts/histogram_lambda_min.py
python scripts/histogram_lambda_max.py
```

* These scripts display/emit histograms so you can decide the **[λ_min, λ_max]** you’ll use in Step 3.

### Step 3 — Build the flux matrix for PCA

Construct a uniformly sampled flux matrix on the chosen wavelength grid and apply your **continuum subtraction / unsharp-mask** (implemented via `src/pyzeutil.py`):

```bat
python scripts/get_flux_matrix.py
```

* **Configure inside the script** (top constants) if you need to change:

  * data folder / list file path
  * wavelength grid (start, end, step; e.g., 0.3 Å sampling)
  * preprocessing toggles (smoothing, unsharp mask)
* **Output (typical):**

  * `data/flux_matrix.npz` (or similar) — shape `(N_spectra, N_pixels)`
  * optional QA CSVs/PNGs depending on your script settings

> **Alternative:** If you prefer the precomputed matrix, download it from Google Drive and place it under `data/` with the expected file name (update the notebook/script paths accordingly).

### Step 4 — PCA, filtering, and exploration (Jupyter)

Open **`notebooks/main_jupyter_notebook.ipynb`** (select the **Python (phy199)** kernel):

This notebook performs:

* **Elbow plot** of explained variance
* **3D PCA projection** (PC1, PC2, PC3)
* **PCA-based filtering algorithm** (remove artifacts/outliers by PC space rules)
* **Corner plot** of the first 3 PCs

It then **exports** the per-spectrum PCA results:

* **Output:** `data/pc_scores_9.csv` (columns should include at least: `obj_id`/filename, `PC1`, `PC2`, `PC3`; and—if available—line measurements such as `log10_OII_OIII` and `Hbeta_flux` used in Step 5).

> If your line measurements are produced elsewhere, merge them into `pc_scores_9.csv` before Step 5 (keys by file name or object ID).

### Step 5 — Make the paper-ready figure bundle (R, ggplot)

You can create **`plots/plot_collection.pdf`** either via the **R script** or the **R Markdown** gallery.

**Option A — R script (batchable, recommended):**

```bat
Rscript scripts/plotting.R
```

* This reads `data/pc_scores_9.csv` and writes `plots/plot_collection.pdf`.
* If your script expects arguments, run e.g.:

  ```bat
  Rscript scripts/plotting.R data/pc_scores_9.csv plots/plot_collection.pdf
  ```

  (Update the top of `scripts/plotting.R` accordingly.)

**Option B — R Markdown gallery (for previewing):**

```r
rmarkdown::render("notebooks/plot_collection.Rmd",
                  output_file = "../plots/plot_collection.pdf")
```

**What gets plotted (by `plotting.R`):**

* Stacked **PC1/PC2/PC3** eigenspectra (standardized)
* **PC3 vs log10([O II]/[O III])**
* **(α·PC1 + β·PC2 + γ·PC3)** vs **F(Hβ)** (standardized PCs; coefficients defined in the script)
* **(α·PC1 + β·PC2)** vs **F(Hβ)**
* **Example spectra** (Gaussian-smoothed) at representative Hβ levels

### (Optional) Step 6 — Reconstruct spectra with the first *n* PCs

Check how well a small number of PCs capture spectral features:

```bat
python scripts/pc_reconstruction.py
```

* Produces diagnostic plots (e.g., original vs reconstruction at n=3/5/10 PCs).

---

## Re-run the whole thing (TL;DR)

```bat
# 0) one-time setup
conda env create -f environment.yml
conda activate phy199
python -m ipykernel install --user --name phy199 --display-name "Python (phy199)"
R -q -e "install.packages(c('ggplot2','readr','dplyr','scales','patchwork','rmarkdown'))"

# 1) list usable spectra
python scripts/get_list_file.py

# 2) inspect wavelength coverage
python scripts/histogram_lambda_min.py
python scripts/histogram_lambda_max.py
# (adjust wavelength settings in get_flux_matrix.py if needed)

# 3) build flux matrix
python scripts/get_flux_matrix.py

# 4) run PCA & export scores in the notebook
# (open notebooks/main_jupyter_notebook.ipynb → run all)

# 5) make final figure bundle
Rscript scripts/plotting.R
# or:
# R -q -e "rmarkdown::render('notebooks/plot_collection.Rmd', output_file='../plots/plot_collection.pdf')"
```

---

## Tips & troubleshooting

* **Conda env not showing in Jupyter?**
  Run: `python -m ipykernel install --user --name phy199 --display-name "Python (phy199)"`, then re-open the notebook.

* **Missing/NaN wavelength headers:**
  `get_list_file.py` **skips** spectra without both `PYZELMIN` and `PYZELMAX`. Regenerate `data/list_file_9.txt` after any fixes.

* **Large data:**
  The flux matrix is intentionally not committed. Keep large artifacts under `data/` (git-ignored) or download from Google Drive as needed.

* **Reproducible R deps:**
  Consider `renv::init()` and `renv::snapshot()` to pin ggplot/tidyverse versions.

---

## License & citation

* **License:** MIT (or choose another license to match collaborators’ needs).
* **Citation:** If you use this code/figures, please cite the upcoming RNAAS note (citation to be added here once available).

---

## Contact

Questions or bugs: open a GitHub Issue or email me. Happy to help others reuse the PCA/line-ratio workflow.
