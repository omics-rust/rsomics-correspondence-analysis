# rsomics-correspondence-analysis

**Correspondence Analysis (CA)** of a non-negative feature/count table —
the chi-square-distance ordination preferred over PCA when the data has many
zeros (steep or long gradients).

Reads a feature/count table TSV — an empty top-left cell, the feature IDs in the
header row, then one row per sample (sample ID followed by tab-separated counts)
— and writes the eigenvalues, proportion explained, sample (row) scores and
feature (column) scores on the CA axes.

```
rsomics-correspondence-analysis otu_table.tsv
rsomics-correspondence-analysis table.tsv --n-axes 2 -o ca.tsv
```

## Method

Matches `skbio.stats.ordination.ca` with `scaling=1` (preserves chi-square
distances between samples), following Legendre & Legendre 1998 §9.4.1:

1. `Q = X / grand_total`; row marginals `p_i+`, column marginals `p_+j`.
2. The chi-square-transformed matrix (Eq. 9.32):
   `Q̄ = (Q − p_i+·p_+j) / √(p_i+·p_+j)`.
3. SVD of `Q̄ = Û · W · Uᵀ`. Centering leaves at most `min(r,c) − 1` non-zero
   singular values; the rank is taken at numpy's `matrix_rank` tolerance.
4. Eigenvalues `λ = W²`; `proportion_explained = λ / Σλ`.
5. `V = D(p_+j)^-½ · U`, `V̂ = D(p_i+)^-½ · Û`, `F = V̂ · W`.
   With scaling 1, **sample scores = F**, **feature scores = V**.

`--n-axes N` keeps the first `N` axes (default: full rank).

### Eigenvector sign

The sign of a singular vector — and therefore of a CA axis — is arbitrary;
flipping a whole axis is an equally valid solution. Eigenvalues and proportion
explained are sign-independent. Output is given a deterministic orientation (the
largest-magnitude feature loading on each axis is made positive); the compat
differential additionally compares scores up to a per-axis sign flip.

### Output

A flat TSV: an `# eigenvalues` block (eigenvalue + proportion_explained rows over
`CA1 … CAk`), a `# samples` block (one row per sample), and a `# features` block
(one row per feature). Floats use Python's shortest round-trip `repr`.

## Origin

This crate is an independent Rust reimplementation of the correspondence-analysis
operation provided by `scikit-bio` (`skbio.stats.ordination.ca`, which delegates
the SVD to `scipy.linalg.svd`), based on:

- Legendre, P. & Legendre, L. (1998), *Numerical Ecology* (2nd ed.), Elsevier,
  §9.4 (the row/column-profile chi-square ordination, Eq. 9.32–9.44),
  <https://doi.org/10.1016/B978-0-444-53868-0.50009-5>.
- Greenacre, M. (1984), *Theory and Applications of Correspondence Analysis*,
  Academic Press; ter Braak, C. J. F. (1985), *Correspondence analysis of
  incidence and abundance data*, Biometrics 41:859-873.
- The black-box behaviour of `skbio.stats.ordination.ca`: the chi-square
  transform, the `min(r,c) − 1` rank cap at the `matrix_rank` tolerance, and the
  scaling-1 choice of `F`/`V` for the sample/feature scores (Eqs. 9.43a / 9.44a,
  which match vegan).

scikit-bio is BSD-3-Clause and was read and cited. The SVD uses
[`faer`](https://crates.io/crates/faer) (pure Rust, SIMD + rayon —
external-dependency quadrant ①). Test fixtures are deterministically generated
count tables.

License: MIT OR Apache-2.0.
Upstream credit: scikit-bio <https://scikit-bio.org> (BSD-3-Clause).

## Compatibility & performance

`tests/compat.rs` runs this binary and the scikit-bio oracle
(`tests/oracle_skbio.py`) and asserts the eigenvalues (`1e-9`) and proportion
explained match directly, and the sample/feature scores match up to a per-axis
sign flip (`1e-6`). It checks a committed skbio-captured golden (always runs) and
a live skbio differential (skipped loudly when scikit-bio is not importable).

On a 1000 × 1200 count table, both pinned to one core, this binary runs the CA
**3.31× faster** than scikit-bio's `ca()` (2.77 s vs 9.18 s wall; single-thread
BLAS on the upstream side), with eigenvalues agreeing to ~14 significant digits.
The hot path is the dense `O(min(r,c)²·max(r,c))` SVD; `faer` carries it.
