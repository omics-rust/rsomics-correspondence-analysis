use std::io::{BufRead, Write};

use faer::Mat;
use faer::linalg::solvers::Svd;
use rsomics_common::{Result, RsomicsError};

mod fmt;
use fmt::push_pyrepr;

/// A non-negative feature/count table: samples are rows, features are columns.
/// TSV form is an empty top-left cell then feature IDs as the header, then one
/// row per sample (sample ID + tab-separated counts).
pub struct FeatureTable {
    pub sample_ids: Vec<String>,
    pub feature_ids: Vec<String>,
    /// Row-major `n_samples × n_features`.
    pub data: Vec<f64>,
}

impl FeatureTable {
    /// # Errors
    /// Errors on a missing header, a ragged body, a non-numeric cell, or a
    /// negative value (CA's chi-square transform is undefined for those).
    pub fn parse<R: BufRead>(reader: R, delim: char) -> Result<FeatureTable> {
        let mut lines = reader.lines();
        let header = loop {
            match lines.next() {
                Some(line) => {
                    let line = line.map_err(RsomicsError::Io)?;
                    if line.trim().is_empty() || line.starts_with('#') {
                        continue;
                    }
                    break line;
                }
                None => return Err(RsomicsError::InvalidInput("empty feature table".into())),
            }
        };
        let feature_ids: Vec<String> = header
            .split(delim)
            .skip(1)
            .map(|s| s.trim().to_string())
            .collect();
        let p = feature_ids.len();
        if p == 0 {
            return Err(RsomicsError::InvalidInput(
                "header has no feature columns (need an empty top-left cell + ≥1 feature)".into(),
            ));
        }

        let mut sample_ids = Vec::new();
        let mut data = Vec::new();
        for line in lines {
            let line = line.map_err(RsomicsError::Io)?;
            if line.trim().is_empty() || line.starts_with('#') {
                continue;
            }
            let mut fields = line.split(delim);
            let label = fields.next().unwrap_or("").trim().to_string();
            let row_start = data.len();
            for field in fields {
                let v: f64 = field.trim().parse().map_err(|_| {
                    RsomicsError::InvalidInput(format!(
                        "sample '{label}', column {}: '{}' is not numeric",
                        data.len() - row_start + 1,
                        field.trim()
                    ))
                })?;
                if v < 0.0 {
                    return Err(RsomicsError::InvalidInput(
                        "input matrix elements must be non-negative".into(),
                    ));
                }
                data.push(v);
            }
            let got = data.len() - row_start;
            if got != p {
                return Err(RsomicsError::InvalidInput(format!(
                    "sample '{label}' has {got} values, expected {p}"
                )));
            }
            sample_ids.push(label);
        }
        if sample_ids.is_empty() {
            return Err(RsomicsError::InvalidInput("no data rows".into()));
        }
        Ok(FeatureTable {
            sample_ids,
            feature_ids,
            data,
        })
    }

    #[must_use]
    pub fn n_samples(&self) -> usize {
        self.sample_ids.len()
    }

    #[must_use]
    pub fn n_features(&self) -> usize {
        self.feature_ids.len()
    }
}

/// Result of a Correspondence Analysis: eigenvalues, proportion explained, and
/// the sample (row) and feature (column) scores on each CA axis. Matches
/// `skbio.stats.ordination.ca` with `scaling=1` (preserves chi-square distances
/// between samples): sample scores are F, feature scores are V.
pub struct Ordination {
    pub sample_ids: Vec<String>,
    pub feature_ids: Vec<String>,
    pub eigvals: Vec<f64>,
    pub proportion_explained: Vec<f64>,
    /// Row-major `n_samples × n_axes`.
    pub sample_scores: Vec<f64>,
    /// Row-major `n_features × n_axes`.
    pub feature_scores: Vec<f64>,
}

impl Ordination {
    /// CA via the chi-square transform (Legendre & Legendre 9.32) then SVD.
    /// `n_axes` caps the number of retained axes; `None` keeps the full rank.
    #[must_use]
    pub fn compute(table: &FeatureTable, n_axes: Option<usize>) -> Ordination {
        let r = table.n_samples();
        let c = table.n_features();

        let grand_total: f64 = table.data.iter().sum();
        let q: Vec<f64> = table.data.iter().map(|&x| x / grand_total).collect();

        let mut row_marginals = vec![0.0_f64; r];
        let mut col_marginals = vec![0.0_f64; c];
        for i in 0..r {
            for j in 0..c {
                let v = q[i * c + j];
                row_marginals[i] += v;
                col_marginals[j] += v;
            }
        }

        // Q_bar = (Q - E) / sqrt(E), E = row_marginal ⊗ col_marginal  (Eq. 9.32)
        let q_bar = Mat::from_fn(r, c, |i, j| {
            let e = row_marginals[i] * col_marginals[j];
            (q[i * c + j] - e) / e.sqrt()
        });

        let svd: Svd<f64> = q_bar.svd().unwrap();
        let s = svd.S().column_vector();
        let u_hat = svd.U();
        let v_right = svd.V();

        // Centering leaves at most min(r,c)-1 non-zero singular values.
        let k = s.nrows();
        let tol = {
            let smax = (0..k).fold(0.0_f64, |m, i| m.max(s[i]));
            smax * r.max(c) as f64 * f64::EPSILON
        };
        let rank = (0..k).filter(|&i| s[i] > tol).count();
        let n_keep = match n_axes {
            Some(n) => n.min(rank),
            None => rank,
        };

        let w: Vec<f64> = (0..n_keep).map(|a| s[a]).collect();

        // Deterministic sign per axis: largest-abs component of the right
        // singular vector positive (matches a stable orientation; the compat
        // test still sign-aligns since the eigensolver's sign is arbitrary).
        let sign: Vec<f64> = (0..n_keep)
            .map(|a| {
                let mut best = 0.0_f64;
                let mut sgn = 1.0_f64;
                for j in 0..c {
                    let val = v_right[(j, a)];
                    if val.abs() > best {
                        best = val.abs();
                        sgn = if val < 0.0 { -1.0 } else { 1.0 };
                    }
                }
                sgn
            })
            .collect();

        // V = D(col_marginal)^-1/2 · U,  V_hat = D(row_marginal)^-1/2 · U_hat
        // F = V_hat · W  (sample scores, scaling 1);  feature scores = V.
        let mut feature_scores = vec![0.0_f64; c * n_keep];
        for j in 0..c {
            let scale = col_marginals[j].powf(-0.5);
            for a in 0..n_keep {
                feature_scores[j * n_keep + a] = scale * v_right[(j, a)] * sign[a];
            }
        }
        let mut sample_scores = vec![0.0_f64; r * n_keep];
        for i in 0..r {
            let scale = row_marginals[i].powf(-0.5);
            for a in 0..n_keep {
                sample_scores[i * n_keep + a] = scale * u_hat[(i, a)] * sign[a] * w[a];
            }
        }

        let eigvals: Vec<f64> = w.iter().map(|&x| x * x).collect();
        let total: f64 = (0..rank).map(|a| s[a] * s[a]).sum();
        let proportion_explained: Vec<f64> = eigvals.iter().map(|&v| v / total).collect();

        Ordination {
            sample_ids: table.sample_ids.clone(),
            feature_ids: table.feature_ids.clone(),
            eigvals,
            proportion_explained,
            sample_scores,
            feature_scores,
        }
    }

    /// Write a flat ordination TSV: an `# eigenvalues` block, a sample-score
    /// table, then a feature-score table, axes labelled `CA1..CAk`.
    ///
    /// # Errors
    /// Propagates write errors.
    pub fn write_tsv<W: Write>(&self, mut out: W) -> Result<()> {
        let n_axes = self.eigvals.len();
        let mut line = String::new();

        writeln!(out, "# eigenvalues").map_err(RsomicsError::Io)?;
        write_axis_header(&mut out, n_axes)?;
        line.push_str("eigval");
        for &v in &self.eigvals {
            line.push('\t');
            push_pyrepr(&mut line, v);
        }
        writeln!(out, "{line}").map_err(RsomicsError::Io)?;

        line.clear();
        line.push_str("proportion_explained");
        for &v in &self.proportion_explained {
            line.push('\t');
            push_pyrepr(&mut line, v);
        }
        writeln!(out, "{line}").map_err(RsomicsError::Io)?;

        writeln!(out, "# samples").map_err(RsomicsError::Io)?;
        write_axis_header(&mut out, n_axes)?;
        for (i, id) in self.sample_ids.iter().enumerate() {
            line.clear();
            line.push_str(id);
            for a in 0..n_axes {
                line.push('\t');
                push_pyrepr(&mut line, self.sample_scores[i * n_axes + a]);
            }
            writeln!(out, "{line}").map_err(RsomicsError::Io)?;
        }

        writeln!(out, "# features").map_err(RsomicsError::Io)?;
        write_axis_header(&mut out, n_axes)?;
        for (j, id) in self.feature_ids.iter().enumerate() {
            line.clear();
            line.push_str(id);
            for a in 0..n_axes {
                line.push('\t');
                push_pyrepr(&mut line, self.feature_scores[j * n_axes + a]);
            }
            writeln!(out, "{line}").map_err(RsomicsError::Io)?;
        }
        Ok(())
    }
}

fn write_axis_header<W: Write>(out: &mut W, n_axes: usize) -> Result<()> {
    let mut header = String::new();
    for a in 1..=n_axes {
        header.push('\t');
        header.push_str("CA");
        header.push_str(&a.to_string());
    }
    writeln!(out, "{header}").map_err(RsomicsError::Io)
}

/// # Errors
/// Propagates parse and write errors.
pub fn run<R: BufRead, W: Write>(
    reader: R,
    out: W,
    delim: char,
    n_axes: Option<usize>,
) -> Result<()> {
    let table = FeatureTable::parse(reader, delim)?;
    let ord = Ordination::compute(&table, n_axes);
    ord.write_tsv(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ll_table() -> &'static str {
        "\tSp1\tSp2\tSp3\n\
         Site1\t10\t10\t20\n\
         Site2\t10\t15\t10\n\
         Site3\t15\t5\t5\n"
    }

    #[test]
    fn parses_table() {
        let t = FeatureTable::parse(ll_table().as_bytes(), '\t').unwrap();
        assert_eq!(t.sample_ids, ["Site1", "Site2", "Site3"]);
        assert_eq!(t.feature_ids, ["Sp1", "Sp2", "Sp3"]);
        assert_eq!(t.data[2 * 3], 15.0);
    }

    #[test]
    fn negative_value_errors() {
        let bad = "\tA\tB\nS1\t1\t-2\n";
        assert!(FeatureTable::parse(bad.as_bytes(), '\t').is_err());
    }

    #[test]
    fn ragged_row_errors() {
        let bad = "\tA\tB\nS1\t1\nS2\t1\t2\n";
        assert!(FeatureTable::parse(bad.as_bytes(), '\t').is_err());
    }

    /// L&L 1998 table 9.11; values verified against skbio ca().
    #[test]
    fn matches_known_eigenvalues() {
        let t = FeatureTable::parse(ll_table().as_bytes(), '\t').unwrap();
        let o = Ordination::compute(&t, None);
        assert_eq!(o.eigvals.len(), 2);
        assert!((o.eigvals[0] - 0.09613302).abs() < 1e-7, "{}", o.eigvals[0]);
        assert!((o.eigvals[1] - 0.04094181).abs() < 1e-7, "{}", o.eigvals[1]);
        let p: f64 = o.proportion_explained.iter().sum();
        assert!((p - 1.0).abs() < 1e-12);
        assert!((o.proportion_explained[0] - 0.70131778).abs() < 1e-7);
    }

    #[test]
    fn n_axes_caps_output() {
        let t = FeatureTable::parse(ll_table().as_bytes(), '\t').unwrap();
        let o = Ordination::compute(&t, Some(1));
        assert_eq!(o.eigvals.len(), 1);
        assert_eq!(o.sample_scores.len(), 3);
        assert_eq!(o.feature_scores.len(), 3);
    }

    #[test]
    fn sample_scores_recover_chi_square_distance() {
        // Scaling-1 euclidean distance between sample scores equals the
        // chi-square distance between rows of the original table.
        let t = FeatureTable::parse(ll_table().as_bytes(), '\t').unwrap();
        let o = Ordination::compute(&t, None);
        let r = t.n_samples();
        let a = o.eigvals.len();
        let chi = chi_square_rows(&t);
        for i in 0..r {
            for j in 0..r {
                let mut s = 0.0;
                for k in 0..a {
                    let d = o.sample_scores[i * a + k] - o.sample_scores[j * a + k];
                    s += d * d;
                }
                assert!(
                    (s.sqrt() - chi[i * r + j]).abs() < 1e-6,
                    "chi[{i}][{j}] {} vs {}",
                    s.sqrt(),
                    chi[i * r + j]
                );
            }
        }
    }

    fn chi_square_rows(t: &FeatureTable) -> Vec<f64> {
        let r = t.n_samples();
        let c = t.n_features();
        let grand: f64 = t.data.iter().sum();
        let mut col = vec![0.0_f64; c];
        for i in 0..r {
            for (j, cj) in col.iter_mut().enumerate() {
                *cj += t.data[i * c + j] / grand;
            }
        }
        let row: Vec<f64> = (0..r)
            .map(|i| (0..c).map(|j| t.data[i * c + j] / grand).sum())
            .collect();
        let mut d = vec![0.0_f64; r * r];
        for i in 0..r {
            for j in 0..r {
                let mut s = 0.0;
                for (k, &ck) in col.iter().enumerate() {
                    let pi = t.data[i * c + k] / grand / row[i];
                    let pj = t.data[j * c + k] / grand / row[j];
                    let diff = pi - pj;
                    s += diff * diff / ck;
                }
                d[i * r + j] = s.sqrt();
            }
        }
        d
    }
}
