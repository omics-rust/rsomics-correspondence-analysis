use std::collections::HashMap;
use std::path::PathBuf;
use std::process::Command;

const EPS_EIG: f64 = 1e-9;
const EPS_SCORE: f64 = 1e-6;

fn ours_bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-correspondence-analysis"))
}

fn fixture(name: &str) -> String {
    format!("{}/tests/golden/{}", env!("CARGO_MANIFEST_DIR"), name)
}

fn oracle_script() -> String {
    format!("{}/tests/oracle_skbio.py", env!("CARGO_MANIFEST_DIR"))
}

/// scikit-bio is the named oracle; skip loudly if it (or python) is unavailable.
/// `RSOMICS_SKBIO_PYTHON` overrides the interpreter (e.g. an isolated venv).
fn skbio_python() -> Option<String> {
    let mut candidates = Vec::new();
    if let Ok(p) = std::env::var("RSOMICS_SKBIO_PYTHON") {
        candidates.push(p);
    }
    candidates.push("python3".into());
    candidates.push("python".into());
    for py in candidates {
        let probe = Command::new(&py)
            .args(["-c", "import skbio.stats.ordination"])
            .output();
        if let Ok(out) = probe
            && out.status.success()
        {
            return Some(py);
        }
    }
    eprintln!(
        "SKIP: scikit-bio not importable — install `scikit-bio` to run the live differential"
    );
    None
}

struct Ca {
    eigvals: Vec<f64>,
    proportion: Vec<f64>,
    samples: HashMap<String, Vec<f64>>,
    features: HashMap<String, Vec<f64>>,
}

fn parse_ours(text: &str) -> Ca {
    let mut eigvals = Vec::new();
    let mut proportion = Vec::new();
    let mut samples = HashMap::new();
    let mut features = HashMap::new();
    let mut section = 0u8;
    for line in text.lines() {
        match line {
            "# eigenvalues" => {
                section = 1;
                continue;
            }
            "# samples" => {
                section = 2;
                continue;
            }
            "# features" => {
                section = 3;
                continue;
            }
            _ => {}
        }
        if line.starts_with('\t') {
            continue;
        }
        let mut it = line.split('\t');
        let label = it.next().unwrap();
        let vals: Vec<f64> = it.map(|s| s.parse().unwrap()).collect();
        match (section, label) {
            (1, "eigval") => eigvals = vals,
            (1, "proportion_explained") => proportion = vals,
            (2, _) => {
                samples.insert(label.to_string(), vals);
            }
            (3, _) => {
                features.insert(label.to_string(), vals);
            }
            _ => {}
        }
    }
    Ca {
        eigvals,
        proportion,
        samples,
        features,
    }
}

fn parse_oracle(text: &str) -> Ca {
    let mut eigvals = Vec::new();
    let mut proportion = Vec::new();
    let mut samples = HashMap::new();
    let mut features = HashMap::new();
    for line in text.lines() {
        let mut it = line.split('\t');
        let label = it.next().unwrap();
        let vals: Vec<f64> = it.map(|s| s.parse().unwrap()).collect();
        if let Some(id) = label.strip_prefix("S:") {
            samples.insert(id.to_string(), vals);
        } else if let Some(id) = label.strip_prefix("F:") {
            features.insert(id.to_string(), vals);
        } else if label == "eigvals" {
            eigvals = vals;
        } else if label == "proportion_explained" {
            proportion = vals;
        }
    }
    Ca {
        eigvals,
        proportion,
        samples,
        features,
    }
}

fn ours_output(table: &str) -> String {
    let out = Command::new(ours_bin())
        .arg(fixture(table))
        .output()
        .expect("run ours");
    assert!(
        out.status.success(),
        "ours failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).unwrap()
}

fn approx(a: f64, b: f64, eps: f64) -> bool {
    (a - b).abs() <= eps + eps * b.abs()
}

/// Per-axis sign from the largest-magnitude score, giving a stable orientation
/// independent of the eigensolver's arbitrary sign.
fn axis_sign(scores: &HashMap<String, Vec<f64>>, ids: &[&String], a: usize) -> f64 {
    let mut best = 0.0_f64;
    let mut sign = 1.0_f64;
    for id in ids {
        let v = scores[*id][a];
        if v.abs() > best {
            best = v.abs();
            sign = if v < 0.0 { -1.0 } else { 1.0 };
        }
    }
    sign
}

fn compare(ours: &Ca, theirs: &Ca, tag: &str) {
    assert_eq!(ours.eigvals.len(), theirs.eigvals.len(), "{tag} axis count");
    for (a, &o) in ours.eigvals.iter().enumerate() {
        assert!(
            approx(o, theirs.eigvals[a], EPS_EIG),
            "{tag} eigval CA{} {o} vs {}",
            a + 1,
            theirs.eigvals[a]
        );
    }
    for (a, &o) in ours.proportion.iter().enumerate() {
        assert!(
            approx(o, theirs.proportion[a], EPS_SCORE),
            "{tag} proportion CA{} {o} vs {}",
            a + 1,
            theirs.proportion[a]
        );
    }
    let n = ours.eigvals.len();
    for (label, ours_map, theirs_map) in [
        ("sample", &ours.samples, &theirs.samples),
        ("feature", &ours.features, &theirs.features),
    ] {
        let ids: Vec<&String> = ours_map.keys().collect();
        for a in 0..n {
            let so = axis_sign(ours_map, &ids, a);
            let st = axis_sign(theirs_map, &ids, a);
            for id in &ids {
                let o = ours_map[*id][a] * so;
                let t = theirs_map[*id][a] * st;
                assert!(
                    approx(o, t, EPS_SCORE),
                    "{tag} {label} {id} CA{} {o} vs {t} (sign-aligned)",
                    a + 1
                );
            }
        }
    }
}

/// Always-on gate: ours vs a committed skbio-captured golden.
#[test]
fn matches_committed_golden() {
    let ours = parse_ours(&ours_output("otu_small.tsv"));
    let golden = std::fs::read_to_string(fixture("otu_small.ca.golden")).unwrap();
    let theirs = parse_oracle(&golden);
    compare(&ours, &theirs, "golden");
}

/// Live differential against an installed scikit-bio, loud-skip if absent.
#[test]
fn matches_live_skbio() {
    let Some(py) = skbio_python() else { return };
    for table in ["otu_small.tsv", "otu_mid.tsv"] {
        let out = Command::new(&py)
            .arg(oracle_script())
            .arg(fixture(table))
            .output()
            .expect("run oracle");
        assert!(
            out.status.success(),
            "oracle failed: {}",
            String::from_utf8_lossy(&out.stderr)
        );
        let ours = parse_ours(&ours_output(table));
        let theirs = parse_oracle(&String::from_utf8(out.stdout).unwrap());
        compare(&ours, &theirs, table);
    }
}
