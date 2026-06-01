use std::hint::black_box;
use std::path::PathBuf;
use std::process::Command;

use criterion::{Criterion, criterion_group, criterion_main};

fn bench_ca(c: &mut Criterion) {
    let bin = env!("CARGO_BIN_EXE_rsomics-correspondence-analysis");
    let table = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/golden/otu_mid.tsv");
    c.bench_function("rsomics-correspondence-analysis otu_mid", |b| {
        b.iter(|| {
            let out = Command::new(black_box(bin))
                .arg(&table)
                .args(["-t", "1"])
                .output()
                .unwrap();
            assert!(out.status.success());
        });
    });
}

criterion_group!(benches, bench_ca);
criterion_main!(benches);
