"""Emit skbio.stats.ordination.ca() results for a feature-table TSV.

Usage: python oracle_skbio.py table.tsv  -> TSV on stdout:
  eigvals\t...
  proportion_explained\t...
  S:<sample_id>\t<CA1>\t<CA2>...      (sample scores)
  F:<feature_id>\t<CA1>\t<CA2>...     (feature scores)
"""
import sys
import numpy as np
from skbio.stats.ordination import ca

path = sys.argv[1]
ids = []
feat = None
rows = []
with open(path) as fh:
    for line in fh:
        line = line.rstrip("\n")
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.split("\t")
        if feat is None:
            feat = [p.strip() for p in parts[1:]]
            continue
        ids.append(parts[0].strip())
        rows.append([float(x) for x in parts[1:]])

X = np.array(rows, dtype=np.float64)
res = ca(X, scaling=1, sample_ids=ids, feature_ids=feat)

ev = np.asarray(res.eigvals.values, dtype=np.float64)
pe = np.asarray(res.proportion_explained.values, dtype=np.float64)
samp = np.asarray(res.samples.values, dtype=np.float64)
feats = np.asarray(res.features.values, dtype=np.float64)

print("eigvals\t" + "\t".join(repr(float(x)) for x in ev))
print("proportion_explained\t" + "\t".join(repr(float(x)) for x in pe))
for sid, row in zip(ids, samp):
    print("S:" + sid + "\t" + "\t".join(repr(float(x)) for x in row))
for fid, row in zip(feat, feats):
    print("F:" + fid + "\t" + "\t".join(repr(float(x)) for x in row))
