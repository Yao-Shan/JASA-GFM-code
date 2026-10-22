# Seeds and repeated networks

The completed with-intercept studies used design seed 238687 and intercept seed
817285. For each fixed n and Lambda scenario, these define one latent design.

- `design_seed` generates points in a K-dimensional unit ball; they are centered
  and whitened to obtain U with column means zero and U' U = n I.
- `alpha_seed` generates raw iid Uniform(-2,0) values a, followed by alpha = a / ||a||2.
  The normalized entries are dependent and are no longer iid Uniform(-2,0).
  Only truth is normalized; the fitted alpha vector has no unit-norm constraint.
- `replication_seed_base` supplies each network seed: replicate r uses base+r,
  for r=1,...,1000. Base 2000 therefore means seeds 2001-3000, not 1000 repeated
  uses of seed 2000. The base controls network observations conditional on truth.

A fixed n, design seed, intercept seed and Lambda scenario hold truth constant
across replications. Changing the network seed base probes variability under another
batch of random network observations. Shared seeds across models or Lambda settings
are common-random-number designs and do not establish independence.
Changing n changes the generated/whitened U and the intercept normalization denominator.
Identical seed labels across n do not imply the same network or independent repeats.

The five report groups are:

| Group | n | Replication seed base | Actual network seeds |
|---|---:|---:|---|
| G1 | 1000 | 2000 | 2001-3000 |
| G2 | 1000 | 12000 | 12001-13000 |
| G3 | 1000 | 32000 | 32001-33000 |
| G4 | 2000 | 32000 | 32001-33000 |
| G5 | 2000 | 2000 | 2001-3000 |

Each group has four models and two Lambda scenarios, diag(0.5,0.5) and diag(0.7,-0.5),
with 1000 attempts per scenario/model. T2 reuses the corresponding T1 fit.
The groups are presented separately, including results that compare less favorably.
These computations are simulation diagnostics, not formal independent validation
of the proofs or a license to discard or replace unfavorable replications.
