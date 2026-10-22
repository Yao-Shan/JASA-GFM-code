# Checks

## Run from the repository root

```sh
python3 tools/verify_package.py
Rscript --vanilla without_intercept/tests/test_gfm.R
Rscript --vanilla without_intercept/tests/test_simulation_config.R
Rscript --vanilla without_intercept/tests/test_plotting.R
Rscript --vanilla with_intercept/tests/test_package.R
Rscript --vanilla with_intercept/tests/test_plots.R
```

The package verifier checks file hashes, portable paths, directory structure,
documentation links and deferred-artifact exclusions. The R checks cover estimation formulas,
configuration and plotting data handling. They do not run the paper-scale
Monte Carlo study.

For small fitting and sequential/parallel execution checks, also run:

```sh
Rscript --vanilla with_intercept/tests/test_small_run.R
Rscript --vanilla with_intercept/tests/test_runtime_free.R
```

`test_runtime_free.R` verifies ordinary R execution without a runtime folder.
These small-network checks are implementation checks, not evidence of
statistical calibration or a proof of the asymptotic results.

## Code-only release: 2026-10-22

Both models contain code, guides and tests. Saved `results/`, paper `reports/`
and local `runtime/` directories are excluded; all remain in the local research
workspace. The verifier expects five model directories, and ignore rules prevent
accidental inclusion of deferred artifacts.

The with-intercept plotting notebooks were consolidated into
`with_intercept/plots/plot_figures.Rmd`; its executable chunks and default
parameters are unchanged. Numerical source and simulation settings are preserved.
The plotting tests and the consolidated notebook's inert default render passed
before this packaging change. No paper-scale simulation was rerun.

Current file checks are recorded in `VALIDATION.json`. Code provenance remains
in each model's documentation. Historical preparation and report checks are
retained locally with the deferred materials.

For this code-only release, the following checks passed on 2026-10-22:

| Check | Result |
| --- | --- |
| No-intercept estimation | 132 checks passed |
| No-intercept simulation configuration | 120 checks passed |
| No-intercept plotting | Passed |
| Intercept package and settings | 202 checks passed |
| Intercept plotting | Passed |

Deferred report files were verified byte-for-byte against the preceding
package. Numerical modules and simulation settings were unchanged by packaging.
