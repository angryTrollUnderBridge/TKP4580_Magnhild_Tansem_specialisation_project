
# NTNU TKP4580 Code for parameter estimation, cross validation, practical identifiability analysis and plotting functions
Author: Magnhild Tansem, magnhirt@stud.ntnu.no
Code assistant: Claude code

MATLAB code. Six kinetic growth models are fitted to bioreactor data
(biomass, flow-cytometry strain fractions, HPLC sugars, off-gas CO2) and
compared by how well they predict an experiment they were not trained on.

**The experimental data is not included in this repository**. Placeholders mark where the data goes; see
[data/README.md](data/README.md) for the expected file names and formats. If you want to ask for access to the data, contact the author :)

The code repository was cleaned up by claude code for readability. It
contains exactly what produced the final model-comparison sweep of
**7 November 2026**. The saved results (fitted parameters, logs and figures)
are included and can be plotted as needed.

## Reproducing the results

Requires MATLAB (tested with R2025b) with the Optimization, Global
Optimization, Statistics and Machine Learning, and Parallel Computing
toolboxes.

```matlab
setup_paths          % once per session
runAllFolds          % 18 fits (6 models x 3 folds), ~40 min on a local parallel pool
makeFigures          % held-out scores and every figure, from the saved fits
```

Refitting requires the data files described in `data/README.md`.
`runAllFolds` writes new, time-stamped runs next to the stored ones; use
`runAllFolds('OutDir', 'rerun')` to keep them apart, and then
`makeFigures('ResultsDir', 'rerun/results', 'OutDir', 'rerun/figures')`.
A single fit: `fitModel('D', [5 6], 4, 20260903)` (model, training
experiments, held-out experiment, seed).

`makeFigures` needs no refitting. The sensitivity heatmaps and correlation
plots are drawn from the stored `results/sensitivity/` files; deleting those
(or `makeFigures('RecomputeSensitivity', true)`) recomputes them, which takes
a few minutes per experiment.

## Repository layout

```
runAllFolds.m        the sweep: every model on every leave-one-experiment-out fold
fitModel.m           one fit (data preparation, stage 1, stage 2, saving); all fit settings are here
makeFigures.m        held-out scoring and all figures
setup_paths.m

src/data/            loadExperimentData   raw files -> one struct per experiment
                     prepareFitData       initial condition and scored observations
                     readFcsExpressReport, groupAndAverageTrips
src/models/          model0, modelA ... modelD   ODE right-hand sides (see below)
                     getParamConfig       parameter names, bounds and units of every model
                     getModelSpec         model name -> model functions
src/fitting/         geneticSearch, paramLeastSquares   stage 1 (GA, then lsqnonlin)
                     estimateCO2Yields                  stage 2 (lsqlin)
                     residuals, residualsSingle, scoredObservations, noiseModel
                     logRunResult
src/analysis/        evalHeldOut, computeHeldOutMetrics  held-out prediction scores
                     computeSensitivities                local sensitivities + Fisher information
                     reduceFIM, loadRun
src/figs/            one function per figure type

data/                placeholder: where the (not included) data files go, see data/README.md
results/             the 18 fitted runs of the 2026-11-07 sweep (+ logs/, figures/)
```

## Data

Not included in this repository. The code expects three triplicate batch
cultivations (experiments 4, 5 and 6) on a glucose/mannose/xylose mixture:

- **Dense logger data**: biomass probe (g/L), off-gas CO2 (vol%), air flow.
- **Flow cytometry (FCM)**: fractions of the four strains (r, g, rg, e), in
  replicate bursts of 2-3 samples. Rows flagged in the export's
  `exclude_flag` column are dropped.
- **HPLC sugars**: glucose, mannose and xylose, typed into the placeholders in
  `loadExperimentData.m` (mM, converted to g/L).

`loadExperimentData` puts everything on one time axis (hours since
inoculation) and applies the biomass corrections specified by the
experimenters (constants at the top of `loadExperimentData.m`). The first
probe readings of experiments 5 and 6 (probe stabilisation) are removed, as
are samples drawn before inoculation. Experiment 4's biomass probe stops
early; after that its manual dry-weight samples are used.

## Models

All models share the state vector `[Xtot Xr Xg Xgr Sglc Sman Sxyl CO2]`
(total biomass, three of the four strains, three sugars, off-gas CO2; the
wild type is `Xe = Xtot - Xr - Xg - Xgr`). Strains: r grows on glucose and
mannose, g on glucose and xylose, rg on all three, e on glucose only.
Growth is Monod-type per (strain, sugar) pair.

| model | growth rate | glucose preference | parameters (kinetic + CO2 yields) |
|---|---|---|---|
| 0  | one lumped population, per sugar | no | 10 + 3 |
| 0B | one lumped population, per sugar | yes (`KI`) | 11 + 3 |
| A  | per strain, independent of the sugar | no | 24 + 8 |
| B  | per strain, independent of the sugar | yes (`KI`) | 25 + 8 |
| C  | per (strain, sugar) pair | no | 28 + 8 |
| D  | per (strain, sugar) pair | yes (`KI`) | 29 + 8 |

The glucose preference multiplies mannose and xylose uptake by
`1/(1 + Sglc/KI)`. Models 0/0B do not resolve strains: their strain states
only carry the initial composition forward, and their FCM data is not fitted.
Each model is split into files `model<X>.m` (full, 8 states),
`model<X>_kinetic.m` (stage 1, no CO2), `model<X>_co2Basis.m` (stage 2) and
`model<X>_pairRates.m` / `model0_rates.m` (the growth rates, shared by the
other three).

Parameter bounds (`getParamConfig.m`) are the same for every model:

| class | bounds | unit |
|---|---|---|
| death rate `kd*` | 1e-5 – 0.01 | 1/h |
| max. growth rate `mumax*`, `mu_*` | 0.01 – 0.6 | 1/h |
| Monod constant `K*` | 1e-3 – 1 | g/L |
| biomass yield `Y*` | 0.1 – 1 | g/g |
| CO2 yield `Y*co2` | 1e-4 – 1 | – |
| glucose preference `KI` | 0.1 – 10 | g/L |

## Fitting

**Folds.** Leave-one-experiment-out over experiments 4, 5, 6: train on two,
hold out the third (`runAllFolds`).

**Initial condition and scored data** (`prepareFitData`). Each experiment is
simulated from an *anchor*: the first FCM burst that has a biomass reading
and an HPLC sample within 0.5 h. The initial biomass and strain split are the burst means, the
initial sugars are interpolated from the HPLC series. Data before the anchor
and the anchor burst itself are not scored. Scored are all HPLC samples, all
FCM samples (as strain biomass, fraction × total biomass), and the biomass at
the times of those samples (plus exp 4's manual samples). For exp 4, whose
probe stops early, the FCM fractions are converted with a smoothed biomass
curve.

**Objective.** Weighted residuals `(simulated - measured)/sigma` with the
assumed measurement errors of `noiseModel.m` (biomass 0.005 + 1 %, strain
fractions 0.017/0.045/0.046, sugars 0.3 g/L, CO2 0.02 + 2 %). These are
reasoned estimates, not calibrated against replicates.

**Stage 1: kinetic parameters**, in log10 space, both training experiments
together:
1. genetic algorithm (`ga`): population 400, at most 40 generations, seed
   20260903; ODE tolerances 1e-4 / 1e-6;
2. `lsqnonlin` from the GA result: central finite differences (step 1e-4),
   at most 200 iterations, function tolerance 1e-6; ODE tolerances
   1e-9 / 1e-11.

The model integrated in stage 1 has no CO2 state and no air-flow input; CO2
does not feed back into the other states, and leaving it out makes each
solve much cheaper.

**Stage 2: CO2 yields.** With the kinetics fixed the CO2 balance is linear in
`1/Yco2`, so the yields are one bounded linear least-squares problem
(`lsqlin`) on the full off-gas record.

Every run saves `results/<runTag>.mat` (fitted parameters, Jacobian,
confidence intervals, prepared data, settings), a console log in `logs/`, and
a row in `results/run_log.csv`; `runAllFolds` adds
`results/sweep_summary_<timestamp>.csv`.

## Evaluation and figures (`makeFigures`)

- **Held-out prediction** (`evalHeldOut`, `computeHeldOutMetrics` →
  `results/heldout_metrics.mat`, `figures/heldout/`): each fit simulates its
  held-out experiment from that experiment's own measured anchor and is scored
  by RMSE per channel, and by the RMSE of the predicted strain *fractions*
  against FCM (compared with the FCM replicate SD).
- **Fits** (`figures/fit_train/`, `figures/fit_held_out/`).
- **Local sensitivities** (`computeSensitivities`, `figures/sensitivity/`):
  forward sensitivity analysis of the full model at the fitted parameters,
  with respect to all parameters (log10) and the measured initial values;
  time-averaged and peak magnitude.
- **Identifiability** (`figures/corr/`): parameter correlations from the
  Fisher information of the kinetic parameters over the scored data.
- **Reduced parameterisations** (`reduceFIM`, `figures/corr_reduced/`): the
  same for models 0, 0B and D with Monod constants fixed, yields pooled per
  sugar and one death rate (plus, for D, fixed rg growth-rate ratio and death
  rate), and the resulting confidence intervals for model D
  (`modelD_forest_reduced`). Computed from the Fisher information, not refitted.

## Notes on the stored results

- Every GA stopped on its 40-generation cap, so stage 1 refines a truncated
  global search; all `lsqnonlin` runs ended on the function tolerance
  (exit flag 3). One fit exists per (model, fold).
- The GA's fitness evaluations run in parallel; results are reproducible with
  the same seed on the same setup, but bit-identical results across machines
  are not guaranteed.
- The `AIC`/`BIC` columns of the stored `results/run_log.csv` and the `AIC`,
  `BIC` fields in the stored `.mat` files were computed with a parameter count
  of 1 instead of the number of fitted kinetic parameters (a bug since fixed in
  `logRunResult.m`). Correct values: `AIC = n ln(SSE/n) + 2k`,
  `BIC = n ln(SSE/n) + k ln n`, with `n = size(sol.J,1)`, `SSE = sol.resnorm`,
  `k = height(P29)`.
- The stored `.mat` files also contain a `boundsSet` field (`'v2'`), the name
  the bounds above had at the time.
