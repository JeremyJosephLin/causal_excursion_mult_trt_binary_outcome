# Consistency and asymptotic normality, Supplementary Section C

Three tables compare EMEE-catA with working-independence GEE (GEE-ind) and exchangeable
GEE (GEE-exch) at n = 20, 30, 40, 50, 100, under a nuisance model that is misspecified for
all three:

- Scenario 1, moderated CEE, randomization free of the state;
- Scenario 2, marginal CEE, randomization depending on the state;
- Scenario 1, marginal CEE, randomization free of the state.

## What to run

```sh
Rscript run_consistency_sim.R [nsim] [ncores]   # defaults: 1000, cores - 2
Rscript summarize_consistency.R [nsim]          # defaults: 1000
```

`run_consistency_sim.R` writes one `.RDS` per setting into `results/` (gitignored);
`summarize_consistency.R` reads them and prints the metrics table for each setting
followed by the LaTeX body of the corresponding table in the supplement. Both take the
replicate count as their first argument, so `Rscript run_consistency_sim.R 20 4` is a quick
check that everything is wired up.

`dgm_consistency.R` holds the generative model, which is the one displayed in the
supplement, together with `true_cee()`, the true causal excursion effects it implies.
Scenario 1 randomizes at 1/3 per level and the estimator is fit at the default
`p~ = 1/3`, so `J_t = 1`; Scenario 2 keeps `p~ = 1/3` while the randomization depends on
the state, which is the `J_t != 1` case that setting exists to illustrate. Neither should
be changed to the Drink Less convention of setting `p~` to the randomization probability:
these two settings are what makes the comparison in the supplement informative.

## The two subfolders

`sim_consistency_moderated/` and
`sim_consistency_marginal_model_with control_non_const_pt/` are the earlier per-scenario
pipelines, written to run on a cluster. They are kept for the record and are **not** what
produced the published tables: they use two different generative models, neither of which
matches the generative model printed in the supplement, and the moderated pipeline fits
its second GEE comparator with working independence rather than an exchangeable working
correlation. `run_consistency_sim.R` replaces both.
