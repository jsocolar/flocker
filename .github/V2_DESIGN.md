# flocker 2.0 design

This document records the design decisions for the version 2.0 development
cycle. It is intended to be sufficient context for a contributor beginning
work without access to the conversation in which the decisions were made.

The immediate release preceding this work is flocker 1.1, whose principal
change is the general `twolevel_single` occupancy family and the refactoring of
the augmented model as a special case of that family. Version 2.0 is the place
to make a deliberate compatibility break in the packed-data and
post-processing contracts before adding more ecological model classes.

## Priorities and scope

The first objective is a stable, dimension-independent packed-data schema and
a simpler internal architecture. New model classes should not be added on top
of the old dimension-dependent schema and then immediately migrated.

After the 2.0 foundation is in place, the highest-priority new model class is
N-mixture. Other plausible additions include two-level multiseason occupancy,
waiting-time likelihoods, false-positive occupancy models, and eventually CJS
and other wildlife models. These need not all enter the first 2.x release.

The design should make subsequent model families additive. In particular,
adding a family should not require another packed-data compatibility break or
new meanings for existing auxiliary Stan variables.

Threading is not part of the 2.0 deliverable. Known-state information is also
deferred. Both should remain possible without another schema redesign.

## Vocabulary

Use the following ecological levels consistently:

- **event**: one observation occasion or visit;
- **unit**: one closure unit, such as a site-season;
- **series**: an ordered collection of units in a dynamic model;
- **level-two group**: a collection of units, or potentially complete series,
  sharing a structural-presence state.

Posterior draws are not an ecological level. They always occupy the final
dimension of returned arrays.

## Public fitting functions

Keep `flock()` as the fitting function for occupancy models it naturally
supports. Give other ecological model classes clear literal names, beginning
with `fit_nmixture()` and potentially including `fit_cjs()`.

Do not invent parallel acronym puns for every model class. The connection
between names such as `flam()` and N-mixture models would impose unnecessary
cognitive work on users.

Some occupancy subclasses may eventually receive their own fitting functions,
for example `flock_fp()` for false-positive models or `flock_wt()` for
waiting-time models. This is especially appropriate when model classes are not
expected to cross and when separating them keeps most arguments to `flock()`
relevant to most calls.

Do not require users to construct a bespoke model-specification object.
Infer structural properties such as one versus two levels and single-season
versus multiseason from formatted data wherever possible. Retain explicit
arguments only for genuine choices that data cannot determine, such as the
form of dynamic occupancy or the initialization assumption.

Remove redundant fitting flags when the formatted-data type is authoritative.
For example, an `augmented` flag should not duplicate information already
encoded by the data object. Continue to reject incompatible formula/data
combinations early.

## Formula interface

Retain explicit `f_*` arguments rather than formula lists or public formula
specification objects. Keep all formula arguments together in each fitting
function signature and expose only formulas relevant to that function's model
classes. In `flock()`, move `f_meta` alongside the other formulas.

A combined `brmsformula` may remain available where advanced brms features
genuinely require it. Continue to reject `mvbrmsformula`: flocker has no
defined semantics for a genuinely multivariate response model.

The existing occupancy component vocabulary should remain stable where the
components apply:

- `det`: event-level detection;
- `occ`: unit-level occupancy;
- `colo`: unit transition colonization;
- `ex`: unit transition extinction;
- `auto`: unit transition autologistic effect;
- `Omega`: level-two group occupancy.

N-mixture should receive similarly clear abundance-specific terms rather than
forcing abundance parameters into occupancy terminology.

## Public data formatters

Keep `make_flocker_data()` as the occupancy formatter. Add separate ecological
formatters such as `make_nmixture_data()` and potentially `make_cjs_data()`.
Keep analogous shared arguments parallel across these functions and add only
genuinely model-specific arguments.

The initial N-mixture formatter is expected to be a thin wrapper around shared
packing internals. Its most obvious distinction from occupancy formatting is
that `obs` accepts nonnegative integer counts greater than one.

Make specialized occupancy formatters internal implementation helpers. This
includes the current static, dynamic, augmented, and two-level single-season
formatters. Their historical public role was largely to give context to errors
that bubbled out of early formatter implementations. Once internal, remove the
`make_flocker_data()` messages saying that all errors must be interpreted in
the context of a named specialized formatter.

### Data role

Replace `newdata_checks` with:

```r
data_role = c("model_data", "new_data")
```

Store the selected role in metadata. Downstream functions must not infer it.

`model_data` means the object obeys fitting-compatible rules. It can still be
used for fitted values, state inference, posterior prediction, likelihoods,
and other post-processing. The name describes the formatting contract, not an
exclusive use.

`new_data` may retain prediction-only positions, including trailing
multiseason units with zero visits and visits whose outcomes are unknown.
Fitting functions must reject `new_data`; post-processing functions accept
both roles.

Likelihood calculations use only observed outcomes. Prediction-only units and
events are neutral to `log_lik_flocker()`.

## Formatted-object and metadata contract

Standardize formatted objects as containers with exactly two conceptual
parts:

```r
list(
  data = <data frame used by brms and/or Stan>,
  metadata = <information needed to interpret and reconstruct packed rows>
)
```

`data` contains only covariates and auxiliary columns consumed by brms or the
Stan likelihood. Do not add unused columns merely to make post-processing
convenient. In particular, never mutate or augment `brmsfit$data` after fitting.

Metadata should contain, where applicable:

- the ecological data/model type;
- `data_role`;
- the creating flocker package version;
- a packing-schema version independent of the package version;
- original dimensions and dimension names;
- packed-to-original and original-to-packed permutations;
- mappings among events, units, series, groups, and sites;
- group, species, series, unit, and event names;
- covariate names classified by ecological level;
- the set of conceptual positions retained for fitting or prediction.

Attach the same metadata object unchanged to a fitted model. Use
`get_flocker_metadata()` as the sole internal access path for metadata from
either formatted data or fitted models. Keep that helper internal unless a
specific user-facing need emerges.

Accept the breaking change from top-level fields such as
`flocker_data$n_rep` to `flocker_data$metadata$n_rep`.

## Packed-data schema

The old schema creates a variable number of columns such as
`ff_rep_index1`, ..., `ff_rep_indexM` and, for dynamic models,
`ff_unit_index1`, ..., `ff_unit_indexT`. Version 2 should replace these with a
fixed-width schema whose number of columns does not depend on the maximum
number of visits or seasons.

Design packing around the outer level at which the likelihood factorizes, not
around minimizing the number of auxiliary columns in isolation:

- ordinary single-season occupancy: unit;
- ordinary multiseason occupancy: series;
- two-level occupancy: level-two group;
- N-mixture: unit or series, according to model type.

Pack each outer likelihood block contiguously. Represent nested children with
segment counts. Add precomputed starts only where they materially simplify or
speed random access. Avoid full child-index vectors unless a model genuinely
requires noncontiguous membership.

Use domain-specific names for counts and offsets rather than abstract
parent/child terminology. The exact columns remain a prototype decision; the
principles above are settled.

Explicitly represent zero-child segments. The prototype must settle how
zero-visit units are represented while retaining unit predictors and avoiding
arbitrary event-level values that would alter brms centering, spline bases,
random-effect levels, or other preprocessing. Synthetic rows are not inert to
brms merely because the custom likelihood ignores them.

For a possible future two-level multiseason model, allow a level-two group to
contain multiple complete series. Require each series to belong wholly to one
level-two group. Do not initially allow groups to cut across series: although
an ID mapping could represent that structure, the joint latent likelihood
generally would not factorize cleanly by either series or group.

Keep outer blocks clean enough that they could later be used as `reduce_sum`
slices, but do not implement or prototype threading during the 2.0 rewrite.

## Compatibility boundary

Version 2 establishes a hard boundary. Do not retain the old packing or
post-processing implementations and do not provide an automatic migration
path.

Store both the creating package version, for provenance, and a distinct
packing-schema version, for compatibility. Compatibility checks use the schema
version directly rather than parsing package versions.

Support exactly the current schema. Reject formatted objects with missing or
incompatible schema tags before fitting. Reject fitted models with incompatible
metadata in flocker-specific post-processing. Error messages should explain
that pre-2.0 objects require a compatible 1.x installation or refitting.
Native brms operations that do not depend on flocker metadata may continue to
work naturally.

There is intentionally no mitigation path that keeps old indexing code alive.
Removing that code is a primary benefit of the major release.

## Stable family definitions

Define one internal descriptor for every retained or newly added model family.
This is an implementation registry, not an object users construct. Each entry
should identify:

- required distributional parameters and links;
- required and optional formulas;
- supported formatted-data schemas;
- auxiliary Stan inputs and their meanings;
- outer likelihood factorization;
- pointwise log-likelihood unit;
- post-processing implementation;
- latent hierarchy;
- zero-visit and prediction-only semantics.

Custom Stan likelihood functions should be stable for a family and independent
of observed maximum visits, maximum seasons, or other dataset dimensions.
brms will still generate formula-specific model code; the objective is to stop
generating structurally different custom likelihoods and auxiliary schemas for
different data dimensions.

Preserve stable ecological ordering within each packed block and stable
meanings for all auxiliary inputs. Do not use Stan generated quantities for
latent-state reconstruction by default. Shared R post-processing can produce
states and coherent conditional samples without increasing saved fit size.
Add generated quantities only for a concrete capability that R
post-processing cannot supply reliably.

## Output contract

Return each fitted component only at its intrinsic ecological level:

- `det` by event;
- `occ`, `colo`, `ex`, and `auto` by unit or applicable transition;
- `Omega` once per level-two group.

Remove the current `unit_level` argument. Do not replace it with one hierarchy
argument that expands or aggregates components living at different levels.
An explicit expansion helper can be considered later if a real workflow needs
ancestor-level values repeated over descendants.

Keep `predict_flocker()` event-shaped. Keep `get_Z()` as an array for one-level
models and as a list with `unit` and `level2` elements for two-level models.
The latter guarantees that jointly sampled unit and group states are coherent.

Place posterior draws in the final dimension and never implicitly drop any
ecological or draw dimension, including for singleton dimensions and
one-draw requests. Summaries replace the draw dimension with a consistently
named summary dimension.

Preserve conceptual input dimensions and user-supplied dimension names.
Provide stable names for package-created entities such as augmented
pseudospecies. Values, names, and metadata mappings must remain synchronized
through all packing permutations.

### Unobserved multiseason units

Unit-level process output uses the conceptual unit grid of the input rather
than silently compressing to observed units.

For `model_data`, units that have no visits but are followed by later observed
units are internal parts of the latent process and receive values. Trailing
zero-visit units should not be retained in fitting rows merely to produce
output: their covariates could affect brms preprocessing despite contributing
no likelihood. Preserve their conceptual positions in metadata and return
`NA` there.

For `new_data`, retain trailing zero-visit units so process quantities can be
forecast. History-conditioned forecasts propagate from preceding observed
history. Hypothetical observations require a supplied event template and the
necessary event covariates; a unit with no supplied visits has no event-level
prediction.

Under explicit initialization, `occ` applies at the first unit and transition
components begin afterward. Under equilibrium initialization, first-unit
`colo` and `ex`, or `colo` and `auto`, are meaningful: they determine the
equilibrium initial state under the interpretation that first-unit covariates
were stable before the observed series began.

## Prediction availability

Relevant post-processing functions should use:

```r
output_scope = c("covariate_complete", "observed")
```

The default is `"covariate_complete"`.

`"observed"` forces output to `NA` wherever the corresponding observation is
`NA`. `"covariate_complete"` returns predictions wherever all quantities
required by the fitted model are available. An `NA` outcome provides no
conditioning information but does not itself suppress prediction.

Observation prediction requires all model-required unit and event covariates.
Latent-state prediction requires all model-required unit covariates.

For a history-conditioned dynamic prediction, covariates must be available
along the latent path leading to the target. Ordinarily that means all
preceding units in the series. An observed detection anchors the state as
present, after which only the path from the last such detection is needed.
Initialization rules determine the required path when no detection provides an
anchor.

Do not implement arbitrary bespoke masks in flocker. Users can construct and
apply any other desired mask to returned values themselves.

## Known states: deferred but anticipated

Do not add a public known-presence or known-absence API during the initial 2.0
packing rewrite.

The schema must nevertheless give every latent-state entity a stable canonical
index and retain relevant zero-visit units so future constraints align directly
with packed entities. Leave a clear extension point for tri-state constraints:
unknown, known present, and known absent. These are state constraints, not
ordinary covariates.

During schema prototyping, decide whether unknown-valued constraint vectors
should be reserved in the initial Stan/data contract or whether they can be
added later without compromising stable family programs.

When implemented, known states should be coherent across applicable one-level,
two-level, and multiseason models. Detections contradicting known absence are
errors. Documentation must discuss the data-generating or ascertainment
assumptions behind external presence information. Dynamic models can enforce a
known state by excluding the contradictory state from marginalization in both
the forward likelihood and forward-backward state reconstruction.

## Threading: explicitly deferred

Do not expose or implement threading in version 2.0.

Useful flocker threading would slice complete ecological likelihood blocks.
Simply asking brms to compile or run its ordinary threaded path is not
obviously correct because brms also slices linear-predictor/design-matrix
work. A possible future implementation could use the cmdstanr backend, enable
threading at C++ compilation (for example through `make/local`), and call a
threaded custom likelihood while leaving brms's own slicing disabled. A cleaner
supported brms mechanism, potentially an upstream option for custom families,
would be preferable.

Preserving contiguous outer blocks makes later threading an additive
optimization rather than a public API or schema break. If implemented, serial
and threaded likelihoods must have dedicated numerical-equivalence tests.

## Required verification

Every retained or added family should receive tests appropriate to its risk:

- exact likelihood comparisons against an independent R calculation;
- event/unit/series/group packing and reconstruction round trips;
- maximally ragged data, including singleton children and empty allowable
  segments;
- preservation of names through nontrivial packing permutations;
- singleton ecological dimensions;
- fits containing one draw and post-processing requests selecting one draw;
- extreme-probability numerical-stability tests where marginal mixtures or
  forward algorithms are involved;
- targeted simulation-based calibration for major new or refactored
  likelihoods;
- explicit rejection of incompatible schema versions and `new_data` supplied
  to fitting functions.

Where practical, retain characterization tests showing that the 2.0 serial
likelihood agrees with 1.1 for equivalent supported inputs even though the
serialized object format is intentionally incompatible.

Tests that require fitted Stan models should use the local fixture cache during
development and continue to exercise fresh fixture construction in the full CI
or release-validation path.

## Recommended implementation sequence

1. Freeze representative 1.1 behavior with characterization tests, especially
   packing round trips, exact likelihoods, output names and dimensions, and
   dynamic missing-unit behavior.
2. Prototype the fixed-width schema separately for ordinary single-season,
   multiseason, and two-level data. Resolve zero-visit and ghost-row behavior
   before broad code changes.
3. Introduce schema-version checks and the standardized `$data`/`$metadata`
   object contract, then make `get_flocker_metadata()` the only internal
   metadata access path.
4. Move shared packing into internal helpers and make specialized occupancy
   formatters internal. Keep `make_flocker_data()` as the stable public
   occupancy boundary.
5. Convert custom family definitions to the internal registry and stable,
   dimension-independent Stan likelihood functions.
6. Rework post-processing around intrinsic-level outputs, final draw
   dimensions, `data_role`, and `output_scope`. Remove old compatibility and
   replicated-output branches rather than carrying them forward.
7. Run the full verification matrix and deliberately test rejection of 1.x
   formatted objects and fits.
8. Once the 2.0 foundation is stable, introduce `make_nmixture_data()` and
   `fit_nmixture()` using the same schema, metadata, naming, and testing
   contracts.

Do not treat this sequence as a requirement to implement all steps in one pull
request. The schema prototypes and invariants should be reviewed before the
large mechanical migration begins.

## Additional technical debt and future work

The following work is useful but is not a prerequisite for settling the 2.0
schema unless it directly intersects changed code:

- Audit log-likelihood calculations across every family and work on the log
  scale wherever appropriate.
- Audit all package documentation for accuracy and consistency after the 2.0
  public APIs settle.
- Improve numerical stability in dynamic forward-backward post-processing and
  add a dedicated stress test with intentionally extreme predictors.
- Extend the shared two-level post-processing preparation used by `get_Z()` and
  `log_lik_flocker()` to one-level and multiseason models. During that work,
  inspect the current autologistic `get_Z()` path where
  `history_condition = FALSE`: historical code can pass `lps1$linpred_det`
  even though `lps1` is not defined on that branch.
- Add known-state constraints in a later coherent PR.
- Explore two-level multiseason, waiting-time, false-positive, and CJS models
  only after their public boundaries and likelihood factorizations are clear.

## Non-goals

- Backward-compatible processing of pre-2.0 flocker data or fits.
- A universal fitting function covering occupancy, abundance, capture-recapture,
  and every future wildlife model.
- Public model-specification or formula-list objects.
- Arbitrary user-supplied prediction masks.
- Threaded likelihoods in the initial 2.0 implementation.
- A public known-state API in the initial packing rewrite.
