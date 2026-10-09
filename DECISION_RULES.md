# Pre-registered decision rules — GNoME PBE/r2SCAN audit
Written 2026-10-09, before any stratified analysis.

Phase 0 (done): 139,541 rows. Flip rate (decomp energy > 0):
MP hull 13.0% (18,198) | MP+OQMD 15.9% (22,222) | All 25.8% (36,021).
Reproduces published 86.8% quaternary survival.

PRIMARY RULE. Flip rate by n_elements x chemistry, strata n>=500,
bootstrap 95% CI, 1000 resamples.
- GO standalone: spread >= 5pp, non-overlapping CIs, at fixed n_elements.
- FOLD into Paper A: spread < 3pp or all CIs overlap.
- Grey 3-5pp: proceed only if AI-vs-known shows separate effect.

SECONDARY. Pre-registered claim: lanthanide flip rate <= global rate
at fixed n_elements. Report confirmed or refuted either way.

ROBUSTNESS (all required): threshold curves 0-50 meV/atom; both
reference hulls; recomputation coverage by family; Li/Na/Mg/Ge/Ga
potential confound as separate stratum; +U vs non-+U split;
stratify by PBE hull-distance bin.

KILL. Abandon if a chemistry-resolved GNoME PBE/r2SCAN flip paper
is found, or if spread vanishes AND AI-vs-known is null.
