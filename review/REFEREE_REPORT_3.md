# Referee report, round 3: hostile read of Draft 5

**Stance.** I read this as the reviewer who wants to reject it. I ran the checks the paper
should have run. Where the paper survives I say so; where it does not I say what changes.

**Verdict: major revision, but the direction is upward.** Three analyses the paper
needed were missing. Two came back strongly in the paper's favour and one overturns the
framing of a central section. The corrected paper is stronger than the submitted one,
but the submitted one contained a claim that was incomplete in a way a domain expert
would have caught.

---

## 1. A pre-registered robustness gate was never run

The pre-registration committed at `82d49db` lists six robustness gates. Gate 6 reads:
*"Near-hull concentration. Stratify by PBE decomposition energy bin and test whether
chemistry effects survive within bins."* Draft 5 cites Supplementary Note 2's observation
that disagreement concentrates near the hull, uses hull distance as a covariate in the
propensity model, and never performs the stratified test.

This is not a technicality. It is the single most obvious alternative explanation for the
headline result. If 3*d* oxides sit closer to the PBE hull than intermetallics, they will
flip more often for no reason connected to the functional. A paper that claims
pre-registration as a strength and silently skips a pre-registered gate invites a referee
to check the repository and find the gap.

**I ran it.** Quaternary materials, matched to PBE energies (117,018 after excluding 25
that were not PBE-stable to begin with, see §5):

| PBE hull distance (eV/atom) | oxide+3d | oxide | pnictide+3d | pnictide | interm+3d | interm |
|---|---|---|---|---|---|---|
| −0.1 to −0.05 | 11.1% | 6.1% | 14.6% | 1.3% | 3.2% | 0.5% |
| −0.05 to −0.02 | 23.1% | 15.5% | 19.8% | 8.2% | 7.9% | 1.2% |
| −0.02 to −0.01 | 29.4% | 24.3% | 28.3% | 20.8% | 16.0% | 4.3% |
| −0.01 to 0 | 40.1% | 27.8% | 37.7% | 27.0% | 27.2% | 15.9% |

The gradient survives inside every bin. Within the −0.05 to −0.02 band, 3*d* oxides flip
at 23.1% and plain intermetallics at 1.2%, a 21.8-point spread. And the cells are not
differently positioned: median PBE hull distance is −20.8 meV/atom for 3*d* oxides and
−22.4 for plain intermetallics. Hull proximity and chemistry act roughly independently.

**Outcome: the paper passes, decisively. But it must report this, and it must say the
gate was pre-registered and run late.**

## 2. The paper reports flip rates when the underlying quantity is a directional shift

A flip is a binary summary of a continuous change: the r²SCAN decomposition energy minus
the PBE one. The paper never examines that quantity. It should, because the sign is the
story.

Quaternary materials, Δ = E_dec(r²SCAN) − E_dec(PBE), in meV/atom:

| Cell | n | mean Δ | median Δ | s.d. | share shifted up |
|---|---|---|---|---|---|
| oxide + 3d | 3,178 | **+11.3** | +4.8 | 52.6 | 60% |
| pnictide + 3d | 9,095 | **+9.6** | +3.6 | 41.8 | 60% |
| oxide | 6,059 | +3.6 | +1.5 | 45.1 | 50% |
| pnictide | 10,380 | −1.8 | −0.8 | 39.0 | 50% |
| intermetallic + 3d | 16,308 | −1.8 | −1.0 | 33.1 | 50% |
| intermetallic | 20,866 | **−12.4** | −4.8 | 36.0 | **30%** |

The spread of Δ (standard deviation 33–53 meV/atom) is comparable across cells. What
differs is the *mean*. r²SCAN systematically destabilises 3*d* oxides and 3*d* pnictides
relative to PBE and systematically *stabilises* plain intermetallics. Only 30% of
intermetallics shift upward at all. The flip-rate gradient is not "some chemistries are
noisier"; it is "r²SCAN disagrees with PBE in opposite directions depending on chemistry,
and flips are where the upward shift crosses zero."

This is a better paper than the one submitted. It is also a direct connection to the
Kothakonda and Bartel results already cited: those papers found GGA and meta-GGA
"comparable" for intermetallics; this analysis shows they are comparable in *scatter*
while differing in *mean*, and that the mean difference points toward greater stability
under r²SCAN. Whether that is an improvement is an experimental question the paper
correctly declines to answer.

**Outcome: add this as a primary result. It belongs in the abstract.**

## 3. Section 4.3's conclusion is incomplete in a way that misleads

Draft 5 argues: the effect survives outside oxides and fluorides, where no *U* is applied,
"therefore the effect is not the correction scheme." True as far as it goes. But I asked
the sharper question: within each anion class, how much does the presence of Fe raise the
flip rate?

| Anion class | n (with Fe) | flip with Fe | flip without Fe | excess |
|---|---|---|---|---|
| boride/carbide/silicide | 1,268 | 54.0% | 13.0% | **+41.0** |
| pnictide | 1,264 | 50.3% | 15.0% | +35.3 |
| intermetallic | 1,904 | 42.2% | 8.0% | +34.2 |
| chalcogenide | 154 | 41.6% | 11.5% | +30.0 |
| halide | 476 | 24.4% | 10.9% | +13.5 |
| **oxide** | 610 | 24.1% | 19.1% | **+5.0** |

Cobalt shows the same pattern (+34 in borides, +22 in intermetallics, +4.6 in oxides).

The iron excess is *smallest* in oxides — the one class where PBE+*U* is applied — and
largest in borides, pnictides and intermetallics, where PBE runs without *U*. The
correction scheme does not *create* the disagreement. It *partially suppresses* it, in
the subset of chemistry where it is applied. The largest PBE/r²SCAN disagreements in the
whole dataset are open-shell 3*d* elements in non-oxide environments, which is exactly the
regime where PBE has no empirical correction at all.

This inverts the current framing. "Not the correction scheme" should become: *"The
correction scheme does not generate the effect; where applied, it partially mitigates it.
Disagreement is largest where an open-shell d element sits in chemistry that receives no
U."* That is both more accurate and more useful. It also connects to Warford et al.
differently than the paper currently does — selective *U* does not merely create two
incompatible surfaces, it leaves the uncorrected surface much further from r²SCAN.

Nickel, for the record, shows no excess in intermetallics (+2.0) and a small one in oxides
(+10.3). Ni is often low-spin or non-magnetic in metallic environments. This is consistent
with *magnetic* open-shell *d* being the operative variable, not list membership, and the
paper should say so rather than calling nickel "unexplained."

**Outcome: rewrite §4.3. The current version is not wrong but a specialist will see
it stopped one question short.**

## 4. The reference-set choice is under-justified and the paper buries the alternative

The paper takes the Materials-Project-referenced decomposition energy as primary because
it reproduces DeepMind's 86.8%. It mentions that the "All" reference gives 25.8% and
moves on. That number is double the headline, and the choice between them is not
innocent.

Anion-class flip rates under each:

| Anion | MP reference | "All" reference |
|---|---|---|
| oxide | 19.3% | 29.5% |
| pnictide | 16.7% | 29.6% |
| boride/carbide/silicide | 15.0% | 29.6% |
| chalcogenide | 12.1% | 21.6% |
| halide | 11.5% | 24.2% |
| intermetallic | 9.2% | 22.3% |

The ordering is roughly preserved, but the spread compresses from 2.1× to 1.4×, and
intermetallics gain the most (+13 points) from the richer reference set. The plausible
reading — and the paper should say it is a reading, since DeepMind's documentation does
not define "All" — is that "All" includes GNoME's own structures as competing phases,
making the hull more complete, and that intermetallics have the least complete MP hulls
because that is where GNoME found the most new material.

If so, this is a *partial measurement* of the hull-completeness bias the limitations
section calls untestable. It is testable, the paper has the data, and the answer is: the
gradient survives but shrinks. That should replace the current "we cannot quantify this."

**Outcome: promote the "All" comparison to a robustness result. Soften the chemistry
spread claim accordingly.**

## 5. Minor but must fix

- **25 recomputed materials were not PBE-stable** (PBE decomposition energy > 0). They
  cannot "flip." Exclude them and say so. 0.02%, immaterial to any number, but a referee
  who finds it will wonder what else was not checked.
- **Per-element rates are not disjoint.** A Fe-oxide counts toward both Fe and O. The
  paper should say this, and should report the within-anion excess (§3 above) rather than
  the raw per-element rate as the headline element-level result.
- **"Thirtyfold" per element** rests on actinium (1.4%, n = 2,181) at the bottom.
  Actinides carry the same frozen-*f* pseudopotential caveat as lanthanides. State that
  the bottom of the element ranking is actinide-dominated and that the "thirtyfold" figure
  should not be over-read.
- **The GNoME-vs-MP comparison has an asymmetric bias.** MP builds r²SCAN hulls only where
  every reference phase has an r²SCAN calculation; GNoME's hulls are not so constrained.
  GNoME flip rates are therefore biased *downward* relative to MP's across the board. Where
  GNoME looks worse (3*d* oxides, +17 points) the true gap is larger. Where GNoME looks
  better (intermetallics, −7.7) part of the advantage may be hull incompleteness. The
  paper states the comparison's limits but not this directional asymmetry.
- **Section 4.3 title** "The driver is not the correction scheme" should go.

## 6. What survived

The clustering correction, the propensity analysis, the repeated cluster-wise budget
evaluation, the multiple-comparison correction and the sandwich regression errors all
held up. I re-derived none of them from scratch but the reported numbers are internally
consistent with the raw tables in `results/`.

## 7. Recommendation

Revise as follows, in order:
1. Add the hull-distance-stratified table and state that it is pre-registered gate 6.
2. Add the Δ analysis as a primary result; move it into the abstract.
3. Rewrite §4.3 around the within-anion Fe/Co excess.
4. Promote the "All" reference comparison to robustness; revise the hull-completeness
   limitation to "partially measured" and give the compressed spread.
5. Exclude the 25 non-PBE-stable materials; add the asymmetric-bias note to §4.8.

With those changes I would recommend posting.
