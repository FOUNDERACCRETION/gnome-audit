#!/usr/bin/env python3
"""Fail if a superseded number is still anywhere in the manuscript source.

check_paper.py parses the LaTeX tables and diffs them against results/. It cannot see
numbers in prose, and across this paper's revisions the same failure has recurred five
times: a quantity is recomputed, the table is regenerated, and a sentence somewhere still
quotes the old value. Appendix A kept the pre-convergence headline intervals for a whole
draft that way.

This is the cheap complement. Every time a number changes, its superseded literal goes in
the table below with what replaced it. The check then greps the whole source --- prose,
captions, appendices --- and fails on any hit. It cannot catch a number that was never
right, only one that used to be right and no longer is, which is precisely the class the
table checker misses.

    python scripts/check_stale.py [paper/gnome_audit_paper.tex]

Exit status 0 if no superseded literal is present, 1 otherwise.
"""
import re, sys
TEX = sys.argv[1] if len(sys.argv) > 1 else "paper/gnome_audit_paper.tex"

# (superseded literal, what replaced it, when and why it changed, allowed occurrences)
#
# The allowance exists because a few superseded values are quoted deliberately, to show
# what a correction changed. Those are historical citations, not stale text, and the count
# is pinned so that an accidental reappearance still fails.
STALE = [
    # Draft 19: bootstrap intervals recomputed at converged replicate counts
    ("[27.6, 33.8]",  "[27.7, 33.8]",  "D19: cell intervals 600 -> 20,000 replicates"),
    ("[3.2, 4.8]",    "[3.3, 4.7]",    "D19: cell intervals 600 -> 20,000 replicates"),
    ("[40.2, 46.2]",  "[39.8, 46.0]",  "D19: per-element 300 -> 20,000 replicates"),
    ("[30.0, 35.6]",  "[30.5, 35.7]",  "D19: per-element 300 -> 20,000 replicates"),
    ("[25.6, 36.8]",  "[26.2, 35.8]",  "D19: per-element 300 -> 20,000 replicates"),
    ("[22.8, 27.9]",  "[22.6, 27.9]",  "D19: per-element 300 -> 20,000 replicates"),
    ("[+0.5$, $+9.6", "[+0.5$, $+10.1","D19: within-anion excesses 400 -> 100,000"),
    ("[+35.6$, $+46.5","[+35.4$, $+47.0","D19: within-anion excesses 400 -> 100,000"),
    # Draft 20: average marginal effects at 1,000 refits
    ("[+14.3$, $+16.6","[+14.2$, $+16.6","D20: AME 150 -> 1,000 refits"),
    ("[$-11.2$, $-4.9$]","[$-11.4$, $-4.8$]","D20: AME 150 -> 1,000 refits"),
    ("[$-10.9$, $-4.6$]","[$-10.7$, $-4.9$]","D20: AME 150 -> 1,000 refits"),
    ("[$+6.4$, $+11.2$]","[$+6.5$, $+11.6$]","D20: AME 150 -> 1,000 refits"),
    # Draft 21: the matched comparison became a GNoME-side cluster bootstrap
    ("[$+11.9$, $+22.5$]","[$+11.4$, $+23.0$]","D21: matched comparison clustered"),
    ("[$-8.1$, $-0.3$]", "[$-9.3$, $+1.9$]",  "D21: matched comparison clustered"),
    ("[$+0.49$, $+1.95$]","[$+0.22$, $+2.28$]","D21: matched comparison clustered", 1),
    ("[$+0.76$, $+2.87$]","[$+0.57$, $+3.11$]","D21: matched comparison clustered", 1),
    ("[$+0.48$, $+2.14$]","[$+0.18$, $+2.51$]","D21: matched comparison clustered"),
    ("[$+1.07$, $+3.36$]","[$+0.83$, $+3.66$]","D21: matched comparison clustered"),
    # phrasings that were corrected and must not come back
    ("positivity failure", "severe estimated lack of overlap",
     "D21: fitted probabilities cannot establish population positivity"),
    ("600 replicates for cells", "20,000 replicates for cells",
     "D19: stale replicate count"),
    ("51 remain so under every bootstrap", "the same 53 survive under every",
     "D21: stale seed-stability count"),
    ("no such choice would close it", "remains large under the sensitivity checks performed",
     "D21: untested clustering guarantee"),
    ("we lose nothing by using it as one", "did not improve held-out catch materially",
     "D21: overstated flexible-hull baseline"),
    ("the eleven of Table", "the ten of Table",
     "D21: clustered matched comparison, 11 -> 10 intervals excluding zero"),
    ("[$+0.3$, $+3.8$]", "[$+0.3$, $+3.7$]",
     "D23: Ni intermetallic within-anion excess, prose quoted an upper endpoint of 3.8 while "
     "tab:feexcess and within_anion_excess_ci.csv both give 3.7"),
]

src = open(TEX).read()
flat = src.replace("{,}", ",")
bad = []
for entry in STALE:
    stale, replacement, why = entry[0], entry[1], entry[2]
    allowed = entry[3] if len(entry) > 3 else 0
    n = max(src.count(stale), flat.count(stale))
    if n > allowed:
        hay = src if src.count(stale) >= flat.count(stale) else flat
        ln = hay[:hay.index(stale)].count("\n") + 1
        bad.append((stale, replacement, why, ln, n, allowed))
print("checked %d superseded literals against %s" % (len(STALE), TEX))
if bad:
    print("\n%d STALE VALUE(S) STILL IN THE SOURCE:" % len(bad))
    for stale, repl, why, ln, n, allowed in bad:
        extra = "" if allowed == 0 else "  (%d found, %d allowed as history)" % (n, allowed)
        print("  line %-5d %-34s should now be %-24s  [%s]%s"
              % (ln, stale, repl, why, extra))
    sys.exit(1)
print("no superseded value found in prose, captions or tables")
