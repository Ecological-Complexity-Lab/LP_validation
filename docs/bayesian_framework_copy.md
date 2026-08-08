# Interactive guide — review log

Working document for `docs/bayesian_framework_interactive.html`, checked against
`SI.tex` (revision of 2026-08-07). Settle text here first, implement in one pass.

**Status:** revision pass done and in the file. Awaiting section-by-section
review — comments go in Part 4 below.

---

## Part 1 — Implemented in this pass

### Previously agreed text, now in the file (removed from the to-do)

| item | where it went |
|---|---|
| Title + three-paragraph intro (*What this is / Scope / How to use*) | page header |
| "What makes a prior non-uniform?" | (i) modal on the prior card; inline `.explain` box removed |
| Short one-line slider descriptions for the three rates | error-rates card |
| "The three error rates" extended explanation | new (i) modal on the error-rates card |

Two of these needed changes before they could go in — see **F10** and **F11**.

### Errors found in the page and corrected

| # | issue | fix |
|---|---|---|
| F1 | **Ŷ used for the evidence.** The revised SI reassigns it: Ŷ is now a *truth* (the prediction a model fitted to the ground truth would give, never available) and the evidence is **Y**, the model we actually have. | Y throughout for evidence; Ŷ kept only where the SI uses it, as a component of C |
| F2 | Category named "Weak support" | "Weakly-supported" (SI Table S2) |
| F3 | Two "Box 2" references | removed; the main text has no Boxes |
| F4 | ρ defaulted to 0.5 | 0.15, the SI's illustrative value — no SI number reproduced at 0.5 |
| F5 | **Local axis symmetric in the cumulative panels.** SI S5 now says explicitly "using the directional rates of Section S4". | local axis directional below: ε_l on a link that is there, f on one that is not |
| F6 | κ asserted as (1−ε_Y)(1−ε_l) in the cumulative panels | that formula holds only under symmetric rates; κ is now computed as the actual limiting posterior |
| F7 | "Hamming distance" | "mismatched sources", the SI's phrase; added as a table column |
| F8 | Effort panel labelled "SI Section 4" | it reads Section 5; Section 4 is the grounding section |
| F9 | φ given as "the chance the interaction can occur at all", defined over O_l / O_r | φ is a **lower bound** (a link may be feasible where we never sampled), defined over the truths z_c,l and z_c,r |
| F10 | "Cross-validation returns a classifier's false-positive and false-negative rates separately" — applied to ε_Y | **wrong under the revised SI.** Cross-validation returns δ, the error of Y against the *observations*. ε_Y is the rate at which Y differs from Ŷ. Rewritten, with δ / ε_Y / δ̂ distinguished |
| F11 | κ(q) = q(1−ε_l), from the old "Model probabilities in place of a threshold" section | that section is **gone**. The score now enters as the *prior*, π_Y = q, and the two treatments are mutually exclusive. κ(q) removed; exclusivity explained in the prior modal |
| F12 | Prior box gave a generic account (sparsity, precision/recall, sampling completeness) | replaced with the SI's: π_d product form, meet-and-match conditions, connectance anchoring π_l and π_r, phylogeny into π_r, and the **double-counting** warning (SI Table S6), which is new |

### Numbers this moved

All verified against the SI text, symmetric panels and directional panels alike.

| quantity | before | now | SI |
|---|---|---|---|
| P(possibly missing), worked example | 50.4% | 50.4% | 50.4% |
| φ | 93.0% | 93.0% | 93.0% |
| κ, symmetric (upper panels) | 56% | 56% | 56% |
| κ, directional (lower panels) | 56% | **60.8%** | 60.8% |
| accumulation, n = R at R = 1 / 3 / 5 | 37.9 / 50.5 / 54.7 | **41.2 / 54.9 / 59.3** | 41.2 / 54.9 / 59.3 |
| bit summary at R = 1 / 5 / 20 | not shown | **41.2 / 39.7 / 35.4** | 41.2 / 39.7 / 35.4 |
| phantom, n = 0 at R = 10 | 36.1% | **39.2%** | 39.2% |
| R*, detections / absences | 3 / 36 | **3 / 26** | 3 / 26 |
| R* at f = 0.10 | 44 | **33** | 33 |

---

## Part 2 — Judgement calls, open to veto

**J1 — Dashed "bit" lines added to the cumulative plot.** The revised SI figure is
titled *"Confidence accumulates with replicates only when the count is
retained"*, and the count-against-bit contrast is its whole point. I added
dashed lines for the two emphasised categories only, to keep the plot readable.
This is a second series on an existing plot, not a new plot, but it is the one
thing here that goes beyond re-wording. Easy to remove.

**J2 — Two different κ values on one page.** 56% in the upper panels, 60.8% in
the lower. This is faithful: the SI defines κ as the ceiling the model and local
errors impose, and gives (1−ε_Y)(1−ε_l) only "under symmetric rates". But on a
single page it can read as a bug. Currently handled with a flagged callout in
the κ modal and a sentence in the accumulation notes. The alternative is to make
the whole page directional, which would cost the SI S2 worked example (50.4%,
φ = 93.0%) that the upper half exists to show.

**J3 — f now also sets the local false positive** in the lower panels, so moving
the f slider changes κ. This is required by the SI: R* = 33 at f = 0.10 only
comes out right this way. It is non-obvious in the UI; currently noted in the
slider description.

**J4 — ε_r slider kept.** The SI now says ε_r "is not a primitive but follows
from ε_l, f and ρ". It is still a free rate in the S1–S2 version that the upper
panels implement, so the slider stays, with a line saying the lower panels
replace it with ρ and f.

---

## Part 3 — In the SI, not on the page

Listed for a decision, not implemented.

- **Fig. S1, the error-mapping square** (δ, ε_Y, δ̂ and the ground-truth
  diagonal). Now central to SI S1. Carried in prose in the error-rates modal; a
  small inline SVG would carry it far better.
- **Table S3, the confusion matrices** composing δ with ε_l and f to reach the
  model's error against the ground truth. Absent.
- **S6 heterogeneous links** and **S7 correlated errors.** No controls, per the
  standing decision that this page does not relax assumptions. The accumulation
  modal now points to them.
- **Table S6, the factor map.** Compressed into prose in the prior modal; the
  table itself may deserve rendering.
- **S9, measuring the errors.** Absent.

**A possible inconsistency in the taxonomy tree.** Its column headers read
"Observed elsewhere" / "Not observed elsewhere" and its inner nodes read
TP/FP/TN/FN with "observed" / "not observed" — all evidence language. The leaves
are categories, which the revised SI defines as *truths*, C = (Ŷ, L_l, L_r).
The old SI blurred this; the new one does not. Worth deciding whether the tree
should be relabelled.

---

## Part 4 — Section-by-section review

Comments go here as you work through the page. Nothing implemented until we
agree the whole set and do one pass.

| # | section | comment | status |
|---|---|---|---|
| R1 | Cumulative replicate evidence | Remove the dashed "bit" lines (J1 vetoed). | **done** |
| R14 | Cumulative replicate evidence | **TODO — revisit when we work on the second half.** With the dashed lines gone, the count-vs-bit contrast survives only as prose: the notes and the modal still quote 39.7% at R = 5, 35.4% at R = 20 and "converging on half of κ", but nothing on the page shows it. Decide then whether to show it again, and how. | todo |
| R12 | New section divider | Heading before the Cumulative card marking the move to the extension. | done, then **superseded by R17** |
| R2 | κ modal + accumulation notes | Keep both κ values (J2 resolved); they just need explaining properly. Draft in Part 6. | drafting |
| R3 | Taxonomy card + modal | Tree is correct as it stands. **The eight are the categories we divide links into given our evidence — not the true categories.** The modal currently leans on the truth framing and must be reworded. Draft in Part 6. | drafting |
| R5 | Taxonomy modal | Delete the closing paragraph ("Note the distinction the framework rests on: Y is evidence, Ŷ is truth … every panel below recomputes for that category."). Folded into the R3 draft. | agreed, pending |
| R4 | New panel | SI Fig. S1 made interactive. Built as a sibling collapsible card to "Start here", titled **Mapping observations to truth**. Standalone rates local to the panel; the diagonal brightens and thickens as the model closes on the trophy. Notes in Part 7. | **done** |
|  | Title + intro |  |  |
| R6 | Prior card, hint text | Replace with: *"A uniform prior treats all eight categories as equally probable before any evidence arrives, so confidence in the chosen label depends only on the error rates — not on which category it is. Making some categories more probable than others would make the posterior no longer equal to the likelihood."* | agreed, pending — see query below |
| R7 | Prior modal | Remove the paragraph beginning "**The prior and the likelihood are not independent sources**". | agreed, pending — see note below |
| R8 | Error rates card, ε_Y slider | **Definition is wrong.** "the probability that the model's prediction Y is wrong" is δ, not ε_Y. Draft in Part 6. | agreed, pending — pick a variant |
| R9 | Error rates card, ε_r slider | "How often" → "the probability that". | done, then **superseded by R15** |
| R15 | Error rates modal **and slider**, ε_r | Drop "is wrong" from both. Use the SI's positive definition (Table S1). Modal takes variant (b), with the gloss on "available"; slider takes the short form. Text in Part 6. | **done** |
| R16 | Error rates modal, orange box | `p₁ = ρ(1−ε_l) · p₀ = f` — the `·` reads as multiplication. Comma in the rates box, semicolon in the cumulative modal (both halves already contain commas). Two other `·` kept. | **done** |
| R10 | Error rates modal | Split ε_l and f into two separate definitions instead of one joint block. Draft in Part 6. | agreed, pending |
| R11 | Error rates modal | Cut the ε_r block right down. State that the simple version uses ε_r, and that the SI relaxes that via ρ and f. Draft in Part 6. | agreed, pending — see note |
|  | Posterior bar chart + modal |  |  |
|  | Key quantities (confidence, φ) + modals |  |  |
| R17 | Page structure | Replace the two section headings with a **segmented pill**, as in the sankey explorer, switching between the simple version and the extension. Setup cards stay visible in both — see note in Part 6. | **done** |
| R18 | Page structure | Three setup cards replaced by a **sticky bar of three compact tiles** (category · prior · error rates), always visible while scrolling; clicking one expands its panel. Notes in Part 10. | **done** |
| R19 | Intro | Moved into a collapsible card titled **Start here — what this is, what it assumes, and how to use it**. Collapsed by default. Note in Part 10. | **done** |
|  | Likelihood table + modal |  |  |
|  | Effort to confidence + modal |  |  |

---

## Part 5 — Agreed text, as implemented

Kept for reference; edits to these go in Part 4.

### Title

> **How confident are we in a *possibly missing* link?**

The italicised category updates with the selection.

### Intro

> **What this is:** an interactive guide to the Bayesian framework developed in
> the Supplementary Information of *Contextual evidence for classifying
> interactions and guiding their discovery*. The taxonomy assigns each candidate
> link to one of eight categories from three pieces of evidence — the model
> prediction Y, the local observation O_l, and the replicate observation O_r.
> That assignment is deterministic, but every source is error-prone, so it
> carries uncertainty. The framework therefore treats the true category C as
> unknown and returns a posterior P(C | E) over the eight.
>
> **Scope:** the upper panels follow the simplest version of the framework
> (SI S1–S2) — binary evidence, one symmetric error rate per source, a uniform
> prior, and independent errors. The two lower panels follow the grounded and
> cumulative version (SI S4–S5), where the local axis errs at different rates in
> each direction and the replicate count is retained. Throughout, the
> illustrative rates are those of the SI (ε_Y = 0.2, ε_l = 0.3, ε_r = 0.1,
> ρ = 0.15, f = 0.05), chosen to show the behaviour of the framework rather than
> to describe any particular system.
>
> **How to use:** pick a category on the taxonomy tree. Its error-free signature
> z_c becomes the observed evidence E, and the confidence in that label is
> (1−ε_Y)(1−ε_l)(1−ε_r). Drag the sliders to see how each error rate moves it,
> and which rival categories absorb the rest.

Changed from the version we agreed: "the model prediction Ŷ" → "Y" (F1), and
the *Scope* paragraph now names the two rate regimes rather than one (F5).

### Slider descriptions

> **ε_Y model error** — the probability that the model's prediction Y is wrong.
>
> **ε_l local error** — the miss rate: the probability that a link realised
> where we sampled goes unrecorded in O_l.
>
> **ε_r replicate error** — how often the replicate evidence O_r — recorded in
> at least one replicate — is wrong.

The sentence we had discussed adding to the ε_l description, naming f, was
dropped: f does not exist at that point on the page. It is introduced properly
in the error-rates modal and again on the cumulative card.

---

## Part 6 — Drafts pending agreement

### R6 — prior card hint (agreed text)

> A uniform prior treats all eight categories as equally probable before any
> evidence arrives, so confidence in the chosen label depends only on the error
> rates — not on which category it is. Making some categories more probable than
> others would make the posterior no longer equal to the likelihood.

**Open query.** This drops the operational sentence "Weights are relative,
rescaled to sum to 1 (shown as %)". Without it the eight sliders are unexplained:
dragging one changes every percentage, which reads as a bug. Options — append it
to the paragraph, demote it to small print under the slider grid, or drop it and
let the behaviour speak for itself.

### R7 — note on what the removed paragraph carried

The deleted paragraph was the SI's **double-counting** argument (S8 and Table
S6): abundance, co-occurrence and generalism raise the prior *and* lower
ε_l / raise ρ, so multiplying prior by likelihood counts the same ecology twice
and the posterior is overconfident; relatedness runs the other way, raising π_r
while also raising f. Nothing else on the page carries this, so with R7 applied
it leaves the guide entirely. Flagged, not contested — the modal is long and
this is the most technical part of it.


### R16 — the `·` separator between p₁ and p₀

**Fix (rates modal, orange box):**

> p₁ = ρ(1−ε_l),  p₀ = f

**Same fault, cumulative modal.** It separates the same two definitions the same
way:

> p₁ = ρ(1−ε_l) (realised, then detected) · p₀ = f (a false detection of a link
> that is not there)

Both halves already contain commas, so a semicolon reads better there:

> p₁ = ρ(1−ε_l) (realised, then detected);  p₀ = f (a false detection of a link
> that is not there)

**The effort modal already uses a comma** for the same pair, so this is an
internal inconsistency as well as a misreading risk.

**Two other `·` uses, proposed to keep.** Neither sits between equations, so
neither can be read as multiplication:

- rates modal, separating three labelled definitions — *δ: Y against the
  observations · ε_Y: Y against Ŷ · δ̂: Ŷ against the ground truth*
- the live readout under the confidence figure — *likelihood 50.4% · prior,
  normalised ⇒ P(C|E) = 50.4%* — where the `·` genuinely does mean "times"

### R15 — ε_r without "is wrong"

SI Table S1 defines it positively, and directionally: *"Probability that an
available link is recorded in none of the replicates."* Using the SI's own
wording removes "is wrong" and makes all three blocks parallel — ε_l already
says "goes unrecorded", not "is wrong".

**Modal block (agreed — variant b):**

> **ε_r — the replicates.** The probability that an available link — one that
> could be realised in the replicate contexts — is recorded in none of them.

The gloss on "available" earns its place because L_r is never spelt out
elsewhere on the page: the taxonomy modal now stops before it, and the rates
modal defines L_l only.

**Slider caption (agreed — "is wrong" removed here too):**

> the probability that an available link is recorded in none of the replicates

Short form, since the (i) button beside it opens the modal that glosses
"available". Note this changes what the caption describes: the old one was about
O_r = 1 ("recorded in at least one replicate"), the new one is about the error
on a link that is there. The second is the SI's definition. This supersedes R9.

### R8 — ε_Y slider caption

The current caption states δ, not ε_Y. Per SI Table S1, ε_Y is the rate at which
Y differs from Ŷ — the prediction a model fitted to the ground truth would give.

> **(a)** how often the prediction would change had the model been fitted to the
> ground truth
>
> **(b)** how often the prediction would change had the model been fitted to the
> ground truth — not the error cross-validation reports

(b) recommended: the whole point of the modal is that readers will assume this
number is the cross-validated one, and the caption is where that assumption
forms.

### R9 — ε_r slider caption

> the probability that the replicate evidence O_r — recorded in at least one
> replicate — is wrong

### R10 + R11 — error-rates modal, restructured

Body covers only the three rates of the simple version, one short block each;
everything grounded goes into an orange `.flag` box at the end. ε_l and f are
separated by being in different places, which is the cleanest split available.

> Every source of evidence is error-prone, and the three rates are where that
> enters the framework. They carry two kinds of uncertainty: the **model**
> (ε_Y) and the **observations** (ε_l and ε_r). In this simplest version each
> source has a single rate, applied whichever way the error runs.
>
> **ε_Y — the model.** This is *not* what cross-validation reports.
> Cross-validation returns **δ = 1 − accuracy**, the rate at which the
> prediction Y departs from the observations it was fitted to. ε_Y is a
> different quantity: the rate at which Y differs from Ŷ, the prediction that
> fitting to the ground truth would have given.
>
> > δ: Y against the observations (measurable) · ε_Y: Y against Ŷ ·
> > δ̂: Ŷ against the ground truth (never available)
>
> Had the observations been error-free, the fitted model would *be* the
> truth-fitted model and ε_Y would be zero. So what makes ε_Y non-zero is the
> same observation error the framework already carries in ε_l — which no
> held-out set can reveal, because the withheld labels are observations too. In
> practice δ is the stand-in, and how it is obtained matters: a random hold-out
> gives the error on pairs of the same kind, whereas blocking by site or by
> clade gives the error to expect when the model is carried to a new context,
> which is nearer to what the framework needs.
>
> **ε_l — the local observation.** The miss rate: the probability that a link
> realised where we sampled goes unrecorded. It concerns detection alone — local
> realisation is not something marginalised over here, it *is* the latent state
> L_l. In this version the same number also serves as the rate at which a link
> that is not there gets recorded anyway.
>
> **ε_r — the replicates.** The probability that the replicate evidence —
> recorded in at least one replicate — is wrong.
>
> ---
> **[orange box]**
>
> **ε_r is not a primitive.** Unlike ε_Y and ε_l it is not a property of a
> method but a consequence of three other quantities. A replicate is only the
> same method applied elsewhere, so the replicate axis adds no sampling
> parameter of its own. What it adds is the **realisation rate ρ**, the
> probability that a feasible link is actually realised in a given replicate — an
> ecological property rather than a methodological one. A record then needs both
> steps, realisation and then detection:
>
> > p₁ = ρ(1−ε_l)  ·  p₀ = f
>
> where **f** is the **false-detection rate**, the probability that a link that
> is not there is recorded anyway through misidentification or contamination. f
> is the other direction of the local error, which this simplest version equates
> with ε_l. From p₁ and p₀, ε_r follows for any number of replicates R. The two
> cumulative panels below take that route, replacing ε_r with ρ and f.

**Still homeless** (flagged once before, not requested back): SI Eq. S3 and
Fig. S2 — ε_r⁺(R) = (1−p₁)^R falls as replicates accumulate while
ε_r⁻(R) = 1−(1−f)^R rises, so no R makes both small; at the illustrative values
they cross at R ≈ 9 where both equal 0.37, far above the 0.1 the slider starts
at. One sentence at the foot of the orange box would carry it if wanted.

### R12 — section divider before the cumulative card

Sits between cards, on the page background rather than inside a card, so it
reads as a break in the page rather than as another panel. Covers both lower
cards, not just the first.

**Title (agreed — pair D):**

> ## Extension: many replicates, and errors that run two ways

Lede beneath it:

> The two panels below relax two of the assumptions above. The replicate
> evidence keeps its full count n rather than collapsing to a single bit, and
> the local axis errs at a different rate in each direction: ε_l on a link that
> is there, f on one that is not. The numbers therefore differ from the worked
> example above, and so does κ.

**Decided.** Keep both the lede and the orange flag box in the cumulative modal.
Add a matching heading above the simple half.

### R13 — matching heading above the simple half

**Placement problem.** The heading cannot simply go after the intro. The first
three cards feed *both* halves: the taxonomy sets the evidence for everything,
the prior is used by every posterior on the page, and ε_Y and ε_l drive both
halves — only ε_r is confined to the simple one. So the page is really three
sections, not two:

| section | cards |
|---|---|
| setup | taxonomy · prior · error rates |
| the simple version | bar chart + key quantities · likelihood table |
| the extension | cumulative replicate evidence · effort to confidence |

So the "simple version" heading belongs **before the bar-chart row**, not before
the taxonomy card. Optionally a third heading over the setup cards.

**Title (agreed — pair D):**

> ## The simple version: one rate per source

Placed immediately before the bar-chart row. The three setup cards above it stay
unlabelled, so the two headings mark exactly the analytical contrast — say if a
third "Setup" heading is wanted after all.

Lede for the simple half:

> Each source is reduced to a single bit and errs at one rate whichever way the
> error runs, the prior is uniform, and the three sources err independently.
> Under those assumptions the posterior equals the likelihood, and confidence in
> the matching category is just (1−ε_Y)(1−ε_l)(1−ε_r).

Replaces the flag box in the κ modal.

> **Why κ is larger in the two panels below.** Those panels let the local axis
> err at a different rate in each direction: a link that is there is missed at
> ε_l = 0.3, but a link that is not there is recorded at only f = 0.05. Misses
> are common; spurious records are rare.
>
> That asymmetry changes what a *non-record* is worth. Under symmetric rates,
> "the link really is absent" and "the link is there but was missed" are treated
> as equally plausible readings of O_l = 0, so the matching category keeps a
> factor of 1−ε_l = 0.7 against a rival's 0.3. Under directional rates it keeps
> 1−f = 0.95 against the same 0.3, so it takes a larger share once the replicate
> axis is resolved: κ rises from 56% to 60.8%.
>
> A second consequence: **κ is no longer the same for every category.** Under
> symmetric rates every category has the same ceiling, which is why confidence
> depends only on the error rates and not on which category you picked. Under
> directional rates it depends on the local bit — for the four categories
> observed locally the ceiling is 74.7%, because a record is strong evidence
> when false detections are rare, while for the four not observed locally it is
> 60.8%. Click between leaves and the dashed line moves.
>
> κ is the same idea throughout: the limit the model and local errors impose on
> the posterior once replication has done all it can. Only the rate structure
> differs.

Verified: κ_dir is 74.7% for the four categories with z_c,l = 1 and 60.8% for
the four with z_c,l = 0; κ_sym is 56% for all eight.

### R3 — taxonomy modal, reworded

The tree is the division we make **given the evidence**, not a claim about the
true categories. Replaces the current taxonomy modal.

> Every candidate link is sorted by three binary pieces of evidence, collected
> as **E = (Y, O_l, O_r)**: the prediction Y of the model we have, the local
> observation O_l (1 if the link was recorded where we sampled), and the
> replicate observation O_r (1 if it was recorded in at least one within-system
> replicate). Every quantity refers to a single candidate link, whose taxon and
> context indices are suppressed throughout.
>
> The first two define the confusion matrix — true/false positive/negative. The
> contextual axis then splits each cell in two ("not observed elsewhere" vs
> "observed elsewhere"), giving the eight categories of the tree.
>
> **This is the division we make given our evidence, not the true state of the
> links.** The tree sorts on what was observed. The bins carry the eight
> category names because there are exactly eight error-free evidence patterns,
> and they correspond one to one with the eight categories. We write
> **z_c = (z_c,Y, z_c,l, z_c,r)** for the pattern category c would produce if no
> source erred — its **signature**.
>
> Because every source is error-prone, the bin a link falls into need not be its
> true category: a *possibly missing* link may be a *recurrent* one that local
> sampling missed. That is what the rest of the page computes — the posterior
> over the eight true categories, given the evidence that put the link where it
> is.
>
*(R5: the closing paragraph on Y-vs-Ŷ, and the "click any leaf" sentence with
it, is cut. The distinction survives in the error-rates modal, which is where
δ, ε_Y and δ̂ are set out; the click instruction survives in the card's summary
line and in the intro's "How to use".)*

---

## Part 8 — TODO: showing count against bit again (R14)

Parked at the same time as R1, to revisit with the second half of the page.

**What the prose still claims, unillustrated.** Retaining n carries *possibly
missing* to κ = 60.8%; collapsing to the bit O_r = 1[n ≥ 1] sends the same
replicates the other way, to 39.7% at R = 5 and 35.4% at R = 20, converging on
κ/2 = 30.4%. The two agree at R = 1 and coincide throughout under n = 0.

**What it would take to show it again.** The removed machinery was:

- `qPostBit(yh, ol, obsR, R)` — same as `qPost` but with the replicate factor
  replaced by `seen = 1 − (1−p_c)^R`, taken as `seen` when O_r = 1 and
  `1 − seen` when O_r = 0. Deleted; recoverable from git history.
- a `bitSeries` array filled alongside `series` in `drawAccumAll`, using
  `qPost(y, l, 0, 0)` at R = 0 where the bit is undefined
- a `drawBit(i)` dashed stroke, drawn for the two emphasised categories only
- legend and caption text for solid-vs-dashed

**Options when we return to it.** Restore the dashed overlay as before; or give
the bit its own small panel rather than overlaying; or drop the claim from the
prose so the page only asserts what it shows.

---

## Part 7 — R4: the error-mapping square, as built

**Placement (A).** Its own collapsible card, sibling to "Start here", collapsed
by default:

> ▸ **Mapping observations to truth** — why cross-validation flatters the model

**Rates are standalone**, local to the panel, changing nothing else on the page.
δ is not a model parameter anywhere else, and the square's diagonal (Y against
the ground truth) is a different quantity from ε_Y (Y against Ŷ), so coupling it
to the page would have implied a link that does not exist.

**Sliders sit beside their edges**: δ above the top edge (Y → O_l), ε_l and f to
the right of the right edge (O_l → L_l). The left and bottom edges — ε_Y and δ̂ —
are drawn faint and dashed with "not measurable" labels, since neither follows
from δ, ε_l and f.

**Each value appears once.** The first build printed every rate twice, on the
slider and again on the diagram. Now the control carries the symbol and its
name, and the diagram's edge label is the only readout — so reading the value
and reading what it does happen in the same place. The δ slider is short and
centred over its edge rather than stretched across the column.

**The diagonal responds.** Its colour runs grey → green and its width 1.7 → 4.6
with `t = (bal − 0.5)/0.5`, where `bal = 1 − (FN + FP)/2`; the trophy's opacity
runs 0.2 → 1 on the same scale. So dragging any rate visibly moves the model
toward or away from the truth.


An interactive version of SI Fig. S1, driven by SI Table S3. The point it makes:
**what cross-validation reports is not the error that matters**, and no amount
of improving either the model or the detection alone will close the gap.

### The composition

Table S3's lower block gives the model's error against the ground truth:

```
FN of Y against the truth  =  (1−ε_l)·δ + ε_l·(1−δ)
FP of Y against the truth  =  f·(1−δ) + (1−f)·δ
```

### Layout

```
  What cross-validation actually tells you                            (i)

     Y ───────────── δ ─────────────▶ O_l          δ   model vs observations
     │  ╲                              │               ▸────────────  0.20
     │    ╲                            │               what CV reports
    ε_Y     ╲  the error that          ε_l , f
     │        ╲   matters  🏆          │           ε_l miss rate
     │          ╲                      │               ▸──────  0.30
     ▼            ╲                    ▼
     Ŷ ───────────── δ̂ ────────────▶ L_l           f   false detection
                                                       ▸──  0.05

   Y against the ground truth        apparent accuracy   80.0%   (= 1−δ)
     misses  38.0%   invents  23.0%  true accuracy       69.5%

   [ perfect model δ=0 ]  [ perfect detection ε_l=f=0 ]  [ both ]
```

Edge weight scales with its rate, so a thick edge reads as error leaking
through. The diagonal is the trophy: it lights up only when both blocks are
clean.

### What the three presets show

| preset | FN | FP | reading |
|---|---|---|---|
| perfect model, δ = 0 | ε_l = 0.30 | f = 0.05 | recovers the *observations*, inheriting detection error in full |
| perfect detection, ε_l = f = 0 | δ = 0.20 | δ = 0.20 | leaves the model error untouched |
| both | 0 | 0 | the identity — the trophy |

That is the SI caption's argument made draggable: neither block alone reaches
the truth.

### Verified against SI Table S3

| preset | misses | invents | real accuracy | SI's claim |
|---|---|---|---|---|
| δ = 0.20, ε_l = 0.30, f = 0.05 | 0.380 | 0.230 | 69.5% | vs 80% apparent |
| perfect model, δ = 0 | 0.300 | 0.050 | 82.5% | recovers the upper block: FN = ε_l, FP = f |
| perfect detection, ε_l = f = 0 | 0.200 | 0.200 | 80.0% | leaves the model error alone: FN = FP = δ |
| both | 0 | 0 | 100% | the identity |

The middle two rows are the caption's argument, and they come out exactly as it
says: neither block alone reaches the truth.


---

## Part 9 — R17: the version pill, and the setup question

Implemented as a `.seg-wrap` segmented control, the same component the scenario
toggle already uses and the same one the sankey explorer uses for "Sampling
method". It replaces both section headings from R12/R13; the two ledes survive
as the pill's caption, swapped on switch.

    [ The simple version ] [ Extension: many replicates ]

**The setup was deliberately not separated.** The taxonomy, the prior and
ε_Y / ε_l all feed *both* versions — only ε_r is confined to the simple one. Put
the setup behind its own tab and you cannot watch a posterior move while
dragging the rate that moves it, which is the point of the page. The error-rates
bar is sticky for exactly that reason.

If a three-tab split is still wanted, the workable version is different: keep
the rates bar pinned outside the tabs and put only the taxonomy and prior cards
behind a "Setup" tab, since those are set once and then left alone. Both are
already collapsible, which achieves most of the same tidying without the cost.

**A follow-up this opens up.** On the extension tab the ε_r slider is inert —
the extension replaces it with ρ and f. It could be dimmed or hidden there. Not
done, since it changes behaviour rather than wording.


---

## Part 10 — R18: the compact setup bar

The taxonomy, prior and error-rate cards are now three tiles in one sticky bar.
Each tile reports its own current state, so the setup is legible at a glance
without opening anything:

    · CATEGORY            · PRIOR        · ERROR RATES
      Possibly missing      uniform        ε_Y 0.20  ε_l 0.30  ε_r 0.10

Tile values update live — the category name follows the tree selection, the
prior reads "uniform" or "custom", and the rates track their sliders.

### Choices made

- **One panel at a time**, and clicking the open tile closes it. An accordion
  keeps the sticky region bounded; independent toggles could stack all three and
  fill the viewport.
- **Panels capped at 47vh with internal scroll**, so an open panel can never
  swallow the page, whatever the window height.
- **Category open on arrival**, so the taxonomy tree is still the first thing
  seen. It collapses on click and stays collapsed.
- **Taxonomy SVG capped at 840px wide**, which keeps the open panel about 245px
  tall instead of 335px.
- The old `<details class="card coll">` cards are gone, and the CSS that styled
  them was removed with them.

### Consequences

- The rates panel is no longer permanently on screen. Adjusting a rate while
  watching a plot now takes one click to reopen it, though the tile keeps the
  current values visible throughout.
- The old sticky `ratesbar` is replaced by `setupbar`; `layoutBars()` and
  `barShadows()` now track that instead.

### R19 — the intro card

The three intro paragraphs now sit in a collapsible card:

> ▸ **Start here** — what this is, what it assumes, and how to use it

**Collapsed by default**, so the taxonomy tree and the first plot are reachable
without scrolling past ~25 lines of prose. "Start here" is meant to carry the
invitation for a first-time reader; if that feels too much of a gamble, one
character flips it (`<details ... open>`).

The lean `details.foldcard` CSS replaces the `details.coll` rules the setup bar
made redundant, so this cost no net stylesheet weight.

Alternative titles considered: "Read me first", "How to read this page", "The
short version", "What you're looking at".
